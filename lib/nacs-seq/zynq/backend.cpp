/*************************************************************************
 *   Copyright (c) 2021 - 2021 Yichao Yu <yyc1992@gmail.com>             *
 *                                                                       *
 *   This library is free software; you can redistribute it and/or       *
 *   modify it under the terms of the GNU Lesser General Public          *
 *   License as published by the Free Software Foundation; either        *
 *   version 3.0 of the License, or (at your option) any later version.  *
 *                                                                       *
 *   This library is distributed in the hope that it will be useful,     *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of      *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU    *
 *   Lesser General Public License for more details.                     *
 *                                                                       *
 *   You should have received a copy of the GNU Lesser General Public    *
 *   License along with this library. If not,                            *
 *   see <http://www.gnu.org/licenses/>.                                 *
 *************************************************************************/

#include "backend.h"

#include "../../nacs-utils/llvm/utils.h"
#include "../../nacs-utils/llvm/codegen.h"
#include "../../nacs-utils/llvm/compile.h"
#include "../../nacs-utils/llvm/execute.h"
#include "../../nacs-utils/llvm/global_rename.h"
#include "../../nacs-utils/processor.h"

#include <llvm/ADT/StringRef.h>
#include <llvm/IR/LegacyPassManager.h>
#include <llvm/IR/Verifier.h>
#include <llvm/Transforms/IPO.h>
#include <llvm/Transforms/IPO/AlwaysInliner.h>

namespace NaCs::Seq::Zynq {

static Device::Register<Backend> register_backend("zynq");

struct Backend::BasicSeq {
    std::unique_ptr<BCGen> bc_gen;
    std::vector<uint8_t> bytecode;
    std::vector<uint8_t> bytecode_first;
    bool clock_set = false; // Whether clock is set at generation time.
    bool all_known = false;
};

struct Backend::TTLManager {
    // All time in sequence time units
    uint32_t off_delay;
    uint32_t on_delay;
    uint32_t skip_time;
    uint32_t min_time;

    bool off_val;
};

// This is the FPGA clock frequency (100 MHz), which is hard coded for now.
static constexpr uint64_t fpga_tick_per_sec = 100ul * 1000ul * 1000ul;

Backend::Backend(Manager &mgr, std::string name)
    : Device(mgr, std::move(name))
{
    m_seq_tick_per_sec = mgr.tick_per_sec();
    if (m_seq_tick_per_sec % fpga_tick_per_sec != 0) {
        throw std::runtime_error(name + ": Sequence time tick must be a fraction of 10ns.");
    }
}

Backend::~Backend()
{
    if (m_obj_id) {
        mgr().exe_engine().free(m_obj_id);
    }
}

NACS_EXPORT() void Backend::set_clock_active_time(
    uint32_t bseq_idx, uint64_t half_period,
    llvm::ArrayRef<std::pair<int64_t,int64_t>> active_times)
{
    if (!m_has_clock)
        throw std::runtime_error(name() + ": Clock output not enabled");
    const auto fpga_clock_div = uint32_t(m_seq_tick_per_sec / fpga_tick_per_sec);
    if (half_period % fpga_clock_div != 0)
        throw std::runtime_error(name() + ": Half clock period tick (" +
                                 std::to_string(half_period * 1e9 / m_seq_tick_per_sec) +
                                 " ns) must be a fraction of 10ns.");
    if (half_period == 0 && !active_times.empty())
        throw std::runtime_error(name() + ": Clock period cannot be zero.");
    if (half_period > fpga_clock_div * 255)
        throw std::runtime_error(name() + ": Clock period (" +
                                 std::to_string(half_period * 2e9 / m_seq_tick_per_sec) +
                                 " ns) exceeds 5.1 us.");
    auto clock_div = uint8_t(half_period / fpga_clock_div);
    std::vector<BCGen::Clock> clocks;
    for (auto [first, last]: active_times) {
        // Start the clock half cycle earlier so that the falling edge is on time.
        first = (first * 2 - 1) * half_period;
        // Stop at the expected time of the falling edge to make sure we turn it off
        // before any raising edges.
        last = last * 2 * half_period;
        clocks.push_back({ .time = first, .period = clock_div });
        clocks.push_back({ .time = last, .period = 0 });
    }
    if (!m_generated) {
        if (!m_clocks.emplace(bseq_idx, std::move(clocks)).second) {
            throw std::runtime_error(name() + ": Clock for basic sequence [index " +
                                     std::to_string(bseq_idx) + "] already set");
        }
    }
    else {
        if (bseq_idx >= m_seqs.size())
            throw std::runtime_error(name() + ": Invalid basic sequence index " +
                                     std::to_string(bseq_idx));
        auto &be_bseq = m_seqs[bseq_idx];
        if (!be_bseq.bc_gen || be_bseq.clock_set)
            throw std::runtime_error(name() + ": Clock for basic sequence [index " +
                                     std::to_string(bseq_idx) + "] already set");
        if (unlikely(mgr().dump_enabled()))
            m_clocks.emplace(bseq_idx, clocks);
        be_bseq.bc_gen->clocks = std::move(clocks);
    }
}

void Backend::add_channel(uint32_t chn_id, const std::string &_chn_name)
{
    llvm::StringRef chn_name(_chn_name);
    if (chn_name.startswith("TTL")) {
        auto chn_num_str = chn_name.substr(3);
        uint8_t chn_num;
        if (chn_num_str.getAsInteger(10, chn_num) || chn_num >= 32)
            throw std::runtime_error(name() + ": Invalid TTL channel " + _chn_name);
        if (chn_num == m_start_ttl_chn)
            throw std::runtime_error(name() + ": Cannot output on start trigger TTL channel");
        m_chn_map.emplace(chn_id, std::make_pair(BCGen::ChnType::TTL, chn_num));
        m_ttl_mask = setBit(m_ttl_mask, chn_num, true);
    }
    else if (chn_name.startswith("DAC")) {
        auto chn_num_str = chn_name.substr(3);
        uint8_t chn_num;
        if (chn_num_str.getAsInteger(10, chn_num) || chn_num >= 4)
            throw std::runtime_error(name() + ": Invalid DAC channel " + _chn_name);
        m_chn_map.emplace(chn_id, std::make_pair(BCGen::ChnType::DAC, chn_num));
    }
    else if (chn_name.startswith("DDS")) {
        auto subname = chn_name.substr(3);
        uint8_t chn_num;
        if (subname.consumeInteger(10, chn_num) || chn_num >= 22)
            throw std::runtime_error(name() + ": Invalid DDS channel " + _chn_name);
        if (subname == "/FREQ") {
            m_chn_map.emplace(chn_id, std::make_pair(BCGen::ChnType::Freq, chn_num));
        }
        else if (subname == "/AMP") {
            m_chn_map.emplace(chn_id, std::make_pair(BCGen::ChnType::Amp, chn_num));
        }
        else {
            throw std::runtime_error(name() + ": Invalid DDS parameter name " + subname.str());
        }
    }
    else {
        throw std::runtime_error(name() + ": Unknown channel name (" +
                                 _chn_name + ") for Zynq backend.");
    }
}

bool Backend::check_noramp(uint32_t chn_id, const std::string &chn_name)
{
    return llvm::StringRef(chn_name).startswith("TTL");
}

static unsigned get_vector_size()
{
    // This needs to be consistent with the bc_gen data stream.
#if !defined(ENABLE_SIMD)
    return 1;
#elif NACS_CPU_X86 || NACS_CPU_X86_64
    static unsigned vector_size = [&] {
        auto &host = CPUInfo::get_host();
        if (host.test_feature(X86::Feature::avx512f) &&
            host.test_feature(X86::Feature::avx512dq)) {
            return 8;
        }
        else if (host.test_feature(X86::Feature::avx)) {
            return 4;
        }
        else {
            return 2;
        }
    }();
    return vector_size;
#elif NACS_CPU_AARCH64
    // TODO SVE ...
    return 2;
#else
    return 1;
#endif
}

void Backend::generate(Manager::ExpSeq &expseq, Compiler &compiler)
{
    mgr().add_debug_printf("%s: compiling\n", name().c_str());
    // Here we create (and later run) the sequence on the FPGA even if no pulses are needed.
    // This is because we also provides the start trigger function for other devices.
    // If we allow multiple FPGA device (where there's one master providing start trigger
    // and maybe clock for other ones),
    // we can expose more control on start trigger and clock output
    // so that each device can know more precisely if they are needed in a basic sequence
    // and may skip the ones where they are completely unused.
    const auto fpga_clock_div = uint32_t(m_seq_tick_per_sec / fpga_tick_per_sec);
    const auto bseqs = compiler.basic_seqs();
    auto nbseqs = bseqs.size();
    m_seqs.resize(nbseqs);
    for (auto &[idx, clocks]: m_clocks) {
        if (idx >= nbseqs) {
            throw std::runtime_error(name() + ": Invalid basic sequence index " +
                                     std::to_string(idx));
        }
    }

    auto &env = compiler.seq()->env();
    auto env_mod = env.llvm_module();
    auto &llvm_ctx = env_mod->getContext();
    llvm::Module mod("", llvm_ctx);
    // For copying code into the new module
    LLVM::FunctionMover mover(&mod);
    // A new cgctx that emits code in the new module
    LLVM::Codegen::Context cgctx(&mod);

    // Do a scan to collect and fix all functions so that we can run
    // always inline pass before trying to vectorize all the functions.
    std::map<llvm::Function*,llvm::Function*> func_map;

    for (uint32_t idx = 0; idx < nbseqs; idx++) {
        auto &bseq = bseqs[idx];
        auto &chn_infos = bseq->get_channel_infos();
        for (auto [chn_id, fpga_chn]: m_chn_map) {
            auto it = chn_infos.find(chn_id);
            assert(it != chn_infos.end());
            auto &chn_info = it->second;
            for (auto &pulse: chn_info.pulses) {
                if (pulse.is_measure())
                    continue;
                auto len = pulse.len();
                if (len) {
                    auto val = pulse.val();
                    assert(val->is_call());
                    auto f0 = val->get_callee().llvm;
                    auto &f = func_map[f0];
                    // Reuse functions
                    if (f)
                        continue;
                    f = mover.clone_function(f0);
                    f = Compiler::get_ramp_func(val, f, cgctx.T_f64, cgctx.T_f64, false);
                }
            }
        }
    }

    // Clean up function exports
    for (auto &go: mod.global_objects()) {
        if (!go.isDeclaration()) {
            go.setLinkage(llvm::GlobalValue::PrivateLinkage);
        }
    }
    for (auto [f0, f]: func_map) {
        f->setLinkage(llvm::GlobalValue::ExternalLinkage);
        f->setVisibility(llvm::GlobalValue::ProtectedVisibility);
        f->setUnnamedAddr(llvm::GlobalValue::UnnamedAddr::Global);
    }
    {
        llvm::legacy::PassManager PM0;
#ifndef NDEBUG
        PM0.add(llvm::createVerifierPass());
#endif
        PM0.add(llvm::createAlwaysInlinerLegacyPass());
#ifndef NDEBUG
        PM0.add(llvm::createVerifierPass());
#endif
        PM0.run(mod);
    }

    struct RampFunc {
        uint32_t bseq_idx;
        uint32_t pulse_idx;
        llvm::Function *f;
        std::string name{};
    };
    std::vector<RampFunc> ramp_funcs;
    auto vector_size = get_vector_size();

    for (uint32_t idx = 0; idx < nbseqs; idx++) {
        auto &bseq = bseqs[idx];
        auto &be_bseq = m_seqs[idx];
        auto bc_gen = new BCGen;
        be_bseq.bc_gen.reset(bc_gen);
        bool all_known = true;
        if (m_has_clock) {
            if (auto it = m_clocks.find(idx); it != m_clocks.end()) {
                be_bseq.clock_set = true;
                if (unlikely(mgr().dump_enabled())) {
                    // copy
                    bc_gen->clocks = it->second;
                }
                else {
                    bc_gen->clocks = std::move(it->second);
                }
            }
            else {
                all_known = false;
            }
        }
        bc_gen->seq_idx = idx;
        bc_gen->fpga_clock_div = fpga_clock_div;
        bc_gen->seq_delay = m_seq_delay;
        bc_gen->start_ttl_chn = m_start_ttl_chn;
        bc_gen->ttl_mask = m_ttl_mask;
        bc_gen->first_bseq = idx == 0;
        for (auto [chn_id, ttl_mgr]: m_ttl_managers) {
            auto it = m_chn_map.find(chn_id);
            if (it == m_chn_map.end())
                throw std::runtime_error(name() + ": Unknown channel " +
                                         std::to_string(chn_id) + " for TTL manager");
            if (it->second.first != BCGen::ChnType::TTL)
                throw std::runtime_error(name() + ": TTL manager can only be applied "
                                         "on TTL channels");
            auto chn = it->second.second;
            bc_gen->add_ttl_manager(chn, ttl_mgr.off_delay, ttl_mgr.on_delay,
                                    ttl_mgr.skip_time, ttl_mgr.min_time, ttl_mgr.off_val);
        }
        auto &chn_infos = bseq->get_channel_infos();
        for (auto [chn_id, fpga_chn]: m_chn_map) {
            auto it = chn_infos.find(chn_id);
            assert(it != chn_infos.end());
            auto &chn_info = it->second;
            for (auto &pulse: chn_info.pulses) {
                if (pulse.is_measure())
                    continue;
                auto len = pulse.len();
                if (len && !len->is_const())
                    all_known = false;
                uint32_t len_id = len ? compiler.var_slot(len) : uint32_t(-1);
                auto cond = pulse.cond();
                if (cond) {
                    assert(!cond->is_const());
                    all_known = false;
                }
                uint32_t cond_id = cond ? compiler.var_slot(cond) : uint32_t(-1);
                uint32_t time_id = compiler.time_slot(&pulse.start());
                if (time_id >= expseq.host_seq.nconsts)
                    all_known = false;
                uint32_t endvalue_id = compiler.var_slot(pulse.endval());
                if (endvalue_id >= expseq.host_seq.nconsts)
                    all_known = false;

                auto pulse_type = BCGen::PulseType::Value;
                if (len) {
                    pulse_type = BCGen::PulseType::Scalar;
                    auto val = pulse.val();
                    uint32_t pulse_idx = bc_gen->seq_pulses.size();
                    assert(val->is_call());

                    auto f0 = val->get_callee().llvm;
                    auto it_f = func_map.find(f0);
                    assert(it_f != func_map.end());
                    auto f = it_f->second;
                    LLVM::Codegen::Wrapper wrap{true};
                    LLVM::Codegen::Wrapper wrap_vector{true};
                    auto args = val->args();
                    uint32_t nargs = args.size();
                    for (uint32_t argi = 0; argi < nargs; argi++) {
                        auto arg = args[argi];
                        if (arg.is_arg()) {
                            assert(arg.get_arg() == 0);
                            if (vector_size > 1) {
                                wrap_vector.add_vector(argi);
                                wrap_vector.add_byref(argi);
                            }
                            continue;
                        }
                        assert(arg.is_var());
                        auto slot = compiler.var_slot(arg.get_var());
                        assert(slot >= expseq.host_seq.nconsts);
                        all_known = false;
                        assert(arg.get_var()->llvm_type() == f->getArg(argi)->getType());
                        if (vector_size > 1)
                            wrap_vector.add_closure(argi, slot);
                        wrap.add_closure(argi, slot);
                    }
                    llvm::Function *wrapf = nullptr;
                    if (vector_size > 1) {
                        wrap_vector.vector_size = vector_size;
                        wrap_vector.add_ret_ref();
                        wrapf = cgctx.emit_wrapper(f, "0", wrap_vector);
                        if (wrapf) {
                            pulse_type = BCGen::PulseType::Vector;
                        }
                    }
                    if (!wrapf) {
                        wrapf = cgctx.emit_wrapper(f, "0", wrap);
                        if (!wrapf) {
                            throw std::runtime_error(name() + ": Failed to compile function.");
                        }
                    }
                    ramp_funcs.push_back(RampFunc{idx, pulse_idx, wrapf});
                }
                bc_gen->seq_pulses.push_back({
                        .chn_type = fpga_chn.first,
                        .pulse_type = pulse_type,
                        .chn = fpga_chn.second,

                        .id = pulse.id(),
                        .time_id = time_id,
                        .len_id = len_id,
                        .endvalue_id = endvalue_id,
                        .cond_id = cond_id,
                        .ramp_func = nullptr
                    });
            }
        }
        be_bseq.all_known = all_known;
    }
    // Compile all the code so that we can do another pass that determines
    // if we can pre-generate the sequence.

    // Clean up function exports
    for (auto &go: mod.global_objects()) {
        if (!go.isDeclaration()) {
            go.setLinkage(llvm::GlobalValue::PrivateLinkage);
        }
    }
    for (auto &ramp_func: ramp_funcs) {
        ramp_func.f->setLinkage(llvm::GlobalValue::ExternalLinkage);
        ramp_func.f->setVisibility(llvm::GlobalValue::ProtectedVisibility);
        ramp_func.f->setUnnamedAddr(llvm::GlobalValue::UnnamedAddr::Global);
    }
    llvm::legacy::PassManager PM;
#ifndef NDEBUG
    PM.add(llvm::createVerifierPass());
#endif
    PM.add(llvm::createGlobalDCEPass());
    // Shorten all global names since we don't care what they are
    // and this should slightly reduce the compiled binary size.
    PM.add(LLVM::createGlobalRenamePass());
#ifndef NDEBUG
    PM.add(llvm::createVerifierPass());
#endif
    PM.run(mod);
    // The optimization passes below (in `emit_objfile`) may recreate the functions
    // so the function handle may not be valid anymore.
    // We are done with the function renaming so we can simply get the names now.
    for (auto &ramp_func: ramp_funcs) {
        ramp_func.name = ramp_func.f->getName();
        assert(!ramp_func.name.empty());
        ramp_func.f = nullptr;
    }

    llvm::SmallVector<char,0> vec;
    auto res = LLVM::Compile::emit_objfile(vec, LLVM::Compile::get_native_target(), &mod, true);
    if (!res)
        throw std::runtime_error(name() + ": Sequence function compilation failed.");
    auto &engine = mgr().exe_engine();
    m_obj_id = engine.load(&vec[0], vec.size(),
                           expseq.cgctx()->get_extern_resolver());
    if (!m_obj_id)
        throw std::runtime_error(name() + ": Loading sequence functions failed.");
    for (auto &ramp_func: ramp_funcs) {
        auto sym = (void(*)(void))engine.get_symbol(ramp_func.name);
        assert(sym);
        m_seqs[ramp_func.bseq_idx].bc_gen->seq_pulses[ramp_func.pulse_idx].ramp_func = sym;
    }
    engine.reset_dyld();
    bool all_pregenerated = true;
    for (uint32_t idx = 0; idx < nbseqs; idx++)
        all_pregenerated &= pregenerate(expseq, compiler, idx);
    if (all_pregenerated) {
        // If everything is pregenerated, free the object file.
        mgr().exe_engine().free(m_obj_id);
        m_obj_id = 0;
    }
    m_generated = true;
    if (likely(!mgr().dump_enabled())) {
        m_clocks.clear();
    }
}

inline bool Backend::pregenerate(Manager::ExpSeq &expseq, Compiler &compiler, uint32_t seq_idx)
{
    auto &be_bseq = m_seqs[seq_idx];
    if (!be_bseq.all_known) {
        mgr().add_debug_printf("%s: unable to pregenerate basic sequence [index %d]\n",
                               name().c_str(), seq_idx);
        return false;
    }
    mgr().add_debug_printf("%s: pregeneratig basic sequence [index %d]\n",
                           name().c_str(), seq_idx);
    const auto &bseq = compiler.basic_seqs()[seq_idx];
    auto &host_seq = expseq.host_seq;
    int64_t len = 0;
    for (auto &time: bseq->get_endtimes()) {
        uint32_t time_id = compiler.time_slot(time.get());
        if (time_id >= host_seq.nconsts)
            return false;
        len = max(len, host_seq.get_time(time_id));
    }

    auto bc_gen = be_bseq.bc_gen.get();
    std::sort(bc_gen->seq_pulses.begin(), bc_gen->seq_pulses.end(), [&] (auto &p1, auto &p2) {
            auto t1 = host_seq.get_time(p1.time_id);
            auto t2 = host_seq.get_time(p2.time_id);
            if (t1 < t2)
                return true;
            if (t1 > t2)
                return false;
            return p1.id < p2.id;
        });
    bc_gen->len_ns = uint64_t(len * 1e9 / m_seq_tick_per_sec);

    if (seq_idx == 0) {
        // We can simply use default values for the first sequence as the starting point.
        bc_gen->start_vals.clear();
        auto defvals = compiler.seq()->get_defvals();
        for (auto [chn_id, fpga_chn]: m_chn_map) {
            auto it = defvals.find(chn_id);
            auto value = it == defvals.end() ? 0 : it->second;
            for (auto &pulse: bc_gen->seq_pulses) {
                auto t = host_seq.get_time(pulse.time_id);
                if (t > 0)
                    break;
                if (pulse.chn_type != fpga_chn.first || pulse.chn != fpga_chn.second)
                    continue;
                if (!pulse.ramp_func) {
                    value = host_seq.get_value(pulse.endvalue_id, seq_idx);
                    continue;
                }
                auto data = host_seq.values.data();
                if (pulse.pulse_type == BCGen::PulseType::Scalar) {
                    auto func = (double (*)(double, const void*))pulse.ramp_func;
                    value = func(0, data);
                }
                else {
                    // This needs to be consistent with the bc_gen data stream.
                    // Note that for simplicity
                    // we simply use the largest size on all platforms.
                    // TODO: SVE
                    double buffer[8];
                    double time_buffer[8] = {0};
                    auto func = (void (*)(double*, const double*, const void*))pulse.ramp_func;
                    func(buffer, time_buffer, data);
                    value = buffer[0];
                }
            }
            bc_gen->start_vals.emplace(fpga_chn, BCGen::convert_value(fpga_chn.first, value));
        }
        bc_gen->first_bseq = true;
        bc_gen->generate(host_seq);
        be_bseq.bytecode_first = std::move(bc_gen->bytecode);
        if (!bseq->has_branchin()) {
            be_bseq.bc_gen.reset(nullptr);
            return true;
        }
    }
    bc_gen->start_vals.clear();
    auto &chn_infos = bseq->get_channel_infos();
    for (auto [chn_id, fpga_chn]: m_chn_map) {
        // We only care about TTL for non-first-bseq.
        if (fpga_chn.first != BCGen::ChnType::TTL)
            continue;
        auto it = chn_infos.find(chn_id);
        assert(it != chn_infos.end());
        auto &chn_info = it->second;
        bool hasval = false;
        bool value = false;
        if (chn_info.startval->is_const()) {
            hasval = true;
            value = chn_info.startval->get_const().get<bool>();
        }
        for (auto &pulse: bc_gen->seq_pulses) {
            auto t = host_seq.get_time(pulse.time_id);
            if (t > 0)
                break;
            if (pulse.chn_type != BCGen::ChnType::TTL || pulse.chn != fpga_chn.second)
                continue;
            assert(!pulse.ramp_func);
            hasval = true;
            value = host_seq.get_value(pulse.endvalue_id, seq_idx) != 0;
        }
        if (!hasval)
            return false;
        bc_gen->start_vals.emplace(fpga_chn, uint32_t(value));
    }
    bc_gen->first_bseq = false;
    bc_gen->generate(host_seq);
    be_bseq.bytecode = std::move(bc_gen->bytecode);
    be_bseq.bc_gen.reset(nullptr);
    return true;
}

void Backend::pre_run(HostSeq &host_seq)
{
    auto idx = host_seq.cur_seq_idx();
    auto &bseq = m_seqs[idx];
    const auto &bc = host_seq.first_bseq ? bseq.bytecode_first : bseq.bytecode;
    if (!bc.empty()) {
        mgr().add_debug_printf("%s: basic sequence [index %d] pregenerated\n",
                               name().c_str(), idx);
        return;
    }
    mgr().add_debug_printf("%s: runtime generation for basic sequence [index %d]\n",
                           name().c_str(), idx);
    auto bc_gen = bseq.bc_gen.get();
    assert(bc_gen);
    auto &host_bseq = host_seq.seqs[idx];
    auto chn_value_offset = host_seq.nconsts + host_seq.nglobals + host_seq.nglobal_vals;
    bc_gen->start_vals.clear();
    for (auto [chn_id, fpga_chn]: m_chn_map) {
        if (!host_seq.first_bseq && fpga_chn.first != BCGen::ChnType::TTL)
            continue;
        auto value = host_seq.get_value(chn_id + chn_value_offset - 1, idx);
        bc_gen->start_vals.emplace(fpga_chn, BCGen::convert_value(fpga_chn.first, value));
    }
    bc_gen->first_bseq = host_seq.first_bseq;
    bc_gen->len_ns = uint64_t(host_bseq.length * 1e9 / m_seq_tick_per_sec);
    bc_gen->generate(host_seq);
}

void Backend::start(HostSeq &host_seq)
{
    auto idx = host_seq.cur_seq_idx();
    auto &bseq = m_seqs[idx];
    const auto &bc = host_seq.first_bseq ? bseq.bytecode_first : bseq.bytecode;
    mgr().add_debug_printf("%s: starting basic sequence [index %d], bytecode length: %zu\n",
                           name().c_str(), idx, bc.size());
    run_bytecode(bc.empty() ? bseq.bc_gen->bytecode : bc);
}

void Backend::cancel(HostSeq &host_seq)
{
    if (!m_sock) {
        if (mgr().dummy_mode)
            return;
        throw std::runtime_error(name() + ": Backend not configured.");
    }
    mgr().add_debug_printf("%s: cancelling\n", name().c_str());
    m_sock->send_msg([&] (auto &sock) {
        ZMQ::send_more(sock, ZMQ::str_msg("cancel_seq"));
        ZMQ::send(sock, zmq::message_t(m_seq_id, sizeof(uint64_t) * 2));
    }).wait();
}

void Backend::wait(HostSeq &host_seq)
{
    if (!m_sock) {
        if (mgr().dummy_mode)
            return;
        throw std::runtime_error(name() + ": Backend not configured.");
    }
    auto reply = m_sock->send_msg([&] (auto &sock) {
        ZMQ::send_more(sock, ZMQ::str_msg("wait_seq"));
        zmq::message_t msg(sizeof(uint64_t) * 2 + 1);
        memcpy(msg.data(), m_seq_id, sizeof(uint64_t) * 2);
        ((char*)msg.data())[sizeof(uint64_t) * 2] = 2; // Wait for finish
        ZMQ::send(sock, std::move(msg));
    });
    reply.wait();
    mgr().add_debug_printf("%s: finished\n", name().c_str());
    auto msgs = reply.get();
    if (msgs.size() < 1) {
        mgr().add_warning(name() + " wait_seq: return empty result\n");
        return;
    }
    auto &msg = msgs[0];
    if (msg.size() < 1) {
        mgr().add_warning(name() + " wait_seq: return empty status code\n");
        return;
    }
    int status = *(uint8_t*)msg.data();
    if (status == 0)
        return;
    if (status == 1) {
        mgr().add_error(name() + ": sequence cancelled");
    }
    else {
        mgr().add_warning(name() + ": unknown return status: " + std::to_string(status));
    }
}

void Backend::config(const YAML::Node &config)
{
    // Currently this is a global setting.
    // If needed, we could make this a number that different backends can negotiate
    // so that we don't need to waste extra time when not needed.
    if (auto seq_delay_ms_node = config["seq_delay_ms"]) {
        auto seq_delay_ms = seq_delay_ms_node.as<double>();
        if (seq_delay_ms <= 0)
            seq_delay_ms = 0;
        if (seq_delay_ms > 3600 * 1000)
            throw std::runtime_error(name() + ": Sequence delay (" +
                                     std::to_string(seq_delay_ms) +
                                     "ms) too long. (Max 1 hour)");
        m_seq_delay = uint32_t(seq_delay_ms * fpga_tick_per_sec / 1000);
    }
    else {
        m_seq_delay = 0;
    }
    if (auto start_ttl_chn_node = config["start_ttl_chn"]) {
        m_start_ttl_chn = start_ttl_chn_node.as<int>();
        if (m_start_ttl_chn < 0 || m_start_ttl_chn >= 32) {
            throw std::runtime_error(name() + ": Invalid start_ttl_chn " +
                                     std::to_string(m_start_ttl_chn));
        }
    }
    else {
        throw std::runtime_error(name() + ": Missing start_ttl_chn in config");
    }
    if (mgr().dummy_mode) {
        m_sock.reset(nullptr);
    }
    else if (auto url_node = config["url"]) {
        auto new_url = url_node.as<std::string>();
        if (new_url != m_url) {
            m_url = std::move(new_url);
            auto &client = ZMQ::MultiClient::global();
            m_sock.reset(new ZMQ::MultiClient::SockRef(client.get_socket(m_url)));
        }
    }
    else {
        throw std::runtime_error(name() + ": Missing url in config");
    }
    m_ttl_ovr_ignore = 0;
    m_dds_ovr_ignore.clear();
    if (auto override_ignore_node = config["override_ignore"]) {
        if (!override_ignore_node.IsSequence())
            throw std::runtime_error(name() + ": override_ignore must be an array");
        for (auto chn_node: override_ignore_node) {
            auto chn_string = chn_node.as<std::string>();
            llvm::StringRef chn_name(chn_string);
            if (chn_name.startswith("TTL")) {
                auto chn_num_str = chn_name.substr(3);
                uint8_t chn_num;
                if (chn_num_str.getAsInteger(10, chn_num) || chn_num >= 32)
                    throw std::runtime_error(name() + ": Invalid TTL channel number " +
                                             chn_string + " in override_ignore");
                if (chn_num == m_start_ttl_chn)
                    throw std::runtime_error(name() + ": Cannot ignore override "
                                             "on start trigger TTL channel");
                m_ttl_ovr_ignore = setBit(m_ttl_ovr_ignore, chn_num, true);
            }
            else if (chn_name.startswith("DDS")) {
                auto subname = chn_name.substr(3);
                uint8_t chn_num;
                if (subname.consumeInteger(10, chn_num) || chn_num >= 22)
                    throw std::runtime_error(name() + ": Invalid DDS channel " +
                                             chn_string + " in override_ignore");
                if (subname == "/FREQ") {
                    m_dds_ovr_ignore.insert(chn_num);
                }
                else if (subname == "/AMP") {
                    m_dds_ovr_ignore.insert(uint8_t((1 << 6) | chn_num));
                }
                else if (subname == "/PHASE") {
                    m_dds_ovr_ignore.insert(uint8_t((2 << 6) | chn_num));
                }
                else if (subname == "") {
                    m_dds_ovr_ignore.insert(chn_num);
                    m_dds_ovr_ignore.insert(uint8_t((1 << 6) | chn_num));
                    m_dds_ovr_ignore.insert(uint8_t((2 << 6) | chn_num));
                }
                else {
                    throw std::runtime_error(name() + ": Invalid DDS parameter name " +
                                             subname.str());
                }
            }
            else {
                throw std::runtime_error(name() + ": Unknown channel name " + chn_string +
                                         " in override_ignore.");
            }
        }
    }
}

inline void Backend::run_bytecode(const std::vector<uint8_t> &bc)
{
    if (!m_sock) {
        if (mgr().dummy_mode)
            return;
        throw std::runtime_error(name() + ": Backend not configured.");
    }
    auto reply = m_sock->send_msg([&] (auto &sock) {
        ZMQ::send_more(sock, ZMQ::str_msg("run_seq"));
        ZMQ::send_more(sock, ZMQ::bits_msg(BCGen::version()));
        ZMQ::send(sock, zmq::message_t(bc.data(), bc.size()));
    }).get();
    if (reply.empty() || reply[0].size() < 16)
        throw std::runtime_error(name() + ": Failed to start FPGA sequence.");
    auto &reply0 = reply[0];
    auto rep_data = (const uint8_t*)reply[0].data();
    memcpy(m_seq_id, rep_data, sizeof(m_seq_id));
    if (reply0.size() < 24)
        return;
    std::string msg;
    auto get_msg = [&] () -> std::string& {
        if (msg.empty()) {
            msg = name() + " channel overridden:";
        }
        else {
            msg += ", ";
        }
        return msg;
    };
    uint32_t ttl_ovr = (Mem::load_unalign<uint32_t>(rep_data + 16) |
                        Mem::load_unalign<uint32_t>(rep_data + 20)) & ~m_ttl_ovr_ignore;
    auto dds_ovr_start = rep_data + 24;
    auto dds_ovr_end = rep_data + reply0.size();
    if (ttl_ovr) {
        for (int i = 0; i < 32; i++) {
            if (ttl_ovr & (uint32_t(1) << i)) {
                get_msg() += "TTL" + std::to_string(i);
            }
        }
    }
    // The format is 1 byte of channel identifier followed by 4 bytes of overriden values.
    for (auto p = dds_ovr_start; p < dds_ovr_end; p += 5) {
        auto chn = *p;
        if (m_dds_ovr_ignore.find(chn) != m_dds_ovr_ignore.end())
            continue;
        auto type = chn >> 6;
        int chn_num = chn & 0x3f;
        if (type == 0) {
            get_msg() += "FREQ" + std::to_string(chn_num);
        }
        else if (type == 1) {
            get_msg() += "AMP" + std::to_string(chn_num);
        }
        else if (type == 2) {
            get_msg() += "PHASE" + std::to_string(chn_num);
        }
        else {
            get_msg() += "unknown(" + std::to_string(int(chn)) + ")";
        }
    }
    if (!msg.empty()) {
        mgr().add_warning(msg + "\n");
    }
}

void Backend::parse_data(const uint8_t *data, size_t len)
{
    mgr().add_debug_printf("%s: parsing backend data\n", name().c_str());
    // Format:
    //   [magic <"ZYNQZYNQ">: 8B]
    //   [version <0>: 1B]
    //   [nttl_mgrs: 1B][[chn_id: 4B][off_delay: 4B][on_delay: 4B]
    //                   [skip_time: 4B][min_time: 4B][off_val: 1B] x nttl_mgrs]
    Mem::Reader reader(data, len);
    // [magic <"ZYNQZYNQ">: 8B]
    // The magic byte is here to make sure that the backend data
    // has reached the desired backend.
    // The caller otherwise only know the device name and may not know if the backend
    // has been configured correctly or not.
    if (memcmp(reader.read_array<char>(8), "ZYNQZYNQ", 8) != 0)
        throw std::runtime_error(name() + ": Incorrect backend data target.");
    // [version <0>: 1B]
    int ver = reader.read<int8_t>();
    if (ver != 0)
        throw std::runtime_error(name() + ": Unknown FPGA backend data version " +
                                 std::to_string(ver));
    // [nttl_mgrs: 1B][[chn_id: 4B][off_delay: 4B][on_delay: 4B]
    //                 [skip_time: 4B][min_time: 4B][off_val: 1B] x nttl_mgrs]
    auto nttl_mgrs = reader.read<uint8_t>();
    for (auto i = 0; i < nttl_mgrs; i++) {
        auto chn_id = reader.read<uint32_t>();
        auto off_delay = reader.read<uint32_t>();
        auto on_delay = reader.read<uint32_t>();
        auto skip_time = reader.read<uint32_t>();
        auto min_time = reader.read<uint32_t>();
        auto off_val = reader.read<uint8_t>() != 0;
        m_ttl_managers.emplace(chn_id, TTLManager{off_delay, on_delay,
                skip_time, min_time, off_val});
    }
}

NACS_EXPORT() Backend *Backend::cast(Device *dev)
{
    return dynamic_cast<Backend*>(dev);
}

NACS_EXPORT() const Backend *Backend::cast(const Device *dev)
{
    return dynamic_cast<const Backend*>(dev);
}

NACS_EXPORT() bool Backend::has_generator(uint32_t bseq_idx) const
{
    return bool(m_seqs[bseq_idx].bc_gen);
}

NACS_EXPORT() bool Backend::pregenerated(uint32_t bseq_idx, bool first_bseq) const
{
    auto &bseq = m_seqs[bseq_idx];
    const auto &bc = first_bseq ? bseq.bytecode_first : bseq.bytecode;
    return !bc.empty();
}

NACS_EXPORT() llvm::ArrayRef<uint8_t> Backend::get_bytecode(uint32_t bseq_idx,
                                                            bool first_bseq) const
{
    auto &bseq = m_seqs[bseq_idx];
    const auto &bc = first_bseq ? bseq.bytecode_first : bseq.bytecode;
    return bc.empty() ? bseq.bc_gen->bytecode : bc;
}

NACS_EXPORT() llvm::ArrayRef<BCGen::Clock> Backend::get_clock(uint32_t bseq_idx) const
{
    auto it = m_clocks.find(bseq_idx);
    if (it == m_clocks.end())
        return llvm::ArrayRef<BCGen::Clock>();
    return it->second;
}

}
