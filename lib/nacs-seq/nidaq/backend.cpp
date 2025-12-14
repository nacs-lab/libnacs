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

#include "../zynq/backend.h"

#include "data_gen.h"

#include "../../nacs-utils/llvm/utils.h"
#include "../../nacs-utils/llvm/codegen.h"
#include "../../nacs-utils/llvm/compile.h"
#include "../../nacs-utils/llvm/execute.h"
#include "../../nacs-utils/llvm/passes.h"
#include "../../nacs-utils/log.h"
#include "../../nacs-utils/processor.h"
#include "../../nacs-utils/streams.h"

#include <llvm/ADT/StringRef.h>

#include <memory>

namespace NaCs::Seq::NiDAQ {

static Device::Register<Backend> register_backend("nidaq");

struct Backend::ChannelInfo {
    std::string dev_name;
    uint32_t chn_num;
    uint32_t linear_idx = 0; // Linear channel index used in data generator

    ChannelInfo(std::string &&dev_name, uint32_t chn_num)
        : dev_name(std::move(dev_name)),
          chn_num(chn_num)
    {}
};

struct Backend::BasicSeq {
    struct DataGenExt : DataGen {
        // Data generator with some extra book keeping information
        // that are only useful during generation...
        // Within `prepare` this is at first a map from basic sequence specific
        // sequence value index to the basic sequence specific backend value index.
        // Afterwards and in all other functions,
        // this is the mapping from non-constant backend index
        // to sequence value index.
        std::vector<uint32_t> val_map;
    };
    // Note that empty `data` may be valid as well.
    // Whether we need to compute the data at runtime
    // is signaled by the existance of `data_gen`.
    std::unique_ptr<DataGenExt> data_gen;
    std::vector<std::pair<uint32_t,size_t>> fill_init;
    std::vector<double> data;
    size_t nsamples;
    // In general, we don't have an initialization sequence so the first basic sequence
    // isn't very special. However, if the a basic sequence contains no pulse,
    // we can usually simply not output anything except when it's the first basic sequence
    // since we still need to apply the default values for all the used channels.
    // If `first_only` is `true`, the data should only be used for the first basic sequence
    // and not for anything else.
    bool first_only = false;
    bool all_time_known;
    bool all_val_known;
    uint32_t nspecific_map;
};

Backend::Backend(Manager &mgr, std::string name)
    : Device(mgr, std::move(name))
{
}

Backend::~Backend()
{
    if (m_obj_id) {
        mgr().exe_engine().free(m_obj_id);
    }
}

NACS_EXPORT() const double *Backend::get_data(uint32_t cur_seq_id, size_t *sz) const
{
    auto &bseq = m_seqs[cur_seq_id];
    if (bseq.data_gen) {
        *sz = bseq.data_gen->data.size();
        return bseq.data_gen->data.data();
    }
    if ((!m_first_bseq && bseq.first_only) || bseq.nsamples == 0) {
        *sz = 0;
        return nullptr;
    }
    *sz = bseq.data.size();
    return bseq.data.data();
}

void Backend::add_channel(uint32_t chn_id, const std::string &_chn_name)
{
    llvm::StringRef chn_name(_chn_name);
    auto [dev_name, chn_num_str] = chn_name.split("/");
    if (chn_num_str.empty())
        throw std::runtime_error(name() + ": No device channel ID in Ni DAQ channel name.");
    if (dev_name.empty())
        throw std::runtime_error(name() + ": Empty device name in Ni DAQ channel name.");
    uint32_t chn_num;
    if (chn_num_str.getAsInteger(10, chn_num))
        throw std::runtime_error(name() + ": Device channel ID must be a number.");
    m_chn_map.try_emplace(chn_id, dev_name.str(), chn_num);
}

static unsigned get_vector_size()
{
    // This needs to be consistent with the data_gen data stream.
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

void Backend::sort_channels()
{
    auto nchn = (uint32_t)m_chn_map.size();
    m_linear_chns.resize(nchn);
    auto it = m_chn_map.begin();
    for (uint32_t i = 0; i < nchn; i++)
        m_linear_chns[i] = it++;
    // Sort the channels based on their name so that the ordering is stable/predictable
    std::sort(m_linear_chns.begin(), m_linear_chns.end(), [&] (auto it1, auto it2) {
        auto &info1 = it1->second;
        auto &info2 = it2->second;
        if (info1.dev_name == info2.dev_name)
            return info1.chn_num < info2.chn_num;
        return info1.dev_name < info2.dev_name;
    });
    for (uint32_t i = 0; i < nchn; i++) {
        m_linear_chns[i]->second.linear_idx = i;
    }
}

NACS_EXPORT() const std::pair<uint32_t,const char*> *Backend::get_channel_info(uint32_t *sz)
{
    if (m_channel_info.empty()) {
        for (auto it: m_linear_chns) {
            m_channel_info.push_back({it->second.chn_num, it->second.dev_name.c_str()});
        }
    }
    *sz = uint32_t(m_channel_info.size());
    return m_channel_info.data();
}

static llvm::Function *clamp_ramp_range(LLVM::Codegen::Context &cgctx, llvm::Function *f)
{
    // TODO: nan check? directly use min/max instruction when available instead of calling sleef
    auto fty = f->getFunctionType();
    auto newf = llvm::Function::Create(fty,  llvm::GlobalValue::ExternalLinkage,
                                       f->getName() + ".c", f->getParent());
    newf->setVisibility(llvm::GlobalValue::ProtectedVisibility);
    newf->setAttributes(llvm::AttributeList::get(
                            f->getContext(), LLVM::getFnAttrs(*f), {}, {}));
    auto *b0 = llvm::BasicBlock::Create(f->getContext(), "top", newf);
    auto nargs = f->arg_size();
    llvm::IRBuilder<> builder(b0);
    llvm::SmallVector<llvm::Value*,8> call_args(nargs);
    for (uint32_t i = 0; i < nargs; i++)
        call_args[i] = newf->getArg(i);
    llvm::Value *res = builder.CreateCall(f, call_args);
    auto max_f = LLVM::get_intrinsic(f->getParent(), llvm::Intrinsic::maxnum,
                                     {cgctx.T_f64});
    auto min_f = LLVM::get_intrinsic(f->getParent(), llvm::Intrinsic::minnum,
                                     {cgctx.T_f64});
    res = builder.CreateCall(min_f, {res, llvm::ConstantFP::get(cgctx.T_f64, 10)});
    res = builder.CreateCall(max_f, {res, llvm::ConstantFP::get(cgctx.T_f64, -10)});
    builder.CreateRet(res);
    return newf;
}

void Backend::prepare(Manager::ExpSeq &expseq, Compiler &compiler)
{
    mgr().add_debug_printf("%s: compiling\n", name().c_str());
    sort_channels();
    m_zynq_dev = dynamic_cast<Zynq::Backend*>(expseq.get_device(m_clock_dev, true));
    if (!m_zynq_dev)
        throw std::runtime_error(name() + ": Cannot find clock generator device.");
    m_zynq_dev->enable_clock();
    const auto bseqs = compiler.basic_seqs();
    auto &host_seq = expseq.host_seq;
    auto nbseqs = bseqs.size();
    m_seqs.resize(nbseqs);
    std::vector<uint32_t> shared_val_map(host_seq.nshared, uint32_t(-1));

    uint32_t nconst_map = 0;
    uint32_t nshared_map = 0;
    uint32_t nspecific_map_max = 0;
    // First pass: map the value IDs from the global sequence one to the one
    // specific to our backend.
    // We map all the basic sequence specific values after all the shared ones
    // so we need to count all the mapped shared ones before we know the final internal ID
    // of backend specific values.
    for (uint32_t idx = 0; idx < nbseqs; idx++) {
        auto &bseq = bseqs[idx];
        auto &chn_infos = bseq->get_channel_infos();
        auto &host_bseq = host_seq.seqs[idx];
        auto &be_bseq = m_seqs[idx];
        auto data_gen = new BasicSeq::DataGenExt;
        be_bseq.data_gen.reset(data_gen);
        bool all_time_known = true;
        bool all_val_known = true;
        bool has_pulse = false;
        data_gen->val_map.resize(host_bseq.nmeasure + host_bseq.ndirect +
                                 host_bseq.nneed_order, uint32_t(-1));
        uint32_t nspecific_map = 0;
        auto map_slot = [&] (uint32_t slot) {
            if (slot < host_seq.nconsts) {
                if (shared_val_map[slot] == uint32_t(-1))
                    shared_val_map[slot] = nconst_map++;
                return true;
            }
            else if (slot < host_seq.nshared) {
                if (shared_val_map[slot] == uint32_t(-1))
                    shared_val_map[slot] = nshared_map++;
                return false;
            }
            slot -= host_seq.nshared;
            assert(slot < data_gen->val_map.size());
            if (data_gen->val_map[slot] == uint32_t(-1))
                data_gen->val_map[slot] = nspecific_map++;
            return false;
        };
        for (auto &[chn_id, be_chn_info]: m_chn_map) {
            auto it = chn_infos.find(chn_id);
            assert(it != chn_infos.end());
            auto &chn_info = it->second;
            for (auto &pulse: chn_info.pulses) {
                if (pulse.is_measure())
                    continue;
                has_pulse = true;
                auto len = pulse.len();
                if (len)
                    all_time_known &= map_slot(compiler.var_slot(len));
                if (auto cond = pulse.cond()) {
                    assert(!map_slot(compiler.var_slot(cond)));
                    all_time_known &= map_slot(compiler.var_slot(cond));
                }
                all_time_known &= map_slot(compiler.time_slot(&pulse.start()));
                all_val_known &= map_slot(compiler.var_slot(pulse.endval()));

                if (len) {
                    auto val = pulse.val();
                    assert(val->is_call());
                    for (auto arg: val->args()) {
                        if (arg.is_arg())
                            continue;
                        assert(arg.is_var());
                        auto slot = compiler.var_slot(arg.get_var());
                        all_val_known &= map_slot(slot);
                        assert(slot >= expseq.host_seq.nconsts);
                        assert(!all_val_known);
                    }
                }
            }
        }
        be_bseq.all_time_known = all_time_known;
        be_bseq.all_val_known = all_val_known;
        be_bseq.nspecific_map = nspecific_map;
        if (!has_pulse) {
            be_bseq.data_gen.reset(nullptr);
            assert(nspecific_map == 0);
            continue;
        }
        data_gen->nchns = (uint32_t)m_linear_chns.size();
        data_gen->step_size = m_step_size;
        data_gen->start_values.resize(data_gen->nchns);
        nspecific_map_max = max(nspecific_map_max, nspecific_map);
    }

    m_nconst_vals = nconst_map;
    // Second pass: map all the values and create LLVM functions.
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
        for (auto &[chn_id, be_chn_info]: m_chn_map) {
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
                    f = clamp_ramp_range(cgctx, f);
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
    LLVM::runAlwaysInlinerPasses(mod);

    auto nvalues = nconst_map + nshared_map + nspecific_map_max;
    auto values_ty = llvm::ArrayType::get(cgctx.T_i8, nvalues * sizeof(double));
    auto values_gv = new llvm::GlobalVariable(mod, values_ty, false,
                                              llvm::GlobalValue::ExternalLinkage,
                                              llvm::ConstantAggregateZero::get(values_ty),
                                              "values");
    struct RampInfo {
        uint32_t bseq_idx;
        uint32_t linear_chn;
        uint32_t pulse_idx;
        llvm::Function *ramp_func;
        std::string name{};
    };
    std::vector<RampInfo> ramp_funcs;
    auto vector_size = get_vector_size();

    for (uint32_t idx = 0; idx < nbseqs; idx++) {
        auto &bseq = bseqs[idx];
        auto &chn_infos = bseq->get_channel_infos();
        auto &be_bseq = m_seqs[idx];
        auto data_gen = be_bseq.data_gen.get();
        // Empty
        if (!data_gen)
            continue;
        auto map_slot = [&] (uint32_t slot) {
            if (slot < host_seq.nconsts) {
                assert(shared_val_map[slot] != uint32_t(-1));
                return shared_val_map[slot];
            }
            else if (slot < host_seq.nshared) {
                assert(shared_val_map[slot] != uint32_t(-1));
                return shared_val_map[slot]+ nconst_map;
            }
            slot -= host_seq.nshared;
            assert(slot < data_gen->val_map.size());
            return data_gen->val_map[slot] + nconst_map + nshared_map;
        };
        data_gen->types.resize(nconst_map + nshared_map + be_bseq.nspecific_map);
        auto nconsts = host_seq.nconsts;
        for (uint32_t slot = 0; slot < nconsts; slot++) {
            if (shared_val_map[slot] == uint32_t(-1))
                continue;
            data_gen->types[shared_val_map[slot]] = host_seq.get_type(slot, idx);
        }
        auto nshared = host_seq.nshared;
        for (uint32_t slot = nconsts; slot < nshared; slot++) {
            if (shared_val_map[slot] == uint32_t(-1))
                continue;
            data_gen->types[shared_val_map[slot] + nconst_map] = host_seq.get_type(slot, idx);
        }
        auto nspecific = data_gen->val_map.size();
        for (uint32_t slot = 0; slot < nspecific; slot++) {
            if (data_gen->val_map[slot] == uint32_t(-1))
                continue;
            data_gen->types[data_gen->val_map[slot] + nconst_map + nshared_map] =
                host_seq.get_type(slot + nshared, idx);
        }

        for (auto &[chn_id, be_chn_info]: m_chn_map) {
            auto it = chn_infos.find(chn_id);
            assert(it != chn_infos.end());
            auto &chn_info = it->second;
            auto &pulses = data_gen->get_pulses(be_chn_info.linear_idx);
            for (auto &pulse: chn_info.pulses) {
                if (pulse.is_measure())
                    continue;
                auto len = pulse.len();
                uint32_t len_id = len ? map_slot(compiler.var_slot(len)) : uint32_t(-1);
                auto cond = pulse.cond();
                uint32_t cond_id = cond ? map_slot(compiler.var_slot(cond)) : uint32_t(-1);
                uint32_t time_id = map_slot(compiler.time_slot(&pulse.start()));
                uint32_t endvalue_id = map_slot(compiler.var_slot(pulse.endval()));

                auto pulse_type = DataGen::PulseType::Value;
                if (len) {
                    pulse_type = DataGen::PulseType::Scalar;
                    auto val = pulse.val();
                    auto pulse_idx = (uint32_t)pulses.size();
                    assert(val->is_call());

                    auto f0 = val->get_callee().llvm;
                    auto it_f = func_map.find(f0);
                    assert(it_f != func_map.end());
                    auto f = it_f->second;
                    LLVM::Codegen::Wrapper wrap;
                    LLVM::Codegen::Wrapper wrap_vector;
                    wrap.closure_ptr = values_gv;
                    wrap_vector.closure_ptr = values_gv;
                    bool has_closure = false;
                    auto args = val->args();
                    uint32_t nargs = args.size();
                    for (uint32_t argi = 0; argi < nargs; argi++) {
                        auto arg = args[argi];
                        if (arg.is_arg()) {
                            assert(arg.get_arg() == 0);
                            if (vector_size > 1) {
                                wrap_vector.add_vector(argi);
                                // Argument pointer is guaranteed to be aligned properly.
                                wrap_vector.add_byref(argi);
                            }
                            continue;
                        }
                        assert(arg.is_var());
                        auto slot = map_slot(compiler.var_slot(arg.get_var()));
                        assert(arg.get_var()->llvm_type() == f->getArg(argi)->getType());
                        if (vector_size > 1)
                            wrap_vector.add_closure(argi, slot);
                        wrap.add_closure(argi, slot);
                        has_closure = true;
                    }
                    wrap.closure = has_closure;
                    wrap_vector.closure = has_closure;
                    llvm::Function *wrapf = nullptr;
                    if (vector_size > 1) {
                        wrap_vector.vector_size = vector_size;
                        // The return value may not be aligned
                        wrap_vector.add_ret_ref(sizeof(double));
                        wrapf = cgctx.emit_wrapper(f, "0", wrap_vector);
                        if (wrapf) {
                            pulse_type = DataGen::PulseType::Vector;
                        }
                    }
                    if (!wrapf) {
                        wrapf = cgctx.emit_wrapper(f, "0", wrap);
                        if (!wrapf) {
                            throw std::runtime_error(name() + ": Failed to compile function.");
                        }
                    }
                    ramp_funcs.push_back({
                            .bseq_idx = idx, .linear_chn = be_chn_info.linear_idx,
                            .pulse_idx = pulse_idx, .ramp_func = wrapf});
                }
                pulses.emplace_back(pulse_type, pulse.id(), time_id,
                                    len_id, endvalue_id, cond_id, nullptr);
            }
        }
        // Convert the value map to backend idx -> sequence idx
        // This is what we care about later
        auto val_map = std::move(data_gen->val_map);
        data_gen->val_map.resize(nshared_map + be_bseq.nspecific_map);
        for (uint32_t slot = nconsts; slot < nshared; slot++) {
            if (shared_val_map[slot] == uint32_t(-1))
                continue;
            data_gen->val_map[shared_val_map[slot]] = slot;
        }
        for (uint32_t slot = 0; slot < nspecific; slot++) {
            if (val_map[slot] == uint32_t(-1))
                continue;
            data_gen->val_map[val_map[slot] + nshared_map] = slot + nshared;
        }
    }

    // Compile all the code so that we can do another pass that determines
    // if we can pre-generate the sequence.

    // Clean up function exports
    for (auto &go: mod.global_objects()) {
        if (!go.isDeclaration()) {
            go.setLinkage(llvm::GlobalValue::PrivateLinkage);
        }
    }
    for (auto &ramp: ramp_funcs) {
        ramp.ramp_func->setLinkage(llvm::GlobalValue::ExternalLinkage);
        ramp.ramp_func->setVisibility(llvm::GlobalValue::ProtectedVisibility);
        ramp.ramp_func->setUnnamedAddr(llvm::GlobalValue::UnnamedAddr::Global);
    }
    values_gv->setLinkage(llvm::GlobalValue::ExternalLinkage);
    values_gv->setVisibility(llvm::GlobalValue::ProtectedVisibility);

    LLVM::runGlobalRenamePasses(mod);

    // The optimization passes below (in `emit_objfile`) may recreate the functions
    // so the function handle may not be valid anymore.
    // We are done with the function renaming so we can simply get the names now.
    for (auto &ramp: ramp_funcs) {
        ramp.name = ramp.ramp_func->getName();
        assert(!ramp.name.empty());
        ramp.ramp_func = nullptr;
    }
    auto values_gv_name = values_gv->getName().str();

    llvm::SmallVector<char,0> vec;
    auto res = LLVM::Compile::emit_objfile(vec, LLVM::Compile::get_native_target(), &mod, true);
    if (!res)
        throw std::runtime_error(name() + ": Sequence function compilation failed.");
    auto &engine = mgr().exe_engine();
    m_obj_id = engine.load(&vec[0], vec.size(), expseq.cgctx()->get_extern_resolver());
    if (!m_obj_id)
        throw std::runtime_error(name() + ": Loading sequence functions failed.");
    HostSeq::Value *values_ptr = nullptr;
    if (nvalues != 0) {
        // The windows loader returned a null address when the size is 0.
        // We don't really care in this case so simply use NULL.
        values_ptr = (HostSeq::Value*)engine.get_symbol(values_gv_name);
        assert(values_ptr);
    }
    // Fill in the constant values here. This saves us from copying these later.
    // Note that these indices are not included in `val_map`.
    for (uint32_t i = 0; i < host_seq.nconsts; i++) {
        if (shared_val_map[i] == uint32_t(-1))
            continue;
        values_ptr[shared_val_map[i]] = host_seq.values[i];
    }
    for (auto &ramp: ramp_funcs) {
        auto sym = (void(*)(void))engine.get_symbol(ramp.name);
        assert(sym);
        auto &pulses = m_seqs[ramp.bseq_idx].data_gen->get_pulses(ramp.linear_chn);
        pulses[ramp.pulse_idx].ramp_func = sym;
    }
    for (uint32_t idx = 0; idx < nbseqs; idx++) {
        auto &be_bseq = m_seqs[idx];
        if (auto data_gen = be_bseq.data_gen.get()) {
            data_gen->values = values_ptr;
        }
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
}

void Backend::generate(Manager::ExpSeq&, Compiler&)
{
    // All the work has been done in `prepare`, nothing to do here...
}

inline bool Backend::pregenerate(Manager::ExpSeq &expseq, Compiler &compiler, uint32_t seq_idx)
{
    auto &be_bseq = m_seqs[seq_idx];
    auto data_gen = be_bseq.data_gen.get();
    assert(m_zynq_dev);
    // Not needed
    if (!data_gen) {
        if (seq_idx == 0) {
            mgr().add_debug_printf("%s: generate init data for empty first basic sequence\n",
                                   name().c_str());
            // For the first sequence, we always need to output something
            // to initialize the channels
            be_bseq.first_only = true;
            std::pair<int64_t,int64_t> active_time[] = {{0, 0}};
            be_bseq.data.resize(m_linear_chns.size());
            be_bseq.nsamples = 1;
            auto defvals = compiler.seq()->get_defvals();
            for (auto [chn_id, be_chn_info]: m_chn_map) {
                auto it = defvals.find(chn_id);
                auto value = it == defvals.end() ? 0 : it->second;
                be_bseq.data[be_chn_info.linear_idx] = value;
            }
            m_zynq_dev->set_clock_active_time(seq_idx, m_step_size / 2, active_time);
        }
        else {
            mgr().add_debug_printf("%s: no output needed for basic sequence [index %d]\n",
                                   name().c_str(), seq_idx);
            be_bseq.nsamples = 0;
            m_zynq_dev->set_clock_active_time(seq_idx, m_step_size / 2, {});
        }
        return true;
    }
    if (!be_bseq.all_time_known)
        return false;
    // All the constant values are already in place
    // so we just need to let the data generator generate the times.
    data_gen->compute_times();
    m_zynq_dev->set_clock_active_time(seq_idx, m_step_size / 2, data_gen->active_times);
    if (!be_bseq.all_val_known) {
        mgr().add_debug_printf("%s: basic sequence [index %d] has known output times "
                               "but undetermined output values.\n",
                               name().c_str(), seq_idx);
        return false;
    }
    mgr().add_debug_printf("%s: basic sequence [index %d] has known output times and values.\n",
                           name().c_str(), seq_idx);
    const auto &bseq = compiler.basic_seqs()[seq_idx];
    auto &chn_infos = bseq->get_channel_infos();
    // All the pulses are know but we still need to know the start values.
    for (auto &[chn_id, be_chn_info]: m_chn_map) {
        auto it = chn_infos.find(chn_id);
        assert(it != chn_infos.end());
        auto &chn_info = it->second;
        double value = 0;
        if (chn_info.startval->is_const()) {
            value = chn_info.startval->get_const().get<double>();
        }
        else {
            auto &pulses = data_gen->get_pulses(be_chn_info.linear_idx);
            if (!pulses.empty()) {
                auto &pulse = pulses.front();
                if (pulse.start_step > 0) {
                    be_bseq.fill_init.emplace_back(be_chn_info.linear_idx,
                                                   size_t(pulse.start_step));
                }
            }
            else {
                be_bseq.fill_init.emplace_back(be_chn_info.linear_idx, size_t(-1));
            }
        }
        data_gen->start_values[be_chn_info.linear_idx] = value;
    }
    data_gen->generate_data();
    be_bseq.data = std::move(data_gen->data);
    auto nsamples = data_gen->nsamples;
    be_bseq.nsamples = nsamples;
    // Clamp fill_init value to the sample numbers.
    // Note that we need this even for ones that has a start_step
    // since if there's only disabled pulses, the start_step can overflow the buffer.
    for (auto &[linear_chn, nele]: be_bseq.fill_init)
        nele = min(nele, nsamples);
    be_bseq.data_gen.reset(nullptr);
    return true;
}

void Backend::populate_values(HostSeq &host_seq, BasicSeq &bseq)
{
    auto data_gen = bseq.data_gen.get();
    assert(data_gen);
    auto nvals = uint32_t(data_gen->val_map.size());
    for (uint32_t i = 0; i < nvals; i++) {
        data_gen->values[m_nconst_vals + i] = host_seq.values[data_gen->val_map[i]];
    }
}

// We need to generate the active times before the clock backend.
void Backend::prepare_run(HostSeq &host_seq)
{
    auto idx = host_seq.cur_seq_idx();
    auto &bseq = m_seqs[idx];
    auto data_gen = bseq.data_gen.get();
    // Already generated or not needed
    if (bseq.all_time_known)
        return;
    assert(data_gen);
    mgr().add_debug_printf("%s: set output clock for basic sequence [index %d] at runtime.\n",
                           name().c_str(), idx);
    populate_values(host_seq, bseq);
    data_gen->compute_times();
    m_zynq_dev->set_clock_active_time(idx, m_step_size / 2, data_gen->active_times);
}

void Backend::pre_run(HostSeq &host_seq)
{
    auto idx = host_seq.cur_seq_idx();
    m_first_bseq = host_seq.first_bseq;
    auto &bseq = m_seqs[idx];
    auto data_gen = bseq.data_gen.get();
    // Already generated or not needed
    if (!data_gen) {
        for (auto [linear_chn, nele]: bseq.fill_init) {
            auto data = &bseq.data[bseq.nsamples * linear_chn];
            auto value = host_seq.start_values[m_linear_chns[linear_chn]->first - 1].f64;
            value = std::isnan(value) ? 0 : bound(-10, value, 10);
            std::fill(data, data + nele, value);
        }
        return;
    }
    mgr().add_debug_printf("%s: generating data for basic sequence [index %d] at runtime.\n",
                           name().c_str(), idx);
    // When `all_time_known`, we have not called `populate_values` in `prepare_run` yet
    // so we need to do it here. Note that the `populate_values` in `prepare_run`
    // is needed for `compute_times`.
    if (bseq.all_time_known)
        populate_values(host_seq, bseq);
    auto nchns = (uint32_t)m_linear_chns.size();
    for (uint32_t linear_chn = 0; linear_chn < nchns; linear_chn++)
        data_gen->start_values[linear_chn] =
            host_seq.start_values[m_linear_chns[linear_chn]->first - 1].f64;
    data_gen->generate_data();
}

void Backend::start(HostSeq&)
{
    // Right now the run and wait is done by the caller in MATLAB
    // so we don't need to do anything here.
}

void Backend::cancel(HostSeq&)
{
    // Ignore for now.
}

void Backend::wait(HostSeq&)
{
    // Right now the run and wait is done by the caller in MATLAB
    // so we don't need to do anything here.
}

void Backend::config(const YAML::Node &config)
{
    if (auto step_size_node = config["step_size"]) {
        m_step_size = step_size_node.as<uint32_t>();
        if (m_step_size % 2 != 0) {
            throw std::runtime_error(name() + ": step_size must be an even number");
        }
    }
    else {
        throw std::runtime_error(name() + ": Missing step_size in config");
    }
    if (auto clock_device_node = config["clock_device"]) {
        m_clock_dev = clock_device_node.as<std::string>();
    }
    else {
        throw std::runtime_error(name() + ": Missing clock_device in config");
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

NACS_EXPORT() Backend::GenerateStatus Backend::get_generate_status(uint32_t cur_seq_id) const
{
    auto &bseq = m_seqs[cur_seq_id];
    GenerateStatus res;
    res.all_time_known = bseq.all_time_known;
    res.all_val_known = bseq.all_val_known;
    if (!bseq.all_time_known || !bseq.all_val_known) {
        assert(bseq.data_gen);
        return res;
    }
    assert(!bseq.data_gen);
    res.nsamples = bseq.nsamples;
    res.first_only = bseq.first_only;
    res.fill_init.resize(m_linear_chns.size(), 0);
    for (auto [linear_chn, nele]: bseq.fill_init)
        res.fill_init[linear_chn] = nele;
    return res;
}

}

extern "C" {

using namespace NaCs;
using namespace NaCs::Seq;
using namespace NaCs::Seq::NiDAQ;

NACS_EXPORT() const double *nacs_seq_manager_expseq_get_nidaq_data(
    Manager::ExpSeq *expseq, const char *name, size_t *sz)
{
    return expseq->mgr().call_guarded([&] () -> const double* {
        auto dev = expseq->get_device(name, false);
        if (!dev) {
            Log::error("Device %s cannot be found.", name);
            return nullptr;
        }
        auto nidaq_dev = dynamic_cast<Backend*>(dev);
        if (!nidaq_dev) {
            Log::error("Device %s is not a Ni DAQ.", name);
            return nullptr;
        }
        auto idx = expseq->host_seq.cur_seq_idx();
        if (idx == uint32_t(-1)) {
            Log::error("Sequence must be started before getting data.");
            return nullptr;
        }
        return nidaq_dev->get_data(idx, sz);
    }, nullptr);
}

NACS_EXPORT() const std::pair<uint32_t,const char*>*
nacs_seq_manager_expseq_get_nidaq_channel_info(
    Manager::ExpSeq *expseq, const char *name, uint32_t *sz)
{
    static const std::pair<uint32_t,const char*> empty[1] = {};
    return expseq->mgr().call_guarded([&] () -> const std::pair<uint32_t,const char*>* {
            auto dev = expseq->get_device(name, false);
            if (!dev) {
                // Return empty array when the backend is not used.
                *sz = 0;
                return empty;
            }
            auto nidaq_dev = dynamic_cast<Backend*>(dev);
            if (!nidaq_dev) {
                Log::error("Device %s is not a Ni DAQ.", name);
                return nullptr;
            }
            return nidaq_dev->get_channel_info(sz);
        }, nullptr);
}

}
