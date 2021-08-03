//

#include "backend.h"

#include "../zynq/backend.h"

#include "../../nacs-utils/llvm/utils.h"
#include "../../nacs-utils/llvm/codegen.h"
#include "../../nacs-utils/llvm/compile.h"
#include "../../nacs-utils/llvm/execute.h"
#include "../../nacs-utils/llvm/global_rename.h"
#include "../../nacs-utils/processor.h"
#include "../../nacs-utils/streams.h"

#include <llvm/ADT/StringRef.h>
#include <llvm/IR/LegacyPassManager.h>
#include <llvm/IR/Verifier.h>
#include <llvm/Transforms/IPO.h>

#include <memory>

namespace NaCs::Seq::AWG {

static Device::Register<Backend> register_backend("awg");

enum class Backend::ChnType : uint8_t
{
    Freq,
    Amp,
    Phase
};

struct Backend::ChannelInfo {
    uint8_t m_phys_chn;
    uint32_t m_chn_num;
    uint32_t linear_idx = 0; // Linear channel number for BasicSeq's list of pulses
    Backend::ChnType m_chn_type;
    ChannelInfo(uint8_t phys_chn, uint32_t chn_num, Backend::ChnType chn_type)
        : m_phys_chn(phys_chn),
          m_chn_num(chn_num),
          m_chn_type(chn_type)
    {
    }
};

struct Backend::BasicSeq{
    struct Pulse {
        Pulse() = default;
        Pulse(uint8_t pulse_type, uint32_t id, uint32_t time_id,
              uint32_t len_id, uint32_t endvalue_id, uint32_t cond_id,
              std::string ramp_func_name)
            : pulse_type(pulse_type), id(id), time_id(time_id),
              len_id(len_id), endvalue_id(endvalue_id),
              cond_id(cond_id), ramp_func_name(ramp_func_name)
        {}
        uint8_t pulse_type;
        uint32_t id;
        uint32_t time_id;
        uint32_t len_id; // may be -1, `len == -1` <=> no-ramp pulse
        uint32_t endvalue_id;
        uint32_t cond_id;
        // double (*)(double in) for scalar function
        // void (*)(double *out, const double *in) for vector function
        //void (*ramp_func)(void);
        std::string ramp_func_name;
        //int64_t time = 0; // in sequence time unit
        //int64_t len = 0; // in sequence time unit
        //int64_t start_step = 0;
        //int64_t end_step = 0;
        friend class BasicSeq; // might not be needed
    };
    std::vector<std::vector<Pulse>> m_pulses;
    // maps from value IDs from global sequence to the one specific to this basic sequence
    std::vector<uint32_t> val_map;
    std::vector<HostSeq::Type> types;
    bool obj_file_generated = false;
    llvm::SmallVector<char,0> obj_file;
    std::string val_array_name;
    std::vector<HostSeq::Value> vals;
    // interpolated datas
    std::map<uint64_t, std::string> iData_name_map;
    // All stuff for network
    std::vector<uint8_t> interp_bc;
    std::vector<uint8_t> obj_file_bc;
    std::vector<uint8_t> types_bc;
    std::vector<uint8_t> pulses_bc;
    bool first_only = false;
    bool obj_file_sent = false;
    bool iData_sent = false;
    bool all_time_known;
    bool all_val_known;
    uint32_t n_const_map;
    uint32_t n_nonconst_map;
};

template<typename T>
void write(std::vector<uint8_t> &vec, uint32_t loc, T obj)
{
    memcpy(&vec[loc], &obj, sizeof(obj));
}

template<typename T>
uint32_t write(std::vector<uint8_t> &vec, T obj)
{
    static_assert(std::is_trivial_v<T>);
    auto len = (uint32_t)vec.size();
    vec.resize(len + sizeof(obj));
    memcpy(&vec[len], &obj, sizeof(obj));
    return len;
}

template<typename T>
size_t vector_sizeof(const typename std::vector<T>& vec)
{
    return sizeof(T) * vec.size();
}

uint32_t write(std::vector<uint8_t> &vec, const uint8_t *src, size_t sz)
{
    auto len = (uint32_t)vec.size();
    vec.resize(len + sz);
    memcpy(&vec[len], src, sz);
    return len;
}

uint32_t write(std::vector<uint8_t> &vec, char *src, size_t sz)
{
    auto len = (uint32_t)vec.size();
    vec.resize(len + sz);
    memcpy(&vec[len], src, sz);
    return len;
}


Backend::Backend(Manager &mgr, std::string name)
    : Device(mgr, std::move(name))
{
    m_seq_cnt_offset = seqcount; //seqcount is a global
}

Backend::~Backend()
{
    // object files already freed...
}

void Backend::add_channel(uint32_t chn_id, const std::string &_chn_name)
{
    // OUT1/CHN1/AMP
    llvm::StringRef chn_name(_chn_name);
    llvm::SmallVector<llvm::StringRef, 3> chn_strs;
    chn_name.split(chn_strs, "/");
    if (chn_strs.size() > 3) {
        throw std::runtime_error("Only 2 / allowed in AWG channel name");
    }
    auto phys_chn_str = chn_strs[0];
    auto chn_num_str = chn_strs[1];
    auto type_str = chn_strs[2];
    if (chn_num_str.empty())
        throw std::runtime_error("No channel number in AWG channel name.");
    if (phys_chn_str.empty())
        throw std::runtime_error("No physical output number in AWG channel name.");
    if (type_str.empty())
        throw std::runtime_error("No type for channel specified in AWG channel name");
    uint8_t phys_chn_num;
    uint32_t chn_num;
    Backend::ChnType chn_type;
    if (!phys_chn_str.startswith("OUT"))
        throw std::runtime_error("Physical channel should be specified with OUT");
    if (!chn_num_str.startswith("CHN"))
        throw std::runtime_error("Virtual channel should be specified with CHN");
    phys_chn_str = phys_chn_str.substr(3);
    chn_num_str = chn_num_str.substr(3);
    if (phys_chn_str.getAsInteger(10, phys_chn_num))
        throw std::runtime_error("Physical output for AWG must be a number.");
    if (chn_num_str.getAsInteger(10, chn_num))
        throw std::runtime_error("Channel for AWG must be a number.");
    if (type_str == "FREQ")
        chn_type = Backend::ChnType::Freq;
    else if (type_str == "AMP")
        chn_type = Backend::ChnType::Amp;
    else if (type_str == "PHASE")
        chn_type = Backend::ChnType::Phase;
    else {
        throw std::runtime_error("Unknown name for channel. Use FREQ, AMP, PHASE");
    }
    m_chn_map.try_emplace(chn_id, phys_chn_num, chn_num, chn_type);
    auto it = out_chns.begin();
    for (; it != out_chns.end(); ++it) {
        if (*it == phys_chn_num) {
            break;
        }
    }
    if (it == out_chns.end()) {
        out_chns.push_back(phys_chn_num);
    }
}

static unsigned get_vector_size()
{
    return 8;
}


void Backend::sort_channels()
{
    // assignment of linear channel id
    auto nchn = (uint32_t)m_chn_map.size();
    m_linear_chns.resize(nchn);
    auto it = m_chn_map.begin();
    for (uint32_t i = 0; i < nchn; i++)
        m_linear_chns[i] = it++;
    std::sort(m_linear_chns.begin(), m_linear_chns.end(), [&] (auto it1, auto it2) {
        auto &info1 = it1->second;
        auto &info2 = it2->second;
        if (info1.m_phys_chn == info2.m_phys_chn) {
            if (info1.m_chn_num == info2.m_chn_num) {
                return info1.m_chn_type < info2.m_chn_type;
            }
            return info1.m_chn_num < info2.m_chn_num;
        }
        return info1.m_phys_chn < info2.m_phys_chn;
    });
    std::sort(out_chns.begin(), out_chns.end());
    for (uint32_t i = 0; i < nchn; i++) {
        m_linear_chns[i]->second.linear_idx = i;
    }
}

const uint8_t *Backend::get_channel_info(uint32_t *sz)
{
    // might not be needed...
    return nullptr;
}

void Backend::prepare (Manager::ExpSeq &expseq, Compiler &compiler)
{
    reqServerInfo();
    sort_channels();
    //m_zynq_dev = dynamic_cast<Zynq::Backend*>(expseq.get_device(m_clock_dev, true));
    //if (!m_zynq_dev)
    //    throw std::runtime_error("Cannot find clock generator device");
    //m_zynq_dev->enable_clock();
    const auto bseqs = compiler.basic_seqs();
    auto &host_seq = expseq.host_seq;
    auto nbseqs = bseqs.size();
    m_seqs.resize(nbseqs);
    uint32_t nconst_map = 0;
    uint32_t n_nonconst_map_max = 0;

    // First pass: map IDs from global sequences to ones specific for the backend, in particular for each basic sequence. The constants are mapped first.
    
    //throw std::runtime_error("made it to first pass");
    for (uint32_t idx = 0; idx < nbseqs; idx++) {
        auto &bseq = bseqs[idx];
        auto &chn_infos = bseq->get_channel_infos();
        auto &host_bseq = host_seq.seqs[idx];
        auto &be_bseq = m_seqs[idx]; //backend basic seq
        bool all_time_known = true;
        bool all_val_known = true;
        bool has_pulse = false;

        be_bseq.val_map.resize(host_seq.nshared + host_bseq.nmeasure + host_bseq.ndirect + host_bseq.nneed_order, uint32_t(-1));
        uint32_t nspecific_map = 0;
        auto map_slot = [&] (uint32_t slot) {
            if (slot < host_seq.nconsts) {
                if (be_bseq.val_map[slot] == uint32_t(-1)) {
                    be_bseq.val_map[slot] = nconst_map++;
                }
                return true;
            }
            // non constant slot
            assert(slot < be_bseq.val_map.size());
            if (be_bseq.val_map[slot] == uint32_t(-1))
                be_bseq.val_map[slot] = nspecific_map++;
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
        be_bseq.n_const_map = nconst_map;
        be_bseq.n_nonconst_map = nspecific_map;
        n_nonconst_map_max = max(n_nonconst_map_max, nspecific_map);
    }

    m_nconst_vals = nconst_map;
    // Second pass: map all values and create LLVM functions
    struct RampInfo {
        uint32_t bseq_idx;
        uint32_t linear_chn;
        uint32_t pulse_idx;
        llvm::Function *ramp_func;
        std::string name {};
    };
    auto vector_size = get_vector_size();

    for (uint32_t idx = 0; idx < nbseqs; idx++) {
        // LLVM for this basic seq
        auto &env = compiler.seq()->env();
        auto env_mod = env.llvm_module();
        auto &llvm_ctx = env_mod->getContext();
        llvm::Module mod("", llvm_ctx);
        // For copying code into the new module
        LLVM::FunctionMover mover(&mod);
        // A new cgctx that emits code in the new module
        LLVM::Codegen::Context cgctx(&mod);

        
        std::vector<RampInfo> ramp_funcs;
        auto &bseq = bseqs[idx];
        auto &chn_infos = bseq->get_channel_infos();
        auto &be_bseq = m_seqs[idx];
        auto values_ty = llvm::ArrayType::get(cgctx.T_i8,
                                              (nconst_map + be_bseq.n_nonconst_map)
                                              * sizeof(double));
        auto values_gv = new llvm::GlobalVariable(mod, values_ty, false,
                                                  llvm::GlobalValue::ExternalLinkage,
                                                  llvm::ConstantAggregateZero::get(values_ty),
                                                  "values");
        auto map_slot = [&] (uint32_t slot) {
            if (slot < host_seq.nconsts) {
                assert(be_bseq.val_map[slot] != uint32_t(-1));
                return be_bseq.val_map[slot];
            }
            assert(slot < be_bseq.val_map.size());
            return be_bseq.val_map[slot] + nconst_map;
        };
        // get types
        be_bseq.types.resize(nconst_map + be_bseq.n_nonconst_map);
        auto nconsts = host_seq.nconsts;
        for (uint32_t slot = 0; slot < be_bseq.val_map.size(); slot++) {
            if (be_bseq.val_map[slot] == uint32_t(-1))
                continue;
            if (slot < nconsts)
                be_bseq.types[be_bseq.val_map[slot]] = host_seq.get_type(slot, idx);
            else {
                be_bseq.types[be_bseq.val_map[slot] + nconst_map] = host_seq.get_type(slot, idx);
            }
        }
        be_bseq.m_pulses.resize(m_chn_map.size());
        for (auto &[chn_id, be_chn_info]: m_chn_map) {
            auto it = chn_infos.find(chn_id);
            assert(it != chn_infos.end());
            auto &chn_info = it->second;
            auto &pulses = be_bseq.m_pulses[be_chn_info.linear_idx];
            for (auto &pulse: chn_info.pulses) {
                if (pulse.is_measure())
                    continue;
                auto len = pulse.len();
                uint32_t len_id = len ? map_slot(compiler.var_slot(len)) : uint32_t(-1);
                auto cond = pulse.cond();
                uint32_t cond_id = cond ? map_slot(compiler.var_slot(cond)) : uint32_t(-1);
                uint32_t time_id = map_slot(compiler.time_slot(&pulse.start()));
                uint32_t endvalue_id = map_slot(compiler.var_slot(pulse.endval()));
                bool is_fn = 0;
                bool is_vector = 0;
                if (len) {
                    is_fn = 1;
                    auto val = pulse.val();
                    auto pulse_idx = (uint32_t)pulses.size();
                    assert(val->is_call());

                    auto f = mover.clone_function(val->get_callee().llvm);
                    f = Compiler::get_ramp_func(val, f, cgctx.T_f64, cgctx.T_f64, false);
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
                        wrap_vector.add_ret_ref();
                        wrapf = cgctx.emit_wrapper(f, "0", wrap_vector);
                        is_vector = 1;
                    }
                    if (!wrapf) {
                        wrapf = cgctx.emit_wrapper(f, "0", wrap);
                        if (!wrapf) {
                            throw std::runtime_error("Failed to compile function");
                        }
                    }
                    ramp_funcs.push_back({.bseq_idx = idx, .linear_chn = be_chn_info.linear_idx,
                                          .pulse_idx = pulse_idx, .ramp_func = wrapf});
                }
                uint8_t pulse_type = get_pulse_type(be_chn_info.m_chn_type, is_fn, is_vector);
                pulses.emplace_back(pulse_type, pulse.id(), time_id,
                                    len_id, endvalue_id, cond_id, "");
            }
        }
        // Convert the value map to backend idx -> sequence idx
        // This is what we care about when we want to send the value array
        auto val_map_copy = std::move(be_bseq.val_map);
        be_bseq.val_map.resize(nconst_map + be_bseq.n_nonconst_map);
        be_bseq.vals.resize(nconst_map + be_bseq.n_nonconst_map);
        for (uint32_t slot = 0; slot < val_map_copy.size(); slot++) {
            if (val_map_copy[slot] == uint32_t(-1))
                continue;
            if (slot < host_seq.nconsts)
            {
                be_bseq.val_map[val_map_copy[slot]] = slot;
                be_bseq.vals[val_map_copy[slot]] = host_seq.values[slot];
            }
            else {
                be_bseq.val_map[val_map_copy[slot] + nconst_map] = slot;
                be_bseq.vals[val_map_copy[slot] + nconst_map].f64 = 1;
            }
        }

        // Compile all the code

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

        // INTERPOLATED DATAS HERE?
        auto it = mod.global_begin();
        for (;it != mod.global_end(); it++) {
            auto iDataName = it->getName();
            uint64_t iData_id = expseq.get_dataid(iDataName);
            if (iData_id != 0) {
                be_bseq.iData_name_map.try_emplace(iData_id, iDataName);
            }
        }
        llvm::legacy::PassManager PM;

#ifndef NDEBUG
        PM.add(llvm::createVerifierPass());
#endif
        PM.add(llvm::createGlobalDCEPass());
        // Shorten all global names since we don't care what they are
        // and this should slightly reduce the compiled binary size.
        PM.add(LLVM::createGlobalRenamePass());
        PM.run(mod);
        // The optimization passes below (in `emit_objfile`) may recreate the functions
        // so the function handle may not be valid anymore.
        // We are done with the function renaming so we can simply get the names now.
        for (auto &ramp: ramp_funcs) {
            ramp.name = ramp.ramp_func->getName();
            assert(!ramp.name.empty());
            ramp.ramp_func = nullptr;
        }

        llvm::SmallVector<char,0> vec;
        // Use my triplestr
        //printf("%s\n", m_triple_str.c_str());
        auto tgtMachine = LLVM::Compile::create_target(m_triple_str, m_cpu_str, m_feature_str);
        auto res = LLVM::Compile::emit_objfile(vec, tgtMachine.get(), &mod, true);
        if (!res)
            throw std::runtime_error("Sequence function compilation failed.");
        be_bseq.obj_file = vec;

        // Put in ramp function names
        for (auto &ramp: ramp_funcs) {
            auto &pulses = be_bseq.m_pulses[ramp.linear_chn];
            pulses[ramp.pulse_idx].ramp_func_name = ramp.name;
        }
        be_bseq.val_array_name = (values_gv->getName()).str();
        // Fill in bytecode after
        //     std::vector<uint8_t> interp_bc; non data version
        //std::vector<uint8_t> obj_file_bc;
        //std::vector<uint8_t> types_bc;
        //std::vector<uint8_t> pulses_bc;
        // obj_file_bc
        // [objfile_size: 4B][object_file: objfile_size][value_array_name: NUL-term-string]
        uint32_t sz = vec.size();
        write(be_bseq.obj_file_bc, sz);
        write(be_bseq.obj_file_bc, vec.data(), sz);
        write(be_bseq.obj_file_bc, be_bseq.val_array_name.data(), be_bseq.val_array_name.size() + 1);
        // types_bc
        // [[types:1B] x nvalues]
        write(be_bseq.types_bc, (uint8_t*) be_bseq.types.data(), be_bseq.types.size());
        // pulses_bc
        // [n_pulses: 4B][[start_time: 4B][len: 4B][endvalue: 4B][pulse_type: 1B][phys_chn_id: 1B][channel_id: 4B][funcname: NUL-term-string] x n_pulses]; Use m_linear_chns to get phys chn and chn num.
        std::vector<uint8_t> pulse_info;
        uint32_t npulses = 0;
        for (uint32_t be_chn_linear = 0; be_chn_linear < m_chn_map.size(); be_chn_linear++)
        {
            auto &pulses = be_bseq.m_pulses[be_chn_linear];
            if (pulses.empty())
                continue;
            uint32_t chn_num;
            uint8_t phys_chn_num;
            Backend::ChnType chn_type;
            get_chn_from_lin(be_chn_linear, chn_type, phys_chn_num, chn_num);
            // uint8_t pulse_type;
            //uint32_t id;
            //uint32_t time_id;
            //uint32_t len_id; // may be -1, `len == -1` <=> no-ramp pulse
            //uint32_t endvalue_id;
            //uint32_t cond_id;
            //std::string ramp_func_name;
            for (uint32_t i = 0; i < pulses.size(); i++)
            {
                npulses++;
                write(pulse_info, pulses[i].cond_id);
                write(pulse_info, pulses[i].id);
                write(pulse_info, pulses[i].time_id);
                write(pulse_info, pulses[i].len_id);
                write(pulse_info, pulses[i].endvalue_id);
                write(pulse_info, pulses[i].pulse_type);
                write(pulse_info, phys_chn_num);
                write(pulse_info, chn_num);
                write(pulse_info, pulses[i].ramp_func_name.data(), pulses[i].ramp_func_name.size() + 1);
            }
        }
        write(be_bseq.pulses_bc, npulses);
        write(be_bseq.pulses_bc, pulse_info.data(), vector_sizeof(pulse_info));
    }
    seqcount += m_seqs.size();
}
void Backend::generate (Manager::ExpSeq &expseq, Compiler &compiler)
{
}
void Backend::init_run(HostSeq &host_seq)
{
    reqServerInfo();
}
void Backend::populate_values(HostSeq &host_seq, BasicSeq &bseq)
{
    auto nvals = uint32_t(bseq.vals.size());
    for (uint32_t i = bseq.n_const_map; i < nvals; i++) {
        bseq.vals[i] = host_seq.values[bseq.val_map[i]];
    }
}
void Backend::prepare_run(HostSeq &host_seq)
{
}
void Backend::pre_run(HostSeq &host_seq)
{
    // send config if first bseq
    if (host_seq.first_bseq) {
        std::vector<uint8_t> config_msg;
        write(config_msg, (uint32_t) out_chns.size());
        write(config_msg, out_chns.data(), out_chns.size());
        auto conf_reply = m_sock->send_msg([&] (auto &sock) {
            ZMQ::send_more(sock, ZMQ::str_msg("set_out_config"));
            ZMQ::send(sock, zmq::message_t(config_msg.data(), config_msg.size()));
        }).get();
        if (conf_reply.empty())
            throw std::runtime_error("did not get reply from server");
        std::string conf_reply_str((char*) conf_reply[0].data(), conf_reply[0].size());
        if (!(conf_reply_str.compare("ok") == 0)) {
            throw std::runtime_error("config of AWG failed");
        }
    }
    // prepare code protocol, value array
    auto idx = host_seq.cur_seq_idx();
    m_cur_seq_id = idx;
    auto &bseq = m_seqs[idx];
    populate_values(host_seq, bseq);
    // send over everything, and handle any missing stuff
    bool success = false;
    while (!success)
    {
        uint32_t version = 0;
        // decide whether to send over object file and iDatas
        std::vector<uint8_t> next_msg;
        uint32_t start_trigger = 1;
        write(next_msg, start_trigger);
        if (!bseq.obj_file_sent || !bseq.iData_sent)
        {
            uint8_t is_seq_sent = 1;
            write(next_msg, is_seq_sent);
            uint32_t n_interp_data = bseq.iData_name_map.size();
            write(next_msg, n_interp_data);
            if (n_interp_data == 0) {
                bseq.iData_sent = true;
            }
            else {
                for (auto it = bseq.iData_name_map.begin(); it != bseq.iData_name_map.end(); it++)
                {
                    auto id = it->first;
                    std::string name = it->second;
                    write(next_msg, id);
                    write(next_msg, name.data(), name.size() + 1);
                    write(next_msg, (uint8_t) !bseq.iData_sent);
                    if (!bseq.iData_sent)
                    {
                        auto data_pair = mgr().get_data(id);
                        write(next_msg, (uint32_t) data_pair.second);
                        write(next_msg, (uint8_t*) data_pair.first, data_pair.second);
                    }
                }
            }
            // send object file
            write(next_msg, bseq.obj_file_bc.data(), bseq.obj_file_bc.size());
            // [nconsts: 4B][nvalues:4B]([data:8B] x nvalues)
            write(next_msg, bseq.n_const_map);
            write(next_msg, bseq.n_const_map + bseq.n_nonconst_map);
            write(next_msg, (uint8_t*) bseq.vals.data(), vector_sizeof(bseq.vals));
            // types
            write(next_msg, bseq.types_bc.data(), bseq.types_bc.size());
            // pulses
            write(next_msg, bseq.pulses_bc.data(), bseq.pulses_bc.size());
        }
        else {
            printf("Not sending object code\n");
            uint8_t is_seq_sent = 0;
            write(next_msg, is_seq_sent);
            write(next_msg, bseq.n_nonconst_map);
            write(next_msg, ((uint8_t*) bseq.vals.data()) + bseq.n_const_map * sizeof(bseq.vals[0]), bseq.n_nonconst_map * sizeof(bseq.vals[0]));
        }
        auto reply = m_sock->send_msg([&] (auto &sock) {
            ZMQ::send_more(sock, ZMQ::str_msg("run_seq"));
            ZMQ::send_more(sock, ZMQ::bits_msg(version));
            write(id_bc, 16, m_seq_cnt_offset + m_cur_seq_id);
            write(id_bc, 24, (uint8_t) host_seq.first_bseq);
            ZMQ::send_more(sock, zmq::message_t(id_bc.data(), id_bc.size()));
            ZMQ::send(sock, zmq::message_t(next_msg.data(), next_msg.size()));
        }).get();
        // handle response
        if (reply.empty())
            throw std::runtime_error("did not get reply from server");
        std::string reply_str((char*) reply[0].data(), reply[0].size());
        if (reply_str.compare("need_seq") == 0)
        {
            bseq.obj_file_sent = false;
        }
        else if (reply_str.compare("need_data") == 0)
        {
            bseq.iData_sent = false;
        }
        else if (reply_str.compare("ok") == 0)
        {
            if (reply.size() < 2) {
                throw std::runtime_error("expecting id after ok");
                // can also request id...
            }
            auto rep_data = (const uint8_t*)reply[1].data();
            memcpy(&m_cur_wait_id, rep_data, sizeof(m_cur_wait_id));
            bseq.obj_file_sent = true;
            bseq.iData_sent = true;
            success = true;
        }
        else {
            throw std::runtime_error("unknown server side error");
        }
    }
    // Sequence has been sent!
}
void Backend::start(HostSeq &host_seq)
{
    // nothing here, server should be ready for trigger...
}
void Backend::cancel(HostSeq &host_seq)
{
    // no cancel for now....
}
void Backend::wait(HostSeq &host_seq)
{
    m_sock->send_msg([&] (auto &sock) {
            ZMQ::send_more(sock, ZMQ::str_msg("wait_seq"));
            ZMQ::send(sock, ZMQ::bits_msg(m_cur_wait_id));
    }).wait();
}
void Backend::finish_run(HostSeq &host_seq)
{
}
uint8_t Backend::get_pulse_type(Backend::ChnType type, bool is_fn, bool is_vector){
    // map from chn_info and whether it's a vector to a uint8_t describing the pulse type
    // mapping from Stream.h for AWG server
    // enum class CmdType : uint8_t
//  {
    // CmdType is a enumerated class that holds all possible commands.
    // They are represented by a uint8_t
    // Meta, // Meta command types are in CmdMeta
    // AmpSet, 1
    // AmpFn, 2
    // AmpVecFn, 3
    // FreqSet, 4
    // FreqFn, 5
    // FreqVecFn, 6
    // ModChn, // add or delete channels
    // Phase,
    // _MAX = Phase // keeps track of how many CmdType options there are
// };
    //auto type = chn_info.m_chn_type;
    if (type == Backend::ChnType::Phase) {
        if (is_fn) {
            throw std::runtime_error("Cannot program a pulse on a phase");
        }
        return 8;
    }
    uint8_t base = 0;
    if (type == Backend::ChnType::Freq) {
        base = 4;
    }
    else if (type == Backend::ChnType::Amp) {
        base = 1;
    }
    if (is_fn) {
        if (is_vector) {
            return base + 2;
        }
        return base + 1;
    }
    else {
        return base;
    }
}
void Backend::get_chn_from_lin(uint32_t lin_idx, Backend::ChnType &type, uint8_t &phys_chn_id, uint32_t &chn_id)
{
    auto chn_info = m_linear_chns[lin_idx]->second;
    type = chn_info.m_chn_type;
    phys_chn_id = chn_info.m_phys_chn;
    chn_id = chn_info.m_chn_num;
}

void Backend::config(const YAML::Node &config)
{
    if (auto url_node = config["url"]) {
        auto new_url = url_node.as<std::string>();
        if (new_url != m_url) {
            m_url = std::move(new_url);
            auto &client = ZMQ::MultiClient::global();
            m_sock.reset(new ZMQ::MultiClient::SockRef(client.get_socket(m_url)));
        }
    }
    else {
        throw std::runtime_error("Missing url in config");
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

NACS_EXPORT() void Backend::reqServerInfo()
{
    // TO DO:
    // Ask server for client_id, if don't have one
    // Ask server for triple if not stored
    if (!m_sock)
        throw std::runtime_error("Backend socket not configured");
    auto it = info_map.find(m_url);
    if (it != info_map.end()) {
        auto global_info = it->second;
        m_client_id = global_info.client_id;
        m_server_id = global_info.server_id;
        m_triple_str = global_info.triple_str;
        m_cpu_str = global_info.cpu_str;
        m_feature_str = global_info.feature_str;
    }
    else {
        auto reply = m_sock->send_msg([&] (auto &sock) {
            ZMQ::send(sock, ZMQ::str_msg("req_client_id"));
        }).get();
        if (reply.empty())
            throw std::runtime_error("client id not obtained");
        auto rep_data = (const uint8_t*)reply[0].data();
        memcpy(&m_client_id, rep_data, sizeof(m_client_id));
        reply = m_sock->send_msg([&] (auto &sock) {
            ZMQ::send(sock, ZMQ::str_msg("req_server_id"));
        }).get();
        if (reply.empty())
            throw std::runtime_error("server id not obtained");
        rep_data = (const uint8_t*)reply[0].data();
        memcpy(&m_server_id, rep_data, sizeof(m_server_id));
        reply = m_sock->send_msg([&] (auto &sock) {
            ZMQ::send(sock, ZMQ::str_msg("req_triple"));
        }).get();
        if (reply.empty())
            throw std::runtime_error("triple not obtained");
        std::string triple_str((char*) reply[0].data(), reply[0].size());
        std::string cpu_str((char*) reply[1].data(), reply[1].size());
        std::string feature_str((char*) reply[2].data(), reply[2].size());
        m_triple_str = triple_str;
        m_cpu_str = cpu_str;
        m_feature_str = feature_str;
        info_map.try_emplace(m_url, m_client_id, m_server_id, m_triple_str, m_cpu_str, m_feature_str);
        // fill in id_bc
        id_bc.resize(25); // 8 bytes server_id, 8 bytes_client_id, 8 bytes, 8 bytes seq_id, 1 byte is_first_seq, Fill in first 16 bytes
        write(id_bc, (uint32_t) 0, m_server_id);
        write(id_bc, (uint32_t) 8, m_client_id);
    }
}

}
