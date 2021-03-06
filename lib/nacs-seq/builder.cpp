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

#include "builder.h"
#include "error.h"
#include "seq.h"

#include "../nacs-utils/llvm/codegen.h"
#include "../nacs-utils/log.h"
#include "../nacs-utils/mem.h"

namespace NaCs::Seq {

NACS_EXPORT_ Builder::Builder(Seq &seq)
    : m_seq(seq)
{
}

NACS_EXPORT() Builder::~Builder()
{
}

NACS_EXPORT() Var *Builder::arg_to_var(ArgType argtype, NodeArg arg)
{
    switch (argtype) {
    case ArgType::ConstBool:
        return m_seq.get_const(IR::TagVal(arg.b));
    case ArgType::ConstInt32:
        return m_seq.get_const(IR::TagVal(arg.i32));
    case ArgType::ConstFloat64:
        return m_seq.get_const(IR::TagVal(arg.f64));
    case ArgType::Node: {
        auto node = get_node(arg.id);
        (void)node;
        assert(node->flag == NodeFlag::Free);
        return m_varmap[arg.id - 1];
    }
    case ArgType::Measure: {
        auto id = arg.id;
        auto it = m_measure_vars.find(id);
        if (it == m_measure_vars.end()) {
            if (m_disabled_measures.insert(id).second)
                Log::warn("Measurement %d taken on disabled channel.\n", id);
            return m_seq.get_const(IR::TagVal(false));
        }
        return it->second;
    }
    case ArgType::Global: {
        auto id = arg.id;
        if (id >= slots.size())
            throw std::runtime_error("Invalid global ID");
        return m_seq.get_slot(slots[id], id);
    }
    case ArgType::Arg:
    default:
        return nullptr;
    }
}

NACS_EXPORT() Arg Builder::arg_to_callarg(ArgType argtype, NodeArg arg)
{
    switch (argtype) {
    case ArgType::ConstBool:
        return Arg::create_const(IR::TagVal(arg.b));
    case ArgType::ConstInt32:
        return Arg::create_const(IR::TagVal(arg.i32));
    case ArgType::ConstFloat64:
        return Arg::create_const(IR::TagVal(arg.f64));
    default:
        return Arg::create_var(arg_to_var(argtype, arg));
    }
}

void Builder::prescan_seq()
{
    auto nchn = chnnames.size();
    for (size_t i = 0; i < nchn; i++) {
        auto id = m_seq.get_chn_id(chnnames[i], true);
        if (id != i + 1) {
            throw std::runtime_error("Duplicated channel name");
        }
    }

    auto check_chnid = [&] (uint32_t chn) {
        if (chn == 0 || chn > nchn) {
            throw std::runtime_error("Invalid channel ID");
        }
    };
    auto check_branch = [&] (uint32_t branch) {
        if (branch > seqs.size()) {
            // 0 or 1-n are valid branch target ID.
            throw std::runtime_error("Invalid branch target");
        }
    };

    for (auto &[chn, val]: defvals) {
        check_chnid(chn);
        m_seq.set_defval(chn, val.get<double>());
    }
    for (auto chn: noramp_chns) {
        check_chnid(chn);
        m_mayramp[chn - 1] = false;
    }
    // Create slot variables.
    auto nslots = slots.size();
    for (size_t i = 0; i < nslots; i++)
        m_seq.get_slot(slots[i], i);

    auto nseq = seqs.size();
    for (uint32_t sid = 1; sid <= nseq; sid++) {
        auto &seq = seqs[sid - 1];
        auto &bseq = seq.m_basicseq;
        bseq = m_seq.add_basicseq(sid);
        for (auto &[oid, output]: seq.outputs)
            check_chnid(output.chn);
        for (auto &[mid, measure]: seq.measures) {
            check_chnid(measure.chn);
            auto &mvar = m_measure_vars[mid];
            if (mvar)
                throw std::runtime_error("Duplicated measure ID");
            mvar = bseq->new_measure(m_seq.env(), mid);
        }
        for (auto &br: seq.branches)
            check_branch(br.target_id);
        check_branch(seq.default_target);
    }
}

namespace {

struct NodeCacheKey {
    Builder::OpCode op;
    uint32_t data_id;
    Arg args[3];
    bool operator<(const NodeCacheKey &other) const
    {
        if (uint8_t(op) != uint8_t(other.op))
            return uint8_t(op) < uint8_t(other.op);
        auto narg = Builder::node_narg(op);
        for (int8_t i = 0; i < narg; i++) {
            auto cmp = args[i].compare(other.args[i]);
            if (cmp != 0) {
                return cmp < 0;
            }
        }
        if (op == Builder::OpCode::Interp && data_id != other.data_id)
            return data_id < other.data_id;
        return false;
    }
};

}

void Builder::scan_node()
{
    // For each node,
    // 1. Check if their arguments are valid.
    // 2. Check if they refer to any `Arg`.
    // 3. For ones that doesn't refer to `Arg`, create a corresponding `Var`.

    struct Visitor : NodeVisitor {
        using NodeVisitor::NodeVisitor;

        bool previsit(Node *node)
        {
            if (node->flag == NodeFlag::Scanning)
                throw std::runtime_error("Dependency loop between values");
            if (node->flag != NodeFlag::None)
                return false;
            node->flag = NodeFlag::Scanning;
            if (node->op == OpCode::Interp && node->data_id >= builder.datas.size())
                throw std::runtime_error("Invalid interp data ID");
            return true;
        }
        void postvisit(Node *node)
        {
            auto op = node->op;
            int8_t narg = node_narg(op);
            assert(narg <= 3);
            bool refarg = false;
            for (int8_t i = 0; i < narg; i++) {
                if (node->argtypes[i] == ArgType::Arg) {
                    refarg = true;
                    break;
                }
                else if (auto child = builder.get_argref(node->argtypes[i], node->args[i])) {
                    assert(child->flag == NodeFlag::Free || child->flag == NodeFlag::RefArg);
                    if (child->flag == NodeFlag::RefArg) {
                        refarg = true;
                        break;
                    }
                }
            }
            node->flag = refarg ? NodeFlag::RefArg : NodeFlag::Free;
            if (refarg)
                return;
            auto nodeid = builder.get_node_id(node);
            builder.m_varmap[nodeid - 1] = create_var(node);
        }

        Var *create_var(Node *node)
        {
            if (node->op == OpCode::Identity)
                return builder.arg_to_var(node->argtypes[0], node->args[0]);
            int8_t narg = node_narg(node->op);
            NodeCacheKey key{node->op, node->data_id, {}};
            std::vector<Arg> call_args(narg);
            for (int8_t i = 0; i < narg; i++)
                key.args[i] = call_args[i] =
                    builder.arg_to_callarg(node->argtypes[i], node->args[i]);
            auto &var = cache[key];
            if (var)
                return var;
            int32_t args[] = {0, 1, 2};
            std::vector<IR::Type> types(narg, IR::Type::Float64);
            IR::Builder irbuilder(IR::Type::Float64, types);
            irbuilder.createRet(builder.create_inst(irbuilder, node, args));
            var = builder.m_seq.get_call(irbuilder.get(), call_args, 0);

            switch (node->op) {
            case OpCode::Add:
            case OpCode::Mul:
            case OpCode::CmpNE:
            case OpCode::CmpEQ:
            case OpCode::And:
            case OpCode::Or:
            case OpCode::Xor:
            case OpCode::Hypot:
            case OpCode::Max:
            case OpCode::Min:
                std::swap(key.args[0], key.args[1]);
                cache[key] = var;
                break;
            case OpCode::CmpLT:
                key.op = OpCode::CmpGT;
                std::swap(key.args[0], key.args[1]);
                cache[key] = var;
                break;
            case OpCode::CmpGT:
                key.op = OpCode::CmpLT;
                std::swap(key.args[0], key.args[1]);
                cache[key] = var;
                break;
            case OpCode::CmpLE:
                key.op = OpCode::CmpGE;
                std::swap(key.args[0], key.args[1]);
                cache[key] = var;
                break;
            case OpCode::CmpGE:
                key.op = OpCode::CmpLE;
                std::swap(key.args[0], key.args[1]);
                cache[key] = var;
                break;
            default:
                break;
            }
            return var;
        }

        // Use the cache to perform CSE
        std::map<NodeCacheKey,Var*> cache;
    };

    Visitor visitor(*this);
    for (auto &node: nodes) {
        if (!visitor.previsit(&node))
            continue;
        visit_dfs(&node, visitor);
    }
}

void Builder::scan_seq()
{
    // Verify type of node
    // Create assumptions
    // Create times
    // Create pulses and measures
    auto nseq = seqs.size();
    for (uint32_t sid = 1; sid <= nseq; sid++) {
        auto &seq = seqs[sid - 1];
        NaCs::Seq::BasicSeq *bseq = seq.m_basicseq;
        uint32_t ntimes = seq.times.size();
        for (uint32_t time_id = 1; time_id <= ntimes; time_id++) {
            auto et = get_eventtime(time_id, seq);
            assert(!et->terms.empty());
            auto &term = et->terms.back();
            if (term.sign == Sign::Unknown)
                continue;
            bseq->add_assume(term.sign, term.var.get(), term.id);
        }
        for (auto time: seq.endtimes)
            bseq->add_endtime(*get_eventtime(time, seq));
        for (auto timeorder: seq.timeorders)
            bseq->add_time_order(get_eventtime(timeorder.before_time_id, seq),
                                 get_eventtime(timeorder.after_time_id, seq),
                                 timeorder.sign != Sign::Pos);
        for (auto &assign: seq.assignments) {
            if (get_node(assign.val_node)->flag != NodeFlag::Free)
                throw std::runtime_error("Invalid assignment value");
            if (assign.global_id >= slots.size())
                throw std::runtime_error("Invalid global ID");
            auto val = m_varmap[assign.val_node - 1];
            assert(val);
            bseq->assign_global(assign.global_id, val, assign.assign_id);
        }
        for (auto &[mid, measure]: seq.measures) {
            auto mvar = m_measure_vars.find(mid)->second;
            bseq->add_measure(measure.chn, mid,
                              *get_eventtime(measure.time_id, seq), mvar);
        }
        for (auto &[oid, output]: seq.outputs) {
            auto &time = *get_eventtime(output.time_id, seq);
            Var *len = nullptr;
            if (output.len_node) {
                if (get_node(output.len_node)->flag != NodeFlag::Free)
                    throw std::runtime_error("Invalid pulse length value");
                len = m_varmap[output.len_node - 1];
            }
            Var *cond = nullptr;
            if (output.cond_node) {
                if (get_node(output.cond_node)->flag != NodeFlag::Free)
                    throw std::runtime_error("Invalid pulse condition value");
                cond = m_varmap[output.cond_node - 1];
            }
            auto val_node = get_node(output.val_node);
            Var *val;
            if (val_node->flag == NodeFlag::Free) {
                val = m_varmap[output.val_node - 1];
                len = nullptr;
            }
            else {
                bool isramp;
                std::tie(val, isramp) = node_to_rampvar(val_node);
                if (!isramp) {
                    len = nullptr;
                }
                else if (!m_mayramp[output.chn - 1]) {
                    throw Error(Error::Type::Pulse, Error::Pulse::NoRamp,
                                Error::Type::Pulse, oid, Error::Type::Channel, output.chn,
                                "Ramp not supported on channel");
                }
                else if (!len) {
                    throw std::runtime_error("Invalid output value");
                }
            }
            bseq->add_pulse(output.chn, oid, time, len, val, cond);
        }

        for (auto &br: seq.branches) {
            if (get_node(br.cond_node)->flag != NodeFlag::Free)
                throw std::runtime_error("Invalid branch condition");
            auto val = m_varmap[br.cond_node - 1];
            auto target = br.target_id ? seqs[br.target_id - 1].m_basicseq : nullptr;
            bseq->add_branch(val, target, br.branch_id);
        }
        auto default_target = seq.default_target ?
            seqs[seq.default_target - 1].m_basicseq : nullptr;
        bseq->set_default_branch(default_target);
    }
}

NACS_EXPORT() void Builder::buildseq()
{
    if (seqs.empty()) {
        auto &seq = seqs.emplace_back();
        seq.default_target = 0;
    }
    m_varmap.resize(nodes.size(), nullptr);
    m_mayramp.resize(chnnames.size(), true);
    for (auto &seq: seqs)
        seq.m_eventtimes.resize(seq.times.size(), nullptr);
    prescan_seq();
    scan_node();
    scan_seq();
}

namespace {
EventTime invalid_time;
}

NACS_EXPORT() EventTime *Builder::get_eventtime(uint32_t time_id, BasicSeq &seq)
{
    auto time = seq.get_time(time_id);
    auto &et = seq.m_eventtimes[time_id - 1];
    assert(et != &invalid_time);
    if (et)
        return et;

    struct TimeFieldIterator : FieldIterator<Time*> {
        TimeFieldIterator(Time *time, BasicSeq &seq)
            : m_time(time),
              m_seq(seq)
        {
        }
        TimeFieldIterator(const TimeFieldIterator &other)
            : m_time(other.m_time),
              m_seq(other.m_seq),
              m_end(other.m_end)
        {
        }
        Time *get() const
        {
            return m_time->prev_id ? m_seq.get_time(m_time->prev_id) : nullptr;
        }
        Time *parent() const
        {
            return m_time;
        }
        bool is_end() const
        {
            return m_end;
        }
        FieldIterator &operator++()
        {
            m_end = true;
            return *this;
        }
        TimeFieldIterator &operator=(const TimeFieldIterator &other)
        {
            assert(&m_seq == &other.m_seq);
            m_time = other.m_time;
            m_end = other.m_end;
            return *this;
        }

    private:
        Time *m_time;
        BasicSeq &m_seq;
        bool m_end = false;
    };
    struct Visitor : DFSVisitor<Time*,TimeFieldIterator,StackVector> {
        Visitor(Builder &builder, BasicSeq &seq)
            : builder(builder),
              seq(seq)
        {
        }

        bool previsit(Time *time)
        {
            auto &et = seq.m_eventtimes[seq.get_time_id(time) - 1];
            if (!et) {
                // Mark this as being processed.
                et = &invalid_time;
                return true;
            }
            if (et == &invalid_time)
                throw std::runtime_error("Dependency loop between times");
            return false;
        }
        void postvisit(Time *time)
        {
            if (time->prev_id)
                assert(seq.m_eventtimes[time->prev_id - 1]);
            EventTime *prev = time->prev_id ? seq.m_eventtimes[time->prev_id - 1] : nullptr;
            assert(prev != &invalid_time);
            if (builder.get_node(time->delta_node)->flag != NodeFlag::Free)
                throw std::runtime_error("Invalid time node");
            // `prev ? *prev : EventTime(0)` doesn't seem to work with clang
            // even though GCC accepts it. Not sure which one is correct.
            auto &et = prev ? seq.m_basicseq->track_time(*prev) :
                seq.m_basicseq->track_time(EventTime(0));
            auto delta = builder.m_varmap[time->delta_node - 1];
            assert(delta);
            et.add_term(time->sign, delta, time->id);
            seq.m_eventtimes[seq.get_time_id(time) - 1] = &et;
        }
        TimeFieldIterator begin(Time *time) const
        {
            return TimeFieldIterator(time, seq);
        }

        Builder &builder;
        BasicSeq &seq;
    };

    visit_dfs(time, Visitor(*this, seq));
    assert(et);
    return et;
}

NACS_EXPORT() int32_t Builder::create_inst(IR::Builder &builder, Node *node, int32_t *args)
{
    auto op = node->op;
    switch (op) {
    case OpCode::Add: return builder.createAdd(args[0], args[1]);
    case OpCode::Sub: return builder.createSub(args[0], args[1]);
    case OpCode::Mul: return builder.createMul(args[0], args[1]);
    case OpCode::Div: return builder.createFDiv(args[0], args[1]);
    case OpCode::CmpLT: return builder.createCmp(IR::CmpType::lt, args[0], args[1]);
    case OpCode::CmpGT: return builder.createCmp(IR::CmpType::gt, args[0], args[1]);
    case OpCode::CmpLE: return builder.createCmp(IR::CmpType::le, args[0], args[1]);
    case OpCode::CmpGE: return builder.createCmp(IR::CmpType::ge, args[0], args[1]);
    case OpCode::CmpNE: return builder.createCmp(IR::CmpType::ne, args[0], args[1]);
    case OpCode::CmpEQ: return builder.createCmp(IR::CmpType::eq, args[0], args[1]);
    case OpCode::And: return builder.createAnd(args[0], args[1]);
    case OpCode::Or: return builder.createOr(args[0], args[1]);
    case OpCode::Xor: return builder.createXor(args[0], args[1]);
    case OpCode::Not: return builder.createNot(args[0]);
    case OpCode::Abs: return builder.createCall(IR::Builtins::abs, {args[0]});
    case OpCode::Ceil: return builder.createCall(IR::Builtins::ceil, {args[0]});
    case OpCode::Exp: return builder.createCall(IR::Builtins::exp, {args[0]});
    case OpCode::Expm1: return builder.createCall(IR::Builtins::expm1, {args[0]});
    case OpCode::Floor: return builder.createCall(IR::Builtins::floor, {args[0]});
    case OpCode::Log: return builder.createCall(IR::Builtins::log, {args[0]});
    case OpCode::Log1p: return builder.createCall(IR::Builtins::log1p, {args[0]});
    case OpCode::Log2: return builder.createCall(IR::Builtins::log2, {args[0]});
    case OpCode::Log10: return builder.createCall(IR::Builtins::log10, {args[0]});
    case OpCode::Pow: return builder.createCall(IR::Builtins::pow, {args[0], args[1]});
    case OpCode::Sqrt: return builder.createCall(IR::Builtins::sqrt, {args[0]});
    case OpCode::Asin: return builder.createCall(IR::Builtins::asin, {args[0]});
    case OpCode::Acos: return builder.createCall(IR::Builtins::acos, {args[0]});
    case OpCode::Atan: return builder.createCall(IR::Builtins::atan, {args[0]});
    case OpCode::Atan2: return builder.createCall(IR::Builtins::atan2, {args[0], args[1]});
    case OpCode::Asinh: return builder.createCall(IR::Builtins::asinh, {args[0]});
    case OpCode::Acosh: return builder.createCall(IR::Builtins::acosh, {args[0]});
    case OpCode::Atanh: return builder.createCall(IR::Builtins::atanh, {args[0]});
    case OpCode::Sin: return builder.createCall(IR::Builtins::sin, {args[0]});
    case OpCode::Cos: return builder.createCall(IR::Builtins::cos, {args[0]});
    case OpCode::Tan: return builder.createCall(IR::Builtins::tan, {args[0]});
    case OpCode::Sinh: return builder.createCall(IR::Builtins::sinh, {args[0]});
    case OpCode::Cosh: return builder.createCall(IR::Builtins::cosh, {args[0]});
    case OpCode::Tanh: return builder.createCall(IR::Builtins::tanh, {args[0]});
    case OpCode::Hypot: return builder.createCall(IR::Builtins::hypot, {args[0], args[1]});
    case OpCode::Erf: return builder.createCall(IR::Builtins::erf, {args[0]});
    case OpCode::Erfc: return builder.createCall(IR::Builtins::erfc, {args[0]});
    case OpCode::Gamma: return builder.createCall(IR::Builtins::gamma, {args[0]});
    case OpCode::Lgamma: return builder.createCall(IR::Builtins::lgamma, {args[0]});
    case OpCode::Rint: return builder.createCall(IR::Builtins::rint, {args[0]});
    case OpCode::Max: return builder.createCall(IR::Builtins::max, {args[0], args[1]});
    case OpCode::Min: return builder.createCall(IR::Builtins::min, {args[0], args[1]});
    case OpCode::Mod: return builder.createCall(IR::Builtins::mod, {args[0], args[1]});
    case OpCode::Interp: {
        assert(node->data_id < datas.size());
        auto &data = datas[node->data_id];
        return builder.createInterp(args[0], args[1], args[2], data.size(), data.data());
    }
    case OpCode::Select: return builder.createSelect(args[0], args[1], args[2]);
    case OpCode::Identity: return args[0];
    default: return 0;
    }
}

std::pair<Var*,bool> Builder::node_to_rampvar(Node *node)
{
    assert(node->flag == NodeFlag::RefArg);

    llvm::SmallVector<Arg,6> callargs;
    callargs.push_back(Arg::create_arg(0));
    callargs.push_back(Arg::create_arg(1));
    std::map<Node*,int32_t> irmap;
    std::map<uint32_t,int32_t> measuremap;
    std::map<uint32_t,int32_t> globalmap;

    struct PreVisitor : NodeVisitor {
        PreVisitor(Builder &builder, std::map<Node*,int32_t> &irmap,
                   std::map<uint32_t,int32_t> &measuremap,
                   std::map<uint32_t,int32_t> &globalmap,
                   llvm::SmallVector<Arg,6> &callargs)
            : NodeVisitor(builder),
              irmap(irmap),
              measuremap(measuremap),
              globalmap(globalmap),
              callargs(callargs)
        {
        }

        Arg &alloc_callarg(int32_t &id)
        {
            id = callargs.size();
            return callargs.emplace_back();
        }

        bool previsit(Node *node)
        {
            if (node->flag == NodeFlag::Free) {
                auto &id = irmap[node];
                if (!id) {
                    auto &carg = alloc_callarg(id);
                    carg = Arg::create_var(builder.m_varmap[builder.get_node_id(node) - 1]);
                }
                return false;
            }
            assert(node->flag == NodeFlag::RefArg);
            auto op = node->op;
            int8_t narg = node_narg(op);
            for (int8_t i = 0; i < narg; i++) {
                auto &arg = node->args[i];
                std::map<uint32_t,int32_t> *map;
                Var *var;
                switch (node->argtypes[i]) {
                case ArgType::Measure: {
                    map = &measuremap;
                    auto measure_it = builder.m_measure_vars.find(arg.id);
                    if (measure_it == builder.m_measure_vars.end()) {
                        if (builder.m_disabled_measures.insert(arg.id).second)
                            Log::warn("Measurement %d taken on disabled channel.\n", arg.id);
                        continue;
                    }
                    var = measure_it->second;
                    break;
                }
                case ArgType::Global:
                    map = &globalmap;
                    if (arg.id >= builder.slots.size())
                        throw std::runtime_error("Invalid global ID");
                    var = builder.m_seq.get_slot(builder.slots[arg.id], arg.id);
                    break;
                case ArgType::Arg:
                    if (arg.id >= 2)
                        throw std::runtime_error("Too many arguments for ramp");
                    argused[arg.id] = true;
                    [[fallthrough]];
                default:
                    continue;
                }
                assert(map);
                auto &id = (*map)[arg.id];
                if (!id) {
                    auto &carg = alloc_callarg(id);
                    carg = Arg::create_var(var);
                }
            }
            return true;
        }

        std::map<Node*,int32_t> &irmap;
        std::map<uint32_t,int32_t> &measuremap;
        std::map<uint32_t,int32_t> &globalmap;
        llvm::SmallVector<Arg,6> &callargs;
        bool argused[2] = {false, false};
    };

    PreVisitor previsitor(*this, irmap, measuremap, globalmap, callargs);
    if (previsitor.previsit(node))
        visit_dfs(node, previsitor);
    if (!previsitor.argused[0])
        callargs[0] = Arg::create_const(false);
    if (!previsitor.argused[1])
        callargs[1] = Arg::create_const(false);

    struct Visitor : NodeVisitor {
        Visitor(Builder &builder, std::map<Node*,int32_t> &irmap,
                std::map<uint32_t,int32_t> &measuremap,
                std::map<uint32_t,int32_t> &globalmap,
                llvm::SmallVector<Arg,6> &callargs)
            : NodeVisitor(builder),
              irmap(irmap),
              measuremap(measuremap),
              globalmap(globalmap),
              callargs(callargs),
              irbuilder{IR::Type::Float64, std::vector<IR::Type>(callargs.size(),
                                                                 IR::Type::Float64)}
        {
        }

        bool previsit(Node *node)
        {
            if (node->flag == NodeFlag::Free)
                return false;
            assert(node->flag == NodeFlag::RefArg);
            auto it = irmap.find(node);
            if (it != irmap.end())
                return false;
            return true;
        }

        int32_t arg_to_ir(ArgType argtype, NodeArg arg)
        {
            switch (argtype) {
            case ArgType::ConstBool:
                return irbuilder.getConst(IR::TagVal(arg.b));
            case ArgType::ConstInt32:
                return irbuilder.getConst(IR::TagVal(arg.i32));
            case ArgType::ConstFloat64:
                return irbuilder.getConstFloat(arg.f64);
            case ArgType::Node:
                return irmap.find(builder.get_node(arg.id))->second;
            case ArgType::Measure: {
                auto measure_it = measuremap.find(arg.id);
                if (measure_it == measuremap.end())
                    return irbuilder.getConst(IR::TagVal(false));
                return measure_it->second;
            }
            case ArgType::Global:
                return globalmap.find(arg.id)->second;
            case ArgType::Arg:
                assert(arg.id < 2);
                return arg.id;
            default:
                throw std::runtime_error("Invalid node argument type");
            }
        }

        void postvisit(Node *node)
        {
            assert(node->flag == NodeFlag::RefArg);
            auto op = node->op;
            int8_t narg = node_narg(op);
            std::vector<int32_t> args;
            for (int8_t i = 0; i < narg; i++)
                args.push_back(arg_to_ir(node->argtypes[i], node->args[i]));
            irmap.emplace(node, builder.create_inst(irbuilder, node, args.data()));
        }

        std::map<Node*,int32_t> &irmap;
        std::map<uint32_t,int32_t> &measuremap;
        std::map<uint32_t,int32_t> &globalmap;
        llvm::SmallVector<Arg,6> &callargs;
        IR::Builder irbuilder;
    };

    Visitor visitor(*this, irmap, measuremap, globalmap, callargs);
    assert(visitor.previsit(node));
    visit_dfs(node, visitor);
    auto ir_it = irmap.find(node);
    assert(ir_it != irmap.end());
    visitor.irbuilder.createRet(ir_it->second);
    return {m_seq.get_call(visitor.irbuilder.get(), callargs, 2), previsitor.argused[0]};
}

NACS_EXPORT() void Builder::deserialize(const uint8_t *data, size_t size)
{
    // Assume little endian...
    Mem::Reader reader(data, size);

    // NodeArg:
    //   Node ID:
    //     1-based index into the node array. (for use in sequence)
    //   Global slot ID: 0-based index into the slot types array.
    //   Measure ID: globally unique ID
    //   Arg: 0-based index of argument
    // Data ID: 0-based index into the data array
    // Channel ID: 1-based index into the channel name array
    // Time ID: 1-based index into the time array
    // Basic Sequence ID: 1-based index into the bseq array. 0 means end of sequence.

    // Format:
    //   [version <0>: 1B]
    //   [nnodes: 4B]
    //   [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //     [OpCode <Interp>: 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    //   [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    //   [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    //   [nslots: 4B][[Type: 1B] x nslots]
    //   [nnoramp: 4B][[chnid: 4B] x nnoramp]
    //   [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    //   [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
    //   [nbackenddatas: 4B][[device_name: NUL-terminated string]
    //                       [size: 4B][data: size B] x nbackenddatas]
    //
    // Basic Sequence format:
    //   [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
    //   [nendtimes: 4B][[time_id: 4B] x nendtimes]
    //   [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
    //   [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
    //   [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
    //   [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
    //   [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
    //   [default_target: 4B]

    // [version <0>: 1B]
    if (reader.read<int8_t>() != 0)
        throw std::runtime_error("Unknown sequence serialization version");
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    nodes.resize(reader.read<uint32_t>());
    for (auto &node: nodes) {
        auto _op = reader.read<uint8_t>();
        if (!_op || _op > uint8_t(OpCode::_MaxOP))
            throw std::runtime_error("Invalid Node OpCode");
        auto op = OpCode(_op);
        auto narg = node_narg(op);
        node.op = op;
        for (int i = 0; i < narg; i++) {
            auto _atype = reader.read<uint8_t>();
            if (_atype > uint8_t(ArgType::_MaxType))
                throw std::runtime_error("Invalid Node argument type");
            auto atype = ArgType(_atype);
            node.argtypes[i] = atype;
            auto &arg = node.args[i];
            switch (atype) {
            case ArgType::ConstBool:
                arg.b = reader.read<uint8_t>() != 0;
                break;
            case ArgType::ConstInt32:
                arg.i32 = reader.read<int32_t>();
                break;
            case ArgType::ConstFloat64:
                arg.f64 = reader.read<double>();
                break;
            case ArgType::Measure:
            case ArgType::Node:
            case ArgType::Global:
            case ArgType::Arg:
                arg.id = reader.read<uint32_t>();
                break;
            default:
                throw std::runtime_error("Invalid Node argument type");
            }
        }
        if (op == OpCode::Interp) {
            node.data_id = reader.read<uint32_t>();
        }
    }
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    chnnames.resize(reader.read<uint32_t>());
    for (auto &chnname: chnnames) {
        auto [name, sz] = reader.read_string();
        chnname.append(name, sz);
    }
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    auto ndefvals = reader.read<uint32_t>();
    for (uint32_t i = 0; i < ndefvals; i++) {
        auto chnid = reader.read<uint32_t>();
        auto _type = reader.read<uint8_t>();
        switch (_type) {
        case uint8_t(IR::Type::Bool):
            defvals[chnid] = IR::TagVal(reader.read<uint8_t>() != 0);
            break;
        case uint8_t(IR::Type::Int32):
            defvals[chnid] = IR::TagVal(reader.read<int32_t>());
            break;
        case uint8_t(IR::Type::Float64):
            defvals[chnid] = IR::TagVal(reader.read<double>());
            break;
        }
    }
    // [nslots: 4B][[Type: 1B] x nslots]
    slots.resize(reader.read<uint32_t>());
    for (auto &slot: slots) {
        auto _type = reader.read<uint8_t>();
        if (_type != uint8_t(IR::Type::Bool) && _type != uint8_t(IR::Type::Int32) &&
            _type != uint8_t(IR::Type::Float64))
            throw std::runtime_error("Invalid value type");
        slot = IR::Type(_type);
    }
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    noramp_chns.resize(reader.read<uint32_t>());
    for (auto &chn: noramp_chns)
        chn = reader.read<uint32_t>();
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    auto nbasicseqs = reader.read<uint32_t>();
    seqs.resize(nbasicseqs);
    for (auto &seq: seqs) {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        seq.times.resize(reader.read<uint32_t>());
        for (auto &time: seq.times) {
            auto sign = reader.read<uint8_t>();
            if (sign != uint8_t(Sign::Unknown) && sign != uint8_t(Sign::NonNeg) &&
                sign != uint8_t(Sign::Pos))
                throw std::runtime_error("Invalid sign specification");
            time.sign = Sign(sign);
            time.id = reader.read<uint32_t>();
            time.delta_node = reader.read<uint32_t>();
            time.prev_id = reader.read<uint32_t>();
        }
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        seq.endtimes.resize(reader.read<uint32_t>());
        for (auto &endtime: seq.endtimes)
            endtime = reader.read<uint32_t>();
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        seq.timeorders.resize(reader.read<uint32_t>());
        for (auto &timeorder: seq.timeorders) {
            auto sign = reader.read<uint8_t>();
            if (sign != uint8_t(Sign::NonNeg) && sign != uint8_t(Sign::Pos))
                throw std::runtime_error("Invalid sign specification");
            timeorder.sign = Sign(sign);
            timeorder.id = reader.read<uint32_t>();
            timeorder.before_time_id = reader.read<uint32_t>();
            timeorder.after_time_id = reader.read<uint32_t>();
        }
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        auto noutputs = reader.read<uint32_t>();
        for (uint32_t i = 0; i < noutputs; i++) {
            auto &output = seq.outputs[reader.read<uint32_t>()];
            output.time_id = reader.read<uint32_t>();
            output.len_node = reader.read<uint32_t>();
            output.val_node = reader.read<uint32_t>();
            output.cond_node = reader.read<uint32_t>();
            output.chn = reader.read<uint32_t>();
        }
        // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
        auto nmeasures = reader.read<uint32_t>();
        for (uint32_t i = 0; i < nmeasures; i++) {
            auto &measure = seq.measures[reader.read<uint32_t>()];
            measure.time_id = reader.read<uint32_t>();
            measure.chn = reader.read<uint32_t>();
        }
        // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
        seq.assignments.resize(reader.read<uint32_t>());
        for (auto &assign: seq.assignments) {
            assign.assign_id = reader.read<uint32_t>();
            assign.global_id = reader.read<uint32_t>();
            assign.val_node = reader.read<uint32_t>();
        }
        // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
        seq.branches.resize(reader.read<uint32_t>());
        for (auto &branch: seq.branches) {
            branch.branch_id = reader.read<uint32_t>();
            branch.target_id = reader.read<uint32_t>();
            branch.cond_node = reader.read<uint32_t>();
        }
        // [default_target: 4B]
        seq.default_target = reader.read<uint32_t>();
    }
    // [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
    datas.resize(reader.read<uint32_t>());
    for (auto &data: datas) {
        auto sz = reader.read<uint32_t>();
        auto p = reader.read_array<double>(sz);
        data.assign(p, p + sz);
    }
    // [nbackenddatas: 4B][[device_name: NUL-terminated string]
    //                     [size: 4B][data: size B] x nbackenddatas]
    auto nbackenddatas = reader.read<uint32_t>();
    for (uint32_t i = 0; i < nbackenddatas; i++) {
        auto [name, name_sz] = reader.read_string();
        std::string chnname(name, name_sz);
        auto sz = reader.read<uint32_t>();
        auto data = reader.read_array<uint8_t>(sz);
        backend_datas.try_emplace(std::move(chnname), data, data + sz);
    }
}

void Builder::print_node_arg(std::ostream &stm, ArgType argtype, NodeArg arg) const
{
    switch (argtype) {
    case ArgType::ConstBool:
        stm << (arg.b ? "true" : "false");
        return;
    case ArgType::ConstInt32:
        stm << arg.i32;
        return;
    case ArgType::ConstFloat64:
        stm << arg.f64;
        return;
    case ArgType::Node:
        stm << "node(" << arg.id << ")";
        return;
    case ArgType::Measure:
        stm << "m(" << arg.id << ")";
        return;
    case ArgType::Global:
        stm << "g(" << arg.id << ")";
        return;
    case ArgType::Arg:
        stm << "arg(" << arg.id << ")";
        return;
    default:
        stm << "<unknown>";
    }
}

static const char *node_op_name(Builder::OpCode opcode)
{
    switch (opcode) {
    case Builder::OpCode::Add:
        return "add";
    case Builder::OpCode::Sub:
        return "sub";
    case Builder::OpCode::Mul:
        return "mul";
    case Builder::OpCode::Div:
        return "div";
    case Builder::OpCode::CmpLT:
        return "cmplt";
    case Builder::OpCode::CmpGT:
        return "cmpgt";
    case Builder::OpCode::CmpLE:
        return "cmple";
    case Builder::OpCode::CmpGE:
        return "cmpge";
    case Builder::OpCode::CmpNE:
        return "cmpne";
    case Builder::OpCode::CmpEQ:
        return "cmpeq";
    case Builder::OpCode::And:
        return "and";
    case Builder::OpCode::Or:
        return "or";
    case Builder::OpCode::Xor:
        return "xor";
    case Builder::OpCode::Not:
        return "not";
    case Builder::OpCode::Abs:
        return "abs";
    case Builder::OpCode::Ceil:
        return "ceil";
    case Builder::OpCode::Exp:
        return "exp";
    case Builder::OpCode::Expm1:
        return "expm1";
    case Builder::OpCode::Floor:
        return "floor";
    case Builder::OpCode::Log:
        return "log";
    case Builder::OpCode::Log1p:
        return "log1p";
    case Builder::OpCode::Log2:
        return "log2";
    case Builder::OpCode::Log10:
        return "log10";
    case Builder::OpCode::Pow:
        return "pow";
    case Builder::OpCode::Sqrt:
        return "sqrt";
    case Builder::OpCode::Asin:
        return "asin";
    case Builder::OpCode::Acos:
        return "acos";
    case Builder::OpCode::Atan:
        return "atan";
    case Builder::OpCode::Atan2:
        return "atan2";
    case Builder::OpCode::Asinh:
        return "asinh";
    case Builder::OpCode::Acosh:
        return "acosh";
    case Builder::OpCode::Atanh:
        return "atanh";
    case Builder::OpCode::Sin:
        return "sin";
    case Builder::OpCode::Cos:
        return "cos";
    case Builder::OpCode::Tan:
        return "tan";
    case Builder::OpCode::Sinh:
        return "sinh";
    case Builder::OpCode::Cosh:
        return "cosh";
    case Builder::OpCode::Tanh:
        return "tanh";
    case Builder::OpCode::Hypot:
        return "hypot";
    case Builder::OpCode::Erf:
        return "erf";
    case Builder::OpCode::Erfc:
        return "erfc";
    case Builder::OpCode::Gamma:
        return "gamma";
    case Builder::OpCode::Lgamma:
        return "lgamma";
    case Builder::OpCode::Rint:
        return "rint";
    case Builder::OpCode::Max:
        return "max";
    case Builder::OpCode::Min:
        return "min";
    case Builder::OpCode::Mod:
        return "mod";
    case Builder::OpCode::Interp:
        return "interp";
    case Builder::OpCode::Select:
        return "select";
    case Builder::OpCode::Identity:
        return "identity";
    }
    return nullptr;
}

void Builder::print_node(std::ostream &stm, const Node &node) const
{
    switch (node.op) {
    case OpCode::Identity:
        print_node_arg(stm, node.argtypes[0], node.args[0]);
        return;
    case OpCode::Add:
        print_node_arg(stm, node.argtypes[0], node.args[0]);
        stm << " + ";
        print_node_arg(stm, node.argtypes[1], node.args[1]);
        return;
    case OpCode::Sub:
        print_node_arg(stm, node.argtypes[0], node.args[0]);
        stm << " - ";
        print_node_arg(stm, node.argtypes[1], node.args[1]);
        return;
    case OpCode::Mul:
        print_node_arg(stm, node.argtypes[0], node.args[0]);
        stm << " * ";
        print_node_arg(stm, node.argtypes[1], node.args[1]);
        return;
    case OpCode::Div:
        print_node_arg(stm, node.argtypes[0], node.args[0]);
        stm << " / ";
        print_node_arg(stm, node.argtypes[1], node.args[1]);
        return;
    case OpCode::CmpLT:
        print_node_arg(stm, node.argtypes[0], node.args[0]);
        stm << " < ";
        print_node_arg(stm, node.argtypes[1], node.args[1]);
        return;
    case OpCode::CmpGT:
        print_node_arg(stm, node.argtypes[0], node.args[0]);
        stm << " > ";
        print_node_arg(stm, node.argtypes[1], node.args[1]);
        return;
    case OpCode::CmpLE:
        print_node_arg(stm, node.argtypes[0], node.args[0]);
        stm << " <= ";
        print_node_arg(stm, node.argtypes[1], node.args[1]);
        return;
    case OpCode::CmpGE:
        print_node_arg(stm, node.argtypes[0], node.args[0]);
        stm << " >= ";
        print_node_arg(stm, node.argtypes[1], node.args[1]);
        return;
    case OpCode::CmpNE:
        print_node_arg(stm, node.argtypes[0], node.args[0]);
        stm << " != ";
        print_node_arg(stm, node.argtypes[1], node.args[1]);
        return;
    case OpCode::CmpEQ:
        print_node_arg(stm, node.argtypes[0], node.args[0]);
        stm << " == ";
        print_node_arg(stm, node.argtypes[1], node.args[1]);
        return;
    case OpCode::And:
        print_node_arg(stm, node.argtypes[0], node.args[0]);
        stm << " & ";
        print_node_arg(stm, node.argtypes[1], node.args[1]);
        return;
    case OpCode::Or:
        print_node_arg(stm, node.argtypes[0], node.args[0]);
        stm << " | ";
        print_node_arg(stm, node.argtypes[1], node.args[1]);
        return;
    case OpCode::Xor:
        print_node_arg(stm, node.argtypes[0], node.args[0]);
        stm << " ^ ";
        print_node_arg(stm, node.argtypes[1], node.args[1]);
        return;
    case OpCode::Not:
        stm << "! ";
        print_node_arg(stm, node.argtypes[0], node.args[0]);
        return;
    case OpCode::Interp:
        stm << "interp(";
        print_node_arg(stm, node.argtypes[0], node.args[0]);
        stm << ", ";
        print_node_arg(stm, node.argtypes[1], node.args[1]);
        stm << ", ";
        print_node_arg(stm, node.argtypes[2], node.args[2]);
        stm << ", " << "data(" << node.data_id << "))";
        return;
    default:
        break;
    }
    switch (node_narg(node.op)) {
    case 1:
        stm << node_op_name(node.op) << "(";
        print_node_arg(stm, node.argtypes[0], node.args[0]);
        stm << ")";
        return;
    case 2:
        stm << node_op_name(node.op) << "(";
        print_node_arg(stm, node.argtypes[0], node.args[0]);
        stm << ", ";
        print_node_arg(stm, node.argtypes[1], node.args[1]);
        stm << ")";
        return;
    case 3:
        stm << node_op_name(node.op) << "(";
        print_node_arg(stm, node.argtypes[0], node.args[0]);
        stm << ", ";
        print_node_arg(stm, node.argtypes[1], node.args[1]);
        stm << ", ";
        print_node_arg(stm, node.argtypes[2], node.args[2]);
        stm << ")";
        return;
    default:
        stm << "<unknown node type>";
    }
}

static void print_node_id(std::ostream &stm, uint32_t id)
{
    if (id) {
        stm << "node(" << id << ")";
    }
    else {
        stm << "null";
    }
}

static void print_time_id(std::ostream &stm, uint32_t id)
{
    if (id) {
        stm << "t(" << id << ")";
    }
    else {
        stm << "null";
    }
}

static void print_bseq_id(std::ostream &stm, uint32_t id)
{
    if (id) {
        stm << "bseq(" << id << ")";
    }
    else {
        stm << "end";
    }
}

NACS_EXPORT() void Builder::print(std::ostream &stm) const
{
    auto nnodes = nodes.size();
    stm << "Nodes[" << nnodes << "]:" << std::endl;
    for (size_t i = 0; i < nnodes; i++) {
        stm << "  node(" << i + 1 << "): ";
        print_node(stm, nodes[i]);
        stm << std::endl;
    }

    auto nchns = chnnames.size();
    stm << "Channels[" << nchns << "]:" << std::endl;
    for (size_t i = 0; i < nchns; i++)
        stm << "  " << i + 1 << ": " << chnnames[i] << std::endl;

    stm << "Default values[" << defvals.size() << "]:" << std::endl;
    for (auto [id, val]: defvals)
        stm << "  " << id << ": " << val << std::endl;

    auto nslots = slots.size();
    stm << "Global variables[" << nslots << "]:" << std::endl;
    for (size_t i = 0; i < nslots; i++)
        stm << "  g(" << i << "): " << slots[i] << std::endl;

    auto nnoramp = noramp_chns.size();
    if (nnoramp > 0) {
        stm << "No ramp channels[" << nnoramp << "]:" << std::endl;
        stm << "  " << noramp_chns[0];
        for (size_t i = 1; i < nnoramp; i++)
            stm << ", " << noramp_chns[i];
        stm << std::endl;
    }

    auto nbseqs = seqs.size();
    stm << "Basic sequences[" << nbseqs << "]:" << std::endl;
    for (size_t i = 0; i < nbseqs; i++) {
        stm << "  bseq(" << i + 1 << "):" << std::endl;
        auto &seq = seqs[i];

        auto ntimes = seq.times.size();
        stm << "    Times[" << ntimes << "]:" << std::endl;
        for (size_t i = 0; i < ntimes; i++) {
            auto time = seq.times[i];
            stm << "      t(" << i + 1 << "): sign=";
            if (time.sign == Sign::Pos) {
                stm << "pos";
            }
            else if (time.sign == Sign::NonNeg) {
                stm << "nonneg";
            }
            else {
                stm << "unknown";
            }
            stm << ", id=" << time.id << ", delta=node(" << time.delta_node << "), prev=";
            print_time_id(stm, time.prev_id);
            stm << std::endl;
        }

        auto nendtimes = seq.endtimes.size();
        if (nendtimes > 0) {
            stm << "    End times[" << nendtimes << "]:" << std::endl;
            stm << "      t(" << seq.endtimes[0] << ")";
            for (size_t i = 1; i < nendtimes; i++)
                stm << ", t(" << seq.endtimes[i] << ")";
            stm << std::endl;
        }

        auto ntimeorders = seq.timeorders.size();
        if (ntimeorders > 0) {
            stm << "    Time orders[" << ntimeorders << "]:" << std::endl;
            for (auto to: seq.timeorders) {
                const char *sign = "?";
                if (to.sign == Sign::Pos) {
                    sign = "<";
                }
                else if (to.sign == Sign::NonNeg) {
                    sign = "<=";
                }
                stm << "      t(" << to.before_time_id << ") " << sign
                    << " t(" << to.after_time_id << ")" << std::endl;
            }
        }

        stm << "    Outputs[" << seq.outputs.size() << "]:" << std::endl;
        for (auto &[id, output]: seq.outputs) {
            stm << "      " << id << ": time=t(" << output.time_id
                << "), len=";
            print_node_id(stm, output.len_node);
            stm << ", val=" << output.val_node << ", cond=";
            print_node_id(stm, output.cond_node);
            stm << ", chn=" << output.chn << std::endl;
        }

        stm << "    Measures[" << seq.measures.size() << "]:" << std::endl;
        for (auto &[id, measure]: seq.measures)
            stm << "      " << id << ": time=t(" << measure.time_id
                << "), chn=" << measure.chn << std::endl;

        stm << "    Assignments[" << seq.assignments.size() << "]:" << std::endl;
        for (auto &assign: seq.assignments) {
            stm << "      " << assign.assign_id << ": g(" << assign.global_id
                << ") <- node(" << assign.val_node << ")" << std::endl;
        }

        auto nbranches = seq.branches.size();
        stm << "    Branches[" << nbranches << "]:" << std::endl;
        for (size_t i = 0; i < nbranches; i++) {
            auto branch = seq.branches[i];
            stm << "      " << i << ": id=" << branch.branch_id << ", target=";
            print_bseq_id(stm, branch.target_id);
            stm << ", cond=node(" << branch.cond_node << ")" << std::endl;
        }

        stm << "    Default target: ";
        print_bseq_id(stm, seq.default_target);
        stm << std::endl;
    }

    auto ndata = datas.size();
    stm << "Datas[" << ndata << "]:" << std::endl;
    for (size_t i = 0; i < ndata; i++) {
        auto &data = datas[i];
        auto data_sz = data.size();
        stm << "  data(" << i << "): [" << data_sz << "] {";
        if (data_sz) {
            stm << data[0];
            for (size_t i = 1; i < data_sz; i++) {
                stm << ", " << data[i];
            }
        }
        stm << "}" << std::endl;
    }

    stm << "Backend datas[" << backend_datas.size() << "]:" << std::endl;
    for (auto &[name, data]: backend_datas) {
        auto data_sz = data.size();
        stm << "  " << name << ": [" << data_sz << "] {";
        if (data_sz) {
            stm << int(data[0]);
            for (size_t i = 1; i < data_sz; i++) {
                stm << ", " << int(data[i]);
            }
        }
        stm << "}" << std::endl;
    }
}

}
