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

#ifndef __NACS_SEQ_BUILDER_H__
#define __NACS_SEQ_BUILDER_H__

#include "sign.h"
#include "env.h"

#include <nacs-utils/ir.h>

#include <list>
#include <map>
#include <ostream>
#include <set>

#include <llvm/ADT/SmallVector.h>

namespace NaCs::Seq {

class EventTime;
class BasicSeq;
class Seq;

class Builder {
    enum class NodeFlag : uint8_t {
        None,
        Scanning,
        RefArg,
        Free
    };
    template<typename T> using StackVector = llvm::SmallVector<T,16>;
    struct NodeFieldIterator;
public:
    enum class OpCode : uint8_t {
        Add = 1,
        Sub = 2,
        Mul = 3,
        Div = 4,
        CmpLT = 5,
        CmpGT = 6,
        CmpLE = 7,
        CmpGE = 8,
        CmpNE = 9,
        CmpEQ = 10,
        And = 11,
        Or = 12,
        Xor = 13,
        Not = 14,
        Abs = 15,
        Ceil = 16,
        Exp = 17,
        Expm1 = 18,
        Floor = 19,
        Log = 20,
        Log1p = 21,
        Log2 = 22,
        Log10 = 23,
        Pow = 24,
        Sqrt = 25,
        Asin = 26,
        Acos = 27,
        Atan = 28,
        Atan2 = 29,
        Asinh = 30,
        Acosh = 31,
        Atanh = 32,
        Sin = 33,
        Cos = 34,
        Tan = 35,
        Sinh = 36,
        Cosh = 37,
        Tanh = 38,
        Hypot = 39,
        Erf = 40,
        Erfc = 41,
        Gamma = 42,
        Lgamma = 43,
        Rint = 44,
        Max = 45,
        Min = 46,
        Mod = 47,
        Interp = 48,
        Select = 49,
        Identity = 50,
        _MaxOP = 50,
    };
    enum class ArgType : uint8_t {
        ConstBool = 0,
        ConstInt32 = 1,
        ConstFloat64 = 2,
        Node = 3,
        Measure = 4,
        Global = 5,
        Arg = 6,
        _MaxType = 6,
    };
    struct Node;
    union NodeArg {
        bool b;
        int32_t i32;
        double f64;
        uint32_t id;
    };
    struct Time {
        Sign sign;
        uint32_t id;
        uint32_t delta_node;
        uint32_t prev_id;
    };
    struct Output {
        // uint32_t id;
        uint32_t time_id;
        uint32_t len_node;
        uint32_t val_node;
        uint32_t cond_node;
        uint32_t chn;
    };
    struct Measure {
        // uint32_t id;
        uint32_t time_id;
        uint32_t chn;
    };
    struct Assignment {
        uint32_t assign_id;
        uint32_t global_id;
        uint32_t val_node;
    };
    struct Branch {
        uint32_t branch_id;
        uint32_t target_id;
        uint32_t cond_node;
    };
    struct TimeOrder {
        Sign sign;
        uint32_t id; // Unused for now
        uint32_t before_time_id;
        uint32_t after_time_id;
    };
    struct BasicSeq;

    struct NodeVisitor;

    Builder(Seq &seq);
    ~Builder();

    static int8_t node_narg(OpCode op);

    inline Node *get_argref(ArgType argtype, NodeArg arg)
    {
        return argtype == ArgType::Node ? get_node(arg.id) : nullptr;
    }

    void buildseq();
    void print(std::ostream &stm) const;

    std::vector<Node> nodes;
    std::vector<BasicSeq> seqs;
    std::vector<std::string> chnnames;
    std::map<uint32_t,IR::TagVal> defvals;
    std::vector<IR::Type> slots;
    std::vector<std::vector<double>> datas;
    std::vector<uint32_t> noramp_chns;
    std::map<std::string,std::vector<uint8_t>> backend_datas;

private:
    // Verify sequence
    // Set default values
    // Create slots
    // Create channels
    // Create basic seqs
    // Create measure vars
    void prescan_seq();
    // Verify and create variables from nodes.
    // `m_measure_vars` must be populated first.
    void scan_node();
    // Verify type of node
    // Create assumptions
    // Create times
    // Create pulses and measures
    void scan_seq();

    EventTime *get_eventtime(uint32_t time_id, BasicSeq &seq);
    Var *arg_to_var(ArgType argtype, NodeArg arg);
    Arg arg_to_callarg(ArgType argtype, NodeArg arg);
    int32_t create_inst(IR::Builder &builder, Node *node, int32_t *args);
    std::pair<Var*,bool> node_to_rampvar(Node *node);
    Node *get_node(uint32_t id);
    const Node *get_node(uint32_t id) const;
    uint32_t get_node_id(const Node *node) const;

    void print_node_arg(std::ostream &stm, ArgType argtype, NodeArg arg) const;
    void print_node(std::ostream &stm, const Node &node) const;

    Seq &m_seq;
    std::vector<bool> m_mayramp;
    std::vector<Var*> m_varmap;
    std::map<uint32_t,Var*> m_measure_vars;
    std::set<uint32_t> m_disabled_measures;
};

struct Builder::Node {
    NodeArg args[3];
    ArgType argtypes[3];
    OpCode op;
    uint32_t data_id = 0;
private:
    NodeFlag flag = NodeFlag::None;
    friend class Builder;
};

struct Builder::NodeFieldIterator : FieldIterator<Node*> {
    NodeFieldIterator(Node *node, Builder &builder)
        : builder(builder),
          m_node(node),
          m_narg(node_narg(node->op))
    {
    }
    NodeFieldIterator(const NodeFieldIterator &other)
        : builder(other.builder),
          m_node(other.m_node),
          m_narg(other.m_narg),
          m_idx(other.m_idx)
    {
    }
    Node *get() const
    {
        return builder.get_argref(m_node->argtypes[m_idx], m_node->args[m_idx]);
    }
    Node *parent() const
    {
        return m_node;
    }
    bool is_end() const
    {
        return m_idx >= m_narg;
    }
    NodeFieldIterator &operator++()
    {
        m_idx++;
        return *this;
    }
    NodeFieldIterator &operator=(const NodeFieldIterator &other)
    {
        assert(&builder == &other.builder);
        m_node = other.m_node;
        m_narg = other.m_narg;
        m_idx = other.m_idx;
        return *this;
    }

private:
    Builder &builder;
    Node *m_node;
    int8_t m_narg;
    int8_t m_idx = 0;
};

struct Builder::BasicSeq {
    std::map<uint32_t,Output> outputs;
    std::map<uint32_t,Measure> measures;
    std::vector<Assignment> assignments;
    std::vector<Time> times;
    std::vector<TimeOrder> timeorders;
    std::vector<uint32_t> endtimes;
    std::vector<Branch> branches;
    uint32_t default_target;
private:
    Time *get_time(uint32_t id)
    {
        if (id == 0 || id > times.size())
            throw std::runtime_error("Out of bound time ID");
        return &times[id - 1];
    }
    const Time *get_time(uint32_t id) const
    {
        if (id == 0 || id > times.size())
            throw std::runtime_error("Out of bound time ID");
        return &times[id - 1];
    }
    uint32_t get_time_id(const Time *time) const
    {
        return uint32_t(time - &times[0]) + 1;
    }

    ::NaCs::Seq::BasicSeq* m_basicseq = nullptr;
    std::vector<EventTime*> m_eventtimes;
    friend class Builder;
};

struct Builder::NodeVisitor: DFSVisitor<Node*,NodeFieldIterator,StackVector> {
    NodeVisitor(Builder &builder)
        : builder(builder)
    {
    }
    NodeFieldIterator begin(Node *node)
    {
        return NodeFieldIterator(node, builder);
    }

    Builder &builder;
};

inline int8_t Builder::node_narg(OpCode op)
{
    switch (op) {
    case OpCode::Add: return 2;
    case OpCode::Sub: return 2;
    case OpCode::Mul: return 2;
    case OpCode::Div: return 2;
    case OpCode::CmpLT: return 2;
    case OpCode::CmpGT: return 2;
    case OpCode::CmpLE: return 2;
    case OpCode::CmpGE: return 2;
    case OpCode::CmpNE: return 2;
    case OpCode::CmpEQ: return 2;
    case OpCode::And: return 2;
    case OpCode::Or: return 2;
    case OpCode::Xor: return 2;
    case OpCode::Not: return 1;
    case OpCode::Abs: return 1;
    case OpCode::Ceil: return 1;
    case OpCode::Exp: return 1;
    case OpCode::Expm1: return 1;
    case OpCode::Floor: return 1;
    case OpCode::Log: return 1;
    case OpCode::Log1p: return 1;
    case OpCode::Log2: return 1;
    case OpCode::Log10: return 1;
    case OpCode::Pow: return 2;
    case OpCode::Sqrt: return 1;
    case OpCode::Asin: return 1;
    case OpCode::Acos: return 1;
    case OpCode::Atan: return 1;
    case OpCode::Atan2: return 2;
    case OpCode::Asinh: return 1;
    case OpCode::Acosh: return 1;
    case OpCode::Atanh: return 1;
    case OpCode::Sin: return 1;
    case OpCode::Cos: return 1;
    case OpCode::Tan: return 1;
    case OpCode::Sinh: return 1;
    case OpCode::Cosh: return 1;
    case OpCode::Tanh: return 1;
    case OpCode::Hypot: return 2;
    case OpCode::Erf: return 1;
    case OpCode::Erfc: return 1;
    case OpCode::Gamma: return 1;
    case OpCode::Lgamma: return 1;
    case OpCode::Rint: return 1;
    case OpCode::Max: return 2;
    case OpCode::Min: return 2;
    case OpCode::Mod: return 2;
    case OpCode::Interp: return 3;
    case OpCode::Select: return 3;
    case OpCode::Identity: return 1;
    default: return 0;
    }
}

inline Builder::Node *Builder::get_node(uint32_t id)
{
    if (id == 0 || id > nodes.size())
        throw std::runtime_error("Out of bound node ID");
    return &nodes[id - 1];
}

inline const Builder::Node *Builder::get_node(uint32_t id) const
{
    if (id == 0 || id > nodes.size())
        throw std::runtime_error("Out of bound node ID");
    return &nodes[id - 1];
}

inline uint32_t Builder::get_node_id(const Node *node) const
{
    return uint32_t(node - &nodes[0]) + 1;
}

static inline std::ostream &operator<<(std::ostream &stm, const Builder &builder)
{
    builder.print(stm);
    return stm;
}

}

#endif
