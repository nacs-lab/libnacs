/*************************************************************************
 *   Copyright (c) 2016 - 2018 Yichao Yu <yyc1992@gmail.com>             *
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

#include "utils.h"

#include <assert.h>
#include <vector>
#include <map>
#include <memory>
#include <iostream>

#ifndef __NACS_UTILS_IR_H__
#define __NACS_UTILS_IR_H__

namespace NaCs {
namespace IR {

enum class Type : uint8_t {
    _Min = 0,
    Bool,
    Int32,
    Float64,
    _Max,
    };

constexpr static inline bool validate(Type typ)
{
    return typ < Type::_Max && typ > Type::_Min;
}

enum class Opcode : uint8_t {
    _Min = 0,
    Ret,
    Br,
    Add,
    Sub,
    Mul,
    FDiv,
    Cmp,
    Phi,
    Call,
    Interp,
    _Max,
    };

constexpr static inline bool validate(Opcode op)
{
    return op < Opcode::_Max && op > Opcode::_Min;
}

enum class CmpType : uint8_t {
    eq,
    gt,
    ge,
    lt,
    le,
    ne,
    };

enum class Builtins : uint16_t {
    // f(f)
    acos,
    acosh,
    asin,
    asinh,
    atan,
    atanh,
    cbrt,
    ceil,
    cos,
    cosh,
    erf,
    erfc,
    exp,
    exp10,
    exp2,
    expm1,
    abs,
    floor,
    gamma, // tgamma
    j0,
    j1,
    lgamma,
    log,
    log10,
    log1p,
    log2,
    pow10,
    rint,
    round,
    sin,
    sinh,
    sqrt,
    tan,
    tanh,
    y0,
    y1,

    // f(f, f)
    atan2,
    copysign,
    fdim,
    max, // fmax
    min, // fmax
    mod, // fmod
    hypot,
    pow,
    remainder,

    // f(f, f, f)
    fma,

    // f(f, i)
    ldexp,

    // f(i, f)
    jn,
    yn,
    };

enum Consts : int32_t {
    False = -1,
    True = -2,
    _Offset = -3,
};

typedef union {
    bool b;
    int32_t i32;
    double f64;
} GenVal;

struct TagVal {
    Type typ;
    GenVal val;
    TagVal() = default;
    TagVal(Type _typ, GenVal _val=GenVal{})
        : typ(_typ),
          val(_val)
    {}
    TagVal(const TagVal&) = default;
    TagVal(bool b)
        : typ(Type::Bool)
    {
        val.b = b;
    }
    template<typename T>
    TagVal(T i, std::enable_if_t<std::is_integral<T>::value>* = nullptr)
        : typ(Type::Int32)
    {
        val.i32 = int32_t(i);
    }
    template<typename T>
    TagVal(T f, std::enable_if_t<std::is_floating_point<T>::value>* = nullptr)
        : typ(Type::Float64)
    {
        val.f64 = double(f);
    }
    template<typename T> T get(void) const
    {
        switch (typ) {
        case Type::Bool:
            return T(val.b);
        case Type::Int32:
            return T(val.i32);
        case Type::Float64:
            return T(val.f64);
        default:
            return T(false);
        }
    }
    TagVal convert(Type new_typ) const
    {
        switch (new_typ) {
        case Type::Bool:
            return get<bool>();
        case Type::Int32:
            return get<int32_t>();
        case Type::Float64:
            return get<double>();
        default:
            return TagVal(new_typ);
        }
    }
    void dump(void) const;
} __attribute__((aligned(8)));

std::ostream &operator<<(std::ostream &stm, const TagVal &val);

inline void TagVal::dump(void) const
{
    std::cerr << *this << std::endl;
}

bool checkBuiltinType(Builtins id, Type *args, size_t narg);
double evalBuiltin(Builtins id, TagVal *args);

struct Function {
    typedef std::vector<int32_t> BB;
    typedef std::pair<int32_t, int32_t> InstRef;
    Function(Type _ret, const std::vector<Type> &args)
        : ret(_ret),
          nargs((int)args.size()),
          vals(args),
          code{BB{}},
          consts{},
          float_table{}
    {}
    Function(const uint32_t*, size_t);
    Function(const uint8_t*, size_t);
    Function(const std::vector<uint32_t> &data)
        : Function(data.data(), data.size())
    {}
    Function(Function&&) = default;
    Function(const Function&) = default;
    void dump(void) const;
    Type valType(int32_t id) const;
    TagVal evalConst(int32_t id) const
    {
        assert(id < 0);
        if (id == Consts::False) {
            return false;
        } else if (id == Consts::True) {
            return true;
        } else {
            return consts[Consts::_Offset - id];
        }
    }
    std::vector<uint32_t> serialize(void) const;
    const Type ret;
    const int nargs;
    std::vector<Type> vals;
    std::vector<BB> code;
    std::vector<TagVal> consts;
    std::vector<double> float_table;
private:
    void printValName(std::ostream &stm, int32_t id) const;
    void printVal(std::ostream &stm, int32_t id) const;
    void printBB(std::ostream &stm, const BB&) const;
    friend std::ostream &operator<<(std::ostream &stm, const Function &f);
};

std::ostream &operator<<(std::ostream &stm, const Function &f);

inline void Function::dump(void) const
{
    std::cerr << *this << std::endl;
}

class NACS_EXPORT Builder {
public:
    Builder(Type ret, std::vector<Type> args)
        : m_f(ret, args),
          m_cur_bb(0),
          const_ints{},
          const_floats{}
    {}
    Function &get(void)
    {
        return m_f;
    }
    int32_t getConst(TagVal val);
    int32_t getConstInt(int32_t val);
    int32_t getConstFloat(double val);

    int32_t newBB(void);
    int32_t &curBB(void);

    void createRet(int32_t val);
    void createBr(int32_t br);
    void createBr(int32_t cond, int32_t bb1, int32_t bb2);
    int32_t createAdd(int32_t val1, int32_t val2);
    int32_t createSub(int32_t val1, int32_t val2);
    int32_t createMul(int32_t val1, int32_t val2);
    int32_t createFDiv(int32_t val1, int32_t val2);
    int32_t createCmp(CmpType cmptyp, int32_t val1, int32_t val2);
    std::pair<int32_t, Function::InstRef> createPhi(Type typ, int ninputs);
    void addPhiInput(Function::InstRef phi, int32_t bb, int32_t val);
    int32_t createCall(Builtins id, int32_t nargs, const int32_t *args);
    int32_t createCall(Builtins id, const std::vector<int32_t> &args)
    {
        return createCall(id, (int32_t)args.size(), args.data());
    }
    int32_t createInterp(int32_t v, double x0, double dx, uint32_t npoints,
                         const double *points);
private:
    int32_t addFloatData(const double *data, uint32_t ndata);
    int32_t *addInst(Opcode op, size_t nop);
    int32_t *addInst(Opcode op, size_t nop, Function::InstRef &inst);
    int32_t newSSA(Type typ);
    int32_t createPromoteOP(Opcode op, int32_t val1, int32_t val2);
    Function m_f;
    int32_t m_cur_bb;
    std::map<int32_t, int> const_ints;
    std::map<double, int> const_floats;
};

struct NACS_EXPORT EvalContext {
    EvalContext(const Function &f)
        : m_f(f),
          m_vals(f.vals.size())
    {}
    void reset(const GenVal *args)
    {
        memcpy(m_vals.data(), args, m_f.nargs * sizeof(GenVal));
    }
    void reset(std::vector<TagVal> tagvals)
    {
        std::vector<GenVal> args(m_f.nargs);
        for (int i = 0;i < m_f.nargs;i++)
            args[i] = tagvals[i].convert(m_f.vals[i]).val;
        reset(args.data());
    }
    void reset(int idx, const GenVal &arg)
    {
        m_vals[idx] = arg;
    }
    TagVal evalVal(int32_t id) const;
    TagVal eval(void);

private:
    const Function &m_f;
    std::vector<GenVal> m_vals;
};

class ExeFunc {
    template<typename FT> struct FuncType;
    template<typename R, typename... Args> struct FuncType<R(Args...)> {
        using ret = R;
        using fptr = R(*)(void*, Args...);
    };
public:
    ExeFunc(void (*cb)(void), void *data, void (*free)(void*))
        : m_cb(cb), m_data(data), m_free(free)
    {}
    ExeFunc(ExeFunc &&other)
        : m_cb(other.m_cb),
          m_data(other.m_data),
          m_free(other.m_free)
    {
        other.m_cb = nullptr;
        other.m_data = nullptr;
        other.m_free = nullptr;
    }
    ExeFunc(const ExeFunc&) = delete;
    ExeFunc &operator=(ExeFunc &&other)
    {
        std::swap(other.m_cb, m_cb);
        std::swap(other.m_data, m_data);
        std::swap(other.m_free, m_free);
        return *this;
    }
    ExeFunc &operator=(const ExeFunc&) = delete;
    ~ExeFunc()
    {
        if (m_free) {
            m_free(m_data);
        }
    }
    template<typename FT, typename... Args> typename FuncType<FT>::ret call(Args&&... args)
    {
        return typename FuncType<FT>::fptr(m_cb)(m_data, std::forward<Args>(args)...);
    }

private:
    void (*m_cb)(void);
    void *m_data;
    void (*m_free)(void*);
};

class ExeContext {
public:
    static std::unique_ptr<ExeContext> get();
    virtual ExeFunc getFunc(const Function &) = 0;
    virtual ExeFunc getFunc(Function&&) = 0;
};

}
}

#endif
