/*************************************************************************
 *   Copyright (c) 2016 - 2016 Yichao Yu <yyc1992@gmail.com>             *
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
#include <vector>
#include <map>

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
    IDiv,
    FDiv,
    Rem,
    Cast,
    Cmp,
    Phi,
    Call,
    _Max,
    };

constexpr static inline bool validate(Opcode op)
{
    return op < Opcode::_Max && op > Opcode::_Min;
}

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

NACS_EXPORT bool checkBuiltinType(Builtins id, Type *args, size_t narg);
NACS_EXPORT double evalBuiltin(Builtins id, uint64_t *args);

enum Consts : int32_t {
    False = -1,
    True = -2,
    _Offset = -3,
};

struct Function {
    // TODO const pool
    typedef std::vector<uint8_t> BB;
    typedef union {
        int32_t i32;
        double f64;
    } ConstVal;
    Function(Type _ret, const std::vector<Type> &args)
        : ret(_ret),
          nargs(args.size()),
          vals(args),
          code{BB{}},
          consts{},
          const_ints{},
          const_floats{}
    {}
    void dump(void);
    Type valType(int32_t id);
    const Type ret;
    const int nargs;
    std::vector<Type> vals;
    std::vector<BB> code;
    std::vector<std::pair<Type, ConstVal>> consts;
    std::map<int32_t, int> const_ints;
    std::map<double, int> const_floats;
private:
    void dumpValName(int32_t id);
    void dumpVal(int32_t id);
    void dumpBB(BB&);
};

class NACS_EXPORT Builder {
public:
    Builder(Type ret, std::vector<Type> args)
        : m_f(ret, args),
          m_cur_bb(0)
    {}
    Function &get(void)
    {
        return m_f;
    }
    void createRet(int32_t val);
    int32_t getConstInt(int32_t val);
    int32_t getConstFloat(double val);
    // Br,
    // Add,
    // Sub,
    // Mul,
    // IDiv,
    // FDiv,
    // Rem,
    // Cast,
    // Cmp,
    // Phi,
    // Call,
private:
    uint8_t *addInst(Opcode op, size_t nbytes);
    Function m_f;
    int m_cur_bb;
};

}
}

#endif
