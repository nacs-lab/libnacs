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

#include "ir.h"

#include <math.h>
#include <iostream>

namespace NaCs {
namespace IR {

typedef double (*fptr_t)(double);

static fptr_t getBuiltinPtr(Builtins id)
{
    switch (id) {
    // f(f)
    case Builtins::acos:
        return ::acos;
    case Builtins::acosh:
        return ::acosh;
    case Builtins::asin:
        return ::asin;
    case Builtins::asinh:
        return ::asinh;
    case Builtins::atan:
        return ::atan;
    case Builtins::atanh:
        return ::atanh;
    case Builtins::cbrt:
        return ::cbrt;
    case Builtins::ceil:
        return ::ceil;
    case Builtins::cos:
        return ::cos;
    case Builtins::cosh:
        return ::cosh;
    case Builtins::erf:
        return ::erf;
    case Builtins::erfc:
        return ::erfc;
    case Builtins::exp:
        return ::exp;
    case Builtins::exp10:
        return ::exp10;
    case Builtins::exp2:
        return ::exp2;
    case Builtins::expm1:
        return ::expm1;
    case Builtins::abs:
        return ::abs;
    case Builtins::floor:
        return ::floor;
    case Builtins::gamma: // tgamma
        return ::tgamma;
    case Builtins::j0:
        return ::j0;
    case Builtins::j1:
        return ::j1;
    case Builtins::lgamma:
        return ::lgamma;
    case Builtins::log:
        return ::log;
    case Builtins::log10:
        return ::log10;
    case Builtins::log1p:
        return ::log1p;
    case Builtins::log2:
        return ::log2;
    case Builtins::pow10:
        return ::pow10;
    case Builtins::rint:
        return ::rint;
    case Builtins::round:
        return ::round;
    case Builtins::sin:
        return ::sin;
    case Builtins::sinh:
        return ::sinh;
    case Builtins::sqrt:
        return ::sqrt;
    case Builtins::tan:
        return ::tan;
    case Builtins::tanh:
        return ::tanh;
    case Builtins::y0:
        return ::y0;
    case Builtins::y1:
        return ::y1;

    // f(f, f)
    case Builtins::atan2:
        return fptr_t(static_cast<double(*)(double, double)>(::atan2));
    case Builtins::copysign:
        return fptr_t(static_cast<double(*)(double, double)>(::copysign));
    case Builtins::fdim:
        return fptr_t(static_cast<double(*)(double, double)>(::fdim));
    case Builtins::max: // fmax
        return fptr_t(static_cast<double(*)(double, double)>(::fmax));
    case Builtins::min: // fmax
        return fptr_t(static_cast<double(*)(double, double)>(::fmin));
    case Builtins::mod: // fmod
        return fptr_t(static_cast<double(*)(double, double)>(::fmod));
    case Builtins::hypot:
        return fptr_t(static_cast<double(*)(double, double)>(::hypot));
    case Builtins::pow:
        return fptr_t(static_cast<double(*)(double, double)>(::pow));
    case Builtins::remainder:
        return fptr_t(static_cast<double(*)(double, double)>(::remainder));

    // f(f, f, f)
    case Builtins::fma:
        return fptr_t(static_cast<double(*)(double, double, double)>(::fma));

    // f(f, i)
    case Builtins::ldexp:
        return fptr_t(static_cast<double(*)(double, int)>(::ldexp));

    // f(i, f)
    case Builtins::jn:
        return fptr_t(static_cast<double(*)(int, double)>(::jn));
    case Builtins::yn:
        return fptr_t(static_cast<double(*)(int, double)>(::yn));
    default:
        return nullptr;
    }
}

enum class BuiltinType : uint8_t {
    Invalid,
    F64_F64,
    F64_F64F64,
    F64_F64F64F64,
    F64_F64I32,
    F64_I32F64,
    };

static BuiltinType getBuiltinType(Builtins id)
{
    switch (id) {
    // f(f)
    case Builtins::acos:
    case Builtins::acosh:
    case Builtins::asin:
    case Builtins::asinh:
    case Builtins::atan:
    case Builtins::atanh:
    case Builtins::cbrt:
    case Builtins::ceil:
    case Builtins::cos:
    case Builtins::cosh:
    case Builtins::erf:
    case Builtins::erfc:
    case Builtins::exp:
    case Builtins::exp10:
    case Builtins::exp2:
    case Builtins::expm1:
    case Builtins::abs:
    case Builtins::floor:
    case Builtins::gamma: // tgamma
    case Builtins::j0:
    case Builtins::j1:
    case Builtins::lgamma:
    case Builtins::log:
    case Builtins::log10:
    case Builtins::log1p:
    case Builtins::log2:
    case Builtins::pow10:
    case Builtins::rint:
    case Builtins::round:
    case Builtins::sin:
    case Builtins::sinh:
    case Builtins::sqrt:
    case Builtins::tan:
    case Builtins::tanh:
    case Builtins::y0:
    case Builtins::y1:
        return BuiltinType::F64_F64;

    // f(f, f)
    case Builtins::atan2:
    case Builtins::copysign:
    case Builtins::fdim:
    case Builtins::max: // fmax
    case Builtins::min: // fmax
    case Builtins::mod: // fmod
    case Builtins::hypot:
    case Builtins::pow:
    case Builtins::remainder:
        return BuiltinType::F64_F64F64;

    // f(f, f, f)
    case Builtins::fma:
        return BuiltinType::F64_F64F64F64;

    // f(f, i)
    case Builtins::ldexp:
        return BuiltinType::F64_F64I32;

    // f(i, f)
    case Builtins::jn:
    case Builtins::yn:
        return BuiltinType::F64_I32F64;
    default:
        return BuiltinType::Invalid;
    }
}

NACS_EXPORT bool checkBuiltinType(Builtins id, Type *args, size_t narg)
{
    switch (getBuiltinType(id)) {
    case BuiltinType::F64_F64:
        return narg == 1 && args[0] == Type::Float64;
    case BuiltinType::F64_F64F64:
        return narg == 2 && args[0] == Type::Float64 && args[1] == Type::Float64;
    case BuiltinType::F64_F64F64F64:
        return (narg == 3 && args[0] == Type::Float64 &&
                args[1] == Type::Float64 && args[2] == Type::Float64);
    case BuiltinType::F64_F64I32:
        return narg == 2 && args[0] == Type::Float64 && args[1] == Type::Int32;
    case BuiltinType::F64_I32F64:
        return narg == 2 && args[0] == Type::Int32 && args[1] == Type::Float64;
    default:
        return false;
    }
}

NACS_EXPORT double evalBuiltin(Builtins id, uint64_t *args)
{
    auto _fptr = getBuiltinPtr(id);
    switch (getBuiltinType(id)) {
    case BuiltinType::F64_F64:
        return _fptr(*(double*)args);
    case BuiltinType::F64_F64F64: {
        auto fptr = (double(*)(double, double))_fptr;
        auto fargs = (double*)args;
        return fptr(fargs[0], fargs[1]);
    }
    case BuiltinType::F64_F64F64F64: {
        auto fptr = (double(*)(double, double, double))_fptr;
        auto fargs = (double*)args;
        return fptr(fargs[0], fargs[1], fargs[2]);
    }
    case BuiltinType::F64_F64I32: {
        auto fptr = (double(*)(double, int))_fptr;
        return fptr(*(double*)args, *(uint32_t*)(args + 1));
    }
    case BuiltinType::F64_I32F64: {
        auto fptr = (double(*)(int, double))_fptr;
        return fptr(*(uint32_t*)args, *(double*)(args + 1));
    }
    default:
        return 0.0;
    }
}

template<typename T, typename T2>
static inline void writeBuff(uint8_t *ptr, T2 _val)
{
    T val(_val);
    unsigned char *src = (unsigned char*)&val;
    for (int i = 0;i < sizeof(T);i++) {
        ptr[i] = src[i];
    }
}

template<typename T>
static inline T readBuff(uint8_t *ptr)
{
    T val;
    unsigned char *dest = (unsigned char*)&val;
    for (int i = 0;i < sizeof(T);i++)
        dest[i] = ptr[i];
    return val;
}

static inline const char *typeName(Type typ)
{
    switch (typ) {
    case Type::Bool:
        return "Bool";
    case Type::Int32:
        return "Int32";
    case Type::Float64:
        return "Float64";
    default:
        return "Bottom";
    }
}

static inline const char *opName(Opcode op)
{
    switch (op) {
    case Opcode::Ret:
        return "ret";
    case Opcode::Br:
        return "br";
    case Opcode::Add:
        return "add";
    case Opcode::Sub:
        return "sub";
    case Opcode::Mul:
        return "mul";
    case Opcode::IDiv:
        return "idiv";
    case Opcode::FDiv:
        return "fdiv";
    case Opcode::Rem:
        return "rem";
    case Opcode::Cast:
        return "cast";
    case Opcode::Cmp:
        return "cmp";
    case Opcode::Phi:
        return "phi";
    case Opcode::Call:
        return "call";
    default:
        return "unknown";
    }
}

void Function::dumpValName(int32_t id)
{
    if (id >= 0) {
        std::cout << "%" << id;
    } else if (id == Consts::False) {
        std::cout << "false";
    } else if (id == Consts::True) {
        std::cout << "true";
    } else {
        auto &constval = consts[Consts::_Offset - id];
        auto val = constval.second;
        switch (constval.first) {
        case Type::Int32:
            std::cout << val.i32;
            break;
        case Type::Float64:
            std::cout << val.f64;
            break;
        default:
            std::cout << "undef";
        }
    }
}

NACS_EXPORT Type Function::valType(int32_t id)
{
    if (id >= 0) {
        return vals[id];
    } else if (id == Consts::False || id == Consts::True) {
        return Type::Bool;
    } else {
        return consts[Consts::_Offset - id].first;
    }
}

void Function::dumpVal(int32_t id)
{
    std::cout << typeName(valType(id)) << " ";
    dumpValName(id);
}

void Function::dumpBB(BB &bb)
{
    uint8_t *pc = bb.data();
    uint8_t *end = pc + bb.size();
    while (end > pc) {
        auto op = Opcode(*pc);
        pc++;
        switch (op) {
        case Opcode::Ret:
            std::cout << "  ret ";
            dumpVal(readBuff<int32_t>(pc));
            std::cout << std::endl;
            pc += 4;
            break;
        case Opcode::Br: {
            auto cond = readBuff<int32_t>(pc);
            pc += 4;
            auto bb1 = readBuff<int32_t>(pc);
            pc += 4;
            if (cond == Consts::True) {
                std::cout << "  br L" << bb1;
            } else {
                auto bb2 = readBuff<int32_t>(pc);
                pc += 4;
                std::cout << "  br ";
                dumpVal(cond);
                std::cout << ", L" << bb1 << ", L" << bb2;
            }
            std::cout << std::endl;
            break;
        }
        case Opcode::Add:
        case Opcode::Sub:
        case Opcode::Mul: {
            auto res = readBuff<int32_t>(pc);
            pc += 4;
            auto val1 = readBuff<int32_t>(pc);
            pc += 4;
            auto val2 = readBuff<int32_t>(pc);
            pc += 4;
            std::cout << "  ";
            dumpVal(res);
            std::cout << " = " << opName(op) << " ";
            dumpVal(val1);
            std::cout << ", ";
            dumpVal(val2);
            std::cout << std::endl;
            break;
        }
        default:
            std::cout << "  unknown op: " << uint8_t(op) << std::endl;
            break;
        }
    }
}

NACS_EXPORT void Function::dump(void)
{
    std::cout << typeName(ret) << " (";
    for (int i = 0;i < nargs;i++) {
        if (i != 0)
            std::cout << ", ";
        std::cout << typeName(valType(i)) << " ";
        dumpValName(i);
    }
    std::cout << ") {" << std::endl;
    for (size_t i = 0;i < code.size();i++) {
        std::cout << "L" << i << ":" << std::endl;
        dumpBB(code[i]);
    }
    std::cout << "}" << std::endl;
}

uint8_t *Builder::addInst(Opcode op, size_t nbytes)
{
    auto &bb = m_f.code[m_cur_bb];
    auto oldlen = bb.size();
    bb.resize(oldlen + nbytes + 1);
    bb[oldlen] = uint8_t(op);
    return &bb[oldlen + 1];
}

void Builder::createRet(int32_t val)
{
    writeBuff<uint32_t>(addInst(Opcode::Ret, 4), val);
}

int32_t Builder::getConstInt(int32_t val)
{
    auto map = m_f.const_ints;
    auto it = map.find(val);
    if (it != map.end())
        return it->second;
    int32_t oldlen = m_f.consts.size();
    Function::ConstVal constval;
    constval.i32 = val;
    m_f.consts.push_back(std::make_pair(Type::Int32, constval));
    int32_t id = Consts::_Offset - oldlen;
    map[val] = id;
    return id;
}

int32_t Builder::getConstFloat(double val)
{
    auto map = m_f.const_floats;
    auto it = map.find(val);
    if (it != map.end())
        return it->second;
    int32_t oldlen = m_f.consts.size();
    Function::ConstVal constval;
    constval.f64 = val;
    m_f.consts.push_back(std::make_pair(Type::Float64, constval));
    int32_t id = Consts::_Offset - oldlen;
    map[val] = id;
    return id;
}

int32_t Builder::newSSA(Type typ)
{
    int32_t id = m_f.vals.size();
    m_f.vals.push_back(typ);
    return id;
}

int32_t Builder::newBB(void)
{
    int32_t id = m_f.code.size();
    m_f.code.push_back({});
    return id;
}

int32_t &Builder::curBB()
{
    return m_cur_bb;
}

void Builder::createBr(int32_t br)
{
    createBr(Consts::True, br, 0);
}

void Builder::createBr(int32_t cond, int32_t bb1, int32_t bb2)
{
    if (cond == Consts::True) {
        uint8_t *ptr = addInst(Opcode::Br, 8);
        writeBuff<uint32_t>(ptr, Consts::True);
        writeBuff<uint32_t>(ptr + 4, bb1);
    }
    else {
        uint8_t *ptr = addInst(Opcode::Br, 12);
        writeBuff<uint32_t>(ptr, cond);
        writeBuff<uint32_t>(ptr + 4, bb1);
        writeBuff<uint32_t>(ptr + 8, bb2);
    }
}

int32_t Builder::createPromoteOP(Opcode op, int32_t val1, int32_t val2)
{
    // TODO optimize for constants
    uint8_t *ptr = addInst(op, 12);
    auto ty1 = m_f.valType(val1);
    auto ty2 = m_f.valType(val2);
    auto resty = std::max(std::max(ty1, ty2), Type::Int32);
    auto res = newSSA(resty);
    writeBuff<uint32_t>(ptr, res);
    writeBuff<uint32_t>(ptr + 4, val1);
    writeBuff<uint32_t>(ptr + 8, val2);
    return res;
}

int32_t Builder::createAdd(int32_t val1, int32_t val2)
{
    return createPromoteOP(Opcode::Add, val1, val2);
}

int32_t Builder::createSub(int32_t val1, int32_t val2)
{
    return createPromoteOP(Opcode::Sub, val1, val2);
}

int32_t Builder::createMul(int32_t val1, int32_t val2)
{
    return createPromoteOP(Opcode::Mul, val1, val2);
}

}
}
