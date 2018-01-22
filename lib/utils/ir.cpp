/*************************************************************************
 *   Copyright (c) 2016 - 2017 Yichao Yu <yyc1992@gmail.com>             *
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
#include "number.h"
#include "mem.h"

#include <math.h>
#include <iostream>

extern "C" void nacs_exefunc_cb(void);

namespace NaCs {
namespace IR {

typedef double (*fptr_t)(double);

#ifndef NACS_HAS_EXP10
double exp10(double x)
{
    return pow(10, x);
}
#endif

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
#ifdef NACS_HAS_EXP10
        return ::exp10;
#else
        return exp10;
#endif
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
#ifdef NACS_HAS_EXP10
        return ::exp10;
#else
        return exp10;
#endif
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

static const char *builtinName(Builtins id)
{
    switch (id) {
        // f(f)
    case Builtins::acos:
        return "acos";
    case Builtins::acosh:
        return "acosh";
    case Builtins::asin:
        return "asin";
    case Builtins::asinh:
        return "asinh";
    case Builtins::atan:
        return "atan";
    case Builtins::atanh:
        return "atanh";
    case Builtins::cbrt:
        return "cbrt";
    case Builtins::ceil:
        return "ceil";
    case Builtins::cos:
        return "cos";
    case Builtins::cosh:
        return "cosh";
    case Builtins::erf:
        return "erf";
    case Builtins::erfc:
        return "erfc";
    case Builtins::exp:
        return "exp";
    case Builtins::exp10:
        return "exp10";
    case Builtins::exp2:
        return "exp2";
    case Builtins::expm1:
        return "expm1";
    case Builtins::abs:
        return "abs";
    case Builtins::floor:
        return "floor";
    case Builtins::gamma:
        return "gamma";
    case Builtins::j0:
        return "j0";
    case Builtins::j1:
        return "j1";
    case Builtins::lgamma:
        return "lgamma";
    case Builtins::log:
        return "log";
    case Builtins::log10:
        return "log10";
    case Builtins::log1p:
        return "log1p";
    case Builtins::log2:
        return "log2";
    case Builtins::pow10:
        return "exp10";
    case Builtins::rint:
        return "rint";
    case Builtins::round:
        return "round";
    case Builtins::sin:
        return "sin";
    case Builtins::sinh:
        return "sinh";
    case Builtins::sqrt:
        return "sqrt";
    case Builtins::tan:
        return "tan";
    case Builtins::tanh:
        return "tanh";
    case Builtins::y0:
        return "y0";
    case Builtins::y1:
        return "y1";

        // f(f, f)
    case Builtins::atan2:
        return "atan2";
    case Builtins::copysign:
        return "copysign";
    case Builtins::fdim:
        return "fdim";
    case Builtins::max: // fmax
        return "max";
    case Builtins::min: // fmax
        return "min";
    case Builtins::mod: // fmod
        return "mod";
    case Builtins::hypot:
        return "hypot";
    case Builtins::pow:
        return "pow";
    case Builtins::remainder:
        return "remainder";

        // f(f, f, f)
    case Builtins::fma:
        return "fma";

        // f(f, i)
    case Builtins::ldexp:
        return "ldexp";

        // f(i, f)
    case Builtins::jn:
        return "jn";
    case Builtins::yn:
        return "yn";
    default:
        return "unknown";
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

NACS_EXPORT() bool checkBuiltinType(Builtins id, Type *args, size_t narg)
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

NACS_EXPORT() double evalBuiltin(Builtins id, TagVal *args)
{
    auto _fptr = getBuiltinPtr(id);
    switch (getBuiltinType(id)) {
    case BuiltinType::F64_F64:
        return _fptr(args[0].get<double>());
    case BuiltinType::F64_F64F64: {
        auto fptr = (double(*)(double, double))_fptr;
        return fptr(args[0].get<double>(), args[1].get<double>());
    }
    case BuiltinType::F64_F64F64F64: {
        auto fptr = (double(*)(double, double, double))_fptr;
        return fptr(args[0].get<double>(), args[1].get<double>(),
                    args[2].get<double>());
    }
    case BuiltinType::F64_F64I32: {
        auto fptr = (double(*)(double, int))_fptr;
        return fptr(args[0].get<double>(), args[1].get<int32_t>());
    }
    case BuiltinType::F64_I32F64: {
        auto fptr = (double(*)(int, double))_fptr;
        return fptr(args[0].get<int32_t>(), args[1].get<double>());
    }
    default:
        return 0.0;
    }
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
    case Opcode::FDiv:
        return "fdiv";
    case Opcode::Cmp:
        return "cmp";
    case Opcode::Phi:
        return "phi";
    case Opcode::Call:
        return "call";
    case Opcode::Interp:
        return "interp";
    default:
        return "unknown";
    }
}

static inline const char *cmpName(CmpType cmp)
{
    switch (cmp) {
    case CmpType::eq:
        return "eq";
    case CmpType::gt:
        return "gt";
    case CmpType::ge:
        return "ge";
    case CmpType::lt:
        return "lt";
    case CmpType::le:
        return "le";
    case CmpType::ne:
        return "ne";
    default:
        return "unknown";
    }
}

NACS_EXPORT() std::ostream &operator<<(std::ostream &stm, const TagVal &val)
{
    stm << typeName(val.typ) << " ";
    switch (val.typ) {
    case Type::Bool:
        stm << (val.val.b ? "true" : "false");
        break;
    case Type::Int32:
        stm << val.val.i32;
        break;
    case Type::Float64:
        stm << val.val.f64;
        break;
    default:
        stm << "undef";
    }
    return stm;
}

void Function::printValName(std::ostream &stm, int32_t id) const
{
    if (id >= 0) {
        stm << "%" << id;
    } else if (id == Consts::False) {
        stm << "false";
    } else if (id == Consts::True) {
        stm << "true";
    } else {
        auto &constval = consts[Consts::_Offset - id];
        auto val = constval.val;
        switch (constval.typ) {
        case Type::Int32:
            stm << val.i32;
            break;
        case Type::Float64:
            stm << val.f64;
            break;
        default:
            stm << "undef";
        }
    }
}

NACS_EXPORT() Type Function::valType(int32_t id) const
{
    if (id >= 0) {
        return vals[id];
    } else if (id == Consts::False || id == Consts::True) {
        return Type::Bool;
    } else {
        return consts[Consts::_Offset - id].typ;
    }
}

void Function::printVal(std::ostream &stm, int32_t id) const
{
    stm << typeName(valType(id)) << " ";
    printValName(stm, id);
}

void Function::printBB(std::ostream &stm, const BB &bb) const
{
    const int32_t *pc = bb.data();
    const int32_t *end = pc + bb.size();
    while (end > pc) {
        auto op = Opcode(*pc);
        pc++;
        switch (op) {
        case Opcode::Ret:
            stm << "  ret ";
            printVal(stm, *pc);
            stm << std::endl;
            pc++;
            break;
        case Opcode::Br: {
            auto cond = *pc;
            pc++;
            auto bb1 = *pc;
            pc++;
            if (cond == Consts::True) {
                stm << "  br L" << bb1;
            } else {
                auto bb2 = *pc;
                pc++;
                stm << "  br ";
                printVal(stm, cond);
                stm << ", L" << bb1 << ", L" << bb2;
            }
            stm << std::endl;
            break;
        }
        case Opcode::Add:
        case Opcode::Sub:
        case Opcode::Mul:
        case Opcode::FDiv: {
            auto res = *pc;
            pc++;
            auto val1 = *pc;
            pc++;
            auto val2 = *pc;
            pc++;
            stm << "  ";
            printVal(stm, res);
            stm << " = " << opName(op) << " ";
            printVal(stm, val1);
            stm << ", ";
            printVal(stm, val2);
            stm << std::endl;
            break;
        }
        case Opcode::Cmp: {
            auto res = *pc;
            pc++;
            auto cmptyp = CmpType(*pc);
            pc++;
            auto val1 = *pc;
            pc++;
            auto val2 = *pc;
            pc++;
            stm << "  ";
            printVal(stm, res);
            stm << " = " << opName(op) << " " << cmpName(cmptyp) << " ";
            printVal(stm, val1);
            stm << ", ";
            printVal(stm, val2);
            stm << std::endl;
            break;
        }
        case Opcode::Phi: {
            auto res = *pc;
            pc++;
            auto nargs = *pc;
            pc++;
            stm << "  ";
            printVal(stm, res);
            stm << " = " << opName(op) << " ";
            for (int i = 0;i < nargs;i++) {
                if (i != 0) {
                    stm << ", [ L";
                } else {
                    stm << "[ L";
                }
                stm << pc[2 * i] << ": ";
                printVal(stm, pc[2 * i + 1]);
                stm << " ]";
            }
            pc += 2 * nargs;
            stm << std::endl;
            break;
        }
        case Opcode::Call: {
            auto res = *pc;
            pc++;
            auto id = Builtins(*pc);
            pc++;
            auto nargs = *pc;
            pc++;
            stm << "  ";
            printVal(stm, res);
            stm << " = " << opName(op) << " " << builtinName(id) << "(";
            for (int i = 0;i < nargs;i++) {
                if (i != 0) {
                    stm << ", ";
                }
                printVal(stm, pc[i]);
            }
            pc += nargs;
            stm << ")" << std::endl;
            break;
        }
        case Opcode::Interp: {
            auto res = *pc;
            pc++;
            auto input = *pc;
            pc++;
            auto x0 = evalConst(*pc).get<double>();
            pc++;
            auto dx = evalConst(*pc).get<double>();
            pc++;
            auto data = *pc;
            pc++;
            auto ndata = *pc;
            pc++;
            stm << "  ";
            printVal(stm, res);
            stm << " = " << opName(op)
                << " [" << x0 << ", (" << ndata << ") +" << dx << "] (";
            printVal(stm, input);
            stm << ") {";
            auto datap = &float_table[data];
            for (int32_t i = 0; i < ndata; i++) {
                if (i != 0)
                    stm << ", ";
                stm << datap[i];
            }
            stm << "}" << std::endl;
            break;
        }
        default:
            stm << "  unknown op: " << uint8_t(op) << std::endl;
            break;
        }
    }
}

NACS_EXPORT() std::ostream &operator<<(std::ostream &stm, const Function &f)
{
    stm << typeName(f.ret) << " (";
    for (int i = 0;i < f.nargs;i++) {
        if (i != 0)
            stm << ", ";
        stm << typeName(f.valType(i)) << " ";
        f.printValName(stm, i);
    }
    stm << ") {" << std::endl;
    for (size_t i = 0;i < f.code.size();i++) {
        stm << "L" << i << ":" << std::endl;
        f.printBB(stm, f.code[i]);
    }
    stm << "}";
    return stm;
}

NACS_EXPORT() std::vector<uint32_t> Function::serialize(void) const
{
    // [ret][nargs][nvals][vals x nvals]
    // [nconsts][consts x nconsts]
    // [nbb][[nword][code x nword] x nbb]
    // <optional>[nfloat][double_data x nfloat]
    std::vector<uint32_t> res{uint32_t(ret), uint32_t(nargs)};
    auto copy_vector = [&res] (auto vec) {
        res.push_back(uint32_t(vec.size()));
        uint32_t idx = (uint32_t)res.size();
        size_t elsz = sizeof(typename decltype(vec)::value_type);
        res.resize(idx + (vec.size() * elsz + 3) / 4);
        memcpy(&res[idx], vec.data(), vec.size() * elsz);
    };
    copy_vector(vals);
    {
        res.push_back(uint32_t(consts.size()));
        uint32_t idx = (uint32_t)res.size();
        res.resize(idx + consts.size() * 3);
        for (size_t i = 0;i < consts.size();i++) {
            res[idx + i * 3] = uint32_t(consts[i].typ);
            memcpy(&res[idx + i * 3 + 1], &(consts[i].val), sizeof(GenVal));
        }
    }
    res.push_back(uint32_t(code.size()));
    for (size_t i = 0;i < code.size();i++)
        copy_vector(code[i]);
    copy_vector(float_table);
    return res;
}

NACS_EXPORT() Function::Function(const uint32_t *data, size_t sz)
    : ret(Type(data[0])),
      nargs(data[1]),
      vals{},
      code{},
      consts{},
      float_table{}
{
    uint32_t cursor = 2;
    auto read_vector = [data, &cursor] (auto &vec) {
        uint32_t size = data[cursor];
        cursor++;
        vec.resize(size);
        int32_t elsz = sizeof(typename std::remove_reference_t<
                              decltype(vec)>::value_type);
        memcpy(vec.data(), &data[cursor], size * elsz);
        cursor += (size * elsz + 3) / 4;
    };
    read_vector(vals);
    {
        uint32_t size = data[cursor];
        cursor++;
        consts.resize(size);
        for (size_t i = 0;i < size;i++) {
            consts[i].typ = Type(data[cursor + i * 3]);
            memcpy(&(consts[i].val), &data[cursor + i * 3 + 1], 8);
        }
        cursor += size * 3;
    }
    code.resize(data[cursor]);
    cursor++;
    for (size_t i = 0;i < code.size();i++)
        read_vector(code[i]);
    if (cursor < sz)
        read_vector(float_table);
    return;
}

NACS_EXPORT() Function::Function(const uint8_t *data, size_t sz)
    : ret(Type(Mem::load_unalign<uint32_t>(data, 0))),
      nargs(Mem::load_unalign<uint32_t>(data, 1)),
      vals{},
      code{},
      consts{},
      float_table{}
{
    sz = sz / sizeof(uint32_t);
    uint32_t cursor = 2;
    auto read_vector = [data, &cursor] (auto &vec) {
        uint32_t size = Mem::load_unalign<uint32_t>(data, cursor);
        cursor++;
        vec.resize(size);
        int32_t elsz = sizeof(typename std::remove_reference_t<
                              decltype(vec)>::value_type);
        memcpy(vec.data(), data + cursor * sizeof(uint32_t), size * elsz);
        cursor += (size * elsz + 3) / 4;
    };
    read_vector(vals);
    {
        uint32_t size = Mem::load_unalign<uint32_t>(data, cursor);
        cursor++;
        consts.resize(size);
        for (size_t i = 0;i < size;i++) {
            consts[i].typ = Type(Mem::load_unalign<uint32_t>(data, cursor + i * 3));
            memcpy(&(consts[i].val), data + (cursor + i * 3 + 1) * sizeof(uint32_t), 8);
        }
        cursor += size * 3;
    }
    code.resize(Mem::load_unalign<uint32_t>(data, cursor));
    cursor++;
    for (size_t i = 0;i < code.size();i++)
        read_vector(code[i]);
    if (cursor < sz)
        read_vector(float_table);
    return;
}

int32_t *Builder::addInst(Opcode op, size_t nop)
{
    Function::InstRef inst;
    return addInst(op, nop, inst);
}

int32_t *Builder::addInst(Opcode op, size_t nop, Function::InstRef &inst)
{
    auto &bb = m_f.code[m_cur_bb];
    auto oldlen = (int32_t)bb.size();
    inst.first = m_cur_bb;
    inst.second = oldlen + 1;
    bb.resize(oldlen + nop + 1);
    bb[oldlen] = uint32_t(op);
    return &bb[oldlen + 1];
}

void Builder::createRet(int32_t val)
{
    *addInst(Opcode::Ret, 1) = val;
}

int32_t Builder::getConstInt(int32_t val)
{
    auto map = const_ints;
    auto it = map.find(val);
    if (it != map.end())
        return it->second;
    int32_t oldlen = (int32_t)m_f.consts.size();
    m_f.consts.emplace_back(val);
    int32_t id = Consts::_Offset - oldlen;
    map[val] = id;
    return id;
}

int32_t Builder::getConstFloat(double val)
{
    auto map = const_floats;
    auto it = map.find(val);
    if (it != map.end())
        return it->second;
    int32_t oldlen = (int32_t)m_f.consts.size();
    m_f.consts.emplace_back(val);
    int32_t id = Consts::_Offset - oldlen;
    map[val] = id;
    return id;
}

int32_t Builder::getConst(TagVal val)
{
    switch (val.typ) {
    case Type::Bool:
        return val.val.b ? Consts::True : Consts::False;
    case Type::Int32:
        return getConstInt(val.val.i32);
    case Type::Float64:
        return getConstFloat(val.val.f64);
    default:
        return Consts::False;
    }
}

int32_t Builder::newSSA(Type typ)
{
    int32_t id = (int32_t)m_f.vals.size();
    m_f.vals.push_back(typ);
    return id;
}

int32_t Builder::newBB(void)
{
    int32_t id = (int32_t)m_f.code.size();
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
        int32_t *ptr = addInst(Opcode::Br, 2);
        ptr[0] = Consts::True;
        ptr[1] = bb1;
    }
    else {
        int32_t *ptr = addInst(Opcode::Br, 3);
        ptr[0] = cond;
        ptr[1] = bb1;
        ptr[2] = bb2;
    }
}

static TagVal evalAdd(Type typ, TagVal val1, TagVal val2)
{
    switch (typ) {
    case Type::Int32:
        return val1.get<int32_t>() + val2.get<int32_t>();
    case Type::Float64:
        return val1.get<double>() + val2.get<double>();
    default:
        return TagVal(typ);
    }
}

static TagVal evalSub(Type typ, TagVal val1, TagVal val2)
{
    switch (typ) {
    case Type::Int32:
        return val1.get<int32_t>() - val2.get<int32_t>();
    case Type::Float64:
        return val1.get<double>() - val2.get<double>();
    default:
        return TagVal(typ);
    }
}

static TagVal evalMul(Type typ, TagVal val1, TagVal val2)
{
    switch (typ) {
    case Type::Int32:
        return val1.get<int32_t>() * val2.get<int32_t>();
    case Type::Float64:
        return val1.get<double>() * val2.get<double>();
    default:
        return TagVal(typ);
    }
}

static TagVal evalFDiv(TagVal val1, TagVal val2)
{
    return val1.get<double>() / val2.get<double>();
}

static TagVal evalCmp(CmpType cmptyp, TagVal val1, TagVal val2)
{
    double v1 = val1.get<double>();
    double v2 = val2.get<double>();
    switch (cmptyp) {
    case CmpType::eq:
        return v1 == v2;
    case CmpType::gt:
        return v1 > v2;
    case CmpType::ge:
        return v1 >= v2;
    case CmpType::lt:
        return v1 < v2;
    case CmpType::le:
        return v1 <= v2;
    case CmpType::ne:
        return v1 != v2;
    default:
        return false;
    }
}

int32_t Builder::createPromoteOP(Opcode op, int32_t val1, int32_t val2)
{
    auto ty1 = m_f.valType(val1);
    auto ty2 = m_f.valType(val2);
    auto resty = std::max(std::max(ty1, ty2), Type::Int32);
    if (val1 < 0 && val2 < 0) {
        switch (op) {
        case Opcode::Add:
            return getConst(evalAdd(resty, m_f.evalConst(val1),
                                    m_f.evalConst(val2)));
        case Opcode::Sub:
            return getConst(evalSub(resty, m_f.evalConst(val1),
                                    m_f.evalConst(val2)));
        case Opcode::Mul:
            return getConst(evalMul(resty, m_f.evalConst(val1),
                                    m_f.evalConst(val2)));
        default:
            break;
        }
    }
    int32_t *ptr = addInst(op, 3);
    auto res = newSSA(resty);
    ptr[0] = res;
    ptr[1] = val1;
    ptr[2] = val2;
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

int32_t Builder::createFDiv(int32_t val1, int32_t val2)
{
    if (val1 < 0 && val2 < 0)
        return getConst(evalFDiv(m_f.evalConst(val1), m_f.evalConst(val2)));
    int32_t *ptr = addInst(Opcode::FDiv, 3);
    auto res = newSSA(Type::Float64);
    ptr[0] = res;
    ptr[1] = val1;
    ptr[2] = val2;
    return res;
}

int32_t Builder::createCmp(CmpType cmptyp, int32_t val1, int32_t val2)
{
    if (val1 < 0 && val2 < 0)
        return getConst(evalCmp(cmptyp, m_f.evalConst(val1),
                                m_f.evalConst(val2)));
    int32_t *ptr = addInst(Opcode::Cmp, 4);
    auto res = newSSA(Type::Bool);
    ptr[0] = res;
    ptr[1] = uint32_t(cmptyp);
    ptr[2] = val1;
    ptr[3] = val2;
    return res;
}

std::pair<int32_t, Function::InstRef> Builder::createPhi(Type typ, int ninputs)
{
    Function::InstRef inst;
    int32_t *ptr = addInst(Opcode::Phi, ninputs * 2 + 2, inst);
    auto res = newSSA(typ);
    ptr[0] = res;
    ptr[1] = ninputs;
    memset(&ptr[2], 0xff, ninputs * 2 * 4);
    return std::make_pair(res, inst);
}

int32_t Builder::createCall(Builtins id, int32_t nargs, const int32_t *args)
{
    switch (getBuiltinType(id)) {
    case BuiltinType::F64_F64:
        if (nargs != 1)
            return getConstFloat(0);
        break;
    case BuiltinType::F64_F64F64:
        if (nargs != 2)
            return getConstFloat(0);
        break;
    case BuiltinType::F64_F64F64F64:
        if (nargs != 3)
            return getConstFloat(0);
        break;
    case BuiltinType::F64_F64I32:
        if (nargs != 2)
            return getConstFloat(0);
        break;
    case BuiltinType::F64_I32F64:
        if (nargs != 2)
            return getConstFloat(0);
        break;
    default:
        return getConstFloat(0);
    }
    bool allconst = true;
    for (int i = 0;i < nargs;i++) {
        if (args[i] >= 0) {
            allconst = false;
            break;
        }
    }
    if (allconst) {
        TagVal constargs[3];
        assert(nargs <= 3);
        for (int i = 0;i < nargs;i++)
            constargs[i] = m_f.evalConst(args[i]);
        return getConstFloat(evalBuiltin(id, constargs));
    }
    int32_t *ptr = addInst(Opcode::Call, nargs + 3);
    auto res = newSSA(Type::Float64);
    ptr[0] = res;
    ptr[1] = uint32_t(id);
    ptr[2] = nargs;
    memcpy(&ptr[3], args, nargs * 4);
    return res;
}

int32_t Builder::createInterp(int32_t v, double x0, double dx, uint32_t npoints,
                              const double *points)
{
    auto x0id = getConstFloat(x0);
    auto dxid = getConstFloat(dx);
    auto data_id = addFloatData(points, npoints);
    auto *ptr = addInst(Opcode::Interp, 6);
    auto res = newSSA(Type::Float64);
    ptr[0] = res;
    ptr[1] = v;
    ptr[2] = x0id;
    ptr[3] = dxid;
    ptr[4] = data_id;
    ptr[5] = (int32_t)npoints;
    return res;
}

int32_t Builder::addFloatData(const double *data, uint32_t ndata)
{
    auto &float_table = m_f.float_table;
    auto res = (int32_t)float_table.size();
    float_table.insert(float_table.end(), data, data + ndata);
    return res;
}

void Builder::addPhiInput(Function::InstRef phi, int32_t bb, int32_t val)
{
    int32_t *inst = &m_f.code[phi.first][phi.second];
    int32_t nargs = inst[1];
    for (int32_t i = 0;i < nargs;i++) {
        if (inst[2 + 2 * i] == bb || inst[2 + 2 * i] == -1) {
            inst[2 + 2 * i] = bb;
            inst[2 + 2 * i + 1] = val;
            break;
        }
    }
}

static inline TagVal eval_val(const Function &f, const GenVal *vals, int32_t id)
{
    if (id >= 0) {
        return TagVal(f.vals[id], vals[id]);
    } else {
        return f.evalConst(id);
    }
}

static inline TagVal eval_func(const Function &f, GenVal *vals)
{
    int32_t bb_num = -1;
    int32_t prev_bb_num;
    const int32_t *pc;
    const int32_t *end;
    auto enter_bb = [&] (int32_t i) {
        prev_bb_num = bb_num;
        bb_num = i;
        auto &bb = f.code[i];
        pc = bb.data();
        end = pc + bb.size();
    };
    enter_bb(0);

    while (end > pc) {
        auto op = Opcode(*pc);
        pc++;
        auto res = *pc;
        pc++;
        auto &res_slot = vals[res];
        switch (op) {
        case Opcode::Ret:
            return eval_val(f, vals, res).convert(f.ret);
        case Opcode::Br:
            if (eval_val(f, vals, res).get<bool>()) {
                enter_bb(pc[0]);
            } else {
                enter_bb(pc[1]);
            }
            continue;
        case Opcode::Add:
        case Opcode::Sub:
        case Opcode::Mul:
        case Opcode::FDiv: {
            auto val1 = eval_val(f, vals, *pc);
            pc++;
            auto val2 = eval_val(f, vals, *pc);
            pc++;
            switch (op) {
            case Opcode::Add:
                res_slot = evalAdd(f.vals[res], val1, val2).val;
                break;
            case Opcode::Sub:
                res_slot = evalSub(f.vals[res], val1, val2).val;
                break;
            case Opcode::Mul:
                res_slot = evalMul(f.vals[res], val1, val2).val;
                break;
            case Opcode::FDiv:
                res_slot = evalFDiv(val1, val2).val;
                break;
            default:
                break;
            }
            break;
        }
        case Opcode::Cmp: {
            auto cmptyp = CmpType(*pc);
            pc++;
            auto val1 = eval_val(f, vals, *pc);
            pc++;
            auto val2 = eval_val(f, vals, *pc);
            pc++;
            res_slot = evalCmp(cmptyp, val1, val2).val;
            break;
        }
        case Opcode::Phi: {
            auto nargs = *pc;
            pc++;
            auto args = pc;
            pc += 2 * nargs;
            for (int i = 0;i < nargs;i++) {
                if (args[2 * i] == prev_bb_num) {
                    auto val = eval_val(f, vals, args[2 * i + 1]);
                    res_slot = val.convert(f.vals[res]).val;
                    break;
                }
            }
            break;
        }
        case Opcode::Call: {
            auto id = Builtins(*pc);
            pc++;
            auto nargs = *pc;
            pc++;
            TagVal argvals[3];
            assert(nargs <= 3);
            for (int i = 0;i < nargs;i++)
                argvals[i] = eval_val(f, vals, pc[i]);
            pc += nargs;
            res_slot = TagVal(evalBuiltin(id, argvals)).val;
            break;
        }
        case Opcode::Interp: {
            auto input = eval_val(f, vals, *pc).get<double>();
            pc++;
            auto x0 = eval_val(f, vals, *pc).get<double>();
            pc++;
            auto dx = eval_val(f, vals, *pc).get<double>();
            pc++;
            auto datap = &f.float_table[*pc];
            pc++;
            auto ndata = *pc;
            pc++;
            res_slot = TagVal(linearInterpolate(input, x0, dx, ndata, datap)).val;
            break;
        }
        default:
            break;
        }
    }
    return TagVal(f.ret);
}

TagVal EvalContext::evalVal(int32_t id) const
{
    return eval_val(m_f, &m_vals[0], id);
}

TagVal EvalContext::eval(void)
{
    return eval_func(m_f, &m_vals[0]);
}

static inline Function *get_interp_func(uint32_t *data)
{
    auto nargs = data[1];
    uint32_t offset = nargs + 2;
    offset = (offset + 3) & ~(uint32_t)3; // align to 16bytes
    return (Function*)&data[offset];
}

struct InterpExeContext : public ExeContext {
    static void exefunc_free(void *data)
    {
        auto func = get_interp_func((uint32_t*)data);
        func->~Function();
        ::free(data);
    }

    static inline ExeFuncBase getFuncIntern(Function f)
    {
        uint32_t nargs = f.nargs;
        uint32_t offset = nargs + 2;
        offset = (offset + 3) & ~(uint32_t)3; // align to 16bytes
        size_t sz = offset * 4 + sizeof(f);
        auto ptr = (uint32_t*)malloc(sz);
        ptr[0] = (f.vals.size() * 8 + 8) & ~(uint32_t)15;
        ptr[1] = nargs;
#if defined(__x86_64__) || defined(__x86_64)
#  ifdef NACS_OS_WINDOWS
        uint32_t nregarg = 0;
        uint32_t nstack = 0;
        for (uint32_t i = 0; i < nargs; i++) {
            auto ty = f.vals[i];
            if (ty == Type::Float64) {
                if (nregarg < 3) {
                    ptr[2 + i] = 3 + nregarg;
                    nregarg += 1;
                    continue;
                }
            }
            else {
                if (nregarg < 3) {
                    ptr[2 + i] = nregarg;
                    nregarg += 1;
                    continue;
                }
            }
            ptr[2 + i] = 48 + nstack * 8;
            nstack += 1;
        }
#else
        uint32_t nintarg = 0;
        uint32_t nfloatarg = 0;
        uint32_t nstack = 0;
        for (uint32_t i = 0; i < nargs; i++) {
            auto ty = f.vals[i];
            if (ty == Type::Float64) {
                if (nfloatarg < 6) {
                    ptr[2 + i] = 5 + nfloatarg;
                    nfloatarg += 1;
                    continue;
                }
            }
            else {
                if (nintarg < 5) {
                    ptr[2 + i] = nintarg;
                    nintarg += 1;
                    continue;
                }
            }
            ptr[2 + i] = 16 + nstack * 8;
            nstack += 1;
        }
#endif
#elif defined(__i386) || defined(__i386__)
        // 4bytes data pointer, 4 bytes return address and 4 bytes base pointer spill
        uint32_t stack_offset = 12;
        for (uint32_t i = 0; i < nargs; i++) {
            auto ty = f.vals[i];
            if (ty == Type::Float64) {
                ptr[2 + i] = stack_offset | 0x80000000;
                stack_offset += 8;
            }
            else {
                ptr[2 + i] = stack_offset;
                stack_offset += 4;
            }
        }
#elif defined(__aarch64__)
        uint32_t nintarg = 0;
        uint32_t nfloatarg = 0;
        uint32_t nstack = 0;
        for (uint32_t i = 0; i < nargs; i++) {
            auto ty = f.vals[i];
            if (ty == Type::Float64) {
                if (nfloatarg < 8) {
                    ptr[2 + i] = 7 + nfloatarg;
                    nfloatarg += 1;
                    continue;
                }
            }
            else {
                if (nintarg < 7) {
                    ptr[2 + i] = nintarg;
                    nintarg += 1;
                    continue;
                }
            }
            ptr[2 + i] = 16 + nstack * 8;
            nstack += 1;
        }
#elif defined(__arm__)
        uint32_t nintarg = 0;
        uint32_t nfloatarg = 0;
        uint32_t stack_offset = 16;
        for (uint32_t i = 0; i < nargs; i++) {
            auto ty = f.vals[i];
            if (ty == Type::Float64) {
                if (nfloatarg < 8) {
                    ptr[2 + i] = 3 + nfloatarg;
                    nfloatarg += 1;
                }
                else {
                    ptr[2 + i] = stack_offset | 0x80000000;
                    stack_offset += 8;
                }
            }
            else {
                if (nintarg < 3) {
                    ptr[2 + i] = nintarg;
                    nintarg += 1;
                }
                else {
                    ptr[2 + i] = stack_offset;
                    stack_offset += 4;
                }
            }
        }
#else
#  error "Unsupported architecture"
#endif
        new (get_interp_func(ptr)) Function(std::move(f));
        return ExeFuncBase(nacs_exefunc_cb, ptr, exefunc_free);
    }

    ExeFuncBase getFuncBase(const Function &f) override
    {
        return getFuncIntern(f);
    }
    ExeFuncBase getFuncBase(Function &&f) override
    {
        return getFuncIntern(std::move(f));
    }
};

NACS_EXPORT() std::unique_ptr<ExeContext> ExeContext::get()
{
    return std::unique_ptr<ExeContext>(new InterpExeContext());
}

}
}

using namespace NaCs::IR;

extern "C" TagVal nacs_exefunc_real(uint32_t *data, GenVal *vals)
{
    auto v = eval_func(*get_interp_func(data), vals);
    if (v.typ == Type::Float64)
        return v;
    if (v.typ == Type::Int32) {
        memset((char*)&v.val + 4, 0, 4);
    }
    else {
        memset((char*)&v.val + 1, 0, 7);
    }
    return v;
}
