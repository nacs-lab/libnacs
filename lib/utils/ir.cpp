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

#include "ir_p.h"
#include "interp_p.h"
#include "mem.h"
#include "dlload.h"

#include "nacs_lib_names.h"

#include <math.h>
#include <iostream>

#ifndef NACS_HAS_EXP10
// This is apparently missing on windows.
NACS_EXPORT() double nacs_exp10(double x)
{
    return pow(10, x);
}
#endif

namespace NaCs {
namespace IR {

void *get_openlibm_handle()
{
    static void *hdl = DL::open(NACS_OPENLIBM_NAME, DL::GLOBAL | DL::LAZY);
    return hdl;
}

namespace {

class LazyLibmFPtr {
public:
    LazyLibmFPtr(const char *name, fptr_t fallback)
        : m_name(name),
          m_fallback(fallback),
          m_fptr(nullptr)
    {
        assert(fallback);
    }
    operator fptr_t()
    {
        if (m_fptr)
            return m_fptr;
        fptr_t fptr = nullptr;
        if (auto hdl = get_openlibm_handle())
            fptr = (fptr_t)DL::sym(hdl, m_name);
        if (!fptr)
            fptr = m_fallback;
        m_fptr = fptr;
        return fptr;
    }

private:
    const char *m_name;
    const fptr_t m_fallback;
    fptr_t m_fptr;
};

#define DEF_FPTR(name) static LazyLibmFPtr p##name(#name, ::name)
DEF_FPTR(acos);
DEF_FPTR(acosh);
DEF_FPTR(asin);
DEF_FPTR(asinh);
DEF_FPTR(atan);
DEF_FPTR(atanh);
DEF_FPTR(cbrt);
DEF_FPTR(ceil);
DEF_FPTR(cos);
DEF_FPTR(cosh);
DEF_FPTR(erf);
DEF_FPTR(erfc);
DEF_FPTR(exp);
#ifdef NACS_HAS_EXP10
DEF_FPTR(exp10);
#else
static LazyLibmFPtr pexp10("exp10", ::nacs_exp10);
#endif
DEF_FPTR(exp2);
DEF_FPTR(expm1);
DEF_FPTR(fabs);
DEF_FPTR(floor);
DEF_FPTR(tgamma);
DEF_FPTR(j0);
DEF_FPTR(j1);
DEF_FPTR(lgamma);
DEF_FPTR(log);
DEF_FPTR(log10);
DEF_FPTR(log1p);
DEF_FPTR(log2);
DEF_FPTR(rint);
DEF_FPTR(round);
DEF_FPTR(sin);
DEF_FPTR(sinh);
DEF_FPTR(sqrt);
DEF_FPTR(tan);
DEF_FPTR(tanh);
DEF_FPTR(y0);
DEF_FPTR(y1);
#undef DEF_FPTR

#define DEF_FPTR(name)                                                  \
    static LazyLibmFPtr p##name(#name,                                  \
                                (fptr_t)(uintptr_t)(double(*)(double, double))(::name))
DEF_FPTR(atan2);
DEF_FPTR(copysign);
DEF_FPTR(fdim);
DEF_FPTR(fmax);
DEF_FPTR(fmin);
DEF_FPTR(fmod);
DEF_FPTR(hypot);
DEF_FPTR(pow);
DEF_FPTR(remainder);
#undef DEF_FPTR

static LazyLibmFPtr pfma("fma", (fptr_t)(uintptr_t)(double(*)(double, double, double))(::fma));
static LazyLibmFPtr pldexp("ldexp", (fptr_t)(uintptr_t)(double(*)(double, int))(::ldexp));
static LazyLibmFPtr pjn("jn", (fptr_t)(uintptr_t)(double(*)(int, double))(::jn));
static LazyLibmFPtr pyn("yn", (fptr_t)(uintptr_t)(double(*)(int, double))(::yn));

}

static fptr_t getBuiltinPtr(Builtins id)
{
    switch (id) {
        // f(f)
    case Builtins::acos:
        return pacos;
    case Builtins::acosh:
        return pacosh;
    case Builtins::asin:
        return pasin;
    case Builtins::asinh:
        return pasinh;
    case Builtins::atan:
        return patan;
    case Builtins::atanh:
        return patanh;
    case Builtins::cbrt:
        return pcbrt;
    case Builtins::ceil:
        return pceil;
    case Builtins::cos:
        return pcos;
    case Builtins::cosh:
        return pcosh;
    case Builtins::erf:
        return perf;
    case Builtins::erfc:
        return perfc;
    case Builtins::exp:
        return pexp;
    case Builtins::exp10:
        return pexp10;
    case Builtins::exp2:
        return pexp2;
    case Builtins::expm1:
        return pexpm1;
    case Builtins::abs:
        return pfabs;
    case Builtins::floor:
        return pfloor;
    case Builtins::gamma: // tgamma
        return ptgamma;
    case Builtins::j0:
        return pj0;
    case Builtins::j1:
        return pj1;
    case Builtins::lgamma:
        return plgamma;
    case Builtins::log:
        return plog;
    case Builtins::log10:
        return plog10;
    case Builtins::log1p:
        return plog1p;
    case Builtins::log2:
        return plog2;
    case Builtins::pow10:
        return pexp10;
    case Builtins::rint:
        return print;
    case Builtins::round:
        return pround;
    case Builtins::sin:
        return psin;
    case Builtins::sinh:
        return psinh;
    case Builtins::sqrt:
        return psqrt;
    case Builtins::tan:
        return ptan;
    case Builtins::tanh:
        return ptanh;
    case Builtins::y0:
        return py0;
    case Builtins::y1:
        return py1;

        // f(f, f)
    case Builtins::atan2:
        return patan2;
    case Builtins::copysign:
        return pcopysign;
    case Builtins::fdim:
        return pfdim;
    case Builtins::max: // fmax
        return pfmax;
    case Builtins::min: // fmax
        return pfmin;
    case Builtins::mod: // fmod
        return pfmod;
    case Builtins::hypot:
        return phypot;
    case Builtins::pow:
        return ppow;
    case Builtins::remainder:
        return premainder;

        // f(f, f, f)
    case Builtins::fma:
        return pfma;

        // f(f, i)
    case Builtins::ldexp:
        return pldexp;

        // f(i, f)
    case Builtins::jn:
        return pjn;
    case Builtins::yn:
        return pyn;
    default:
        return nullptr;
    }
}

const char *getBuiltinSymbol(Builtins id)
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
    case Builtins::gamma: // tgamma
        return "tgamma";
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
        return "fmax";
    case Builtins::min: // fmax
        return "fmin";
    case Builtins::mod: // fmod
        return "fmod";
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
        return nullptr;
    }
}

static const char *builtinName(Builtins id)
{
    switch (id) {
    case Builtins::gamma:
        return "gamma";
    case Builtins::max:
        return "max";
    case Builtins::min:
        return "min";
    case Builtins::mod:
        return "mod";
    default:
        if (auto name = getBuiltinSymbol(id))
            return name;
        return "unknown";
    }
}

BuiltinType getBuiltinType(Builtins id)
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
        auto fptr = (double(*)(double, double))(uintptr_t)_fptr;
        return fptr(args[0].get<double>(), args[1].get<double>());
    }
    case BuiltinType::F64_F64F64F64: {
        auto fptr = (double(*)(double, double, double))(uintptr_t)_fptr;
        return fptr(args[0].get<double>(), args[1].get<double>(),
                    args[2].get<double>());
    }
    case BuiltinType::F64_F64I32: {
        auto fptr = (double(*)(double, int))(uintptr_t)_fptr;
        return fptr(args[0].get<double>(), args[1].get<int32_t>());
    }
    case BuiltinType::F64_I32F64: {
        auto fptr = (double(*)(int, double))(uintptr_t)_fptr;
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
    case Opcode::Convert:
        return "convert";
    case Opcode::Select:
        return "select";
    case Opcode::And:
        return "and";
    case Opcode::Or:
        return "or";
    case Opcode::Xor:
        return "xor";
    case Opcode::Not:
        return "not";
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
        case Opcode::Convert: {
            auto res = *pc;
            pc++;
            auto input = *pc;
            pc++;
            stm << "  ";
            printVal(stm, res);
            stm << " = " << opName(op) << "(";
            printVal(stm, input);
            stm << ")" << std::endl;
            break;
        }
        case Opcode::Select: {
            auto res = *pc;
            pc++;
            auto cond = *pc;
            pc++;
            auto v1 = *pc;
            pc++;
            auto v2 = *pc;
            pc++;
            stm << "  ";
            printVal(stm, res);
            stm << " = " << opName(op) << " ";
            printVal(stm, cond);
            stm << ", ";
            printVal(stm, v1);
            stm << ", ";
            printVal(stm, v2);
            stm << std::endl;
            break;
        }
        case Opcode::And:
        case Opcode::Or:
        case Opcode::Xor: {
            auto res = *pc;
            pc++;
            auto v1 = *pc;
            pc++;
            auto v2 = *pc;
            pc++;
            stm << "  ";
            printVal(stm, res);
            stm << " = " << opName(op) << " ";
            printVal(stm, v1);
            stm << ", ";
            printVal(stm, v2);
            stm << std::endl;
            break;
        }
        case Opcode::Not: {
            auto res = *pc;
            pc++;
            auto v = *pc;
            pc++;
            stm << "  ";
            printVal(stm, res);
            stm << " = " << opName(op) << " ";
            printVal(stm, v);
            stm << std::endl;
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

NACS_EXPORT() void Builder::createRet(int32_t val)
{
    *addInst(Opcode::Ret, 1) = val;
}

NACS_EXPORT() int32_t Builder::getConstInt(int32_t val)
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

NACS_EXPORT() int32_t Builder::getConstFloat(double val)
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

NACS_EXPORT() int32_t Builder::getConst(TagVal val)
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

NACS_EXPORT() int32_t Builder::newBB(void)
{
    int32_t id = (int32_t)m_f.code.size();
    m_f.code.push_back({});
    return id;
}

NACS_EXPORT() int32_t &Builder::curBB()
{
    return m_cur_bb;
}

NACS_EXPORT() void Builder::createBr(int32_t br)
{
    createBr(Consts::True, br, 0);
}

NACS_EXPORT() void Builder::createBr(int32_t cond, int32_t bb1, int32_t bb2)
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

NACS_EXPORT() int32_t Builder::createAdd(int32_t val1, int32_t val2)
{
    return createPromoteOP(Opcode::Add, val1, val2);
}

NACS_EXPORT() int32_t Builder::createSub(int32_t val1, int32_t val2)
{
    return createPromoteOP(Opcode::Sub, val1, val2);
}

NACS_EXPORT() int32_t Builder::createMul(int32_t val1, int32_t val2)
{
    return createPromoteOP(Opcode::Mul, val1, val2);
}

NACS_EXPORT() int32_t Builder::createFDiv(int32_t val1, int32_t val2)
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

NACS_EXPORT() int32_t Builder::createCmp(CmpType cmptyp, int32_t val1, int32_t val2)
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

NACS_EXPORT() std::pair<int32_t, Function::InstRef>
Builder::createPhi(Type typ, int ninputs)
{
    Function::InstRef inst;
    int32_t *ptr = addInst(Opcode::Phi, ninputs * 2 + 2, inst);
    auto res = newSSA(typ);
    ptr[0] = res;
    ptr[1] = ninputs;
    memset(&ptr[2], 0xff, ninputs * 2 * 4);
    return std::make_pair(res, inst);
}

NACS_EXPORT() int32_t Builder::createCall(Builtins id, int32_t nargs, const int32_t *args)
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

NACS_EXPORT() int32_t Builder::createInterp(int32_t v, double x0, double dx,
                                            uint32_t npoints, const double *points)
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

NACS_EXPORT() int32_t Builder::createConvert(Type typ, int32_t v)
{
    if (m_f.valType(v) == typ)
        return v;
    if (v < 0)
        return getConst(m_f.evalConst(v).convert(typ));
    auto *ptr = addInst(Opcode::Convert, 2);
    auto res = newSSA(typ);
    ptr[0] = res;
    ptr[1] = v;
    return res;
}

NACS_EXPORT() int32_t Builder::createSelect(int32_t cond, int32_t val1, int32_t val2)
{
    auto ty1 = m_f.valType(val1);
    auto ty2 = m_f.valType(val2);
    auto resty = std::max(std::max(ty1, ty2), Type::Int32);
    if (cond < 0)
        return createConvert(resty, m_f.evalConst(cond).get<bool>() ? val1 : val2);
    auto *ptr = addInst(Opcode::Select, 4);
    auto res = newSSA(resty);
    ptr[0] = res;
    ptr[1] = cond;
    ptr[2] = val1;
    ptr[3] = val2;
    return res;
}

NACS_EXPORT() int32_t Builder::createAnd(int32_t val1, int32_t val2)
{
    if (val1 < 0) {
        auto c1 = m_f.evalConst(val1).get<bool>();
        if (!c1)
            return Consts::False;
        return createConvert(Type::Bool, val2);
    }
    else if (val2 < 0) {
        auto c2 = m_f.evalConst(val2).get<bool>();
        if (!c2)
            return Consts::False;
        return createConvert(Type::Bool, val1);
    }
    auto *ptr = addInst(Opcode::And, 3);
    auto res = newSSA(Type::Bool);
    ptr[0] = res;
    ptr[1] = val1;
    ptr[2] = val2;
    return res;
}

NACS_EXPORT() int32_t Builder::createOr(int32_t val1, int32_t val2)
{
    if (val1 < 0) {
        auto c1 = m_f.evalConst(val1).get<bool>();
        if (c1)
            return Consts::True;
        return createConvert(Type::Bool, val2);
    }
    else if (val2 < 0) {
        auto c2 = m_f.evalConst(val2).get<bool>();
        if (c2)
            return Consts::True;
        return createConvert(Type::Bool, val1);
    }
    auto *ptr = addInst(Opcode::Or, 3);
    auto res = newSSA(Type::Bool);
    ptr[0] = res;
    ptr[1] = val1;
    ptr[2] = val2;
    return res;
}

NACS_EXPORT() int32_t Builder::createXor(int32_t val1, int32_t val2)
{
    if (val1 < 0) {
        auto c1 = m_f.evalConst(val1).get<bool>();
        if (c1)
            return createConvert(Type::Bool, val2);
        return createNot(val2);
    }
    else if (val2 < 0) {
        auto c2 = m_f.evalConst(val2).get<bool>();
        if (c2)
            return createConvert(Type::Bool, val1);
        return createNot(val1);
    }
    auto *ptr = addInst(Opcode::Xor, 3);
    auto res = newSSA(Type::Bool);
    ptr[0] = res;
    ptr[1] = val1;
    ptr[2] = val2;
    return res;
}

NACS_EXPORT() int32_t Builder::createNot(int32_t val)
{
    if (val < 0) {
        auto c1 = m_f.evalConst(val).get<bool>();
        return c1 ? Consts::True : Consts::False;
    }
    auto *ptr = addInst(Opcode::Not, 2);
    auto res = newSSA(Type::Bool);
    ptr[0] = res;
    ptr[1] = val;
    return res;
}

int32_t Builder::addFloatData(const double *data, uint32_t ndata)
{
    auto &float_table = m_f.float_table;
    auto res = (int32_t)float_table.size();
    float_table.insert(float_table.end(), data, data + ndata);
    return res;
}

NACS_EXPORT() void Builder::addPhiInput(Function::InstRef phi, int32_t bb, int32_t val)
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

NACS_EXPORT() std::unique_ptr<ExeContext> ExeContext::get()
{
    return std::unique_ptr<ExeContext>(get_interp_context());
}

}
}
