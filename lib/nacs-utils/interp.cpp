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

#include "interp_p.h"

// This is the assembly function defined in `ir-interp-shim.S` that grabs
// function arguments from the correct location and rearange them into a format
// that can be understand by the interpreter.
// This will be used as the function pointer returned by the interpreter execution context.
extern "C" void nacs_exefunc_cb(void);

namespace NaCs::IR {

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
        case Opcode::Convert: {
            auto input = eval_val(f, vals, *pc);
            pc++;
            res_slot = input.convert(f.vals[res]).val;
            break;
        }
        case Opcode::Select: {
            auto cond = eval_val(f, vals, *pc).get<bool>();
            pc++;
            auto val = eval_val(f, vals, cond ? pc[0] : pc[1]);
            pc += 2;
            res_slot = val.convert(f.vals[res]).val;
            break;
        }
        case Opcode::And:
        case Opcode::Or:
        case Opcode::Xor: {
            auto v1 = eval_val(f, vals, *pc).get<bool>();
            pc++;
            auto v2 = eval_val(f, vals, *pc).get<bool>();
            pc++;
            bool v = false;
            if (op == Opcode::And) {
                v = v1 && v2;
            }
            else if (op == Opcode::Or) {
                v = v1 || v2;
            }
            else if (op == Opcode::Xor) {
                v = v1 xor v2;
            }
            res_slot = TagVal(v).val;
            break;
        }
        case Opcode::Not: {
            auto v = eval_val(f, vals, *pc).get<bool>();
            pc++;
            res_slot = TagVal(!v).val;
            break;
        }
        default:
            break;
        }
    }
    return TagVal(f.ret);
}

namespace {

static inline Function *get_interp_func(uint32_t *data, uint32_t nargs)
{
    uint32_t offset = nargs + 2;
    offset = (offset + 3) & ~(uint32_t)3; // align to 16bytes
    return (Function*)&data[offset];
}

static inline Function *get_interp_func(uint32_t *data)
{
    return get_interp_func(data, data[1]);
}

template<typename T> static T getV(GenVal);
template<> inline bool getV<bool>(GenVal v)
{
    return v.b;
}
template<> inline int32_t getV<int32_t>(GenVal v)
{
    return v.i32;
}
template<> inline double getV<double>(GenVal v)
{
    return v.f64;
}

static inline void set_args(GenVal*)
{
}

template<typename Arg1, typename... Args>
static inline void set_args(GenVal *vals, Arg1 arg1, Args... args)
{
    *vals = TagVal(arg1).val;
    set_args(vals + 1, args...);
}

// Template to generate a number of specialized callback functions at compile time.
template<int n, int max, typename... Args>
struct Dispatch {
    static_assert(n + sizeof...(Args) == max, "");
    static inline void (*check(Function &f))(void)
    {
        if (f.nargs > max)
            return nullptr;
        static constexpr int i = max - n;
        if (i >= f.nargs)
            return Dispatch<-1,max,Args...>::check(f);
        switch (f.vals[i]) {
        case Type::Bool:
            return Dispatch<n - 1,max,Args...,bool>::check(f);
        case Type::Int32:
            return Dispatch<n - 1,max,Args...,int32_t>::check(f);
        case Type::Float64:
            return Dispatch<n - 1,max,Args...,double>::check(f);
        default:
            return nullptr;
        }
    }
};

template<int max, typename... Args>
struct Dispatch<-1,max,Args...> {
    template<typename Ret>
    static Ret spec_execfunc_cb(void *data, Args... args)
    {
        uint32_t *data32 = (uint32_t*)data;
        GenVal vals[data32[0]];
        set_args(vals, args...);
        auto v = eval_func(*get_interp_func(data32, uint32_t(sizeof...(Args))), vals).val;
        return getV<Ret>(v);
    }
    static inline void (*check(Function &f))(void)
    {
        switch (f.ret) {
        case Type::Bool:
            return (void(*)())spec_execfunc_cb<bool>;
        case Type::Int32:
            return (void(*)())spec_execfunc_cb<int32_t>;
        case Type::Float64:
            return (void(*)())spec_execfunc_cb<double>;
        default:
            return nullptr;
        }
    }
};

struct InterpExeContext : public ExeContext {
    static void exefunc_free(void *data)
    {
        auto func = get_interp_func((uint32_t*)data);
        func->~Function();
        ::free(data);
    }

    static inline FuncBase getFuncIntern(Function f)
    {
        uint32_t nargs = f.nargs;
        uint32_t offset = nargs + 2;
        offset = (offset + 3) & ~(uint32_t)3; // align to 16bytes
        size_t sz = offset * 4 + sizeof(f);
        auto ptr = (uint32_t*)malloc(sz);
        ptr[0] = (f.vals.size() * 8 + 8) & ~(uint32_t)15;
        ptr[1] = nargs;
        if (auto fptr = Dispatch<3,3>::check(f)) {
            new (get_interp_func(ptr)) Function(std::move(f));
            return FuncBase(fptr, ptr, exefunc_free);
        }
        // Now use the generic version.
        // The format used below has to be synchronized with the assembly function.
#if defined(__x86_64__) || defined(__x86_64)
#  if NACS_OS_WINDOWS
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
                if (nfloatarg < 8) {
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
        uint32_t stack_offset = 0;
        // AAPCS64 requires 8 bytes alignment for all stack slots for arguments.
        // However, Apple deviates from this and do not require extra alignment.
        for (uint32_t i = 0; i < nargs; i++) {
            auto ty = f.vals[i];
            if (ty == Type::Float64) {
                if (nfloatarg < 8) {
                    ptr[2 + i] = 7 + nfloatarg;
                    nfloatarg += 1;
                    continue;
                }
#  if NACS_OS_DARWIN
                else {
                    stack_offset = alignTo(stack_offset, 8);
                    ptr[2 + i] = 16 + stack_offset;
                    stack_offset += 8;
                }
                continue;
#  endif
            }
            else {
                if (nintarg < 7) {
                    ptr[2 + i] = nintarg;
                    nintarg += 1;
                    continue;
                }
#  if NACS_OS_DARWIN
                else if (ty == Type::Int32) {
                    stack_offset = alignTo(stack_offset, 4);
                    ptr[2 + i] = 16 + stack_offset;
                    stack_offset += 4;
                }
                else {
                    ptr[2 + i] = 16 + stack_offset;
                    stack_offset += 1;
                }
                continue;
#  endif
            }
#  if !NACS_OS_DARWIN
            ptr[2 + i] = 16 + stack_offset;
            stack_offset += 8;
#  endif
        }
#elif defined(__arm__)
        uint32_t nintarg = 0;
        uint32_t nfloatarg = 0;
        uint32_t stack_offset = 24;
        for (uint32_t i = 0; i < nargs; i++) {
            auto ty = f.vals[i];
            if (ty == Type::Float64) {
                if (nfloatarg < 8) {
                    ptr[2 + i] = 3 + nfloatarg;
                    nfloatarg += 1;
                }
                else {
                    stack_offset = (stack_offset + 7) & ~7;
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
        return FuncBase(nacs_exefunc_cb, ptr, exefunc_free);
    }

    FuncBase getFuncBase(const Function &f) override
    {
        return getFuncIntern(f);
    }
    FuncBase getFuncBase(Function &&f) override
    {
        return getFuncIntern(std::move(f));
    }
};

}

ExeContext *get_interp_context()
{
    return new InterpExeContext();
}

} // NaCs::IR

using namespace NaCs::IR;
// This is called from the assembly callback function after rearranging the argument format.
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
