/*************************************************************************
 *   Copyright (c) 2018 - 2018 Yichao Yu <yyc1992@gmail.com>             *
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

#include "codegen.h"
#include "utils.h"
#include "compile_p.h"

#include "../ir_p.h"
#include "../number.h"

#include <llvm/ADT/SetVector.h>
#include <llvm/ADT/SmallVector.h>

namespace NaCs {
namespace LLVM {
namespace Codegen {

Context::Context(Module *mod)
    : m_mod(mod),
      m_ctx(mod->getContext()),
      T_bool(Type::getInt1Ty(m_ctx)),
      T_i8(Type::getInt8Ty(m_ctx)),
      T_i32(Type::getInt32Ty(m_ctx)),
      T_f64(Type::getDoubleTy(m_ctx)),
      F_f64_f64(FunctionType::get(T_f64, {T_f64}, false)),
      F_f64_f64f64(FunctionType::get(T_f64, {T_f64, T_f64}, false)),
      F_f64_f64f64f64(FunctionType::get(T_f64, {T_f64, T_f64, T_f64}, false)),
      F_f64_f64i32(FunctionType::get(T_f64, {T_f64, T_i32}, false)),
      F_f64_i32f64(FunctionType::get(T_f64, {T_i32, T_f64}, false)),
      F_f64_f64f64f64i32pf64(FunctionType::get(T_f64, {T_f64, T_f64, T_f64,
                      T_i32, T_f64->getPointerTo()}, false)),
      m_mdbuilder(m_ctx),
      tbaa_root(m_mdbuilder.createTBAARoot("nacs_tbaa"))
{
    MDNode *const_s = m_mdbuilder.createTBAAScalarTypeNode("nacs_tbaa_const", tbaa_root);
    tbaa_const = m_mdbuilder.createTBAAStructTagNode(const_s, const_s, 0, true);
}

Constant *Context::ensurePureFunc(StringRef name, FunctionType *ft, bool canread) const
{
    if (auto f = m_mod->getNamedValue(name)) {
        auto pft = ft->getPointerTo();
        if (f->getType() != pft)
            return ConstantExpr::getBitCast(f, pft);
        return f;
    }
    Function *f = Function::Create(ft, GlobalValue::ExternalLinkage, name, m_mod);
    f->addFnAttr(Attribute::NoRecurse);
    f->addFnAttr(Attribute::NoUnwind);
    if (canread) {
        f->addFnAttr(Attribute::ReadOnly);
        f->addFnAttr(Attribute::ArgMemOnly);
    }
    else {
        f->addFnAttr(Attribute::ReadNone);
    }
    return f;
}

Type *Context::llvm_ty(IR::Type ty) const
{
    switch (ty) {
    case IR::Type::Bool:
        return T_bool;
    case IR::Type::Int32:
        return T_i32;
    case IR::Type::Float64:
        return T_f64;
    default:
        abort();
    }
}

Type *Context::llvm_argty(IR::Type ty) const
{
    switch (ty) {
    case IR::Type::Bool:
        return T_i8;
    case IR::Type::Int32:
        return T_i32;
    case IR::Type::Float64:
        return T_f64;
    default:
        abort();
    }
}

Value *Context::emit_const(IR::TagVal c) const
{
    switch (c.typ) {
    case IR::Type::Bool:
        return ConstantInt::get(T_bool, c.val.b);
    case IR::Type::Int32:
        return ConstantInt::get(T_i32, c.val.i32);
    case IR::Type::Float64:
        return ConstantFP::get(T_f64, c.val.f64);
    default:
        abort();
    }
}

Value *Context::emit_convert(IRBuilder<> &builder, IR::Type ty, Value *val) const
{
    auto lty = val->getType();
    switch (ty) {
    case IR::Type::Bool:
    case IR::Type::Int32:
        if (lty->isIntegerTy())
            return builder.CreateZExtOrTrunc(val, llvm_ty(ty));
        assert(lty->isFloatingPointTy());
        if (ty == IR::Type::Bool)
            return builder.CreateFCmpOEQ(val, ConstantFP::get(lty, 0));
        return builder.CreateFPToSI(val, T_i32);
    case IR::Type::Float64:
        if (lty == T_bool)
            return builder.CreateUIToFP(val, T_f64);
        if (lty->isIntegerTy())
            return builder.CreateSIToFP(val, T_f64);
        assert(lty->isFloatingPointTy());
        return builder.CreateFPCast(val, T_f64);
    default:
        abort();
    }
}

Value *Context::emit_add(IRBuilder<> &builder, IR::Type ty, Value *val1, Value *val2) const
{
    switch (ty) {
    case IR::Type::Int32:
        return builder.CreateAdd(emit_convert(builder, ty, val1),
                                 emit_convert(builder, ty, val2));
    case IR::Type::Float64:
        return builder.CreateFAdd(emit_convert(builder, ty, val1),
                                  emit_convert(builder, ty, val2));
    default:
        return UndefValue::get(llvm_ty(ty));
    }
}

Value *Context::emit_sub(IRBuilder<> &builder, IR::Type ty, Value *val1, Value *val2) const
{
    switch (ty) {
    case IR::Type::Int32:
        return builder.CreateSub(emit_convert(builder, ty, val1),
                                 emit_convert(builder, ty, val2));
    case IR::Type::Float64:
        return builder.CreateFSub(emit_convert(builder, ty, val1),
                                  emit_convert(builder, ty, val2));
    default:
        return UndefValue::get(llvm_ty(ty));
    }
}

Value *Context::emit_mul(IRBuilder<> &builder, IR::Type ty, Value *val1, Value *val2) const
{
    switch (ty) {
    case IR::Type::Int32:
        return builder.CreateMul(emit_convert(builder, ty, val1),
                                 emit_convert(builder, ty, val2));
    case IR::Type::Float64:
        return builder.CreateFMul(emit_convert(builder, ty, val1),
                                  emit_convert(builder, ty, val2));
    default:
        return UndefValue::get(llvm_ty(ty));
    }
}

Value *Context::emit_fdiv(IRBuilder<> &builder, Value *val1, Value *val2) const
{
    return builder.CreateFDiv(emit_convert(builder, IR::Type::Float64, val1),
                              emit_convert(builder, IR::Type::Float64, val2));
}

Value *Context::emit_cmp(IRBuilder<> &builder, IR::CmpType cmptyp,
                         Value *val1, Value *val2) const
{
    if (val1->getType()->isIntegerTy() && val2->getType()->isIntegerTy()) {
        val1 = emit_convert(builder, IR::Type::Int32, val1);
        val2 = emit_convert(builder, IR::Type::Int32, val2);
        switch (cmptyp) {
        case IR::CmpType::eq:
            return builder.CreateICmpEQ(val1, val2);
        case IR::CmpType::gt:
            return builder.CreateICmpSGT(val1, val2);
        case IR::CmpType::ge:
            return builder.CreateICmpSGE(val1, val2);
        case IR::CmpType::lt:
            return builder.CreateICmpSLT(val1, val2);
        case IR::CmpType::le:
            return builder.CreateICmpSLE(val1, val2);
        case IR::CmpType::ne:
            return builder.CreateICmpNE(val1, val2);
        default:
            abort();
        }
    }
    val1 = emit_convert(builder, IR::Type::Float64, val1);
    val2 = emit_convert(builder, IR::Type::Float64, val2);
    switch (cmptyp) {
    case IR::CmpType::eq:
        return builder.CreateFCmpOEQ(val1, val2);
    case IR::CmpType::gt:
        return builder.CreateFCmpOGT(val1, val2);
    case IR::CmpType::ge:
        return builder.CreateFCmpOGE(val1, val2);
    case IR::CmpType::lt:
        return builder.CreateFCmpOLT(val1, val2);
    case IR::CmpType::le:
        return builder.CreateFCmpOLE(val1, val2);
    case IR::CmpType::ne:
        return builder.CreateFCmpONE(val1, val2);
    default:
        abort();
    }
}

NACS_EXPORT()
Function *Context::emit_wrapper(Function *func, StringRef name, const Wrapper &spec)
{
    // 1. Create function signature
    auto fty = func->getFunctionType();
    assert(!fty->isVarArg());
    auto nargs = fty->getNumParams();
    auto rt = fty->getReturnType();
    SmallVector<Type*, 8> fsig;
    for (int i = 0; i < nargs; i++) {
        auto arg_spec = spec.arg_map.find(i);
        auto argt = fty->getParamType(i);
        if (arg_spec == spec.arg_map.end()) {
            fsig.push_back(argt);
            continue;
        }
        if (arg_spec->second.type == Wrapper::Closure) {
            if (!spec.closure)
                return nullptr;
            // Nothing to do with this argument
        }
        else {
            // Unknown argument type
            return nullptr;
        }
    }
    if (spec.closure)
        fsig.push_back(T_i8->getPointerTo());
    auto wrapf_type = FunctionType::get(rt, fsig, false);

    // 2. Create function
    auto wrapf = Function::Create(wrapf_type, GlobalValue::ExternalLinkage, name, m_mod);
    wrapf->setVisibility(GlobalValue::ProtectedVisibility);
    wrapf->addFnAttr(Attribute::NoRecurse);
    wrapf->addFnAttr(Attribute::NoUnwind);
    wrapf->addFnAttr(Attribute::ReadOnly);

    BasicBlock *b0 = BasicBlock::Create(m_ctx, "top", wrapf);
    IRBuilder<> builder(b0);

    // 3. Create arguments
    SmallVector<Value*, 8> call_args(nargs);
    auto argit = wrapf->arg_end();
    auto clarg = spec.closure ? &*(--argit) : nullptr;
    argit = wrapf->arg_begin();
    uint32_t max_offset = 0;
    for (int i = 0, j = 0; i < nargs; i++) {
        auto arg_spec = spec.arg_map.find(i);
        if (arg_spec == spec.arg_map.end()) {
            call_args[i] = &*(argit + j);
            j++;
            continue;
        }
        // Only type for now.
        assert(arg_spec->second.type == Wrapper::Closure);
        auto offset = arg_spec->second.idx;
        auto argt = fty->getParamType(i);
        auto ptr = builder.CreateBitCast(builder.CreateConstGEP1_32(T_i8, clarg, offset * 8),
                                         argt->getPointerTo());
        max_offset = max(max_offset, offset + 1);
        auto load = builder.CreateLoad(argt, ptr);
        load->setAlignment(8);
        load->setMetadata(LLVMContext::MD_tbaa, tbaa_const);
        call_args[i] = load;
    }
    if (max_offset) {
        assert(spec.closure);
        clarg->addAttr(Attribute::NonNull);
        wrapf->addDereferenceableParamAttr(fsig.size() - 1, max_offset * 8);
    }
    auto res = builder.CreateCall(func, call_args);
    builder.CreateRet(res);
    return wrapf;
}

NACS_EXPORT()
Function *Context::emit_function(const IR::Function &func, StringRef name, bool _export)
{
    auto nargs = func.nargs;
    auto nslots = func.vals.size();
    // 1. Create function signature
    auto rt = llvm_argty(func.ret);
    SmallVector<Type*, 8> fsig(nargs);
    for (int i = 0; i < nargs; i++)
        fsig[i] = llvm_argty(func.vals[i]);
    auto ftype = FunctionType::get(rt, fsig, false);

    // 2. Create function
    auto f = Function::Create(ftype, _export ? GlobalValue::ExternalLinkage :
                              GlobalValue::InternalLinkage, name, m_mod);
    if (_export)
        f->setVisibility(GlobalValue::ProtectedVisibility);
    f->addFnAttr(Attribute::AlwaysInline);
    f->addFnAttr(Attribute::NoRecurse);
    f->addFnAttr(Attribute::NoUnwind);
    f->addFnAttr(Attribute::ReadOnly);

    BasicBlock *b0 = BasicBlock::Create(m_ctx, "top", f);
    IRBuilder<> builder(b0);
    FastMathFlags fmf;
    fmf.setAllowContract(true);
    builder.setFastMathFlags(fmf);

    // 3. Create variable slots
    SmallVector<AllocaInst*, 16> slots(nslots);
    for (unsigned i = 0; i < nslots; i++)
        slots[i] = builder.CreateAlloca(llvm_ty(func.vals[i]));
    auto argit = f->arg_begin();
    for (int i = 0; i < nargs; i++) {
        Value *argv = &*argit++;
        if (func.vals[i] == IR::Type::Bool)
            argv = builder.CreateTrunc(argv, T_bool);
        builder.CreateStore(argv, slots[i]);
    }
    auto prev_bb_var = builder.CreateAlloca(T_i32);

    // 4. Initialize BB info
    auto nbbs = func.code.size();
    SmallVector<BasicBlock*, 8> bbs(nbbs);
    // We allow branching to the entry BB so don't use it as the LLVM entry BB.
    for (unsigned i = 0; i < nbbs; i++)
        bbs[i] = BasicBlock::Create(m_ctx, "B" + std::to_string(i), f);
    builder.CreateBr(bbs[0]);

    SetVector<int, SmallVector<int, 8>> workstack;
    auto push_work = [&] (int i) {
        if (bbs[i]->getTerminator())
            return;
        workstack.insert(i);
    };

    const int32_t *pc;
    const int32_t *end;
    Value *prev_bb_num;
    auto enter_bb = [&] (int i) {
        auto &bb = func.code[i];
        pc = bb.data();
        end = pc + bb.size();
        builder.SetInsertPoint(bbs[i]);
        prev_bb_num = builder.CreateLoad(T_i32, prev_bb_var);
        builder.CreateStore(ConstantInt::get(T_i32, i), prev_bb_var);
    };
    auto enter_or_pop = [&] (int i) {
        if (i < 0 || bbs[i]->getTerminator()) {
            if (workstack.empty()) {
                pc = end = nullptr;
                return -1;
            }
            i = workstack.pop_back_val();
        }
        enter_bb(i);
        return i;
    };
    enter_bb(0);

    // 5. Start emitting bbs
    auto emit_val = [&] (int id) -> Value* {
        if (id >= 0) {
            return builder.CreateLoad(llvm_ty(func.vals[id]), slots[id]);
        }
        else {
            return emit_const(func.evalConst(id));
        }
    };
    while (pc && end > pc) {
        const auto op = IR::Opcode(*pc);
        pc++;
        const auto res = *pc;
        pc++;
        Value *lres;
        switch (op) {
        case IR::Opcode::Ret: {
            auto val = emit_val(res);
            val = emit_convert(builder, func.ret, val);
            if (func.ret == IR::Type::Bool)
                val = builder.CreateZExt(val, T_i8);
            builder.CreateRet(val);
            enter_or_pop(-1);
            continue;
        }
        case IR::Opcode::Br: {
            auto true_br = *pc;
            pc++;
            if (res == IR::Consts::True) {
                builder.CreateBr(bbs[true_br]);
                enter_or_pop(true_br);
                continue;
            }
            auto false_br = *pc;
            pc++;
            auto val = emit_val(res);
            val = emit_convert(builder, IR::Type::Bool, val);
            builder.CreateCondBr(val, bbs[true_br], bbs[false_br]);
            push_work(true_br);
            enter_or_pop(false_br);
            continue;
        }
        case IR::Opcode::Add:
        case IR::Opcode::Sub:
        case IR::Opcode::Mul:
        case IR::Opcode::FDiv: {
            auto val1 = emit_val(*pc);
            pc++;
            auto val2 = emit_val(*pc);
            pc++;
            switch (op) {
            case IR::Opcode::Add:
                lres = emit_add(builder, func.vals[res], val1, val2);
                break;
            case IR::Opcode::Sub:
                lres = emit_sub(builder, func.vals[res], val1, val2);
                break;
            case IR::Opcode::Mul:
                lres = emit_mul(builder, func.vals[res], val1, val2);
                break;
            case IR::Opcode::FDiv:
                lres = emit_fdiv(builder, val1, val2);
                break;
            default:
                lres = UndefValue::get(llvm_ty(func.vals[res]));
                break;
            }
            break;
        }
        case IR::Opcode::Cmp: {
            auto cmptyp = IR::CmpType(*pc);
            pc++;
            auto val1 = emit_val(*pc);
            pc++;
            auto val2 = emit_val(*pc);
            pc++;
            lres = emit_cmp(builder, cmptyp, val1, val2);
            break;
        }
        case IR::Opcode::Phi: {
            auto nargs = *pc;
            pc++;
            const auto args = pc;
            pc += 2 * nargs;
            if (nargs == 0) {
                lres = UndefValue::get(llvm_ty(func.vals[res]));
                break;
            }
            lres = emit_val(args[1]);
            lres = emit_convert(builder, func.vals[res], lres);
            for (int i = 1;i < nargs;i++) {
                auto cmp = builder.CreateICmpEQ(prev_bb_num,
                                                ConstantInt::get(T_i32, args[2 * i]));
                auto newval = emit_val(args[2 * i + 1]);
                lres = builder.CreateSelect(cmp, emit_convert(builder, func.vals[res],
                                                              newval), lres);
            }
            break;
        }
        case IR::Opcode::Call: {
            auto id = IR::Builtins(*pc);
            pc++;
            auto nargs = *pc;
            pc++;
            const auto args = pc;
            pc += nargs;
            auto emit_intrinsic_f64 = [&] (Intrinsic::ID intrinsic) {
                assert(nargs == 1);
                auto intrin = Intrinsic::getDeclaration(m_mod, intrinsic, {T_f64});
                auto arg = emit_convert(builder, IR::Type::Float64, emit_val(args[0]));
                lres = builder.CreateCall(intrin, arg);
            };
            auto emit_intrinsic_f64f64 = [&] (Intrinsic::ID intrinsic) {
                assert(nargs == 2);
                auto intrin = Intrinsic::getDeclaration(m_mod, intrinsic, {T_f64});
                auto arg1 = emit_convert(builder, IR::Type::Float64, emit_val(args[0]));
                auto arg2 = emit_convert(builder, IR::Type::Float64, emit_val(args[1]));
                lres = builder.CreateCall(intrin, {arg1, arg2});
            };
            auto emit_intrinsic_f64f64f64 = [&] (Intrinsic::ID intrinsic) {
                assert(nargs == 3);
                auto intrin = Intrinsic::getDeclaration(m_mod, intrinsic, {T_f64});
                auto arg1 = emit_convert(builder, IR::Type::Float64, emit_val(args[0]));
                auto arg2 = emit_convert(builder, IR::Type::Float64, emit_val(args[1]));
                auto arg3 = emit_convert(builder, IR::Type::Float64, emit_val(args[2]));
                lres = builder.CreateCall(intrin, {arg1, arg2, arg3});
            };
            switch (id) {
            case IR::Builtins::sqrt:
                emit_intrinsic_f64(Intrinsic::sqrt);
                break;
            case IR::Builtins::sin:
                emit_intrinsic_f64(Intrinsic::sin);
                break;
            case IR::Builtins::cos:
                emit_intrinsic_f64(Intrinsic::cos);
                break;
            case IR::Builtins::pow:
                emit_intrinsic_f64f64(Intrinsic::pow);
                break;
            case IR::Builtins::exp:
                emit_intrinsic_f64(Intrinsic::exp);
                break;
            case IR::Builtins::exp2:
                emit_intrinsic_f64(Intrinsic::exp2);
                break;
            case IR::Builtins::log:
                emit_intrinsic_f64(Intrinsic::log);
                break;
            case IR::Builtins::log2:
                emit_intrinsic_f64(Intrinsic::log2);
                break;
            case IR::Builtins::log10:
                emit_intrinsic_f64(Intrinsic::log10);
                break;
            case IR::Builtins::fma:
                emit_intrinsic_f64f64f64(Intrinsic::fma);
                break;
            case IR::Builtins::abs:
                emit_intrinsic_f64(Intrinsic::fabs);
                break;
            case IR::Builtins::min:
                emit_intrinsic_f64f64(Intrinsic::minnum);
                break;
            case IR::Builtins::max:
                emit_intrinsic_f64f64(Intrinsic::maxnum);
                break;
            case IR::Builtins::copysign:
                emit_intrinsic_f64f64(Intrinsic::copysign);
                break;
            case IR::Builtins::floor:
                emit_intrinsic_f64(Intrinsic::floor);
                break;
            case IR::Builtins::ceil:
                emit_intrinsic_f64(Intrinsic::ceil);
                break;
            case IR::Builtins::rint:
                emit_intrinsic_f64(Intrinsic::rint);
                break;
            case IR::Builtins::round:
                emit_intrinsic_f64(Intrinsic::round);
                break;
            default: {
                auto sym = IR::getBuiltinSymbol(id);
                if (!sym) {
                    lres = UndefValue::get(llvm_ty(func.vals[res]));
                    break;
                }
                switch (IR::getBuiltinType(id)) {
                case IR::BuiltinType::F64_F64: {
                    assert(nargs == 1);
                    auto callee = ensurePureFunc(sym, F_f64_f64);
                    auto arg = emit_convert(builder, IR::Type::Float64, emit_val(args[0]));
                    lres = builder.CreateCall(callee, arg);
                    break;
                }
                case IR::BuiltinType::F64_F64F64: {
                    assert(nargs == 2);
                    auto callee = ensurePureFunc(sym, F_f64_f64f64);
                    auto arg1 = emit_convert(builder, IR::Type::Float64, emit_val(args[0]));
                    auto arg2 = emit_convert(builder, IR::Type::Float64, emit_val(args[1]));
                    lres = builder.CreateCall(callee, {arg1, arg2});
                    break;
                }
                case IR::BuiltinType::F64_F64F64F64: {
                    assert(nargs == 3);
                    auto callee = ensurePureFunc(sym, F_f64_f64f64f64);
                    auto arg1 = emit_convert(builder, IR::Type::Float64, emit_val(args[0]));
                    auto arg2 = emit_convert(builder, IR::Type::Float64, emit_val(args[1]));
                    auto arg3 = emit_convert(builder, IR::Type::Float64, emit_val(args[2]));
                    lres = builder.CreateCall(callee, {arg1, arg2, arg3});
                    break;
                }
                case IR::BuiltinType::F64_F64I32: {
                    assert(nargs == 2);
                    auto callee = ensurePureFunc(sym, F_f64_f64i32);
                    auto arg1 = emit_convert(builder, IR::Type::Float64, emit_val(args[0]));
                    auto arg2 = emit_convert(builder, IR::Type::Int32, emit_val(args[1]));
                    lres = builder.CreateCall(callee, {arg1, arg2});
                    break;
                }
                case IR::BuiltinType::F64_I32F64: {
                    assert(nargs == 2);
                    auto callee = ensurePureFunc(sym, F_f64_i32f64);
                    auto arg1 = emit_convert(builder, IR::Type::Int32, emit_val(args[0]));
                    auto arg2 = emit_convert(builder, IR::Type::Float64, emit_val(args[1]));
                    lres = builder.CreateCall(callee, {arg1, arg2});
                    break;
                }
                default:
                    lres = UndefValue::get(llvm_ty(func.vals[res]));
                }
            }
                break;
            }
            break;
        }
        case IR::Opcode::Interp: {
            auto input = emit_convert(builder, IR::Type::Float64, emit_val(*pc));
            pc++;
            auto x0 = emit_convert(builder, IR::Type::Float64, emit_val(*pc));
            pc++;
            auto dx = emit_convert(builder, IR::Type::Float64, emit_val(*pc));
            pc++;
            auto data_offset = *pc;
            pc++;
            auto ndata = *pc;
            pc++;

            ArrayRef<double> dataref(&func.float_table[data_offset], ndata);
            auto table = ConstantDataArray::get(m_ctx, dataref);
            Constant *datap = new GlobalVariable(*m_mod, table->getType(), true,
                                                 GlobalValue::InternalLinkage, table,
                                                 ".L.nacs." + std::to_string(m_counter++));
            cast<GlobalVariable>(datap)->setUnnamedAddr(GlobalValue::UnnamedAddr::Global);
            datap = ConstantExpr::getBitCast(datap, T_f64->getPointerTo());

            auto interp_f = ensurePureFunc("interp", F_f64_f64f64f64i32pf64, true);
            lres = builder.CreateCall(interp_f, {input, x0, dx,
                        ConstantInt::get(T_i32, ndata), datap});
            break;
        }
        case IR::Opcode::Convert: {
            lres = emit_val(*pc);
            pc++;
            break;
        }
        case IR::Opcode::Select: {
            auto cond = emit_convert(builder, IR::Type::Bool, emit_val(*pc));
            pc++;
            auto ty = func.vals[res];
            auto v1 = emit_convert(builder, ty, emit_val(*pc));
            pc++;
            auto v2 = emit_convert(builder, ty, emit_val(*pc));
            pc++;
            lres = builder.CreateSelect(cond, v1, v2);
            break;
        }
        case IR::Opcode::And:
        case IR::Opcode::Or:
        case IR::Opcode::Xor: {
            auto v1 = emit_convert(builder, IR::Type::Bool, emit_val(*pc));
            pc++;
            auto v2 = emit_convert(builder, IR::Type::Bool, emit_val(*pc));
            pc++;
            if (op == IR::Opcode::And) {
                lres = builder.CreateAnd(v1, v2);
            }
            else if (op == IR::Opcode::Or) {
                lres = builder.CreateOr(v1, v2);
            }
            else if (op == IR::Opcode::Xor) {
                lres = builder.CreateXor(v1, v2);
            }
            break;
        }
        case IR::Opcode::Not: {
            auto v = emit_convert(builder, IR::Type::Bool, emit_val(*pc));
            pc++;
            lres = builder.CreateNot(v);
            break;
        }
        default:
            lres = UndefValue::get(llvm_ty(func.vals[res]));
            break;
        }
        builder.CreateStore(emit_convert(builder, func.vals[res], lres), slots[res]);
    }

    // 6. Fix up unreachable bbs
    for (auto bb: bbs) {
        if (!bb->getTerminator()) {
            builder.SetInsertPoint(bb);
            builder.CreateUnreachable();
        }
    }
    return f;
}

}
}
}
