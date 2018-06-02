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

namespace NaCs {
namespace LLVM {
namespace Codegen {

Context::Context(Module *mod)
    : m_mod(mod),
      m_ctx(mod->getContext()),
      T_bool(Type::getInt1Ty(m_ctx)),
      T_i32(Type::getInt32Ty(m_ctx)),
      T_isz(sizeof(void*) == 8 ? Type::getInt64Ty(m_ctx) : T_i32),
      T_f64(Type::getDoubleTy(m_ctx)),
      F_f64_f64(FunctionType::get(T_f64, {T_f64}, false)),
      F_f64_f64f64(FunctionType::get(T_f64, {T_f64, T_f64}, false)),
      F_f64_f64f64f64(FunctionType::get(T_f64, {T_f64, T_f64, T_f64}, false)),
      F_f64_f64i32(FunctionType::get(T_f64, {T_f64, T_i32}, false)),
      F_f64_i32f64(FunctionType::get(T_f64, {T_i32, T_f64}, false))
{
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

}
}
}
