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

#define LLVM_DISABLE_ABI_BREAKING_CHECKS_ENFORCING 1

#include "../../lib/seq/seq.h"
#include "../../lib/seq/error.h"

#include "../../lib/utils/llvm/codegen.h"
#include "../../lib/utils/llvm/utils.h"

#include <llvm/IR/LLVMContext.h>

#include <assert.h>
#include <iostream>
#include <memory>

using namespace NaCs;

static void test_channels(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn") == 0);
    assert(seq.get_chn_id("test_chn", true) == 1);
    assert(seq.get_chn_id("test_chn", true) == 1);
    assert(seq.get_chn_id("test_chn") == 1);
    assert(seq.get_chn_name(1) == "test_chn");
    assert(seq.get_chn_id("test_chn2", true) == 2);
    assert(seq.get_chn_name(1) == "test_chn");
    assert(seq.get_chn_name(2) == "test_chn2");
    assert(seq.get_chn_names().size() == 2);
    assert(seq.get_chn_names()[0] == "test_chn");
    assert(seq.get_chn_names()[1] == "test_chn2");
}

static void test_default(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    assert(seq.get_chn_id("test_chn", true) == 1);
    auto defval = seq.defval(1);
    assert(defval == 0);
    // Reuse of constant variables
    assert(seq.get_const(IR::TagVal(2.3)) == seq.get_const(IR::TagVal(2.3)));
    assert(seq.get_const(IR::TagVal(2.3)) != seq.get_const(IR::TagVal(2.0)));
    seq.set_defval(1, 2.3);
    defval = seq.defval(1);
    assert(defval == 2.3);
}

static void test_global(llvm::LLVMContext &llvm_ctx)
{
    Seq::Seq seq(llvm_ctx);
    auto slot1 = seq.get_slot(IR::Type::Float64, 0);
    auto slot2 = seq.get_slot(IR::Type::Bool, 1);
    auto slot3 = seq.get_slot(IR::Type::Int32, 2);
    assert(slot1 == seq.get_slot(IR::Type::Bool, 0));
    assert(slot2 == seq.get_slot(IR::Type::Int32, 1));
    assert(slot3 == seq.get_slot(IR::Type::Int32, 2));
    assert(slot1->type() == IR::Type::Float64);
    assert(slot2->type() == IR::Type::Bool);
    assert(slot3->type() == IR::Type::Int32);
    auto bs = seq.add_basicseq(1);
    // This corresponds to swapping the values
    // All the assignment uses the value of the variable before the assignment happens.
    bs->assign_global(0, slot2, 41);
    bs->assign_global(1, slot1, 52);
    bs->assign_global(2, seq.get_const(IR::TagVal(3.4)), 13);
    auto &assigns = bs->get_assigns();
    assert(assigns.size() == 3);
    assert(assigns.find(0)->second.val.get() == slot2);
    assert(assigns.find(0)->second.id == 41);
    assert(assigns.find(1)->second.val.get() == slot1);
    assert(assigns.find(1)->second.id == 52);
    // get_const provides pointer identity to Float64 constants.
    assert(assigns.find(2)->second.val.get() == seq.get_const(IR::TagVal(3.4)));
    assert(assigns.find(2)->second.id == 13);
}

int main()
{
    auto llvm_ctx = LLVM::new_context();

    test_channels(*llvm_ctx);
    test_default(*llvm_ctx);
    test_global(*llvm_ctx);

    return 0;
}
