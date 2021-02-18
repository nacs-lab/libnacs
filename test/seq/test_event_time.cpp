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

#include "../error_helper.h"

#include "../../lib/seq/event_time.h"
#include "../../lib/seq/error.h"
#include "../../lib/utils/llvm/utils.h"

#include <llvm/IR/LLVMContext.h>

#include <assert.h>
#include <iostream>

using namespace NaCs;

int main()
{
    auto llvm_ctx = LLVM::new_context();
    Seq::Env env(*llvm_ctx);

    Seq::EventTime t;
    t.add_term(Seq::EventTime::Pos, env.new_const(IR::TagVal(1.2)), 1);
    t.normalize();
    assert(t.tconst == 1);
    assert(t.terms.empty());
    t.add_term(Seq::EventTime::Pos, env.new_extern({IR::Type::Float64, 1234}), 2);
    t.normalize();
    assert(t.tconst == 1);
    assert(t.terms.size() == 1);
    assert(t.terms[0].var->is_extern());
    t.tconst += 5;

    auto tv = env.new_extern({IR::Type::Float64, 1111});
    Seq::EventTime t2;
    t2.tconst = 10;
    t.add_term(Seq::EventTime::Unknown, tv, 3);
    t2.add_term(Seq::EventTime::Unknown, tv, 3);
    t.normalize();
    t2.normalize();
    // We sort by varid, so `tv` should remain as the second term.
    assert(t.terms.size() == 2);
    assert(t.terms[1].var.get() == tv);
    assert(t2.terms.size() == 1);
    assert(t2.terms[0].var.get() == tv);
    assert(t.isless_terms(t2) == Seq::EventTime::Unknown);
    assert(t2.isless_terms(t) == Seq::EventTime::Pos);
    assert(t.isless(t2) == Seq::EventTime::Unknown);
    assert(t2.isless(t) == Seq::EventTime::Unknown);

    t2.tconst = 1;
    assert(t.isless_terms(t2) == Seq::EventTime::Unknown);
    assert(t2.isless_terms(t) == Seq::EventTime::Pos);
    assert(t.isless(t2) == Seq::EventTime::Unknown);
    assert(t2.isless(t) == Seq::EventTime::Pos);

    t2.tconst = 6;
    assert(t.isless_terms(t2) == Seq::EventTime::Unknown);
    assert(t2.isless_terms(t) == Seq::EventTime::Pos);
    assert(t.isless(t2) == Seq::EventTime::Unknown);
    assert(t2.isless(t) == Seq::EventTime::Pos);

    assert(t.terms[0].var.get() != tv);
    assert(t.terms[0].sign == Seq::EventTime::Pos);
    t.terms[0].sign = Seq::EventTime::NonNeg;
    assert(t.isless_terms(t2) == Seq::EventTime::Unknown);
    assert(t2.isless_terms(t) == Seq::EventTime::NonNeg);
    assert(t.isless(t2) == Seq::EventTime::Unknown);
    assert(t2.isless(t) == Seq::EventTime::NonNeg);

    t2.tconst = 4;
    assert(t.isless_terms(t2) == Seq::EventTime::Unknown);
    assert(t2.isless_terms(t) == Seq::EventTime::NonNeg);
    assert(t.isless(t2) == Seq::EventTime::Unknown);
    assert(t2.isless(t) == Seq::EventTime::Pos);

    auto tdiff = t - t2;
    assert(tdiff.tconst == 2);
    assert(tdiff.terms.size() == 1);
    assert(tdiff.terms[0].sign == Seq::EventTime::NonNeg);
    assert(tdiff.terms[0].var == t.terms[0].var);

    {
        // Test copy constructor
        Seq::EventTime t3(t);
        assert(t == t3);
        // Make sure id is also copied.
        // `==` does not check this since it is insignificant to the sequence semantics
        auto nterms = t.terms.size();
        for (size_t i = 0; i < nterms; i++) {
            assert(t.terms[i].id == t3.terms[i].id);
            assert(t.terms[i].id != 0);
        }
    }

    {
        // Test error
        Seq::EventTime t;
        t.add_term(Seq::EventTime::Pos, env.new_const(IR::TagVal(0)), 42);
        auto err = expect_error<Seq::Error>([&] {
            t.normalize();
        });
        assert(err.code >> 16 == Seq::Error::EventTime);
        assert((err.code & ((1u << 16) - 1)) == Seq::EventTime::NonPosTime);
        assert(err.type1 == Seq::Error::EventTime);
        assert(err.id1 == 42);
        assert(strcmp(err.what(), "Positive time expected.") == 0);
    }

    {
        // Test error
        Seq::EventTime t;
        t.add_term(Seq::EventTime::NonNeg, env.new_const(IR::TagVal(-2)), 49);
        auto err = expect_error<Seq::Error>([&] {
            t.normalize();
        });
        assert(err.code >> 16 == Seq::Error::EventTime);
        assert((err.code & ((1u << 16) - 1)) == Seq::EventTime::NegTime);
        assert(err.type1 == Seq::Error::EventTime);
        assert(err.id1 == 49);
        assert(strcmp(err.what(), "Non-negative time expected.") == 0);
    }

    return 0;
}
