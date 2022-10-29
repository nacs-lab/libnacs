/*************************************************************************
 *   Copyright (c) 2021 - 2022 Yichao Yu <yyc1992@gmail.com>             *
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
#define CATCH_CONFIG_MAIN

#include "../error_helper.h"

#include "../../lib/nacs-seq/event_time.h"
#include "../../lib/nacs-seq/error.h"
#include "../../lib/nacs-utils/llvm/utils.h"

#include <llvm/IR/LLVMContext.h>

#include <iostream>

using namespace NaCs;
using namespace std::literals::string_literals;

static Seq::Env env = [] {
    static auto llvm_ctx = LLVM::new_context();
    return Seq::Env(*llvm_ctx);
}();

TEST_CASE("EventTime") {
    Seq::EventTime t;
    t.add_term(Seq::Sign::Pos, env.new_const(IR::TagVal(1.2)), 1);
    t.normalize();
    REQUIRE(t.tconst == 1);
    REQUIRE(t.terms.empty());
    t.add_term(Seq::Sign::Pos, env.new_extern({IR::Type::Float64, 1234}), 2);
    t.normalize();
    REQUIRE(t.tconst == 1);
    REQUIRE(t.terms.size() == 1);
    REQUIRE(t.terms[0].var->is_extern());
    REQUIRE(t.min_const() == 1);
    t.tconst += 5;
    REQUIRE(t.min_const() == 6);

    auto tv = env.new_extern({IR::Type::Float64, 1111});
    Seq::EventTime t2;
    t2.tconst = 10;
    REQUIRE(t2.min_const() == 10);
    t.add_term(Seq::Sign::Unknown, tv, 3);
    REQUIRE(t.min_const() == 0);
    t2.add_term(Seq::Sign::Unknown, tv, 3);
    REQUIRE(t2.min_const() == 0);
    t.normalize();
    t2.normalize();
    // We sort by varid, so `tv` should remain as the second term.
    REQUIRE(t.terms.size() == 2);
    REQUIRE(t.terms[1].var.get() == tv);
    REQUIRE(t2.terms.size() == 1);
    REQUIRE(t2.terms[0].var.get() == tv);
    REQUIRE(t.isless_terms(t2) == Seq::Sign::Unknown);
    REQUIRE(t2.isless_terms(t) == Seq::Sign::Pos);
    REQUIRE(t.isless(t2) == Seq::Sign::Unknown);
    REQUIRE(t2.isless(t) == Seq::Sign::Unknown);

    t2.tconst = 1;
    REQUIRE(t.isless_terms(t2) == Seq::Sign::Unknown);
    REQUIRE(t2.isless_terms(t) == Seq::Sign::Pos);
    REQUIRE(t.isless(t2) == Seq::Sign::Unknown);
    REQUIRE(t2.isless(t) == Seq::Sign::Pos);

    t2.tconst = 6;
    REQUIRE(t.isless_terms(t2) == Seq::Sign::Unknown);
    REQUIRE(t2.isless_terms(t) == Seq::Sign::Pos);
    REQUIRE(t.isless(t2) == Seq::Sign::Unknown);
    REQUIRE(t2.isless(t) == Seq::Sign::Pos);

    REQUIRE(t.terms[0].var.get() != tv);
    REQUIRE(t.terms[0].sign == Seq::Sign::Pos);
    t.terms[0].sign = Seq::Sign::NonNeg;
    REQUIRE(t.isless_terms(t2) == Seq::Sign::Unknown);
    REQUIRE(t2.isless_terms(t) == Seq::Sign::NonNeg);
    REQUIRE(t.isless(t2) == Seq::Sign::Unknown);
    REQUIRE(t2.isless(t) == Seq::Sign::NonNeg);

    t2.tconst = 4;
    REQUIRE(t.isless_terms(t2) == Seq::Sign::Unknown);
    REQUIRE(t2.isless_terms(t) == Seq::Sign::NonNeg);
    REQUIRE(t.isless(t2) == Seq::Sign::Unknown);
    REQUIRE(t2.isless(t) == Seq::Sign::Pos);

    auto tdiff = t - t2;
    REQUIRE(tdiff.tconst == 2);
    REQUIRE(tdiff.terms.size() == 1);
    REQUIRE(tdiff.terms[0].sign == Seq::Sign::NonNeg);
    REQUIRE(tdiff.terms[0].var == t.terms[0].var);

    SECTION("Copy") {
        Seq::EventTime t3(t);
        REQUIRE(t == t3);
        // Make sure id is also copied.
        // `==` does not check this since it is insignificant to the sequence semantics
        auto nterms = t.terms.size();
        for (size_t i = 0; i < nterms; i++) {
            REQUIRE(t.terms[i].id == t3.terms[i].id);
            REQUIRE(t.terms[i].id != 0);
        }
    }

    SECTION("Compare extra terms") {
        auto term1 = env.new_extern({IR::Type::Float64, 10});
        auto term2 = env.new_extern({IR::Type::Float64, 20});

        Seq::EventTime t;
        t.add_term(Seq::Sign::Pos, term1, 1);
        Seq::EventTime t2;
        t2.add_term(Seq::Sign::Pos, term1, 1);
        t2.add_term(Seq::Sign::Unknown, term2, 2);
        Seq::EventTime t3;
        t3.add_term(Seq::Sign::Unknown, term2, 2);
        t3.add_term(Seq::Sign::Pos, term1, 1);

        REQUIRE(t.isless_terms(t2) == Seq::Sign::Unknown);
        REQUIRE(t2.isless_terms(t) == Seq::Sign::Unknown);

        REQUIRE(t.isless_terms(t3) == Seq::Sign::Unknown);
        REQUIRE(t3.isless_terms(t) == Seq::Sign::Unknown);
    }
}

TEST_CASE("Error") {
    SECTION("Pos") {
        Seq::EventTime t;
        t.add_term(Seq::Sign::Pos, env.new_const(IR::TagVal(0)), 42);
        auto err = expect_error<Seq::Error>([&] {
            t.normalize();
        });
        REQUIRE(err.type == Seq::Error::Type::EventTime);
        REQUIRE(err.code == uint16_t(Seq::Error::EventTime::NonPosTime));
        REQUIRE(err.type1 == Seq::Error::Type::EventTime);
        REQUIRE(err.id1 == 42);
        REQUIRE(err.what() == "Positive time expected."s);
    }

    SECTION("NonNeg") {
        Seq::EventTime t;
        t.add_term(Seq::Sign::NonNeg, env.new_const(IR::TagVal(-2)), 49);
        auto err = expect_error<Seq::Error>([&] {
            t.normalize();
        });
        REQUIRE(err.type == Seq::Error::Type::EventTime);
        REQUIRE(err.code == uint16_t(Seq::Error::EventTime::NegTime));
        REQUIRE(err.type1 == Seq::Error::Type::EventTime);
        REQUIRE(err.id1 == 49);
        REQUIRE(err.what() == "Non-negative time expected."s);
    }
}
