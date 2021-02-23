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
#define CATCH_CONFIG_MAIN

#include "../error_helper.h"

#include "../../lib/nacs-seq/pulse.h"
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

TEST_CASE("Pulse") {
    Seq::EventTime t;
    t.add_term(Seq::Sign::Pos, env.new_const(IR::TagVal(1.45)));
    t.add_term(Seq::Sign::Pos, env.new_extern({IR::Type::Float64, 1234}));
    t.add_term(Seq::Sign::Pos, env.new_const(IR::TagVal(5.45)));
    t.add_term(Seq::Sign::NonNeg, env.new_extern({IR::Type::Float64, 2345}));

    {
        Seq::Pulse p1(1, t, env.new_const(IR::TagVal(4.5)),
                      env.new_const(IR::TagVal(8.9)), nullptr, false);
        REQUIRE(p1.id() == 1);
        REQUIRE(!p1.is_measure());
        REQUIRE(p1.start() == t);
        REQUIRE(p1.len());
        REQUIRE(p1.val());
        REQUIRE(!p1.needs_oldval());
        auto ev = p1.endval();
        REQUIRE(ev);
        REQUIRE(ev->is_const());
        REQUIRE(ev->get_const().is(IR::TagVal(8.9)));
        REQUIRE(p1.len());

        p1.clear_unused_args();
        p1.optimize();

        REQUIRE(!p1.len());
        REQUIRE(p1.val());
        REQUIRE(p1.endval() == ev);
    }

    {
        auto v = [&] {
            IR::Builder builder(IR::Type::Float64, {IR::Type::Float64, IR::Type::Float64});
            builder.createRet(builder.createAdd(0, 1));
            return env.new_call(builder.get(), {Seq::Arg::create_const(0.5),
                    Seq::Arg::create_arg(1)}, 2);
        }();
        Seq::Pulse p1(1, t, env.new_const(IR::TagVal(4.5)), v, nullptr, false);
        REQUIRE(p1.id() == 1);
        REQUIRE(!p1.is_measure());
        REQUIRE(p1.start() == t);
        REQUIRE(p1.len());
        REQUIRE(p1.val());
        REQUIRE(p1.needs_oldval());
        auto ev = p1.endval();
        REQUIRE(!ev); // Old value needed but unspecified.
        REQUIRE(p1.len());

        p1.set_oldval(env.new_extern({IR::Type::Float64, 1234}));

        ev = p1.endval();
        REQUIRE(ev);
        // set_oldval should also check for unused time arguments
        // and optimize out len if possible.
        REQUIRE(!p1.len());

        Seq::Pulse p2(2, t, nullptr, v, nullptr, false);

        REQUIRE(p1.start() == p2.start());
        REQUIRE(p1.start() == t);

        p1.optimize();
        p2.optimize();
    }

    {
        Seq::Pulse m(1, t, nullptr,
                     env.new_extern({IR::Type::Float64, 1234}), nullptr, true);
        REQUIRE(m.is_measure());
    }

    {
        Seq::Pulse p(1, t, env.new_const(IR::TagVal(4.5)), [&] {
            IR::Builder builder(IR::Type::Float64,
                                {IR::Type::Float64, IR::Type::Float64});
            builder.createRet(
                builder.createAdd(1, builder.getConst(0.0008)));
            return env.new_call(builder.get(), {Seq::Arg::create_arg(0),
                    Seq::Arg::create_arg(1)}, 2);
        }(), nullptr, false);
        REQUIRE(p.len());
        env.optimize();
        p.clear_unused_args();
        REQUIRE(p.id() == 1);
        REQUIRE(p.needs_oldval());
        REQUIRE(!p.endval());
        REQUIRE(!p.len());
    }

    {
        // No error for negative time offset error
        Seq::EventTime t(-10);
        t.add_term(Seq::Sign::NonNeg, env.new_extern({IR::Type::Float64, 1234}), 49);
        Seq::Pulse p(1, t, nullptr, env.new_const(IR::TagVal(4.5)), nullptr, false);
        p.optimize();
        p.check_start();

        Seq::EventTime t2(-10);
        Seq::Pulse p2(13, t2, nullptr, env.new_const(IR::TagVal(4.5)), nullptr, false);
        auto err = expect_error<Seq::Error>([&] {
            p2.check_start();
        });
        REQUIRE(err.type == Seq::Error::Type::Pulse);
        REQUIRE(err.code == uint16_t(Seq::Error::Pulse::NegTime));
        REQUIRE(err.type1 == Seq::Error::Type::Pulse);
        REQUIRE(err.id1 == 13);
        REQUIRE(err.what() == "Pulse time must be positive."s);
    }

    // Conditional pulse + old value handling
    {
        Seq::Pulse p1(1, t, env.new_const(IR::TagVal(4.5)),
                      env.new_const(IR::TagVal(8.9)),
                      env.new_extern({IR::Type::Bool, 4234}), false);
        REQUIRE(p1.id() == 1);
        REQUIRE(!p1.is_measure());
        REQUIRE(p1.start() == t);
        REQUIRE(p1.len());
        REQUIRE(p1.val());
        REQUIRE(p1.val()->is_const());
        REQUIRE(!p1.needs_oldval());
        auto ev = p1.endval();
        REQUIRE(ev);
        REQUIRE(ev->is_const());
        REQUIRE(p1.len());

        p1.clear_unused_args();
        p1.optimize();

        REQUIRE(!p1.len());
    }
}
