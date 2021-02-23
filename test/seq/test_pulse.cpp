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

#include "../../lib/seq/pulse.h"
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
    t.add_term(Seq::EventTime::Pos, env.new_const(IR::TagVal(1.45)));
    t.add_term(Seq::EventTime::Pos, env.new_extern({IR::Type::Float64, 1234}));
    t.add_term(Seq::EventTime::Pos, env.new_const(IR::TagVal(5.45)));
    t.add_term(Seq::EventTime::NonNeg, env.new_extern({IR::Type::Float64, 2345}));

    {
        Seq::Pulse p1(1, t, env.new_const(IR::TagVal(4.5)),
                      env.new_const(IR::TagVal(8.9)), false);
        assert(p1.id() == 1);
        assert(!p1.is_measure());
        assert(p1.start() == t);
        assert(p1.len());
        assert(p1.val());
        assert(!p1.needs_oldval());
        auto ev = p1.endval();
        assert(ev);
        assert(ev->is_const());
        assert(ev->get_const().is(IR::TagVal(8.9)));
        assert(p1.len());

        p1.clear_unused_args();
        p1.optimize();

        assert(!p1.len());
        assert(p1.val());
        assert(p1.endval() == ev);
    }

    {
        auto v = [&] {
            IR::Builder builder(IR::Type::Float64, {IR::Type::Float64, IR::Type::Float64});
            builder.createRet(builder.createAdd(0, 1));
            return env.new_call(builder.get(), {Seq::Arg::create_const(0.5),
                    Seq::Arg::create_arg(1)}, 2);
        }();
        Seq::Pulse p1(1, t, env.new_const(IR::TagVal(4.5)), v, false);
        assert(p1.id() == 1);
        assert(!p1.is_measure());
        assert(p1.start() == t);
        assert(p1.len());
        assert(p1.val());
        assert(p1.needs_oldval());
        auto ev = p1.endval();
        assert(!ev); // Old value needed but unspecified.
        assert(p1.len());

        p1.set_oldval(env.new_extern({IR::Type::Float64, 1234}));

        ev = p1.endval();
        assert(ev);
        // set_oldval should also check for unused time arguments
        // and optimize out len if possible.
        assert(!p1.len());

        Seq::Pulse p2(2, t, nullptr, v, false);

        assert(p1.start() == p2.start());
        assert(p1.start() == t);
        assert(p1.known_before(p2) == Seq::EventTime::Pos);
        assert(p2.known_before(p1) == Seq::EventTime::Unknown);

        p1.optimize();
        p2.optimize();

        assert(p1.known_before(p2) == Seq::EventTime::Pos);
        assert(p2.known_before(p1) == Seq::EventTime::Unknown);
    }

    {
        Seq::Pulse m(1, t, nullptr,
                     env.new_extern({IR::Type::Float64, 1234}), true);
        assert(m.is_measure());
    }

    {
        Seq::Pulse p(1, t, env.new_const(IR::TagVal(4.5)), [&] {
            IR::Builder builder(IR::Type::Float64,
                                {IR::Type::Float64, IR::Type::Float64});
            builder.createRet(
                builder.createAdd(1, builder.getConst(0.0008)));
            return env.new_call(builder.get(), {Seq::Arg::create_arg(0),
                    Seq::Arg::create_arg(1)}, 2);
        }(), false);
        assert(p.len());
        env.optimize();
        p.clear_unused_args();
        assert(p.id() == 1);
        assert(p.needs_oldval());
        assert(!p.endval());
        assert(!p.len());
    }

    {
        // No error for negative time offset error
        Seq::EventTime t(-10);
        t.add_term(Seq::EventTime::NonNeg, env.new_extern({IR::Type::Float64, 1234}), 49);
        Seq::Pulse p(1, t, nullptr, env.new_const(IR::TagVal(4.5)), false);
        p.optimize();
        p.check_start();

        Seq::EventTime t2(-10);
        Seq::Pulse p2(13, t2, nullptr, env.new_const(IR::TagVal(4.5)), false);
        auto err = expect_error<Seq::Error>([&] {
            p2.check_start();
        });
        assert(err.code >> 16 == Seq::Error::Pulse);
        assert((err.code & ((1u << 16) - 1)) == Seq::Pulse::NegTime);
        assert(err.type1 == Seq::Error::Pulse);
        assert(err.id1 == 13);
        assert(strcmp(err.what(), "Pulse time must be positive.") == 0);
    }

    return 0;
}
