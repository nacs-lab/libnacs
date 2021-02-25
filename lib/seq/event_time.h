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

#ifndef __NACS_SEQ_EVENT_TIME_H__
#define __NACS_SEQ_EVENT_TIME_H__

#include "env.h"

namespace NaCs::Seq {

struct EventTime {
    enum Errno : uint8_t {
        NegTime,
        NonPosTime
    };
    enum Sign : uint8_t {
        Unknown = 0,
        NonNeg,
        Pos
    };
    struct Term {
        // This layout is not ideal for memory consumption
        // but makes it much easier to sort without c++20
        Sign sign;
        Var::Ref var;
        // This ID is only used for error reporting.
        uint32_t id = 0;
        bool operator==(const Term &other) const
        {
            return sign == other.sign && var.get() == other.var.get();
        }
        bool operator!=(const Term &other) const
        {
            return !operator==(other);
        }
    };

    EventTime() = default;
    EventTime(uint64_t tconst, std::vector<Term> &&terms={})
        : tconst(tconst),
          terms(std::move(terms))
    {}
    EventTime(EventTime&&) = default;
    EventTime(const EventTime&);
    bool operator==(const EventTime &other) const
    {
        return tconst == other.tconst && terms == other.terms;
    }
    bool operator!=(const EventTime &other) const
    {
        return !operator==(other);
    }

    void add_term(Sign sign, Var *var, uint32_t id=0)
    {
        terms.push_back({sign, var->ref(), id});
    }

    bool normalize();
    Var *to_var(Env &env) const;
    void print(std::ostream &stm, bool newline=false) const;

    // Compare the terms.
    // If `this` is known to be less than `other`, return `Pos`.
    // If `this` is known to be no greater than `other`, return `NonNeg`.
    // Otherwise, return `Unknown`.
    Sign isless_terms(const EventTime &other) const;
    Sign isless(const EventTime &other) const;
    // other must be known to be smaller than/before us
    // in a way that's consistent with `isless_terms`
    EventTime operator-(const EventTime &other) const;

    uint64_t tconst = 0;
    // I'd like to use `llvm::SmallVector` here
    // but I also don't want to expose any LLVM linkage directly
    // (in order to support statically linking LLVM).
    std::vector<Term> terms;
};

static inline std::ostream &operator<<(std::ostream &stm, const EventTime &t)
{
    t.print(stm);
    return stm;
}

}

#endif
