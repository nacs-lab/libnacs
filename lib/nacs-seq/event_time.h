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

#include "sign.h"
#include "env.h"

namespace NaCs::Seq {

class BasicSeq;

class EventTime {
    struct Unrefer {
        void operator()(EventTime *time) const
        {
            time->unref();
        }
    };
public:
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
    // Use unique pointer since I'm too lazy to write my own.
    // The only limitation is that copy is not allowed
    // but we could/should use `->ref()` most of the time in place of copy anyway.
    class Ref : public std::unique_ptr<EventTime,Unrefer> {
    public:
        // Unlike `Var::Ref`, we do not allow keeping the `Ref` alive
        // past the lifetime of the allocator of `EventTime` itself.
        // This should be fine since all the `Ref`'s should be stored in the sequence
        // that owns the `EventTime`.
        // This way the use of the reference counting to manage the lifetime of `EventTime`
        // is optional and we can still easily use our own `EventTime` however we want
        // when constructing the sequence.
        using std::unique_ptr<EventTime,Unrefer>::unique_ptr;
        void reset(EventTime *ptr) noexcept
        {
            if (ptr)
                ptr->incref();
            std::unique_ptr<EventTime,Unrefer>::reset(ptr);
        }
        EventTime *release() noexcept = delete;
    };

    Ref ref()
    {
        incref();
        return Ref(this);
    }

    EventTime() = default;
    EventTime(int64_t tconst, std::vector<Term> &&terms={})
        : tconst(tconst),
          terms(std::move(terms))
    {}
    // We want pointer identity and this will be kept track of by the parent (BasicSeq)
    // so we don't want or need moving.
    EventTime(EventTime&&) = delete;
    explicit EventTime(const EventTime&); // We don't want to accidentally copy either
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
    int32_t order_id() const
    {
        return m_order_id;
    }

    bool normalize();
    Var *to_var(Env &env) const;
    void print(std::ostream &stm, bool newline=false) const;

    // Compare the terms.
    // If `this` is known to be less than `other`, return `Sign::Pos`.
    // If `this` is known to be no greater than `other`, return `Sign::NonNeg`.
    // Otherwise, return `Sign::Unknown`.
    Sign isless_terms(const EventTime &other) const;
    Sign isless(const EventTime &other) const;
    // other must be known to be smaller than/before us
    // in a way that's consistent with `isless_terms`
    EventTime operator-(const EventTime &other) const;

    int64_t min_const() const;

    int64_t tconst = 0;
    // I'd like to use `llvm::SmallVector` here
    // but I also don't want to expose any LLVM linkage directly
    // (in order to support statically linking LLVM).
    std::vector<Term> terms;

private:
    void incref()
    {
        m_ref_count++;
    }
    void unref()
    {
        assert(m_ref_count > 0);
        --m_ref_count;
    }
    uint32_t m_ref_count = 0;
    int32_t m_order_id = -1;
    int32_t m_scc_index = 0;
    int32_t m_scc_index_lowest = 0;
    bool m_on_stack = false;
    friend class BasicSeq;
};

static inline std::ostream &operator<<(std::ostream &stm, const EventTime &t)
{
    t.print(stm);
    return stm;
}

}

#endif
