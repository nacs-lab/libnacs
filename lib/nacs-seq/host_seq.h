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

#ifndef __NACS_SEQ_HOST_SEQ_H__
#define __NACS_SEQ_HOST_SEQ_H__

#include "sign.h"

#include <vector>

namespace NaCs::Seq {

/**
 * This hosts information for finalizing the sequence (all the variables) at runtime.
 * This may not contain all the information for the backends.
 */
class HostSeq {
public:
    enum class Type {
        Bool = 1,
        Int32 = 2,
        Float64 = 3,
        Int64 = 4, // Used only for times
    };
    union Value {
        bool b;
        int32_t i32;
        double f64;
        int64_t i64;
    };
    struct Pulse {
        uint32_t id;

        uint32_t time;
        uint32_t measure; // may be -1
        uint32_t len; // may be -1 or -2, `len == -1` <=> no-ramp pulse, `len == -2` <=> measure
        // For ramps, this is a dummy slot tracking the dependency of the ramp function.
        uint32_t value; // -1 for measures
        uint32_t chn; // 1-based indexing
        double (*ramp_func)(double, void*);

        // This value is generally useful for all backends.
        uint32_t endvalue;
        uint32_t cond; // may be -1

        bool is_measure() const
        {
            return len == uint32_t(-2);
        }
    };
    struct Assignment {
        uint32_t global_id;
        uint32_t value;
    };
    struct Assumption {
        Sign sign;
        uint32_t value;
        uint32_t id;
    };
    struct Branch {
        uint32_t cond;
        uint32_t target; // -1 means exit
    };
    struct BasicSeq {
        uint32_t id;

        /**
         * Values info
         */
        // Value counts
        // The basic sequence specific values are also stored after the shared ones
        // in this order in the `values` array.
        uint32_t nmeasure; // including start values of pulses, all measures have type Float64
        // Ones that can directly be computed
        uint32_t ndirect;
        // Ones that has to evaluate by figuring out the pulse order.
        uint32_t nneed_order;

        // Types, no need for measures since those are all `Float64`
        std::vector<Type> types; // `.size() == ndirect + nneed_order`
        // Functions to populate each slots
        // (again, no need for measures since we have to fill those in ourselves)
        std::vector<void (*)(void*)> evals; // `.size() == ndirect + nneed_order`

        /**
         * Values dependency info
         */
        // Global values that are used by this basic sequence,
        // used to minimize computation. In range [0, nglobal_vals).
        // We don't need to recoard what global variables (in range [0, nglobals))
        // are used since those are supplied by the user and we never need to compute those.
        std::vector<uint32_t> global_refs;
        // The values (in [0, nneed_order)) that depends on each measures
        std::vector<std::vector<uint32_t>> reverse_depends; // `.size() == nmeasure`
        // Number of measures each value depends on.
        // Combined with `reverse_depends` we can figure out the set of values
        // that can be evaluated when a measure is fullfilled.
        std::vector<uint32_t> deps_count; // `.size() == nneed_order`

        /**
         * Sequence info
         */
        mutable std::vector<Pulse> pulses;
        std::vector<Assignment> assignments;
        std::vector<Assumption> assumptions;
        std::vector<Branch> branches;
        uint32_t default_branch; // -1 means exit
        std::vector<uint32_t> endtimes;

        /**
         * Assumption dependency info
         */
        // Number of assumes that doesn't depend on pulse ordering
        // This corresponds to `ndirect` values
        uint32_t ndirect_assumes;
        // The index of the corresponding assumption from measure or computed values
        // `-1` means there's no corresponding assumptions
        // Since all assumptions has type `Int64` there's no overlap with measures
        std::vector<uint32_t> assumptions_idx; // `.size() == nneed_order`

        /**
         * Condition dependency info
         */
        std::vector<uint32_t> cond_global_refs;

        // Output
        mutable uint64_t length;

    private:
        void init() const;
        void init_run() const;

        // Workspace variables
        // These are the main ones that are mutated in each run
        // (`pulses` is being resorted but the set of pulses doesn't change).
        mutable std::vector<bool> m_measure_filled; // `.size() == nmeasure`
        mutable std::vector<uint32_t> m_deps_count_copy; // `.size() == nneed_order`
        friend class HostSeq;
    };

    // Call before all runs
    void init();

    // Call before starting this sequence
    void init_run();
    // Call before running a basic sequence within this sequence
    void pre_run();
    // Call after running a basic sequence within this sequence
    // Returns the user supplied ID of the next basic sequence.
    uint32_t post_run();

    double get_global(uint32_t i) const;
    void set_global(uint32_t i, double val);
    uint32_t cur_seq_idx() const
    {
        // For backends
        return m_cur_seq_idx;
    }
    const std::vector<uint64_t> &global_ages() const
    {
        return m_ages;
    }
    double get_value(uint32_t i, uint32_t seq_idx) const
    {
        return get_value(seqs[seq_idx], i);
    }
    Type get_type(uint32_t i, uint32_t seq_idx) const;
    int64_t get_time(uint32_t i) const
    {
        return values[i].i64;
    }

    // Public states (for construction of sequence)

    /**
     * Value info
     */
    // Value counts
    // The shared values are stored at the beginning of the `values` array in this order
    uint32_t nconsts;
    // User read/write-able global variables
    uint32_t nglobals;
    uint32_t nglobal_vals;
    // Current values on each channels
    uint32_t nchannels; // corresponding values must be `double`
    uint32_t nshared; // `== nconsts + nglobals + nglobal_vals + nchannels`

    uint32_t npublic_globals;
    bool first_bseq;

    // The values here are valid after `pre_run` and `post_run`
    // for the current run. Automatic update for the next sequence,
    // e.g. updating of the start values for each channel must wait until
    // before the next `pre_run` so that backends can read the correct numbers.
    std::vector<Value> values;
    // Backend should look at `start_values` instead.
    std::vector<Value> default_values;
    // Start values corrected for any t=0 pulses
    std::vector<Value> start_values;
    std::vector<Type> types; // `.size() == nconst + nglobals + nglobal_vals`

    // Only includes dependency on global variables.
    // Dependency to other computed global values is encoded in the order of the values
    // in the `values` array.
    std::vector<std::vector<uint32_t>> depends; // `.size() == nglobal_vals`
    std::vector<void (*)(void*)> global_evals; // `.size() == nglobal_vals`

    std::vector<BasicSeq> seqs;

private:
    static double convert_value(Type type, Value value);
    static bool assign_value(Type type, Value &slot, double val);
    // Input is index of value in the range
    // [0, nshared + bseq.nmeasure + bseq.ndirect + bseq.nneed_order)
    bool value_computed(const BasicSeq &bseq, uint32_t i) const;
    double get_value(const BasicSeq &bseq, uint32_t i) const;

    /**
     * Functions to iteratively evaluate all the values
     */
    // Input is index of global values in the range [0, nglobal_vals)
    // This should be called in increasing order of `i` each time
    // and all the values that the value depend on must be updated first.
    void update_global_value(uint32_t i);
    // Input is index of measure in the range [0, bseq.nmeasures)
    // Returns whether any time value has been evaluated.
    // The caller may need to resort the pulses when that happens.
    bool propagate_measure(const BasicSeq &bseq, uint32_t i);
    void check_assumption(const BasicSeq &bseq, const Assumption &assumption) const;
    void update_values(const BasicSeq &bseq);
    void _set_global(uint32_t i, double val);

    uint32_t eval_branch(const BasicSeq &bseq); // -1 means exit

    // Global age tracking
    uint64_t m_cur_age;
    std::vector<uint64_t> m_ages; // `.size() == nglobals + nglobal_vals`
    bool m_global_dirty;
    bool m_first_bseq = false;
    bool m_initialized = false;
    bool m_running = false;

    // Sequence tracking
    uint32_t m_cur_seq_idx = -1;

    // Workspace variables
    // In additional to be the buffer we use to compute the values,
    // this holds the values of each channel before and after each basic sequence.
    std::vector<double> m_chn_vals;
    std::vector<const Pulse*> m_chn_pulses;
};

}

#endif
