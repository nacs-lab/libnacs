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

#define CATCH_CONFIG_MAIN

#include "../../lib/nacs-seq/host_seq.h"
#include "../../lib/nacs-seq/error.h"

#include "../../lib/nacs-utils/number.h"

#include <iostream>
#include <memory>

#include <catch2/catch.hpp>

using namespace NaCs;

TEST_CASE("default_empty") {
    Seq::HostSeq host_seq;
    host_seq.nconsts = 0;
    host_seq.nglobals = 0;
    host_seq.nglobal_vals = 0;
    host_seq.nchannels = 2;
    host_seq.nshared = (host_seq.nconsts + host_seq.nglobals +
                        host_seq.nglobal_vals + host_seq.nchannels);
    host_seq.npublic_globals = 0;

    host_seq.values.resize(2);
    host_seq.default_values.resize(2);
    host_seq.default_values[0].f64 = 1.2;
    host_seq.default_values[1].f64 = 3.4;
    host_seq.types.resize(0); // nconst + nglobals + nglobal_vals

    host_seq.depends.resize(0); // nglobal_vals
    host_seq.global_evals.resize(0); // nglobal_vals

    host_seq.seqs.resize(1);

    {
        auto &host_bseq = host_seq.seqs[0];
        host_bseq.id = 1;

        host_bseq.nmeasure = 0;
        host_bseq.ndirect = 0;
        host_bseq.nneed_order = 0;

        host_bseq.types.resize(0); // ndirect + nneed_order
        host_bseq.evals.resize(0); // ndirect + nneed_order

        host_bseq.global_refs.resize(0);
        host_bseq.reverse_depends.resize(0); // nmeasure
        host_bseq.deps_count.resize(0); // nneed_order

        host_bseq.pulses.resize(0);
        host_bseq.assignments.resize(0);
        host_bseq.assumptions.resize(0);
        host_bseq.branches.resize(0);
        host_bseq.default_branch = -1;
        host_bseq.endtimes.resize(0);

        host_bseq.ndirect_assumes = 0;
        host_bseq.assumptions_idx.resize(0); // nneed_order

        host_bseq.cond_global_refs.resize(0);
    }

    host_seq.init();
    REQUIRE(host_seq.global_ages().size() == 0);

    for (uint32_t i = 0; i < 3; i++) {
        REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
        host_seq.init_run();
        REQUIRE(host_seq.cur_seq_idx() == 0);
        host_seq.pre_run();
        REQUIRE(host_seq.cur_seq_idx() == 0);
        REQUIRE(host_seq.first_bseq);

        REQUIRE(host_seq.seqs[0].length == 0);

        REQUIRE(host_seq.start_values[0].f64 == 1.2);
        REQUIRE(host_seq.start_values[1].f64 == 3.4);

        REQUIRE(host_seq.values[0].f64 == 1.2);
        REQUIRE(host_seq.values[1].f64 == 3.4);
        REQUIRE(host_seq.post_run() == 0);
        REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
    }
}

TEST_CASE("increment_loop") {
    // Ch1: increase by 0.2 per cycle
    // Ch2: decrease by 0.3 per cycle
    // Branch condition: stop when Ch1 > Ch2
    Seq::HostSeq host_seq;
    host_seq.nconsts = 1;
    host_seq.nglobals = 2;
    host_seq.nglobal_vals = 1;
    host_seq.nchannels = 2;
    host_seq.nshared = (host_seq.nconsts + host_seq.nglobals +
                        host_seq.nglobal_vals + host_seq.nchannels);
    host_seq.npublic_globals = 0;

    host_seq.values.resize(8);
    host_seq.values[0].i64 = 125;
    host_seq.default_values.resize(2);
    host_seq.default_values[0].f64 = 1.2;
    host_seq.default_values[1].f64 = 3.4;
    host_seq.types.resize(4); // nconst + nglobals + nglobal_vals
    host_seq.types[0] = Seq::HostSeq::Type::Int64;
    host_seq.types[1] = Seq::HostSeq::Type::Float64;
    host_seq.types[2] = Seq::HostSeq::Type::Float64;
    host_seq.types[3] = Seq::HostSeq::Type::Bool;

    host_seq.depends.resize(1); // nglobal_vals
    host_seq.depends[0].push_back(0);
    host_seq.depends[0].push_back(1);
    host_seq.global_evals.resize(1); // nglobal_vals
    host_seq.global_evals[0] = [] (void *_p) {
        auto p = (Seq::HostSeq::Value*)_p;
        p[3].b = p[1].f64 <= p[2].f64;
    };

    host_seq.seqs.resize(1);

    {
        auto &host_bseq = host_seq.seqs[0];
        host_bseq.id = 1;

        host_bseq.nmeasure = 0;
        host_bseq.ndirect = 2;
        host_bseq.nneed_order = 0;

        host_bseq.types.resize(2); // ndirect + nneed_order
        host_bseq.types[0] = Seq::HostSeq::Type::Float64;
        host_bseq.types[1] = Seq::HostSeq::Type::Float64;
        host_bseq.evals.resize(2); // ndirect + nneed_order
        host_bseq.evals[0] = [] (void *_p) {
            auto p = (Seq::HostSeq::Value*)_p;
            p[6].f64 = p[4].f64 + 0.2;
        };
        host_bseq.evals[1] = [] (void *_p) {
            auto p = (Seq::HostSeq::Value*)_p;
            p[7].f64 = p[5].f64 - 0.3;
        };

        host_bseq.global_refs.resize(0);
        host_bseq.reverse_depends.resize(0); // nmeasure
        host_bseq.deps_count.resize(0); // nneed_order

        host_bseq.pulses.resize(2);
        host_bseq.pulses[0] = Seq::HostSeq::Pulse{
            .id = 1, .time = 0, .measure = uint32_t(-1), .len = uint32_t(-1),
            .value = 6, .chn = 1, .ramp_func = nullptr, .endvalue = 6, .cond = uint32_t(-1)
        };
        host_bseq.pulses[1] = Seq::HostSeq::Pulse{
            .id = 2, .time = 0, .measure = uint32_t(-1), .len = uint32_t(-1),
            .value = 7, .chn = 2, .ramp_func = nullptr, .endvalue = 7, .cond = uint32_t(-1)
        };
        host_bseq.assignments.resize(2);
        host_bseq.assignments[0] = Seq::HostSeq::Assignment{ .global_id = 0, .value = 6 };
        host_bseq.assignments[1] = Seq::HostSeq::Assignment{ .global_id = 1, .value = 7 };
        host_bseq.assumptions.resize(0);
        host_bseq.branches.resize(1);
        host_bseq.branches[0] = Seq::HostSeq::Branch{ .cond = 3, .target = 0 };
        host_bseq.default_branch = -1;
        host_bseq.endtimes.resize(1);
        host_bseq.endtimes[0] = 0;

        host_bseq.ndirect_assumes = 0;
        host_bseq.assumptions_idx.resize(0); // nneed_order

        host_bseq.cond_global_refs.resize(1);
        host_bseq.cond_global_refs[0] = 0;
    }

    host_seq.init();
    auto &ages = host_seq.global_ages();
    REQUIRE(ages.size() == 3);
    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 1);
    REQUIRE(ages[2] == 0);

    uint32_t age = 1;

    for (uint32_t i = 0; i < 3; i++) {
        REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
        host_seq.init_run();
        REQUIRE(host_seq.cur_seq_idx() == uint32_t(0));

        double v1 = 1.2;
        double v2 = 3.4;

        for (int j = 0; ; j++) {
            REQUIRE(host_seq.cur_seq_idx() == 0);
            host_seq.pre_run();

            REQUIRE(host_seq.start_values[0].f64 == v1);
            REQUIRE(host_seq.start_values[1].f64 == v2);
            REQUIRE(host_seq.first_bseq == (j == 0));
            REQUIRE(host_seq.cur_seq_idx() == 0);

            auto oldv1 = v1;
            auto oldv2 = v2;

            v1 += 0.2;
            v2 -= 0.3;

            // Assigned globals
            REQUIRE(host_seq.values[1].f64 == v1);
            REQUIRE(host_seq.values[2].f64 == v2);
            // Channel initial values
            REQUIRE(host_seq.values[4].f64 == oldv1);
            REQUIRE(host_seq.values[5].f64 == oldv2);
            // Pulse values
            REQUIRE(host_seq.values[6].f64 == v1);
            REQUIRE(host_seq.values[7].f64 == v2);

            auto host_bseq = &host_seq.seqs[0];
            REQUIRE(host_bseq->pulses[0].id == 1);
            REQUIRE(host_bseq->pulses[0].chn == 1);
            REQUIRE(host_seq.values[host_bseq->pulses[0].endvalue].f64 == v1);
            REQUIRE(host_bseq->pulses[1].id == 2);
            REQUIRE(host_bseq->pulses[1].chn == 2);
            REQUIRE(host_seq.values[host_bseq->pulses[1].endvalue].f64 == v2);

            REQUIRE(host_bseq->length == 125);

            age++;
            REQUIRE(ages[0] == age);
            REQUIRE(ages[1] == age);
            if (age == 2) {
                REQUIRE(ages[2] == 0);
            }
            else {
                REQUIRE(ages[2] == age - 1);
            }
            bool cond = v1 <= v2;
            REQUIRE(host_seq.post_run() == (cond ? 1 : 0));
            REQUIRE(host_seq.values[3].b == cond);

            REQUIRE(ages[0] == age);
            REQUIRE(ages[1] == age);
            REQUIRE(ages[2] == age);

            if (!cond) {
                REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
                break;
            }
            REQUIRE(host_seq.cur_seq_idx() == 0);
        }
    }
}

TEST_CASE("increment_loop_t0") {
    // Ch1: increase by 0.2 per cycle
    // Ch2: decrease by 0.3 per cycle
    // Branch condition: stop when Ch1 > Ch2
    Seq::HostSeq host_seq;
    host_seq.nconsts = 1;
    host_seq.nglobals = 2;
    host_seq.nglobal_vals = 1;
    host_seq.nchannels = 2;
    host_seq.nshared = (host_seq.nconsts + host_seq.nglobals +
                        host_seq.nglobal_vals + host_seq.nchannels);
    host_seq.npublic_globals = 0;

    host_seq.values.resize(8);
    host_seq.values[0].i64 = 0;
    host_seq.default_values.resize(2);
    host_seq.default_values[0].f64 = 1.2;
    host_seq.default_values[1].f64 = 3.4;
    host_seq.types.resize(4); // nconst + nglobals + nglobal_vals
    host_seq.types[0] = Seq::HostSeq::Type::Int64;
    host_seq.types[1] = Seq::HostSeq::Type::Float64;
    host_seq.types[2] = Seq::HostSeq::Type::Float64;
    host_seq.types[3] = Seq::HostSeq::Type::Bool;

    host_seq.depends.resize(1); // nglobal_vals
    host_seq.depends[0].push_back(0);
    host_seq.depends[0].push_back(1);
    host_seq.global_evals.resize(1); // nglobal_vals
    host_seq.global_evals[0] = [] (void *_p) {
        auto p = (Seq::HostSeq::Value*)_p;
        p[3].b = p[1].f64 <= p[2].f64;
    };

    host_seq.seqs.resize(1);

    {
        auto &host_bseq = host_seq.seqs[0];
        host_bseq.id = 1;

        host_bseq.nmeasure = 0;
        host_bseq.ndirect = 2;
        host_bseq.nneed_order = 0;

        host_bseq.types.resize(2); // ndirect + nneed_order
        host_bseq.types[0] = Seq::HostSeq::Type::Float64;
        host_bseq.types[1] = Seq::HostSeq::Type::Float64;
        host_bseq.evals.resize(2); // ndirect + nneed_order
        host_bseq.evals[0] = [] (void *_p) {
            auto p = (Seq::HostSeq::Value*)_p;
            p[6].f64 = p[4].f64 + 0.2;
        };
        host_bseq.evals[1] = [] (void *_p) {
            auto p = (Seq::HostSeq::Value*)_p;
            p[7].f64 = p[5].f64 - 0.3;
        };

        host_bseq.global_refs.resize(0);
        host_bseq.reverse_depends.resize(0); // nmeasure
        host_bseq.deps_count.resize(0); // nneed_order

        host_bseq.pulses.resize(2);
        host_bseq.pulses[0] = Seq::HostSeq::Pulse{
            .id = 1, .time = 0, .measure = uint32_t(-1), .len = uint32_t(-1),
            .value = 6, .chn = 1, .ramp_func = nullptr, .endvalue = 6, .cond = uint32_t(-1)
        };
        host_bseq.pulses[1] = Seq::HostSeq::Pulse{
            .id = 2, .time = 0, .measure = uint32_t(-1), .len = uint32_t(-1),
            .value = 7, .chn = 2, .ramp_func = nullptr, .endvalue = 7, .cond = uint32_t(-1)
        };
        host_bseq.assignments.resize(2);
        host_bseq.assignments[0] = Seq::HostSeq::Assignment{ .global_id = 0, .value = 6 };
        host_bseq.assignments[1] = Seq::HostSeq::Assignment{ .global_id = 1, .value = 7 };
        host_bseq.assumptions.resize(0);
        host_bseq.branches.resize(1);
        host_bseq.branches[0] = Seq::HostSeq::Branch{ .cond = 3, .target = 0 };
        host_bseq.default_branch = -1;
        host_bseq.endtimes.resize(1);
        host_bseq.endtimes[0] = 0;

        host_bseq.ndirect_assumes = 0;
        host_bseq.assumptions_idx.resize(0); // nneed_order

        host_bseq.cond_global_refs.resize(1);
        host_bseq.cond_global_refs[0] = 0;
    }

    host_seq.init();
    auto &ages = host_seq.global_ages();
    REQUIRE(ages.size() == 3);
    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 1);
    REQUIRE(ages[2] == 0);

    uint32_t age = 1;

    for (uint32_t i = 0; i < 3; i++) {
        REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
        host_seq.init_run();
        REQUIRE(host_seq.cur_seq_idx() == uint32_t(0));

        double v1 = 1.2;
        double v2 = 3.4;

        for (int j = 0; ; j++) {
            REQUIRE(host_seq.cur_seq_idx() == 0);
            host_seq.pre_run();

            auto oldv1 = v1;
            auto oldv2 = v2;

            v1 += 0.2;
            v2 -= 0.3;

            // t=0 pulses will cause the start values to be updated
            REQUIRE(host_seq.start_values[0].f64 == v1);
            REQUIRE(host_seq.start_values[1].f64 == v2);
            REQUIRE(host_seq.first_bseq == (j == 0));
            REQUIRE(host_seq.cur_seq_idx() == 0);

            // Assigned globals
            REQUIRE(host_seq.values[1].f64 == v1);
            REQUIRE(host_seq.values[2].f64 == v2);
            // Channel final values
            REQUIRE(host_seq.values[4].f64 == oldv1);
            REQUIRE(host_seq.values[5].f64 == oldv2);
            // Pulse values
            REQUIRE(host_seq.values[6].f64 == v1);
            REQUIRE(host_seq.values[7].f64 == v2);

            auto host_bseq = &host_seq.seqs[0];
            REQUIRE(host_bseq->pulses[0].id == 1);
            REQUIRE(host_bseq->pulses[0].chn == 1);
            REQUIRE(host_seq.values[host_bseq->pulses[0].endvalue].f64 == v1);
            REQUIRE(host_bseq->pulses[1].id == 2);
            REQUIRE(host_bseq->pulses[1].chn == 2);
            REQUIRE(host_seq.values[host_bseq->pulses[1].endvalue].f64 == v2);

            REQUIRE(host_bseq->length == 0);

            age++;
            REQUIRE(ages[0] == age);
            REQUIRE(ages[1] == age);
            if (age == 2) {
                REQUIRE(ages[2] == 0);
            }
            else {
                REQUIRE(ages[2] == age - 1);
            }
            bool cond = v1 <= v2;
            REQUIRE(host_seq.post_run() == (cond ? 1 : 0));
            REQUIRE(host_seq.values[3].b == cond);

            REQUIRE(ages[0] == age);
            REQUIRE(ages[1] == age);
            REQUIRE(ages[2] == age);

            if (!cond) {
                REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
                break;
            }
            REQUIRE(host_seq.cur_seq_idx() == 0);
        }
    }
}

TEST_CASE("global") {
    // Const:
    //   C0[V0]: time
    // Slot:
    //   S0[V1]: set before BS1
    //   S1[V2]: set before BS2
    //   S2[V3]: stashed value of CH1, set by BS1/BS2
    // Global Value:
    //   G0[V4]: S0 + 0.5
    //   G1[V5]: S0 + S1
    //   G2[V6]: S1 > S2
    //   G3[V7]: S0 > S2
    // Channel:
    //   CV1[V8]: Channel 1 value
    // BS1:
    //   Values:
    //     D0[V9]: CV1 + 0.2 + G0
    //   Assignments:
    //     S2 <- D0
    //   Branch:
    //     G2: BS1
    //     Default: BS2
    // BS2:
    //   Values:
    //     D0[V9]: CV1 + 3.4 - G1
    //   Assignments:
    //     S2 <- D0
    //   Branch:
    //     G3: BS2
    //     Default: end
    Seq::HostSeq host_seq;
    host_seq.nconsts = 1;
    host_seq.nglobals = 3;
    host_seq.nglobal_vals = 4;
    host_seq.nchannels = 1;
    host_seq.nshared = (host_seq.nconsts + host_seq.nglobals +
                        host_seq.nglobal_vals + host_seq.nchannels);
    host_seq.npublic_globals = 2;

    host_seq.values.resize(10);
    host_seq.values[0].i64 = 123;
    host_seq.default_values.resize(1);
    host_seq.default_values[0].f64 = 1.2;
    host_seq.types.resize(9); // nconst + nglobals + nglobal_vals
    host_seq.types[0] = Seq::HostSeq::Type::Float64;
    host_seq.types[1] = Seq::HostSeq::Type::Float64;
    host_seq.types[2] = Seq::HostSeq::Type::Float64;
    host_seq.types[3] = Seq::HostSeq::Type::Float64;
    host_seq.types[4] = Seq::HostSeq::Type::Float64;
    host_seq.types[5] = Seq::HostSeq::Type::Float64;
    host_seq.types[6] = Seq::HostSeq::Type::Bool; // Condition 1
    host_seq.types[7] = Seq::HostSeq::Type::Bool; // Condition 2

    host_seq.depends.resize(4); // nglobal_vals
    host_seq.depends[0].push_back(0);
    host_seq.depends[1].push_back(0);
    host_seq.depends[1].push_back(1);
    host_seq.depends[2].push_back(1);
    host_seq.depends[2].push_back(2);
    host_seq.depends[3].push_back(0);
    host_seq.depends[3].push_back(2);
    host_seq.global_evals.resize(4); // nglobal_vals
    host_seq.global_evals[0] = [] (void *_p) {
        auto p = (Seq::HostSeq::Value*)_p;
        p[4].f64 = p[1].f64 + 0.5;
    };
    host_seq.global_evals[1] = [] (void *_p) {
        auto p = (Seq::HostSeq::Value*)_p;
        p[5].f64 = p[1].f64 + p[2].f64;
    };
    host_seq.global_evals[2] = [] (void *_p) {
        auto p = (Seq::HostSeq::Value*)_p;
        p[6].b = p[2].f64 > p[3].f64;
    };
    host_seq.global_evals[3] = [] (void *_p) {
        auto p = (Seq::HostSeq::Value*)_p;
        p[7].b = p[1].f64 > p[3].f64;
    };

    host_seq.seqs.resize(2);

    {
        auto &host_bseq = host_seq.seqs[0];
        host_bseq.id = 1;

        host_bseq.nmeasure = 0;
        host_bseq.ndirect = 1;
        host_bseq.nneed_order = 0;

        host_bseq.types.resize(1); // ndirect + nneed_order
        host_bseq.types[0] = Seq::HostSeq::Type::Float64;
        host_bseq.evals.resize(1); // ndirect + nneed_order
        host_bseq.evals[0] = [] (void *_p) {
            auto p = (Seq::HostSeq::Value*)_p;
            p[9].f64 = p[8].f64 + 0.2 + p[4].f64;
        };

        host_bseq.global_refs.resize(1);
        host_bseq.global_refs[0] = 0;
        host_bseq.reverse_depends.resize(0); // nmeasure
        host_bseq.deps_count.resize(0); // nneed_order

        host_bseq.pulses.resize(1);
        host_bseq.pulses[0] = Seq::HostSeq::Pulse{
            .id = 1, .time = 0, .measure = uint32_t(-1), .len = uint32_t(-1),
            .value = 9, .chn = 1, .ramp_func = nullptr, .endvalue = 9, .cond = uint32_t(-1)
        };
        host_bseq.assignments.resize(1);
        host_bseq.assignments[0] = Seq::HostSeq::Assignment{ .global_id = 2, .value = 9 };
        host_bseq.assumptions.resize(0);
        host_bseq.branches.resize(1);
        host_bseq.branches[0] = Seq::HostSeq::Branch{ .cond = 6, .target = 0 };
        host_bseq.default_branch = 1;
        host_bseq.endtimes.resize(1);
        host_bseq.endtimes[0] = 0;

        host_bseq.ndirect_assumes = 0;
        host_bseq.assumptions_idx.resize(0); // nneed_order

        host_bseq.cond_global_refs.resize(1);
        host_bseq.cond_global_refs[0] = 2;
    }

    {
        auto &host_bseq = host_seq.seqs[1];
        host_bseq.id = 2;

        host_bseq.nmeasure = 0;
        host_bseq.ndirect = 1;
        host_bseq.nneed_order = 0;

        host_bseq.types.resize(1); // ndirect + nneed_order
        host_bseq.types[0] = Seq::HostSeq::Type::Float64;
        host_bseq.evals.resize(1); // ndirect + nneed_order
        host_bseq.evals[0] = [] (void *_p) {
            auto p = (Seq::HostSeq::Value*)_p;
            p[9].f64 = p[8].f64 + 3.4 - p[5].f64;
        };

        host_bseq.global_refs.resize(1);
        host_bseq.global_refs[0] = 1;
        host_bseq.reverse_depends.resize(0); // nmeasure
        host_bseq.deps_count.resize(0); // nneed_order

        host_bseq.pulses.resize(1);
        host_bseq.pulses[0] = Seq::HostSeq::Pulse{
            .id = 2, .time = 0, .measure = uint32_t(-1), .len = uint32_t(-1),
            .value = 9, .chn = 1, .ramp_func = nullptr, .endvalue = 9, .cond = uint32_t(-1)
        };
        host_bseq.assignments.resize(1);
        host_bseq.assignments[0] = Seq::HostSeq::Assignment{ .global_id = 2, .value = 9 };
        host_bseq.assumptions.resize(0);
        host_bseq.branches.resize(1);
        host_bseq.branches[0] = Seq::HostSeq::Branch{ .cond = 7, .target = 1 };
        host_bseq.default_branch = -1;
        host_bseq.endtimes.resize(1);
        host_bseq.endtimes[0] = 0;

        host_bseq.ndirect_assumes = 0;
        host_bseq.assumptions_idx.resize(0); // nneed_order

        host_bseq.cond_global_refs.resize(1);
        host_bseq.cond_global_refs[0] = 3;
    }

    host_seq.init();
    auto &ages = host_seq.global_ages();
    REQUIRE(ages.size() == 7);
    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 1);
    REQUIRE(ages[2] == 1);
    REQUIRE(ages[3] == 0);
    REQUIRE(ages[4] == 0);
    REQUIRE(ages[5] == 0);
    REQUIRE(ages[6] == 0);

    // Set and get global values
    REQUIRE(host_seq.values[1].f64 == 0);
    REQUIRE(host_seq.values[2].f64 == 0);
    REQUIRE(host_seq.values[3].f64 == 0);
    REQUIRE(host_seq.get_global(0) == 0);
    REQUIRE(host_seq.get_global(1) == 0);
    REQUIRE(host_seq.get_global(2) == 0);

    host_seq.set_global(0, 3.4);
    host_seq.set_global(1, 5.6);
    host_seq.set_global(2, 6.9);

    REQUIRE(host_seq.values[1].f64 == 3.4);
    REQUIRE(host_seq.values[2].f64 == 5.6);
    REQUIRE(host_seq.values[3].f64 == 0);
    REQUIRE(host_seq.get_global(0) == 3.4);
    REQUIRE(host_seq.get_global(1) == 5.6);
    REQUIRE(host_seq.get_global(2) == 0);

    host_seq.values[3].f64 = 9.8;
    REQUIRE(host_seq.values[3].f64 == 9.8);
    REQUIRE(host_seq.get_global(2) == 0);

    host_seq.values[3].f64 = 0;

    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 1);
    REQUIRE(ages[2] == 1);
    REQUIRE(ages[3] == 0);
    REQUIRE(ages[4] == 0);
    REQUIRE(ages[5] == 0);
    REQUIRE(ages[6] == 0);

    /**
     * Run 1
     */
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
    host_seq.init_run();
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(0));

    /**
     * Run 1 SubRun 1: BSeq 1
     */
    REQUIRE(host_seq.cur_seq_idx() == 0);
    host_seq.pre_run();
    auto host_bseq = &host_seq.seqs[0];
    REQUIRE(host_bseq->length == 123);
    REQUIRE(host_seq.start_values[0].f64 == 1.2);

    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 1);
    REQUIRE(ages[2] == 2);
    REQUIRE(ages[3] == 1);
    REQUIRE(ages[4] == 0);
    REQUIRE(ages[5] == 0);
    REQUIRE(ages[6] == 0);

    REQUIRE(host_seq.values[0].i64 == 123);
    REQUIRE(host_seq.values[1].f64 == 3.4);
    REQUIRE(host_seq.values[2].f64 == 5.6);
    REQUIRE(host_seq.values[3].f64 == Approx(5.3));
    REQUIRE(host_seq.values[4].f64 == Approx(3.9)); // 3.4 + 0.5
    // REQUIRE(host_seq.values[5].f64 == 0); // Unused by BS1
    // REQUIRE(host_seq.values[6].b == false); // Unused by BS1
    // REQUIRE(host_seq.values[7].b == false); // Unused by BS1
    REQUIRE(host_seq.values[8].f64 == 1.2);
    REQUIRE(host_seq.values[9].f64 == Approx(5.3)); // 1.2 + 0.2 + 3.9

    REQUIRE(host_bseq->pulses.size() == 1);
    REQUIRE(host_bseq->pulses[0].id == 1);
    REQUIRE(host_bseq->pulses[0].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[0].endvalue].f64 == Approx(5.3));

    REQUIRE(host_seq.post_run() == 1);
    REQUIRE(host_seq.cur_seq_idx() == 0);

    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 1);
    REQUIRE(ages[2] == 2);
    REQUIRE(ages[3] == 1);
    REQUIRE(ages[4] == 0);
    REQUIRE(ages[5] == 2);
    REQUIRE(ages[6] == 0);

    REQUIRE(host_seq.values[0].i64 == 123);
    REQUIRE(host_seq.values[1].f64 == 3.4);
    REQUIRE(host_seq.values[2].f64 == 5.6);
    REQUIRE(host_seq.values[3].f64 == Approx(5.3));
    REQUIRE(host_seq.values[4].f64 == Approx(3.9));
    // REQUIRE(host_seq.values[5].f64 == 0); // Unused by BS1
    REQUIRE(host_seq.values[6].b == true);
    // REQUIRE(host_seq.values[7].b == false); // Unused by BS1
    REQUIRE(host_seq.values[8].f64 == 1.2);
    REQUIRE(host_seq.values[9].f64 == Approx(5.3));

    /**
     * Run 1 SubRun 2: BSeq 1
     */
    host_seq.set_global(0, -0.3);

    REQUIRE(ages[0] == 3);
    REQUIRE(ages[1] == 1);
    REQUIRE(ages[2] == 2);
    REQUIRE(ages[3] == 1);
    REQUIRE(ages[4] == 0);
    REQUIRE(ages[5] == 2);
    REQUIRE(ages[6] == 0);

    REQUIRE(host_seq.cur_seq_idx() == 0);
    host_seq.pre_run();
    host_bseq = &host_seq.seqs[0];
    REQUIRE(host_bseq->length == 123);
    REQUIRE(host_seq.start_values[0].f64 == Approx(5.3));

    REQUIRE(ages[0] == 3);
    REQUIRE(ages[1] == 1);
    REQUIRE(ages[2] == 4);
    REQUIRE(ages[3] == 3);
    REQUIRE(ages[4] == 0);
    REQUIRE(ages[5] == 2);
    REQUIRE(ages[6] == 0);

    REQUIRE(host_seq.values[0].i64 == 123);
    REQUIRE(host_seq.values[1].f64 == -0.3);
    REQUIRE(host_seq.values[2].f64 == 5.6);
    REQUIRE(host_seq.values[3].f64 == Approx(5.7));
    REQUIRE(host_seq.values[4].f64 == Approx(0.2)); // -0.3 + 0.5
    // REQUIRE(host_seq.values[5].f64 == 0); // Unused by BS1
    REQUIRE(host_seq.values[6].b == true);
    // REQUIRE(host_seq.values[7].b == false); // Unused by BS1
    REQUIRE(host_seq.values[8].f64 == Approx(5.3));
    REQUIRE(host_seq.values[9].f64 == Approx(5.7)); // 5.3 + 0.2 + 0.2

    REQUIRE(host_bseq->pulses.size() == 1);
    REQUIRE(host_bseq->pulses[0].id == 1);
    REQUIRE(host_bseq->pulses[0].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[0].endvalue].f64 == Approx(5.7));

    host_seq.set_global(1, 5.9);

    REQUIRE(ages[0] == 3);
    REQUIRE(ages[1] == 4);
    REQUIRE(ages[2] == 4);
    REQUIRE(ages[3] == 3);
    REQUIRE(ages[4] == 0);
    REQUIRE(ages[5] == 2);
    REQUIRE(ages[6] == 0);

    REQUIRE(host_seq.values[0].i64 == 123);
    REQUIRE(host_seq.values[1].f64 == -0.3);
    REQUIRE(host_seq.values[2].f64 == 5.9);
    REQUIRE(host_seq.values[3].f64 == Approx(5.7));
    REQUIRE(host_seq.values[4].f64 == Approx(0.2));
    // REQUIRE(host_seq.values[5].f64 == 0); // Unused by BS1
    REQUIRE(host_seq.values[6].b == true);
    // REQUIRE(host_seq.values[7].b == false); // Unused by BS1
    REQUIRE(host_seq.values[8].f64 == Approx(5.3));
    REQUIRE(host_seq.values[9].f64 == Approx(5.7));

    REQUIRE(host_seq.post_run() == 1);
    REQUIRE(host_seq.cur_seq_idx() == 0);

    REQUIRE(ages[0] == 3);
    REQUIRE(ages[1] == 4);
    REQUIRE(ages[2] == 4);
    REQUIRE(ages[3] == 3);
    REQUIRE(ages[4] == 0);
    REQUIRE(ages[5] == 4);
    REQUIRE(ages[6] == 0);

    REQUIRE(host_seq.values[0].i64 == 123);
    REQUIRE(host_seq.values[1].f64 == -0.3);
    REQUIRE(host_seq.values[2].f64 == 5.9);
    REQUIRE(host_seq.values[3].f64 == Approx(5.7));
    REQUIRE(host_seq.values[4].f64 == Approx(0.2));
    // REQUIRE(host_seq.values[5].f64 == 0); // Unused by BS1
    REQUIRE(host_seq.values[6].b == true);
    // REQUIRE(host_seq.values[7].b == false); // Unused by BS1
    REQUIRE(host_seq.values[8].f64 == Approx(5.3));
    REQUIRE(host_seq.values[9].f64 == Approx(5.7));

    // Setting global variables to the same value should not increase age
    host_seq.set_global(0, -0.3);
    host_seq.set_global(1, 5.9);

    REQUIRE(ages[0] == 3);
    REQUIRE(ages[1] == 4);
    REQUIRE(ages[2] == 4);
    REQUIRE(ages[3] == 3);
    REQUIRE(ages[4] == 0);
    REQUIRE(ages[5] == 4);
    REQUIRE(ages[6] == 0);

    REQUIRE(host_seq.values[0].i64 == 123);
    REQUIRE(host_seq.values[1].f64 == -0.3);
    REQUIRE(host_seq.values[2].f64 == 5.9);
    REQUIRE(host_seq.values[3].f64 == Approx(5.7));
    REQUIRE(host_seq.values[4].f64 == Approx(0.2));
    // REQUIRE(host_seq.values[5].f64 == 0); // Unused by BS1
    REQUIRE(host_seq.values[6].b == true);
    // REQUIRE(host_seq.values[7].b == false); // Unused by BS1
    REQUIRE(host_seq.values[8].f64 == Approx(5.3));
    REQUIRE(host_seq.values[9].f64 == Approx(5.7));

    /**
     * Run 1 SubRun 3: BSeq 1
     */
    REQUIRE(host_seq.cur_seq_idx() == 0);
    host_seq.pre_run();
    host_bseq = &host_seq.seqs[0];
    REQUIRE(host_bseq->length == 123);
    REQUIRE(host_seq.start_values[0].f64 == Approx(5.7));

    REQUIRE(ages[0] == 3);
    REQUIRE(ages[1] == 4);
    REQUIRE(ages[2] == 5);
    REQUIRE(ages[3] == 3);
    REQUIRE(ages[4] == 0);
    REQUIRE(ages[5] == 4);
    REQUIRE(ages[6] == 0);

    REQUIRE(host_seq.values[0].i64 == 123);
    REQUIRE(host_seq.values[1].f64 == -0.3);
    REQUIRE(host_seq.values[2].f64 == 5.9);
    REQUIRE(host_seq.values[3].f64 == Approx(6.1));
    REQUIRE(host_seq.values[4].f64 == Approx(0.2));
    // REQUIRE(host_seq.values[5].f64 == 0); // Unused by BS1
    REQUIRE(host_seq.values[6].b == true);
    // REQUIRE(host_seq.values[7].b == false); // Unused by BS1
    REQUIRE(host_seq.values[8].f64 == Approx(5.7));
    REQUIRE(host_seq.values[9].f64 == Approx(6.1)); // 5.7 + 0.2 + 0.2

    REQUIRE(host_bseq->pulses.size() == 1);
    REQUIRE(host_bseq->pulses[0].id == 1);
    REQUIRE(host_bseq->pulses[0].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[0].endvalue].f64 == Approx(6.1));

    // Should have no effect
    host_seq.set_global(1, 5.9);

    REQUIRE(ages[0] == 3);
    REQUIRE(ages[1] == 4);
    REQUIRE(ages[2] == 5);
    REQUIRE(ages[3] == 3);
    REQUIRE(ages[4] == 0);
    REQUIRE(ages[5] == 4);
    REQUIRE(ages[6] == 0);

    REQUIRE(host_seq.post_run() == 2);
    REQUIRE(host_seq.cur_seq_idx() == 1);

    REQUIRE(ages[0] == 3);
    REQUIRE(ages[1] == 4);
    REQUIRE(ages[2] == 5);
    REQUIRE(ages[3] == 3);
    REQUIRE(ages[4] == 0);
    REQUIRE(ages[5] == 5);
    REQUIRE(ages[6] == 0);

    REQUIRE(host_seq.values[0].i64 == 123);
    REQUIRE(host_seq.values[1].f64 == -0.3);
    REQUIRE(host_seq.values[2].f64 == 5.9);
    REQUIRE(host_seq.values[3].f64 == Approx(6.1));
    REQUIRE(host_seq.values[4].f64 == Approx(0.2));
    // REQUIRE(host_seq.values[5].f64 == 0); // Unused by BS1
    REQUIRE(host_seq.values[6].b == false);
    // REQUIRE(host_seq.values[7].b == false); // Unused by BS1
    REQUIRE(host_seq.values[8].f64 == Approx(5.7));
    REQUIRE(host_seq.values[9].f64 == Approx(6.1));

    /**
     * Run 1 SubRun 4: BSeq 2
     */
    REQUIRE(host_seq.cur_seq_idx() == 1);
    host_seq.pre_run();
    host_bseq = &host_seq.seqs[1];
    REQUIRE(host_bseq->length == 123);
    REQUIRE(host_seq.start_values[0].f64 == Approx(6.1));

    REQUIRE(ages[0] == 3);
    REQUIRE(ages[1] == 4);
    REQUIRE(ages[2] == 6);
    REQUIRE(ages[3] == 3);
    REQUIRE(ages[4] == 5);
    REQUIRE(ages[5] == 5);
    REQUIRE(ages[6] == 0);

    REQUIRE(host_seq.values[0].i64 == 123);
    REQUIRE(host_seq.values[1].f64 == -0.3);
    REQUIRE(host_seq.values[2].f64 == 5.9);
    REQUIRE(host_seq.values[3].f64 == Approx(3.9));
    REQUIRE(host_seq.values[4].f64 == Approx(0.2));
    REQUIRE(host_seq.values[5].f64 == Approx(5.6)); // 5.9 + -0.3
    REQUIRE(host_seq.values[6].b == false);
    // REQUIRE(host_seq.values[7].b == false); // Unused by BS2
    REQUIRE(host_seq.values[8].f64 == Approx(6.1));
    REQUIRE(host_seq.values[9].f64 == Approx(3.9)); // 6.1 + 3.4 - 5.6

    REQUIRE(host_bseq->pulses.size() == 1);
    REQUIRE(host_bseq->pulses[0].id == 2);
    REQUIRE(host_bseq->pulses[0].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[0].endvalue].f64 == Approx(3.9));

    host_seq.set_global(0, 4.0);

    REQUIRE(ages[0] == 6);
    REQUIRE(ages[1] == 4);
    REQUIRE(ages[2] == 6);
    REQUIRE(ages[3] == 3);
    REQUIRE(ages[4] == 5);
    REQUIRE(ages[5] == 5);
    REQUIRE(ages[6] == 0);

    REQUIRE(host_seq.post_run() == 2);
    REQUIRE(host_seq.cur_seq_idx() == 1);

    REQUIRE(ages[0] == 6);
    REQUIRE(ages[1] == 4);
    REQUIRE(ages[2] == 6);
    REQUIRE(ages[3] == 3);
    REQUIRE(ages[4] == 5);
    REQUIRE(ages[5] == 5);
    REQUIRE(ages[6] == 6);

    REQUIRE(host_seq.values[0].i64 == 123);
    REQUIRE(host_seq.values[1].f64 == 4.0);
    REQUIRE(host_seq.values[2].f64 == 5.9);
    REQUIRE(host_seq.values[3].f64 == Approx(3.9));
    REQUIRE(host_seq.values[4].f64 == Approx(0.2));
    REQUIRE(host_seq.values[5].f64 == Approx(5.6));
    REQUIRE(host_seq.values[6].b == false);
    REQUIRE(host_seq.values[7].b == true);
    REQUIRE(host_seq.values[8].f64 == Approx(6.1));
    REQUIRE(host_seq.values[9].f64 == Approx(3.9));

    /**
     * Run 1 SubRun 5: BSeq 2
     */
    REQUIRE(host_seq.cur_seq_idx() == 1);
    host_seq.pre_run();
    host_bseq = &host_seq.seqs[1];
    REQUIRE(host_bseq->length == 123);
    REQUIRE(host_seq.start_values[0].f64 == Approx(3.9));

    REQUIRE(ages[0] == 6);
    REQUIRE(ages[1] == 4);
    REQUIRE(ages[2] == 7);
    REQUIRE(ages[3] == 3);
    REQUIRE(ages[4] == 6);
    REQUIRE(ages[5] == 5);
    REQUIRE(ages[6] == 6);

    REQUIRE(host_seq.values[0].i64 == 123);
    REQUIRE(host_seq.values[1].f64 == 4.0);
    REQUIRE(host_seq.values[2].f64 == 5.9);
    REQUIRE(host_seq.values[3].f64 == Approx(-2.6));
    REQUIRE(host_seq.values[4].f64 == Approx(0.2));
    REQUIRE(host_seq.values[5].f64 == Approx(9.9)); // 5.9 + 4.0
    REQUIRE(host_seq.values[6].b == false);
    REQUIRE(host_seq.values[7].b == true);
    REQUIRE(host_seq.values[8].f64 == Approx(3.9));
    REQUIRE(host_seq.values[9].f64 == Approx(-2.6)); // 3.9 + 3.4 - 9.9

    REQUIRE(host_bseq->pulses.size() == 1);
    REQUIRE(host_bseq->pulses[0].id == 2);
    REQUIRE(host_bseq->pulses[0].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[0].endvalue].f64 == Approx(-2.6));

    host_seq.set_global(0, -3.0);

    REQUIRE(ages[0] == 7);
    REQUIRE(ages[1] == 4);
    REQUIRE(ages[2] == 7);
    REQUIRE(ages[3] == 3);
    REQUIRE(ages[4] == 6);
    REQUIRE(ages[5] == 5);
    REQUIRE(ages[6] == 6);

    REQUIRE(host_seq.post_run() == 0);
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));

    REQUIRE(ages[0] == 7);
    REQUIRE(ages[1] == 4);
    REQUIRE(ages[2] == 7);
    REQUIRE(ages[3] == 3);
    REQUIRE(ages[4] == 6);
    REQUIRE(ages[5] == 5);
    REQUIRE(ages[6] == 7);

    REQUIRE(host_seq.values[0].i64 == 123);
    REQUIRE(host_seq.values[1].f64 == -3.0);
    REQUIRE(host_seq.values[2].f64 == 5.9);
    REQUIRE(host_seq.values[3].f64 == Approx(-2.6));
    REQUIRE(host_seq.values[4].f64 == Approx(0.2));
    REQUIRE(host_seq.values[5].f64 == Approx(9.9));
    REQUIRE(host_seq.values[6].b == false);
    REQUIRE(host_seq.values[7].b == false);
    REQUIRE(host_seq.values[8].f64 == Approx(3.9));
    REQUIRE(host_seq.values[9].f64 == Approx(-2.6));

    /**
     * Run 2
     */
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
    host_seq.init_run();
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(0));

    REQUIRE(ages[0] == 7);
    REQUIRE(ages[1] == 4);
    REQUIRE(ages[2] == 7);
    REQUIRE(ages[3] == 3);
    REQUIRE(ages[4] == 6);
    REQUIRE(ages[5] == 5);
    REQUIRE(ages[6] == 7);

    REQUIRE(host_seq.values[0].i64 == 123);
    // Global variable values remains unchanged
    REQUIRE(host_seq.values[1].f64 == -3.0);
    REQUIRE(host_seq.values[2].f64 == 5.9);
    REQUIRE(host_seq.values[3].f64 == Approx(-2.6));
    REQUIRE(host_seq.values[4].f64 == Approx(0.2));
    REQUIRE(host_seq.values[5].f64 == Approx(9.9));
    REQUIRE(host_seq.values[6].b == false);
    REQUIRE(host_seq.values[7].b == false);
    REQUIRE(host_seq.values[8].f64 == Approx(3.9));
    REQUIRE(host_seq.values[9].f64 == Approx(-2.6));

    host_seq.set_global(1, -2.6);

    REQUIRE(ages[0] == 7);
    REQUIRE(ages[1] == 8);
    REQUIRE(ages[2] == 7);
    REQUIRE(ages[3] == 3);
    REQUIRE(ages[4] == 6);
    REQUIRE(ages[5] == 5);
    REQUIRE(ages[6] == 7);

    REQUIRE(host_seq.values[0].i64 == 123);
    // Global variable values remains unchanged
    REQUIRE(host_seq.values[1].f64 == -3.0);
    REQUIRE(host_seq.values[2].f64 == -2.6);
    REQUIRE(host_seq.values[3].f64 == Approx(-2.6));
    REQUIRE(host_seq.values[4].f64 == Approx(0.2));
    REQUIRE(host_seq.values[5].f64 == Approx(9.9));
    REQUIRE(host_seq.values[6].b == false);
    REQUIRE(host_seq.values[7].b == false);
    REQUIRE(host_seq.values[8].f64 == Approx(3.9));
    REQUIRE(host_seq.values[9].f64 == Approx(-2.6));

    /**
     * Run 2 SubRun 1: BSeq 1
     */
    REQUIRE(host_seq.cur_seq_idx() == 0);
    host_seq.pre_run();
    host_bseq = &host_seq.seqs[0];
    REQUIRE(host_bseq->length == 123);
    REQUIRE(host_seq.start_values[0].f64 == 1.2);

    REQUIRE(ages[0] == 7);
    REQUIRE(ages[1] == 8);
    REQUIRE(ages[2] == 9);
    REQUIRE(ages[3] == 8);
    REQUIRE(ages[4] == 6);
    REQUIRE(ages[5] == 5);
    REQUIRE(ages[6] == 7);

    REQUIRE(host_seq.values[0].i64 == 123);
    REQUIRE(host_seq.values[1].f64 == -3.0);
    REQUIRE(host_seq.values[2].f64 == -2.6);
    REQUIRE(host_seq.values[3].f64 == Approx(-1.1));
    REQUIRE(host_seq.values[4].f64 == Approx(-2.5)); // -3.0 + 0.5
    REQUIRE(host_seq.values[5].f64 == Approx(9.9));
    REQUIRE(host_seq.values[6].b == false);
    REQUIRE(host_seq.values[7].b == false);
    REQUIRE(host_seq.values[8].f64 == Approx(1.2));
    REQUIRE(host_seq.values[9].f64 == Approx(-1.1)); // 1.2 + 0.2 + -2.5

    REQUIRE(host_bseq->pulses.size() == 1);
    REQUIRE(host_bseq->pulses[0].id == 1);
    REQUIRE(host_bseq->pulses[0].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[0].endvalue].f64 == Approx(-1.1));

    REQUIRE(host_seq.post_run() == 2);
    REQUIRE(host_seq.cur_seq_idx() == 1);

    REQUIRE(ages[0] == 7);
    REQUIRE(ages[1] == 8);
    REQUIRE(ages[2] == 9);
    REQUIRE(ages[3] == 8);
    REQUIRE(ages[4] == 6);
    REQUIRE(ages[5] == 9);
    REQUIRE(ages[6] == 7);

    REQUIRE(host_seq.values[0].i64 == 123);
    REQUIRE(host_seq.values[1].f64 == -3.0);
    REQUIRE(host_seq.values[2].f64 == -2.6);
    REQUIRE(host_seq.values[3].f64 == Approx(-1.1));
    REQUIRE(host_seq.values[4].f64 == Approx(-2.5)); // -3.0 + 0.5
    REQUIRE(host_seq.values[5].f64 == Approx(9.9));
    REQUIRE(host_seq.values[6].b == false);
    REQUIRE(host_seq.values[7].b == false);
    REQUIRE(host_seq.values[8].f64 == Approx(1.2));
    REQUIRE(host_seq.values[9].f64 == Approx(-1.1));

    /**
     * Run 2 SubRun 2: BSeq 2
     */
    REQUIRE(host_seq.cur_seq_idx() == 1);
    host_seq.pre_run();
    host_bseq = &host_seq.seqs[1];
    REQUIRE(host_bseq->length == 123);
    REQUIRE(host_seq.start_values[0].f64 == Approx(-1.1));

    REQUIRE(ages[0] == 7);
    REQUIRE(ages[1] == 8);
    REQUIRE(ages[2] == 10);
    REQUIRE(ages[3] == 8);
    // Age of G1 should jump to 9 (current age) even though max dependency age is 8
    REQUIRE(ages[4] == 9);
    REQUIRE(ages[5] == 9);
    REQUIRE(ages[6] == 7);

    REQUIRE(host_seq.values[0].i64 == 123);
    REQUIRE(host_seq.values[1].f64 == -3.0);
    REQUIRE(host_seq.values[2].f64 == -2.6);
    REQUIRE(host_seq.values[3].f64 == Approx(7.9));
    REQUIRE(host_seq.values[4].f64 == Approx(-2.5));
    REQUIRE(host_seq.values[5].f64 == Approx(-5.6)); // -3.0 + -2.6
    REQUIRE(host_seq.values[6].b == false);
    REQUIRE(host_seq.values[7].b == false);
    REQUIRE(host_seq.values[8].f64 == Approx(-1.1));
    REQUIRE(host_seq.values[9].f64 == Approx(7.9)); // -1.1 + 3.4 - -5.6

    REQUIRE(host_bseq->pulses.size() == 1);
    REQUIRE(host_bseq->pulses[0].id == 2);
    REQUIRE(host_bseq->pulses[0].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[0].endvalue].f64 == Approx(7.9));

    REQUIRE(host_seq.post_run() == 0);
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));

    REQUIRE(ages[0] == 7);
    REQUIRE(ages[1] == 8);
    REQUIRE(ages[2] == 10);
    REQUIRE(ages[3] == 8);
    REQUIRE(ages[4] == 9);
    REQUIRE(ages[5] == 9);
    REQUIRE(ages[6] == 10); // Unchanged value but age should update

    REQUIRE(host_seq.values[0].i64 == 123);
    REQUIRE(host_seq.values[1].f64 == -3.0);
    REQUIRE(host_seq.values[2].f64 == -2.6);
    REQUIRE(host_seq.values[3].f64 == Approx(7.9));
    REQUIRE(host_seq.values[4].f64 == Approx(-2.5));
    REQUIRE(host_seq.values[5].f64 == Approx(-5.6)); // -3.0 + -2.6
    REQUIRE(host_seq.values[6].b == false);
    REQUIRE(host_seq.values[7].b == false);
    REQUIRE(host_seq.values[8].f64 == Approx(-1.1));
    REQUIRE(host_seq.values[9].f64 == Approx(7.9));
}

TEST_CASE("neg_endtime") {
    Seq::HostSeq host_seq;
    host_seq.nconsts = 2;
    host_seq.nglobals = 0;
    host_seq.nglobal_vals = 0;
    host_seq.nchannels = 1;
    host_seq.nshared = (host_seq.nconsts + host_seq.nglobals +
                        host_seq.nglobal_vals + host_seq.nchannels);
    host_seq.npublic_globals = 0;

    host_seq.values.resize(4);
    host_seq.values[0].i64 = -1;
    host_seq.values[1].f64 = 2.3;
    host_seq.default_values.resize(1);
    host_seq.default_values[0].f64 = 1.2;
    host_seq.types.resize(2); // nconst + nglobals + nglobal_vals
    host_seq.types[0] = Seq::HostSeq::Type::Int64;
    host_seq.types[1] = Seq::HostSeq::Type::Float64;

    host_seq.depends.resize(0); // nglobal_vals
    host_seq.global_evals.resize(0); // nglobal_vals

    host_seq.seqs.resize(1);

    {
        auto &host_bseq = host_seq.seqs[0];
        host_bseq.id = 1;

        host_bseq.nmeasure = 0;
        host_bseq.ndirect = 0;
        host_bseq.nneed_order = 0;

        host_bseq.types.resize(0); // ndirect + nneed_order
        host_bseq.evals.resize(0); // ndirect + nneed_order

        host_bseq.global_refs.resize(0);
        host_bseq.reverse_depends.resize(0); // nmeasure
        host_bseq.deps_count.resize(0); // nneed_order

        host_bseq.pulses.resize(0);
        host_bseq.assignments.resize(0);
        host_bseq.assumptions.resize(0);
        host_bseq.branches.resize(0);
        host_bseq.default_branch = -1;
        host_bseq.endtimes.resize(1);
        host_bseq.endtimes[0] = 0;

        host_bseq.ndirect_assumes = 0;
        host_bseq.assumptions_idx.resize(0); // nneed_order

        host_bseq.cond_global_refs.resize(0);
    }

    host_seq.init();
    REQUIRE(host_seq.global_ages().empty());

    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
    host_seq.init_run();
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(0));
    host_seq.pre_run();
    auto host_bseq = &host_seq.seqs[0];
    REQUIRE(host_bseq->length == 0);
    REQUIRE(host_seq.post_run() == 0);
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
}

TEST_CASE("need_order_assign") {
    // To get a need_order value, we need two pulses with known order;
    // a measure that depends on the order of the two,
    // and a value that depends on that measure... Sigh...
    // Const:
    //   C0[V0]: time 1 (20)
    //   C1[V1]: val 1 (1.2)
    //   C2[V2]: val 2 (2.3)
    // Slot:
    //   S0[V3]: for time 2
    //   S1[V4]: assignment
    // Global Value:
    //   G0[V5]: round S0 (time 2)
    //   G1[V6]: C0 + G0
    // Channel:
    //   CV1[V7]: Channel 1 value
    // BS1:
    //   Values:
    //     MV0[V8]: measure value of M1
    //     N0[V9]: MV0 + 3.4
    //   Pulse:
    //     P1: @C0 val=C1
    //     P2: @G0 val=C2
    //     M1: @G1 measure=MV0
    //   Assignments:
    //     S1 <- N0
    //   Branch:
    //     Default: end
    Seq::HostSeq host_seq;
    host_seq.nconsts = 3;
    host_seq.nglobals = 2;
    host_seq.nglobal_vals = 2;
    host_seq.nchannels = 1;
    host_seq.nshared = (host_seq.nconsts + host_seq.nglobals +
                        host_seq.nglobal_vals + host_seq.nchannels);
    host_seq.npublic_globals = 2;

    host_seq.values.resize(10);
    host_seq.values[0].i64 = 20;
    host_seq.values[1].f64 = 1.2;
    host_seq.values[2].f64 = 2.3;
    host_seq.default_values.resize(1);
    host_seq.default_values[0].f64 = 0.2;
    host_seq.types.resize(7); // nconst + nglobals + nglobal_vals
    // Const
    host_seq.types[0] = Seq::HostSeq::Type::Int64;
    host_seq.types[1] = Seq::HostSeq::Type::Float64;
    host_seq.types[2] = Seq::HostSeq::Type::Float64;
    // Slot
    host_seq.types[3] = Seq::HostSeq::Type::Float64;
    host_seq.types[4] = Seq::HostSeq::Type::Float64;
    // Global Value
    host_seq.types[5] = Seq::HostSeq::Type::Int64;
    host_seq.types[6] = Seq::HostSeq::Type::Int64;

    host_seq.depends.resize(2); // nglobal_vals
    host_seq.depends[0].push_back(0);
    host_seq.depends[1].push_back(0);
    host_seq.global_evals.resize(2); // nglobal_vals
    host_seq.global_evals[0] = [] (void *_p) {
        auto p = (Seq::HostSeq::Value*)_p;
        p[5].i64 = round<int64_t>(p[3].f64);
    };
    host_seq.global_evals[1] = [] (void *_p) {
        auto p = (Seq::HostSeq::Value*)_p;
        p[6].i64 = p[5].i64 + 20;
    };

    host_seq.seqs.resize(1);

    {
        auto &host_bseq = host_seq.seqs[0];
        host_bseq.id = 1;

        host_bseq.nmeasure = 1;
        host_bseq.ndirect = 0;
        host_bseq.nneed_order = 1;

        host_bseq.types.resize(1); // ndirect + nneed_order
        host_bseq.types[0] = Seq::HostSeq::Type::Float64;
        host_bseq.evals.resize(1); // ndirect + nneed_order
        host_bseq.evals[0] = [] (void *_p) {
            auto p = (Seq::HostSeq::Value*)_p;
            p[9].f64 = p[8].f64 + 3.4;
        };

        host_bseq.global_refs.resize(2);
        host_bseq.global_refs[0] = 0;
        host_bseq.global_refs[1] = 1;
        host_bseq.reverse_depends.resize(1); // nmeasure
        host_bseq.reverse_depends[0].push_back(0);
        host_bseq.deps_count.resize(1); // nneed_order
        host_bseq.deps_count[0] = 1;

        host_bseq.pulses.resize(3);
        host_bseq.pulses[0] = Seq::HostSeq::Pulse{
            .id = 1, .time = 0, .measure = uint32_t(-1), .len = uint32_t(-1),
            .value = 1, .chn = 1, .ramp_func = nullptr, .endvalue = 1, .cond = uint32_t(-1)
        };
        host_bseq.pulses[1] = Seq::HostSeq::Pulse{
            .id = 2, .time = 5, .measure = uint32_t(-1), .len = uint32_t(-1),
            .value = 2, .chn = 1, .ramp_func = nullptr, .endvalue = 2, .cond = uint32_t(-1)
        };
        host_bseq.pulses[2] = Seq::HostSeq::Pulse{
            .id = 3, .time = 6, .measure = 0, .len = uint32_t(-2),
            .value = uint32_t(-1), .chn = 1, .ramp_func = nullptr,
            .endvalue = uint32_t(-1), .cond = uint32_t(-1)
        };
        host_bseq.assignments.resize(1);
        host_bseq.assignments[0] = Seq::HostSeq::Assignment{ .global_id = 1, .value = 9 };
        host_bseq.assumptions.resize(0);
        host_bseq.branches.resize(0);
        host_bseq.default_branch = -1;
        host_bseq.endtimes.resize(1);
        host_bseq.endtimes[0] = 6;

        host_bseq.ndirect_assumes = 0;
        host_bseq.assumptions_idx.resize(1); // nneed_order
        host_bseq.assumptions_idx[0] = uint32_t(-1);

        host_bseq.cond_global_refs.resize(0);
    }

    host_seq.init();
    auto &ages = host_seq.global_ages();
    REQUIRE(ages.size() == 4);
    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 1);
    REQUIRE(ages[2] == 0);
    REQUIRE(ages[3] == 0);
    host_seq.set_global(0, 0.4); // this will make G0 == 0
    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 1);
    REQUIRE(ages[2] == 0);
    REQUIRE(ages[3] == 0);

    // Run 1
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
    host_seq.init_run();
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(0));

    host_seq.pre_run();
    auto host_bseq = &host_seq.seqs[0];
    REQUIRE(host_bseq->length == 20);
    REQUIRE(host_seq.start_values[0].f64 == Approx(2.3));

    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 2);
    REQUIRE(ages[2] == 1);
    REQUIRE(ages[3] == 1);

    REQUIRE(host_seq.values[0].i64 == 20);
    REQUIRE(host_seq.values[1].f64 == 1.2);
    REQUIRE(host_seq.values[2].f64 == 2.3);
    REQUIRE(host_seq.values[3].f64 == 0.4);
    REQUIRE(host_seq.values[4].f64 == Approx(4.6));
    REQUIRE(host_seq.values[5].i64 == 0);
    REQUIRE(host_seq.values[6].i64 == 20);
    REQUIRE(host_seq.values[7].f64 == 0.2);
    REQUIRE(host_seq.values[8].f64 == 1.2);
    REQUIRE(host_seq.values[9].f64 == Approx(4.6));

    REQUIRE(host_bseq->pulses.size() == 3);
    REQUIRE(host_bseq->pulses[0].id == 2);
    REQUIRE(host_bseq->pulses[0].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[0].endvalue].f64 == 2.3);
    REQUIRE(host_bseq->pulses[1].id == 1);
    REQUIRE(host_bseq->pulses[1].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[1].endvalue].f64 == 1.2);
    REQUIRE(host_bseq->pulses[2].id == 3); // measure
    REQUIRE(host_bseq->pulses[2].chn == 1);
    REQUIRE(host_bseq->pulses[2].len == uint32_t(-2));

    REQUIRE(host_seq.post_run() == 0);
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));

    // Run 2
    host_seq.init_run();
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(0));
    host_seq.set_global(0, 25.6);

    REQUIRE(ages[0] == 3);
    REQUIRE(ages[1] == 2);
    REQUIRE(ages[2] == 1);
    REQUIRE(ages[3] == 1);

    host_seq.pre_run();
    host_bseq = &host_seq.seqs[0];
    REQUIRE(host_bseq->length == 46);
    REQUIRE(host_seq.start_values[0].f64 == Approx(0.2));

    REQUIRE(ages[0] == 3);
    REQUIRE(ages[1] == 4);
    REQUIRE(ages[2] == 3);
    REQUIRE(ages[3] == 3);

    REQUIRE(host_seq.values[0].i64 == 20);
    REQUIRE(host_seq.values[1].f64 == 1.2);
    REQUIRE(host_seq.values[2].f64 == 2.3);
    REQUIRE(host_seq.values[3].f64 == 25.6);
    REQUIRE(host_seq.values[4].f64 == Approx(5.7));
    REQUIRE(host_seq.values[5].i64 == 26);
    REQUIRE(host_seq.values[6].i64 == 46);
    REQUIRE(host_seq.values[7].f64 == 0.2);
    REQUIRE(host_seq.values[8].f64 == 2.3);
    REQUIRE(host_seq.values[9].f64 == Approx(5.7));

    REQUIRE(host_bseq->pulses.size() == 3);
    REQUIRE(host_bseq->pulses[0].id == 1);
    REQUIRE(host_bseq->pulses[0].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[0].endvalue].f64 == 1.2);
    REQUIRE(host_bseq->pulses[1].id == 2);
    REQUIRE(host_bseq->pulses[1].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[1].endvalue].f64 == 2.3);
    REQUIRE(host_bseq->pulses[2].id == 3); // measure
    REQUIRE(host_bseq->pulses[2].chn == 1);
    REQUIRE(host_bseq->pulses[2].len == uint32_t(-2));

    REQUIRE(host_seq.post_run() == 0);
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
}

TEST_CASE("measure_order_id") {
    // Const:
    //   C0[V0]: time (20)
    //   C1[V1]: val 1 (1.2)
    //   C2[V2]: val 2 (2.3)
    // Slot:
    //   S0[V3]: assignment
    // Global Value:
    // Channel:
    //   CV1[V4]: Channel 1 value
    // BS1:
    //   Values:
    //     MV0[V5]: measure value of M1
    //   Pulse:
    //     P1: @C0 val=C1
    //     M1: @C0 measure=MV0
    //     P2: @C0 val=C2
    //   Assignments:
    //     S0 <- MV0
    //   Branch:
    //     Default: end
    Seq::HostSeq host_seq;
    host_seq.nconsts = 3;
    host_seq.nglobals = 1;
    host_seq.nglobal_vals = 0;
    host_seq.nchannels = 1;
    host_seq.nshared = (host_seq.nconsts + host_seq.nglobals +
                        host_seq.nglobal_vals + host_seq.nchannels);
    host_seq.npublic_globals = 1;

    host_seq.values.resize(6);
    host_seq.values[0].i64 = 20;
    host_seq.values[1].f64 = 1.2;
    host_seq.values[2].f64 = 2.3;
    host_seq.default_values.resize(1);
    host_seq.default_values[0].f64 = 0.2;
    host_seq.types.resize(4); // nconst + nglobals + nglobal_vals
    // Const
    host_seq.types[0] = Seq::HostSeq::Type::Int64;
    host_seq.types[1] = Seq::HostSeq::Type::Float64;
    host_seq.types[2] = Seq::HostSeq::Type::Float64;
    // Slot
    host_seq.types[3] = Seq::HostSeq::Type::Float64;

    host_seq.depends.resize(0); // nglobal_vals
    host_seq.global_evals.resize(0); // nglobal_vals

    host_seq.seqs.resize(1);

    {
        auto &host_bseq = host_seq.seqs[0];
        host_bseq.id = 1;

        host_bseq.nmeasure = 1;
        host_bseq.ndirect = 0;
        host_bseq.nneed_order = 0;

        host_bseq.types.resize(0); // ndirect + nneed_order
        host_bseq.evals.resize(0); // ndirect + nneed_order

        host_bseq.global_refs.resize(0);
        host_bseq.reverse_depends.resize(1); // nmeasure
        host_bseq.deps_count.resize(0); // nneed_order

        host_bseq.pulses.resize(3);
        host_bseq.pulses[0] = Seq::HostSeq::Pulse{
            .id = 1, .time = 0, .measure = uint32_t(-1), .len = uint32_t(-1),
            .value = 1, .chn = 1, .ramp_func = nullptr, .endvalue = 1, .cond = uint32_t(-1)
        };
        host_bseq.pulses[1] = Seq::HostSeq::Pulse{
            .id = 3, .time = 0, .measure = uint32_t(-1), .len = uint32_t(-1),
            .value = 2, .chn = 1, .ramp_func = nullptr, .endvalue = 2, .cond = uint32_t(-1)
        };
        host_bseq.pulses[2] = Seq::HostSeq::Pulse{
            .id = 2, .time = 0, .measure = 0, .len = uint32_t(-2),
            .value = uint32_t(-1), .chn = 1, .ramp_func = nullptr,
            .endvalue = uint32_t(-1), .cond = uint32_t(-1)
        };
        host_bseq.assignments.resize(1);
        host_bseq.assignments[0] = Seq::HostSeq::Assignment{ .global_id = 0, .value = 5 };
        host_bseq.assumptions.resize(0);
        host_bseq.branches.resize(0);
        host_bseq.default_branch = -1;
        host_bseq.endtimes.resize(1);
        host_bseq.endtimes[0] = 0;

        host_bseq.ndirect_assumes = 0;
        host_bseq.assumptions_idx.resize(0); // nneed_order

        host_bseq.cond_global_refs.resize(0);
    }

    host_seq.init();
    auto &ages = host_seq.global_ages();
    REQUIRE(ages.size() == 1);
    REQUIRE(ages[0] == 1);

    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
    host_seq.init_run();
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(0));

    host_seq.pre_run();
    auto host_bseq = &host_seq.seqs[0];
    REQUIRE(host_bseq->length == 20);
    REQUIRE(host_seq.start_values[0].f64 == 0.2);

    REQUIRE(ages[0] == 2);

    REQUIRE(host_seq.values[0].i64 == 20);
    REQUIRE(host_seq.values[1].f64 == 1.2);
    REQUIRE(host_seq.values[2].f64 == 2.3);
    REQUIRE(host_seq.values[3].f64 == 1.2);
    REQUIRE(host_seq.values[4].f64 == 0.2);
    REQUIRE(host_seq.values[5].f64 == 1.2);

    REQUIRE(host_bseq->pulses.size() == 3);
    REQUIRE(host_bseq->pulses[0].id == 1);
    REQUIRE(host_bseq->pulses[0].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[0].endvalue].f64 == 1.2);
    REQUIRE(host_bseq->pulses[1].id == 2); // measure
    REQUIRE(host_bseq->pulses[1].chn == 1);
    REQUIRE(host_bseq->pulses[1].len == uint32_t(-2));
    REQUIRE(host_bseq->pulses[2].id == 3);
    REQUIRE(host_bseq->pulses[2].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[2].endvalue].f64 == 2.3);

    REQUIRE(host_seq.post_run() == 0);
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
}

TEST_CASE("ramp_measure") {
    // Test subjects:
    // * Need order value with more than one dependencies
    // * Start value of pulse capturing end value of previous pulse with overlap
    // * Start value of pulse capturing end value of previous pulse without overlap
    // * Measure capture mid value of a pulse.
    // * Measure capture end value of a pulse.
    // * Resort pulse when new time is determined
    //   (pulse at later time having known time before pulse at a earlier time does)

    // Const:
    //   C0[V0]: time for P1 (200)
    //   C1[V1]: length for P1 (1000)
    // Slot:
    //   S0[V2]: for time for P2
    //   S1[V3]: for time for M2+P3
    //   S2[V4]: for time for P4
    // Global Value:
    //   G0[V5]: round(S0) (time for P2)
    //   G1[V6]: G0 + round(S1) (time for M2+P3)
    //   G2[V7]: G1 + 5000 (time for P5)
    // Channel:
    //   CV1[V8]: Channel 1 value
    // BS1:
    //   Values:
    //     MV0[V9]: measure value of P1
    //     MV1[V10]: measure value of P2
    //     MV2[V11]: measure value of M1
    //     MV3[V12]: measure value of P3
    //     MV4[V13]: measure value of P4
    //     MV5[V14]: measure value of P5
    //     N0[V15]: MV0 + t * 0.001
    //     N1[V16]: MV1 * 1.3
    //     N2[V17]: MV3 * 0.8 + 0.7
    //     N3[V18]: round(500 + MV2 * 1000 + N2 * 900) + round(S2)
    //     N4[V19]: MV4 - 5.7
    //     N5[V20]: MV5 * 2.8
    //     N6[V21]: MV0 + 1 // End value of P1
    //   Pulse:
    //     P1(id=1): @C0 measure=MV0 ramp=[N0] MV0 + t * 0.001 len=1000 endvalue = N6
    //     P2(id=2): @G0 measure=MV1 val=N1
    //     M1(id=3): @G1 measure=MV2
    //     P3(id=4): @G1 measure=MV3 val=N2
    //     P4(id=5): @N3 measure=MV4 val=N4
    //     P5(id=6): @G2 measure=MV5 val=N5
    //   Branch:
    //     Default: end

    // Conditions:
    //   C0 > G0 / C0 < G0 < C0 + 1000 / C0 + 1000 < G0
    //   G1 > G0, G1 > C0
    //       (round(S1) > 0, round(S0) + round(S1) > 200)
    //   G1 > C0 + 1000 / G1 < C0 + 1000
    //   N3 > G1, G2 > G1
    //       (round(500 + MV2 * 1000 + N2 * 900) + round(S2) > round(S0) + round(S1))
    //   N3 <= G2 / G2 < N3

    Seq::HostSeq host_seq;
    host_seq.nconsts = 2;
    host_seq.nglobals = 3;
    host_seq.nglobal_vals = 3;
    host_seq.nchannels = 1;
    host_seq.nshared = (host_seq.nconsts + host_seq.nglobals +
                        host_seq.nglobal_vals + host_seq.nchannels);
    host_seq.npublic_globals = 3;

    host_seq.values.resize(22);
    host_seq.values[0].i64 = 200;
    host_seq.values[1].f64 = 1000;
    host_seq.default_values.resize(1);
    host_seq.default_values[0].f64 = -0.4;
    host_seq.types.resize(8); // nconst + nglobals + nglobal_vals
    // Const
    host_seq.types[0] = Seq::HostSeq::Type::Int64;
    host_seq.types[1] = Seq::HostSeq::Type::Float64;
    // Slot
    host_seq.types[2] = Seq::HostSeq::Type::Float64;
    host_seq.types[3] = Seq::HostSeq::Type::Float64;
    host_seq.types[4] = Seq::HostSeq::Type::Float64;
    // Global Value
    host_seq.types[5] = Seq::HostSeq::Type::Int64;
    host_seq.types[6] = Seq::HostSeq::Type::Int64;
    host_seq.types[7] = Seq::HostSeq::Type::Int64;

    host_seq.depends.resize(3); // nglobal_vals
    host_seq.depends[0].push_back(0);
    host_seq.depends[1].push_back(0);
    host_seq.depends[1].push_back(1);
    host_seq.depends[2].push_back(0);
    host_seq.depends[2].push_back(1);
    host_seq.global_evals.resize(3); // nglobal_vals
    host_seq.global_evals[0] = [] (void *_p) {
        auto p = (Seq::HostSeq::Value*)_p;
        p[5].i64 = round<int64_t>(p[2].f64);
    };
    host_seq.global_evals[1] = [] (void *_p) {
        auto p = (Seq::HostSeq::Value*)_p;
        p[6].i64 = p[5].i64 + round<int64_t>(p[3].f64);
    };
    host_seq.global_evals[2] = [] (void *_p) {
        auto p = (Seq::HostSeq::Value*)_p;
        p[7].i64 = p[6].i64 + 5000;
    };

    host_seq.seqs.resize(1);

    {
        auto &host_bseq = host_seq.seqs[0];
        host_bseq.id = 1;

        host_bseq.nmeasure = 6;
        host_bseq.ndirect = 0;
        host_bseq.nneed_order = 7;

        host_bseq.types.resize(7); // ndirect + nneed_order
        host_bseq.types[0] = Seq::HostSeq::Type::Float64;
        host_bseq.types[1] = Seq::HostSeq::Type::Float64;
        host_bseq.types[2] = Seq::HostSeq::Type::Float64;
        host_bseq.types[3] = Seq::HostSeq::Type::Int64;
        host_bseq.types[4] = Seq::HostSeq::Type::Float64;
        host_bseq.types[5] = Seq::HostSeq::Type::Float64;
        host_bseq.types[6] = Seq::HostSeq::Type::Float64;
        host_bseq.evals.resize(7); // ndirect + nneed_order
        host_bseq.evals[0] = [] (void*) {};
        host_bseq.evals[1] = [] (void *_p) {
            auto p = (Seq::HostSeq::Value*)_p;
            p[16].f64 = p[10].f64 * 1.3;
        };
        host_bseq.evals[2] = [] (void *_p) {
            auto p = (Seq::HostSeq::Value*)_p;
            p[17].f64 = p[12].f64 * 0.8 + 0.7;
        };
        host_bseq.evals[3] = [] (void *_p) {
            auto p = (Seq::HostSeq::Value*)_p;
            p[18].i64 = 500 + round<int64_t>(p[11].f64 * 1000 + p[17].f64 * 900) +
                round<int64_t>(p[4].f64);
        };
        host_bseq.evals[4] = [] (void *_p) {
            auto p = (Seq::HostSeq::Value*)_p;
            p[19].f64 = p[13].f64 - 5.7;
        };
        host_bseq.evals[5] = [] (void *_p) {
            auto p = (Seq::HostSeq::Value*)_p;
            p[20].f64 = p[14].f64 * 2.8;
        };
        host_bseq.evals[6] = [] (void *_p) {
            auto p = (Seq::HostSeq::Value*)_p;
            p[21].f64 = p[9].f64 + 1;
        };

        host_bseq.global_refs.resize(3);
        host_bseq.global_refs[0] = 0;
        host_bseq.global_refs[1] = 1;
        host_bseq.global_refs[2] = 2;
        host_bseq.reverse_depends.resize(6); // nmeasure
        host_bseq.reverse_depends[0].push_back(0);
        host_bseq.reverse_depends[0].push_back(6);
        host_bseq.reverse_depends[1].push_back(1);
        host_bseq.reverse_depends[2].push_back(3);
        host_bseq.reverse_depends[3].push_back(2);
        host_bseq.reverse_depends[3].push_back(3);
        host_bseq.reverse_depends[4].push_back(4);
        host_bseq.reverse_depends[5].push_back(5);
        host_bseq.deps_count.resize(7); // nneed_order
        host_bseq.deps_count[0] = 1;
        host_bseq.deps_count[1] = 1;
        host_bseq.deps_count[2] = 1;
        host_bseq.deps_count[3] = 2;
        host_bseq.deps_count[4] = 1;
        host_bseq.deps_count[5] = 1;
        host_bseq.deps_count[6] = 1;

        host_bseq.pulses.resize(6);
        host_bseq.pulses[0] = Seq::HostSeq::Pulse{
            .id = 1, .time = 0, .measure = 0, .len = 1,
            .value = 15, .chn = 1, .ramp_func = [] (double t, void *_p) {
                auto p = (Seq::HostSeq::Value*)_p;
                return p[9].f64 + t * 0.001;
            }, .endvalue = 21, .cond = uint32_t(-1)
        };
        host_bseq.pulses[1] = Seq::HostSeq::Pulse{
            .id = 2, .time = 5, .measure = 1, .len = uint32_t(-1),
            .value = 16, .chn = 1, .ramp_func = nullptr, .endvalue = 16, .cond = uint32_t(-1)
        };
        host_bseq.pulses[2] = Seq::HostSeq::Pulse{
            .id = 3, .time = 6, .measure = 2, .len = uint32_t(-2),
            .value = uint32_t(-1), .chn = 1, .ramp_func = nullptr,
            .endvalue = uint32_t(-1), .cond = uint32_t(-1)
        };
        host_bseq.pulses[3] = Seq::HostSeq::Pulse{
            .id = 4, .time = 6, .measure = 3, .len = uint32_t(-1),
            .value = 17, .chn = 1, .ramp_func = nullptr, .endvalue = 17, .cond = uint32_t(-1)
        };
        host_bseq.pulses[4] = Seq::HostSeq::Pulse{
            .id = 5, .time = 18, .measure = 4, .len = uint32_t(-1),
            .value = 19, .chn = 1, .ramp_func = nullptr, .endvalue = 19, .cond = uint32_t(-1)
        };
        host_bseq.pulses[5] = Seq::HostSeq::Pulse{
            .id = 6, .time = 7, .measure = 5, .len = uint32_t(-1),
            .value = 20, .chn = 1, .ramp_func = nullptr, .endvalue = 20, .cond = uint32_t(-1)
        };
        host_bseq.assignments.resize(0);
        host_bseq.assumptions.resize(0);
        host_bseq.branches.resize(0);
        host_bseq.default_branch = -1;
        host_bseq.endtimes.resize(2);
        host_bseq.endtimes[0] = 7;
        host_bseq.endtimes[1] = 18;

        host_bseq.ndirect_assumes = 0;
        host_bseq.assumptions_idx.resize(7); // nneed_order
        host_bseq.assumptions_idx[0] = uint32_t(-1);
        host_bseq.assumptions_idx[1] = uint32_t(-1);
        host_bseq.assumptions_idx[2] = uint32_t(-1);
        host_bseq.assumptions_idx[3] = uint32_t(-1);
        host_bseq.assumptions_idx[4] = uint32_t(-1);
        host_bseq.assumptions_idx[5] = uint32_t(-1);
        host_bseq.assumptions_idx[6] = uint32_t(-1);

        host_bseq.cond_global_refs.resize(0);
    }

    host_seq.init();
    auto &ages = host_seq.global_ages();
    REQUIRE(ages.size() == 6);
    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 1);
    REQUIRE(ages[2] == 1);
    REQUIRE(ages[3] == 0);
    REQUIRE(ages[4] == 0);
    REQUIRE(ages[5] == 0);
    auto host_bseq = &host_seq.seqs[0];

    // Run 1
    //   C0 > G0
    //   G1 > C0 + 1000
    //   N2 < G2
    host_seq.set_global(0, 100);
    host_seq.set_global(1, 1300);
    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 1);
    REQUIRE(ages[2] == 1);
    REQUIRE(ages[3] == 0);
    REQUIRE(ages[4] == 0);
    REQUIRE(ages[5] == 0);
    // start: -0.4
    // P2(id=2): @G0=100 measure=(MV1=-0.4) val=(N1=-0.4 * 1.3=-0.52)
    // P1(id=1): @C0=200 measure=(MV0=-0.52) ramp=[N0] (MV0=-0.52) + t * 0.001 len=1000
    // M1(id=3): @G1=1300 measure=(MV2=-0.52 + 1000 * 0.001=0.48)
    // P3(id=4): @G1=1300 measure=(MV3=-0.52 + 1000 * 0.001=0.48) val=(N2=0.48 * 0.8 + 0.7=1.084)
    // P4(id=5): @N3=1956 measure=(MV4=1.084) val=(N4=1.084 - 5.7=-4.616)
    // P5(id=6): @G2=6400 measure=(MV5=-4.616) val=(N5=-4.616 * 2.8=-12.9248)

    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
    host_seq.init_run();
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(0));

    host_seq.pre_run();
    REQUIRE(host_bseq->length == 6400);
    REQUIRE(host_seq.start_values[0].f64 == Approx(-0.4));

    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 1);
    REQUIRE(ages[2] == 1);
    REQUIRE(ages[3] == 1);
    REQUIRE(ages[4] == 1);
    REQUIRE(ages[5] == 1);

    // Consts
    REQUIRE(host_seq.values[0].i64 == 200);
    REQUIRE(host_seq.values[1].f64 == 1000);
    // Slots
    REQUIRE(host_seq.values[2].f64 == 100);
    REQUIRE(host_seq.values[3].f64 == 1300);
    REQUIRE(host_seq.values[4].f64 == 0);
    // Global Values
    REQUIRE(host_seq.values[5].i64 == 100);
    REQUIRE(host_seq.values[6].i64 == 1400);
    REQUIRE(host_seq.values[7].i64 == 6400);
    // Channels
    REQUIRE(host_seq.values[8].f64 == -0.4);
    // Measures
    REQUIRE(host_seq.values[9].f64 == Approx(-0.52));
    REQUIRE(host_seq.values[10].f64 == -0.4);
    REQUIRE(host_seq.values[11].f64 == Approx(0.48));
    REQUIRE(host_seq.values[12].f64 == Approx(0.48));
    REQUIRE(host_seq.values[13].f64 == Approx(1.084));
    REQUIRE(host_seq.values[14].f64 == Approx(-4.616));
    // Need Orders
    // V15 is a ramp
    REQUIRE(host_seq.values[16].f64 == Approx(-0.52));
    REQUIRE(host_seq.values[17].f64 == Approx(1.084));
    REQUIRE(host_seq.values[18].i64 == 1956);
    REQUIRE(host_seq.values[19].f64 == Approx(-4.616));
    REQUIRE(host_seq.values[20].f64 == Approx(-12.9248));

    REQUIRE(host_bseq->pulses.size() == 6);
    REQUIRE(host_bseq->pulses[0].id == 2);
    REQUIRE(host_bseq->pulses[0].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[0].endvalue].f64 == Approx(-0.52));
    REQUIRE(host_bseq->pulses[1].id == 1);
    REQUIRE(host_bseq->pulses[1].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[1].endvalue].f64 == Approx(0.48));
    REQUIRE(host_bseq->pulses[2].id == 3); // measure
    REQUIRE(host_bseq->pulses[2].chn == 1);
    REQUIRE(host_bseq->pulses[2].len == uint32_t(-2));
    REQUIRE(host_bseq->pulses[3].id == 4);
    REQUIRE(host_bseq->pulses[3].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[3].endvalue].f64 == Approx(1.084));
    REQUIRE(host_bseq->pulses[4].id == 5);
    REQUIRE(host_bseq->pulses[4].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[4].endvalue].f64 == Approx(-4.616));
    REQUIRE(host_bseq->pulses[5].id == 6);
    REQUIRE(host_bseq->pulses[5].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[5].endvalue].f64 == Approx(-12.9248));

    REQUIRE(host_seq.post_run() == 0);
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));

    // Run 2
    //   C0 > G0
    //   G1 > C0 + 1000
    //   N2 > G2
    host_seq.set_global(2, 5000);
    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 1);
    REQUIRE(ages[2] == 2);
    REQUIRE(ages[3] == 1);
    REQUIRE(ages[4] == 1);
    REQUIRE(ages[5] == 1);
    // start: -0.4
    // P2(id=2): @G0=100 measure=(MV1=-0.4) val=(N1=-0.4 * 1.3=-0.52)
    // P1(id=1): @C0=200 measure=(MV0=-0.52) ramp=[N0] (MV0=-0.52) + t * 0.001 len=1000
    // M1(id=3): @G1=1300 measure=(MV2=-0.52 + 1000 * 0.001=0.48)
    // P3(id=4): @G1=1300 measure=(MV3=-0.52 + 1000 * 0.001=0.48) val=(N2=0.48 * 0.8 + 0.7=1.084)
    // P5(id=6): @G2=6400 measure=(MV5=1.084) val=(N5=1.084 * 2.8=3.0352)
    // P4(id=5): @N3=6956 measure=(MV4=3.0352) val=(N4=3.0352 - 5.7=-2.6648)

    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
    host_seq.init_run();
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(0));

    host_seq.pre_run();
    REQUIRE(host_bseq->length == 6956);
    REQUIRE(host_seq.start_values[0].f64 == Approx(-0.4));

    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 1);
    REQUIRE(ages[2] == 2);
    REQUIRE(ages[3] == 1);
    REQUIRE(ages[4] == 1);
    REQUIRE(ages[5] == 1);

    // Consts
    REQUIRE(host_seq.values[0].i64 == 200);
    REQUIRE(host_seq.values[1].f64 == 1000);
    // Slots
    REQUIRE(host_seq.values[2].f64 == 100);
    REQUIRE(host_seq.values[3].f64 == 1300);
    REQUIRE(host_seq.values[4].f64 == 5000);
    // Global Values
    REQUIRE(host_seq.values[5].i64 == 100);
    REQUIRE(host_seq.values[6].i64 == 1400);
    REQUIRE(host_seq.values[7].i64 == 6400);
    // Channels
    REQUIRE(host_seq.values[8].f64 == -0.4);
    // Measures
    REQUIRE(host_seq.values[9].f64 == Approx(-0.52));
    REQUIRE(host_seq.values[10].f64 == -0.4);
    REQUIRE(host_seq.values[11].f64 == Approx(0.48));
    REQUIRE(host_seq.values[12].f64 == Approx(0.48));
    REQUIRE(host_seq.values[13].f64 == Approx(3.0352));
    REQUIRE(host_seq.values[14].f64 == Approx(1.084));
    // Need Orders
    // V15 is a ramp
    REQUIRE(host_seq.values[16].f64 == Approx(-0.52));
    REQUIRE(host_seq.values[17].f64 == Approx(1.084));
    REQUIRE(host_seq.values[18].i64 == 6956);
    REQUIRE(host_seq.values[19].f64 == Approx(-2.6648));
    REQUIRE(host_seq.values[20].f64 == Approx(3.0352));

    REQUIRE(host_bseq->pulses.size() == 6);
    REQUIRE(host_bseq->pulses[0].id == 2);
    REQUIRE(host_bseq->pulses[0].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[0].endvalue].f64 == Approx(-0.52));
    REQUIRE(host_bseq->pulses[1].id == 1);
    REQUIRE(host_bseq->pulses[1].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[1].endvalue].f64 == Approx(0.48));
    REQUIRE(host_bseq->pulses[2].id == 3); // measure
    REQUIRE(host_bseq->pulses[2].chn == 1);
    REQUIRE(host_bseq->pulses[2].len == uint32_t(-2));
    REQUIRE(host_bseq->pulses[3].id == 4);
    REQUIRE(host_bseq->pulses[3].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[3].endvalue].f64 == Approx(1.084));
    REQUIRE(host_bseq->pulses[4].id == 6);
    REQUIRE(host_bseq->pulses[4].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[4].endvalue].f64 == Approx(3.0352));
    REQUIRE(host_bseq->pulses[5].id == 5);
    REQUIRE(host_bseq->pulses[5].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[5].endvalue].f64 == Approx(-2.6648));

    REQUIRE(host_seq.post_run() == 0);
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));

    // Run 3
    //   C0 > G0
    //   G1 < C0 + 1000
    //   N3 == G2
    host_seq.set_global(1, 700);
    host_seq.set_global(2, 4244);
    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 3);
    REQUIRE(ages[2] == 3);
    REQUIRE(ages[3] == 1);
    REQUIRE(ages[4] == 1);
    REQUIRE(ages[5] == 1);
    // start: -0.4
    // P2(id=2): @G0=100 measure=(MV1=-0.4) val=(N1=-0.4 * 1.3=-0.52)
    // P1(id=1): @C0=200 measure=(MV0=-0.52) ramp=[N0] (MV0=-0.52) + t * 0.001 len=1000
    // M1(id=3): @G1=800 measure=(MV2=-0.52 + 600 * 0.001=0.08)
    // P3(id=4): @G1=800 measure=(MV3=-0.52 + 1000 * 0.001=0.48) val=(N2=0.48 * 0.8 + 0.7=1.084)
    // P4(id=5): @N3=5800 measure=(MV4=1.084) val=(N4=-4.616)
    // P5(id=6): @G2=5800 measure=(MV5=-4.616) val=(N5=-12.9248)

    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
    host_seq.init_run();
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(0));

    host_seq.pre_run();
    REQUIRE(host_bseq->length == 5800);
    REQUIRE(host_seq.start_values[0].f64 == Approx(-0.4));

    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 3);
    REQUIRE(ages[2] == 3);
    REQUIRE(ages[3] == 1);
    REQUIRE(ages[4] == 3);
    REQUIRE(ages[5] == 3);

    // Consts
    REQUIRE(host_seq.values[0].i64 == 200);
    REQUIRE(host_seq.values[1].f64 == 1000);
    // Slots
    REQUIRE(host_seq.values[2].f64 == 100);
    REQUIRE(host_seq.values[3].f64 == 700);
    REQUIRE(host_seq.values[4].f64 == 4244);
    // Global Values
    REQUIRE(host_seq.values[5].i64 == 100);
    REQUIRE(host_seq.values[6].i64 == 800);
    REQUIRE(host_seq.values[7].i64 == 5800);
    // Channels
    REQUIRE(host_seq.values[8].f64 == -0.4);
    // Measures
    REQUIRE(host_seq.values[9].f64 == Approx(-0.52));
    REQUIRE(host_seq.values[10].f64 == -0.4);
    REQUIRE(host_seq.values[11].f64 == Approx(0.08));
    REQUIRE(host_seq.values[12].f64 == Approx(0.48));
    REQUIRE(host_seq.values[13].f64 == Approx(1.084));
    REQUIRE(host_seq.values[14].f64 == Approx(-4.616));
    // Need Orders
    // V15 is a ramp
    REQUIRE(host_seq.values[16].f64 == Approx(-0.52));
    REQUIRE(host_seq.values[17].f64 == Approx(1.084));
    REQUIRE(host_seq.values[18].i64 == 5800);
    REQUIRE(host_seq.values[19].f64 == Approx(-4.616));
    REQUIRE(host_seq.values[20].f64 == Approx(-12.9248));

    REQUIRE(host_bseq->pulses.size() == 6);
    REQUIRE(host_bseq->pulses[0].id == 2);
    REQUIRE(host_bseq->pulses[0].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[0].endvalue].f64 == Approx(-0.52));
    REQUIRE(host_bseq->pulses[1].id == 1);
    REQUIRE(host_bseq->pulses[1].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[1].endvalue].f64 == Approx(0.48));
    REQUIRE(host_bseq->pulses[2].id == 3); // measure
    REQUIRE(host_bseq->pulses[2].chn == 1);
    REQUIRE(host_bseq->pulses[2].len == uint32_t(-2));
    REQUIRE(host_bseq->pulses[3].id == 4);
    REQUIRE(host_bseq->pulses[3].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[3].endvalue].f64 == Approx(1.084));
    REQUIRE(host_bseq->pulses[4].id == 5);
    REQUIRE(host_bseq->pulses[4].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[4].endvalue].f64 == Approx(-4.616));
    REQUIRE(host_bseq->pulses[5].id == 6);
    REQUIRE(host_bseq->pulses[5].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[5].endvalue].f64 == Approx(-12.9248));

    REQUIRE(host_seq.post_run() == 0);
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));

    // Run 4
    //   C0 < G0 < C0 + 1000
    //   G1 < C0 + 1000
    //   G2 < N3
    host_seq.set_global(0, 300);
    host_seq.set_global(2, 3529);
    REQUIRE(ages[0] == 4);
    REQUIRE(ages[1] == 3);
    REQUIRE(ages[2] == 4);
    REQUIRE(ages[3] == 1);
    REQUIRE(ages[4] == 3);
    REQUIRE(ages[5] == 3);
    // start: -0.4
    // P1(id=1): @C0=200 measure=(MV0=-0.4) ramp=[N0] (MV0=-0.4) + t * 0.001 len=1000
    // P2(id=2): @G0=300 measure=(MV1=0.6) val=(N1=0.78)
    // M1(id=3): @G1=1000 measure=(MV2=0.78)
    // P3(id=4): @G1=1000 measure=(MV3=0.78) val=(N2=0.78 * 0.8 + 0.7=1.324)
    // P5(id=6): @G2=6000 measure=(MV5=1.324) val=(N5=3.7072)
    // P4(id=5): @N3=6001 measure=(MV4=3.7072) val=(N4=-1.9928)

    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
    host_seq.init_run();
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(0));

    host_seq.pre_run();
    REQUIRE(host_bseq->length == 6001);
    REQUIRE(host_seq.start_values[0].f64 == Approx(-0.4));

    REQUIRE(ages[0] == 4);
    REQUIRE(ages[1] == 3);
    REQUIRE(ages[2] == 4);
    REQUIRE(ages[3] == 4);
    REQUIRE(ages[4] == 4);
    REQUIRE(ages[5] == 4);

    // Consts
    REQUIRE(host_seq.values[0].i64 == 200);
    REQUIRE(host_seq.values[1].f64 == 1000);
    // Slots
    REQUIRE(host_seq.values[2].f64 == 300);
    REQUIRE(host_seq.values[3].f64 == 700);
    REQUIRE(host_seq.values[4].f64 == 3529);
    // Global Values
    REQUIRE(host_seq.values[5].i64 == 300);
    REQUIRE(host_seq.values[6].i64 == 1000);
    REQUIRE(host_seq.values[7].i64 == 6000);
    // Channels
    REQUIRE(host_seq.values[8].f64 == -0.4);
    // Measures
    REQUIRE(host_seq.values[9].f64 == -0.4);
    REQUIRE(host_seq.values[10].f64 == Approx(0.6));
    REQUIRE(host_seq.values[11].f64 == Approx(0.78));
    REQUIRE(host_seq.values[12].f64 == Approx(0.78));
    REQUIRE(host_seq.values[13].f64 == Approx(3.7072));
    REQUIRE(host_seq.values[14].f64 == Approx(1.324));
    // Need Orders
    // V15 is a ramp
    REQUIRE(host_seq.values[16].f64 == Approx(0.78));
    REQUIRE(host_seq.values[17].f64 == Approx(1.324));
    REQUIRE(host_seq.values[18].i64 == 6001);
    REQUIRE(host_seq.values[19].f64 == Approx(-1.9928));
    REQUIRE(host_seq.values[20].f64 == Approx(3.7072));

    REQUIRE(host_bseq->pulses.size() == 6);
    REQUIRE(host_bseq->pulses[0].id == 1);
    REQUIRE(host_bseq->pulses[0].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[0].endvalue].f64 == Approx(0.6));
    REQUIRE(host_bseq->pulses[1].id == 2);
    REQUIRE(host_bseq->pulses[1].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[1].endvalue].f64 == Approx(0.78));
    REQUIRE(host_bseq->pulses[2].id == 3); // measure
    REQUIRE(host_bseq->pulses[2].chn == 1);
    REQUIRE(host_bseq->pulses[2].len == uint32_t(-2));
    REQUIRE(host_bseq->pulses[3].id == 4);
    REQUIRE(host_bseq->pulses[3].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[3].endvalue].f64 == Approx(1.324));
    REQUIRE(host_bseq->pulses[4].id == 6);
    REQUIRE(host_bseq->pulses[4].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[4].endvalue].f64 == Approx(3.7072));
    REQUIRE(host_bseq->pulses[5].id == 5);
    REQUIRE(host_bseq->pulses[5].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[5].endvalue].f64 == Approx(-1.9928));

    REQUIRE(host_seq.post_run() == 0);
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));

    // Run 5
    //   C0 + 1000 < G0
    //   G1 > C0 + 1000
    //   N3 == G2
    host_seq.set_global(0, 1300);
    host_seq.set_global(2, 4528);
    REQUIRE(ages[0] == 5);
    REQUIRE(ages[1] == 3);
    REQUIRE(ages[2] == 5);
    REQUIRE(ages[3] == 4);
    REQUIRE(ages[4] == 4);
    REQUIRE(ages[5] == 4);
    // start: -0.4
    // P1(id=1): @C0=200 measure=(MV0=-0.4) ramp=[N0] (MV0=-0.4) + t * 0.001 len=1000
    // P2(id=2): @G0=1300 measure=(MV1=0.6) val=(N1=0.78)
    // M1(id=3): @G1=2000 measure=(MV2=0.78)
    // P3(id=4): @G1=2000 measure=(MV3=0.78) val=(N2=0.78 * 0.8 + 0.7=1.324)
    // P4(id=5): @N3=7000 measure=(MV4=1.324) val=(N4=-4.376)
    // P5(id=6): @G2=7000 measure=(MV5=-4.376) val=(N5=-12.2528)

    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
    host_seq.init_run();
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(0));

    host_seq.pre_run();
    REQUIRE(host_bseq->length == 7000);
    REQUIRE(host_seq.start_values[0].f64 == Approx(-0.4));

    REQUIRE(ages[0] == 5);
    REQUIRE(ages[1] == 3);
    REQUIRE(ages[2] == 5);
    REQUIRE(ages[3] == 5);
    REQUIRE(ages[4] == 5);
    REQUIRE(ages[5] == 5);

    // Consts
    REQUIRE(host_seq.values[0].i64 == 200);
    REQUIRE(host_seq.values[1].f64 == 1000);
    // Slots
    REQUIRE(host_seq.values[2].f64 == 1300);
    REQUIRE(host_seq.values[3].f64 == 700);
    REQUIRE(host_seq.values[4].f64 == 4528);
    // Global Values
    REQUIRE(host_seq.values[5].i64 == 1300);
    REQUIRE(host_seq.values[6].i64 == 2000);
    REQUIRE(host_seq.values[7].i64 == 7000);
    // Channels
    REQUIRE(host_seq.values[8].f64 == -0.4);
    // Measures
    REQUIRE(host_seq.values[9].f64 == -0.4);
    REQUIRE(host_seq.values[10].f64 == Approx(0.6));
    REQUIRE(host_seq.values[11].f64 == Approx(0.78));
    REQUIRE(host_seq.values[12].f64 == Approx(0.78));
    REQUIRE(host_seq.values[13].f64 == Approx(1.324));
    REQUIRE(host_seq.values[14].f64 == Approx(-4.376));
    // Need Orders
    // V15 is a ramp
    REQUIRE(host_seq.values[16].f64 == Approx(0.78));
    REQUIRE(host_seq.values[17].f64 == Approx(1.324));
    REQUIRE(host_seq.values[18].i64 == 7000);
    REQUIRE(host_seq.values[19].f64 == Approx(-4.376));
    REQUIRE(host_seq.values[20].f64 == Approx(-12.2528));

    REQUIRE(host_bseq->pulses.size() == 6);
    REQUIRE(host_bseq->pulses[0].id == 1);
    REQUIRE(host_bseq->pulses[0].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[0].endvalue].f64 == Approx(0.6));
    REQUIRE(host_bseq->pulses[1].id == 2);
    REQUIRE(host_bseq->pulses[1].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[1].endvalue].f64 == Approx(0.78));
    REQUIRE(host_bseq->pulses[2].id == 3); // measure
    REQUIRE(host_bseq->pulses[2].chn == 1);
    REQUIRE(host_bseq->pulses[2].len == uint32_t(-2));
    REQUIRE(host_bseq->pulses[3].id == 4);
    REQUIRE(host_bseq->pulses[3].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[3].endvalue].f64 == Approx(1.324));
    REQUIRE(host_bseq->pulses[4].id == 5);
    REQUIRE(host_bseq->pulses[4].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[4].endvalue].f64 == Approx(-4.376));
    REQUIRE(host_bseq->pulses[5].id == 6);
    REQUIRE(host_bseq->pulses[5].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[5].endvalue].f64 == Approx(-12.2528));

    REQUIRE(host_seq.post_run() == 0);
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
}

TEST_CASE("cond") {
    // Test subjects:
    // * Evaluation of condition (true/false)
    // * Measure value (with conditional true/false)
    // * Pulse start value (with conditional true/false)

    // Const:
    //   C0[V0]: time for P1 (100)
    //   C1[V1]: time for P2 (200)
    //   C2[V2]: time for M1 (300)
    //   C3[V3]: time for P3 (400)
    //   C4[V4]: length for P1 (1000)
    //   C5[V5]: -0.1 + t * 0.02
    //   C6[V6]: -0.1 + 1000 * 0.02
    //   C7[V7]: -10
    // Slot:
    //   S0[V8]: condition for P2
    // Global Value:
    // Channel:
    //   CV1[V9]: Channel 1 value
    // BS1:
    //   Values:
    //     MV0[V10]: measure value of M1
    //     MV1[V11]: measure value of M3
    //     N0[V12]: -2 * MV0 + MV1
    //   Pulse:
    //     P1(id=1): @C0=100 ramp=[C5] -0.1 + t * 0.02 len=C4=1000 endvalue=C6=-0.1 + 1000 * 0.02
    //     P2(id=2): @C1=200 val=C7=-10 cond=S0
    //     M1(id=3): @C2=300 measure=MV0
    //     P3(id=4): @C3=400 measure=MV1 val=N0
    //   Branch:
    //     Default: end

    Seq::HostSeq host_seq;
    host_seq.nconsts = 8;
    host_seq.nglobals = 1;
    host_seq.nglobal_vals = 0;
    host_seq.nchannels = 1;
    host_seq.nshared = (host_seq.nconsts + host_seq.nglobals +
                        host_seq.nglobal_vals + host_seq.nchannels);
    host_seq.npublic_globals = 1;

    host_seq.values.resize(13);
    host_seq.values[0].i64 = 100;
    host_seq.values[1].i64 = 200;
    host_seq.values[2].i64 = 300;
    host_seq.values[3].i64 = 400;
    host_seq.values[4].f64 = 1000;
    host_seq.values[5].f64 = 0; // Ramp
    host_seq.values[6].f64 = -0.1 + 1000 * 0.02;
    host_seq.values[7].f64 = -10;
    host_seq.default_values.resize(1);
    host_seq.default_values[0].f64 = -0.2;
    host_seq.types.resize(9); // nconst + nglobals + nglobal_vals
    // Const
    host_seq.types[0] = Seq::HostSeq::Type::Int64;
    host_seq.types[1] = Seq::HostSeq::Type::Int64;
    host_seq.types[2] = Seq::HostSeq::Type::Int64;
    host_seq.types[3] = Seq::HostSeq::Type::Int64;
    host_seq.types[4] = Seq::HostSeq::Type::Float64;
    host_seq.types[5] = Seq::HostSeq::Type::Float64;
    host_seq.types[6] = Seq::HostSeq::Type::Float64;
    host_seq.types[7] = Seq::HostSeq::Type::Float64;
    // Slot
    host_seq.types[8] = Seq::HostSeq::Type::Bool;
    // Global Value

    host_seq.depends.resize(0); // nglobal_vals
    host_seq.global_evals.resize(0); // nglobal_vals

    host_seq.seqs.resize(1);

    {
        auto &host_bseq = host_seq.seqs[0];
        host_bseq.id = 1;

        host_bseq.nmeasure = 2;
        host_bseq.ndirect = 0;
        host_bseq.nneed_order = 1;

        host_bseq.types.resize(1); // ndirect + nneed_order
        host_bseq.types[0] = Seq::HostSeq::Type::Float64;
        host_bseq.evals.resize(1); // ndirect + nneed_order
        host_bseq.evals[0] = [] (void *_p) {
            auto p = (Seq::HostSeq::Value*)_p;
            p[12].f64 = -2 * p[10].f64 + p[11].f64;
        };

        host_bseq.global_refs.resize(0);
        host_bseq.reverse_depends.resize(2); // nmeasure
        host_bseq.reverse_depends[0].push_back(0);
        host_bseq.reverse_depends[1].push_back(0);
        host_bseq.deps_count.resize(1); // nneed_order
        host_bseq.deps_count[0] = 2;

        host_bseq.pulses.resize(4);
        host_bseq.pulses[0] = Seq::HostSeq::Pulse{
            .id = 1, .time = 0, .measure = uint32_t(-1), .len = 4,
            .value = 5, .chn = 1, .ramp_func = [] (double t, void*) {
                return -0.1 + t * 0.02;
            }, .endvalue = 6, .cond = uint32_t(-1)
        };
        host_bseq.pulses[1] = Seq::HostSeq::Pulse{
            .id = 2, .time = 1, .measure = uint32_t(-1), .len = uint32_t(-1),
            .value = 7, .chn = 1, .ramp_func = nullptr, .endvalue = 7, .cond = 8
        };
        host_bseq.pulses[2] = Seq::HostSeq::Pulse{
            .id = 3, .time = 2, .measure = 0, .len = uint32_t(-2),
            .value = uint32_t(-1), .chn = 1, .ramp_func = nullptr,
            .endvalue = uint32_t(-1), .cond = uint32_t(-1)
        };
        host_bseq.pulses[3] = Seq::HostSeq::Pulse{
            .id = 4, .time = 3, .measure = 1, .len = uint32_t(-1),
            .value = 12, .chn = 1, .ramp_func = nullptr, .endvalue = 12, .cond = uint32_t(-1)
        };
        host_bseq.assignments.resize(0);
        host_bseq.assumptions.resize(0);
        host_bseq.branches.resize(0);
        host_bseq.default_branch = -1;
        host_bseq.endtimes.resize(1);
        host_bseq.endtimes[0] = 3;

        host_bseq.ndirect_assumes = 0;
        host_bseq.assumptions_idx.resize(1); // nneed_order
        host_bseq.assumptions_idx[0] = uint32_t(-1);

        host_bseq.cond_global_refs.resize(0);
    }

    host_seq.init();
    auto &ages = host_seq.global_ages();
    REQUIRE(ages.size() == 1);
    REQUIRE(ages[0] == 1);
    auto host_bseq = &host_seq.seqs[0];

    // Run 1
    //   S0 == false
    REQUIRE(host_seq.get_global(0) == 0);
    // start: -0.2
    // P1(id=1): @C0=100 ramp=[C5] -0.1 + t * 0.02 len=C4=1000 endvalue=C6=-0.1 + 1000 * 0.02=19.9
    // P2(id=2): @C1=200 val=C7=-10 cond=S0=false
    // M1(id=3): @C2=300 measure=MV0=-0.1 + 200 * 0.02=3.9
    // P3(id=4): @C3=400 measure=MV1=-0.1 + 1000 * 0.02=19.9 val=-2 * 3.9 + 19.9=12.1

    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
    host_seq.init_run();
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(0));

    host_seq.pre_run();
    REQUIRE(host_bseq->length == 400);
    REQUIRE(host_seq.start_values[0].f64 == Approx(-0.2));

    REQUIRE(ages[0] == 1);

    // Consts
    REQUIRE(host_seq.values[0].i64 == 100);
    REQUIRE(host_seq.values[1].i64 == 200);
    REQUIRE(host_seq.values[2].i64 == 300);
    REQUIRE(host_seq.values[3].i64 == 400);
    REQUIRE(host_seq.values[4].f64 == 1000);
    // REQUIRE(host_seq.values[5].f64 == 0); // Ramp
    REQUIRE(host_seq.values[6].f64 == -0.1 + 1000 * 0.02);
    REQUIRE(host_seq.values[7].f64 == -10);
    // Slots
    REQUIRE(host_seq.values[8].b == false);
    // Global Values
    // Channels
    REQUIRE(host_seq.values[9].f64 == -0.2);
    // Measures
    REQUIRE(host_seq.values[10].f64 == Approx(3.9));
    REQUIRE(host_seq.values[11].f64 == Approx(19.9));
    // Need Orders
    REQUIRE(host_seq.values[12].f64 == Approx(12.1));

    REQUIRE(host_bseq->pulses.size() == 4);
    REQUIRE(host_bseq->pulses[0].id == 1);
    REQUIRE(host_bseq->pulses[0].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[0].endvalue].f64 == Approx(19.9));
    REQUIRE(host_bseq->pulses[1].id == 2);
    REQUIRE(host_bseq->pulses[1].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[1].cond].b == false);
    REQUIRE(host_bseq->pulses[2].id == 3); // measure
    REQUIRE(host_bseq->pulses[2].chn == 1);
    REQUIRE(host_bseq->pulses[2].len == uint32_t(-2));
    REQUIRE(host_bseq->pulses[3].id == 4);
    REQUIRE(host_bseq->pulses[3].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[3].endvalue].f64 == Approx(12.1));

    REQUIRE(host_seq.post_run() == 0);
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));

    // Run 2
    //   S0 == true
    host_seq.set_global(0, 1);
    REQUIRE(ages[0] == 2);
    // start: -0.2
    // P1(id=1): @C0=100 ramp=[C5] -0.1 + t * 0.02 len=C4=1000 endvalue=C6=-0.1 + 1000 * 0.02=19.9
    // P2(id=2): @C1=200 val=C7=-10 cond=S0=true
    // M1(id=3): @C2=300 measure=MV0=-10
    // P3(id=4): @C3=400 measure=MV1=-10 val=-2 * -10 + -10=10

    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
    host_seq.init_run();
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(0));

    host_seq.pre_run();
    REQUIRE(host_bseq->length == 400);
    REQUIRE(host_seq.start_values[0].f64 == Approx(-0.2));

    REQUIRE(ages[0] == 2);

    // Consts
    REQUIRE(host_seq.values[0].i64 == 100);
    REQUIRE(host_seq.values[1].i64 == 200);
    REQUIRE(host_seq.values[2].i64 == 300);
    REQUIRE(host_seq.values[3].i64 == 400);
    REQUIRE(host_seq.values[4].f64 == 1000);
    // REQUIRE(host_seq.values[5].f64 == 0); // Ramp
    REQUIRE(host_seq.values[6].f64 == -0.1 + 1000 * 0.02);
    REQUIRE(host_seq.values[7].f64 == -10);
    // Slots
    REQUIRE(host_seq.values[8].b == true);
    // Global Values
    // Channels
    REQUIRE(host_seq.values[9].f64 == -0.2);
    // Measures
    REQUIRE(host_seq.values[10].f64 == Approx(-10));
    REQUIRE(host_seq.values[11].f64 == Approx(-10));
    // Need Orders
    REQUIRE(host_seq.values[12].f64 == Approx(10));

    REQUIRE(host_bseq->pulses.size() == 4);
    REQUIRE(host_bseq->pulses[0].id == 1);
    REQUIRE(host_bseq->pulses[0].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[0].endvalue].f64 == Approx(19.9));
    REQUIRE(host_bseq->pulses[1].id == 2);
    REQUIRE(host_bseq->pulses[1].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[1].cond].b == true);
    REQUIRE(host_bseq->pulses[2].id == 3); // measure
    REQUIRE(host_bseq->pulses[2].chn == 1);
    REQUIRE(host_bseq->pulses[2].len == uint32_t(-2));
    REQUIRE(host_bseq->pulses[3].id == 4);
    REQUIRE(host_bseq->pulses[3].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[3].endvalue].f64 == Approx(10));

    REQUIRE(host_seq.post_run() == 0);
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
}
