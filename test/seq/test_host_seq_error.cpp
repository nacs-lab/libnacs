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

#include "../error_helper.h"

#include "../../lib/nacs-seq/host_seq.h"
#include "../../lib/nacs-seq/error.h"

#include "../../lib/nacs-utils/number.h"

#include <iostream>
#include <memory>

using namespace NaCs;
using namespace std::literals::string_literals;

TEST_CASE("run_error") {
    Seq::HostSeq host_seq;
    host_seq.nconsts = 0;
    host_seq.nglobals = 0;
    host_seq.nglobal_vals = 0;
    host_seq.nchannels = 0;
    host_seq.nshared = (host_seq.nconsts + host_seq.nglobals +
                        host_seq.nglobal_vals + host_seq.nchannels);
    host_seq.npublic_globals = 0;

    host_seq.values.resize(0);
    host_seq.default_values.resize(0);
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

    auto err = expect_error<std::runtime_error>([&] {
        host_seq.init_run();
    });
    REQUIRE(err.what() == "Sequence must be initialized first."s);
    err = expect_error<std::runtime_error>([&] {
        host_seq.pre_run();
    });
    REQUIRE(err.what() == "Sequence must be initialized first."s);
    err = expect_error<std::runtime_error>([&] {
        host_seq.post_run();
    });
    REQUIRE(err.what() == "Sequence must be initialized first."s);
    err = expect_error<std::runtime_error>([&] {
        host_seq.get_global(0);
    });
    REQUIRE(err.what() == "Sequence must be initialized first."s);
    err = expect_error<std::runtime_error>([&] {
        host_seq.set_global(0, 0);
    });
    REQUIRE(err.what() == "Sequence must be initialized first."s);

    host_seq.init();
    err = expect_error<std::runtime_error>([&] {
        host_seq.init();
    });
    REQUIRE(err.what() == "Sequence already initialized."s);
    err = expect_error<std::runtime_error>([&] {
        host_seq.pre_run();
    });
    REQUIRE(err.what() == "Sequence must be started before running."s);
    err = expect_error<std::runtime_error>([&] {
        host_seq.post_run();
    });
    REQUIRE(err.what() == "Sequence must be started before running."s);

    for (uint32_t i = 0; i < 3; i++) {
        REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
        host_seq.init_run();
        REQUIRE(host_seq.cur_seq_idx() == 0);
        auto err = expect_error<std::runtime_error>([&] {
            host_seq.init();
        });
        REQUIRE(err.what() == "Sequence already initialized."s);
        err = expect_error<std::runtime_error>([&] {
            host_seq.init_run();
        });
        REQUIRE(err.what() == "Unfinished sequence cannot be restarted again."s);
        err = expect_error<std::runtime_error>([&] {
            host_seq.post_run();
        });
        REQUIRE(err.what() == "Sequence not running."s);


        host_seq.pre_run();
        REQUIRE(host_seq.cur_seq_idx() == 0);
        REQUIRE(host_seq.first_bseq);
        err = expect_error<std::runtime_error>([&] {
            host_seq.init();
        });
        REQUIRE(err.what() == "Sequence already initialized."s);
        err = expect_error<std::runtime_error>([&] {
            host_seq.init_run();
        });
        REQUIRE(err.what() == "Unfinished sequence cannot be restarted again."s);
        err = expect_error<std::runtime_error>([&] {
            host_seq.pre_run();
        });
        REQUIRE(err.what() == "Sequence already running."s);

        REQUIRE(host_seq.post_run() == 0);
        REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
        err = expect_error<std::runtime_error>([&] {
            host_seq.init();
        });
        REQUIRE(err.what() == "Sequence already initialized."s);
        err = expect_error<std::runtime_error>([&] {
            host_seq.pre_run();
        });
        REQUIRE(err.what() == "Sequence must be started before running."s);
        err = expect_error<std::runtime_error>([&] {
            host_seq.post_run();
        });
        REQUIRE(err.what() == "Sequence must be started before running."s);
    }
}

TEST_CASE("neg_const_pulse_time") {
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

        host_bseq.pulses.resize(1);
        host_bseq.pulses[0] = Seq::HostSeq::Pulse{
            .id = 17, .time = 0, .measure = uint32_t(-1), .len = uint32_t(-1),
            .value = 1, .chn = 1, .ramp_func = nullptr, .endvalue = 1, .cond = uint32_t(-1)
        };
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

    auto err = expect_error<Seq::Error>([&] {
        host_seq.pre_run();
    });
    REQUIRE(err.type == Seq::Error::Type::Pulse);
    REQUIRE(err.code == uint16_t(Seq::Error::Pulse::NegTime));
    REQUIRE(err.type1 == Seq::Error::Type::Pulse);
    REQUIRE(err.id1 == 17);
    REQUIRE(err.what() == "Pulse time must be positive."s);
}

TEST_CASE("neg_const_measure_time") {
    Seq::HostSeq host_seq;
    host_seq.nconsts = 2;
    host_seq.nglobals = 0;
    host_seq.nglobal_vals = 0;
    host_seq.nchannels = 1;
    host_seq.nshared = (host_seq.nconsts + host_seq.nglobals +
                        host_seq.nglobal_vals + host_seq.nchannels);
    host_seq.npublic_globals = 0;

    host_seq.values.resize(5);
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

        host_bseq.nmeasure = 1;
        host_bseq.ndirect = 0;
        host_bseq.nneed_order = 0;

        host_bseq.types.resize(0); // ndirect + nneed_order
        host_bseq.evals.resize(0); // ndirect + nneed_order

        host_bseq.global_refs.resize(0);
        host_bseq.reverse_depends.resize(1); // nmeasure
        host_bseq.deps_count.resize(0); // nneed_order

        host_bseq.pulses.resize(1);
        host_bseq.pulses[0] = Seq::HostSeq::Pulse{
            .id = 13, .time = 0, .measure = 0, .len = uint32_t(-2),
            .value = 1, .chn = 1, .ramp_func = nullptr, .endvalue = 1, .cond = uint32_t(-1)
        };
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

    auto err = expect_error<Seq::Error>([&] {
        host_seq.pre_run();
    });
    REQUIRE(err.type == Seq::Error::Type::Pulse);
    REQUIRE(err.code == uint16_t(Seq::Error::Pulse::NegTime));
    REQUIRE(err.type1 == Seq::Error::Type::Pulse);
    REQUIRE(err.id1 == 13);
    REQUIRE(err.what() == "Pulse time must be positive."s);
}

TEST_CASE("neg_computed_pulse_time") {
    Seq::HostSeq host_seq;
    host_seq.nconsts = 1;
    host_seq.nglobals = 1;
    host_seq.nglobal_vals = 1;
    host_seq.nchannels = 1;
    host_seq.nshared = (host_seq.nconsts + host_seq.nglobals +
                        host_seq.nglobal_vals + host_seq.nchannels);
    host_seq.npublic_globals = 1;

    host_seq.values.resize(5);
    host_seq.values[0].f64 = 2.3;
    host_seq.default_values.resize(1);
    host_seq.default_values[0].f64 = 1.2;
    host_seq.types.resize(3); // nconst + nglobals + nglobal_vals
    host_seq.types[0] = Seq::HostSeq::Type::Float64;
    host_seq.types[1] = Seq::HostSeq::Type::Float64;
    host_seq.types[2] = Seq::HostSeq::Type::Int64;

    host_seq.depends.resize(1); // nglobal_vals
    host_seq.depends[0].push_back(0); // nglobal_vals
    host_seq.global_evals.resize(1); // nglobal_vals
    host_seq.global_evals[0] = [] (void *_p) {
        auto p = (Seq::HostSeq::Value*)_p;
        p[2].i64 = round<int64_t>(p[1].f64);
    };

    host_seq.seqs.resize(1);

    {
        auto &host_bseq = host_seq.seqs[0];
        host_bseq.id = 1;

        host_bseq.nmeasure = 0;
        host_bseq.ndirect = 0;
        host_bseq.nneed_order = 0;

        host_bseq.types.resize(0); // ndirect + nneed_order
        host_bseq.evals.resize(0); // ndirect + nneed_order

        host_bseq.global_refs.resize(1);
        host_bseq.global_refs[0] = 0;
        host_bseq.reverse_depends.resize(0); // nmeasure
        host_bseq.deps_count.resize(0); // nneed_order

        host_bseq.pulses.resize(1);
        host_bseq.pulses[0] = Seq::HostSeq::Pulse{
            .id = 19, .time = 2, .measure = uint32_t(-1), .len = uint32_t(-1),
            .value = 0, .chn = 1, .ramp_func = nullptr, .endvalue = 0, .cond = uint32_t(-1)
        };
        host_bseq.assignments.resize(0);
        host_bseq.assumptions.resize(0);
        host_bseq.branches.resize(0);
        host_bseq.default_branch = -1;
        host_bseq.endtimes.resize(1);
        host_bseq.endtimes[0] = 2;

        host_bseq.ndirect_assumes = 0;
        host_bseq.assumptions_idx.resize(0); // nneed_order

        host_bseq.cond_global_refs.resize(0);
    }

    host_seq.init();
    REQUIRE(host_seq.global_ages().size() == 2);
    host_seq.set_global(0, -0.6);

    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
    host_seq.init_run();
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(0));

    auto err = expect_error<Seq::Error>([&] {
        host_seq.pre_run();
    });
    REQUIRE(err.type == Seq::Error::Type::Pulse);
    REQUIRE(err.code == uint16_t(Seq::Error::Pulse::NegTime));
    REQUIRE(err.type1 == Seq::Error::Type::Pulse);
    REQUIRE(err.id1 == 19);
    REQUIRE(err.what() == "Pulse time must be positive."s);
}

TEST_CASE("neg_computed_measure_time") {
    Seq::HostSeq host_seq;
    host_seq.nconsts = 1;
    host_seq.nglobals = 1;
    host_seq.nglobal_vals = 1;
    host_seq.nchannels = 1;
    host_seq.nshared = (host_seq.nconsts + host_seq.nglobals +
                        host_seq.nglobal_vals + host_seq.nchannels);
    host_seq.npublic_globals = 1;

    host_seq.values.resize(6);
    host_seq.values[0].f64 = 2.3;
    host_seq.default_values.resize(1);
    host_seq.default_values[0].f64 = 1.2;
    host_seq.types.resize(3); // nconst + nglobals + nglobal_vals
    host_seq.types[0] = Seq::HostSeq::Type::Float64;
    host_seq.types[1] = Seq::HostSeq::Type::Float64;
    host_seq.types[2] = Seq::HostSeq::Type::Int64;

    host_seq.depends.resize(1); // nglobal_vals
    host_seq.depends[0].push_back(0); // nglobal_vals
    host_seq.global_evals.resize(1); // nglobal_vals
    host_seq.global_evals[0] = [] (void *_p) {
        auto p = (Seq::HostSeq::Value*)_p;
        p[2].i64 = round<int64_t>(p[1].f64);
    };

    host_seq.seqs.resize(1);

    {
        auto &host_bseq = host_seq.seqs[0];
        host_bseq.id = 1;

        host_bseq.nmeasure = 1;
        host_bseq.ndirect = 0;
        host_bseq.nneed_order = 0;

        host_bseq.types.resize(0); // ndirect + nneed_order
        host_bseq.evals.resize(0); // ndirect + nneed_order

        host_bseq.global_refs.resize(1);
        host_bseq.global_refs[0] = 0;
        host_bseq.reverse_depends.resize(1); // nmeasure
        host_bseq.deps_count.resize(0); // nneed_order

        host_bseq.pulses.resize(1);
        host_bseq.pulses[0] = Seq::HostSeq::Pulse{
            .id = 11, .time = 2, .measure = 0, .len = uint32_t(-2),
            .value = 0, .chn = 1, .ramp_func = nullptr, .endvalue = 0, .cond = uint32_t(-1)
        };
        host_bseq.assignments.resize(0);
        host_bseq.assumptions.resize(0);
        host_bseq.branches.resize(0);
        host_bseq.default_branch = -1;
        host_bseq.endtimes.resize(1);
        host_bseq.endtimes[0] = 2;

        host_bseq.ndirect_assumes = 0;
        host_bseq.assumptions_idx.resize(0); // nneed_order

        host_bseq.cond_global_refs.resize(0);
    }

    host_seq.init();
    REQUIRE(host_seq.global_ages().size() == 2);
    host_seq.set_global(0, -0.6);

    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
    host_seq.init_run();
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(0));

    auto err = expect_error<Seq::Error>([&] {
        host_seq.pre_run();
    });
    REQUIRE(err.type == Seq::Error::Type::Pulse);
    REQUIRE(err.code == uint16_t(Seq::Error::Pulse::NegTime));
    REQUIRE(err.type1 == Seq::Error::Type::Pulse);
    REQUIRE(err.id1 == 11);
    REQUIRE(err.what() == "Pulse time must be positive."s);
}

TEST_CASE("direct_assumption") {
    Seq::HostSeq host_seq;
    host_seq.nconsts = 0;
    host_seq.nglobals = 1;
    host_seq.nglobal_vals = 1;
    host_seq.nchannels = 0;
    host_seq.nshared = (host_seq.nconsts + host_seq.nglobals +
                        host_seq.nglobal_vals + host_seq.nchannels);
    host_seq.npublic_globals = 1;

    host_seq.values.resize(2);
    host_seq.default_values.resize(0);
    host_seq.types.resize(2); // nconst + nglobals + nglobal_vals
    host_seq.types[0] = Seq::HostSeq::Type::Float64;
    host_seq.types[1] = Seq::HostSeq::Type::Int64;

    host_seq.depends.resize(1); // nglobal_vals
    host_seq.depends[0].push_back(0); // nglobal_vals
    host_seq.global_evals.resize(1); // nglobal_vals
    host_seq.global_evals[0] = [] (void *_p) {
        auto p = (Seq::HostSeq::Value*)_p;
        p[1].i64 = round<int64_t>(p[0].f64);
    };

    host_seq.seqs.resize(1);

    {
        auto &host_bseq = host_seq.seqs[0];
        host_bseq.id = 1;

        host_bseq.nmeasure = 0;
        host_bseq.ndirect = 0;
        host_bseq.nneed_order = 0;

        host_bseq.types.resize(0); // ndirect + nneed_order
        host_bseq.evals.resize(0); // ndirect + nneed_order

        host_bseq.global_refs.resize(1);
        host_bseq.global_refs[0] = 0;
        host_bseq.reverse_depends.resize(0); // nmeasure
        host_bseq.deps_count.resize(0); // nneed_order

        host_bseq.pulses.resize(0);
        host_bseq.assignments.resize(0);
        host_bseq.assumptions.resize(1);
        host_bseq.assumptions[0] = Seq::HostSeq::Assumption{
            .sign = Seq::Sign::Pos, .value = 1, .id = 23
        };
        host_bseq.branches.resize(0);
        host_bseq.default_branch = -1;
        host_bseq.endtimes.resize(0);

        host_bseq.ndirect_assumes = 1;
        host_bseq.assumptions_idx.resize(0); // nneed_order

        host_bseq.cond_global_refs.resize(0);
    }

    host_seq.init();
    REQUIRE(host_seq.global_ages().size() == 2);
    host_seq.set_global(0, 0.4);

    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
    host_seq.init_run();
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(0));

    auto err = expect_error<Seq::Error>([&] {
        host_seq.pre_run();
    });
    REQUIRE(err.type == Seq::Error::Type::EventTime);
    REQUIRE(err.code == uint16_t(Seq::Error::EventTime::NonPosTime));
    REQUIRE(err.type1 == Seq::Error::Type::EventTime);
    REQUIRE(err.id1 == 23);
    REQUIRE(err.what() == "Positive time expected."s);
}

TEST_CASE("need_order_assumption") {
    // Const:
    //   C0[V0]: time 1 (20)
    //   C1[V1]: val 1 (1.2)
    //   C2[V2]: val 2 (2.3)
    // Slot:
    //   S0[V3]: for time 2
    // Global Value:
    //   G0[V4]: round S0 (time 2)
    //   G1[V5]: C0 + G0
    // Channel:
    //   CV1[V6]: Channel 1 value
    // BS1:
    //   Values:
    //     MV0[V7]: measure value of M1
    //     N0[V8]: round((1.82 - MV0) * 10)
    //   Pulse:
    //     P1: @C0 val=C1
    //     P2: @G0 val=C2
    //     M1: @G1 measure=MV0
    //   Assumption:
    //     N0 > 0
    //   Branch:
    //     Default: end
    Seq::HostSeq host_seq;
    host_seq.nconsts = 3;
    host_seq.nglobals = 1;
    host_seq.nglobal_vals = 2;
    host_seq.nchannels = 1;
    host_seq.nshared = (host_seq.nconsts + host_seq.nglobals +
                        host_seq.nglobal_vals + host_seq.nchannels);
    host_seq.npublic_globals = 1;

    host_seq.values.resize(9);
    host_seq.values[0].i64 = 20;
    host_seq.values[1].f64 = 1.2;
    host_seq.values[2].f64 = 2.3;
    host_seq.default_values.resize(1);
    host_seq.default_values[0].f64 = 0.2;
    host_seq.types.resize(6); // nconst + nglobals + nglobal_vals
    // Const
    host_seq.types[0] = Seq::HostSeq::Type::Int64;
    host_seq.types[1] = Seq::HostSeq::Type::Float64;
    host_seq.types[2] = Seq::HostSeq::Type::Float64;
    // Slot
    host_seq.types[3] = Seq::HostSeq::Type::Float64;
    // Global Value
    host_seq.types[4] = Seq::HostSeq::Type::Int64;
    host_seq.types[5] = Seq::HostSeq::Type::Int64;

    host_seq.depends.resize(2); // nglobal_vals
    host_seq.depends[0].push_back(0);
    host_seq.depends[1].push_back(0);
    host_seq.global_evals.resize(2); // nglobal_vals
    host_seq.global_evals[0] = [] (void *_p) {
        auto p = (Seq::HostSeq::Value*)_p;
        p[4].i64 = round<int64_t>(p[3].f64);
    };
    host_seq.global_evals[1] = [] (void *_p) {
        auto p = (Seq::HostSeq::Value*)_p;
        p[5].i64 = p[4].i64 + 20;
    };

    host_seq.seqs.resize(1);

    {
        auto &host_bseq = host_seq.seqs[0];
        host_bseq.id = 1;

        host_bseq.nmeasure = 1;
        host_bseq.ndirect = 0;
        host_bseq.nneed_order = 1;

        host_bseq.types.resize(1); // ndirect + nneed_order
        host_bseq.types[0] = Seq::HostSeq::Type::Int64;
        host_bseq.evals.resize(1); // ndirect + nneed_order
        host_bseq.evals[0] = [] (void *_p) {
            auto p = (Seq::HostSeq::Value*)_p;
            p[8].i64 = round<int64_t>(10 * (1.82 - p[7].f64));
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
            .id = 2, .time = 4, .measure = uint32_t(-1), .len = uint32_t(-1),
            .value = 2, .chn = 1, .ramp_func = nullptr, .endvalue = 2, .cond = uint32_t(-1)
        };
        host_bseq.pulses[2] = Seq::HostSeq::Pulse{
            .id = 3, .time = 5, .measure = 0, .len = uint32_t(-2),
            .value = uint32_t(-1), .chn = 1, .ramp_func = nullptr,
            .endvalue = uint32_t(-1), .cond = uint32_t(-1)
        };
        host_bseq.assignments.resize(0);
        host_bseq.assumptions.resize(1);
        host_bseq.assumptions[0] = Seq::HostSeq::Assumption{
            .sign = Seq::Sign::NonNeg, .value = 8, .id = 12
        };
        host_bseq.branches.resize(0);
        host_bseq.default_branch = -1;
        host_bseq.endtimes.resize(1);
        host_bseq.endtimes[0] = 5;

        host_bseq.ndirect_assumes = 0;
        host_bseq.assumptions_idx.resize(1); // nneed_order
        host_bseq.assumptions_idx[0] = 0;

        host_bseq.cond_global_refs.resize(0);
    }

    host_seq.init();
    auto &ages = host_seq.global_ages();
    REQUIRE(ages.size() == 3);
    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 0);
    REQUIRE(ages[2] == 0);
    host_seq.set_global(0, 0.4); // this will make G0 == 0
    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 0);
    REQUIRE(ages[2] == 0);

    // Run 1
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
    host_seq.init_run();
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(0));

    host_seq.pre_run();
    auto host_bseq = &host_seq.seqs[0];
    REQUIRE(host_bseq->length == 20);
    REQUIRE(host_seq.start_values[0].f64 == Approx(2.3));

    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 1);
    REQUIRE(ages[2] == 1);

    REQUIRE(host_seq.values[0].i64 == 20);
    REQUIRE(host_seq.values[1].f64 == 1.2);
    REQUIRE(host_seq.values[2].f64 == 2.3);
    REQUIRE(host_seq.values[3].f64 == 0.4);
    REQUIRE(host_seq.values[4].i64 == 0);
    REQUIRE(host_seq.values[5].i64 == 20);
    REQUIRE(host_seq.values[6].f64 == 0.2);
    REQUIRE(host_seq.values[7].f64 == 1.2);
    REQUIRE(host_seq.values[8].i64 == 6);

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

    REQUIRE(ages[0] == 2);
    REQUIRE(ages[1] == 1);
    REQUIRE(ages[2] == 1);

    auto err = expect_error<Seq::Error>([&] {
        host_seq.pre_run();
    });
    REQUIRE(err.type == Seq::Error::Type::EventTime);
    REQUIRE(err.code == uint16_t(Seq::Error::EventTime::NegTime));
    REQUIRE(err.type1 == Seq::Error::Type::EventTime);
    REQUIRE(err.id1 == 12);
    REQUIRE(err.what() == "Non-negative time expected."s);

    REQUIRE(host_seq.values[8].i64 == -5);
}

TEST_CASE("deps_cycle") {
    // Const:
    // Slot:
    // Global Value:
    // Channel:
    //   CV1[V0]: Channel 1 value
    // BS1:
    //   Values:
    //     MV0[V1]: measure value of M1
    //     N0[V2]: round(MV0)
    //   Pulse:
    //     M1: @N0 measure=MV0
    //   Branch:
    //     Default: end
    Seq::HostSeq host_seq;
    host_seq.nconsts = 0;
    host_seq.nglobals = 0;
    host_seq.nglobal_vals = 0;
    host_seq.nchannels = 1;
    host_seq.nshared = (host_seq.nconsts + host_seq.nglobals +
                        host_seq.nglobal_vals + host_seq.nchannels);
    host_seq.npublic_globals = 0;

    host_seq.values.resize(3);
    host_seq.default_values.resize(1);
    host_seq.default_values[0].f64 = 0.2;
    host_seq.types.resize(0); // nconst + nglobals + nglobal_vals

    host_seq.depends.resize(0); // nglobal_vals
    host_seq.global_evals.resize(0); // nglobal_vals

    host_seq.seqs.resize(1);

    {
        auto &host_bseq = host_seq.seqs[0];
        host_bseq.id = 1;

        host_bseq.nmeasure = 1;
        host_bseq.ndirect = 0;
        host_bseq.nneed_order = 1;

        host_bseq.types.resize(1); // ndirect + nneed_order
        host_bseq.types[0] = Seq::HostSeq::Type::Int64;
        host_bseq.evals.resize(1); // ndirect + nneed_order
        host_bseq.evals[0] = [] (void *_p) {
            auto p = (Seq::HostSeq::Value*)_p;
            p[2].i64 = round<int64_t>(p[1].f64);
        };

        host_bseq.global_refs.resize(0);
        host_bseq.reverse_depends.resize(1); // nmeasure
        host_bseq.reverse_depends[0].push_back(0);
        host_bseq.deps_count.resize(1); // nneed_order
        host_bseq.deps_count[0] = 1;

        host_bseq.pulses.resize(1);
        host_bseq.pulses[0] = Seq::HostSeq::Pulse{
            .id = 3, .time = 2, .measure = 0, .len = uint32_t(-2),
            .value = uint32_t(-1), .chn = 1, .ramp_func = nullptr,
            .endvalue = uint32_t(-1), .cond = uint32_t(-1)
        };
        host_bseq.assignments.resize(0);
        host_bseq.assumptions.resize(0);
        host_bseq.branches.resize(0);
        host_bseq.default_branch = -1;
        host_bseq.endtimes.resize(1);
        host_bseq.endtimes[0] = 2;

        host_bseq.ndirect_assumes = 0;
        host_bseq.assumptions_idx.resize(1); // nneed_order
        host_bseq.assumptions_idx[0] = uint32_t(-1);

        host_bseq.cond_global_refs.resize(0);
    }

    host_seq.init();
    REQUIRE(host_seq.global_ages().empty());

    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
    host_seq.init_run();
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(0));

    auto err = expect_error<std::runtime_error>([&] {
        host_seq.pre_run();
    });
    REQUIRE(err.what() == "Cannot compute sequence."s);
}

TEST_CASE("backward_time") {
    // In order to trigger a go back in time error,
    // we need a pulse time that depends on a measure.

    // Const:
    //   C0[V0]: time 1 (20)
    //   C1[V1]: val 1 (1.83)
    //   C2[V2]: val 2 (2.3)
    // Slot:
    //   S0[V3]: for time 2
    // Global Value:
    //   G0[V4]: round S0 (time 2)
    //   G1[V5]: C0 + G0
    // Channel:
    //   CV1[V6]: Channel 1 value
    // BS1:
    //   Values:
    //     MV0[V7]: measure value of M1
    //     N0[V8]: G1 + round((1.82 - MV0) * 10)
    //   Pulse:
    //     P1: @C0 val=C1
    //     P2: @G0 val=C2
    //     M1: @G1 measure=MV0
    //     P3: @N0 val=S0
    //   Branch:
    //     Default: end
    Seq::HostSeq host_seq;
    host_seq.nconsts = 3;
    host_seq.nglobals = 1;
    host_seq.nglobal_vals = 2;
    host_seq.nchannels = 1;
    host_seq.nshared = (host_seq.nconsts + host_seq.nglobals +
                        host_seq.nglobal_vals + host_seq.nchannels);
    host_seq.npublic_globals = 1;

    host_seq.values.resize(9);
    host_seq.values[0].i64 = 20;
    host_seq.values[1].f64 = 1.83;
    host_seq.values[2].f64 = 2.3;
    host_seq.default_values.resize(1);
    host_seq.default_values[0].f64 = 0.2;
    host_seq.types.resize(6); // nconst + nglobals + nglobal_vals
    // Const
    host_seq.types[0] = Seq::HostSeq::Type::Int64;
    host_seq.types[1] = Seq::HostSeq::Type::Float64;
    host_seq.types[2] = Seq::HostSeq::Type::Float64;
    // Slot
    host_seq.types[3] = Seq::HostSeq::Type::Float64;
    // Global Value
    host_seq.types[4] = Seq::HostSeq::Type::Int64;
    host_seq.types[5] = Seq::HostSeq::Type::Int64;

    host_seq.depends.resize(2); // nglobal_vals
    host_seq.depends[0].push_back(0);
    host_seq.depends[1].push_back(0);
    host_seq.global_evals.resize(2); // nglobal_vals
    host_seq.global_evals[0] = [] (void *_p) {
        auto p = (Seq::HostSeq::Value*)_p;
        p[4].i64 = round<int64_t>(p[3].f64);
    };
    host_seq.global_evals[1] = [] (void *_p) {
        auto p = (Seq::HostSeq::Value*)_p;
        p[5].i64 = p[4].i64 + 20;
    };

    host_seq.seqs.resize(1);

    {
        auto &host_bseq = host_seq.seqs[0];
        host_bseq.id = 1;

        host_bseq.nmeasure = 1;
        host_bseq.ndirect = 0;
        host_bseq.nneed_order = 1;

        host_bseq.types.resize(1); // ndirect + nneed_order
        host_bseq.types[0] = Seq::HostSeq::Type::Int64;
        host_bseq.evals.resize(1); // ndirect + nneed_order
        host_bseq.evals[0] = [] (void *_p) {
            auto p = (Seq::HostSeq::Value*)_p;
            p[8].i64 = p[5].i64 + round<int64_t>(10 * (1.82 - p[7].f64));
        };

        host_bseq.global_refs.resize(2);
        host_bseq.global_refs[0] = 0;
        host_bseq.global_refs[1] = 1;
        host_bseq.reverse_depends.resize(1); // nmeasure
        host_bseq.reverse_depends[0].push_back(0);
        host_bseq.deps_count.resize(1); // nneed_order
        host_bseq.deps_count[0] = 1;

        host_bseq.pulses.resize(4);
        host_bseq.pulses[0] = Seq::HostSeq::Pulse{
            .id = 1, .time = 0, .measure = uint32_t(-1), .len = uint32_t(-1),
            .value = 1, .chn = 1, .ramp_func = nullptr, .endvalue = 1, .cond = uint32_t(-1)
        };
        host_bseq.pulses[1] = Seq::HostSeq::Pulse{
            .id = 2, .time = 4, .measure = uint32_t(-1), .len = uint32_t(-1),
            .value = 2, .chn = 1, .ramp_func = nullptr, .endvalue = 2, .cond = uint32_t(-1)
        };
        host_bseq.pulses[2] = Seq::HostSeq::Pulse{
            .id = 3, .time = 5, .measure = 0, .len = uint32_t(-2),
            .value = uint32_t(-1), .chn = 1, .ramp_func = nullptr,
            .endvalue = uint32_t(-1), .cond = uint32_t(-1)
        };
        host_bseq.pulses[3] = Seq::HostSeq::Pulse{
            .id = 7, .time = 8, .measure = uint32_t(-1), .len = uint32_t(-1),
            .value = 3, .chn = 1, .ramp_func = nullptr, .endvalue = 3, .cond = uint32_t(-1)
        };
        host_bseq.assignments.resize(0);
        host_bseq.assumptions.resize(0);
        host_bseq.branches.resize(0);
        host_bseq.default_branch = -1;
        host_bseq.endtimes.resize(1);
        host_bseq.endtimes[0] = 8;

        host_bseq.ndirect_assumes = 0;
        host_bseq.assumptions_idx.resize(1); // nneed_order
        host_bseq.assumptions_idx[0] = uint32_t(-1);

        host_bseq.cond_global_refs.resize(0);
    }

    host_seq.init();
    auto &ages = host_seq.global_ages();
    REQUIRE(ages.size() == 3);
    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 0);
    REQUIRE(ages[2] == 0);
    host_seq.set_global(0, 0.4); // this will make G0 == 0
    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 0);
    REQUIRE(ages[2] == 0);

    // Run 1
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
    host_seq.init_run();
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(0));

    host_seq.pre_run();
    auto host_bseq = &host_seq.seqs[0];
    REQUIRE(host_bseq->length == 20);
    REQUIRE(host_seq.start_values[0].f64 == Approx(2.3));

    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 1);
    REQUIRE(ages[2] == 1);

    REQUIRE(host_seq.values[0].i64 == 20);
    REQUIRE(host_seq.values[1].f64 == 1.83);
    REQUIRE(host_seq.values[2].f64 == 2.3);
    REQUIRE(host_seq.values[3].f64 == 0.4);
    REQUIRE(host_seq.values[4].i64 == 0);
    REQUIRE(host_seq.values[5].i64 == 20);
    REQUIRE(host_seq.values[6].f64 == 0.2);
    REQUIRE(host_seq.values[7].f64 == 1.83);
    REQUIRE(host_seq.values[8].i64 == 20);

    REQUIRE(host_bseq->pulses.size() == 4);
    REQUIRE(host_bseq->pulses[0].id == 2);
    REQUIRE(host_bseq->pulses[0].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[0].endvalue].f64 == 2.3);
    REQUIRE(host_bseq->pulses[1].id == 1);
    REQUIRE(host_bseq->pulses[1].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[1].endvalue].f64 == 1.83);
    REQUIRE(host_bseq->pulses[2].id == 3); // measure
    REQUIRE(host_bseq->pulses[2].chn == 1);
    REQUIRE(host_bseq->pulses[2].len == uint32_t(-2));
    REQUIRE(host_bseq->pulses[3].id == 7);
    REQUIRE(host_bseq->pulses[3].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[3].endvalue].f64 == 0.4);

    REQUIRE(host_seq.post_run() == 0);
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));

    // Run 2
    host_seq.init_run();
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(0));
    host_seq.set_global(0, 25.6);

    REQUIRE(ages[0] == 2);
    REQUIRE(ages[1] == 1);
    REQUIRE(ages[2] == 1);

    auto err = expect_error<std::runtime_error>([&] {
        host_seq.pre_run();
    });
    REQUIRE(err.what() == "Going back in time not allowed."s);

    REQUIRE(host_seq.values[8].i64 == 41);
}

TEST_CASE("backward_time_id") {
    // In order to trigger a go back in time error,
    // we need a pulse time that depends on a measure.

    // Const:
    //   C0[V0]: time 1 (20)
    //   C1[V1]: val 1 (1.2)
    //   C2[V2]: val 2 (1.81)
    // Slot:
    //   S0[V3]: for time 2
    // Global Value:
    //   G0[V4]: round S0 (time 2)
    //   G1[V5]: C0 + G0
    // Channel:
    //   CV1[V6]: Channel 1 value
    // BS1:
    //   Values:
    //     MV0[V7]: measure value of M1
    //     N0[V8]: G1 + round((1.82 - MV0) * 10)
    //   Pulse:
    //     P1: @C0 val=C1
    //     P2: @G0 val=C2
    //     M1: @G1 measure=MV0
    //     P3: @N0 val=S0
    //   Branch:
    //     Default: end
    Seq::HostSeq host_seq;
    host_seq.nconsts = 3;
    host_seq.nglobals = 1;
    host_seq.nglobal_vals = 2;
    host_seq.nchannels = 1;
    host_seq.nshared = (host_seq.nconsts + host_seq.nglobals +
                        host_seq.nglobal_vals + host_seq.nchannels);
    host_seq.npublic_globals = 1;

    host_seq.values.resize(9);
    host_seq.values[0].i64 = 20;
    host_seq.values[1].f64 = 1.2;
    host_seq.values[2].f64 = 1.81;
    host_seq.default_values.resize(1);
    host_seq.default_values[0].f64 = 0.2;
    host_seq.types.resize(6); // nconst + nglobals + nglobal_vals
    // Const
    host_seq.types[0] = Seq::HostSeq::Type::Int64;
    host_seq.types[1] = Seq::HostSeq::Type::Float64;
    host_seq.types[2] = Seq::HostSeq::Type::Float64;
    // Slot
    host_seq.types[3] = Seq::HostSeq::Type::Float64;
    // Global Value
    host_seq.types[4] = Seq::HostSeq::Type::Int64;
    host_seq.types[5] = Seq::HostSeq::Type::Int64;

    host_seq.depends.resize(2); // nglobal_vals
    host_seq.depends[0].push_back(0);
    host_seq.depends[1].push_back(0);
    host_seq.global_evals.resize(2); // nglobal_vals
    host_seq.global_evals[0] = [] (void *_p) {
        auto p = (Seq::HostSeq::Value*)_p;
        p[4].i64 = round<int64_t>(p[3].f64);
    };
    host_seq.global_evals[1] = [] (void *_p) {
        auto p = (Seq::HostSeq::Value*)_p;
        p[5].i64 = p[4].i64 + 20;
    };

    host_seq.seqs.resize(1);

    {
        auto &host_bseq = host_seq.seqs[0];
        host_bseq.id = 1;

        host_bseq.nmeasure = 1;
        host_bseq.ndirect = 0;
        host_bseq.nneed_order = 1;

        host_bseq.types.resize(1); // ndirect + nneed_order
        host_bseq.types[0] = Seq::HostSeq::Type::Int64;
        host_bseq.evals.resize(1); // ndirect + nneed_order
        host_bseq.evals[0] = [] (void *_p) {
            auto p = (Seq::HostSeq::Value*)_p;
            p[8].i64 = p[5].i64 + round<int64_t>(10 * (1.82 - p[7].f64));
        };

        host_bseq.global_refs.resize(2);
        host_bseq.global_refs[0] = 0;
        host_bseq.global_refs[1] = 1;
        host_bseq.reverse_depends.resize(1); // nmeasure
        host_bseq.reverse_depends[0].push_back(0);
        host_bseq.deps_count.resize(1); // nneed_order
        host_bseq.deps_count[0] = 1;

        host_bseq.pulses.resize(4);
        host_bseq.pulses[0] = Seq::HostSeq::Pulse{
            .id = 1, .time = 0, .measure = uint32_t(-1), .len = uint32_t(-1),
            .value = 1, .chn = 1, .ramp_func = nullptr, .endvalue = 1, .cond = uint32_t(-1)
        };
        host_bseq.pulses[1] = Seq::HostSeq::Pulse{
            .id = 2, .time = 4, .measure = uint32_t(-1), .len = uint32_t(-1),
            .value = 2, .chn = 1, .ramp_func = nullptr, .endvalue = 2, .cond = uint32_t(-1)
        };
        host_bseq.pulses[2] = Seq::HostSeq::Pulse{
            .id = 4, .time = 5, .measure = 0, .len = uint32_t(-2),
            .value = uint32_t(-1), .chn = 1, .ramp_func = nullptr,
            .endvalue = uint32_t(-1), .cond = uint32_t(-1)
        };
        host_bseq.pulses[3] = Seq::HostSeq::Pulse{
            .id = 3, .time = 8, .measure = uint32_t(-1), .len = uint32_t(-1),
            .value = 3, .chn = 1, .ramp_func = nullptr, .endvalue = 3, .cond = uint32_t(-1)
        };
        host_bseq.assignments.resize(0);
        host_bseq.assumptions.resize(0);
        host_bseq.branches.resize(0);
        host_bseq.default_branch = -1;
        host_bseq.endtimes.resize(1);
        host_bseq.endtimes[0] = 8;

        host_bseq.ndirect_assumes = 0;
        host_bseq.assumptions_idx.resize(1); // nneed_order
        host_bseq.assumptions_idx[0] = uint32_t(-1);

        host_bseq.cond_global_refs.resize(0);
    }

    host_seq.init();
    auto &ages = host_seq.global_ages();
    REQUIRE(ages.size() == 3);
    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 0);
    REQUIRE(ages[2] == 0);
    host_seq.set_global(0, 0.4); // this will make G0 == 0
    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 0);
    REQUIRE(ages[2] == 0);

    // Run 1
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
    host_seq.init_run();
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(0));

    host_seq.pre_run();
    auto host_bseq = &host_seq.seqs[0];
    REQUIRE(host_bseq->length == 26);
    REQUIRE(host_seq.start_values[0].f64 == Approx(1.81));

    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 1);
    REQUIRE(ages[2] == 1);

    REQUIRE(host_seq.values[0].i64 == 20);
    REQUIRE(host_seq.values[1].f64 == 1.2);
    REQUIRE(host_seq.values[2].f64 == 1.81);
    REQUIRE(host_seq.values[3].f64 == 0.4);
    REQUIRE(host_seq.values[4].i64 == 0);
    REQUIRE(host_seq.values[5].i64 == 20);
    REQUIRE(host_seq.values[6].f64 == 0.2);
    REQUIRE(host_seq.values[7].f64 == 1.2);
    REQUIRE(host_seq.values[8].i64 == 26);

    REQUIRE(host_bseq->pulses.size() == 4);
    REQUIRE(host_bseq->pulses[0].id == 2);
    REQUIRE(host_bseq->pulses[0].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[0].endvalue].f64 == 1.81);
    REQUIRE(host_bseq->pulses[1].id == 1);
    REQUIRE(host_bseq->pulses[1].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[1].endvalue].f64 == 1.2);
    REQUIRE(host_bseq->pulses[2].id == 4); // measure
    REQUIRE(host_bseq->pulses[2].chn == 1);
    REQUIRE(host_bseq->pulses[2].len == uint32_t(-2));
    REQUIRE(host_bseq->pulses[3].id == 3);
    REQUIRE(host_bseq->pulses[3].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[3].endvalue].f64 == 0.4);

    REQUIRE(host_seq.post_run() == 0);
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));

    // Run 2
    host_seq.init_run();
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(0));
    host_seq.set_global(0, 25.6);

    REQUIRE(ages[0] == 2);
    REQUIRE(ages[1] == 1);
    REQUIRE(ages[2] == 1);

    auto err = expect_error<std::runtime_error>([&] {
        host_seq.pre_run();
    });
    REQUIRE(err.what() == "Going back in time not allowed."s);

    REQUIRE(host_seq.values[8].i64 == 46);
}

TEST_CASE("backward_time_assumption") {
    // In order to trigger a go back in time error,
    // we need a pulse time that depends on a measure.

    // Const:
    //   C0[V0]: time 1 (20)
    //   C1[V1]: val 1 (1.83)
    //   C2[V2]: val 2 (2.3)
    // Slot:
    //   S0[V3]: for time 2
    // Global Value:
    //   G0[V4]: round S0 (time 2)
    //   G1[V5]: C0 + G0
    // Channel:
    //   CV1[V6]: Channel 1 value
    // BS1:
    //   Values:
    //     MV0[V7]: measure value of M1
    //     N0[V8]: round((1.82 - MV0) * 10)
    //     N1[V9]: G1 + N0
    //   Pulse:
    //     P1: @C0 val=C1
    //     P2: @G0 val=C2
    //     M1: @G1 measure=MV0
    //     P3: @N1 val=S0
    //   Assumption:
    //     N0 >= 0
    //   Branch:
    //     Default: end
    Seq::HostSeq host_seq;
    host_seq.nconsts = 3;
    host_seq.nglobals = 1;
    host_seq.nglobal_vals = 2;
    host_seq.nchannels = 1;
    host_seq.nshared = (host_seq.nconsts + host_seq.nglobals +
                        host_seq.nglobal_vals + host_seq.nchannels);
    host_seq.npublic_globals = 1;

    host_seq.values.resize(10);
    host_seq.values[0].i64 = 20;
    host_seq.values[1].f64 = 1.83;
    host_seq.values[2].f64 = 2.3;
    host_seq.default_values.resize(1);
    host_seq.default_values[0].f64 = 0.2;
    host_seq.types.resize(6); // nconst + nglobals + nglobal_vals
    // Const
    host_seq.types[0] = Seq::HostSeq::Type::Int64;
    host_seq.types[1] = Seq::HostSeq::Type::Float64;
    host_seq.types[2] = Seq::HostSeq::Type::Float64;
    // Slot
    host_seq.types[3] = Seq::HostSeq::Type::Float64;
    // Global Value
    host_seq.types[4] = Seq::HostSeq::Type::Int64;
    host_seq.types[5] = Seq::HostSeq::Type::Int64;

    host_seq.depends.resize(2); // nglobal_vals
    host_seq.depends[0].push_back(0);
    host_seq.depends[1].push_back(0);
    host_seq.global_evals.resize(2); // nglobal_vals
    host_seq.global_evals[0] = [] (void *_p) {
        auto p = (Seq::HostSeq::Value*)_p;
        p[4].i64 = round<int64_t>(p[3].f64);
    };
    host_seq.global_evals[1] = [] (void *_p) {
        auto p = (Seq::HostSeq::Value*)_p;
        p[5].i64 = p[4].i64 + 20;
    };

    host_seq.seqs.resize(1);

    {
        auto &host_bseq = host_seq.seqs[0];
        host_bseq.id = 1;

        host_bseq.nmeasure = 1;
        host_bseq.ndirect = 0;
        host_bseq.nneed_order = 2;

        host_bseq.types.resize(2); // ndirect + nneed_order
        host_bseq.types[0] = Seq::HostSeq::Type::Int64;
        host_bseq.types[1] = Seq::HostSeq::Type::Int64;
        host_bseq.evals.resize(2); // ndirect + nneed_order
        host_bseq.evals[0] = [] (void *_p) {
            auto p = (Seq::HostSeq::Value*)_p;
            p[8].i64 = round<int64_t>(10 * (1.82 - p[7].f64));
        };
        host_bseq.evals[1] = [] (void *_p) {
            auto p = (Seq::HostSeq::Value*)_p;
            p[9].i64 = p[5].i64 + p[8].i64;
        };

        host_bseq.global_refs.resize(2);
        host_bseq.global_refs[0] = 0;
        host_bseq.global_refs[1] = 1;
        host_bseq.reverse_depends.resize(1); // nmeasure
        host_bseq.reverse_depends[0].push_back(0);
        host_bseq.reverse_depends[0].push_back(1);
        host_bseq.deps_count.resize(2); // nneed_order
        host_bseq.deps_count[0] = 1;
        host_bseq.deps_count[1] = 1;

        host_bseq.pulses.resize(4);
        host_bseq.pulses[0] = Seq::HostSeq::Pulse{
            .id = 1, .time = 0, .measure = uint32_t(-1), .len = uint32_t(-1),
            .value = 1, .chn = 1, .ramp_func = nullptr, .endvalue = 1, .cond = uint32_t(-1)
        };
        host_bseq.pulses[1] = Seq::HostSeq::Pulse{
            .id = 2, .time = 4, .measure = uint32_t(-1), .len = uint32_t(-1),
            .value = 2, .chn = 1, .ramp_func = nullptr, .endvalue = 2, .cond = uint32_t(-1)
        };
        host_bseq.pulses[2] = Seq::HostSeq::Pulse{
            .id = 3, .time = 5, .measure = 0, .len = uint32_t(-2),
            .value = uint32_t(-1), .chn = 1, .ramp_func = nullptr,
            .endvalue = uint32_t(-1), .cond = uint32_t(-1)
        };
        host_bseq.pulses[3] = Seq::HostSeq::Pulse{
            .id = 7, .time = 9, .measure = uint32_t(-1), .len = uint32_t(-1),
            .value = 3, .chn = 1, .ramp_func = nullptr, .endvalue = 3, .cond = uint32_t(-1)
        };
        host_bseq.assignments.resize(0);
        host_bseq.assumptions.resize(1);
        host_bseq.assumptions[0] = Seq::HostSeq::Assumption{
            .sign = Seq::Sign::NonNeg, .value = 8, .id = 13
        };
        host_bseq.branches.resize(0);
        host_bseq.default_branch = -1;
        host_bseq.endtimes.resize(1);
        host_bseq.endtimes[0] = 9;

        host_bseq.ndirect_assumes = 0;
        host_bseq.assumptions_idx.resize(2); // nneed_order
        host_bseq.assumptions_idx[0] = 0;
        host_bseq.assumptions_idx[1] = uint32_t(-1);

        host_bseq.cond_global_refs.resize(0);
    }

    host_seq.init();
    auto &ages = host_seq.global_ages();
    REQUIRE(ages.size() == 3);
    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 0);
    REQUIRE(ages[2] == 0);
    host_seq.set_global(0, 0.4); // this will make G0 == 0
    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 0);
    REQUIRE(ages[2] == 0);

    // Run 1
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
    host_seq.init_run();
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(0));

    host_seq.pre_run();
    auto host_bseq = &host_seq.seqs[0];
    REQUIRE(host_bseq->length == 20);
    REQUIRE(host_seq.start_values[0].f64 == Approx(2.3));

    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 1);
    REQUIRE(ages[2] == 1);

    REQUIRE(host_seq.values[0].i64 == 20);
    REQUIRE(host_seq.values[1].f64 == 1.83);
    REQUIRE(host_seq.values[2].f64 == 2.3);
    REQUIRE(host_seq.values[3].f64 == 0.4);
    REQUIRE(host_seq.values[4].i64 == 0);
    REQUIRE(host_seq.values[5].i64 == 20);
    REQUIRE(host_seq.values[6].f64 == 0.2);
    REQUIRE(host_seq.values[7].f64 == 1.83);
    REQUIRE(host_seq.values[8].i64 == 0);
    REQUIRE(host_seq.values[9].i64 == 20);

    REQUIRE(host_bseq->pulses.size() == 4);
    REQUIRE(host_bseq->pulses[0].id == 2);
    REQUIRE(host_bseq->pulses[0].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[0].endvalue].f64 == 2.3);
    REQUIRE(host_bseq->pulses[1].id == 1);
    REQUIRE(host_bseq->pulses[1].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[1].endvalue].f64 == 1.83);
    REQUIRE(host_bseq->pulses[2].id == 3); // measure
    REQUIRE(host_bseq->pulses[2].chn == 1);
    REQUIRE(host_bseq->pulses[2].len == uint32_t(-2));
    REQUIRE(host_bseq->pulses[3].id == 7);
    REQUIRE(host_bseq->pulses[3].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[3].endvalue].f64 == 0.4);

    REQUIRE(host_seq.post_run() == 0);
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));

    // Run 2
    host_seq.init_run();
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(0));
    host_seq.set_global(0, 25.6);

    REQUIRE(ages[0] == 2);
    REQUIRE(ages[1] == 1);
    REQUIRE(ages[2] == 1);

    auto err = expect_error<Seq::Error>([&] {
        host_seq.pre_run();
    });
    REQUIRE(err.type == Seq::Error::Type::EventTime);
    REQUIRE(err.code == uint16_t(Seq::Error::EventTime::NegTime));
    REQUIRE(err.type1 == Seq::Error::Type::EventTime);
    REQUIRE(err.id1 == 13);
    REQUIRE(err.what() == "Non-negative time expected."s);

    REQUIRE(host_seq.values[8].i64 == -5);
    // This value should not have been evaluated
    REQUIRE(host_seq.values[9].i64 == 20);
}

TEST_CASE("backward_time_id_assumption") {
    // In order to trigger a go back in time error,
    // we need a pulse time that depends on a measure.

    // Const:
    //   C0[V0]: time 1 (20)
    //   C1[V1]: val 1 (1.2)
    //   C2[V2]: val 2 (1.81)
    // Slot:
    //   S0[V3]: for time 2
    // Global Value:
    //   G0[V4]: round S0 (time 2)
    //   G1[V5]: C0 + G0
    // Channel:
    //   CV1[V6]: Channel 1 value
    // BS1:
    //   Values:
    //     MV0[V7]: measure value of M1
    //     N0[V8]: round((1.82 - MV0) * 10)
    //     N1[V9]: G1 + N0
    //   Pulse:
    //     P1: @C0 val=C1
    //     P2: @G0 val=C2
    //     M1: @G1 measure=MV0
    //     P3: @N0 val=S0
    //   Assumption:
    //     N0 > 0
    //   Branch:
    //     Default: end
    Seq::HostSeq host_seq;
    host_seq.nconsts = 3;
    host_seq.nglobals = 1;
    host_seq.nglobal_vals = 2;
    host_seq.nchannels = 1;
    host_seq.nshared = (host_seq.nconsts + host_seq.nglobals +
                        host_seq.nglobal_vals + host_seq.nchannels);
    host_seq.npublic_globals = 1;

    host_seq.values.resize(10);
    host_seq.values[0].i64 = 20;
    host_seq.values[1].f64 = 1.2;
    host_seq.values[2].f64 = 1.81;
    host_seq.default_values.resize(1);
    host_seq.default_values[0].f64 = 0.2;
    host_seq.types.resize(6); // nconst + nglobals + nglobal_vals
    // Const
    host_seq.types[0] = Seq::HostSeq::Type::Int64;
    host_seq.types[1] = Seq::HostSeq::Type::Float64;
    host_seq.types[2] = Seq::HostSeq::Type::Float64;
    // Slot
    host_seq.types[3] = Seq::HostSeq::Type::Float64;
    // Global Value
    host_seq.types[4] = Seq::HostSeq::Type::Int64;
    host_seq.types[5] = Seq::HostSeq::Type::Int64;

    host_seq.depends.resize(2); // nglobal_vals
    host_seq.depends[0].push_back(0);
    host_seq.depends[1].push_back(0);
    host_seq.global_evals.resize(2); // nglobal_vals
    host_seq.global_evals[0] = [] (void *_p) {
        auto p = (Seq::HostSeq::Value*)_p;
        p[4].i64 = round<int64_t>(p[3].f64);
    };
    host_seq.global_evals[1] = [] (void *_p) {
        auto p = (Seq::HostSeq::Value*)_p;
        p[5].i64 = p[4].i64 + 20;
    };

    host_seq.seqs.resize(1);

    {
        auto &host_bseq = host_seq.seqs[0];
        host_bseq.id = 1;

        host_bseq.nmeasure = 1;
        host_bseq.ndirect = 0;
        host_bseq.nneed_order = 2;

        host_bseq.types.resize(2); // ndirect + nneed_order
        host_bseq.types[0] = Seq::HostSeq::Type::Int64;
        host_bseq.types[1] = Seq::HostSeq::Type::Int64;
        host_bseq.evals.resize(2); // ndirect + nneed_order
        host_bseq.evals[0] = [] (void *_p) {
            auto p = (Seq::HostSeq::Value*)_p;
            p[8].i64 = round<int64_t>(10 * (1.82 - p[7].f64));
        };
        host_bseq.evals[1] = [] (void *_p) {
            auto p = (Seq::HostSeq::Value*)_p;
            p[9].i64 = p[5].i64 + p[8].i64;
        };

        host_bseq.global_refs.resize(2);
        host_bseq.global_refs[0] = 0;
        host_bseq.global_refs[1] = 1;
        host_bseq.reverse_depends.resize(1); // nmeasure
        host_bseq.reverse_depends[0].push_back(0);
        host_bseq.reverse_depends[0].push_back(1);
        host_bseq.deps_count.resize(2); // nneed_order
        host_bseq.deps_count[0] = 1;
        host_bseq.deps_count[1] = 1;

        host_bseq.pulses.resize(4);
        host_bseq.pulses[0] = Seq::HostSeq::Pulse{
            .id = 1, .time = 0, .measure = uint32_t(-1), .len = uint32_t(-1),
            .value = 1, .chn = 1, .ramp_func = nullptr, .endvalue = 1, .cond = uint32_t(-1)
        };
        host_bseq.pulses[1] = Seq::HostSeq::Pulse{
            .id = 2, .time = 4, .measure = uint32_t(-1), .len = uint32_t(-1),
            .value = 2, .chn = 1, .ramp_func = nullptr, .endvalue = 2, .cond = uint32_t(-1)
        };
        host_bseq.pulses[2] = Seq::HostSeq::Pulse{
            .id = 4, .time = 5, .measure = 0, .len = uint32_t(-2),
            .value = uint32_t(-1), .chn = 1, .ramp_func = nullptr,
            .endvalue = uint32_t(-1), .cond = uint32_t(-1)
        };
        host_bseq.pulses[3] = Seq::HostSeq::Pulse{
            .id = 3, .time = 9, .measure = uint32_t(-1), .len = uint32_t(-1),
            .value = 3, .chn = 1, .ramp_func = nullptr, .endvalue = 3, .cond = uint32_t(-1)
        };
        host_bseq.assignments.resize(0);
        host_bseq.assumptions.resize(1);
        host_bseq.assumptions[0] = Seq::HostSeq::Assumption{
            .sign = Seq::Sign::Pos, .value = 8, .id = 381
        };
        host_bseq.branches.resize(0);
        host_bseq.default_branch = -1;
        host_bseq.endtimes.resize(1);
        host_bseq.endtimes[0] = 9;

        host_bseq.ndirect_assumes = 0;
        host_bseq.assumptions_idx.resize(2); // nneed_order
        host_bseq.assumptions_idx[0] = 0;
        host_bseq.assumptions_idx[1] = uint32_t(-1);

        host_bseq.cond_global_refs.resize(0);
    }

    host_seq.init();
    auto &ages = host_seq.global_ages();
    REQUIRE(ages.size() == 3);
    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 0);
    REQUIRE(ages[2] == 0);
    host_seq.set_global(0, 0.4); // this will make G0 == 0
    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 0);
    REQUIRE(ages[2] == 0);

    // Run 1
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));
    host_seq.init_run();
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(0));

    host_seq.pre_run();
    auto host_bseq = &host_seq.seqs[0];
    REQUIRE(host_bseq->length == 26);
    REQUIRE(host_seq.start_values[0].f64 == Approx(1.81));

    REQUIRE(ages[0] == 1);
    REQUIRE(ages[1] == 1);
    REQUIRE(ages[2] == 1);

    REQUIRE(host_seq.values[0].i64 == 20);
    REQUIRE(host_seq.values[1].f64 == 1.2);
    REQUIRE(host_seq.values[2].f64 == 1.81);
    REQUIRE(host_seq.values[3].f64 == 0.4);
    REQUIRE(host_seq.values[4].i64 == 0);
    REQUIRE(host_seq.values[5].i64 == 20);
    REQUIRE(host_seq.values[6].f64 == 0.2);
    REQUIRE(host_seq.values[7].f64 == 1.2);
    REQUIRE(host_seq.values[8].i64 == 6);
    REQUIRE(host_seq.values[9].i64 == 26);

    REQUIRE(host_bseq->pulses.size() == 4);
    REQUIRE(host_bseq->pulses[0].id == 2);
    REQUIRE(host_bseq->pulses[0].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[0].endvalue].f64 == 1.81);
    REQUIRE(host_bseq->pulses[1].id == 1);
    REQUIRE(host_bseq->pulses[1].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[1].endvalue].f64 == 1.2);
    REQUIRE(host_bseq->pulses[2].id == 4); // measure
    REQUIRE(host_bseq->pulses[2].chn == 1);
    REQUIRE(host_bseq->pulses[2].len == uint32_t(-2));
    REQUIRE(host_bseq->pulses[3].id == 3);
    REQUIRE(host_bseq->pulses[3].chn == 1);
    REQUIRE(host_seq.values[host_bseq->pulses[3].endvalue].f64 == 0.4);

    REQUIRE(host_seq.post_run() == 0);
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(-1));

    // Run 2
    host_seq.init_run();
    REQUIRE(host_seq.cur_seq_idx() == uint32_t(0));
    host_seq.set_global(0, 25.6);

    REQUIRE(ages[0] == 2);
    REQUIRE(ages[1] == 1);
    REQUIRE(ages[2] == 1);

    auto err = expect_error<Seq::Error>([&] {
        host_seq.pre_run();
    });
    REQUIRE(err.type == Seq::Error::Type::EventTime);
    REQUIRE(err.code == uint16_t(Seq::Error::EventTime::NonPosTime));
    REQUIRE(err.type1 == Seq::Error::Type::EventTime);
    REQUIRE(err.id1 == 381);
    REQUIRE(err.what() == "Positive time expected."s);

    REQUIRE(host_seq.values[8].i64 == 0);
    // This value should not have been evaluated
    REQUIRE(host_seq.values[9].i64 == 26);
}
