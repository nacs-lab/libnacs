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

#ifndef _COMPILER_HELPER
#define _COMPILER_HELPER

#include "../../lib/nacs-seq/seq.h"
#include "../../lib/nacs-seq/error.h"
#include "../../lib/nacs-seq/compiler.h"
#include "../../lib/nacs-seq/host_seq.h"

#include "../../lib/nacs-utils/llvm/codegen.h"
#include "../../lib/nacs-utils/llvm/execute.h"

#include <algorithm>

#include <catch2/catch.hpp>

namespace {

using namespace NaCs;

struct CompileSeq {
    CompileSeq(Seq::Seq &seq)
        : seq(seq),
          compiler(host_seq, seq, engine,
                   seq.env().cg_context()->get_extern_resolver())
    {
        obj_id = compiler.compile();
        REQUIRE(obj_id);
        test_host_seq();
    }
    void test_host_seq()
    {
        REQUIRE(host_seq.nshared == host_seq.nconsts + host_seq.nglobals +
                host_seq.nglobal_vals + host_seq.nchannels);
        REQUIRE(host_seq.npublic_globals <= host_seq.nglobals);
        REQUIRE(host_seq.values.size() >= host_seq.nshared);
        REQUIRE(host_seq.default_values.size() == host_seq.nchannels);
        REQUIRE(host_seq.types.size() == host_seq.nconsts + host_seq.nglobals +
                host_seq.nglobal_vals);
        REQUIRE(host_seq.depends.size() == host_seq.nglobal_vals);
        for (auto &deps: host_seq.depends) {
            REQUIRE(!deps.empty());
            for (auto dep: deps)
                REQUIRE(dep < host_seq.nglobals);
            REQUIRE(std::is_sorted(deps.begin(), deps.end()));
        }
        REQUIRE(host_seq.global_evals.size() == host_seq.nglobal_vals);
        for (auto p: host_seq.global_evals)
            REQUIRE(p);
        REQUIRE(host_seq.seqs.size() == seq.get_basicseqs().size());
        for (auto &host_bseq: host_seq.seqs) {
            uint32_t bseq_nvars = (host_seq.nshared + host_bseq.nmeasure +
                                   host_bseq.ndirect + host_bseq.nneed_order);
            REQUIRE(host_seq.values.size() >= bseq_nvars);
            REQUIRE(host_bseq.types.size() == host_bseq.ndirect + host_bseq.nneed_order);
            REQUIRE(host_bseq.evals.size() == host_bseq.ndirect + host_bseq.nneed_order);
            for (auto p: host_bseq.evals)
                REQUIRE(p);
            REQUIRE(std::is_sorted(host_bseq.global_refs.begin(), host_bseq.global_refs.end()));
            for (auto ref: host_bseq.global_refs)
                REQUIRE(ref < host_seq.nglobal_vals);
            REQUIRE(std::is_sorted(host_bseq.cond_global_refs.begin(),
                                   host_bseq.cond_global_refs.end()));
            for (auto ref: host_bseq.cond_global_refs)
                REQUIRE(ref < host_seq.nglobal_vals);
            std::vector<uint32_t> deps_count(host_bseq.nneed_order, 0);
            REQUIRE(host_bseq.reverse_depends.size() == host_bseq.nmeasure);
            for (auto &rdeps: host_bseq.reverse_depends) {
                REQUIRE(std::is_sorted(rdeps.begin(), rdeps.end()));
                for (auto rdep: rdeps) {
                    REQUIRE(rdep < host_bseq.nneed_order);
                    deps_count[rdep] += 1;
                }
            }
            REQUIRE(host_bseq.deps_count.size() == host_bseq.nneed_order);
            REQUIRE(host_bseq.deps_count == deps_count);
            for (auto &assign: host_bseq.assignments) {
                REQUIRE(assign.global_id < host_seq.nglobals);
                REQUIRE(assign.value < bseq_nvars);
            }
            for (auto &assume: host_bseq.assumptions) {
                REQUIRE(assume.sign != Seq::Sign::Unknown);
                REQUIRE(assume.value < bseq_nvars);
                if (assume.value >= host_seq.nshared + host_bseq.nmeasure + host_bseq.ndirect) {
                    auto need_order_id = assume.value -
                        (host_seq.nshared + host_bseq.nmeasure + host_bseq.ndirect);
                    auto assume_id = host_bseq.assumptions_idx[need_order_id];
                    REQUIRE(assume_id != uint32_t(-1));
                    REQUIRE(&host_bseq.assumptions[assume_id] == &assume);
                }
            }
            for (auto &br: host_bseq.branches) {
                if (br.target != uint32_t(-1))
                    REQUIRE(br.target < host_seq.seqs.size());
                REQUIRE(br.cond < bseq_nvars);
            }
            if (host_bseq.default_branch != uint32_t(-1))
                REQUIRE(host_bseq.default_branch < host_seq.seqs.size());
            for (auto et: host_bseq.endtimes)
                REQUIRE(et < bseq_nvars);
            REQUIRE(host_bseq.ndirect_assumes <= host_bseq.assumptions.size());
            REQUIRE(host_bseq.assumptions_idx.size() == host_bseq.nneed_order);
            uint32_t nneed_order_assume = 0;
            for (auto ai: host_bseq.assumptions_idx) {
                if (ai == uint32_t(-1))
                    continue;
                nneed_order_assume++;
                REQUIRE(ai < host_bseq.assumptions.size());
            }
            REQUIRE(host_bseq.ndirect_assumes + nneed_order_assume ==
                    host_bseq.assumptions.size());
        }
    }
    ~CompileSeq()
    {
        engine.free(obj_id);
    }

    Seq::Seq &seq;
    Seq::HostSeq host_seq;
    LLVM::Exe::Engine engine;
    Seq::Compiler compiler;
    uint64_t obj_id = 0;
};
}

#endif
