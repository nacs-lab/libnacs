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

#include "zynq_helper.h"
#include "../error_helper.h"

#include "../../lib/nacs-seq/builder.h"
#include "../../lib/nacs-seq/manager.h"
#include "../../lib/nacs-seq/zynq/backend.h"
#include "../../lib/nacs-utils/streams.h"

using namespace std::literals::string_literals;

// A dummy backend to enable and set the clock output on the Zynq backend.
// We only need this for the precomputed version since it must be done
// inside the `create_sequence` from `prepare` of another backend.
// For the runtime version, we don't have to do it in `prepare_run` and can simply
// do it before we call `pre_run` because we know the Zynq backend shouldn't care.
class ClockEnabler : public Device {
public:
    using Device::Device;

private:
    void add_channel(uint32_t chn_id, const std::string &chn_name) override {}
    // bool check_noramp(uint32_t chn_id, const std::string &chn_name) override;
    void prepare(Manager::ExpSeq &expseq, Compiler &compiler) override
    {
        auto zynq_dev = Zynq::Backend::cast(expseq.get_device("FPGA1", true));
        if (!zynq_dev)
            throw std::runtime_error("Cannot find clock generator device.");
        zynq_dev->enable_clock();
        for (auto &[idx, times]: m_active_times) {
            zynq_dev->set_clock_active_time(idx, m_half_period, times);
        }
    }
    void generate(Manager::ExpSeq &expseq, Compiler &compiler) override {}

    // void init_run(HostSeq &host_seq) override;
    // void prepare_run(HostSeq &host_seq) override;
    void pre_run(HostSeq &host_seq) override {}
    void start(HostSeq &host_seq) override {}
    void cancel(HostSeq &host_seq) override {}
    void wait(HostSeq &host_seq) override {}
    // void finish_run(HostSeq &host_seq) override;

    void config(const YAML::Node&) override {}
    void parse_data(const uint8_t *data, size_t len) override
    {
        Mem::Reader reader(data, len);
        m_half_period = reader.read<uint64_t>();
        auto ntimes = reader.read<uint32_t>();
        for (uint32_t i = 0; i < ntimes; i++) {
            auto idx = reader.read<uint32_t>();
            auto sz = reader.read<uint32_t>();
            auto p = (std::pair<int64_t,int64_t>*)(reader.data + reader.cursor);
            reader.forward_cursor(sizeof(std::pair<int64_t,int64_t>) * sz);
            m_active_times[idx].assign(p, p + sz);
        }
    }

    uint64_t m_half_period;
    std::map<uint32_t,std::vector<std::pair<int64_t,int64_t>>> m_active_times;
};

static const char *const config_str =
    "tick_per_sec: 1000000000000\n"
    "devices:\n"
    "  Clocker:\n"
    "    backend: clocker\n"
    "  FPGA1:\n"
    "    backend: zynq\n"
    "    config:\n"
    "      start_ttl_chn: 0\n";

static Manager &manager()
{
    static Manager mgr;
    static Device::Register<ClockEnabler> register_backend("clocker");
    static Call init([&] {
        // TODO disable dummy mode if we have a server URL in the ENV.
        mgr.dummy_mode = true;
        mgr.load_config_string(config_str);
    });
    return mgr;
}

struct ExpSeqDeleter {
    ExpSeqDeleter(Manager &mgr)
        : m_mgr(&mgr)
    {
    }
    void operator()(Manager::ExpSeq *expseq) const
    {
        m_mgr->free_sequence(expseq);
    }
private:
    Manager *m_mgr;
};

struct ostream : uvector_ostream {
    using uvector_ostream::write;
    template<typename T>
    void write(T v)
    {
        STATIC_REQUIRE(std::is_trivial_v<T>);
        write((const char*)&v, sizeof(T));
    }
    void write_string(const char *s)
    {
        write(s, strlen(s) + 1);
    }
    std::unique_ptr<Manager::ExpSeq,ExpSeqDeleter> create_sequence()
    {
        auto data = get_buf();
        auto &mgr = manager();
        return {mgr.create_sequence(data.data(), data.size()), ExpSeqDeleter(mgr)};
    }
};

TEST_CASE("empty") {
    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(0);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(2);
    stm.write_string("FPGA1/DDS2/FREQ");
    stm.write_string("FPGA1/TTL3");
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(0);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(0);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(0);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(2);
    {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(0);
        // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
        stm.write<uint32_t>(0);
        // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
        stm.write<uint32_t>(0);
        // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
        stm.write<uint32_t>(0);
        // [default_target: 4B]
        stm.write<uint32_t>(2);
    }
    {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(0);
        // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
        stm.write<uint32_t>(0);
        // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
        stm.write<uint32_t>(0);
        // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
        stm.write<uint32_t>(0);
        // [default_target: 4B]
        stm.write<uint32_t>(0);
    }
    // [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
    stm.write<uint32_t>(0);
    // [nbackenddatas: 4B][[device_name: NUL-terminated string]
    //                     [size: 4B][data: size B] x nbackenddatas]
    stm.write<uint32_t>(0);

    auto expseq = stm.create_sequence();
    REQUIRE(expseq);

    REQUIRE(expseq->devices().size() == 1);
    auto zynq_dev = expseq->get_device("FPGA1", false);
    REQUIRE(zynq_dev);
    auto zynq_backend = Zynq::Backend::cast(zynq_dev);
    REQUIRE(zynq_backend);
    REQUIRE(!zynq_backend->has_generator(0));
    REQUIRE(zynq_backend->pregenerated(0, true));
    REQUIRE(!zynq_backend->pregenerated(0, false));
    REQUIRE(!zynq_backend->has_generator(1));
    REQUIRE(!zynq_backend->pregenerated(1, true));
    REQUIRE(zynq_backend->pregenerated(1, false));

    {
        auto bc = zynq_backend->get_bytecode(0, true);
        auto bc_ver = zynq_backend->get_bytecode_version(0, true);
        Checker checker(bc.data(), bc.size());
        INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
        REQUIRE(checker.cmp<int64_t>(0)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(9)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq2{7, 2, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 97}));
        REQUIRE(checker.is_end());
    }

    {
        auto bc = zynq_backend->get_bytecode(1, false);
        auto bc_ver = zynq_backend->get_bytecode_version(1, false);
        Checker checker(bc.data(), bc.size());
        INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
        REQUIRE(checker.cmp<int64_t>(0)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(9)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 97}));
        REQUIRE(checker.is_end());
    }

    expseq->init_run();
    expseq->pre_run();
    expseq->start();
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 2);
    expseq->pre_run();
    expseq->start();
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 0);
}

TEST_CASE("empty+clock_unset") {
    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(0);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(3);
    stm.write_string("FPGA1/DDS2/FREQ");
    stm.write_string("FPGA1/TTL3");
    stm.write_string("Clocker/Dummy");
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(0);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(0);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(0);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(2);
    {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(0);
        // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
        stm.write<uint32_t>(0);
        // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
        stm.write<uint32_t>(0);
        // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
        stm.write<uint32_t>(0);
        // [default_target: 4B]
        stm.write<uint32_t>(2);
    }
    {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(0);
        // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
        stm.write<uint32_t>(0);
        // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
        stm.write<uint32_t>(0);
        // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
        stm.write<uint32_t>(0);
        // [default_target: 4B]
        stm.write<uint32_t>(0);
    }
    // [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
    stm.write<uint32_t>(0);
    // [nbackenddatas: 4B][[device_name: NUL-terminated string]
    //                     [size: 4B][data: size B] x nbackenddatas]
    stm.write<uint32_t>(0);

    auto expseq = stm.create_sequence();
    REQUIRE(expseq);

    REQUIRE(expseq->devices().size() == 2);
    auto zynq_dev = expseq->get_device("FPGA1", false);
    REQUIRE(zynq_dev);
    auto zynq_backend = Zynq::Backend::cast(zynq_dev);
    REQUIRE(zynq_backend);
    REQUIRE(zynq_backend->has_generator(0));
    REQUIRE(!zynq_backend->pregenerated(0, true));
    REQUIRE(!zynq_backend->pregenerated(0, false));
    REQUIRE(zynq_backend->has_generator(1));
    REQUIRE(!zynq_backend->pregenerated(1, true));
    REQUIRE(!zynq_backend->pregenerated(1, false));

    expseq->init_run();
    zynq_backend->set_clock_active_time(0, 1000000, {});
    expseq->pre_run();
    {
        auto bc = zynq_backend->get_bytecode(0, true);
        auto bc_ver = zynq_backend->get_bytecode_version(0, true);
        Checker checker(bc.data(), bc.size());
        INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
        REQUIRE(checker.cmp<int64_t>(0)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(9)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq2{7, 2, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 296}));
        REQUIRE(checker.is_end());
    }
    expseq->start();
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 2);

    zynq_backend->set_clock_active_time(1, 1000000, {{0, 20}, {30, 100}});
    expseq->pre_run();
    {
        auto bc = zynq_backend->get_bytecode(1, false);
        auto bc_ver = zynq_backend->get_bytecode_version(1, false);
        Checker checker(bc.data(), bc.size());
        INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
        REQUIRE(checker.cmp<int64_t>(0)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(9)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 0, 0, 0}));
        // Sequence
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 99}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 4194}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 255}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 1795}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 99}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 14195}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 255}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 100}));
        REQUIRE(checker.is_end());
    }
    expseq->start();
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 0);

    expseq->init_run();
    zynq_backend->set_clock_active_time(0, 1000000, {});
    expseq->pre_run();
    {
        auto bc = zynq_backend->get_bytecode(0, true);
        auto bc_ver = zynq_backend->get_bytecode_version(0, true);
        Checker checker(bc.data(), bc.size());
        INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
        REQUIRE(checker.cmp<int64_t>(0)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(9)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq2{7, 2, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 296}));
        REQUIRE(checker.is_end());
    }
    expseq->start();
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 2);

    zynq_backend->set_clock_active_time(1, 1000000, {{0, 20}, {30, 100}});
    expseq->pre_run();
    {
        auto bc = zynq_backend->get_bytecode(1, false);
        auto bc_ver = zynq_backend->get_bytecode_version(1, false);
        Checker checker(bc.data(), bc.size());
        INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
        REQUIRE(checker.cmp<int64_t>(0)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(9)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 0, 0, 0}));
        // Sequence
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 99}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 4194}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 255}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 1795}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 99}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 14195}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 255}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 100}));
        REQUIRE(checker.is_end());
    }
    expseq->start();
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 0);
}

TEST_CASE("empty+clock_set") {
    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(0);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(3);
    stm.write_string("FPGA1/DDS2/FREQ");
    stm.write_string("FPGA1/TTL3");
    stm.write_string("Clocker/Dummy");
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(0);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(0);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(0);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(2);
    {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(0);
        // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
        stm.write<uint32_t>(0);
        // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
        stm.write<uint32_t>(0);
        // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
        stm.write<uint32_t>(0);
        // [default_target: 4B]
        stm.write<uint32_t>(2);
    }
    {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(0);
        // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
        stm.write<uint32_t>(0);
        // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
        stm.write<uint32_t>(0);
        // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
        stm.write<uint32_t>(0);
        // [default_target: 4B]
        stm.write<uint32_t>(0);
    }
    // [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
    stm.write<uint32_t>(0);
    // [nbackenddatas: 4B][[device_name: NUL-terminated string]
    //                     [size: 4B][data: size B] x nbackenddatas]
    stm.write<uint32_t>(1);
    stm.write_string("Clocker");
    stm.write<uint32_t>(60);
    // half_period
    stm.write<uint64_t>(1000000);
    // ntimes
    stm.write<uint32_t>(2);
    // bseq_idx
    stm.write<uint32_t>(0);
    // size
    stm.write<uint32_t>(0);
    // bseq_idx
    stm.write<uint32_t>(1);
    // size
    stm.write<uint32_t>(2);
    stm.write<uint64_t>(0);
    stm.write<uint64_t>(20);
    stm.write<uint64_t>(30);
    stm.write<uint64_t>(100);

    auto expseq = stm.create_sequence();
    REQUIRE(expseq);

    REQUIRE(expseq->devices().size() == 2);
    auto zynq_dev = expseq->get_device("FPGA1", false);
    REQUIRE(zynq_dev);
    auto zynq_backend = Zynq::Backend::cast(zynq_dev);
    REQUIRE(zynq_backend);
    REQUIRE(!zynq_backend->has_generator(0));
    REQUIRE(zynq_backend->pregenerated(0, true));
    REQUIRE(!zynq_backend->pregenerated(0, false));
    REQUIRE(!zynq_backend->has_generator(1));
    REQUIRE(!zynq_backend->pregenerated(1, true));
    REQUIRE(zynq_backend->pregenerated(1, false));
    {
        auto bc = zynq_backend->get_bytecode(0, true);
        auto bc_ver = zynq_backend->get_bytecode_version(0, true);
        Checker checker(bc.data(), bc.size());
        INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
        REQUIRE(checker.cmp<int64_t>(0)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(9)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq2{7, 2, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 296}));
        REQUIRE(checker.is_end());
    }
    {
        auto bc = zynq_backend->get_bytecode(1, false);
        auto bc_ver = zynq_backend->get_bytecode_version(1, false);
        Checker checker(bc.data(), bc.size());
        INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
        REQUIRE(checker.cmp<int64_t>(0)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(9)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 0, 0, 0}));
        // Sequence
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 99}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 4194}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 255}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 1795}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 99}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 14195}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 255}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 100}));
        REQUIRE(checker.is_end());
    }

    expseq->init_run();
    expseq->pre_run();
    expseq->start();
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 2);

    expseq->pre_run();
    expseq->start();
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 0);

    expseq->init_run();
    expseq->pre_run();
    expseq->start();
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 2);

    expseq->pre_run();
    expseq->start();
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 0);
}

TEST_CASE("loop") {
    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(0);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(2);
    stm.write_string("FPGA1/DDS2/FREQ");
    stm.write_string("FPGA1/TTL3");
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(0);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(0);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(0);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(2);
    {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(0);
        // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
        stm.write<uint32_t>(0);
        // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
        stm.write<uint32_t>(0);
        // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
        stm.write<uint32_t>(0);
        // [default_target: 4B]
        stm.write<uint32_t>(2);
    }
    {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(0);
        // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
        stm.write<uint32_t>(0);
        // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
        stm.write<uint32_t>(0);
        // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
        stm.write<uint32_t>(0);
        // [default_target: 4B]
        stm.write<uint32_t>(1);
    }
    // [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
    stm.write<uint32_t>(0);
    // [nbackenddatas: 4B][[device_name: NUL-terminated string]
    //                     [size: 4B][data: size B] x nbackenddatas]
    stm.write<uint32_t>(0);

    auto expseq = stm.create_sequence();
    REQUIRE(expseq);

    REQUIRE(expseq->devices().size() == 1);
    auto zynq_dev = expseq->get_device("FPGA1", false);
    REQUIRE(zynq_dev);
    auto zynq_backend = Zynq::Backend::cast(zynq_dev);
    REQUIRE(zynq_backend);
    REQUIRE(!zynq_backend->has_generator(0));
    REQUIRE(zynq_backend->pregenerated(0, true));
    REQUIRE(zynq_backend->pregenerated(0, false));
    REQUIRE(!zynq_backend->has_generator(1));
    REQUIRE(!zynq_backend->pregenerated(1, true));
    REQUIRE(zynq_backend->pregenerated(1, false));

    {
        auto bc = zynq_backend->get_bytecode(0, true);
        auto bc_ver = zynq_backend->get_bytecode_version(0, true);
        Checker checker(bc.data(), bc.size());
        INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
        REQUIRE(checker.cmp<int64_t>(0)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(9)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq2{7, 2, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 97}));
        REQUIRE(checker.is_end());
    }

    {
        auto bc = zynq_backend->get_bytecode(0, false);
        auto bc_ver = zynq_backend->get_bytecode_version(0, false);
        Checker checker(bc.data(), bc.size());
        INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
        REQUIRE(checker.cmp<int64_t>(0)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(9)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 97}));
        REQUIRE(checker.is_end());
    }

    {
        auto bc = zynq_backend->get_bytecode(1, false);
        auto bc_ver = zynq_backend->get_bytecode_version(1, false);
        Checker checker(bc.data(), bc.size());
        INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
        REQUIRE(checker.cmp<int64_t>(0)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(9)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 97}));
        REQUIRE(checker.is_end());
    }

    expseq->init_run();
    expseq->pre_run();
    expseq->start();
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 2);
    expseq->pre_run();
    expseq->start();
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 1);
    expseq->pre_run();
    expseq->start();
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 2);
    expseq->pre_run();
    expseq->start();
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 1);
}

TEST_CASE("setpulse") {
    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(6);
    // 1
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(0.301);
    // 2
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(0);
    // 3
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(120e6);
    // 4
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(12000000);
    // 5
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(5000000);
    // 6
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(15000000);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(3);
    stm.write_string("FPGA1/DDS4/FREQ");
    stm.write_string("FPGA1/DDS20/AMP");
    stm.write_string("FPGA1/TTL9");
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(3);
    stm.write<uint32_t>(1);
    stm.write(IR::Type::Float64);
    stm.write<double>(102e6);
    stm.write<uint32_t>(2);
    stm.write(IR::Type::Float64);
    stm.write<double>(0.4);
    stm.write<uint32_t>(3);
    stm.write(IR::Type::Bool);
    stm.write<bool>(true);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(0);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(0);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(1);
    {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(3);
        stm.write<uint8_t>((uint8_t)Sign::Pos);
        stm.write<uint32_t>(105);
        stm.write<uint32_t>(4);
        stm.write<uint32_t>(0);
        stm.write<uint8_t>((uint8_t)Sign::Pos);
        stm.write<uint32_t>(204);
        stm.write<uint32_t>(5);
        stm.write<uint32_t>(1);
        stm.write<uint8_t>((uint8_t)Sign::Unknown);
        stm.write<uint32_t>(209);
        stm.write<uint32_t>(6);
        stm.write<uint32_t>(2);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(3);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(3);
        // 1
        stm.write<uint32_t>(8); // id
        stm.write<uint32_t>(1); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(3); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(1); // chn
        // 2
        stm.write<uint32_t>(65); // id
        stm.write<uint32_t>(3); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(1); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(2); // chn
        // 3
        stm.write<uint32_t>(32); // id
        stm.write<uint32_t>(2); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(2); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(3); // chn
        // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
        stm.write<uint32_t>(0);
        // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
        stm.write<uint32_t>(0);
        // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
        stm.write<uint32_t>(0);
        // [default_target: 4B]
        stm.write<uint32_t>(0);
    }
    // [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
    stm.write<uint32_t>(0);
    // [nbackenddatas: 4B][[device_name: NUL-terminated string]
    //                     [size: 4B][data: size B] x nbackenddatas]
    stm.write<uint32_t>(0);

    auto expseq = stm.create_sequence();
    REQUIRE(expseq);

    REQUIRE(expseq->devices().size() == 1);
    auto zynq_dev = expseq->get_device("FPGA1", false);
    REQUIRE(zynq_dev);
    auto zynq_backend = Zynq::Backend::cast(zynq_dev);
    REQUIRE(zynq_backend);
    REQUIRE(!zynq_backend->has_generator(0));
    REQUIRE(zynq_backend->pregenerated(0, true));
    REQUIRE(!zynq_backend->pregenerated(0, false));

    {
        auto bc = zynq_backend->get_bytecode(0, true);
        auto bc_ver = zynq_backend->get_bytecode_version(0, true);
        Checker checker(bc.data(), bc.size());
        INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
        REQUIRE(checker.cmp<int64_t>(32000)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(0x201)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSFreq{6, 4, 0x775e802}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSAmp{10, 0, 20, 0x666}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 9}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        // Sequence
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 1196}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSFreq{6, 4, 0x8c6f2d6}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 450}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 9, 9}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 1496}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSAmp{10, 0, 20, 0x4d1}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 100}));
        REQUIRE(checker.is_end());
    }

    expseq->init_run();
    expseq->pre_run();
    expseq->start();
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 0);

    expseq->init_run();
    expseq->pre_run();
    expseq->start();
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 0);
}

TEST_CASE("ramp") {
    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(6);
    // 1
    stm.write(Builder::OpCode::Mul);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(0.8149072527885437);
    stm.write(Builder::ArgType::Arg);
    stm.write<uint32_t>(0);
    // 2
    stm.write(Builder::OpCode::Mul);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(2.9304029304029303e-8);
    stm.write(Builder::ArgType::Arg);
    stm.write<uint32_t>(0);
    // 3
    stm.write(Builder::OpCode::Sub);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(1.0);
    stm.write(Builder::ArgType::Node);
    stm.write<uint32_t>(2);
    // 4
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(12000000);
    // 5
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(15000000);
    // 6
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(10000000);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(2);
    stm.write_string("FPGA1/DDS4/FREQ");
    stm.write_string("FPGA1/DDS20/AMP");
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(0);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(0);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(0);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(1);
    {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(2);
        stm.write<uint8_t>((uint8_t)Sign::Pos);
        stm.write<uint32_t>(105);
        stm.write<uint32_t>(4);
        stm.write<uint32_t>(0);
        stm.write<uint8_t>((uint8_t)Sign::Pos);
        stm.write<uint32_t>(204);
        stm.write<uint32_t>(5);
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(2);
        // 1
        stm.write<uint32_t>(8); // id
        stm.write<uint32_t>(1); // time_id
        stm.write<uint32_t>(6); // len
        stm.write<uint32_t>(1); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(1); // chn
        // 2
        stm.write<uint32_t>(65); // id
        stm.write<uint32_t>(2); // time_id
        stm.write<uint32_t>(6); // len
        stm.write<uint32_t>(3); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(2); // chn
        // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
        stm.write<uint32_t>(0);
        // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
        stm.write<uint32_t>(0);
        // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
        stm.write<uint32_t>(0);
        // [default_target: 4B]
        stm.write<uint32_t>(0);
    }
    // [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
    stm.write<uint32_t>(0);
    // [nbackenddatas: 4B][[device_name: NUL-terminated string]
    //                     [size: 4B][data: size B] x nbackenddatas]
    stm.write<uint32_t>(0);

    auto expseq = stm.create_sequence();
    REQUIRE(expseq);

    REQUIRE(expseq->devices().size() == 1);
    auto zynq_dev = expseq->get_device("FPGA1", false);
    REQUIRE(zynq_dev);
    auto zynq_backend = Zynq::Backend::cast(zynq_dev);
    REQUIRE(zynq_backend);
    REQUIRE(!zynq_backend->has_generator(0));
    REQUIRE(zynq_backend->pregenerated(0, true));
    REQUIRE(!zynq_backend->pregenerated(0, false));

    {
        auto bc = zynq_backend->get_bytecode(0, true);
        auto bc_ver = zynq_backend->get_bytecode_version(0, true);
        Checker checker(bc.data(), bc.size());
        INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
        REQUIRE(checker.cmp<int64_t>(0)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(0x1)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq2{7, 4, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetAmp{11, 20, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        // Sequence
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 1246}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq4{9, 4, 500000}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq4{9, 4, 500000}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq4{9, 4, 500000}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq4{9, 4, 500000}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq4{9, 4, 500000}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSAmp{10, 0, 20, 0xfff}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetAmp{11, 20, 0x44})); // -60
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq4{9, 4, 1500000}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSAmp{10, 0, 20, 0xf4b}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq4{9, 4, 1000000}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSAmp{10, 0, 20, 0xed3}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq4{9, 4, 1000000}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSAmp{10, 0, 20, 0xe5b}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq4{9, 4, 1000000}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSAmp{10, 0, 20, 0xde3}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq4{9, 4, 1000000}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSAmp{10, 0, 20, 0xd6b}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq4{9, 4, 2000000}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSAmp{10, 0, 20, 0xcf3}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetAmp{11, 20, 0x44})); // -60
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetAmp{11, 20, 0x44})); // -60
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetAmp{11, 20, 0x44})); // -60
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetAmp{11, 20, 0x44})); // -60
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetAmp{11, 20, 0x44})); // -60
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSAmp{10, 0, 20, 0xb4f}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 100}));
        REQUIRE(checker.is_end());
    }

    expseq->init_run();
    expseq->pre_run();
    expseq->start();
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 0);

    expseq->init_run();
    expseq->pre_run();
    expseq->start();
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 0);
}

TEST_CASE("noramp") {
    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(3);
    // 1
    stm.write(Builder::OpCode::Mul);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(1e-5);
    stm.write(Builder::ArgType::Arg);
    stm.write<uint32_t>(0);
    // 2
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(12000000);
    // 3
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(10000000);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(1);
    stm.write_string("FPGA1/TTL2");
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(0);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(0);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(0);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(1);
    {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(1);
        stm.write<uint8_t>((uint8_t)Sign::Pos);
        stm.write<uint32_t>(105);
        stm.write<uint32_t>(2);
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(1);
        // 1
        stm.write<uint32_t>(23); // id
        stm.write<uint32_t>(1); // time_id
        stm.write<uint32_t>(3); // len
        stm.write<uint32_t>(1); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(1); // chn
        // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
        stm.write<uint32_t>(0);
        // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
        stm.write<uint32_t>(0);
        // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
        stm.write<uint32_t>(0);
        // [default_target: 4B]
        stm.write<uint32_t>(0);
    }
    // [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
    stm.write<uint32_t>(0);
    // [nbackenddatas: 4B][[device_name: NUL-terminated string]
    //                     [size: 4B][data: size B] x nbackenddatas]
    stm.write<uint32_t>(0);

    auto err = expect_error<Error>([&] {
        stm.create_sequence();
    });

    REQUIRE(err.type == Error::Type::Pulse);
    REQUIRE(err.code == uint16_t(Error::Pulse::NoRamp));
    REQUIRE(err.type1 == Error::Type::Pulse);
    REQUIRE(err.id1 == 23);
    REQUIRE(err.type2 == Error::Type::Channel);
    REQUIRE(err.id2 == 1);
    REQUIRE(err.what() == "Ramp not supported on channel"s);
}

TEST_CASE("setpulse_global") {
    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(6);
    // 1
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::Global);
    stm.write<uint32_t>(0);
    // 2
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::Global);
    stm.write<uint32_t>(1);
    // 3
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::Global);
    stm.write<uint32_t>(2);
    // 4
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(12000000);
    // 5
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(5000000);
    // 6
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(15000000);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(3);
    stm.write_string("FPGA1/DDS4/FREQ");
    stm.write_string("FPGA1/DDS20/AMP");
    stm.write_string("FPGA1/TTL9");
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(3);
    stm.write<uint32_t>(1);
    stm.write(IR::Type::Float64);
    stm.write<double>(102e6);
    stm.write<uint32_t>(2);
    stm.write(IR::Type::Float64);
    stm.write<double>(0.4);
    stm.write<uint32_t>(3);
    stm.write(IR::Type::Bool);
    stm.write<bool>(true);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(3);
    stm.write(IR::Type::Float64);
    stm.write(IR::Type::Float64);
    stm.write(IR::Type::Float64);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(0);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(1);
    {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(3);
        stm.write<uint8_t>((uint8_t)Sign::Pos);
        stm.write<uint32_t>(105);
        stm.write<uint32_t>(4);
        stm.write<uint32_t>(0);
        stm.write<uint8_t>((uint8_t)Sign::Pos);
        stm.write<uint32_t>(204);
        stm.write<uint32_t>(5);
        stm.write<uint32_t>(1);
        stm.write<uint8_t>((uint8_t)Sign::Unknown);
        stm.write<uint32_t>(209);
        stm.write<uint32_t>(6);
        stm.write<uint32_t>(2);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(3);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(3);
        // 1
        stm.write<uint32_t>(8); // id
        stm.write<uint32_t>(1); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(3); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(1); // chn
        // 2
        stm.write<uint32_t>(65); // id
        stm.write<uint32_t>(3); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(1); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(2); // chn
        // 3
        stm.write<uint32_t>(32); // id
        stm.write<uint32_t>(2); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(2); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(3); // chn
        // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
        stm.write<uint32_t>(0);
        // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
        stm.write<uint32_t>(0);
        // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
        stm.write<uint32_t>(0);
        // [default_target: 4B]
        stm.write<uint32_t>(0);
    }
    // [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
    stm.write<uint32_t>(0);
    // [nbackenddatas: 4B][[device_name: NUL-terminated string]
    //                     [size: 4B][data: size B] x nbackenddatas]
    stm.write<uint32_t>(0);

    auto expseq = stm.create_sequence();
    REQUIRE(expseq);

    REQUIRE(expseq->devices().size() == 1);
    auto zynq_dev = expseq->get_device("FPGA1", false);
    REQUIRE(zynq_dev);
    auto zynq_backend = Zynq::Backend::cast(zynq_dev);
    REQUIRE(zynq_backend);
    REQUIRE(zynq_backend->has_generator(0));
    REQUIRE(!zynq_backend->pregenerated(0, true));
    REQUIRE(!zynq_backend->pregenerated(0, false));

    expseq->init_run();
    expseq->pre_run();
    {
        auto bc = zynq_backend->get_bytecode(0, true);
        auto bc_ver = zynq_backend->get_bytecode_version(0, true);
        Checker checker(bc.data(), bc.size());
        INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
        REQUIRE(checker.cmp<int64_t>(32000)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(0x201)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSFreq{6, 4, 0x775e802}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSAmp{10, 0, 20, 0x666}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 9}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        // Sequence
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 1196}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSFreq{6, 4, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 450}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 9, 9}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 1496}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSAmp{10, 0, 20, 0}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 100}));
        REQUIRE(checker.is_end());
    }
    expseq->start();
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 0);

    expseq->init_run();
    expseq->host_seq.set_global(0, 0.301);
    expseq->host_seq.set_global(1, 1);
    expseq->host_seq.set_global(2, 120e6);
    expseq->pre_run();
    {
        auto bc = zynq_backend->get_bytecode(0, true);
        auto bc_ver = zynq_backend->get_bytecode_version(0, true);
        Checker checker(bc.data(), bc.size());
        INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
        REQUIRE(checker.cmp<int64_t>(32000)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(0x201)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSFreq{6, 4, 0x775e802}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSAmp{10, 0, 20, 0x666}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 9}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        // Sequence
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 1196}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSFreq{6, 4, 0x8c6f2d6}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 1950}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSAmp{10, 0, 20, 0x4d1}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 100}));
        REQUIRE(checker.is_end());
    }
    expseq->start();
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 0);
}

TEST_CASE("ramp_global") {
    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(6);
    // 1
    stm.write(Builder::OpCode::Mul);
    stm.write(Builder::ArgType::Global);
    stm.write<uint32_t>(0);
    stm.write(Builder::ArgType::Arg);
    stm.write<uint32_t>(0);
    // 2
    stm.write(Builder::OpCode::Mul);
    stm.write(Builder::ArgType::Global);
    stm.write<uint32_t>(1);
    stm.write(Builder::ArgType::Arg);
    stm.write<uint32_t>(0);
    // 3
    stm.write(Builder::OpCode::Sub);
    stm.write(Builder::ArgType::Global);
    stm.write<uint32_t>(2);
    stm.write(Builder::ArgType::Node);
    stm.write<uint32_t>(2);
    // 4
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(12000000);
    // 5
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(15000000);
    // 6
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(10000000);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(2);
    stm.write_string("FPGA1/DDS4/FREQ");
    stm.write_string("FPGA1/DDS20/AMP");
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(0);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(3);
    stm.write(IR::Type::Float64);
    stm.write(IR::Type::Float64);
    stm.write(IR::Type::Float64);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(0);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(1);
    {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(2);
        stm.write<uint8_t>((uint8_t)Sign::Pos);
        stm.write<uint32_t>(105);
        stm.write<uint32_t>(4);
        stm.write<uint32_t>(0);
        stm.write<uint8_t>((uint8_t)Sign::Pos);
        stm.write<uint32_t>(204);
        stm.write<uint32_t>(5);
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(2);
        // 1
        stm.write<uint32_t>(8); // id
        stm.write<uint32_t>(1); // time_id
        stm.write<uint32_t>(6); // len
        stm.write<uint32_t>(1); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(1); // chn
        // 2
        stm.write<uint32_t>(65); // id
        stm.write<uint32_t>(2); // time_id
        stm.write<uint32_t>(6); // len
        stm.write<uint32_t>(3); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(2); // chn
        // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
        stm.write<uint32_t>(0);
        // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
        stm.write<uint32_t>(0);
        // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
        stm.write<uint32_t>(0);
        // [default_target: 4B]
        stm.write<uint32_t>(0);
    }
    // [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
    stm.write<uint32_t>(0);
    // [nbackenddatas: 4B][[device_name: NUL-terminated string]
    //                     [size: 4B][data: size B] x nbackenddatas]
    stm.write<uint32_t>(0);

    auto expseq = stm.create_sequence();
    REQUIRE(expseq);

    REQUIRE(expseq->devices().size() == 1);
    auto zynq_dev = expseq->get_device("FPGA1", false);
    REQUIRE(zynq_dev);
    auto zynq_backend = Zynq::Backend::cast(zynq_dev);
    REQUIRE(zynq_backend);
    REQUIRE(zynq_backend->has_generator(0));
    REQUIRE(!zynq_backend->pregenerated(0, true));
    REQUIRE(!zynq_backend->pregenerated(0, false));

    expseq->init_run();
    expseq->pre_run();
    {
        auto bc = zynq_backend->get_bytecode(0, true);
        auto bc_ver = zynq_backend->get_bytecode_version(0, true);
        Checker checker(bc.data(), bc.size());
        INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
        REQUIRE(checker.cmp<int64_t>(0)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(0x1)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq2{7, 4, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetAmp{11, 20, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        // Sequence
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 97}));
        // Finish
        REQUIRE(checker.is_end());
    }
    expseq->start();
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 0);

    expseq->init_run();
    expseq->host_seq.set_global(0, 0.8149072527885437);
    expseq->host_seq.set_global(1, 2.9304029304029303e-8);
    expseq->host_seq.set_global(2, 1);
    expseq->pre_run();
    {
        auto bc = zynq_backend->get_bytecode(0, true);
        auto bc_ver = zynq_backend->get_bytecode_version(0, true);
        Checker checker(bc.data(), bc.size());
        INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
        REQUIRE(checker.cmp<int64_t>(0)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(0x1)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq2{7, 4, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetAmp{11, 20, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        // Sequence
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 1246}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq4{9, 4, 500000}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq4{9, 4, 500000}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq4{9, 4, 500000}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq4{9, 4, 500000}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq4{9, 4, 500000}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSAmp{10, 0, 20, 0xfff}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetAmp{11, 20, 0x44})); // -60
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq4{9, 4, 1500000}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSAmp{10, 0, 20, 0xf4b}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq4{9, 4, 1000000}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSAmp{10, 0, 20, 0xed3}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq4{9, 4, 1000000}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSAmp{10, 0, 20, 0xe5b}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq4{9, 4, 1000000}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSAmp{10, 0, 20, 0xde3}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq4{9, 4, 1000000}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSAmp{10, 0, 20, 0xd6b}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq4{9, 4, 2000000}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSAmp{10, 0, 20, 0xcf3}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetAmp{11, 20, 0x44})); // -60
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetAmp{11, 20, 0x44})); // -60
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetAmp{11, 20, 0x44})); // -60
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetAmp{11, 20, 0x44})); // -60
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetAmp{11, 20, 0x44})); // -60
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSAmp{10, 0, 20, 0xb4f}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 100}));
        REQUIRE(checker.is_end());
    }
    expseq->start();
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 0);
}

TEST_CASE("ttl_mgr") {
    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(9);
    // 1
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstBool);
    stm.write<bool>(false);
    // 2
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstBool);
    stm.write<bool>(true);
    // 3
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(100000);
    // 4
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(3000000);
    // 5
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(8000000);
    // 6
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(8100000);
    // 7
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(12000000);
    // 8
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(15000000);
    // 9
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(15500000);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(2);
    stm.write_string("FPGA1/TTL1");
    stm.write_string("FPGA1/TTL2");
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(1);
    stm.write<uint32_t>(2);
    stm.write(IR::Type::Bool);
    stm.write<bool>(true);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(0);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(0);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(1);
    {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(7);
        stm.write<uint8_t>((uint8_t)Sign::Pos);
        stm.write<uint32_t>(105);
        stm.write<uint32_t>(3);
        stm.write<uint32_t>(0);
        stm.write<uint8_t>((uint8_t)Sign::Pos);
        stm.write<uint32_t>(204);
        stm.write<uint32_t>(4);
        stm.write<uint32_t>(0);
        stm.write<uint8_t>((uint8_t)Sign::Pos);
        stm.write<uint32_t>(301);
        stm.write<uint32_t>(5);
        stm.write<uint32_t>(0);
        stm.write<uint8_t>((uint8_t)Sign::Pos);
        stm.write<uint32_t>(404);
        stm.write<uint32_t>(6);
        stm.write<uint32_t>(0);
        stm.write<uint8_t>((uint8_t)Sign::Pos);
        stm.write<uint32_t>(504);
        stm.write<uint32_t>(7);
        stm.write<uint32_t>(0);
        stm.write<uint8_t>((uint8_t)Sign::Pos);
        stm.write<uint32_t>(543);
        stm.write<uint32_t>(8);
        stm.write<uint32_t>(0);
        stm.write<uint8_t>((uint8_t)Sign::Pos);
        stm.write<uint32_t>(104);
        stm.write<uint32_t>(9);
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(7);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(14);
        // 1
        stm.write<uint32_t>(1); // id
        stm.write<uint32_t>(1); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(2); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(1); // chn
        // 2
        stm.write<uint32_t>(2); // id
        stm.write<uint32_t>(2); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(1); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(1); // chn
        // 3
        stm.write<uint32_t>(3); // id
        stm.write<uint32_t>(3); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(2); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(1); // chn
        // 4
        stm.write<uint32_t>(4); // id
        stm.write<uint32_t>(4); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(1); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(1); // chn
        // 5
        stm.write<uint32_t>(5); // id
        stm.write<uint32_t>(5); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(2); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(1); // chn
        // 6
        stm.write<uint32_t>(6); // id
        stm.write<uint32_t>(6); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(1); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(1); // chn
        // 7
        stm.write<uint32_t>(7); // id
        stm.write<uint32_t>(7); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(2); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(1); // chn
        // 8
        stm.write<uint32_t>(8); // id
        stm.write<uint32_t>(1); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(1); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(2); // chn
        // 9
        stm.write<uint32_t>(9); // id
        stm.write<uint32_t>(2); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(2); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(2); // chn
        // 10
        stm.write<uint32_t>(10); // id
        stm.write<uint32_t>(3); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(1); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(2); // chn
        // 11
        stm.write<uint32_t>(11); // id
        stm.write<uint32_t>(4); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(2); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(2); // chn
        // 12
        stm.write<uint32_t>(12); // id
        stm.write<uint32_t>(5); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(1); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(2); // chn
        // 13
        stm.write<uint32_t>(13); // id
        stm.write<uint32_t>(6); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(2); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(2); // chn
        // 14
        stm.write<uint32_t>(14); // id
        stm.write<uint32_t>(7); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(1); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(2); // chn
        // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
        stm.write<uint32_t>(0);
        // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
        stm.write<uint32_t>(0);
        // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
        stm.write<uint32_t>(0);
        // [default_target: 4B]
        stm.write<uint32_t>(0);
    }
    // [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
    stm.write<uint32_t>(0);
    // [nbackenddatas: 4B][[device_name: NUL-terminated string]
    //                     [size: 4B][data: size B] x nbackenddatas]
    stm.write<uint32_t>(1);
    stm.write_string("FPGA1");
    auto szpos = stm.tellp();
    stm.write<uint32_t>(0);
    // [magic <"ZYNQZYNQ">: 8B]
    stm.write("ZYNQZYNQ", 8);
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nttl_mgrs: 1B]
    stm.write<uint8_t>(2);
    // [[chn_id: 4B][off_delay: 4B][on_delay: 4B]
    //  [skip_time: 4B][min_time: 4B][off_val: 1B] x nttl_mgrs]
    stm.write<uint32_t>(1);
    stm.write<uint32_t>(300000);
    stm.write<uint32_t>(200000);
    stm.write<uint32_t>(1000000);
    stm.write<uint32_t>(1300000);
    stm.write<bool>(false);
    stm.write<uint32_t>(2);
    stm.write<uint32_t>(300000);
    stm.write<uint32_t>(200000);
    stm.write<uint32_t>(1000000);
    stm.write<uint32_t>(1300000);
    stm.write<bool>(true);
    auto endpos = stm.tellp();
    stm.seekp(szpos);
    stm.write(uint32_t(endpos - szpos - 4));
    stm.seekp(endpos);

    auto expseq = stm.create_sequence();
    REQUIRE(expseq);

    REQUIRE(expseq->devices().size() == 1);
    auto zynq_dev = expseq->get_device("FPGA1", false);
    REQUIRE(zynq_dev);
    auto zynq_backend = Zynq::Backend::cast(zynq_dev);
    REQUIRE(zynq_backend);
    REQUIRE(!zynq_backend->has_generator(0));
    REQUIRE(zynq_backend->pregenerated(0, true));
    REQUIRE(!zynq_backend->pregenerated(0, false));

    {
        auto bc = zynq_backend->get_bytecode(0, true);
        auto bc_ver = zynq_backend->get_bytecode_version(0, true);
        Checker checker(bc.data(), bc.size());
        INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
        REQUIRE(checker.cmp<int64_t>(15500)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(0x7)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 1}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        // Sequence
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 266}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 1, 2}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 506}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 1, 2}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 126}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 1, 2}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 266}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 1, 2}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 97}));
        REQUIRE(checker.is_end());
    }

    expseq->init_run();
    expseq->pre_run();
    expseq->start();
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 0);

    expseq->init_run();
    expseq->pre_run();
    expseq->start();
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 0);
}

TEST_CASE("offset") {
    for (int offset = 300; offset >= 0; offset -= 100) {
        // Make sure the last one is offset == 0
        // to match the expectation on the global config for other tests.
        if (offset == 0) {
            manager().load_config_string(config_str);
        }
        else {
            auto config = config_str + "      seq_delay_ms: "s +
                std::to_string(offset * 1e-5) + "\n"s;
            manager().load_config_string(config.data());
        }
        ostream stm;
        // [version <0>: 1B]
        stm.write<uint8_t>(0);
        // [nnodes: 4B]
        // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
        //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
        stm.write<uint32_t>(0);
        // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
        stm.write<uint32_t>(2);
        stm.write_string("FPGA1/DDS2/FREQ");
        stm.write_string("FPGA1/TTL3");
        // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
        stm.write<uint32_t>(0);
        // [nslots: 4B][[Type: 1B] x nslots]
        stm.write<uint32_t>(0);
        // [nnoramp: 4B][[chnid: 4B] x nnoramp]
        stm.write<uint32_t>(0);
        // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
        stm.write<uint32_t>(2);
        {
            // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
            stm.write<uint32_t>(0);
            // [nendtimes: 4B][[time_id: 4B] x nendtimes]
            stm.write<uint32_t>(0);
            // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
            stm.write<uint32_t>(0);
            // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
            stm.write<uint32_t>(0);
            // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
            stm.write<uint32_t>(0);
            // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
            stm.write<uint32_t>(0);
            // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
            stm.write<uint32_t>(0);
            // [default_target: 4B]
            stm.write<uint32_t>(2);
        }
        {
            // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
            stm.write<uint32_t>(0);
            // [nendtimes: 4B][[time_id: 4B] x nendtimes]
            stm.write<uint32_t>(0);
            // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
            stm.write<uint32_t>(0);
            // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
            stm.write<uint32_t>(0);
            // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
            stm.write<uint32_t>(0);
            // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
            stm.write<uint32_t>(0);
            // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
            stm.write<uint32_t>(0);
            // [default_target: 4B]
            stm.write<uint32_t>(0);
        }
        // [ndatas: 4B][[ndouble: 4B][data: 8B x ndouble] x ndatas]
        stm.write<uint32_t>(0);
        // [nbackenddatas: 4B][[device_name: NUL-terminated string]
        //                     [size: 4B][data: size B] x nbackenddatas]
        stm.write<uint32_t>(0);

        auto expseq = stm.create_sequence();
        REQUIRE(expseq);

        REQUIRE(expseq->devices().size() == 1);
        auto zynq_dev = expseq->get_device("FPGA1", false);
        REQUIRE(zynq_dev);
        auto zynq_backend = Zynq::Backend::cast(zynq_dev);
        REQUIRE(zynq_backend);
        REQUIRE(!zynq_backend->has_generator(0));
        REQUIRE(zynq_backend->pregenerated(0, true));
        REQUIRE(!zynq_backend->pregenerated(0, false));
        REQUIRE(!zynq_backend->has_generator(1));
        REQUIRE(!zynq_backend->pregenerated(1, true));
        REQUIRE(zynq_backend->pregenerated(1, false));

        // Make sure the wait time fits in a `wait2`
        REQUIRE(offset + 97 < 2048);
        {
            auto bc = zynq_backend->get_bytecode(0, true);
            auto bc_ver = zynq_backend->get_bytecode_version(0, true);
            Checker checker(bc.data(), bc.size());
            INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
            REQUIRE(checker.cmp<int64_t>(0)); // len_ns
            REQUIRE(checker.cmp<uint32_t>(9)); // ttl_mask
            // Startup
            REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::DDSDetFreq2{7, 2, 0}));
            REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
            REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
            REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
            REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
            // Finish
            REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, uint16_t(offset == 0 ?
                                                                              97 : offset + 96)}));
            REQUIRE(checker.is_end());
        }

        {
            auto bc = zynq_backend->get_bytecode(1, false);
            auto bc_ver = zynq_backend->get_bytecode_version(1, false);
            Checker checker(bc.data(), bc.size());
            INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
            REQUIRE(checker.cmp<int64_t>(0)); // len_ns
            REQUIRE(checker.cmp<uint32_t>(9)); // ttl_mask
            // Startup
            REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
            REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
            REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
            REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
            // Finish
            REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, uint16_t(offset == 0 ?
                                                                              97 : offset + 96)}));
            REQUIRE(checker.is_end());
        }

        expseq->init_run();
        expseq->pre_run();
        expseq->start();
        REQUIRE(expseq->wait(20000)); // 20s timeout
        REQUIRE(expseq->post_run() == 2);
        expseq->pre_run();
        expseq->start();
        REQUIRE(expseq->wait(20000)); // 20s timeout
        REQUIRE(expseq->post_run() == 0);
    }
}
