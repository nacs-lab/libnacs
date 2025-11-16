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
#include "../../lib/nacs-seq/builder.h"
#include "../../lib/nacs-seq/nidaq/backend.h"
#include "../../lib/nacs-seq/zynq/backend.h"

using namespace std::literals::string_literals;

using active_time_t = std::vector<std::pair<int64_t,int64_t>>;

static void check_channels(NiDAQ::Backend *backend,
                           std::vector<std::pair<std::string,uint32_t>> channels)
{
    uint32_t nchns;
    auto info = backend->get_channel_info(&nchns);
    REQUIRE(channels.size() == nchns);
    for (size_t i = 0; i < nchns; i++) {
        INFO(i);
        REQUIRE(info[i].first == channels[i].second);
        REQUIRE(info[i].second == channels[i].first);
    }
}

static void check_nsamples(NiDAQ::Backend *backend, uint32_t cur_seq_id, size_t nsamples)
{
    uint32_t nchns;
    backend->get_channel_info(&nchns);
    size_t data_sz;
    backend->get_data(cur_seq_id, &data_sz);
    REQUIRE(nsamples * nchns == data_sz);
}

static void check_data_range(NiDAQ::Backend *backend, uint32_t cur_seq_id, uint32_t chn,
                             size_t begin, size_t end, double val)
{
    uint32_t nchns;
    backend->get_channel_info(&nchns);
    size_t data_sz;
    auto data = backend->get_data(cur_seq_id, &data_sz);
    REQUIRE(data_sz % nchns == 0);
    size_t nsamples = data_sz / nchns;
    REQUIRE(begin <= end);
    REQUIRE(end <= nsamples);
    REQUIRE(chn < nchns);
    auto offset = chn * nsamples;
    for (size_t i = begin; i < end; i++) {
        INFO(i);
        REQUIRE(data[i + offset] == Approx(val));
    }
}

template<typename Func>
static void check_data_range(NiDAQ::Backend *backend, uint32_t cur_seq_id, uint32_t chn,
                             size_t begin, size_t end, Func &&func)
{
    uint32_t nchns;
    backend->get_channel_info(&nchns);
    size_t data_sz;
    auto data = backend->get_data(cur_seq_id, &data_sz);
    REQUIRE(data_sz % nchns == 0);
    size_t nsamples = data_sz / nchns;
    REQUIRE(begin <= end);
    REQUIRE(end <= nsamples);
    REQUIRE(chn < nchns);
    auto offset = chn * nsamples;
    for (size_t i = begin; i < end; i++) {
        INFO(i);
        REQUIRE(data[i + offset] == Approx(func(i)));
    }
}

extern "C" {

const std::pair<uint32_t,const char*> *nacs_seq_manager_expseq_get_nidaq_channel_info(
    Manager::ExpSeq *expseq, const char *name, uint32_t *sz);
const double *nacs_seq_manager_expseq_get_nidaq_data(
    Manager::ExpSeq *expseq, const char *name, size_t *sz);

}

static void check_c_api(NiDAQ::Backend *backend, Manager::ExpSeq *expseq)
{
    uint32_t nchns;
    auto chns = backend->get_channel_info(&nchns);
    uint32_t nchns_c;
    auto chns_c = nacs_seq_manager_expseq_get_nidaq_channel_info(expseq, backend->name().data(),
                                                                 &nchns_c);
    REQUIRE(chns);
    REQUIRE(chns_c);
    REQUIRE(nchns == nchns_c);
    for (uint32_t i = 0; i < nchns; i++) {
        INFO(i);
        REQUIRE(chns[i] == chns_c[i]);
    }

    auto cur_seq_id = expseq->host_seq.cur_seq_idx();
    size_t data_sz;
    auto data = backend->get_data(cur_seq_id, &data_sz);
    size_t data_sz_c;
    auto data_c = nacs_seq_manager_expseq_get_nidaq_data(expseq, backend->name().data(),
                                                         &data_sz_c);
    REQUIRE(data_sz == data_sz_c);
    for (size_t i = 0; i < data_sz; i++) {
        INFO(i);
        REQUIRE(data[i] == data_c[i]);
    }
}

static const char *const config_str =
    "tick_per_sec: 1000000000000\n"
    "devices:\n"
    "  NIDAQ:\n"
    "    backend: nidaq\n"
    "    config:\n"
    "      step_size: 2000000\n"
    "      clock_device: FPGA1\n"
    "  FPGA1:\n"
    "    backend: zynq\n"
    "    config:\n"
    "      start_ttl_chn: 0\n";

static Manager &manager()
{
    static Manager mgr;
    static Call init([&] {
        mgr.dummy_mode = true; // For the Zynq backend
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
    stm.write_string("NIDAQ/Dev2/3");
    stm.write_string("NIDAQ/Dev0/1");
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(1);
    stm.write<uint32_t>(1);
    stm.write(IR::Type::Float64);
    stm.write<double>(3.4);
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
    REQUIRE(!zynq_backend->has_generator(0));
    REQUIRE(zynq_backend->pregenerated(0, true));
    REQUIRE(!zynq_backend->pregenerated(0, false));
    REQUIRE(!zynq_backend->has_generator(1));
    REQUIRE(!zynq_backend->pregenerated(1, true));
    REQUIRE(zynq_backend->pregenerated(1, false));

    auto nidaq_dev = expseq->get_device("NIDAQ", false);
    REQUIRE(nidaq_dev);
    auto nidaq_backend = NiDAQ::Backend::cast(nidaq_dev);
    REQUIRE(nidaq_backend);
    check_channels(nidaq_backend, {{"Dev0", 1}, {"Dev2", 3}});
    auto status = nidaq_backend->get_generate_status(0);
    REQUIRE(status.fill_init == std::vector<size_t>({0, 0}));
    REQUIRE(status.nsamples == 1);
    REQUIRE(status.first_only);
    REQUIRE(status.all_time_known);
    REQUIRE(status.all_val_known);
    status = nidaq_backend->get_generate_status(1);
    REQUIRE(status.fill_init == std::vector<size_t>({0, 0}));
    REQUIRE(status.nsamples == 0);
    REQUIRE(!status.first_only);
    REQUIRE(status.all_time_known);
    REQUIRE(status.all_val_known);

    {
        auto bc = zynq_backend->get_bytecode(0, true);
        auto bc_ver = zynq_backend->get_bytecode_version(0, true);
        Checker checker(bc.data(), bc.size());
        INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
        REQUIRE(checker.cmp<int64_t>(0)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 0, 0, 0}));
        // Sequence
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 99}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 194}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 255}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 100}));
        REQUIRE(checker.is_end());
    }

    {
        auto bc = zynq_backend->get_bytecode(1, false);
        auto bc_ver = zynq_backend->get_bytecode_version(1, false);
        Checker checker(bc.data(), bc.size());
        INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
        REQUIRE(checker.cmp<int64_t>(0)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 296}));
        REQUIRE(checker.is_end());
    }

    expseq->init_run();
    expseq->pre_run();
    expseq->start();
    check_nsamples(nidaq_backend, 0, 1);
    check_c_api(nidaq_backend, expseq.get());
    check_data_range(nidaq_backend, 0, 0, 0, 1, 0.0);
    check_data_range(nidaq_backend, 0, 1, 0, 1, 3.4);
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 2);
    expseq->pre_run();
    check_nsamples(nidaq_backend, 1, 0);
    check_c_api(nidaq_backend, expseq.get());
    expseq->start();
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 0);
}

TEST_CASE("Short sequence") {
    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(4);
    // 1
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(400000000);
    // 2
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(1.9);
    // 3
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(600002000);
    // 4
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(-2.3);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(2);
    stm.write_string("NIDAQ/Dev2/3");
    stm.write_string("NIDAQ/Dev0/1");
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(2);
    stm.write<uint32_t>(1);
    stm.write(IR::Type::Float64);
    stm.write<double>(2.6);
    stm.write<uint32_t>(2);
    stm.write(IR::Type::Float64);
    stm.write<double>(1.2);
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
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(0);
        stm.write<uint8_t>((uint8_t)Sign::Pos);
        stm.write<uint32_t>(204);
        stm.write<uint32_t>(3);
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(2);
        // 1
        stm.write<uint32_t>(1); // id
        stm.write<uint32_t>(1); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(2); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(2); // chn
        // 2
        stm.write<uint32_t>(2); // id
        stm.write<uint32_t>(2); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(4); // val
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

    auto nidaq_dev = expseq->get_device("NIDAQ", false);
    REQUIRE(nidaq_dev);
    auto nidaq_backend = NiDAQ::Backend::cast(nidaq_dev);
    REQUIRE(nidaq_backend);
    check_channels(nidaq_backend, {{"Dev0", 1}, {"Dev2", 3}});
    auto status = nidaq_backend->get_generate_status(0);
    REQUIRE(status.fill_init == std::vector<size_t>({0, 0}));
    REQUIRE(status.nsamples == 302);
    REQUIRE(!status.first_only);
    REQUIRE(status.all_time_known);
    REQUIRE(status.all_val_known);

    {
        auto bc = zynq_backend->get_bytecode(0, true);
        auto bc_ver = zynq_backend->get_bytecode_version(0, true);
        Checker checker(bc.data(), bc.size());
        INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
        REQUIRE(checker.cmp<int64_t>(0)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 0, 0, 0}));
        // Sequence
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 99}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 60394}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 255}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 100}));
        REQUIRE(checker.is_end());
    }

    expseq->init_run();
    expseq->pre_run();
    expseq->start();
    check_nsamples(nidaq_backend, 0, 302);
    check_c_api(nidaq_backend, expseq.get());
    check_data_range(nidaq_backend, 0, 0, 0, 200, 1.2);
    check_data_range(nidaq_backend, 0, 0, 201, 302, 1.9);
    check_data_range(nidaq_backend, 0, 1, 0, 301, 2.6);
    check_data_range(nidaq_backend, 0, 1, 301, 302, -2.3);
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 0);
}

TEST_CASE("Short sequence with ramp") {
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
    stm.write<double>(400000000);
    // 2
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(1.9);
    // 3
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(600800000);
    // 4
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(-2.3);
    // 5
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(201600000);
    // 6
    stm.write(Builder::OpCode::Div);
    stm.write(Builder::ArgType::Arg);
    stm.write<uint32_t>(0);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(200000000);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(2);
    stm.write_string("NIDAQ/Dev2/3");
    stm.write_string("NIDAQ/Dev0/1");
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(2);
    stm.write<uint32_t>(1);
    stm.write(IR::Type::Float64);
    stm.write<double>(2.6);
    stm.write<uint32_t>(2);
    stm.write(IR::Type::Float64);
    stm.write<double>(1.2);
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
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(0);
        stm.write<uint8_t>((uint8_t)Sign::Pos);
        stm.write<uint32_t>(204);
        stm.write<uint32_t>(3);
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(2);
        // 1
        stm.write<uint32_t>(1); // id
        stm.write<uint32_t>(1); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(2); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(2); // chn
        // 2
        stm.write<uint32_t>(2); // id
        stm.write<uint32_t>(2); // time_id
        stm.write<uint32_t>(5); // len
        stm.write<uint32_t>(6); // val
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

    auto nidaq_dev = expseq->get_device("NIDAQ", false);
    REQUIRE(nidaq_dev);
    auto nidaq_backend = NiDAQ::Backend::cast(nidaq_dev);
    REQUIRE(nidaq_backend);
    check_channels(nidaq_backend, {{"Dev0", 1}, {"Dev2", 3}});
    auto status = nidaq_backend->get_generate_status(0);
    REQUIRE(status.fill_init == std::vector<size_t>({0, 0}));
    REQUIRE(status.nsamples == 403);
    REQUIRE(!status.first_only);
    REQUIRE(status.all_time_known);
    REQUIRE(status.all_val_known);

    {
        auto bc = zynq_backend->get_bytecode(0, true);
        auto bc_ver = zynq_backend->get_bytecode_version(0, true);
        Checker checker(bc.data(), bc.size());
        INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
        REQUIRE(checker.cmp<int64_t>(0)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 0, 0, 0}));
        // Sequence
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 99}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 1, 0x275a}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 2}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 255}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 100}));
        REQUIRE(checker.is_end());
    }

    expseq->init_run();
    expseq->pre_run();
    expseq->start();
    check_nsamples(nidaq_backend, 0, 403);
    check_c_api(nidaq_backend, expseq.get());
    check_data_range(nidaq_backend, 0, 0, 0, 200, 1.2);
    check_data_range(nidaq_backend, 0, 0, 201, 403, 1.9);
    check_data_range(nidaq_backend, 0, 1, 0, 301, 2.6);
    check_data_range(nidaq_backend, 0, 1, 301, 402, [&] (auto t) {
        return (double(t) * 2000000 - 600800000) / 200000000;
    });
    check_data_range(nidaq_backend, 0, 1, 402, 403, 1.008);
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 0);
}

TEST_CASE("Ramp with clip") {
    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(7);
    // 1
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(400000000);
    // 2
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(1.9);
    // 3
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(600800000);
    // 4
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(-2.3);
    // 5
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(201600000);
    // 6
    stm.write(Builder::OpCode::Div);
    stm.write(Builder::ArgType::Arg);
    stm.write<uint32_t>(0);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(5000000);
    // 7
    stm.write(Builder::OpCode::Sub);
    stm.write(Builder::ArgType::Node);
    stm.write<uint32_t>(6);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(20);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(2);
    stm.write_string("NIDAQ/Dev2/3");
    stm.write_string("NIDAQ/Dev0/1");
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(2);
    stm.write<uint32_t>(1);
    stm.write(IR::Type::Float64);
    stm.write<double>(2.6);
    stm.write<uint32_t>(2);
    stm.write(IR::Type::Float64);
    stm.write<double>(1.2);
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
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(0);
        stm.write<uint8_t>((uint8_t)Sign::Pos);
        stm.write<uint32_t>(204);
        stm.write<uint32_t>(3);
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(2);
        // 1
        stm.write<uint32_t>(1); // id
        stm.write<uint32_t>(1); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(2); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(2); // chn
        // 2
        stm.write<uint32_t>(2); // id
        stm.write<uint32_t>(2); // time_id
        stm.write<uint32_t>(5); // len
        stm.write<uint32_t>(7); // val
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

    auto nidaq_dev = expseq->get_device("NIDAQ", false);
    REQUIRE(nidaq_dev);
    auto nidaq_backend = NiDAQ::Backend::cast(nidaq_dev);
    REQUIRE(nidaq_backend);
    check_channels(nidaq_backend, {{"Dev0", 1}, {"Dev2", 3}});
    auto status = nidaq_backend->get_generate_status(0);
    REQUIRE(status.fill_init == std::vector<size_t>({0, 0}));
    REQUIRE(status.nsamples == 403);
    REQUIRE(!status.first_only);
    REQUIRE(status.all_time_known);
    REQUIRE(status.all_val_known);

    {
        auto bc = zynq_backend->get_bytecode(0, true);
        auto bc_ver = zynq_backend->get_bytecode_version(0, true);
        Checker checker(bc.data(), bc.size());
        INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
        REQUIRE(checker.cmp<int64_t>(0)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 0, 0, 0}));
        // Sequence
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 99}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 1, 0x275a}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 2}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 255}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 100}));
        REQUIRE(checker.is_end());
    }

    expseq->init_run();
    expseq->pre_run();
    expseq->start();
    check_nsamples(nidaq_backend, 0, 403);
    check_c_api(nidaq_backend, expseq.get());
    check_data_range(nidaq_backend, 0, 0, 0, 200, 1.2);
    check_data_range(nidaq_backend, 0, 0, 201, 403, 1.9);
    check_data_range(nidaq_backend, 0, 1, 0, 301, 2.6);
    check_data_range(nidaq_backend, 0, 1, 301, 402, [&] (auto t) {
        return bound(-10, (double(t) * 20 - 6008) / 50 - 20, 10);
    });
    // Check the start and end points of the ramp twice
    // to make sure our ramp is actually clipped.
    check_data_range(nidaq_backend, 0, 1, 301, 303, -10.0);
    check_data_range(nidaq_backend, 0, 1, 400, 402, 10.0);
    check_data_range(nidaq_backend, 0, 1, 402, 403, 10.0);
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 0);
}

TEST_CASE("Unknown time") {
    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(7);
    // 1
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(400000000);
    // 2
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(1.9);
    // 3
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(600800000);
    // 4
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(-2.3);
    // 5
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::Global);
    stm.write<uint32_t>(0);
    // 6
    stm.write(Builder::OpCode::Div);
    stm.write(Builder::ArgType::Arg);
    stm.write<uint32_t>(0);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(5000000);
    // 7
    stm.write(Builder::OpCode::Sub);
    stm.write(Builder::ArgType::Node);
    stm.write<uint32_t>(6);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(20);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(2);
    stm.write_string("NIDAQ/Dev2/3");
    stm.write_string("NIDAQ/Dev0/1");
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(2);
    stm.write<uint32_t>(1);
    stm.write(IR::Type::Float64);
    stm.write<double>(2.6);
    stm.write<uint32_t>(2);
    stm.write(IR::Type::Float64);
    stm.write<double>(1.2);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(1);
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
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(0);
        stm.write<uint8_t>((uint8_t)Sign::Pos);
        stm.write<uint32_t>(204);
        stm.write<uint32_t>(3);
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(2);
        // 1
        stm.write<uint32_t>(1); // id
        stm.write<uint32_t>(1); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(2); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(2); // chn
        // 2
        stm.write<uint32_t>(2); // id
        stm.write<uint32_t>(2); // time_id
        stm.write<uint32_t>(5); // len
        stm.write<uint32_t>(7); // val
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

    auto nidaq_dev = expseq->get_device("NIDAQ", false);
    REQUIRE(nidaq_dev);
    auto nidaq_backend = NiDAQ::Backend::cast(nidaq_dev);
    REQUIRE(nidaq_backend);
    check_channels(nidaq_backend, {{"Dev0", 1}, {"Dev2", 3}});
    auto status = nidaq_backend->get_generate_status(0);
    REQUIRE(status.fill_init == std::vector<size_t>());
    REQUIRE(status.nsamples == size_t(-1));
    REQUIRE(!status.first_only);
    REQUIRE(!status.all_time_known);
    REQUIRE(!status.all_val_known);

    expseq->init_run();
    expseq->host_seq.set_global(0, 201600000);
    expseq->pre_run();
    {
        auto bc = zynq_backend->get_bytecode(0, true);
        auto bc_ver = zynq_backend->get_bytecode_version(0, true);
        Checker checker(bc.data(), bc.size());
        INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
        REQUIRE(checker.cmp<int64_t>(0)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 0, 0, 0}));
        // Sequence
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 99}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 1, 0x275a}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 2}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 255}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 100}));
        REQUIRE(checker.is_end());
    }
    expseq->start();
    check_nsamples(nidaq_backend, 0, 403);
    check_c_api(nidaq_backend, expseq.get());
    check_data_range(nidaq_backend, 0, 0, 0, 200, 1.2);
    check_data_range(nidaq_backend, 0, 0, 201, 403, 1.9);
    check_data_range(nidaq_backend, 0, 1, 0, 301, 2.6);
    check_data_range(nidaq_backend, 0, 1, 301, 402, [&] (auto t) {
        return bound(-10, (double(t) * 20 - 6008) / 50 - 20, 10);
    });
    // Check the start and end points of the ramp twice
    // to make sure our ramp is actually clipped.
    check_data_range(nidaq_backend, 0, 1, 301, 303, -10.0);
    check_data_range(nidaq_backend, 0, 1, 400, 402, 10.0);
    check_data_range(nidaq_backend, 0, 1, 402, 403, 10.0);
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 0);
}

TEST_CASE("Unknown value") {
    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(7);
    // 1
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(400000000);
    // 2
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(1.9);
    // 3
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(600800000);
    // 4
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(-2.3);
    // 5
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(201600000);
    // 6
    stm.write(Builder::OpCode::Div);
    stm.write(Builder::ArgType::Arg);
    stm.write<uint32_t>(0);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(5000000);
    // 7
    stm.write(Builder::OpCode::Sub);
    stm.write(Builder::ArgType::Node);
    stm.write<uint32_t>(6);
    stm.write(Builder::ArgType::Global);
    stm.write<uint32_t>(0);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(2);
    stm.write_string("NIDAQ/Dev2/3");
    stm.write_string("NIDAQ/Dev0/1");
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(2);
    stm.write<uint32_t>(1);
    stm.write(IR::Type::Float64);
    stm.write<double>(2.6);
    stm.write<uint32_t>(2);
    stm.write(IR::Type::Float64);
    stm.write<double>(1.2);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(1);
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
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(0);
        stm.write<uint8_t>((uint8_t)Sign::Pos);
        stm.write<uint32_t>(204);
        stm.write<uint32_t>(3);
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(2);
        // 1
        stm.write<uint32_t>(1); // id
        stm.write<uint32_t>(1); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(2); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(2); // chn
        // 2
        stm.write<uint32_t>(2); // id
        stm.write<uint32_t>(2); // time_id
        stm.write<uint32_t>(5); // len
        stm.write<uint32_t>(7); // val
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

    auto nidaq_dev = expseq->get_device("NIDAQ", false);
    REQUIRE(nidaq_dev);
    auto nidaq_backend = NiDAQ::Backend::cast(nidaq_dev);
    REQUIRE(nidaq_backend);
    check_channels(nidaq_backend, {{"Dev0", 1}, {"Dev2", 3}});
    auto status = nidaq_backend->get_generate_status(0);
    REQUIRE(status.fill_init == std::vector<size_t>());
    REQUIRE(status.nsamples == size_t(-1));
    REQUIRE(!status.first_only);
    REQUIRE(status.all_time_known);
    REQUIRE(!status.all_val_known);

    {
        auto bc = zynq_backend->get_bytecode(0, true);
        auto bc_ver = zynq_backend->get_bytecode_version(0, true);
        Checker checker(bc.data(), bc.size());
        INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
        REQUIRE(checker.cmp<int64_t>(0)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 0, 0, 0}));
        // Sequence
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 99}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 1, 0x275a}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 2}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 255}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 100}));
        REQUIRE(checker.is_end());
    }

    expseq->init_run();
    expseq->host_seq.set_global(0, 20);
    expseq->pre_run();
    expseq->start();
    check_nsamples(nidaq_backend, 0, 403);
    check_c_api(nidaq_backend, expseq.get());
    check_data_range(nidaq_backend, 0, 0, 0, 200, 1.2);
    check_data_range(nidaq_backend, 0, 0, 201, 403, 1.9);
    check_data_range(nidaq_backend, 0, 1, 0, 301, 2.6);
    check_data_range(nidaq_backend, 0, 1, 301, 402, [&] (auto t) {
        return bound(-10, (double(t) * 20 - 6008) / 50 - 20, 10);
    });
    // Check the start and end points of the ramp twice
    // to make sure our ramp is actually clipped.
    check_data_range(nidaq_backend, 0, 1, 301, 303, -10.0);
    check_data_range(nidaq_backend, 0, 1, 400, 402, 10.0);
    check_data_range(nidaq_backend, 0, 1, 402, 403, 10.0);
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 0);
}

TEST_CASE("Unknown init") {
    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(8);
    // 1
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(400000000);
    // 2
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(1.9);
    // 3
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(600800000);
    // 4
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(-2.3);
    // 5
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(201600000);
    // 6
    stm.write(Builder::OpCode::Div);
    stm.write(Builder::ArgType::Arg);
    stm.write<uint32_t>(0);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(5000000);
    // 7
    stm.write(Builder::OpCode::Sub);
    stm.write(Builder::ArgType::Node);
    stm.write<uint32_t>(6);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(20);
    // 8
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::Global);
    stm.write<uint32_t>(0);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(2);
    stm.write_string("NIDAQ/Dev2/3");
    stm.write_string("NIDAQ/Dev0/1");
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(2);
    stm.write<uint32_t>(1);
    stm.write(IR::Type::Float64);
    stm.write<double>(2.6);
    stm.write<uint32_t>(2);
    stm.write(IR::Type::Float64);
    stm.write<double>(1.2);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(1);
    stm.write(IR::Type::Bool);
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
        stm.write<uint32_t>(2);
        stm.write<uint8_t>((uint8_t)Sign::Pos);
        stm.write<uint32_t>(105);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(0);
        stm.write<uint8_t>((uint8_t)Sign::Pos);
        stm.write<uint32_t>(204);
        stm.write<uint32_t>(3);
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(0);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(2);
        // 1
        stm.write<uint32_t>(1); // id
        stm.write<uint32_t>(1); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(2); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(2); // chn
        // 2
        stm.write<uint32_t>(2); // id
        stm.write<uint32_t>(2); // time_id
        stm.write<uint32_t>(5); // len
        stm.write<uint32_t>(7); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(1); // chn
        // [nmeasures: 4B][[id: 4B][time_id: 4B][chn: 4B] x nmeasures]
        stm.write<uint32_t>(0);
        // [nassigns: 4B][[assign_id: 4B][global_id: 4B][val: 4B] x nassigns]
        stm.write<uint32_t>(0);
        // [nbranches: 4B][[branch_id: 4B][target_id: 4B][cond: 4B] x nbranches]
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(123);
        stm.write<uint32_t>(2);
        stm.write<uint32_t>(8);
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
    REQUIRE(!zynq_backend->has_generator(0));
    REQUIRE(zynq_backend->pregenerated(0, true));
    REQUIRE(!zynq_backend->pregenerated(0, false));
    REQUIRE(!zynq_backend->has_generator(1));
    REQUIRE(!zynq_backend->pregenerated(1, true));
    REQUIRE(zynq_backend->pregenerated(1, false));

    auto nidaq_dev = expseq->get_device("NIDAQ", false);
    REQUIRE(nidaq_dev);
    auto nidaq_backend = NiDAQ::Backend::cast(nidaq_dev);
    REQUIRE(nidaq_backend);
    check_channels(nidaq_backend, {{"Dev0", 1}, {"Dev2", 3}});
    auto status = nidaq_backend->get_generate_status(0);
    REQUIRE(status.fill_init == std::vector<size_t>({0, 0}));
    REQUIRE(status.nsamples == 1);
    REQUIRE(status.first_only);
    REQUIRE(status.all_time_known);
    REQUIRE(status.all_val_known);
    status = nidaq_backend->get_generate_status(1);
    REQUIRE(status.fill_init == std::vector<size_t>({200, 301}));
    REQUIRE(status.nsamples == 403);
    REQUIRE(!status.first_only);
    REQUIRE(status.all_time_known);
    REQUIRE(status.all_val_known);

    {
        auto bc = zynq_backend->get_bytecode(0, true);
        auto bc_ver = zynq_backend->get_bytecode_version(0, true);
        Checker checker(bc.data(), bc.size());
        INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
        REQUIRE(checker.cmp<int64_t>(0)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 0, 0, 0}));
        // Sequence
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 99}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 194}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 255}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 100}));
        REQUIRE(checker.is_end());
    }

    {
        auto bc = zynq_backend->get_bytecode(1, false);
        auto bc_ver = zynq_backend->get_bytecode_version(1, false);
        Checker checker(bc.data(), bc.size());
        INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
        REQUIRE(checker.cmp<int64_t>(0)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 0, 0, 0}));
        // Sequence
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 99}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 1, 0x275a}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 2}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 255}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 100}));
        REQUIRE(checker.is_end());
    }

    expseq->init_run();
    expseq->pre_run();
    expseq->start();
    check_nsamples(nidaq_backend, 0, 1);
    check_c_api(nidaq_backend, expseq.get());
    check_data_range(nidaq_backend, 0, 0, 0, 1, 1.2);
    check_data_range(nidaq_backend, 0, 1, 0, 1, 2.6);
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 2);
    expseq->pre_run();
    expseq->start();
    check_nsamples(nidaq_backend, 1, 403);
    check_c_api(nidaq_backend, expseq.get());
    check_data_range(nidaq_backend, 1, 0, 0, 200, 1.2);
    check_data_range(nidaq_backend, 1, 0, 201, 403, 1.9);
    check_data_range(nidaq_backend, 1, 1, 0, 301, 2.6);
    check_data_range(nidaq_backend, 1, 1, 301, 402, [&] (auto t) {
        return bound(-10, (double(t) * 20 - 6008) / 50 - 20, 10);
    });
    // Check the start and end points of the ramp twice
    // to make sure our ramp is actually clipped.
    check_data_range(nidaq_backend, 1, 1, 301, 303, -10.0);
    check_data_range(nidaq_backend, 1, 1, 400, 402, 10.0);
    check_data_range(nidaq_backend, 1, 1, 402, 403, 10.0);
    REQUIRE(expseq->wait(20000)); // 20s timeout
    expseq->host_seq.set_global(0, true);
    REQUIRE(expseq->post_run() == 2);
    expseq->pre_run();
    expseq->start();
    check_nsamples(nidaq_backend, 1, 403);
    check_c_api(nidaq_backend, expseq.get());
    check_data_range(nidaq_backend, 1, 0, 0, 403, 1.9);
    check_data_range(nidaq_backend, 1, 1, 0, 301, 10.0);
    check_data_range(nidaq_backend, 1, 1, 301, 402, [&] (auto t) {
        return bound(-10, (double(t) * 20 - 6008) / 50 - 20, 10);
    });
    // Check the start and end points of the ramp twice
    // to make sure our ramp is actually clipped.
    check_data_range(nidaq_backend, 1, 1, 301, 303, -10.0);
    check_data_range(nidaq_backend, 1, 1, 400, 402, 10.0);
    check_data_range(nidaq_backend, 1, 1, 402, 403, 10.0);
    REQUIRE(expseq->wait(20000)); // 20s timeout
    expseq->host_seq.set_global(0, false);
    REQUIRE(expseq->post_run() == 0);
}

TEST_CASE("Specific map") {
    // Test mapping of basic sequence specific values
    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]
    stm.write<uint32_t>(10);
    // 1
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(3e9);
    // 2
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(0);
    // 3
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(1e9);
    // 4
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::Global);
    stm.write<uint32_t>(0);
    // 5
    stm.write(Builder::OpCode::CmpGT);
    stm.write(Builder::ArgType::Global);
    stm.write<uint32_t>(0);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(0.005);
    // 6
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(-1);
    // 7
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(2e9);
    // 8
    stm.write(Builder::OpCode::CmpEQ);
    stm.write(Builder::ArgType::Global);
    stm.write<uint32_t>(0);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(0.005);
    // 9
    stm.write(Builder::OpCode::Div);
    stm.write(Builder::ArgType::Arg);
    stm.write<uint32_t>(0);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(1e12);
    // 10
    stm.write(Builder::OpCode::Add);
    stm.write(Builder::ArgType::Arg);
    stm.write<uint32_t>(1);
    stm.write(Builder::ArgType::Node);
    stm.write<uint32_t>(9);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(1);
    stm.write_string("NIDAQ/Dev1/7");
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(1);
    stm.write<uint32_t>(1);
    stm.write(IR::Type::Float64);
    stm.write<double>(0.465);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(1);
    stm.write(IR::Type::Float64);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(0);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(1);
    {
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(4);
        stm.write<uint8_t>((uint8_t)Sign::Unknown);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(0);
        stm.write<uint8_t>((uint8_t)Sign::Unknown);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(2);
        stm.write<uint32_t>(0);
        stm.write<uint8_t>((uint8_t)Sign::Unknown);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(3);
        stm.write<uint32_t>(0);
        stm.write<uint8_t>((uint8_t)Sign::Unknown);
        stm.write<uint32_t>(0);
        stm.write<uint32_t>(7);
        stm.write<uint32_t>(0);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(1);
        stm.write<uint32_t>(1);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(3);
        // 1
        stm.write<uint32_t>(0); // id
        stm.write<uint32_t>(2); // time_id
        stm.write<uint32_t>(3); // len
        stm.write<uint32_t>(4); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(1); // chn
        // 2
        stm.write<uint32_t>(1); // id
        stm.write<uint32_t>(3); // time_id
        stm.write<uint32_t>(3); // len
        stm.write<uint32_t>(6); // val
        stm.write<uint32_t>(5); // cond
        stm.write<uint32_t>(1); // chn
        // 3
        stm.write<uint32_t>(2); // id
        stm.write<uint32_t>(4); // time_id
        stm.write<uint32_t>(3); // len
        stm.write<uint32_t>(10); // val
        stm.write<uint32_t>(8); // cond
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

    auto nidaq_dev = expseq->get_device("NIDAQ", false);
    REQUIRE(nidaq_dev);
    auto nidaq_backend = NiDAQ::Backend::cast(nidaq_dev);
    REQUIRE(nidaq_backend);
    check_channels(nidaq_backend, {{"Dev1", 7}});
    auto status = nidaq_backend->get_generate_status(0);
    REQUIRE(status.fill_init == std::vector<size_t>());
    REQUIRE(status.nsamples == size_t(-1));
    REQUIRE(!status.first_only);
    REQUIRE(!status.all_time_known);
    REQUIRE(!status.all_val_known);

    // expseq->init_run();
    // expseq->host_seq.set_global(0, 0.001);
    // expseq->pre_run();
    // expseq->start();
    // {
    //     auto bc = zynq_backend->get_bytecode(0, true);
    //     auto bc_ver = zynq_backend->get_bytecode(0, true);
    //     Checker checker(bc.data(), bc.size());
    //     INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
    //     REQUIRE(checker.cmp<int64_t>(3000000)); // len_ns
    //     REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
    //     // Startup
    //     REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
    //     REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
    //     REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
    //     REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 0, 0, 0}));
    //     // Sequence
    //     REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 99}));
    //     REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 194}));
    //     REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 255}));
    //     // Finish
    //     REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 100}));
    //     REQUIRE(checker.is_end());
    // }
    // check_nsamples(nidaq_backend, 0, 1);
    // check_c_api(nidaq_backend, expseq.get());
    // check_data_range(nidaq_backend, 0, 0, 0, 1, 0.001);
    // REQUIRE(expseq->wait(20000)); // 20s timeout
    // REQUIRE(expseq->post_run() == 0);

    // expseq->init_run();
    // expseq->host_seq.set_global(0, 0.009);
    // expseq->pre_run();
    // expseq->start();
    // {
    //     auto bc = zynq_backend->get_bytecode(0, true);
    //     auto bc_ver = zynq_backend->get_bytecode(0, true);
    //     Checker checker(bc.data(), bc.size());
    //     INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
    //     REQUIRE(checker.cmp<int64_t>(3000000)); // len_ns
    //     REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
    //     // Startup
    //     REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
    //     REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
    //     REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
    //     REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 0, 0, 0}));
    //     // Sequence
    //     REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 99}));
    //     REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 1, 0x30ec})); // total wait 100095
    //     REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 2}));
    //     REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 255}));
    //     // Finish
    //     REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 100}));
    //     REQUIRE(checker.is_end());
    // }
    // check_nsamples(nidaq_backend, 0, 501);
    // check_c_api(nidaq_backend, expseq.get());
    // check_data_range(nidaq_backend, 0, 0, 0, 500, 0.009);
    // check_data_range(nidaq_backend, 0, 0, 500, 501, -1.0);
    // REQUIRE(expseq->wait(20000)); // 20s timeout
    // REQUIRE(expseq->post_run() == 0);

    expseq->init_run();
    expseq->host_seq.set_global(0, 0.005);
    expseq->pre_run();
    expseq->start();
    {
        auto bc = zynq_backend->get_bytecode(0, true);
        auto bc_ver = zynq_backend->get_bytecode_version(0, true);
        Checker checker(bc.data(), bc.size());
        INFO(dump_bytecode(bc.data(), bc.size(), bc_ver));
        REQUIRE(checker.cmp<int64_t>(3000000)); // len_ns
        REQUIRE(checker.cmp<uint32_t>(1)); // ttl_mask
        // Startup
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 0, 2500}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 3, 0, 0}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 96}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::TTL2{1, 0, 0, 0}));
        // Sequence
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 99}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait{4, 1, 0x9294})); // total wait 300095
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 2}));
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Clock{5, 0, 255}));
        // Finish
        REQUIRE(checker.cmp(Zynq::ByteCode::Inst_v1::Wait2{5, 1, 100}));
        REQUIRE(checker.is_end());
    }
    check_nsamples(nidaq_backend, 0, 1501);
    check_c_api(nidaq_backend, expseq.get());
    check_data_range(nidaq_backend, 0, 0, 0, 1000, 0.005);
    check_data_range(nidaq_backend, 0, 0, 1000, 1501, [&] (auto t) {
        return 0.005 + double(t - 1000) * 2e-6;
    });
    REQUIRE(expseq->wait(20000)); // 20s timeout
    REQUIRE(expseq->post_run() == 0);
}
