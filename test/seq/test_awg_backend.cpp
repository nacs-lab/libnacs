#define LLVM_DISABLE_ABI_BREAKING_CHECKS_ENFORCING 1

//#define CATCH_CONFIG_MAIN

//#include "zynq_helper.h" // For Checker

#include "../../lib/seq/builder.h"
#include "../../lib/seq/manager.h"
#include "../../lib/seq/awg/backend.h"
#include "../../lib/utils/streams.h"
using namespace NaCs;
using namespace NaCs::Seq;
using namespace std::literals::string_literals;

static const char *const config_str =
    "tick_per_sec: 1000000000000\n"
    "devices: \n"
    "  AWG1:\n"
    "    backend: awg\n"
    "    config:\n"
    "      url: tcp://localhost:8888\n";

Manager &manager()
{
    static Manager mgr;
    static Call init([&] {
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
    template <typename T>
    void write(T v)
    {
        write((const char*)&v, sizeof(T));
    }
    void write_string(const char *s)
    {
        write(s, strlen(s) + 1);
    }
    std::unique_ptr<Manager::ExpSeq, ExpSeqDeleter> create_sequence()
    {
        auto data = get_buf();
        auto &mgr = manager();
        return {mgr.create_sequence(data.data(), data.size()), ExpSeqDeleter(mgr)};
    }
};

void test_set_pulse() {
    ostream stm;
    // [version <0>: 1B]
    stm.write<uint8_t>(0);
    // [nnodes: 4B]
    // [[[OpCode: 1B][[ArgType: 1B][NodeArg: 1-8B] x narg] /
    //   [OpCode <Interp> : 1B][[ArgType: 1B][NodeArg: 1-8B] x 3][data_id: 4B]] x nnodes]//
    stm.write<uint32_t>(8);
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
    stm.write<double>(12000000000); // 12 ms if used as time
    // 5
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(5000000000);
    // 6
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(15000000000);
    // 7
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(100e6);
    // 8
    stm.write(Builder::OpCode::Identity);
    stm.write(Builder::ArgType::ConstFloat64);
    stm.write<double>(0.7);
    // [nchns: 4B][[chnname: non-empty NUL-terminated string] x nchns]
    stm.write<uint32_t>(4);
    stm.write_string("AWG1/OUT0/CHN1/FREQ");
    stm.write_string("AWG1/OUT0/CHN1/AMP");
    stm.write_string("AWG1/OUT0/CHN2/FREQ");
    stm.write_string("AWG1/OUT0/CHN2/AMP");
    // [ndefvals: 4B][[chnid: 4B][Type: 1B][value: 1-8B] x ndefvals]
    stm.write<uint32_t>(0);
    // [nslots: 4B][[Type: 1B] x nslots]
    stm.write<uint32_t>(0);
    // [nnoramp: 4B][[chnid: 4B] x nnoramp]
    stm.write<uint32_t>(0);
    // [nbasicseqs: 4B][[Basic Sequence] x nbasicseqs]
    stm.write<uint32_t>(1);
    {
        //Basic Seq 1
        // [ntimes: 4B][[sign: 1B][id: 4B][delta_node: 4B][prev_id: 4B] x ntimes]
        stm.write<uint32_t>(4);
        // t1 12 ms after start
        stm.write<uint8_t>((uint8_t)Sign::Pos);
        stm.write<uint32_t>(105);
        stm.write<uint32_t>(4);
        stm.write<uint32_t>(0);
        // t2 5 ms after t1
        stm.write<uint8_t>((uint8_t)Sign::Pos);
        stm.write<uint32_t>(204);
        stm.write<uint32_t>(5);
        stm.write<uint32_t>(1);
        // t3 15 ms after t2
        stm.write<uint8_t>((uint8_t)Sign::Unknown);
        stm.write<uint32_t>(209);
        stm.write<uint32_t>(6);
        stm.write<uint32_t>(2);
        // t4 15 ms after t1
        stm.write<uint8_t>((uint8_t)Sign::Pos);
        stm.write<uint32_t>(209);
        stm.write<uint32_t>(6);
        stm.write<uint32_t>(1);
        // [nendtimes: 4B][[time_id: 4B] x nendtimes]
        stm.write<uint32_t>(2);
        stm.write<uint32_t>(3);
        stm.write<uint32_t>(4);
        // [ntimeorders: 4B][[sign: 1B][id: 4B][before_id: 4B][after_id: 4B] x ntimeorders]
        stm.write<uint32_t>(0);
        // [noutputs: 4B][[id: 4B][time_id: 4B][len: 4B][val: 4B][cond: 4B][chn: 4B] x noutputs]
        stm.write<uint32_t>(5);
        // 1 at 12 ms after start, set chn1 freq to 120 MHz
        stm.write<uint32_t>(8); // id
        stm.write<uint32_t>(1); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(3); // val
        stm.write<uint32_t>(0); // cond // cond of 0 corresponds to default which is on. If it refers to a node, then it can be on or off depending on the evaluation of the node.
        stm.write<uint32_t>(1); // chn
        // 2 at 12 ms after start set chn1 amp to 0.301
        stm.write<uint32_t>(65); // id
        stm.write<uint32_t>(1); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(1); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(2); // chn
        // 3 at 17 ms after start set chn2 freq to 100 MHz
        stm.write<uint32_t>(32); // id
        stm.write<uint32_t>(2); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(7); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(3); // chn
        // 4 at 17 ms after start set chn2 amp to 0.7
        stm.write<uint32_t>(105); // id
        stm.write<uint32_t>(2); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(8); // val
        stm.write<uint32_t>(0); // cond
        stm.write<uint32_t>(4); // chn
        // 5 at 32 ms after start set chn1 amp to 0
        stm.write<uint32_t>(105); // id
        stm.write<uint32_t>(3); // time_id
        stm.write<uint32_t>(0); // len
        stm.write<uint32_t>(2); // val
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

    auto awg_dev = expseq->get_device("AWG1", false);
    auto awg_backend = AWG::Backend::cast(awg_dev);

    expseq->init_run();
    expseq->pre_run();
    expseq->start();
    expseq->wait(60000); // units of ms
}

int main() {

}
