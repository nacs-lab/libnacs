// This file should send sequences that are pretty real
// frequency units are in Hz
// amplitude is scale from 0 to 1
// phase is scale from 0 to 1
// time is in units of ps, required from backend

#define LLVM_DISABLE_ABI_BREAKING_CHECKS_ENFORCING 1

//#include "DummyClient.h"
#include <llvm/IR/LLVMContext.h>
#include "../../lib/nacs-seq/awg/DummyClient.h"

#include <llvm/ADT/StringRef.h>
#include <llvm/IR/LegacyPassManager.h>
#include <llvm/IR/Verifier.h>
#include <llvm/Transforms/IPO.h>


#include <iostream>
#include <stdio.h>

using namespace NaCs;
int main ()
{
    //std::string fname = "C:/msys64/etc/client_config.yml"; for windows
    std::string fname = "/etc/client_config.yml";
    Seq::AWG::Dummy::DummyClient client = Seq::AWG::Dummy::DummyClient();
    client.loadConfig(fname.data());
    // test network
    client.init_run();
    std::cout << "server id: " << client.m_server_id << std::endl;
    std::cout << "client id: " << client.m_client_id << std::endl;
    std::cout << "triple str: " << client.m_triple_str << std::endl;
    std::cout << "cpu str: " << client.m_cpu_str << std::endl;
    std::cout << "feature str " << client.m_feature_str << std::endl;
    
    // channel testing
    client.add_channel(1, "OUT0/CHN1/FREQ");
    client.add_channel(2, "OUT0/CHN1/AMP");
    client.add_channel(3, "OUT0/CHN1/PHASE");
    client.add_channel(4, "OUT0/CHN2/FREQ");
    client.add_channel(5, "OUT0/CHN2/AMP");
    client.add_channel(6, "OUT0/CHN2/PHASE");
    /* client.add_channel(7, "OUT0/CHN1/FREQ"); */
    /* client.add_channel(8, "OUT0/CHN1/AMP"); */
    /* client.add_channel(9, "OUT0/CHN1/PHASE"); */
    client.sort_channels();
    Seq::HostSeq::Value v0, v1, v2, v3, v4, v5, v6, v7, v8;
    v0.i64 = 1e12; // 1 second
    v1.i64 = 6e12; // 6 seconds
    v2.b = true;
    v3.f64 = 0.5;
    v4.f64 = 80e6;
    v5.b = false;
    v6.f64 = 0.1;
    v7.i64 = 4e12; // 4 seconds
    v8.f64 = 85e6;
    // Freq 0 , amp 1, phase 2
    client.addConstVal(0, Seq::HostSeq::Type::Int64, v0); // 0 // 51.2 ns is 1 cycle on awg. This is in units of ps supposedly.
    client.addConstVal(0, Seq::HostSeq::Type::Int64, v1); // 1
    client.addConstVal(0, Seq::HostSeq::Type::Bool, v2); // 2
    client.addConstVal(0, Seq::HostSeq::Type::Float64, v3); // 3
    client.addNConstVal(0, Seq::HostSeq::Type::Float64, v4); // 4
    client.addNConstVal(0, Seq::HostSeq::Type::Bool, v5); // 5
    client.addNConstVal(0, Seq::HostSeq::Type::Float64, v6); // 6
    client.addNConstVal(0, Seq::HostSeq::Type::Int64, v7); // 7
    client.addNConstVal(0, Seq::HostSeq::Type::Float64, v8); // 8
    // pulse type: 1 amp set, 4 freq set,8 phase
    client.addPulse(0, 1, 1, 0, uint32_t(-1), 4, 2, "");
    client.addPulse(0, 2, 2, 0, uint32_t(-1), 3, 2, ""); // set amplitude and frequency of chn1 to 0.5 and 80e6 at 1 second
    client.addPulse(0, 4, 3, 1, uint32_t(-1), 8, 2, ""); 
    client.addPulse(0, 5, 4, 1, uint32_t(-1), 6, 2, ""); // set amplitude and frequency of chn2 to 0.1 and 85e6 at 6 seconds;
    client.addPulse(0, 2, 5, 7, uint32_t(-1), 6, 5, ""); // set amplitude of chn1 to 0.1, but should be off!

    /* // SECOND BASIC SEQ */
    /* Seq::HostSeq::Value w0, w1, w2, w3, w4, w5, w6; */
    /* w0.i64 = 51201 * 2; */
    /* w1.b = true; */
    /* w2.b = false; */
    /* w3.f64 = 12.1; */
    /* w4.f64 = 21.1; */
    /* w5.f64 = 80; */
    /* w6.i64 = 51201 * 3; */
    /* // Freq 0 , amp 1, phase 2 */
    /* client.addConstVal(1, Seq::HostSeq::Type::Int64, w0); // 0 // 51.2 ns is 1 cycle on awg. This is in units of ps supposedly. */
    /* client.addConstVal(1, Seq::HostSeq::Type::Bool, w1); // 1 */
    /* client.addNConstVal(1, Seq::HostSeq::Type::Bool, w2); // 2 */
    /* client.addNConstVal(1, Seq::HostSeq::Type::Float64, w3); // 3 */
    /* client.addNConstVal(1, Seq::HostSeq::Type::Float64, w4); // 4 */
    /* client.addNConstVal(1, Seq::HostSeq::Type::Float64, w5); // 5 */
    /* client.addNConstVal(1, Seq::HostSeq::Type::Int64, w6); // 6 */

    /* // pulse type: 1 amp set, 4 freq set,8 phase */
    /* client.addPulse(1, 7, 1, 0, uint32_t(-1), 3, 1, ""); // set frequency of out0/chn1 to 12.1 at time 2 */
    /* client.addPulse(1, 8, 2, 0, uint32_t(-1), 4, 1, ""); // set amplitude of out0/chn1 to 21.1 at time 2 */
    /* client.addPulse(1, 4, 3, 6, uint32_t(-1), 5, 1, ""); // set frequency of out1/chn2 to 80 at time 2 */
    /* client.addPulse(1, 5, 4, 0, uint32_t(-1), 3, 1, ""); // set amplitude of out1/chn2 to 12.1 at time 2 */
    /* client.addPulse(1, 6, 5, 6, uint32_t(-1), 5, 2, ""); // nothing happens cause pulse is off */

    
    client.prepare(); // compiles
    /* auto bseq = client.m_seqs[0]; */
    // for (int i = 0; i < bseq.pulses_bc.size() ; i++) {
    //     printf(" %X", bseq.pulses_bc[i]);
    // }
    // std::cout << std::endl;
    // for (int i = 0; i < bseq.types.size(); i++) {
    //     printf(" %X", bseq.types[i]);
    // }
    // std::cout << std::endl;
    // for (int i = 0; i < bseq.types_bc.size() ; i++) {
    //     printf(" %X", bseq.types_bc[i]);
    // }
    // std::cout << std::endl;
    client.pre_run();
    client.wait();

    // second time around. change some values
    v4.f64 = 83e6;
    v5.b = true;
    v6.f64 = 0.1;
    v7.i64 = 4e12; // 4 seconds
    v8.f64 = 81e6;
    client.setVals(0, 0, v4);
    client.setVals(0, 1, v5);
    client.setVals(0, 2, v6);
    client.setVals(0, 3, v7);
    client.setVals(0, 4, v8);
    client.pre_run();
    client.wait();

    // client.set_cur_seq_id(1);
    // client.pre_run();
    // client.wait();
}