//

#define CATCH_CONFIG_MAIN

#include <catch2/catch.hpp>

#define LLVM_DISABLE_ABI_BREAKING_CHECKS_ENFORCING 1

//#include "DummyClient.h"
#include <llvm/IR/LLVMContext.h>
#include "../../lib/nacs-seq/awg/DummyClient.h"

#include <iostream>
#include <stdio.h>

using namespace NaCs;

TEST_CASE("client") {
    std::string fname = "/etc/client_config.yml";
    std::cout << "here" << std::endl;
    Seq::AWG::Dummy::DummyClient client;
#if 0
    client.loadConfig(fname.data());
    std::cout << client.m_url << std::endl;
    // test network
    client.init_run();
    std::cout << "server id: " << client.m_server_id << std::endl;
    std::cout << "client id: " << client.m_client_id << std::endl;
    std::cout << "triple str: " << client.m_triple_str << std::endl;
    std::cout << "cpu str: " << client.m_cpu_str << std::endl;
    std::cout << "feature str " << client.m_feature_str << std::endl;
    //std::cout << client.id_bc << std::endl;

    //std::cout << id_bc: << std::endl;
    for (int i = 0; i < client.id_bc.size() ; i++) {
        printf(" %X", client.id_bc[i]);
    }
    std::cout << std::endl;
#endif
    // channel testing
    client.add_channel(15, "OUT1/CHN2/FREQ");
    client.add_channel(10, "OUT2/CHN5/AMP");
    client.add_channel(12, "OUT1/CHN1/PHASE");
    client.sort_channels();
    // Freq 0 , amp 1, phase 2
    Seq::AWG::Dummy::DummyClient::ChnType chn_type;
    uint8_t phys_id;
    uint32_t chn_id;
    size_t nchns = client.m_linear_chns.size();
    REQUIRE(nchns == 3);
    client.get_chn_from_lin(0, chn_type, phys_id, chn_id);
    REQUIRE(chn_type == Seq::AWG::Dummy::DummyClient::ChnType::Phase);
    REQUIRE(phys_id == 1);
    REQUIRE(chn_id == 1);
    client.get_chn_from_lin(1, chn_type, phys_id, chn_id);
    REQUIRE(chn_type == Seq::AWG::Dummy::DummyClient::ChnType::Freq);
    REQUIRE(phys_id == 1);
    REQUIRE(chn_id == 2);
    client.get_chn_from_lin(2, chn_type, phys_id, chn_id);
    REQUIRE(chn_type == Seq::AWG::Dummy::DummyClient::ChnType::Amp);
    REQUIRE(phys_id == 2);
    REQUIRE(chn_id == 5);
}
