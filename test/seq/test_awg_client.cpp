//

#define LLVM_DISABLE_ABI_BREAKING_CHECKS_ENFORCING 1

//#include "DummyClient.h"
#include <llvm/IR/LLVMContext.h>
#include "../../lib/seq/awg/DummyClient.h"

#include <iostream>
#include <stdio.h>

using namespace NaCs;
int main ()
{
    std::string fname = "/etc/client_config.yml";
    std::cout << "here" << std::endl;
    Seq::AWG::Dummy::DummyClient client = Seq::AWG::Dummy::DummyClient();
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
        printf("%X", client.id_bc[i]);
    }
    std::cout << std::endl;
}
