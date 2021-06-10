//

#define LLVM_DISABLE_ABI_BREAKING_CHECKS_ENFORCING 1

//#include "DummyClient.h"
#include <llvm/IR/LLVMContext.h>
#include "../../lib/seq/awg/DummyClient.h"

#include <iostream>

using namespace NaCs;
int main ()
{
    std::string fname = "/etc/client_config.yml";
    std::cout << "here" << std::endl;
    Seq::AWG::Dummy::DummyClient client = Seq::AWG::Dummy::DummyClient();
    //client.loadConfig(fname.data());
    //std::cout << client.m_url << std::endl;
    // test network
    
}
