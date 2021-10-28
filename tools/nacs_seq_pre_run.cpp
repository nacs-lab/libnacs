//

#define LLVM_DISABLE_ABI_BREAKING_CHECKS_ENFORCING 1

#include "../lib/nacs-utils/log.h"
#include "../lib/nacs-seq/manager.h"

#include <iostream>
#include <fstream>

using namespace NaCs;

int main(int argc, char **argv)
{
    if (argc < 2) {
        Log::error("No input file specified.\n");
        return 1;
    }
    if (argc < 3) {
        Log::error("No config file specified.\n");
        return 1;
    }

    Seq::Manager mgr;
    mgr.load_config_file(argv[2]);

    std::ifstream istm(argv[1], std::ios::binary);
    if (!istm) {
        Log::error("Cannot open input file.\n");
        return 1;
    }

    std::string str(std::istreambuf_iterator<char>(istm), {});
    auto seq = mgr.create_sequence((const uint8_t*)str.data(), str.size());

    seq->init_run();
    seq->pre_run();
    return 0;
}
