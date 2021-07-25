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

#include "../lib/nacs-seq/builder.h"
#include "../lib/nacs-seq/compiler.h"
#include "../lib/nacs-seq/host_seq.h"
#include "../lib/nacs-seq/manager.h"
#include "../lib/nacs-seq/seq.h"
#include "../lib/nacs-utils/errors.h"
#include "../lib/nacs-utils/llvm/codegen.h"
#include "../lib/nacs-utils/llvm/execute.h"
#include "../lib/nacs-utils/llvm/utils.h"
#include "../lib/nacs-utils/log.h"
#include "../lib/nacs-utils/streams.h"

#include <iostream>
#include <fstream>

using namespace NaCs;

namespace {

enum class Print {
    Builder,
    Seq,
    SeqOpt
};

}

int print(int argc, char **argv, Print print_type)
{
    std::ostream *stm;
    std::unique_ptr<std::ostream> ustm;
    if (argc < 1) {
        Log::error("No input file specified.\n");
        return 1;
    }
    else if (argc < 2) {
        stm = &std::cout;
    }
    else if (argc == 2) {
        ustm = std::make_unique<std::ofstream>(argv[1]);
        if (!ustm->good()) {
            Log::error("Cannot open output file.\n");
            return 1;
        }
        stm = ustm.get();
    }
    else {
        Log::error("Wrong number of arguments.\n");
        return 1;
    }
    std::ifstream istm(argv[0], std::ios::binary);
    if (!istm) {
        Log::error("Cannot open input file.\n");
        return 1;
    }

    std::string str(std::istreambuf_iterator<char>(istm), {});
    auto llvm_ctx = LLVM::new_context();
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.deserialize((const uint8_t*)str.data(), str.size());

    if (print_type == Print::Builder) {
        builder.print(*stm);
        return 0;
    }

    builder.buildseq();

    if (print_type == Print::SeqOpt) {
        seq.prepare();
        seq.optimize();
    }

    seq.print(*stm);

    return 0;
}

int compile(int argc, char **argv)
{
    if (argc < 1) {
        Log::error("No input file specified.\n");
        return 1;
    }
    std::ifstream istm(argv[0], std::ios::binary);
    if (!istm) {
        Log::error("Cannot open input file.\n");
        return 1;
    }

    std::string str(std::istreambuf_iterator<char>(istm), {});
    auto llvm_ctx = LLVM::new_context();
    Seq::Seq seq(*llvm_ctx);
    Seq::Builder builder(seq);

    builder.deserialize((const uint8_t*)str.data(), str.size());
    builder.buildseq();
    seq.prepare();
    seq.optimize();

    Seq::HostSeq host_seq;
    LLVM::Exe::Engine engine;
    Seq::Compiler compiler(host_seq, seq, engine,
                           seq.env().cg_context()->get_extern_resolver());
    auto obj_id = compiler.compile();
    engine.free(obj_id);

    printf("Compilation completed successfully.\n");

    return 0;
}

int create_sequence(int argc, char **argv)
{
    if (argc < 1) {
        Log::error("No input file specified.\n");
        return 1;
    }
    if (argc < 2) {
        Log::error("No config file specified.\n");
        return 1;
    }
    std::ifstream istm(argv[0], std::ios::binary);
    if (!istm) {
        Log::error("Cannot open input file.\n");
        return 1;
    }

    Seq::Manager manager;
    manager.load_config_file(argv[1]);

    std::string str(std::istreambuf_iterator<char>(istm), {});
    auto seq = manager.create_sequence((const uint8_t*)str.data(), str.size());

    if (!seq) {
        printf("Sequence creation failed.\n");
        return 1;
    }

    printf("Sequence created successfully.\n");
    manager.free_sequence(seq);
    return 0;
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        Log::error("No action specified.\n");
        return 1;
    }
    if (strcmp(argv[1], "print_builder") == 0) {
        return print(argc - 2, argv + 2, Print::Builder);
    }
    else if (strcmp(argv[1], "print_seq") == 0) {
        return print(argc - 2, argv + 2, Print::Seq);
    }
    else if (strcmp(argv[1], "print_seq_opt") == 0) {
        return print(argc - 2, argv + 2, Print::SeqOpt);
    }
    else if (strcmp(argv[1], "compile") == 0) {
        return compile(argc - 2, argv + 2);
    }
    else if (strcmp(argv[1], "create_sequence") == 0) {
        return create_sequence(argc - 2, argv + 2);
    }
    else {
        Log::error("Unknown action: %s.\n", argv[1]);
        return 1;
    }
}
