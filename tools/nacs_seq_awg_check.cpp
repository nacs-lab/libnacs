//

// Offline test driver for the AWG RF power check (lib/nacs-seq/awg_check.cpp).
//
//   nacs-seq-awg-check <seq.bin> <config.yml> [globals.txt] [shots] [trace.json]
//
// seq.bin      the bytes MATLAB passes to create_sequence (ExpSeq.serialize)
// config.yml   MUST set `use_dummy_device: true` (refused otherwise: no hardware)
// globals.txt  lines "<global id> <value>" applied before every shot ("-" = none)
// shots        shots to run (default 20), for the timing statistics
// trace.json   optional: the AWG channels of every bseq of the first shot, from
//              get_nominal_output (the .seq dump path), for an independent recount
//
// Prints one JSON object: the AWG outputs and tones found, the first shot's per-bseq
// report, the first refusal (if any), and the time spent in the check and in pre_run.
// Runs with the check disabled too (`awg_power: {enabled: false}`): then only the
// pre_run time is reported, the baseline for the cost.

#define LLVM_DISABLE_ABI_BREAKING_CHECKS_ENFORCING 1

#include "../lib/nacs-utils/log.h"
#include "../lib/nacs-seq/manager.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace NaCs;

static std::string jstr(const std::string &s)
{
    std::string o = "\"";
    for (char c: s) {
        if (c == '"' || c == '\\') {
            o += '\\';
            o += c;
        }
        else if ((unsigned char)c < 0x20) {
            o += ' ';
        }
        else {
            o += c;
        }
    }
    return o + "\"";
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        Log::error("usage: nacs-seq-awg-check <seq.bin> <config.yml> [globals.txt] [shots]\n");
        return 1;
    }
    Seq::Manager mgr;
    mgr.load_config_file(argv[2]);
    if (!mgr.uses_dummy_device()) {
        Log::error("refused: the config must set use_dummy_device: true\n");
        return 2;
    }
    std::vector<std::pair<uint32_t,double>> globals;
    if (argc > 3 && std::string(argv[3]) != "-") {
        std::ifstream g(argv[3]);
        uint32_t id;
        double val;
        while (g >> id >> val) {
            globals.emplace_back(id, val);
        }
    }
    int shots = argc > 4 ? atoi(argv[4]) : 20;
    if (shots < 1)
        shots = 1;
    std::string trace_path = argc > 5 ? argv[5] : "";
    std::ostringstream trace;
    trace.precision(17);
    trace << "[";
    bool trace_first = true;

    std::ifstream istm(argv[1], std::ios::binary);
    if (!istm) {
        Log::error("Cannot open input file.\n");
        return 1;
    }
    std::string str(std::istreambuf_iterator<char>(istm), {});
    auto t0 = std::chrono::steady_clock::now();
    auto seq = mgr.create_sequence((const uint8_t*)str.data(), str.size());
    double compile_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - t0).count();
    bool enabled = mgr.awg_limits().enabled;

    std::string refusal;
    int refused_shot = -1;
    int refused_bseq = -1;
    // Per bseq index: check / pre_run times over all shots.
    std::vector<std::vector<double>> check_us, pre_run_us;
    Seq::Manager::ExpSeq::AwgReport first;
    for (int k = 0; k < shots && refused_shot < 0; k++) {
        for (auto &[id, val]: globals)
            seq->host_seq.set_global(id, val);
        seq->init_run();
        while (true) {
            auto idx = seq->host_seq.cur_seq_idx();
            if (check_us.size() <= idx) {
                check_us.resize(idx + 1);
                pre_run_us.resize(idx + 1);
            }
            auto p0 = std::chrono::steady_clock::now();
            try {
                seq->pre_run();
            }
            catch (const Seq::Error &err) {
                refusal = err.what();
                refused_shot = k;
                refused_bseq = (int)idx;
                break;
            }
            pre_run_us[idx].push_back(std::chrono::duration<double, std::micro>(
                                          std::chrono::steady_clock::now() - p0).count());
            if (enabled && !seq->awg_report.bseqs.empty())
                check_us[idx].push_back(seq->awg_report.bseqs.back().check_us);
            if (k == 0 && !trace_path.empty()) {
                size_t nout = 0;
                auto outs = seq->get_nominal_output(mgr.awg_limits().pts_per_ramp, &nout);
                trace << (trace_first ? "" : ",") << "{\"bseq_idx\":" << idx
                      << ",\"tick_per_sec\":" << mgr.tick_per_sec() << ",\"chns\":{";
                trace_first = false;
                bool cf = true;
                for (size_t c = 0; c < nout; c++) {
                    std::string name(outs[c].name, outs[c].name_sz);
                    while (!name.empty() && name.back() == '\0')
                        name.pop_back();
                    if (name.compare(0, 4, "AWG1") != 0 && name.compare(0, 4, "AWG2") != 0)
                        continue;
                    trace << (cf ? "" : ",") << jstr(name) << ":[";
                    cf = false;
                    for (size_t j = 0; j < outs[c].npts; j++)
                        trace << (j ? "," : "") << "[" << outs[c].times[j] << ","
                              << outs[c].values[j] << "]";
                    trace << "]";
                }
                trace << "}}";
            }
            seq->start();
            while (!seq->wait(100)) {
            }
            if (seq->post_run() == 0) {
                break;
            }
        }
        if (k == 0)
            first = seq->awg_report;
    }
    if (refused_shot == 0)
        first = seq->awg_report;

    auto stat = [] (std::vector<double> v, double q) {
        if (v.empty())
            return 0.0;
        std::sort(v.begin(), v.end());
        return v[std::min(v.size() - 1, size_t(q * (v.size() - 1) + 0.5))];
    };
    std::ostringstream o;
    o.precision(10);
    o << "{\"enabled\":" << (enabled ? "true" : "false") << ",\"compile_ms\":" << compile_ms
      << ",\"verdict\":" << jstr(refused_shot >= 0 ? "FAIL" :
                                 (first.warns.empty() ? "PASS" : "WARN"))
      << ",\"refused_shot\":" << refused_shot << ",\"refused_bseq_idx\":" << refused_bseq
      << ",\"refusal\":" << jstr(refusal) << ",\"shots\":" << shots << ",\"outs\":[";
    for (size_t i = 0; i < seq->awg_outs.size(); i++) {
        auto &out = seq->awg_outs[i];
        int namp = 0, nfreq = 0;
        for (auto &t: out.tones) {
            namp += t.amp_chn != 0;
            nfreq += t.freq_chn != 0;
        }
        o << (i ? "," : "") << "{\"name\":" << jstr(out.name) << ",\"tones\":"
          << out.tones.size() << ",\"amp_chns\":" << namp << ",\"freq_chns\":" << nfreq << "}";
    }
    o << "],\"warns\":[";
    for (size_t i = 0; i < first.warns.size(); i++)
        o << (i ? "," : "") << jstr(first.warns[i]);
    o << "],\"bseqs\":[";
    for (size_t i = 0; i < first.bseqs.size(); i++) {
        auto &b = first.bseqs[i];
        o << (i ? "," : "") << "{\"id\":" << b.bseq_id << ",\"outs\":[";
        for (size_t j = 0; j < b.outs.size(); j++) {
            auto &r = b.outs[j];
            o << (j ? "," : "") << "{\"name\":" << jstr(r.name) << ",\"bound\":" << r.bound
              << ",\"exact\":" << (r.exact ? "true" : "false") << ",\"sum_amp2\":"
              << r.sum_amp2 << ",\"t_ms\":" << r.t_ms << "}";
        }
        o << "]}";
    }
    o << "],\"timing\":[";
    for (size_t i = 0; i < pre_run_us.size(); i++) {
        o << (i ? "," : "") << "{\"bseq_idx\":" << i << ",\"n\":" << pre_run_us[i].size()
          << ",\"pre_run_us_median\":" << stat(pre_run_us[i], 0.5)
          << ",\"check_us_median\":" << stat(check_us[i], 0.5)
          << ",\"check_us_p99\":" << stat(check_us[i], 0.99)
          << ",\"check_us_max\":" << stat(check_us[i], 1.0)
          << ",\"check_us_first\":" << (check_us[i].empty() ? 0.0 : check_us[i][0])
          << ",\"check_us_rest_max\":"
          << (check_us[i].size() > 1 ? *std::max_element(check_us[i].begin() + 1,
                                                          check_us[i].end()) : 0.0)
          << "}";
    }
    o << "]}";
    std::cout << o.str() << std::endl;
    if (!trace_path.empty()) {
        trace << "]";
        std::ofstream(trace_path) << trace.str();
    }
    mgr.free_sequence(seq);
    return 0;
}
