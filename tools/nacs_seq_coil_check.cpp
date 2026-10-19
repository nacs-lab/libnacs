//

// Offline test driver for the Feshbach coil check (lib/nacs-seq/coil_check.cpp).
//
//   nacs-seq-coil-check <seq.bin> <config.yml> [globals.txt] [shots] [period_s]
//
// seq.bin      the bytes MATLAB passes to create_sequence (ExpSeq.serialize)
// config.yml   MUST set `use_dummy_device: true` (refused otherwise: no hardware)
// globals.txt  lines "<global id> <value>" applied before every shot ("-" = none)
// shots        simulated back-to-back shots (default: enough to fill the duty window)
// period_s     shot period for the simulated clock; "auto" (default) = the compiled
//              shot length + 0.10 s per bseq, the coil_guard.py estimate
//
// Prints one JSON object: the first shot's per-bseq report, warnings, the shot and
// bseq of the first refusal (if any), and the time spent in the check.

#define LLVM_DISABLE_ABI_BREAKING_CHECKS_ENFORCING 1

#include "../lib/nacs-utils/log.h"
#include "../lib/nacs-seq/manager.h"

#include <algorithm>
#include <chrono>
#include <cmath>
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

static std::string jintervals(const std::vector<Seq::Manager::ExpSeq::CoilInterval> &v)
{
    std::ostringstream o;
    o.precision(17);
    o << "[";
    for (size_t i = 0; i < v.size(); i++)
        o << (i ? "," : "") << "[" << v[i].t0_ms << "," << v[i].t1_ms << "]";
    o << "]";
    return o.str();
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        Log::error("usage: nacs-seq-coil-check <seq.bin> <config.yml> [globals.txt] [shots] [period_s]\n");
        return 1;
    }
    Seq::Manager mgr;
    mgr.load_config_file(argv[2]);
    if (!mgr.uses_dummy_device()) {
        Log::error("refused: the config must set use_dummy_device: true\n");
        return 2;
    }
    if (!mgr.coil_limits().enabled) {
        Log::error("refused: the config has no enabled feshbach_coil block\n");
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
    int shots = argc > 4 ? atoi(argv[4]) : 0;
    std::string period_arg = argc > 5 ? argv[5] : "auto";

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
    double s_per_tick = 1.0 / double(mgr.tick_per_sec());

    // One shot on the virtual clock starting at t_shot. Returns false on a refusal.
    std::string refusal;
    int refused_bseq = -1;
    std::vector<double> check_us, pre_run_us;
    auto run_shot = [&] (double t_shot, double &shot_s, int &nbseq) {
        for (auto &[id, val]: globals)
            seq->host_seq.set_global(id, val);
        seq->init_run();
        double t = t_shot;
        shot_s = 0;
        nbseq = 0;
        while (true) {
            mgr.coil_clock_override = t;
            auto p0 = std::chrono::steady_clock::now();
            try {
                seq->pre_run();
            }
            catch (const Seq::Error &err) {
                refusal = err.what();
                refused_bseq = (int)seq->host_seq.cur_seq_idx();
                return false;
            }
            pre_run_us.push_back(std::chrono::duration<double, std::micro>(
                                     std::chrono::steady_clock::now() - p0).count());
            check_us.push_back(seq->coil_report.bseqs.back().check_us);
            seq->start();
            while (!seq->wait(100)) {
            }
            double len_s = double(seq->host_seq.seqs[seq->host_seq.cur_seq_idx()].length) *
                s_per_tick;
            shot_s += len_s;
            nbseq++;
            t += len_s + 0.10;
            if (seq->post_run() == 0) {
                break;
            }
        }
        return true;
    };

    // First shot: its report, and the coil_guard period estimate.
    double shot_s = 0;
    int nbseq = 0;
    mgr.coil_clear_history();
    bool ok = run_shot(0, shot_s, nbseq);
    auto first = seq->coil_report;
    double period = shot_s + 0.10 * nbseq;
    if (period_arg != "auto")
        period = atof(period_arg.c_str());
    double W = mgr.coil_limits().duty_window_s;
    if (shots <= 0)
        shots = int(std::ceil(W / std::max(period, 1e-3))) + 5;
    int refused_shot = ok ? -1 : 0;
    double max_duty = first.duty;
    std::vector<std::string> warns = first.warns;
    if (ok) {
        // The whole stream as a new scan: empty history, no pending refusal.
        mgr.coil_clear_history();
        mgr.coil_new_run();
        for (int k = 0; k < shots; k++) {
            double s;
            int n;
            if (!run_shot(k * period, s, n)) {
                refused_shot = k;
                break;
            }
            max_duty = std::max(max_duty, seq->coil_report.duty);
            for (auto &w: seq->coil_report.warns) {
                if (std::find(warns.begin(), warns.end(), w) == warns.end()) {
                    warns.push_back(w);
                }
            }
        }
    }

    auto stat = [] (std::vector<double> v, double q) {
        if (v.empty())
            return 0.0;
        std::sort(v.begin(), v.end());
        return v[std::min(v.size() - 1, size_t(q * (v.size() - 1) + 0.5))];
    };
    std::ostringstream o;
    o.precision(10);
    o << "{\"coil_check_build\":\"coil-check\",\"compile_ms\":" << compile_ms
      << ",\"verdict\":" << jstr(refused_shot >= 0 ? "FAIL" : (warns.empty() ? "PASS" : "WARN"))
      << ",\"refused_shot\":" << refused_shot << ",\"refused_bseq_idx\":" << refused_bseq
      << ",\"refusal\":" << jstr(refusal)
      << ",\"shots\":" << shots << ",\"period_s\":" << period
      << ",\"shot_on_ms\":" << first.shot_on_ms << ",\"peak_G\":" << first.peak_G
      << ",\"max_duty\":" << max_duty << ",\"pred_duty\":" << first.pred_duty
      << ",\"pred_period_s\":" << first.pred_period_s << ",\"warns\":[";
    for (size_t i = 0; i < warns.size(); i++)
        o << (i ? "," : "") << jstr(warns[i]);
    o << "],\"bseqs\":[";
    for (size_t i = 0; i < first.bseqs.size(); i++) {
        auto &b = first.bseqs[i];
        o << (i ? "," : "") << "{\"id\":" << b.bseq_id << ",\"seq_ms\":" << b.seq_ms
          << ",\"on_ms\":" << b.on_ms << ",\"peak_G\":" << b.peak_G << ",\"ilk_ms\":" << b.ilk_ms
          << ",\"windows\":" << jintervals(b.windows) << ",\"gaps\":" << jintervals(b.gaps)
          << ",\"exempt\":" << jintervals(b.exempt) << "}";
    }
    o << "],\"check_us\":{\"n\":" << check_us.size() << ",\"median\":" << stat(check_us, 0.5)
      << ",\"p99\":" << stat(check_us, 0.99) << ",\"max\":" << stat(check_us, 1.0)
      << "},\"pre_run_us\":{\"median\":" << stat(pre_run_us, 0.5)
      << ",\"max\":" << stat(pre_run_us, 1.0) << "}}";
    std::cout << o.str() << std::endl;
    mgr.free_sequence(seq);
    return 0;
}
