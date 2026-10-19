/*************************************************************************
 *   Copyright (c) 2026 NaCs lab                                         *
 *                                                                       *
 *   This library is free software; you can redistribute it and/or       *
 *   modify it under the terms of the GNU Lesser General Public          *
 *   License as published by the Free Software Foundation; either        *
 *   version 3.0 of the License, or (at your option) any later version.  *
 *************************************************************************/

// Feshbach coil check in `pre_run` (enabled by the `feshbach_coil` config block).
//
// For the basic sequence about to run, rebuild the FeshbachCurr and interlockTTL
// traces exactly as `get_nominal_output` does (the trace the .seq dumps carry and
// web/verify/coil_guard.py checks), then apply the coil_guard rules:
//   energized   field above thr_G (B = b_per_v * V + b_offset), summed over the shot
//   per shot    <= max_on_ms
//   field       <= max_B_G on every sample
//   TTL26       high over every energized interval, except a window that starts from
//               and returns to the park voltage (fbCondSeq conditioning ramp)
//   duty        <= max_duty of any duty_window_s, two ways: predicted at the end of
//               each shot for the scan at full speed (coil_predict; the next shot is
//               refused), and measured from the energized intervals of the bseqs that
//               actually started (manager history) plus this one (this bseq refused)
// A hard limit throws Error(BasicSeq, CoilLimit) before any device sees the bseq.
// Warn levels are recorded in `coil_report` only.

#include "manager.h"
#include "error.h"

#include "../nacs-utils/number.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>

namespace NaCs::Seq {

namespace {

static constexpr double time_tol_ms = 1e-6; // coil_guard TIME_TOL_MS

struct Trace {
    std::vector<double> t; // ms
    std::vector<double> v;
};

// Same walk as Manager::ExpSeq::get_nominal_output, for one channel.
static Trace channel_trace(const HostSeq &host_seq, uint32_t seq_idx, uint32_t chn,
                           uint64_t pts_per_ramp, double ms_per_tick)
{
    Trace out;
    if (chn == 0)
        return out;
    struct P {
        uint32_t id;
        int64_t time;
        int64_t len;
        double (*ramp_func)(double, void*);
        double endvalue;
        bool cond;
    };
    std::vector<P> pulses;
    auto &bseq = host_seq.seqs[seq_idx];
    for (auto &pulse: bseq.pulses) {
        if (pulse.is_measure() || pulse.chn != chn)
            continue;
        int64_t time = host_seq.get_time(pulse.time);
        int64_t len = pulse.len == uint32_t(-1) ? 0 :
            round<int64_t>(host_seq.get_value(pulse.len, seq_idx));
        double endvalue = host_seq.get_value(pulse.endvalue, seq_idx);
        bool cond = pulse.cond == uint32_t(-1) ? true :
            (bool)host_seq.get_value(pulse.cond, seq_idx);
        pulses.push_back({pulse.id, time, len, pulse.ramp_func, endvalue, cond});
    }
    std::sort(pulses.begin(), pulses.end(), [] (auto &p1, auto &p2) {
        if (p1.time != p2.time)
            return p1.time < p2.time;
        return p1.id < p2.id;
    });
    std::vector<int64_t> times{0};
    std::vector<double> values{host_seq.start_values[chn - 1].f64};
    auto npulses = pulses.size();
    for (size_t i = 0; i < npulses;) {
        auto &pulse = pulses[i];
        if (!pulse.cond) {
            i++;
            continue;
        }
        const P *next = nullptr;
        i++;
        while (i < npulses) {
            if (pulses[i].cond) {
                next = &pulses[i];
                break;
            }
            i++;
        }
        if (pulse.len == 0) {
            if (next && next->time <= pulse.time)
                continue;
            if (times.back() != pulse.time) {
                times.push_back(pulse.time);
                values.push_back(values.back());
            }
            times.push_back(pulse.time);
            values.push_back(pulse.endvalue);
        }
        else {
            int64_t end_time = pulse.time + pulse.len;
            if (next)
                end_time = std::min(next->time, end_time);
            int64_t ramp_time = end_time - pulse.time;
            if (ramp_time <= 0)
                continue;
            if (times.back() != pulse.time) {
                times.push_back(pulse.time);
                values.push_back(values.back());
            }
            uint64_t npts = std::max<uint64_t>(2, std::min((uint64_t)ramp_time, pts_per_ramp));
            for (uint64_t j = 0; j < npts; j++) {
                int64_t this_pt = round<int64_t>(j * ramp_time / (npts - 1));
                times.push_back(pulse.time + this_pt);
                values.push_back(pulse.ramp_func(this_pt, (void*)host_seq.values.data()));
            }
        }
    }
    out.t.reserve(times.size());
    for (auto t: times)
        out.t.push_back(double(t) * ms_per_tick);
    out.v = std::move(values);
    return out;
}

using Intervals = std::vector<Manager::ExpSeq::CoilInterval>;

// coil_guard._intervals: held-value pieces where pred(value) holds.
template<typename Pred>
static Intervals intervals(const Trace &tr, double t_end, Pred &&pred)
{
    Intervals out;
    auto n = tr.t.size();
    if (n == 0)
        return out;
    bool open = false;
    double a = 0, b = 0;
    for (size_t i = 0; i < n; i++) {
        double next_t = i + 1 < n ? tr.t[i + 1] : std::max(t_end, tr.t[n - 1]);
        if (pred(tr.v[i])) {
            if (!open) {
                open = true;
                a = tr.t[i];
            }
            b = next_t;
        }
        else if (open) {
            out.push_back({a, b});
            open = false;
        }
    }
    if (open)
        out.push_back({a, b});
    Intervals res;
    for (auto &w: out) {
        if (w.t1_ms - w.t0_ms > 0) {
            res.push_back(w);
        }
    }
    return res;
}

// coil_guard._parked_window
static bool parked_window(const Trace &tr, double a, double b, double park_v, double margin)
{
    bool before = false;
    bool after = false;
    for (size_t i = 0; i < tr.t.size(); i++) {
        if (std::abs(tr.v[i] - park_v) >= 1e-9)
            continue;
        auto t = tr.t[i];
        if (t <= a + time_tol_ms && t >= a - margin)
            before = true;
        if (t >= b - time_tol_ms && t <= b + margin)
            after = true;
    }
    return before && after;
}

// coil_guard._uncovered
static Intervals uncovered(const Intervals &on, Intervals high)
{
    std::sort(high.begin(), high.end(), [] (auto &x, auto &y) {
        if (x.t0_ms != y.t0_ms)
            return x.t0_ms < y.t0_ms;
        return x.t1_ms < y.t1_ms;
    });
    Intervals gaps;
    for (auto &w: on) {
        double x = w.t0_ms;
        for (auto &h: high) {
            if (h.t1_ms <= x + time_tol_ms)
                continue;
            if (h.t0_ms > x + time_tol_ms)
                gaps.push_back({x, std::min(h.t0_ms, w.t1_ms)});
            x = std::max(x, h.t1_ms);
            if (x >= w.t1_ms - time_tol_ms) {
                break;
            }
        }
        if (x < w.t1_ms - time_tol_ms)
            gaps.push_back({x, w.t1_ms});
    }
    Intervals res;
    for (auto &g: gaps) {
        if (g.t1_ms - g.t0_ms > time_tol_ms) {
            res.push_back(g);
        }
    }
    return res;
}

static std::string fmt(const char *f, double a, double b)
{
    char buf[256];
    snprintf(buf, sizeof(buf), f, a, b);
    return buf;
}

}

NACS_EXPORT() double Manager::coil_clock() const
{
    if (coil_clock_override >= 0)
        return coil_clock_override;
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}

void Manager::ExpSeq::coil_check()
{
    auto t_start = std::chrono::steady_clock::now();
    auto &L = m_mgr.m_coil_limits;
    auto seq_idx = host_seq.cur_seq_idx();
    auto &bseq = host_seq.seqs[seq_idx];
    double ms_per_tick = 1e3 / double(m_mgr.tick_per_sec());
    if (!m_mgr.m_coil_predicted_fail.empty()) {
        // An earlier shot of this scan would exceed the duty limit once the scan runs
        // at full speed: refuse before the next shot is sent.
        std::string msg = "Feshbach coil check refused this sequence: " +
            m_mgr.m_coil_predicted_fail;
        m_mgr.m_coil_predicted_fail.clear();
        coil_report.fails.push_back(msg);
        throw Error(Error::Type::BasicSeq, Error::BasicSeq::CoilLimit,
                    Error::Type::BasicSeq, bseq.id, msg);
    }

    CoilBseqReport rep;
    rep.bseq_id = bseq.id;
    auto coil = channel_trace(host_seq, seq_idx, coil_chn, L.pts_per_ramp, ms_per_tick);
    auto ttl = channel_trace(host_seq, seq_idx, ttl_chn, L.pts_per_ramp, ms_per_tick);
    double t_end = double(bseq.length) * ms_per_tick;
    if (!coil.t.empty())
        t_end = std::max(t_end, coil.t.back());
    if (!ttl.t.empty())
        t_end = std::max(t_end, ttl.t.back());
    rep.seq_ms = t_end;

    auto to_B = [&] (double v) { return L.b_per_v * v + L.b_offset; };
    auto on = intervals(coil, t_end, [&] (double v) { return to_B(v) > L.thr_G; });
    auto high = intervals(ttl, t_end, [&] (double v) { return v > L.ttl_high; });
    Intervals checked;
    for (auto &w: on) {
        if (parked_window(coil, w.t0_ms, w.t1_ms, L.park_v, L.park_margin_ms)) {
            rep.exempt.push_back(w);
        }
        else {
            checked.push_back(w);
        }
    }
    rep.gaps = uncovered(checked, high);
    rep.windows = on;
    for (auto &w: on)
        rep.on_ms += w.t1_ms - w.t0_ms;
    for (auto &h: high)
        rep.ilk_ms += h.t1_ms - h.t0_ms;
    for (auto v: coil.v)
        rep.peak_G = std::max(rep.peak_G, to_B(v));

    std::vector<std::string> fails;
    auto &R = coil_report;
    R.shot_on_ms += rep.on_ms;
    R.peak_G = std::max(R.peak_G, rep.peak_G);
    if (R.shot_on_ms > L.max_on_ms) {
        fails.push_back(fmt("energized %.1f ms per shot > %.0f ms", R.shot_on_ms, L.max_on_ms));
    }
    else if (R.shot_on_ms > L.warn_on_ms) {
        R.warns.push_back(fmt("energized %.1f ms per shot > warn %.0f ms",
                              R.shot_on_ms, L.warn_on_ms));
    }
    if (rep.peak_G > L.max_B_G) {
        fails.push_back(fmt("field %.1f G > %.0f G", rep.peak_G, L.max_B_G));
    }
    else if (rep.peak_G > L.warn_B_G) {
        R.warns.push_back(fmt("field %.1f G > warn %.0f G", rep.peak_G, L.warn_B_G));
    }
    if (!rep.gaps.empty()) {
        double tot = 0;
        for (auto &g: rep.gaps)
            tot += g.t1_ms - g.t0_ms;
        fails.push_back(fmt("interlockTTL (TTL26) low while energized: %.3f ms, first at %.3f ms",
                            tot, rep.gaps[0].t0_ms));
    }
    // Duty: the largest energized time in any window of length W that holds part of
    // this bseq's energized time, from the intervals of the bseqs that started
    // (history) plus this one placed at `now` (it starts no earlier). The energy in
    // [a, a + W] is piecewise linear in a, so its maximum is at a breakpoint: a or
    // a + W on an interval edge.
    if (!on.empty()) {
        double W = L.duty_window_s;
        double now = m_mgr.coil_clock();
        auto &hist = m_mgr.m_coil_history;
        while (!hist.empty() && hist.front().second < now - 2 * W)
            hist.pop_front();
        std::vector<std::pair<double,double>> ivs(hist.begin(), hist.end());
        for (auto &w: on)
            ivs.emplace_back(now + w.t0_ms * 1e-3, now + w.t1_ms * 1e-3);
        double lo = now + on.front().t0_ms * 1e-3 - W;
        double hi = now + on.back().t1_ms * 1e-3;
        auto energy = [&] (double a) {
            double e = 0;
            for (auto &[x, y]: ivs) {
                double l = std::max(x, a);
                double r = std::min(y, a + W);
                if (r > l) {
                    e += r - l;
                }
            }
            return e;
        };
        double emax = 0;
        auto try_at = [&] (double a) {
            if (a >= lo && a <= hi) {
                emax = std::max(emax, energy(a));
            }
        };
        try_at(lo);
        try_at(hi);
        for (auto &[x, y]: ivs) {
            try_at(x);
            try_at(y);
            try_at(x - W);
            try_at(y - W);
        }
        double duty = emax / W;
        R.duty = std::max(R.duty, duty);
        if (duty > L.max_duty) {
            fails.push_back(fmt("duty %.1f %% of the last %.0f s", 100 * duty, W) +
                            fmt(" > %.0f %%", 100 * L.max_duty, 0));
        }
        else if (duty > L.warn_duty) {
            R.warns.push_back(fmt("duty %.1f %% > warn %.0f %%", 100 * duty, 100 * L.warn_duty));
        }
    }
    m_coil_pending = on;
    rep.check_us = std::chrono::duration<double, std::micro>(
        std::chrono::steady_clock::now() - t_start).count();
    R.bseqs.push_back(std::move(rep));
    if (!fails.empty()) {
        m_coil_pending.clear();
        std::string msg = "Feshbach coil check refused this sequence:";
        for (auto &f: fails) {
            msg += " ";
            msg += f;
            msg += ";";
            R.fails.push_back(f);
        }
        throw Error(Error::Type::BasicSeq, Error::BasicSeq::CoilLimit,
                    Error::Type::BasicSeq, bseq.id, msg);
    }
}

// End of a shot: the duty if this shot repeated at the fastest period the machine can
// run it (compiled length + host_overhead_s per bseq; the real period is never shorter,
// however long the first shots took to build). Same formula as coil_guard.py
// worst_window_duty: n = floor(W / P) whole shots plus at most min(E, W - n P).
// Above the limit, the next pre_run of this scan refuses (cleared by new_run).
void Manager::ExpSeq::coil_predict()
{
    auto &L = m_mgr.m_coil_limits;
    auto &R = coil_report;
    double seq_s = 0;
    for (auto &b: R.bseqs)
        seq_s += b.seq_ms * 1e-3;
    double P = seq_s + L.host_overhead_s * double(R.bseqs.size());
    double E = R.shot_on_ms * 1e-3;
    double W = L.duty_window_s;
    double duty = 0;
    if (E > 0) {
        if (P > 0) {
            double n = std::floor(W / P);
            duty = std::min(1.0, (n * E + std::min(E, W - n * P)) / W);
        }
        else {
            duty = 1;
        }
    }
    R.pred_duty = duty;
    R.pred_period_s = P;
    if (duty > L.max_duty) {
        m_mgr.m_coil_predicted_fail =
            fmt("duty %.1f %% of the worst %.0f s", 100 * duty, W) +
            fmt(" at the full-speed period %.3f s > %.0f %%", P, 100 * L.max_duty);
    }
    else if (duty > L.warn_duty) {
        R.warns.push_back(fmt("predicted duty %.1f %% > warn %.0f %%",
                              100 * duty, 100 * L.warn_duty));
    }
}

void Manager::ExpSeq::coil_commit_start()
{
    double t0 = m_mgr.coil_clock();
    for (auto &w: m_coil_pending)
        m_mgr.m_coil_history.emplace_back(t0 + w.t0_ms * 1e-3, t0 + w.t1_ms * 1e-3);
    m_coil_pending.clear();
}

}
