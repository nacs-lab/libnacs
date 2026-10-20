/*************************************************************************
 *   Copyright (c) 2026 NaCs lab                                         *
 *                                                                       *
 *   This library is free software; you can redistribute it and/or       *
 *   modify it under the terms of the GNU Lesser General Public          *
 *   License as published by the Free Software Foundation; either        *
 *   version 3.0 of the License, or (at your option) any later version.  *
 *************************************************************************/

// AWG RF power check in `pre_run` (enabled by the `awg_power` config block).
//
// For the basic sequence about to run (its globals already applied, so this includes the
// rearrangement moves), rebuild every AMP / FREQ tone trace of each checked AWG output
// the way `get_nominal_output` does, then:
//   power   sum over the tones of one output of amp^2 <= max_sum_amp2 at every time
//           (awgInterlock.m)
//   freq    while a tone's amp is above amp_on, its frequency stays inside the
//           output's window, if one is configured (awgFreqInterlock.m)
// A violation throws Error(BasicSeq, AwgLimit) before any device sees the bseq.
//
// Cost (this runs on the rearrangement branch, between the image and the moves):
//   - the AWG pulses of each bseq are found once, then only those are read per shot;
//   - buffers are kept between shots;
//   - a quick bound first (each tone at its own peak, all at once); only if it is
//     above the warn level does the time sweep run (sorted change points, running sum);
//   - the frequency walk runs only if some frequency sample is outside its window.
// Ramps are sampled at pts_per_ramp points; between two samples of one ramp the larger
// end is used, so a monotonic ramp is never under-counted.

#include "manager.h"
#include "error.h"

#include "../nacs-utils/number.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>

namespace NaCs::Seq {

namespace awg_detail {

struct Sample {
    int64_t t; // ticks
    double v;
    bool cont; // the value ramps from the previous sample to this one
};

struct P {
    uint32_t id;
    int64_t time;
    int64_t len;
    double (*ramp_func)(double, void*);
    double endvalue;
    bool cond;
};

struct Event {
    int64_t t;
    uint32_t tone;
    double a2; // the tone's amp^2 bound from t until its next event
};

// Same walk as Manager::ExpSeq::get_nominal_output (and coil_check's channel_trace),
// with the ramp-continuation flag kept.
static void build_trace(const HostSeq &host_seq, uint32_t chn, std::vector<P> &pulses,
                        uint64_t pts_per_ramp, std::vector<Sample> &out)
{
    out.clear();
    std::sort(pulses.begin(), pulses.end(), [] (auto &p1, auto &p2) {
        if (p1.time != p2.time)
            return p1.time < p2.time;
        return p1.id < p2.id;
    });
    out.push_back({0, host_seq.start_values[chn - 1].f64, false});
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
            if (out.back().t != pulse.time)
                out.push_back({pulse.time, out.back().v, false});
            out.push_back({pulse.time, pulse.endvalue, false});
        }
        else {
            int64_t end_time = pulse.time + pulse.len;
            if (next)
                end_time = std::min(next->time, end_time);
            int64_t ramp_time = end_time - pulse.time;
            if (ramp_time <= 0)
                continue;
            if (out.back().t != pulse.time)
                out.push_back({pulse.time, out.back().v, false});
            uint64_t npts = std::max<uint64_t>(2, std::min((uint64_t)ramp_time, pts_per_ramp));
            for (uint64_t j = 0; j < npts; j++) {
                int64_t this_pt = round<int64_t>(j * ramp_time / (npts - 1));
                out.push_back({pulse.time + this_pt,
                               pulse.ramp_func(this_pt, (void*)host_seq.values.data()),
                               j > 0});
            }
        }
    }
}

// Walks one trace forward in time. `range` = the values over [t, next change): the held
// sample, and the next one too if a ramp leads to it.
struct Cursor {
    const std::vector<Sample> *tr = nullptr;
    size_t p = 0;
    void advance(int64_t t)
    {
        auto &s = *tr;
        while (p + 1 < s.size() && s[p + 1].t <= t) {
            p++;
        }
    }
    void range(double &lo, double &hi) const
    {
        auto &s = *tr;
        lo = hi = s[p].v;
        if (p + 1 < s.size() && s[p + 1].cont) {
            lo = std::min(lo, s[p + 1].v);
            hi = std::max(hi, s[p + 1].v);
        }
    }
};

static std::string fmt(const char *f, const std::string &a, double b, double c)
{
    char buf[256];
    snprintf(buf, sizeof(buf), f, a.c_str(), b, c);
    return buf;
}

}

using namespace awg_detail;

// Kept between shots so the check does not allocate once warmed up.
struct Manager::ExpSeq::AwgScratch {
    std::vector<std::vector<std::vector<P>>> amp_p; // [out][tone]
    std::vector<std::vector<std::vector<P>>> freq_p;
    std::vector<std::vector<Sample>> amp_tr; // [tone], reused per out
    std::vector<std::vector<Sample>> freq_tr;
    std::vector<double> peak;
    std::vector<char> freq_bad; // [tone] some freq sample outside the window
    std::vector<double> cur;
    std::vector<Event> events;
    std::vector<Event> events2;
    std::vector<size_t> runs; // end of each tone's (sorted) run in `events`
};

void Manager::ExpSeq::awg_add_channel(uint32_t chn_id, const std::string &dev,
                                      const std::string &chn)
{
    auto &devs = m_mgr.m_awg_limits.devices;
    if (std::find(devs.begin(), devs.end(), dev) == devs.end())
        return;
    // OUT<n>/CHN<m>/AMP or OUT<n>/CHN<m>/FREQ; anything else (PHASE, A<k>) is ignored.
    auto s1 = chn.find('/');
    if (s1 == chn.npos)
        return;
    auto s2 = chn.find('/', s1 + 1);
    if (s2 == chn.npos)
        return;
    auto out_str = chn.substr(0, s1);
    auto tone_str = chn.substr(s1 + 1, s2 - s1 - 1);
    auto type_str = chn.substr(s2 + 1);
    bool is_amp = type_str == "AMP";
    if (!is_amp && type_str != "FREQ")
        return;
    if (out_str.compare(0, 3, "OUT") != 0 || tone_str.compare(0, 3, "CHN") != 0)
        return;
    uint32_t chn_num;
    try {
        chn_num = (uint32_t)std::stoul(tone_str.substr(3));
    }
    catch (...) {
        return;
    }
    auto out_name = dev + "/" + out_str;
    int32_t out_idx = -1;
    for (size_t i = 0; i < awg_outs.size(); i++) {
        if (awg_outs[i].name == out_name) {
            out_idx = int32_t(i);
            break;
        }
    }
    if (out_idx < 0) {
        out_idx = int32_t(awg_outs.size());
        awg_outs.push_back({out_name, {}});
    }
    auto &tones = awg_outs[out_idx].tones;
    uint32_t tone_idx = uint32_t(tones.size());
    for (uint32_t i = 0; i < tones.size(); i++) {
        if (tones[i].chn_num == chn_num) {
            tone_idx = i;
            break;
        }
    }
    if (tone_idx == tones.size())
        tones.push_back({chn_num, 0, 0});
    if (is_amp) {
        tones[tone_idx].amp_chn = chn_id;
    }
    else {
        tones[tone_idx].freq_chn = chn_id;
    }
    if (m_awg_chn_ref.size() <= chn_id)
        m_awg_chn_ref.resize(chn_id + 1);
    m_awg_chn_ref[chn_id] = {out_idx, tone_idx, is_amp};
}

void Manager::ExpSeq::awg_check()
{
    auto t_start = std::chrono::steady_clock::now();
    auto &L = m_mgr.m_awg_limits;
    auto seq_idx = host_seq.cur_seq_idx();
    auto &bseq = host_seq.seqs[seq_idx];
    double ms_per_tick = 1e3 / double(m_mgr.tick_per_sec());
    if (!m_awg_scratch)
        m_awg_scratch = std::make_shared<AwgScratch>();
    auto &S = *m_awg_scratch;

    // The AWG pulses of this bseq (fixed at compile time; found once).
    if (m_awg_pulse_idx.size() < host_seq.seqs.size()) {
        m_awg_pulse_idx.resize(host_seq.seqs.size());
        m_awg_pulse_idx_done.resize(host_seq.seqs.size(), false);
    }
    auto &idxs = m_awg_pulse_idx[seq_idx];
    if (!m_awg_pulse_idx_done[seq_idx]) {
        auto npulses = uint32_t(bseq.pulses.size());
        for (uint32_t i = 0; i < npulses; i++) {
            auto &pulse = bseq.pulses[i];
            if (pulse.is_measure() || pulse.chn >= m_awg_chn_ref.size())
                continue;
            if (m_awg_chn_ref[pulse.chn].out < 0)
                continue;
            idxs.push_back(i);
        }
        m_awg_pulse_idx_done[seq_idx] = true;
    }

    // Pulses of each tone channel, with this shot's times and values.
    S.amp_p.resize(awg_outs.size());
    S.freq_p.resize(awg_outs.size());
    for (size_t o = 0; o < awg_outs.size(); o++) {
        S.amp_p[o].resize(awg_outs[o].tones.size());
        S.freq_p[o].resize(awg_outs[o].tones.size());
        for (auto &v: S.amp_p[o])
            v.clear();
        for (auto &v: S.freq_p[o]) {
            v.clear();
        }
    }
    for (auto i: idxs) {
        auto &pulse = bseq.pulses[i];
        auto &ref = m_awg_chn_ref[pulse.chn];
        int64_t time = host_seq.get_time(pulse.time);
        int64_t len = pulse.len == uint32_t(-1) ? 0 :
            round<int64_t>(host_seq.get_value(pulse.len, seq_idx));
        double endvalue = host_seq.get_value(pulse.endvalue, seq_idx);
        bool cond = pulse.cond == uint32_t(-1) ? true :
            (bool)host_seq.get_value(pulse.cond, seq_idx);
        auto &dst = ref.is_amp ? S.amp_p[ref.out][ref.tone] : S.freq_p[ref.out][ref.tone];
        dst.push_back({pulse.id, time, len, pulse.ramp_func, endvalue, cond});
    }

    AwgBseqReport rep;
    rep.bseq_id = bseq.id;
    std::vector<std::string> fails;
    for (size_t o = 0; o < awg_outs.size(); o++) {
        auto &out = awg_outs[o];
        auto ntones = out.tones.size();
        if (S.amp_tr.size() < ntones) {
            S.amp_tr.resize(ntones);
            S.freq_tr.resize(ntones);
        }
        S.peak.assign(ntones, 0);
        const std::pair<double,double> *win = nullptr;
        if (auto it = L.freq_windows.find(out.name); it != L.freq_windows.end())
            win = &it->second;

        // Quick bound: every tone at its own peak at the same time.
        AwgOutReport orep;
        orep.name = out.name;
        S.freq_bad.assign(ntones, 0);
        for (size_t k = 0; k < ntones; k++) {
            auto &tone = out.tones[k];
            S.amp_tr[k].clear();
            S.freq_tr[k].clear();
            if (!tone.amp_chn)
                continue;
            build_trace(host_seq, tone.amp_chn, S.amp_p[o][k], L.pts_per_ramp, S.amp_tr[k]);
            for (auto &s: S.amp_tr[k])
                S.peak[k] = std::max(S.peak[k], std::abs(s.v));
            orep.bound += S.peak[k] * S.peak[k];
            if (win && tone.freq_chn && S.peak[k] > L.amp_on) {
                build_trace(host_seq, tone.freq_chn, S.freq_p[o][k], L.pts_per_ramp,
                            S.freq_tr[k]);
                for (auto &s: S.freq_tr[k]) {
                    if (!(s.v >= win->first && s.v <= win->second)) {
                        S.freq_bad[k] = 1;
                        break;
                    }
                }
            }
        }

        if (orep.bound > L.warn_sum_amp2) {
            // Time sweep: each tone's amp^2 bound changes only at its own samples.
            orep.exact = true;
            S.events.clear();
            S.runs.clear();
            for (size_t k = 0; k < ntones; k++) {
                auto &s = S.amp_tr[k];
                auto n = s.size();
                for (size_t p = 0; p < n; p++) {
                    // A later sample at the same time replaces this one.
                    if (p + 1 < n && s[p + 1].t == s[p].t)
                        continue;
                    double a = std::abs(s[p].v);
                    if (p + 1 < n && s[p + 1].cont)
                        a = std::max(a, std::abs(s[p + 1].v));
                    S.events.push_back({s[p].t, uint32_t(k), a * a});
                }
                S.runs.push_back(S.events.size());
            }
            // Each tone's events are already in time order: merge the runs pairwise
            // (n log(tones)) instead of sorting everything (n log n).
            // Ping-pong between two kept buffers, so nothing is allocated here.
            auto by_time = [] (auto &e1, auto &e2) { return e1.t < e2.t; };
            S.events2.resize(S.events.size());
            while (S.runs.size() > 1) {
                size_t w = 0;
                size_t begin = 0;
                auto src = S.events.begin();
                auto dst = S.events2.begin();
                for (size_t r = 0; r < S.runs.size(); r += 2) {
                    size_t end = r + 1 < S.runs.size() ? S.runs[r + 1] : S.runs[r];
                    std::merge(src + begin, src + S.runs[r], src + S.runs[r], src + end,
                               dst + begin, by_time);
                    begin = end;
                    S.runs[w++] = end;
                }
                S.runs.resize(w);
                std::swap(S.events, S.events2);
            }
            S.cur.assign(ntones, 0);
            double sum = 0;
            auto nev = S.events.size();
            for (size_t i = 0; i < nev;) {
                auto t = S.events[i].t;
                for (; i < nev && S.events[i].t == t; i++) {
                    auto &e = S.events[i];
                    sum += e.a2 - S.cur[e.tone];
                    S.cur[e.tone] = e.a2;
                }
                if (sum > orep.sum_amp2) {
                    // Recount, so the reported value carries no rounding from the
                    // running sum.
                    double exact = 0;
                    for (auto c: S.cur)
                        exact += c;
                    sum = exact;
                    if (exact > orep.sum_amp2) {
                        orep.sum_amp2 = exact;
                        orep.t_ms = double(t) * ms_per_tick;
                    }
                }
            }
            if (orep.sum_amp2 > L.max_sum_amp2) {
                fails.push_back(fmt("%s sum(amp^2) %.4f > %.4f", out.name, orep.sum_amp2,
                                    L.max_sum_amp2) +
                                fmt("%s at %.3f ms", "", orep.t_ms, 0));
            }
            else if (orep.sum_amp2 > L.warn_sum_amp2) {
                awg_report.warns.push_back(fmt("%s sum(amp^2) %.4f > warn %.4f", out.name,
                                               orep.sum_amp2, L.warn_sum_amp2));
            }
        }

        // A tone with a frequency sample outside the window (e.g. parked at 0 Hz with
        // its amp at 0): walk its own amp and freq traces together; is it on there?
        for (size_t k = 0; k < ntones; k++) {
            if (!S.freq_bad[k])
                continue;
            Cursor ac{&S.amp_tr[k]};
            Cursor fc{&S.freq_tr[k]};
            int64_t t = 0;
            while (true) {
                ac.advance(t);
                fc.advance(t);
                double lo, hi;
                ac.range(lo, hi);
                double a = std::max(std::abs(lo), std::abs(hi));
                double flo, fhi;
                fc.range(flo, fhi);
                if (a > L.amp_on && !(flo >= win->first && fhi <= win->second)) {
                    char buf[256];
                    snprintf(buf, sizeof(buf),
                             "%s CHN%u frequency %.6g Hz outside [%.6g, %.6g] Hz "
                             "at %.3f ms (amp %.4g)", out.name.c_str(),
                             out.tones[k].chn_num, flo < win->first ? flo : fhi,
                             win->first, win->second, double(t) * ms_per_tick, a);
                    fails.push_back(buf);
                    break;
                }
                // Next change of either trace.
                auto &as = S.amp_tr[k];
                auto &fs = S.freq_tr[k];
                bool more_a = ac.p + 1 < as.size();
                bool more_f = fc.p + 1 < fs.size();
                if (!more_a && !more_f)
                    break;
                t = !more_f ? as[ac.p + 1].t : !more_a ? fs[fc.p + 1].t :
                    std::min(as[ac.p + 1].t, fs[fc.p + 1].t);
            }
        }
        rep.outs.push_back(std::move(orep));
    }
    rep.check_us = std::chrono::duration<double, std::micro>(
        std::chrono::steady_clock::now() - t_start).count();
    awg_report.bseqs.push_back(std::move(rep));
    if (!fails.empty()) {
        std::string msg = "AWG power check refused this sequence:";
        for (auto &f: fails) {
            msg += " ";
            msg += f;
            msg += ";";
            awg_report.fails.push_back(f);
        }
        throw Error(Error::Type::BasicSeq, Error::BasicSeq::AwgLimit,
                    Error::Type::BasicSeq, bseq.id, msg);
    }
}

}
