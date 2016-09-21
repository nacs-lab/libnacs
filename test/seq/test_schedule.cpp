/*************************************************************************
 *   Copyright (c) 2016 - 2016 Yichao Yu <yyc1992@gmail.com>             *
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

#include <nacs-seq/seq.h>
#include <nacs-seq/pulser.h>
#include <nacs-utils/log.h>
#include <functional>

#include <math.h>
#include <inttypes.h>

using namespace NaCs;

struct cid_t {
    uint32_t id;
};
bool operator<(cid_t id1, cid_t id2)
{
    return id1.id < id2.id;
}

struct val_t {
    double v;
    val_t(double _v=0)
        : v(_v)
    {}
};

struct accum_t {
    uint64_t operator()(Seq::Channel chn, val_t val, uint64_t t)
    {
        // nacsLog("t = %" PRIu64 ", chn = %d, v = %f\n", t, chn.id, val.v);
        return 1;
    }
};

typedef Seq::Pulse<Seq::Channel,std::function<
                                    val_t(uint64_t,val_t,uint64_t)>> pulse_t;

struct filter_t {
    bool operator()(Seq::Channel cid)
    {
        return true;
    }
    bool operator()(Seq::Channel cid, val_t val1, val_t val2)
    {
        return val1.v != val2.v;
    }
};

static const auto seq_cb = [&] (auto &accum, uint64_t cur_t, Seq::Event evt) {
    nacsLog("Start time: %" PRIu64 "\n", cur_t);
    return cur_t;
};

int main()
{
    // uint64_t tlen = 10000000;
    uint64_t tlen = 10000;
    accum_t accum;
    std::vector<pulse_t> seq = {
        {0, tlen, {Seq::Channel::DAC, 0}, [] (auto t, auto start, auto len) {
             return val_t(start.v + (double)t);
         }},
        {0, tlen, {Seq::Channel::DAC, 1}, [] (auto t, auto start, auto len) {
             return val_t(start.v - (double)t);
         }},
        {0, tlen, {Seq::Channel::DAC, 3}, [] (auto t, auto start, auto len) {
             return val_t(start.v + sin((double)t / 1000.0));
         }},
        {5000, 0, {Seq::Channel::DAC, 2}, [] (auto t, auto start, auto len) {
             return val_t(20);
         }},
        {5002, 0, {Seq::Channel::DAC, 2}, [] (auto t, auto start, auto len) {
             return val_t(10);
         }},
        {5202, 0, {Seq::Channel::DAC, 2}, [] (auto t, auto start, auto len) {
             return val_t(10);
         }}
    };
    Seq::Time::Constraints t_cons{100, 40, 4096};
    std::map<Seq::Channel,val_t> defaults;
    filter_t filter;
    sort(seq);
    schedule(accum, seq, t_cons, defaults, filter, seq_cb);
    return 0;
}
