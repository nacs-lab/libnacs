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
#include <nacs-utils/timer.h>
#include <nacs-utils/log.h>
#include <nacs-utils/base64.h>

#include <math.h>
#include <inttypes.h>
#include <iostream>
#include <fstream>
#include <ostream>
#include <iomanip>

using namespace NaCs;

static const auto seq_cb = [&] (auto &, uint64_t cur_t, Seq::Event evt) {
    if (evt == Seq::Event::start) {
        nacsLog("Start time: %" PRIu64 "\n", cur_t);
    } else {
        nacsLog("End time: %" PRIu64 "\n", cur_t);
    }
    return cur_t;
};

int main(int argn, char **argv)
{
    assert(argn == 2);
    uint64_t tlen = 10000000;
    // uint64_t tlen = 10000;
    Seq::PulsesBuilder builder =
        [] (Seq::Channel chn, Seq::Val val, uint64_t t) -> uint64_t {
        (void)chn;
        (void)val;
        (void)t;
        // nacsLog("t = %" PRIu64 ", chn = (%d, %d), v = %f\n",
        //         t, int(chn.typ), chn.id, val.val.f64);
        return 1;
    };
    std::vector<Seq::Pulse> seq;
    seq.push_back(Seq::Pulse{0, tlen, {Seq::Channel::DAC, 0},
                [] (auto t, auto start, auto) {
                    return Seq::Val::get<double>(start.val.f64 + (double)t);
                }});
    seq.push_back(Seq::Pulse{0, tlen, {Seq::Channel::DAC, 1},
                [] (auto t, auto start, auto) {
                    return Seq::Val::get<double>(start.val.f64 - (double)t);
                }});
    seq.push_back(Seq::Pulse{0, tlen, {Seq::Channel::DAC, 3},
                [] (auto t, auto start, auto) {
                    return Seq::Val::get<double>(start.val.f64 +
                                                 sin((double)t / 1000.0));
                }});
    seq.push_back(Seq::Pulse{5000, 0, {Seq::Channel::DAC, 2},
                Seq::Val::get<double>(20)});
    seq.push_back(Seq::Pulse{5002, 0, {Seq::Channel::DAC, 2},
                Seq::Val::get<double>(10)});
    seq.push_back(Seq::Pulse{5202, 0, {Seq::Channel::DAC, 2},
                Seq::Val::get<double>(10)});
    std::map<Seq::Channel,Seq::Val> defaults;
    builder.schedule(seq, defaults, seq_cb);

    std::ofstream ostm;
    Seq::PulsesBuilder file_builder =
        [&] (Seq::Channel chn, Seq::Val val, uint64_t t) -> uint64_t {
        uint64_t min_time = 50;
        const char *channel_name;
        switch (chn.typ) {
        case Seq::Channel::TTL:
            assert(chn.id == 0);
            channel_name = "TTL";
            min_time = 3;
            break;
        case Seq::Channel::DDS_FREQ:
            channel_name = "FREQ";
            break;
        case Seq::Channel::DDS_AMP:
            channel_name = "AMP";
            break;
        case Seq::Channel::DAC:
            channel_name = "DAC";
            break;
        default:
            assert(0 && "Unknown channel type.");
        }
        if (!ostm.is_open())
            return min_time;
        ostm << t << ", " << channel_name;
        if (chn.typ == Seq::Channel::TTL) {
            ostm << "(all) = 0x" << std::hex << val.val.i32 << std::dec
                 << std::endl;
        } else {
            ostm << "(" << chn.id << ") = "
            << std::setprecision(18) << val.val.f64 << std::endl;
        }
        return min_time;
    };

    {
        nacsLog("NaSingleAtom\n");
        std::ifstream istm(std::string(argv[1]) + "/na_single_atom.seq");
        ostm.open("na_single_atom.txt");
        std::string data(std::istreambuf_iterator<char>(istm), {});
        assert(Base64::validate((const uint8_t*)data.data(), data.size()));
        std::vector<Seq::Pulse> seq;
        std::map<Seq::Channel,Seq::Val> defaults;
        tic();
        std::tie(seq, defaults) =
            Seq::PulsesBuilder::fromBase64((const uint8_t*)data.data(),
                                           data.size());
        printToc();
        tic();
        file_builder.schedule(seq, defaults, seq_cb);
        printToc();
        ostm.close();
    }

    {
        nacsLog("CsSingleAtom\n");
        std::ifstream istm(std::string(argv[1]) + "/cs_single_atom.seq");
        ostm.open("cs_single_atom.txt");
        std::string data(std::istreambuf_iterator<char>(istm), {});
        assert(Base64::validate((const uint8_t*)data.data(), data.size()));
        std::vector<Seq::Pulse> seq;
        std::map<Seq::Channel,Seq::Val> defaults;
        tic();
        std::tie(seq, defaults) =
            Seq::PulsesBuilder::fromBase64((const uint8_t*)data.data(),
                                           data.size());
        printToc();
        tic();
        file_builder.schedule(seq, defaults, seq_cb);
        printToc();
        ostm.close();
    }
    return 0;
}
