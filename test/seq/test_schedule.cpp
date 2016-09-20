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
#include <functional>

using namespace NaCs;

typedef uint32_t cid_t;

struct accum_t {
    uint64_t operator()(cid_t, double, uint64_t)
    {
        return 1;
    }
};

typedef Seq::Pulse<cid_t,std::function<
                             double(uint64_t,double,uint64_t)>> pulse_t;

struct filter_t {
    bool operator()(cid_t cid)
    {
        return true;
    }
    bool operator()(cid_t cid, double val1, double val2)
    {
        return true;
    }
};

static const auto seq_cb = [&] (auto &accum, uint64_t cur_t, Seq::Event evt) {
    return cur_t;
};

// template<typename Accum, typename Cid, typename Cb, typename Filter,
//          typename SeqCB, typename ValT>
// static void schedule(Accum &accum, std::vector<Pulse<Cid,Cb>> &seq,
//                      const Time::Constraints &t_cons,
//                      const std::map<Cid,ValT> &defaults,
//                      Filter &&filter, SeqCB &&seq_cb)

int main()
{
    accum_t accum;
    std::vector<pulse_t> seq;
    Seq::Time::Constraints t_cons{1, 1, 1};
    std::map<cid_t,double> defaults;
    filter_t filter;
    schedule(accum, seq, t_cons, defaults, filter, seq_cb);
    return 0;
}
