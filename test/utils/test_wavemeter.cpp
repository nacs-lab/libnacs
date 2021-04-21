/*************************************************************************
 *   Copyright (c) 2018 - 2019 Yichao Yu <yyc1992@gmail.com>             *
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

#include "../../lib/utils/streams.h"
#include "../../lib/utils/timer.h"
#include "../../lib/utils/wavemeter.h"

#include <algorithm>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include <assert.h>

using namespace NaCs;

struct TestFile {
    static constexpr double teps = 1.5e-3;
    TestFile(double t0, double t1, double dt, double _lo, double _hi, int npeaks=5)
        : lo(_lo),
          hi(_hi)
    {
        assert(dt > 10 * teps);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> data_dis(lo, hi);
        std::uniform_real_distribution<> dt_dis(dt * 0.5, dt * 1.5);
        // Populate the data.
        double t = t0;
        while (t < t1) {
            times.push_back(t);
            datas.push_back(data_dis(gen));
            auto dt2 = dt_dis(gen);
            t += dt2;
        }
        // Populate the file
        std::uniform_real_distribution<> strength_dis(-10, 0);
        data_dis = std::uniform_real_distribution<>((3 * lo - hi) / 2, (3 * hi - lo) / 2);

        string_ostream stm;
        // Enough precision for round trip of time and data.
        stm.precision(20);
        // Header
        stm << "Timestamp";
        for (int i = 0; i < npeaks; i++)
            stm << ",Cursor " << i + 1 << " Wavelength,Cursor " << i + 1 << " Intensity";
        stm << std::endl;

        std::vector<double> freqs(npeaks);
        auto gen_line = [&] (double tf, double data, bool valid) {
            double tsecf;
            double tmsf = modf(tf, &tsecf);

            auto tsec = std::time_t(tsecf);
            auto tms = int(tmsf * 1000);
            stm << std::put_time(std::gmtime(&tsec), "%Y-%m-%dT%H:%M:%S")
                << '.' << std::setfill('0') << std::setw(3) << tms;

            if (valid) {
                for (int j = 0; j < npeaks - 1; j++)
                    freqs[j] = data_dis(gen);
                freqs[npeaks - 1] = data;
            }
            else {
                for (int j = 0; j < npeaks; j++) {
                    do {
                        freqs[j] = data_dis(gen);
                    } while (lo <= freqs[j] && freqs[j] <= hi);
                }
            }
            std::sort(freqs.begin(), freqs.end());

            for (int j = 0; j < npeaks; j++) {
                double strength = strength_dis(gen);
                if (valid && freqs[j] == data)
                    strength += 11;
                stm << ',' << freqs[j] << ',' << strength;
            }

            stm << std::endl;
        };

        size_t nline = times.size();
        for (size_t i = 0; ; i++) {
            gen_line(times[i], datas[i], true);
            if (i == nline - 1)
                break;
            auto dt2 = dt_dis(gen);
            auto t2 = times[i] + dt2;
            if (times[i + 1] > t2) {
                gen_line(t2, 0, false);
            }
        }

        file = stm.get_buf();
    }

    void test_parse(Wavemeter &parser, double tstart, double tend)
    {
        const_istream stm(file);
        size_t sz;
        const double *times_ptr;
        const double *datas_ptr;
        const double *heights_ptr;
        std::tie(times_ptr, datas_ptr, heights_ptr) = parser.parse(stm, &sz, tstart, tend);
        if (!times_ptr || *times_ptr <= tstart) {
            // If the time that was right before the requested start time
            // does not have a valid data point we'll omit it from the cache
            // and the second run might expand the data range down.
            // Only check that the cache is taken when we know this doesn't occur.
            size_t sz2;
            auto res = parser.parse(stm, &sz2, tstart, tend);
            assert(times_ptr == nullptr || sz == sz2);
            assert(res == std::make_tuple(times_ptr, datas_ptr, heights_ptr));
        }
        if (times_ptr == nullptr) {
            assert(datas_ptr == nullptr);
            assert(tstart + teps > times.back() ||
                   tend - teps < times.front());
            return;
        }
        assert(datas_ptr != nullptr);
        // The first data must be at least the same as tstart.
        assert(*times_ptr >= tstart);
        // Given the constraint on `dt`, the first data must be `idx_start` or `idx_start + 1`.
        size_t idx_start = std::lower_bound(times.begin(), times.end(),
                                            tstart - teps) - times.begin();
        if (abs(*times_ptr - times[idx_start]) > teps) {
            idx_start++;
            assert(idx_start < times.size());
            assert(abs(*times_ptr - times[idx_start]) <= teps);
        }
        if (idx_start > 0) {
            assert(times[idx_start - 1] < tstart + teps);
        }
        for (size_t i = 0; i < sz; i++) {
            assert(abs(times_ptr[i] - times[idx_start + i]) <= teps);
            assert(datas_ptr[i] == datas[idx_start + i]);
            assert(datas_ptr[i] >= lo);
            assert(datas_ptr[i] <= hi);
        }
        if (idx_start + sz < times.size() - 1) {
            assert(times[idx_start + sz + 1] > tend - teps);
        }
    }

    const double lo;
    const double hi;
    std::string file;
    std::vector<double> times;
    std::vector<double> datas;
};

int main()
{
    auto test0 = [&] (double dt) {
        TestFile test(1472368320, 1472549760, dt, 288000, 289000);

        Wavemeter parser(test.lo, test.hi);
        test.test_parse(parser, 1472368320, 1472541120);
        test.test_parse(parser, 1472368320, 1472541120);
        test.test_parse(parser, 1472550000, 1472580000);

        parser.clear();
        test.test_parse(parser, 1472368320, 1472376960);
        test.test_parse(parser, 1472368320, 1472376960);
        test.test_parse(parser, 1472385600, 1472411520);
        test.test_parse(parser, 1472385600, 1472411520);
        test.test_parse(parser, 1472428800, 1472515200);
        test.test_parse(parser, 1472428800, 1472515200);
        test.test_parse(parser, 1472368320, 1472549760);
        test.test_parse(parser, 1472368320, 1472549760);

        parser.clear();
        test.test_parse(parser, 1472428800, 1472515200);
        test.test_parse(parser, 1472385600, 1472411520);
        test.test_parse(parser, 1472368320, 1472376960);
        test.test_parse(parser, 1472368320, 1472549760);

        parser.clear();
        test.test_parse(parser, 1472428800, 1472515200);
        test.test_parse(parser, 1472385600, 1472411520);
        test.test_parse(parser, 1472402880, 1472437440);
        test.test_parse(parser, 1472368320, 1472549760);
    };

    test0(10);
    test0(1);
    test0(0.3);

    auto test1 = [&] () {
        TestFile test(1539887328, 1543334688, 5, 288000, 289000);

        Wavemeter parser(test.lo, test.hi);
        test.test_parse(parser, 1543190400, 1543363200);
        test.test_parse(parser, 1543190400, 1543363200);
        test.test_parse(parser, 0, std::numeric_limits<double>::max());
    };
    test1();

    return 0;
}
