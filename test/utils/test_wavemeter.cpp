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
    TestFile(double t0, double t1, double dt, double _lo, double _hi, int npeaks=5)
        : lo(_lo),
          hi(_hi)
    {
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
        // Header
        stm << "Timestamp";
        for (int i = 0; i < npeaks; i++)
            stm << ",Cursor " << i + 1 << " Wavelength,Cursor " << i + 1 << " Intensity";
        stm << std::endl;

        std::vector<double> freqs(npeaks);
        size_t nline = times.size();
        for (size_t i = 0; i < nline; i++) {
            double tf = (times[i] - 719529) * 86400;
            double tsecf;
            double tmsf = modf(tf, &tsecf);

            auto tsec = std::time_t(tsecf);
            auto tms = int(tmsf * 1000);
            stm << std::put_time(std::gmtime(&tsec), "%Y-%m-%dT%H:%M:%S")
                << '.' << std::setfill('0') << std::setw(3) << tms;

            for (int j = 0; j < npeaks - 1; j++)
                freqs[j] = data_dis(gen);
            freqs[npeaks - 1] = datas[i];
            std::sort(freqs.begin(), freqs.end());

            for (int j = 0; j < npeaks; j++) {
                double strength = strength_dis(gen);
                if (freqs[j] == datas[i])
                    strength += 11;
                stm << ',' << freqs[j] << ',' << strength;
            }

            stm << std::endl;
        }

        file = stm.get_buf();
    }

    const double lo;
    const double hi;
    std::string file;
    std::vector<double> times;
    std::vector<double> datas;
};

int main()
{
    return 0;
}
