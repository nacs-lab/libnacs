/*************************************************************************
 *   Copyright (c) 2019 - 2019 Yichao Yu <yyc1992@gmail.com>             *
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

#include "tweezer.h"

#include <cmath>

namespace NaCs::Tweezer {

NACS_EXPORT() double sideband(double eta, int n1, int n2)
{
    if (n1 < 0 || n2 < 0)
        return 0;
    if (eta == 0)
        return n1 == n2 ? 1 : 0;
    // Ref http://journals.aps.org/pra/pdf/10.1103/PhysRevA.20.1521
    // Δn ≡ |n1 - n2|
    // n₋ ≡ min(n1, n2)
    // n₊ ≡ max(n1, n2)
    //   ⟨n1|exp(ikx)|n2⟩
    // = ⟨n1|exp(iη(a + a†))|n2⟩
    // = exp(-η^2 / 2) η^Δn √(γ(n₋ + 1) / γ(n₊ + 1)) L^Δn_n₋(η^2)
    // = exp(-η^2 / 2 + Δn log(η) + lγ(n₋ + 1) / 2 - lγ(n₊ + 1) / 2) L^Δn_n₋(η^2)
    auto eta2 = eta * eta;
    double lpre;
    double lag;
    if (n1 == n2) {
        lpre = -0.5 * eta2;
        lag = std::assoc_laguerre(n1, 0, eta2);
    }
    else {
        int n_m = std::min(n1, n2);
        int n_p = std::max(n1, n2);
        int dn = n_p - n_m;
        lpre = (-eta2 + std::lgamma(n_m + 1) - std::lgamma(n_p + 1)) / 2 + std::log(eta) * dn;
        lag = std::assoc_laguerre(n_m, dn, eta2);
    }
    return lag * std::exp(lpre);
}

static double find_max_sideband(double eta, int order)
{
    double mmax = sideband(eta, 0, order);
    for (int i = 1; ; i++) {
        double m = sideband(eta, i, i + order);
        if (m <= mmax)
            return mmax;
        mmax = m;
    }
}

NACS_EXPORT() void max_sideband_ratios(double eta, int orders, double *ratios)
{
    double m0 = sideband(eta, 0, 1);
    for (int i = 0; i < orders; i++) {
        ratios[i] = m0 / find_max_sideband(eta, i + 1);
    }
}

}

extern "C" NACS_EXPORT() double nacs_tweezer_sideband(double eta, int n1, int n2)
{
    return NaCs::Tweezer::sideband(eta, n1, n2);
}

extern "C" NACS_EXPORT() void nacs_tweezer_max_sideband_ratios(
    double eta, int orders, double *ratios)
{
    return NaCs::Tweezer::max_sideband_ratios(eta, orders, ratios);
}
