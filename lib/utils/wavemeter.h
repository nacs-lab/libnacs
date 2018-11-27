/*************************************************************************
 *   Copyright (c) 2018 - 2018 Yichao Yu <yyc1992@gmail.com>             *
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

#ifndef __NACS_UTILS_WAVEMETER_H__
#define __NACS_UTILS_WAVEMETER_H__

#include <istream>

namespace NaCs {

class Wavemeter {
    // Stateless parsing functions
    static bool try_parsetime(std::istream &stm, double *tsf);
    static bool try_parsenumber(std::istream &stm, double *val, bool *eol);
    static inline bool try_parsenumber(std::istream &stm, double *val)
    {
        bool eol;
        return try_parsenumber(stm, val, &eol);
    }
    static bool try_parseval(std::istream &stm, double *val, double lo, double hi);
    bool try_parseline(std::istream &stm, double *tsf, double *val) const;

    using pos_type = std::istream::pos_type;
    static pos_type find_linestart(std::istream &stm, pos_type ub, pos_type lb);
};

}

#endif // __NACS_UTILS_WAVEMETER_H__
