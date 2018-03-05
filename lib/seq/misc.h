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

#ifndef __NACS_SEQ_MISC_H__
#define __NACS_SEQ_MISC_H__

/**
 * These are some functions that are useful for the experiment.
 * They are in this library mainly to make my life easier on windows........
 * They can probably be moved to `nacs-utils` but I'm too lazy to try....
 */

#include <nacs-utils/utils.h>

#include <istream>
#include <locale>
#include <string>
#include <vector>
#include <utility>

namespace NaCs {

class WavemeterParser {
    bool try_parseline(std::istream &stm);
    bool do_parse(std::istream &stm, bool inc);
    using pos_type = std::istream::pos_type;

public:
    std::pair<const double*,const double*> parse(std::istream &stm, size_t *sz,
                                                 bool allow_cache);

private:
    std::vector<double> m_time;
    std::vector<double> m_data;
    std::locale m_locale{"en_US.utf-8"};
    pos_type m_last_pos;
    std::string m_last_line;
};

}

#endif
