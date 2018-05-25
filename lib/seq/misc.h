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

#include <nacs-utils/macros.h>
#include <nacs-utils/utils.h>

#include <istream>
#include <string>
#include <vector>
#include <utility>

namespace NaCs {

class WavemeterParser {
    using pos_type = std::istream::pos_type;

    bool try_parsetime(std::istream &stm, double &tsf);
    bool try_parsenumber(std::istream &stm, double &val, bool &eol);
    inline bool try_parsenumber(std::istream &stm, double &val)
    {
        bool eol;
        return try_parsenumber(stm, val, eol);
    }
    bool try_parseline(std::istream &stm);
    bool do_parse(std::istream &stm, bool inc);
    bool match_cache(std::istream &stm) const;

public:
    std::pair<const double*,const double*> parse(std::istream &stm, size_t *sz,
                                                 bool allow_cache);
    WavemeterParser();

private:
    std::vector<double> m_time;
    std::vector<double> m_data;
    pos_type m_last_pos = 0;
    std::string m_last_line;
};

}

extern "C" {

void *nacs_seq_new_wavemeter_parser(void);
size_t nacs_seq_wavemeter_parse(void *parser, const char *name, const double **ts,
                                const double **data, int cache);
void nacs_seq_free_wavemeter_parser(void *parser);

}

#endif
