/*************************************************************************
 *   Copyright (c) 2021 - 2021 Yichao Yu <yyc1992@gmail.com>             *
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

#ifndef __NACS_SEQ_ERROR_H__
#define __NACS_SEQ_ERROR_H__

#include <stdexcept>

namespace NaCs::Seq {

struct Error : std::runtime_error {
    enum ErrorType : uint8_t {
        None,
        EventTime,
        Pulse,
        Measure,
        BasicSeq,
        Assignment,
        Branch,
        Seq
    };

    Error(uint32_t code, ErrorType type1, uint64_t id1, const char *what);
    Error(uint32_t code, ErrorType type1, uint64_t id1, const std::string &what);
    Error(uint32_t code, ErrorType type1, uint64_t id1,
          ErrorType type2, uint64_t id2, const char *what);
    Error(uint32_t code, ErrorType type1, uint64_t id1,
          ErrorType type2, uint64_t id2, const std::string &what);
    Error(const Error&);
    ~Error() override;

    uint32_t code;
    ErrorType type1;
    ErrorType type2 = None;
    uint64_t id1;
    uint64_t id2 = 0;
};

}

#endif
