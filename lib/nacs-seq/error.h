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

#include <nacs-utils/utils.h>

#include <stdexcept>
#include <type_traits>

namespace NaCs::Seq {

struct NACS_EXPORT_ Error : std::runtime_error {
    enum class Type : uint8_t {
        None,
        EventTime,
        Pulse,
        Measure,
        BasicSeq,
        Assignment,
        Branch,
        Seq,
        Channel,
    };
    enum class EventTime : uint8_t {
        NegTime,
        NonPosTime,
    };
    enum class Pulse : uint8_t {
        NegTime,
    };

    Error(Type type, uint16_t code, Type type1, uint64_t id1, const char *what);
    Error(Type type, uint16_t code, Type type1, uint64_t id1, const std::string &what);
    Error(Type type, uint16_t code, Type type1, uint64_t id1,
          Type type2, uint64_t id2, const char *what);
    Error(Type type, uint16_t code, Type type1, uint64_t id1,
          Type type2, uint64_t id2, const std::string &what);
    template<typename Code, typename... Args,
             typename=std::enable_if_t<!std::is_same_v<
                 std::remove_cv_t<std::remove_reference_t<Code>>,uint16_t>>>
        Error(Type type, Code code, Args&&... args)
        : Error(type, uint16_t(code), std::forward<Args>(args)...)
    {}
    Error(const Error&);
    ~Error() override;

    Type type;
    Type type1;
    Type type2 = Type::None;
    uint16_t code;
    uint64_t id1;
    uint64_t id2 = 0;
};

}

#endif
