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

#ifndef _NACS_UTILS_ERRORS_H_
#define _NACS_UTILS_ERRORS_H_

#include "utils.h"

#include <system_error>

namespace NaCs {

// Helper function for C style system calls.
// (return -1 for error and type of error is in errno)
template<typename T, typename... Arg>
static inline T checkErrno(T res, Arg&&... arg)
{
    if (res == -1)
        throw std::system_error(errno, std::system_category(), std::forward<Arg>(arg)...);
    return res;
}

}

#endif
