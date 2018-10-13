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

#include "utils.h"

#include <zmq.hpp>
#include <string.h>

#ifndef __NACS_UTILS_ZMQ_H__
#define __NACS_UTILS_ZMQ_H__

namespace NaCs {
namespace ZMQ {

static inline bool has_more(zmq::socket_t &sock)
{
    int more = 0;
    size_t more_size = sizeof(more);
    sock.getsockopt(ZMQ_RCVMORE, &more, &more_size);
    return more != 0;
}

static inline bool recv_more(zmq::socket_t &sock, zmq::message_t &msg)
{
    if (has_more(sock)) {
        sock.recv(&msg);
        return true;
    }
    return false;
}

static inline void readall(zmq::socket_t &sock)
{
    while (has_more(sock)) {
        zmq::message_t msg;
        sock.recv(&msg);
    }
}

static inline bool match(const zmq::message_t &msg, const char *str)
{
    auto len = strlen(str);
    if (msg.size() != len)
        return false;
    return memcmp(msg.data(), str, len) == 0;
}

template<typename T>
static inline zmq::message_t bits_msg(T v)
{
    zmq::message_t msg(sizeof(T));
    std::memcpy(msg.data(), &v, sizeof(T));
    return msg;
}

static inline zmq::message_t str_msg(const char *str)
{
    return zmq::message_t(str, strlen(str));
}

}
}

#endif
