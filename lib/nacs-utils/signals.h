/*************************************************************************
 *   Copyright (c) 2014 - 2016 Yichao Yu <yyc1992@gmail.com>             *
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

#include <functional>
#include <list>

#ifndef __NACS_UTILS_SIGNALS_H__
#define __NACS_UTILS_SIGNALS_H__

namespace NaCs {

template<typename... Arg>
struct Signal {
    using callback_t = std::function<void(Arg...)>;
    using connect_id_t = typename std::list<callback_t>::iterator;

    template<typename CB>
    connect_id_t connect(CB &&cb)
    {
        m_handlers.emplace_back(std::forward<CB>(cb));
        return --m_handlers.end();
    }

    void disconnect(connect_id_t id)
    {
        m_handlers.erase(id);
    }

    void emit(Arg&... args) const
    {
        for (auto &cb: m_handlers) {
            try {
                cb(args...);
            }
            catch (...) {
                // Ignore exceptions for now...
            }
        }
    }

private:
    std::list<callback_t> m_handlers;
};

}

#endif
