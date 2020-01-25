/*************************************************************************
 *   Copyright (c) 2020 - 2021 Yichao Yu <yyc1992@gmail.com>             *
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

#include "rtenv.h"

namespace NaCs::Seq {

NACS_EXPORT_ RTEnv::RTEnv(LLVM::Exe::Engine &engine,
                          const std::function<uintptr_t(const std::string&)> &resolver,
                          uint32_t n_rtdata, const char *code, size_t len)
    : m_engine(engine),
      m_resolver(std::move(resolver)),
      m_rtdata(n_rtdata, 0),
      m_objid(engine.load(code, len, std::bind(&RTEnv::resolver, this, std::placeholders::_1)))
{
    if (!m_objid) {
        throw std::runtime_error("Unable to load object.");
    }
}

NACS_EXPORT() RTEnv::~RTEnv()
{
    if (m_objid) {
        m_engine.free(m_objid);
    }
}

uintptr_t RTEnv::resolver(const std::string &name) const
{
    if (name == "rt.data")
        return (uintptr_t)m_rtdata.data();
    return m_resolver(name);
}

}
