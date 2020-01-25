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

#ifndef __NACS_SEQ_RTENV_H__
#define __NACS_SEQ_RTENV_H__

#include <nacs-utils/llvm/execute.h>

#include <vector>

#include <string.h>

namespace NaCs::Seq {

/**
 * This is the runtime version of `Env`, which manages (most) information needed to run
 * code at runtime. Interpolation data may not be included in this class since it'll be
 * managed by another class in order to save memory.
 */
class RTEnv {
public:
    RTEnv(LLVM::Exe::Engine &engine,
          const std::function<uintptr_t(const std::string&)> &resolver,
          uint32_t n_rtdata, const char *code, size_t len);
    ~RTEnv();
    void copy_data(uint32_t idx, const uint64_t *data, uint32_t ndata)
    {
        memcpy(&m_rtdata[idx], data, ndata * 8);
    }
    void reset_data()
    {
        memset(m_rtdata.data(), 0, m_rtdata.size() * 8);
    }
    void *get_symbol(llvm::StringRef name) const
    {
        return m_engine.get_symbol(name);
    }

private:
    uintptr_t resolver(const std::string &name) const;

    LLVM::Exe::Engine &m_engine;
    std::function<uintptr_t(const std::string&)> m_resolver;
    // These are runtime data used by the code.
    // The data will be filled up by the user (from server).
    std::vector<uint64_t> m_rtdata;
    const uint64_t m_objid;
};

}

#endif
