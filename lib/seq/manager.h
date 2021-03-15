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

#ifndef __NACS_SEQ_MANAGER_H__
#define __NACS_SEQ_MANAGER_H__

#include "seq.h"
#include "host_seq.h"

#include <nacs-utils/llvm/execute.h>

#include <algorithm>
#include <map>
#include <vector>

namespace NaCs::Seq {

class Manager {
    using dataset_t = std::map<std::vector<uint8_t>,uint64_t,std::less<>>;
    struct DataInfo {
        dataset_t::iterator data;
        uint32_t ref_count;
    };
    class CGContext;
public:
    struct ExpSeq;

    Manager();
    ~Manager();

    ExpSeq *create_sequence(const uint8_t *data, size_t size);
    void free_sequence(ExpSeq *seq);

private:
    uint64_t add_data(const void *data, size_t size);
    std::pair<const void*,size_t> get_data(uint64_t id) const;
    void unref_data(uint64_t id);

    llvm::LLVMContext m_llvm_ctx;
    LLVM::Exe::Engine m_engine;

    uint64_t m_data_id_cnt = 0;
    dataset_t m_dataset;
    std::map<uint64_t,DataInfo> m_datainfo;
};

}

#endif
