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

#include "manager.h"

#include "builder.h"
#include "compiler.h"

#include "../utils/llvm/codegen.h"

#include <llvm/ADT/StringMap.h>
#include <llvm/ADT/StringSet.h>

namespace NaCs::Seq {

class Manager::CGContext: public LLVM::Codegen::CachedContext {
public:
    CGContext(Manager &mgr)
        : LLVM::Codegen::CachedContext(mgr.m_llvm_ctx),
          mgr(mgr)
    {
    }

    bool use_extern_data() const override
    {
        return true;
    }
    void add_extern_data(llvm::StringRef name, const void *data, size_t size) override
    {
        dataids[name] = mgr.add_data(data, size);
    }
    uintptr_t get_extern_data(llvm::StringRef name) override
    {
        auto it = dataids.find(name);
        if (it == dataids.end())
            return 0;
        return (uintptr_t)mgr.get_data(it->second).first;
    }

    void remove_unused_data(llvm::Module *mod)
    {
        llvm::StringSet to_delete;
        for (auto &entry: dataids) {
            auto key = entry.getKey();
            if (!mod->getGlobalVariable(key)) {
                to_delete.insert(key);
            }
        }
        for (auto &entry: to_delete) {
            auto key = entry.getKey();
            auto it = dataids.find(key);
            mgr.unref_data(it->getValue());
            dataids.erase(it);
        }
    }

    Manager &mgr;
    llvm::StringMap<uint64_t> dataids;
};
struct Manager::ExpSeq {
    HostSeq host_seq;
    uint64_t obj_id;
    // TODO
};
NACS_EXPORT_ Manager::Manager()
{
}

NACS_EXPORT() Manager::~Manager()
{
}

NACS_EXPORT() Manager::ExpSeq *Manager::create_sequence(const uint8_t *data, size_t size)
{
    auto cgctx = new CGContext(*this);
    Seq seq{std::unique_ptr<CGContext>(cgctx)};
    {
        Builder builder(seq);
        builder.deserialize(data, size);
        // TODO? apply noramp sequence config
        builder.buildseq();
    }
    seq.prepare();
    seq.optimize();
    cgctx->remove_unused_data(seq.env().llvm_module());
    // TODO memory management?
    auto expseq = std::make_unique<ExpSeq>();
    Compiler compiler(expseq->host_seq, seq, m_engine, cgctx->get_extern_resolver());
    expseq->obj_id = compiler.compile();
    return expseq.release();
}

NACS_EXPORT() void Manager::free_sequence(Manager::ExpSeq *expseq)
{
    m_engine.free(expseq->obj_id);
    delete expseq;
}

uint64_t Manager::add_data(const void *_data, size_t size)
{
    auto data = (const uint8_t*)_data;
    ArrayKey<uint8_t> key{data, size};
    auto it = m_dataset.find(key);
    uint64_t id;
    if (it != m_dataset.end()) {
        id = it->second;
        m_datainfo.find(id)->second.ref_count++;
        return id;
    }
    id = ++m_data_id_cnt;
    it = m_dataset.emplace(std::vector<uint8_t>(data, data + size), id).first;
    m_datainfo.emplace(id, DataInfo{it, 1});
    return id;
}

std::pair<const void*,size_t> Manager::get_data(uint64_t id) const
{
    auto it = m_datainfo.find(id);
    if (it == m_datainfo.end())
        return {nullptr, 0};
    auto &data = it->second.data->first;
    return {data.data(), data.size()};
}

void Manager::unref_data(uint64_t id)
{
    auto it = m_datainfo.find(id);
    assert(it != m_datainfo.end());
    auto &info = it->second;
    if (--info.ref_count > 0)
        return;
    m_dataset.erase(info.data);
    m_datainfo.erase(it);
}

}
