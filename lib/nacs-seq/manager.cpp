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
#include "error.h"

#include "../nacs-utils/llvm/codegen.h"
#include "../nacs-utils/number.h"

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
    ~CGContext() override
    {
        for (auto &entry: dataids) {
            mgr.unref_data(entry.getValue());
        }
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

Manager::ExpSeq::ExpSeq(Manager *mgr, CGContext *cgctx)
    : m_mgr(*mgr),
      m_cgctx(cgctx)
{
    assert(cgctx);
}

Manager::ExpSeq::~ExpSeq()
{
    for (auto &dataid: m_dataids)
        m_mgr.unref_data(dataid);
    m_mgr.m_engine.free(obj_id);
}

inline Manager::CGContext *Manager::ExpSeq::_cgctx() const
{
    assert(m_cgctx);
    return static_cast<CGContext*>(m_cgctx);
}

inline void Manager::ExpSeq::clear_compiler()
{
    assert(m_cgctx);
    for (auto &entry: _cgctx()->dataids)
        m_dataids.push_back(entry.getValue());
    // Clear the data ID so that we don't release it when the cgctx is destructed
    _cgctx()->dataids.clear();
    m_cgctx = nullptr;
}

NACS_EXPORT() uint64_t Manager::ExpSeq::get_dataid(llvm::StringRef name) const
{
    assert(m_cgctx);
    auto it = _cgctx()->dataids.find(name);
    if (it == _cgctx()->dataids.end())
        return 0;
    return it->second;
}

NACS_EXPORT() void Manager::ExpSeq::init_run()
{
    host_seq.init_run();
    // TODO backend
}

NACS_EXPORT() void Manager::ExpSeq::pre_run()
{
    host_seq.pre_run();
    if (max_seq_length > 0) {
        auto seq_idx = host_seq.cur_seq_idx();
        assert(seq_idx != uint32_t(-1));
        auto &host_bseq = host_seq.seqs[seq_idx];
        if (host_bseq.length > max_seq_length) {
            throw Error(Error::Type::BasicSeq, Error::BasicSeq::LengthLimit,
                        Error::Type::BasicSeq, host_bseq.id, "Sequence exceeds length limit.");
        }
    }
    // TODO backend
}

NACS_EXPORT() void Manager::ExpSeq::start()
{
    // TODO backend
    // TODO backend ordering
}

NACS_EXPORT() void Manager::ExpSeq::cancel()
{
    // TODO backend
}

NACS_EXPORT() bool Manager::ExpSeq::wait(uint64_t timeout_ms)
{
    // TODO backend
    return false;
}

NACS_EXPORT() uint32_t Manager::ExpSeq::post_run()
{
    return host_seq.post_run();
    // TODO backend
}

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
    auto expseq = std::make_unique<ExpSeq>(this, cgctx);
    expseq->max_seq_length = m_max_seq_length;
    {
        Builder builder(seq);
        builder.deserialize(data, size);
        m_backend_datas = std::move(builder.backend_datas);
        // TODO: apply noramp sequence setting from backend
        builder.buildseq();
    }
    seq.prepare();
    seq.optimize();
    cgctx->remove_unused_data(seq.env().llvm_module());
    Compiler compiler(expseq->host_seq, seq, m_engine, cgctx->get_extern_resolver());
    expseq->obj_id = compiler.compile();
    expseq->host_seq.init();
    // TODO: call backends
    expseq->clear_compiler();
    return expseq.release();
}

NACS_EXPORT() void Manager::free_sequence(Manager::ExpSeq *expseq)
{
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

NACS_EXPORT() std::pair<const void*,size_t> Manager::get_data(uint64_t id) const
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

NACS_EXPORT() void Manager::load_config_file(const char *fname)
{
    update_config(YAML::LoadFile(fname));
}

NACS_EXPORT() void Manager::load_config_string(const char *str)
{
    update_config(YAML::Load(str));
}

void Manager::update_config(const YAML::Node &config)
{
    if (auto tick_per_sec_node = config["tick_per_sec"]) {
        auto tick_per_sec = tick_per_sec_node.as<int64_t>();
        if (tick_per_sec <= 0)
            throw std::runtime_error("Invalid tick_per_sec value.");
        // Since this is used at runtime, changing this could affect sequence that are running.
        // making them inconsistent.
        // If necessary, we could keep track of all the sequences
        // and require all of them to be free'd before we can change this value.
        if (m_tick_per_sec > 0 && m_tick_per_sec != tick_per_sec)
            throw std::runtime_error("tick_per_sec value cannot be changed without restart.");
        m_tick_per_sec = tick_per_sec;
    }
    else {
        throw std::runtime_error("tick_per_sec configuration missing.");
    }
    if (auto max_length_sec_node = config["max_length_sec"]) {
        auto max_length_sec = max_length_sec_node.as<double>();
        if (max_length_sec <= 0)
            max_length_sec = 0;
        auto ticks = max_length_sec * (double)m_tick_per_sec;
        // Overflowed, which is fine.
        if (ticks >= (double)UINT64_MAX)
            ticks = 0;
        m_max_seq_length = uint64_t(ticks);
    }
    else {
        m_max_seq_length = 0;
    }
    m_device_info.clear();
    m_device_order.clear();
    if (auto devices_node = config["devices"]) {
        if (!devices_node.IsMap())
            throw std::runtime_error("Invalid value for devices.");
        for (auto node: devices_node) {
            auto dev_name = node.first.as<std::string>();
            auto dev_info = node.second;
            if (!dev_info.IsMap())
                throw std::runtime_error("Invalid device info: "
                                         "expect `backend` and `config` keys.");
            auto backend_node = dev_info["backend"];
            if (!backend_node)
                throw std::runtime_error("Invalid device info: missing `backend`.");
            m_device_info.emplace(dev_name,
                                  DeviceInfo{backend_node.as<std::string>(),
                                      dev_info["config"]});
            m_device_order.push_back(std::move(dev_name));
        }
    }
}

NACS_EXPORT() int64_t Manager::tick_per_sec() const
{
    if (m_tick_per_sec == 0)
        throw std::runtime_error("Sequence time unit not initialized.");
    return m_tick_per_sec;
}

}
