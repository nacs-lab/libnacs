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
#include "device.h"
#include "error.h"

#include "../nacs-utils/llvm/codegen.h"
#include "../nacs-utils/log.h"
#include "../nacs-utils/number.h"
#include "../nacs-utils/timer.h"

#include <llvm/ADT/StringMap.h>
#include <llvm/ADT/StringSet.h>

#include <algorithm>
#include <chrono>
#ifdef __cpp_lib_execution
#  include <execution>
#endif
#include <thread>

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

template<typename Str>
Device *Manager::ExpSeq::_get_device(Str &&name, bool create)
{
    auto it = m_devices.find(name);
    if (it != m_devices.end())
        return it->second.get();
    if (!create)
        return nullptr;
    if (m_mgr.m_use_dummy_device) {
        auto dev = Device::create(m_mgr, "dummy", name);
        assert(dev);
        auto res = dev.get();
        m_devices.emplace(std::forward<Str>(name), std::move(dev));
        return res;
    }
    auto it2 = m_mgr.m_device_info.find(name);
    if (it2 == m_mgr.m_device_info.end())
        return nullptr;
    auto dev = Device::create(m_mgr, it2->second.backend, name);
    auto res = dev.get();
    if (dev) {
        dev->config(it2->second.config);
        auto data_it = m_mgr.m_backend_datas.find(name);
        if (data_it != m_mgr.m_backend_datas.end())
            dev->parse_data(data_it->second.data(), data_it->second.size());
        m_devices.emplace(std::forward<Str>(name), std::move(dev));
    }
    return res;
}

NACS_EXPORT() Device *Manager::ExpSeq::get_device(const std::string &name, bool create)
{
    return _get_device(name, create);
}

NACS_EXPORT() Device *Manager::ExpSeq::get_device(std::string &&name, bool create)
{
    return _get_device(std::move(name), create);
}

NACS_EXPORT() void Manager::ExpSeq::init_run()
{
    mgr().add_debug("Initialize sequence run\n");
    host_seq.init_run();
    for (auto &[name, dev]: m_devices) {
        dev->init_run(host_seq);
    }
}

NACS_EXPORT() void Manager::ExpSeq::pre_run()
{
    mgr().add_debug_printf("Preparing basic sequence [index %d]\n", host_seq.cur_seq_idx());
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
    // Assuming there's no complex dependencies between different `prepare_run`s
    // or that the backend can handle those themselves.
    for (auto &[name, dev]: m_devices)
        dev->prepare_run(host_seq);
#ifdef __cpp_lib_execution
    std::for_each(std::execution::par_unseq, m_devices.begin(), m_devices.end(),
                  [&] (auto &item) { item.second->pre_run(host_seq); });
#else
    std::for_each(m_devices.begin(), m_devices.end(),
                  [&] (auto &item) { item.second->pre_run(host_seq); });
#endif
}

NACS_EXPORT() void Manager::ExpSeq::start()
{
    mgr().add_debug_printf("Starting basic sequence [index %d]\n", host_seq.cur_seq_idx());
    m_cancelled.store(false, std::memory_order_relaxed);
    for (auto &name: m_mgr.m_device_order) {
        // Start in the correct order
        if (auto dev = get_device(name, false)) {
            dev->start(host_seq);
        }
    }
    m_running = true;
    std::thread wait_thread([this] {
        for (auto &[name, dev]: m_devices) {
            if (m_cancelled.load(std::memory_order_relaxed))
                break;
            dev->wait(host_seq);
        }
        {
            std::unique_lock locker(m_run_lock);
            m_running = false;
        }
        m_run_cond.notify_all();
    });
    wait_thread.detach();
}

NACS_EXPORT() void Manager::ExpSeq::cancel()
{
    mgr().add_debug_printf("Cancelling basic sequence [index %d]\n", host_seq.cur_seq_idx());
    m_cancelled.store(true, std::memory_order_relaxed);
    for (auto &[name, dev]: m_devices) {
        dev->cancel(host_seq);
    }
}

NACS_EXPORT() bool Manager::ExpSeq::wait(uint64_t timeout_ms)
{
    mgr().add_debug_printf("Waiting for basic sequence [index %d]\n", host_seq.cur_seq_idx());
    std::unique_lock locker(m_run_lock);
    return m_run_cond.wait_for(locker, std::chrono::milliseconds(timeout_ms), [&] {
        return !m_running || m_cancelled;
    });
}

NACS_EXPORT() uint32_t Manager::ExpSeq::post_run()
{
    mgr().add_debug_printf("Finishing basic sequence [index %d]\n", host_seq.cur_seq_idx());
    auto id = host_seq.post_run();
    mgr().add_debug_printf("Next basic sequence [index %d]\n", host_seq.cur_seq_idx());
    if (id == 0) {
        mgr().add_debug_printf("Sequence ended\n");
        for (auto &[name, dev]: m_devices) {
            dev->finish_run(host_seq);
        }
    }
    return id;
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
        add_debug("Deserializig sequence.\n");
        builder.deserialize(data, size);
        m_backend_datas = std::move(builder.backend_datas);
        auto nchn = builder.chnnames.size();
        for (uint32_t i = 0; i < nchn; i++) {
            const auto &name = builder.chnnames[i];
            auto chn_id = i + 1;
            auto sep = name.find_first_of('/');
            auto name_len = name.size();
            if (sep == name.npos || sep == name_len - 1)
                throw std::runtime_error("Invalid channel name.");
            add_debug_printf("Adding channel: %s\n", name.c_str());
            auto dev = expseq->get_device(name.substr(0, sep), true);
            auto chn_name = name.substr(sep + 1);
            assert(chn_name.size() == name_len - sep - 1);
            dev->add_channel(chn_id, chn_name);
            if (dev->check_noramp(chn_id, chn_name)) {
                builder.noramp_chns.push_back(chn_id);
            }
        }
        add_debug("Building sequence\n");
        builder.buildseq();
    }
    add_debug("Preparing sequence\n");
    seq.prepare();
    add_debug("Optimizing sequence\n");
    seq.optimize();
    cgctx->remove_unused_data(seq.env().llvm_module());
    Compiler compiler(expseq->host_seq, seq, m_engine, cgctx->get_extern_resolver());
    add_debug("Compiling sequence\n");
    expseq->obj_id = compiler.compile();
    add_debug("Initializing runtime sequence\n");
    expseq->host_seq.init();
    // Note that backends may call `get_device` for other backends in `prepare`
    // and the new device may not have `prepare` called on it.
    // This shouldn't be a problem right now since we don't have very deep dependencies.
    // I'd like to wait until we know what actually need before implementing a proper solution.
    // (e.g. maybe record if something is prepared and call prepare on it,
    // or re-prepare if one device triggers update on another one)
    add_debug("Preparing backend drivers\n");
    for (auto &[name, dev]: expseq->devices())
        dev->prepare(*expseq, compiler);
    add_debug("Generating backend sequences\n");
    for (auto &[name, dev]: expseq->devices())
        dev->generate(*expseq, compiler);
    add_debug("Cleaning up compiler\n");
    expseq->clear_compiler();
    return expseq.release();
}

NACS_EXPORT() void Manager::free_sequence(Manager::ExpSeq *expseq)
{
    add_debug("Freeing sequence\n");
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
    add_debug_printf("Loading config file %s\n", fname);
    update_config(YAML::LoadFile(fname));
}

NACS_EXPORT() void Manager::load_config_string(const char *str)
{
    add_debug_printf("Loading config string\n");
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
    m_use_dummy_device = false;
    if (auto dummy_devices_node = config["use_dummy_device"];
        dummy_devices_node && dummy_devices_node.as<bool>()) {
        m_use_dummy_device = true;
    }
    else if (auto devices_node = config["devices"]) {
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

// Error+warning API for python/matlab
NACS_EXPORT() uint8_t *Manager::take_messages(size_t *sz)
{
    std::lock_guard<std::mutex> locker(m_message_lock);
    return (uint8_t*)m_messages.get_buf(*sz);
}

// Error+warning API for C++
NACS_EXPORT() void Manager::add_info(const char *msg)
{
    std::lock_guard<std::mutex> locker(m_message_lock);
    m_messages.put((char)MsgType::Info);
    uint64_t t = getTime();
    m_messages.write((const char*)&t, sizeof(t));
    uint32_t len = (uint32_t)strlen(msg);
    m_messages.write((const char*)&len, sizeof(uint32_t));
    m_messages.write(msg, len);
}

NACS_EXPORT() void Manager::add_warning(const char *msg)
{
    std::lock_guard<std::mutex> locker(m_message_lock);
    m_messages.put((char)MsgType::Warn);
    uint64_t t = getTime();
    m_messages.write((const char*)&t, sizeof(t));
    uint32_t len = (uint32_t)strlen(msg);
    m_messages.write((const char*)&len, sizeof(uint32_t));
    m_messages.write(msg, len);
}

NACS_EXPORT() void Manager::add_error(const char *msg)
{
    std::lock_guard<std::mutex> locker(m_message_lock);
    m_messages.put((char)MsgType::Error);
    uint64_t t = getTime();
    m_messages.write((const char*)&t, sizeof(t));
    uint32_t len = (uint32_t)strlen(msg);
    // Automatically delete a new line from error messages
    // to make `Log::error` and throwing an error more consistent.
    if (len > 0 && msg[len - 1] == '\n')
        len -= 1;
    m_messages.write((const char*)&len, sizeof(uint32_t));
    m_messages.write(msg, len);
}

NACS_EXPORT() void Manager::add_error(const std::exception &error)
{
    add_error(error.what());
}

NACS_EXPORT() void Manager::add_error(const Error &error)
{
    std::lock_guard<std::mutex> locker(m_message_lock);
    m_messages.put((char)MsgType::SeqError);
    uint64_t t = getTime();
    m_messages.write((const char*)&t, sizeof(t));
    write_bits(m_messages, error.type);
    write_bits(m_messages, error.code);
    write_bits(m_messages, error.type1);
    write_bits(m_messages, error.id1);
    write_bits(m_messages, error.type2);
    write_bits(m_messages, error.id2);
    auto msg = error.what();
    uint32_t len = (uint32_t)strlen(msg);
    m_messages.write((const char*)&len, sizeof(uint32_t));
    m_messages.write(msg, len);
}

NACS_EXPORT() void Manager::_add_debug(const char *msg)
{
    std::lock_guard<std::mutex> locker(m_message_lock);
    m_messages.put((char)MsgType::Debug);
    uint64_t t = getTime();
    m_messages.write((const char*)&t, sizeof(t));
    uint32_t len = (uint32_t)strlen(msg);
    m_messages.write((const char*)&len, sizeof(uint32_t));
    m_messages.write(msg, len);
}

NACS_EXPORT() void Manager::_add_debug_printf(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
#if NACS_OS_LINUX
    char *str = nullptr;
    if (vasprintf(&str, fmt, ap) == -1) {
        add_error("Unable to allocate memory for debug message\n");
        va_end(ap);
        return;
    }
#else
    va_list aq;
    va_copy(aq, ap);
    auto size = vsnprintf(nullptr, 0, fmt, aq);
    va_end(aq);
    // size doesn't include the NUL byte at the end.
    char *str = (char*)malloc(size + 1);
    auto size2 = vsnprintf(str, size + 1, fmt, ap);
    assert(size == size2);
    (void)size2;
#endif
    va_end(ap);
    _add_debug(str);
    free(str);
}

NACS_EXPORT() void Manager::enable_debug(bool enable)
{
    if (m_debug_messages == enable)
        return;
    m_debug_messages = enable;
    if (m_debug_messages) {
        Log::level = Log::Debug;
    }
    else {
        Log::level = Log::Info;
    }
}

NACS_EXPORT() void Manager::register_logger()
{
    Log::pushLogger([&] (Log::Level level, const char *, const char *msg) {
        if (level == Log::Info || level == Log::Force) {
            add_info(msg);
        }
        else if (level == Log::Warn) {
            add_warning(msg);
        }
        else if (level == Log::Error) {
            add_error(msg);
        }
        else if (level == Log::Debug) {
            add_debug(msg);
        }
    });
}

NACS_EXPORT() void Manager::unregister_logger()
{
    Log::popLogger();
}

}
