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
#include "error.h"

#include <nacs-utils/llvm/execute.h>
#include <nacs-utils/signals.h>
#include <nacs-utils/streams.h>

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <any>
#include <atomic>
#include <condition_variable>
#include <map>
#include <memory>
#include <mutex>
#include <vector>

namespace NaCs::LLVM::Codegen {
class CachedContext;
}

namespace NaCs::Seq {

class Device;

class Manager {
    using dataset_t = std::map<std::vector<uint8_t>,uint64_t,std::less<>>;
    struct DataInfo {
        dataset_t::iterator data;
        uint32_t ref_count;
    };
    class CGContext;
    struct DeviceInfo {
        std::string backend;
        YAML::Node config;
        uint32_t n_restarts = 0;
        bool marked; // flag for mark-and-sweep when updating config
    };
    struct LoggerRegister {
        LoggerRegister(Manager*);
        ~LoggerRegister();
    private:
        Manager &m_mgr;
    };

public:
    enum class MsgType : uint8_t {
        Info = 0,
        Warn,
        Error,
        SeqError,
        Debug,
    };

    class ExpSeq {
    public:
        struct ChnOutput {
            ChnOutput(char *name, size_t name_sz, int64_t *times, double *values, uint32_t *pulse_ids, size_t npts)
            : name(name), name_sz(name_sz), times(times), values(values), pulse_ids(pulse_ids), npts(npts)
            {
                //
            }

            char *name;
            size_t name_sz;
            int64_t *times;
            double *values;
            uint32_t *pulse_ids;
            size_t npts;
        };
        ExpSeq(Manager *mgr, CGContext *cgctx);
        ~ExpSeq();

        HostSeq host_seq;
        uint64_t obj_id = 0;
        uint64_t max_seq_length;

        Manager &mgr()
        {
            return m_mgr;
        }

        void clear_compiler();

        // For backend during compilation time only.
        LLVM::Codegen::CachedContext *cgctx() const
        {
            return m_cgctx;
        }
        uint64_t get_dataid(llvm::StringRef name) const;
        Device *get_device(const std::string &name, bool create);
        Device *get_device(std::string &&name, bool create);
        const std::map<std::string,std::unique_ptr<Device>> &devices()
        {
            return m_devices;
        }

        void init_run();
        void pre_run();
        void start();
        void cancel();
        bool wait(uint64_t timeout_ms);
        uint32_t post_run();

        const uint8_t *get_builder_dump(size_t *sz) const
        {
            if (!dump)
                return nullptr;
            *sz = dump->builder.size();
            return dump->builder.data();
        }
        const uint8_t *get_seq_dump(size_t *sz) const
        {
            if (!dump)
                return nullptr;
            *sz = dump->seq.size();
            return dump->seq.data();
        }
        const uint8_t *get_seq_opt_dump(size_t *sz) const
        {
            if (!dump)
                return nullptr;
            *sz = dump->seq_opt.size();
            return dump->seq_opt.data();
        }
        struct Dump {
            std::vector<uint8_t> builder;
            std::vector<uint8_t> seq;
            std::vector<uint8_t> seq_opt;
        };
        std::unique_ptr<Dump> dump;
        std::map<std::string,uint32_t> m_chn_map; // Tracks the mapping from channel name to index.

        uint32_t refresh_device_restart(const char *dname);
        ChnOutput* get_nominal_output(uint64_t pts_per_ramp, size_t *out_sz);

    private:
        template<typename Str>
        Device *_get_device(Str &&name, bool create);
        CGContext *_cgctx() const;

        Manager &m_mgr;
        LLVM::Codegen::CachedContext *m_cgctx;
        std::vector<uint64_t> m_dataids;
        std::map<std::string,std::unique_ptr<Device>> m_devices;

        bool m_running = false;
        std::atomic<bool> m_cancelled = false;
        std::mutex m_run_lock;
        std::condition_variable m_run_cond;
    };

    template<typename T>
    struct Storage {
        Storage() = default;
        Storage(const Storage&) = delete;
        Storage(Storage&&) = delete;
        T &get(Manager &mgr) const
        {
            return mgr.get_storage(this);
        }
    };

    Manager();
    ~Manager();

    ExpSeq *create_sequence(const uint8_t *data, size_t size);
    void free_sequence(ExpSeq *seq);
    std::pair<const void*,size_t> get_data(uint64_t id) const;
    void load_config_file(const char *fname);
    void load_config_string(const char *str);
    char *get_config_str(const std::string &name, size_t *sz) const;
    int64_t tick_per_sec() const;
    uint64_t max_seq_len() const;
    LLVM::Exe::Engine &exe_engine()
    {
        return m_engine;
    }

    // Log device restart API
    bool add_device_restart(const std::string &name, uint32_t n = 1);
    uint32_t get_device_restart(const char *dname);
    bool set_device_restart(const std::string &name, uint32_t n = 0);

    // Error+warning API for python/matlab
    uint8_t *take_messages(size_t *sz);

    // Error+warning API for C++
    void add_info(const char *msg);
    void add_warning(const char *msg);
    void add_error(const char *msg);
    void add_error(const std::exception &error);
    void add_error(const Error &error);
    void add_debug(const char *msg)
    {
        if (likely(!m_debug_messages))
            return;
        _add_debug(msg);
    }
    // This version allows the formatting to be skipped when debugging is off.
    template<typename... Args>
    void add_debug_printf(const char *fmt, Args&&... args)
    {
        if (likely(!m_debug_messages))
            return;
        _add_debug_printf(fmt, std::forward<Args>(args)...);
    }
    void enable_debug(bool enable);
    bool debug_enabled()
    {
        return m_debug_messages;
    }
    void register_logger();
    void unregister_logger();

    void add_info(const std::string &msg)
    {
        add_info(msg.c_str());
    }
    void add_warning(const std::string &msg)
    {
        add_warning(msg.c_str());
    }
    void add_error(const std::string &msg)
    {
        add_error(msg.c_str());
    }
    void add_debug(const std::string &msg)
    {
        if (likely(!m_debug_messages))
            return;
        _add_debug(msg.c_str());
    }

    template<typename F>
    void call_guarded(F &&func)
    {
        LoggerRegister reg(this);
        try {
            func();
        }
        catch (const Error &error) {
            add_error(error);
        }
        catch (const std::exception &error) {
            add_error(error);
        }
        catch (...) {
            add_error("Unknown error");
        }
    }

    template<typename F, typename T>
    auto call_guarded(F &&func, T def) -> decltype(func())
    {
        LoggerRegister reg(this);
        try {
            return func();
        }
        catch (const Error &error) {
            add_error(error);
        }
        catch (const std::exception &error) {
            add_error(error);
        }
        catch (...) {
            add_error("Unknown error");
        }
        return def;
    }

    void enable_dump(bool enable)
    {
        m_dump_sequences = enable;
    }
    bool dump_enabled()
    {
        return m_dump_sequences;
    }

    // Not affecting the manager itself very much
    // but let the backend know that it should not send any data for real execution
    bool dummy_mode = false;

    static Signal<Manager*> &new_run();

private:
    uint64_t add_data(const void *data, size_t size);
    void unref_data(uint64_t id);
    void update_config(const YAML::Node &config);
    void _add_debug(const char *msg);
    void _add_debug_printf(const char *fmt, ...);
    template<typename T>
    T &get_storage(const Storage<T> *key)
    {
        auto &val = m_storage[key];
        if (!val.has_value())
            val = T();
        return *std::any_cast<T>(&val);
    }

    llvm::LLVMContext m_llvm_ctx;
    LLVM::Exe::Engine m_engine;

    uint64_t m_data_id_cnt = 0;
    dataset_t m_dataset;
    int64_t m_tick_per_sec = 0;
    uint64_t m_max_seq_length = 0;
    bool m_use_dummy_device = false;
    std::map<std::string,DeviceInfo> m_device_info;
    std::vector<std::string> m_device_order;
    std::map<uint64_t,DataInfo> m_datainfo;
    std::map<std::string,std::vector<uint8_t>> m_backend_datas;
    std::map<const void*,std::any> m_storage;

    malloc_ostream m_messages;
    std::mutex m_message_lock;
    // Whether debug messages should be logged.
    bool m_debug_messages = false;
    // Whether we should save a dump of intermediate results for debugging
    bool m_dump_sequences = false;
};

inline Manager::LoggerRegister::LoggerRegister(Manager *mgr)
    : m_mgr(*mgr)
{
    m_mgr.register_logger();
}

inline Manager::LoggerRegister::~LoggerRegister()
{
    m_mgr.unregister_logger();
}

}

#endif
