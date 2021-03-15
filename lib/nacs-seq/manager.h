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

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <map>
#include <vector>

namespace NaCs::LLVM::Codegen {
class CachedContext;
}

namespace NaCs::Seq {

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
    };

public:
    class ExpSeq {
    public:
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

        void init_run();
        void pre_run();
        void start();
        void cancel();
        bool wait(uint64_t timeout_ms);
        uint32_t post_run();

    private:
        CGContext *_cgctx() const;

        Manager &m_mgr;
        LLVM::Codegen::CachedContext *m_cgctx;
        std::vector<uint64_t> m_dataids;
        // TODO backends
    };

    Manager();
    ~Manager();

    ExpSeq *create_sequence(const uint8_t *data, size_t size);
    void free_sequence(ExpSeq *seq);
    std::pair<const void*,size_t> get_data(uint64_t id) const;
    void load_config_file(const char *fname);
    void load_config_string(const char *str);
    int64_t tick_per_sec() const;
    LLVM::Exe::Engine &exe_engine()
    {
        return m_engine;
    }

    // Not affecting the manager itself very much
    // but let the backend know that it should not send any data for real execution
    bool dummy_mode = false;

private:
    uint64_t add_data(const void *data, size_t size);
    void unref_data(uint64_t id);
    void update_config(const YAML::Node &config);

    llvm::LLVMContext m_llvm_ctx;
    LLVM::Exe::Engine m_engine;

    uint64_t m_data_id_cnt = 0;
    dataset_t m_dataset;
    int64_t m_tick_per_sec = 0;
    uint64_t m_max_seq_length = 0;
    std::map<std::string,DeviceInfo> m_device_info;
    std::vector<std::string> m_device_order;
    std::map<uint64_t,DataInfo> m_datainfo;
    std::map<std::string,std::vector<uint8_t>> m_backend_datas;
};

}

#endif
