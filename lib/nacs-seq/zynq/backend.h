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

#ifndef __NACS_SEQ_ZYNQ_BACKEND_H__
#define __NACS_SEQ_ZYNQ_BACKEND_H__

#include "../device.h"

#include "bc_gen.h"

#include <nacs-utils/zmq_utils.h>

#include <map>
#include <memory>

namespace NaCs::Seq::Zynq {

class Backend : public Device {
    struct BasicSeq;
    struct TTLManager;
public:
    Backend(Manager &mgr, std::string name);
    ~Backend() override;

    void enable_clock()
    {
        m_has_clock = true;
    }

    void set_clock_active_time(uint32_t bseq_idx, uint64_t half_period,
                               llvm::ArrayRef<std::pair<int64_t,int64_t>> active_times);

    // For testing only
    static Backend *cast(Device *dev);
    static const Backend *cast(const Device *dev);
    bool has_generator(uint32_t bseq_idx) const;
    bool pregenerated(uint32_t bseq_idx, bool first_bseq) const;
    llvm::ArrayRef<uint8_t> get_bytecode(uint32_t bseq_idx, bool first_bseq) const;
    llvm::ArrayRef<BCGen::Clock> get_clock(uint32_t bseq_idx) const;

private:
    void add_channel(uint32_t chn_id, const std::string &chn_name) override;
    bool check_noramp(uint32_t chn_id, const std::string &chn_name) override;
    // void prepare(Manager::ExpSeq &expseq, Compiler &compiler) override;
    void generate(Manager::ExpSeq &expseq, Compiler &compiler) override;

    // void init_run(HostSeq &host_seq) override;
    // void prepare_run(HostSeq &host_seq) override;
    void pre_run(HostSeq &host_seq) override;
    void start(HostSeq &host_seq) override;
    void cancel(HostSeq &host_seq) override;
    void wait(HostSeq &host_seq) override;
    // void finish_run(HostSeq &host_seq) override;

    void run_bytecode(const std::vector<uint8_t> &bc);
    // Return whether pregeneration completed
    bool pregenerate(Manager::ExpSeq &expseq, Compiler &compiler, uint32_t seq_idx);

    void config(const YAML::Node&) override;
    void parse_data(const uint8_t *data, size_t len) override;

    std::vector<BasicSeq> m_seqs;
    std::map<uint32_t,std::pair<BCGen::ChnType,uint8_t>> m_chn_map;
    std::map<uint32_t,std::vector<BCGen::Clock>> m_clocks;
    uint32_t m_ttl_mask = 0;
    uint32_t m_seq_delay; // in unit of 10ns.
    int8_t m_start_ttl_chn;
    bool m_has_clock = false;
    bool m_generated = false;
    uint64_t m_seq_tick_per_sec;
    uint64_t m_obj_id = 0;
    std::string m_url;
    std::unique_ptr<ZMQ::MultiClient::SockRef> m_sock;
    uint64_t m_seq_id[2];
    // The keys in this `m_ttl_managers` are the sequence channel ID
    // since we don't know about all the channel mapping yet when we populate this.
    std::map<uint32_t,TTLManager> m_ttl_managers;

    uint32_t m_ttl_ovr_ignore = 0;
    std::set<uint8_t> m_dds_ovr_ignore;
};

}

#endif
