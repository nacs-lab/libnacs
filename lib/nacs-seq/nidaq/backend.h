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

#ifndef __NACS_SEQ_NIDAQ_BACKEND_H__
#define __NACS_SEQ_NIDAQ_BACKEND_H__

#include "../device.h"

#include <map>
#include <vector>

namespace NaCs::Seq {

namespace Zynq {
class Backend;
}

namespace NiDAQ {

class Backend : public Device {
    struct BasicSeq;
    struct ChannelInfo;

public:
    Backend(Manager &mgr, std::string name);
    ~Backend() override;
    // Can only be called after `pre_run`
    const double *get_data(uint32_t cur_seq_id, size_t *nsamples) const;
    // Can only be called after `prepare`
    const std::pair<uint32_t,const char*> *get_channel_info(uint32_t *sz);

    // For testing only
    static Backend *cast(Device *dev);
    static const Backend *cast(const Device *dev);
    struct GenerateStatus {
        std::vector<size_t> fill_init{};
        size_t nsamples = size_t(-1); // size_t(-1) if no static generation.
        bool first_only = false;
        bool all_time_known = false;
        bool all_val_known = false;
    };
    GenerateStatus get_generate_status(uint32_t cur_seq_id) const;

private:
    void add_channel(uint32_t chn_id, const std::string &chn_name) override;
    // bool check_noramp(uint32_t chn_id, const std::string &chn_name) override;
    void prepare(Manager::ExpSeq &expseq, Compiler &compiler) override;
    void generate(Manager::ExpSeq &expseq, Compiler &compiler) override;

    // void init_run(HostSeq &host_seq) override;
    void prepare_run(HostSeq &host_seq) override;
    void pre_run(HostSeq &host_seq) override;
    void start(HostSeq &host_seq) override;
    void cancel(HostSeq &host_seq) override;
    void wait(HostSeq &host_seq) override;
    // void post_run(HostSeq &host_seq) override;

    void config(const YAML::Node&) override;
    // void parse_data(const uint8_t *data, size_t len) override;

    void sort_channels();
    // Return whether pregeneration completed
    bool pregenerate(Manager::ExpSeq &expseq, Compiler &compiler, uint32_t seq_idx);
    void populate_values(HostSeq &host_seq, BasicSeq &bseq);

    std::vector<BasicSeq> m_seqs;
    // Map from sequence channel ID to physical device and channel name/ids etc.
    std::map<uint32_t,ChannelInfo> m_chn_map;
    // From linear channel ID within the backend to more detailed channel info.
    std::vector<std::map<uint32_t,ChannelInfo>::iterator> m_linear_chns;
    uint32_t m_step_size; // How many sequence clock cycle per refresh step
    uint32_t m_nconst_vals;
    bool m_first_bseq = false;
    uint64_t m_obj_id = 0;
    std::string m_clock_dev;
    Zynq::Backend *m_zynq_dev = nullptr;
    // Channel info for python/MATLAB.
    std::vector<std::pair<uint32_t,const char*>> m_channel_info;
};

}
}

#endif
