//

#ifndef __NACS_SEQ_AWG_BACKEND_H__
#define __NACS_SEQ_AWG_BACKEND_H__

#include "../device.h"

#include <map>
#include <vector>

namespace NaCs::Seq {

namespace AWG {

class Backend : public Device {
    struct BasicSeq;
    struct ChannelInfo;

public:
    Backend(Manager &mgr, std::string name);
    ~Backend() override;
    const uint8_t *get_channel_info(uint32_t *sz);

private:
    void add_channel(uint32_t chn_id, const std::string &chn_name) override;
    void prepare (Manager::ExpSeq &expseq, Compiler &compiler) override;
    void generate (Manager::ExpSeq &expseq, Compiler &compiler) override;

    void prepare_run(HostSeq &host_seq) override;
    void pre_run(HostSeq &host_seq) override;
    void start(HostSeq &host_seq) override;
    void cancel(HostSeq &host_seq) override;
    void wait(HostSeq &host_seq) override;

    // void config(const YAML::Node&) override;

    void sort_channels();
    // Return whether pregeneration competed
    bool pregenerate(Manager::ExpSeq &expseq, Compiler &compiler, uint32_t seq_idx);
    void populate_values(HostSeq &host_seq, BasicSeq &bseq);

    std::vector<BasicSeq> m_seqs;
    // Map from sequence channel ID to physical device and channel name/ids etc.
    std::map<uint32_t,ChannelInfo> m_chn_map;
    // From linear channel ID within the backend to more detailed channel info
    std::vector<std::map<uint32_t, ChannelInfo>::iterator> m_linear_chns;
    uint32_t m_step_Size;
    uint32_t m_nconst_vals;
    bool m_first_bseq = false;
    uint32_t m_cur_seq_id = 0;
    uint64_t m_obj_id = 0;
    // Serialized channel info to make passing it to python/MATLAB easier
    std::vector<uint8_t> m_channel_info;

};

}

}

#endif
