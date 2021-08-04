//

#ifndef __NACS_SEQ_AWG_BACKEND_H__
#define __NACS_SEQ_AWG_BACKEND_H__

#include "../device.h"

#include <map>
#include <vector>
#include <memory>

#include <nacs-utils/zmq_utils.h>


namespace NaCs::Seq {

namespace Zynq {
class Backend;
}

namespace AWG {

//globals for all AWG backends
uint64_t seqcount = 0;

struct GlobalInfo {
    uint64_t client_id;
    uint64_t server_id;
    std::string triple_str;
    std::string cpu_str;
    std::string feature_str;
    
    GlobalInfo(uint64_t client, uint64_t server, std::string triple,
               std::string cpu, std::string feature):
        client_id(client),
        server_id(server),
        triple_str(triple),
        cpu_str(cpu),
        feature_str(feature)
    {
    }
};

std::map <std::string, GlobalInfo> info_map;

class Backend : public Device {
    enum class ChnType : uint8_t;
    struct BasicSeq;
    struct ChannelInfo;

public:
    Backend(Manager &mgr, std::string name);
    ~Backend() override;
    const uint8_t *get_channel_info(uint32_t *sz);

    // testing
    static Backend *cast(Device *dev);
    static const Backend *cast(const Device *dev);
    //llvm::ArrayRef<uint8_t> get_bytecode(uint32_t bseq_id);
private:
    void add_channel(uint32_t chn_id, const std::string &chn_name) override;
    void prepare (Manager::ExpSeq &expseq, Compiler &compiler) override;
    void generate (Manager::ExpSeq &expseq, Compiler &compiler) override;

    void init_run(HostSeq &host_seq) override;
    void prepare_run(HostSeq &host_seq) override;
    void pre_run(HostSeq &host_seq) override;
    void start(HostSeq &host_seq) override;
    void cancel(HostSeq &host_seq) override;
    void wait(HostSeq &host_seq) override;
    void finish_run(HostSeq &host_seq) override;
    // void config(const YAML::Node&) override;

    void sort_channels();
    // Return whether pregeneration competed
    bool pregenerate(Manager::ExpSeq &expseq, Compiler &compiler, uint32_t seq_idx);
    void populate_values(HostSeq &host_seq, BasicSeq &bseq);
    void reqServerInfo();

    uint8_t get_pulse_type(Backend::ChnType type, bool is_fn, bool is_vector);
    void get_chn_from_lin(uint32_t lin_idx, Backend::ChnType &type, uint8_t &phys_chn_id, uint32_t &chn_id);
    void config (const YAML::Node&) override;
    std::vector<BasicSeq> m_seqs;
    // Map from sequence channel ID to physical device and channel name/ids etc.
    std::map<uint32_t,ChannelInfo> m_chn_map;
    // From linear channel ID within the backend to more detailed channel info
    std::vector<std::map<uint32_t, ChannelInfo>::iterator> m_linear_chns;
    uint32_t m_step_Size;
    uint32_t m_nconst_vals;
    bool m_first_bseq = false;
    uint32_t m_cur_seq_id = 0;
    uint64_t m_cur_wait_id = 0;
    uint64_t m_obj_id = 0;
    uint64_t m_seq_cnt_offset = 0;
    // Serialized channel info to make passing it to python/MATLAB easier
    std::vector<uint8_t> m_channel_info;
    std::string m_clock_dev;
    Zynq::Backend *m_zynq_dev = nullptr;
    // network related variables...
    std::string m_url;
    std::unique_ptr<ZMQ::MultiClient::SockRef> m_sock;
    uint64_t m_client_id;
    uint64_t m_server_id;
    std::string m_triple_str;
    std::string m_cpu_str;
    std::string m_feature_str;
    std::vector<uint8_t> id_bc;

    std::vector<uint8_t> out_chns;
};

}

}

#endif
