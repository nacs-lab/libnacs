//

#ifndef __NACS_SEQ_AWG_DUMMYCLIENT_H__
#define __NACS_SEQ_AWG_DUMMYCLIENT_H__

// include "../device.h"

#include "../../utils/llvm/utils.h"
#include "../../utils/llvm/codegen.h"
#include "../../utils/llvm/compile.h"
#include "../../utils/llvm/execute.h"
#include "../../utils/llvm/global_rename.h"
#include "../../utils/processor.h"
#include "../../utils/streams.h"

#include <llvm/ADT/StringRef.h>
#include <llvm/IR/LegacyPassManager.h>
#include <llvm/IR/Verifier.h>
#include <llvm/Transforms/IPO.h>

#include <map>
#include <vector>
#include <memory>

#include <nacs-utils/zmq_utils.h>
#include <yaml-cpp/yaml.h>

namespace NaCs::Seq::AWG::Dummy {

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

class DummyClient {
    enum class ChnType : uint8_t;
    struct BasicSeq;
    struct ChannelInfo;

public:
    DummyClient();
    //{
    //    std::cout << "calling constructor" << std::endl;
    //    m_seq_cnt_offset = seqcount; //seqcount is a global
    //}
    ~DummyClient();
    const uint8_t *get_channel_info(uint32_t *sz);

public:
    void add_channel(uint32_t chn_id, const std::string &chn_name);
    //void prepare (Manager::ExpSeq &expseq, Compiler &compiler);
    //void generate (Manager::ExpSeq &expseq, Compiler &compiler);

    void init_run() ;
    //void prepare_run(HostSeq &host_seq) ;
    //void pre_run(HostSeq &host_seq) ;
    void start() ;
    //void cancel(HostSeq &host_seq);
    void wait() ;
    //void finish_run(HostSeq &host_seq);
    // void config(const YAML::Node&);

    void sort_channels();
    // Return whether pregeneration competed
    //bool pregenerate(Manager::ExpSeq &expseq, Compiler &compiler, uint32_t seq_idx);
    //void populate_values(HostSeq &host_seq, BasicSeq &bseq);


    uint8_t get_pulse_type(DummyClient::ChnType type, bool is_fn, bool is_vector);
    void get_chn_from_lin(uint32_t lin_idx, DummyClient::ChnType &type, uint8_t &phys_chn_id, uint32_t &chn_id);
    //void config (const YAML::Node&);
    void loadConfig (const char *fname);
    // std::vector<BasicSeq> m_seqs;
    // Map from sequence channel ID to physical device and channel name/ids etc.
    std::map<uint32_t,ChannelInfo> m_chn_map;
    // From linear channel ID within the backend to more detailed channel info
    std::vector<std::map<uint32_t, ChannelInfo>::iterator> m_linear_chns;
    //uint32_t m_step_Size;
    // uint32_t m_nconst_vals;
    //bool m_first_bseq = false;
    //uint32_t m_cur_seq_id = 0;
    uint64_t m_cur_wait_id = 0;
    //uint64_t m_obj_id = 0;
    uint64_t m_seq_cnt_offset = 0;
    // Serialized channel info to make passing it to python/MATLAB easier
    //std::vector<uint8_t> m_channel_info;
    //std::string m_clock_dev;
    //Zynq::Backend *m_zynq_dev = nullptr;
    // network related variables...
    std::unique_ptr<ZMQ::MultiClient::SockRef> m_sock;
public:
    std::string m_url;
    uint64_t m_client_id;
    uint64_t m_server_id;
    std::string m_triple_str;
    std::string m_cpu_str;
    std::string m_feature_str;
    std::vector<uint8_t> id_bc;
};


}
#endif
