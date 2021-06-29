//

#ifndef __NACS_SEQ_AWG_DUMMYCLIENT_H__
#define __NACS_SEQ_AWG_DUMMYCLIENT_H__

#include "../device.h"

#include "../../nacs-utils/llvm/utils.h"
#include "../../nacs-utils/llvm/codegen.h"
#include "../../nacs-utils/llvm/compile.h"
#include "../../nacs-utils/llvm/execute.h"
#include "../../nacs-utils/processor.h"
#include "../../nacs-utils/streams.h"

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
public:
    enum class ChnType : uint8_t {
        Freq,
        Amp,
        Phase
    };
    struct ChannelInfo {
        uint8_t m_phys_chn;
        uint32_t m_chn_num;
        uint32_t linear_idx = 0; // Linear channel number for BasicSeq's list of pulses
        ChnType m_chn_type;
        ChannelInfo(uint8_t phys_chn, uint32_t chn_num, ChnType chn_type)
            : m_phys_chn(phys_chn),
              m_chn_num(chn_num),
              m_chn_type(chn_type)
        {
        }
    };
    struct BasicSeq{
        struct Pulse {
            Pulse() = default;
            Pulse(uint8_t pulse_type, uint32_t id, uint32_t time_id,
                  uint32_t len_id, uint32_t endvalue_id, uint32_t cond_id,
                  std::string ramp_func_name)
                : pulse_type(pulse_type), id(id), time_id(time_id),
                  len_id(len_id), endvalue_id(endvalue_id),
                  cond_id(cond_id), ramp_func_name(ramp_func_name)
            {}
            uint8_t pulse_type;
            uint32_t id;
            uint32_t time_id;
            uint32_t len_id; // may be -1, `len == -1` <=> no-ramp pulse
            uint32_t endvalue_id;
            uint32_t cond_id;
            // double (*)(double in) for scalar function
            // void (*)(double *out, const double *in) for vector function
            //void (*ramp_func)(void);
            std::string ramp_func_name;
            //int64_t time = 0; // in sequence time unit
            //int64_t len = 0; // in sequence time unit
            //int64_t start_step = 0;
            //int64_t end_step = 0;
            friend class BasicSeq; // might not be needed
        };
        std::vector<std::vector<Pulse>> m_pulses;
        // maps from value IDs from global sequence to the one specific to this basic sequence
        std::vector<uint32_t> val_map;
        std::vector<HostSeq::Type> types;
        bool obj_file_generated = false;
        llvm::SmallVector<char,0> obj_file;
        std::string val_array_name;
        std::vector<HostSeq::Value> vals;
        // interpolated datas
        std::map<uint64_t, std::string> iData_name_map;
        // All stuff for network
        std::vector<uint8_t> interp_bc;
        std::vector<uint8_t> obj_file_bc;
        std::vector<uint8_t> types_bc;
        std::vector<uint8_t> pulses_bc;
        bool first_only = false;
        bool obj_file_sent = false;
        bool iData_sent = false;
        bool all_time_known;
        bool all_val_known;
        uint32_t n_const_map = 0;
        uint32_t n_nonconst_map = 0;
    };

    DummyClient();
    //{
    //    std::cout << "calling constructor" << std::endl;
    //    m_seq_cnt_offset = seqcount; //seqcount is a global
    //}
    ~DummyClient();
    const uint8_t *get_channel_info(uint32_t *sz);

public:
    void add_channel(uint32_t chn_id, const std::string &chn_name);
    void prepare();
    //void generate (Manager::ExpSeq &expseq, Compiler &compiler);

    void init_run();
    //void prepare_run(HostSeq &host_seq) ;
    void pre_run() ;
    void start() ;
    //void cancel(HostSeq &host_seq);
    void wait() ;
    //void finish_run(HostSeq &host_seq);
    // void config(const YAML::Node&);

    void sort_channels();
    // Return whether pregeneration competed
    //bool pregenerate(Manager::ExpSeq &expseq, Compiler &compiler, uint32_t seq_idx);
    //void populate_values(HostSeq &host_seq, BasicSeq &bseq);


    void addConstVal(uint32_t idx, HostSeq::Type type, HostSeq::Value val);
    void addNConstVal(uint32_t idx, HostSeq::Type type, HostSeq::Value first_val);
    void setVals(uint32_t idx, uint32_t val_idx, HostSeq::Value val);
    void addPulse(uint32_t idx, uint32_t chn_idx, uint32_t id,
                  uint32_t time_id,uint32_t len_id, uint32_t endvalue_id,
                  uint32_t cond_id, std::string func_name);
    void set_cur_seq_id(uint32_t idx);
    uint8_t get_pulse_type(DummyClient::ChnType type, bool is_fn, bool is_vector);
    void get_chn_from_lin(uint32_t lin_idx, DummyClient::ChnType &type, uint8_t &phys_chn_id, uint32_t &chn_id);
    //void config (const YAML::Node&);
    void loadConfig (const char *fname);
    std::vector<BasicSeq> m_seqs;
    // Map from sequence channel ID to physical device and channel name/ids etc.
    std::map<uint32_t,ChannelInfo> m_chn_map;
    // From linear channel ID within the backend to more detailed channel info
    std::vector<std::map<uint32_t, ChannelInfo>::iterator> m_linear_chns;
    //uint32_t m_step_Size;
    // uint32_t m_nconst_vals;
    bool m_first_bseq = false;
    uint32_t m_cur_seq_id = 0;
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
