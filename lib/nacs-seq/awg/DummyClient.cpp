//

#include "DummyClient.h"

//#include "../zynq/backend.h"

#include <memory>

namespace NaCs::Seq::AWG::Dummy {

//static Device::Register<Backend> register_backend("awg");

enum class DummyClient::ChnType : uint8_t
{
    Freq,
    Amp,
    Phase
};

struct DummyClient::ChannelInfo {
    uint8_t m_phys_chn;
    uint32_t m_chn_num;
    uint32_t linear_idx = 0; // Linear channel number for BasicSeq's list of pulses
    DummyClient::ChnType m_chn_type;
    ChannelInfo(uint8_t phys_chn, uint32_t chn_num, DummyClient::ChnType chn_type)
        : m_phys_chn(phys_chn),
          m_chn_num(chn_num),
          m_chn_type(chn_type)
    {
    }
};

struct DummyClient::BasicSeq{
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
    //std::vector<HostSeq::Type> types;
    bool obj_file_generated = false;
    llvm::SmallVector<char,0> obj_file;
    std::string val_array_name;
    //std::vector<HostSeq::Value> vals;
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
    uint32_t n_const_map;
    uint32_t n_nonconst_map;
};

template<typename T>
void write(std::vector<uint8_t> &vec, uint32_t loc, T obj)
{
    memcpy(&vec[loc], &obj, sizeof(obj));
}

template<typename T>
uint32_t write(std::vector<uint8_t> &vec, T obj)
{
    static_assert(std::is_trivial_v<T>);
    auto len = (uint32_t)vec.size();
    vec.resize(len + sizeof(obj));
    memcpy(&vec[len], &obj, sizeof(obj));
    return len;
}

template<typename T>
size_t vector_sizeof(const typename std::vector<T>& vec)
{
    return sizeof(T) * vec.size();
}

uint32_t write(std::vector<uint8_t> &vec, const void *src, size_t sz)
{
    auto len = (uint32_t)vec.size();
    vec.resize(len + sz);
    memcpy(&vec[len], src, sz);
    return len;
}

NACS_EXPORT() DummyClient::DummyClient()
{
    std::cout << "calling constructor" << std::endl;
    m_seq_cnt_offset = seqcount; //seqcount is a global
}

NACS_EXPORT() DummyClient::~DummyClient()
{
    // object files already freed...
}

NACS_EXPORT() void DummyClient::add_channel(uint32_t chn_id, const std::string &_chn_name)
{
    // OUT1/CHN1/AMP
    llvm::StringRef chn_name(_chn_name);
    llvm::SmallVector<llvm::StringRef, 3> chn_strs;
    chn_name.split(chn_strs, "/");
    if (chn_strs.size() > 3) {
        throw std::runtime_error("Only 2 / allowed in AWG channel name");
    }
    auto phys_chn_str = chn_strs[0];
    auto chn_num_str = chn_strs[1];
    auto type_str = chn_strs[2];
    if (chn_num_str.empty())
        throw std::runtime_error("No channel number in AWG channel name.");
    if (phys_chn_str.empty())
        throw std::runtime_error("No physical output number in AWG channel name.");
    if (type_str.empty())
        throw std::runtime_error("No type for channel specified in AWG channel name");
    uint8_t phys_chn_num;
    uint32_t chn_num;
    DummyClient::ChnType chn_type;
    if (!phys_chn_str.startswith("OUT"))
        throw std::runtime_error("Physical channel should be specified with OUT");
    if (!chn_num_str.startswith("CHN"))
        throw std::runtime_error("Virtual channel should be specified with CHN");
    phys_chn_str = phys_chn_str.substr(3);
    chn_num_str = chn_num_str.substr(3);
    //std::cout << phys_chn_str.str() << std::endl;
    if (phys_chn_str.getAsInteger(10, phys_chn_num))
        throw std::runtime_error("Physical output for AWG must be a number.");
    if (chn_num_str.getAsInteger(10, chn_num))
        throw std::runtime_error("Channel for AWG must be a number.");
    if (type_str == "FREQ")
    {
        //std::cout << "This is a freq" << std::endl;
        chn_type = DummyClient::ChnType::Freq;
    }
    else if (type_str == "AMP")
        chn_type = DummyClient::ChnType::Amp;
    else if (type_str == "PHASE")
        chn_type = DummyClient::ChnType::Phase;
    else {
        throw std::runtime_error("Unknown name for channel. Use FREQ, AMP, PHASE");
    }
    std::cout << "phys_chn num:" << phys_chn_num << std::endl;
    m_chn_map.try_emplace(chn_id, phys_chn_num, chn_num, chn_type);
}

static unsigned get_vector_size()
{
    return 8;
}


NACS_EXPORT() void DummyClient::sort_channels()
{
    // assignment of linear channel id
    auto nchn = (uint32_t)m_chn_map.size();
    m_linear_chns.resize(nchn);
    auto it = m_chn_map.begin();
    for (uint32_t i = 0; i < nchn; i++)
        m_linear_chns[i] = it++;
    std::sort(m_linear_chns.begin(), m_linear_chns.end(), [&] (auto it1, auto it2) {
        auto &info1 = it1->second;
        auto &info2 = it2->second;
        if (info1.m_phys_chn == info2.m_phys_chn) {
            if (info1.m_chn_num == info2.m_chn_num) {
                return info1.m_chn_type < info2.m_chn_type;
            }
            return info1.m_chn_num < info2.m_chn_num;
        }
        return info1.m_phys_chn < info2.m_phys_chn;
    });
    for (uint32_t i = 0; i < nchn; i++) {
        m_linear_chns[i]->second.linear_idx = i;
    }
}

const uint8_t *DummyClient::get_channel_info(uint32_t *sz)
{
    // might not be needed...
    return nullptr;
}

NACS_EXPORT() void DummyClient::init_run()
{
    // TO DO:
    // Ask server for client_id, if don't have one
    // Ask server for triple if not stored
    if (!m_sock)
        throw std::runtime_error("DummyClient socket not configured");
    auto it = info_map.find(m_url);
    if (it != info_map.end()) {
        auto global_info = it->second;
        m_client_id = global_info.client_id;
        m_server_id = global_info.server_id;
        m_triple_str = global_info.triple_str;
        m_cpu_str = global_info.cpu_str;
        m_feature_str = global_info.feature_str;
    }
    else {
        auto reply = m_sock->send_msg([&] (auto &sock) {
            ZMQ::send(sock, ZMQ::str_msg("req_client_id"));
        }).get();
        if (reply.empty())
            throw std::runtime_error("client id not obtained");
        auto rep_data = (const uint8_t*)reply[0].data();
        memcpy(&m_client_id, rep_data, sizeof(m_client_id));
        reply = m_sock->send_msg([&] (auto &sock) {
            ZMQ::send(sock, ZMQ::str_msg("req_server_id"));
        }).get();
        if (reply.empty())
            throw std::runtime_error("server id not obtained");
        rep_data = (const uint8_t*)reply[0].data();
        memcpy(&m_server_id, rep_data, sizeof(m_server_id));
        reply = m_sock->send_msg([&] (auto &sock) {
            ZMQ::send(sock, ZMQ::str_msg("req_triple"));
        }).get();
        if (reply.empty())
            throw std::runtime_error("triple not obtained");
        std::string triple_str((char*) reply[0].data(), reply[0].size());
        std::string cpu_str((char*) reply[1].data(), reply[1].size());
        std::string feature_str((char*) reply[2].data(), reply[2].size());
        m_triple_str = triple_str;
        m_cpu_str = cpu_str;
        m_feature_str = feature_str;
        info_map.try_emplace(m_url, m_client_id, m_server_id, m_triple_str, m_cpu_str, m_feature_str);
        // fill in id_bc
        id_bc.resize(25); // 8 bytes server_id, 8 bytes_client_id, 8 bytes, 8 bytes seq_id, 1 byte is_first_seq, Fill in first 16 bytes
        write(id_bc, (uint32_t) 0, m_server_id);
        write(id_bc, (uint32_t) 8, m_client_id);
    }
}
void DummyClient::start()
{

}
void DummyClient::wait()
{
    m_sock->send_msg([&] (auto &sock) {
            ZMQ::send_more(sock, ZMQ::str_msg("wait_seq"));
            ZMQ::send_more(sock, ZMQ::bits_msg(m_cur_wait_id));
    }).wait();
}
uint8_t DummyClient::get_pulse_type(DummyClient::ChnType type, bool is_fn, bool is_vector){
    // map from chn_info and whether it's a vector to a uint8_t describing the pulse type
    // mapping from Stream.h for AWG server
    // enum class CmdType : uint8_t
//  {
    // CmdType is a enumerated class that holds all possible commands.
    // They are represented by a uint8_t
    // Meta, // Meta command types are in CmdMeta
    // AmpSet,
    // AmpFn,
    // AmpVecFn,
    // FreqSet,
    // FreqFn,
    // FreqVecFn,
    // ModChn, // add or delete channels
    // Phase,
    // _MAX = Phase // keeps track of how many CmdType options there are
// };
    //auto type = chn_info.m_chn_type;
    if (type == DummyClient::ChnType::Phase) {
        if (is_fn) {
            throw std::runtime_error("Cannot program a pulse on a phase");
        }
        return 6;
    }
    uint8_t base = 0;
    if (type == DummyClient::ChnType::Freq) {
        base = 4;
    }
    else if (type == DummyClient::ChnType::Amp) {
        base = 1;
    }
    if (is_fn) {
        if (is_vector) {
            return base + 2;
        }
        return base + 1;
    }
    else {
        return base;
    }
}
NACS_EXPORT() void DummyClient::get_chn_from_lin(uint32_t lin_idx, DummyClient::ChnType &type, uint8_t &phys_chn_id, uint32_t &chn_id)
{
    auto chn_info = m_linear_chns[lin_idx]->second;
    type = chn_info.m_chn_type;
    //std::cout << "chn type: " << (uint8_t) type << std::endl;
    phys_chn_id = chn_info.m_phys_chn;
    //std::cout << "phys chn_id " << phys_chn_id << std::endl;
    chn_id = chn_info.m_chn_num;
    //std::cout << "chn num: " << chn_id << std::endl;
}

NACS_EXPORT() void DummyClient::loadConfig(const char *fname)
{
    auto config = YAML::LoadFile(fname);
    if (auto url_node = config["url"]) {
        auto new_url = url_node.as<std::string>();
        if (new_url != m_url) {
            m_url = std::move(new_url);
            auto &client = ZMQ::MultiClient::global();
            m_sock.reset(new ZMQ::MultiClient::SockRef(client.get_socket(m_url)));
        }
    }
    else {
        throw std::runtime_error("Missing url in config");
    }
}

}
