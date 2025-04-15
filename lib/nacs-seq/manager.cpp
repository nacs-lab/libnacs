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
    for (auto &[name, dev]: m_devices) {
        mgr().call_guarded([&] {
	    dev->post_run(host_seq);
        });
    }
    if (id == 0) {
        mgr().add_debug_printf("Sequence ended\n");
    }
    return id;
}

NACS_EXPORT() uint32_t Manager::ExpSeq::refresh_device_restart(const char *dname)
{
    auto it = m_devices.find(dname);
    if (it != m_devices.end()) {
        Device* dev = it->second.get();
        return dev->refresh_restart();
    }
    return uint32_t(-1);
}

NACS_EXPORT() Manager::ExpSeq::ChnOutput* Manager::ExpSeq::get_nominal_output(uint64_t pts_per_ramp, size_t *out_sz)
{
    // pts_per_ramp includes the first and last point, so it must be at least 2. It will be set to 2 if not the case.
    // Temporary insofar as keeping the start value, and pulses for this channel. 
    struct FilledPulse {
        FilledPulse(uint32_t id, int64_t time, int64_t len, double (*ramp_func)(double, void*),
                    double endvalue, bool cond)
            : id(id), time(time), len(len), ramp_func(ramp_func),
              endvalue(endvalue), cond(cond)
        {
            // empty
        }
        uint32_t id;

        int64_t time;
        int64_t len; // 0 for a set pulse and nonzero for a ramp.
        double (*ramp_func)(double, void*);

        double endvalue;
        bool cond;
    };
    struct TempChnOutput {
        TempChnOutput(std::string name, double start_val)
            : m_name(name),
              m_start_val(start_val)
            {
                //
            }
        std::string m_name;
        double m_start_val;
        std::vector<FilledPulse> m_pulses;
        std::vector<int64_t> times;
        std::vector<double> values;
        std::vector<uint32_t> pulse_ids;
    };
    uint32_t nchannels = host_seq.nchannels;
    auto seq_idx = host_seq.cur_seq_idx();
    if (nchannels == 0 || seq_idx == uint32_t(-1)) {
        *out_sz = 0;
        return nullptr;
    }

    std::map<uint32_t,TempChnOutput> temp_chn_outputs;
    std::vector<Manager::ExpSeq::ChnOutput> outputs;

    // Initialize the temp output for each channel.
    for (auto &[name, chn_id]: m_chn_map) {
        std::string chn_name = name;
        double start_val = host_seq.start_values[chn_id - 1].f64;
        printf("Channel %s (ID: %d) start value: %f\n", chn_name.c_str(), chn_id, start_val);
        temp_chn_outputs.try_emplace(chn_id - 1, chn_name, start_val);
    }
    // Collect the pulses along with their filled values for each channel.
    auto &bseq = host_seq.seqs[seq_idx];
    // Keep track of the max time
    int64_t max_time = 0;
    for (auto &pulse: bseq.pulses) {
        if (pulse.is_measure())
            continue;
        int64_t time = host_seq.get_time(pulse.time);
        // Using this call double get_value(uint32_t i, uint32_t seq_idx) const
        int64_t len = pulse.len == uint32_t(-1) ? 0 : round<int64_t>(host_seq.get_value(pulse.len, seq_idx));
        double endvalue = host_seq.get_value(pulse.endvalue, seq_idx);
        bool cond = pulse.cond == uint32_t(-1) ? true : (bool) host_seq.get_value(pulse.cond, seq_idx);
        printf("Pulse ID: %d, Channel: %d, Time: %ld, Length: %ld, End Value: %f, Condition: %d\n",
               pulse.id, pulse.chn, time, len, endvalue, cond);
        auto it = temp_chn_outputs.find(pulse.chn - 1);
        if (it == temp_chn_outputs.end()) {
            throw std::runtime_error("Missing channel ID");
        } else {
            auto &temp_output = it->second;
            temp_output.m_pulses.emplace_back(pulse.id, time, len, pulse.ramp_func, endvalue, cond);
        }
    }
    // With all the pulses, we now build our output.
    for (auto &[chn_id, temp_output]: temp_chn_outputs) {
        // First we sort them.
        auto &pulses = temp_output.m_pulses;
        std::sort(pulses.begin(), pulses.end(), [] (auto &p1, auto &p2) {
            if (p1.time < p2.time)
                return true;
            if (p1.time > p2.time)
                return false;
            return p1.id < p2.id;
        });
        // Start building the output with the start value
        auto &times = temp_output.times;
        auto &values = temp_output.values;
        auto &pulse_ids = temp_output.pulse_ids;
        times.push_back(0);
        values.push_back(temp_output.m_start_val);
        pulse_ids.push_back(uint32_t(-1)); // the start value is not identified with any particular pulse.
        auto npulses = pulses.size();
        for (uint32_t i = 0; i < npulses;) {
            auto &pulse = pulses[i];
            // This check should rarely execute, except if the first pulse is disabled.
            if (!pulse.cond) {
                // Skip disabled pulse
                i++;
                continue;
            }
            // Find next pulse which is important for ensuring that this set pulse will take place and
            // determines the end time for a potentially interrupted ramp.
            bool found = false;
            FilledPulse *next_pulse_ptr = nullptr;
            i++; // Move to the next pulse
            while (!found) {
                if (i >= npulses) {
                    found = true;
                }
                else {
                    if (pulses[i].cond) {
                        next_pulse_ptr = &pulses[i];
                        found = true;
                    } else {
                        // Skip disabled pulse
                        i++;
                    }
                }
            }
            // Below here, i does not need to be advanced!
            if (pulse.len == 0) {
                // Set pulse
                // Only valid if next pulse is not at the same time.
                if (next_pulse_ptr && (next_pulse_ptr->time <= pulse.time)) {
                    // Skip current pulse
                    continue;
                }
                // If previous point is not at the current time, then add a point with the previous value.
                if (times.back() != pulse.time) {
                    times.push_back(pulse.time);
                    values.push_back(values.back()); // Use the last value.
                    pulse_ids.push_back(pulse_ids.back()); // No pulse ID for this point.
                }
                times.push_back(pulse.time);
                values.push_back(pulse.endvalue);
                pulse_ids.push_back(pulse.id);
            } else {
                // Ramp pulse
                // Get end time by looking for the next pulse.
                int64_t end_time;
                if (next_pulse_ptr) {
                    end_time = min(next_pulse_ptr->time, pulse.time + pulse.len);
                }
                else {
                    end_time = pulse.time + pulse.len;
                }
                int64_t ramp_time = end_time - pulse.time;
                assert(ramp_time >= 0);
                if (ramp_time == 0) {
                    // If the ramp time is 0, we let the next pulse dictate what happens.
                    continue;
                }
                // If previous point is not at the current time, then add a point with the previous value.
                if (times.back() != pulse.time) {
                    times.push_back(pulse.time);
                    values.push_back(values.back()); // Use the last value.
                    pulse_ids.push_back(pulse_ids.back()); // No pulse ID for this point.
                }
                // Calculate the ramp with at least pts_per_ramp points.
                uint64_t npts = max(2, min((uint64_t) ramp_time, pts_per_ramp));
                for (uint64_t j = 0; j < npts; j++) {
                    int64_t this_pt = round<int64_t>(j * ramp_time / (npts - 1));
                    int64_t this_time = pulse.time + this_pt;
                    times.push_back(this_time);
                    values.push_back(pulse.ramp_func(this_pt, host_seq.values.data()));
                    pulse_ids.push_back(pulse.id);
                }
            }
        }
        int64_t last_time = times.back();
        if (last_time > max_time) {
            max_time = last_time;
        }
    }

    // Now, we have built our output for each channel, we need to allocate some memory to pass it onto the next user.
    // We also append an extra point if not already at the max time
    for (auto &[chn_id, temp_output]: temp_chn_outputs) {
        auto &times = temp_output.times;
        auto &values = temp_output.values;
        auto &pulse_ids = temp_output.pulse_ids;
        // If the last time is not the max time, we add an extra point with the last value.
        if (times.back() != max_time) {
            times.push_back(max_time);
            values.push_back(values.back());
            pulse_ids.push_back(pulse_ids.back());
        }
        auto names_ptr = (char*) malloc(temp_output.m_name.size() + 1);
        auto times_ptr = (int64_t*) malloc(times.size() * sizeof(int64_t));
        auto values_ptr = (double*) malloc(values.size() * sizeof(double));
        auto pulse_ids_ptr = (uint32_t*) malloc(pulse_ids.size() * sizeof(uint32_t));
        size_t out_pts = times.size(); // All sizes should be the same.
        memcpy(names_ptr, temp_output.m_name.data(), temp_output.m_name.size() + 1);
        memcpy(times_ptr, times.data(), out_pts * sizeof(int64_t));
        memcpy(values_ptr, values.data(), out_pts * sizeof(double));
        memcpy(pulse_ids_ptr, pulse_ids.data(), out_pts * sizeof(uint32_t));
        // Create the output object
        outputs.emplace_back(names_ptr, temp_output.m_name.size() + 1, times_ptr, values_ptr, pulse_ids_ptr, out_pts);
    }

    // Last allocation of memory for the outputs array, which carries a bunch of addresses.
    auto outputs_ptr = (Manager::ExpSeq::ChnOutput*) malloc(outputs.size() * sizeof(Manager::ExpSeq::ChnOutput));
    memcpy(outputs_ptr, outputs.data(), outputs.size() * sizeof(Manager::ExpSeq::ChnOutput));
    *out_sz = outputs.size();

    return outputs_ptr;
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
    auto dump_sequences = m_dump_sequences;
    {
        Builder builder(seq);
        add_debug("Deserializig sequence.\n");
        builder.deserialize(data, size);
        if (unlikely(dump_sequences)) {
            expseq->dump.reset(new ExpSeq::Dump);
            basic_vector_ostream<std::vector<uint8_t>> stm;
            builder.print(stm);
            expseq->dump->builder = stm.get_buf();
        }
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
            if (!dev)
                throw std::runtime_error("Invalid device name.");
            auto chn_name = name.substr(sep + 1);
            assert(chn_name.size() == name_len - sep - 1);
            dev->add_channel(chn_id, chn_name);
            if (dev->check_noramp(chn_id, chn_name)) {
                builder.noramp_chns.push_back(chn_id);
            }
            expseq->m_chn_map.emplace(name, chn_id);
        }
        add_debug("Building sequence\n");
        builder.buildseq();
        if (unlikely(dump_sequences)) {
            basic_vector_ostream<std::vector<uint8_t>> stm;
            seq.print(stm);
            assert(expseq->dump);
            expseq->dump->seq = stm.get_buf();
        }
    }
    add_debug("Preparing sequence\n");
    seq.prepare();
    add_debug("Optimizing sequence\n");
    seq.optimize();
    cgctx->remove_unused_data(seq.env().llvm_module());
    if (unlikely(dump_sequences)) {
        basic_vector_ostream<std::vector<uint8_t>> stm;
        seq.print(stm);
        assert(expseq->dump);
        expseq->dump->seq_opt = stm.get_buf();
    }
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

NACS_EXPORT() uint32_t Manager::get_device_restart(const char *dname)
{
    // This function should not throw. Returns uint32_t(-1) if device not found.
    auto it = m_device_info.find(dname);
    if (it != m_device_info.end())
        return it->second.n_restarts;
    return uint32_t(-1);
}

NACS_EXPORT() bool Manager::add_device_restart(const std::string &name, uint32_t n)
{
    auto it = m_device_info.find(name);
    if (it != m_device_info.end()) {
        it->second.n_restarts = it->second.n_restarts + n;
        return true;
    }
    return false;
}

NACS_EXPORT() bool Manager::set_device_restart(const std::string &name, uint32_t n)
{
    auto it = m_device_info.find(name);
    if (it != m_device_info.end()) {
        it->second.n_restarts = n;
        return true;
    }
    return false;
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
            auto &info = m_device_info[dev_name];
            info.backend = backend_node.as<std::string>();
            info.config = dev_info["config"];
            info.marked = true;
            m_device_order.push_back(std::move(dev_name));
        }
    }
    for (auto it = m_device_info.begin(), end = m_device_info.end(); it != end;) {
        auto &this_info = it->second;
        if (!this_info.marked) {
            // not marked, so delete
            it = m_device_info.erase(it);
            continue;
        }
        else {
            this_info.marked = false; // remove the mark
        }
        ++it;
    }
}

NACS_EXPORT() char *Manager::get_config_str(const std::string &name, size_t *sz) const
{
    auto it = m_device_info.find(name);
    if (it == m_device_info.end())
        throw std::runtime_error("Invalid device name");
    YAML::Emitter out;
    out << it->second.config;
    *sz = out.size();
    return strdup(out.c_str());
}

NACS_EXPORT() int64_t Manager::tick_per_sec() const
{
    if (m_tick_per_sec == 0)
        throw std::runtime_error("Sequence time unit not initialized.");
    return m_tick_per_sec;
}

NACS_EXPORT() uint64_t Manager::max_seq_len() const
{
    if (m_max_seq_length == 0)
        throw std::runtime_error("Max sequence length not initialized.");
    return m_max_seq_length;
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

NACS_EXPORT() Signal<Manager*> &Manager::new_run()
{
    // Use local global to avoid global initialization order issues.
    static Signal<Manager*> sig;
    return sig;
}

}

extern "C" {

using namespace NaCs::Seq;

NACS_EXPORT() Manager *nacs_seq_new_manager()
{
    return new Manager;
}

NACS_EXPORT() void nacs_seq_free_manager(Manager *mgr)
{
    delete mgr;
}

NACS_EXPORT() Manager::ExpSeq *nacs_seq_manager_create_sequence(
    Manager *mgr, const uint8_t *data, size_t size)
{
    return mgr->call_guarded([&] {
        return mgr->create_sequence(data, size);
    }, nullptr);
}

NACS_EXPORT() void nacs_seq_manager_free_sequence(Manager *mgr, Manager::ExpSeq *seq)
{
    mgr->call_guarded([&] {
        mgr->free_sequence(seq);
    });
}

NACS_EXPORT() uint32_t nacs_seq_manager_get_device_restart(Manager *mgr, const char *dname)
{
    return mgr->get_device_restart(dname);
}

NACS_EXPORT() void nacs_seq_manager_load_config_file(Manager *mgr, const char *fname)
{
    mgr->call_guarded([&] {
        mgr->load_config_file(fname);
    });
}

NACS_EXPORT() void nacs_seq_manager_load_config_string(Manager *mgr, const char *str)
{
    mgr->call_guarded([&] {
        mgr->load_config_string(str);
    });
}

NACS_EXPORT() char *nacs_seq_manager_get_config_string(Manager *mgr, const char *name, size_t *sz)
{
    return mgr->call_guarded([&] {
        return mgr->get_config_str(name, sz);
    }, nullptr);
}

NACS_EXPORT() int64_t nacs_seq_manager_tick_per_sec(Manager *mgr)
{
    return mgr->call_guarded([&] {
        return mgr->tick_per_sec();
    }, 0);
}

NACS_EXPORT() uint64_t nacs_seq_manager_max_seq_len(Manager *mgr)
{
    return mgr->call_guarded([&] {
        return mgr->max_seq_len();
    }, 0);
}

NACS_EXPORT() uint8_t *nacs_seq_manager_take_messages(Manager *mgr, size_t *sz)
{
    return mgr->take_messages(sz);
}

NACS_EXPORT() void nacs_seq_manager_enable_debug(Manager *mgr, bool enable)
{
    mgr->enable_debug(enable);
}

NACS_EXPORT() bool nacs_seq_manager_debug_enabled(Manager *mgr)
{
    return mgr->debug_enabled();
}

NACS_EXPORT() void nacs_seq_manager_enable_dump(Manager *mgr, bool enable)
{
    mgr->enable_dump(enable);
}

NACS_EXPORT() bool nacs_seq_manager_dump_enabled(Manager *mgr)
{
    return mgr->dump_enabled();
}

NACS_EXPORT() void nacs_seq_manager_new_run(Manager *mgr)
{
    Manager::new_run().emit(mgr);
}

NACS_EXPORT() void nacs_seq_manager_expseq_init_run(Manager::ExpSeq *expseq)
{
    expseq->mgr().call_guarded([&] {
        expseq->init_run();
    });
}

NACS_EXPORT() void nacs_seq_manager_expseq_pre_run(Manager::ExpSeq *expseq)
{
    expseq->mgr().call_guarded([&] {
        expseq->pre_run();
    });
}

NACS_EXPORT() void nacs_seq_manager_expseq_start(Manager::ExpSeq *expseq)
{
    expseq->mgr().call_guarded([&] {
        expseq->start();
    });
}

NACS_EXPORT() void nacs_seq_manager_expseq_cancel(Manager::ExpSeq *expseq)
{
    expseq->mgr().call_guarded([&] {
        expseq->cancel();
    });
}

NACS_EXPORT() bool nacs_seq_manager_expseq_wait(Manager::ExpSeq *expseq, uint64_t timeout_ms)
{
    return expseq->mgr().call_guarded([&] {
        return expseq->wait(timeout_ms);
    }, false);
}

NACS_EXPORT() uint32_t nacs_seq_manager_expseq_post_run(Manager::ExpSeq *expseq)
{
    return expseq->mgr().call_guarded([&] {
        return expseq->post_run();
    }, 0);
}

NACS_EXPORT() double nacs_seq_manager_expseq_get_global(Manager::ExpSeq *expseq, uint32_t i)
{
    return expseq->mgr().call_guarded([&] {
        return expseq->host_seq.get_global(i);
    }, 0);
}

NACS_EXPORT() void nacs_seq_manager_expseq_set_global(Manager::ExpSeq *expseq,
                                                      uint32_t i, double val)
{
    expseq->mgr().call_guarded([&] {
        expseq->host_seq.set_global(i, val);
    });
}

NACS_EXPORT() uint64_t nacs_seq_manager_expseq_cur_bseq_length(Manager::ExpSeq *expseq)
{
    auto &host_seq = expseq->host_seq;
    auto seq_idx = host_seq.cur_seq_idx();
    if (seq_idx == uint32_t(-1)) {
        expseq->mgr().add_error("Sequence not running.\n");
        return 0;
    }
    return host_seq.seqs[seq_idx].length;
}

NACS_EXPORT() const uint8_t *nacs_seq_manager_expseq_get_builder_dump(Manager::ExpSeq *expseq,
                                                                      size_t *sz)
{
    return expseq->get_builder_dump(sz);
}

NACS_EXPORT() const uint8_t *nacs_seq_manager_expseq_get_seq_dump(Manager::ExpSeq *expseq,
                                                                  size_t *sz)
{
    return expseq->get_seq_dump(sz);
}

NACS_EXPORT() const uint8_t *nacs_seq_manager_expseq_get_seq_opt_dump(Manager::ExpSeq *expseq,
                                                                      size_t *sz)
{
    return expseq->get_seq_opt_dump(sz);
}

NACS_EXPORT() uint32_t nacs_seq_manager_expseq_refresh_device_restart(Manager::ExpSeq *expseq, const char *dname)
{
    return expseq->mgr().call_guarded([&] {
        return expseq->refresh_device_restart(dname);
    }, uint32_t(-1));
}

NACS_EXPORT() Manager::ExpSeq::ChnOutput*
nacs_seq_manager_expseq_get_nominal_output(
    Manager::ExpSeq *expseq, uint64_t pts_per_ramp, size_t *sz)
{
    return expseq->mgr().call_guarded([&] () -> Manager::ExpSeq::ChnOutput* {
            return expseq->get_nominal_output(pts_per_ramp, sz);
        }, nullptr);
}


}
