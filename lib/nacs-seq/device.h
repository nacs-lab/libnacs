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

#ifndef __NACS_SEQ_DEVICE_H__
#define __NACS_SEQ_DEVICE_H__

#include "compiler.h"
#include "manager.h"

#include <functional>
#include <memory>
#include <string>

namespace NaCs::Seq {

class NACS_EXPORT() Device {
public:
    static void register_factory(
        std::string backend_name,
        std::function<std::unique_ptr<Device>(Manager&, const std::string&)> factory);
    static std::unique_ptr<Device> create(Manager &mgr, const std::string &backend_name,
                                          const std::string &dev_name);
    template<typename DevClass>
    struct Register : Call {
        Register(std::string name)
            : Call([&] {
                Device::register_factory(std::move(name), [] (Manager &mgr,
                                                              const std::string &name) {
                    return std::make_unique<DevClass>(mgr, name);
                });
            })
        {}
    };

    Device(Manager &mgr, std::string name);
    virtual ~Device();

    // Compile time APIs
    virtual void add_channel(uint32_t chn_id, const std::string &chn_name) = 0;
    virtual bool check_noramp(uint32_t chn_id, const std::string &chn_name);
    virtual void prepare(Manager::ExpSeq &expseq, Compiler &compiler);
    virtual void generate(Manager::ExpSeq &expseq, Compiler &compiler) = 0;

    // Runtime APIs
    virtual void init_run(HostSeq &host_seq);
    // Not allowed to throw
    virtual void prepare_run(HostSeq &host_seq);
    virtual void pre_run(HostSeq &host_seq) = 0;
    virtual void start(HostSeq &host_seq) = 0;
    virtual void cancel(HostSeq &host_seq) = 0;
    // Not allowed to throw
    virtual void wait(HostSeq &host_seq) = 0;
    virtual void finish_run(HostSeq &host_seq);

    virtual void config(const YAML::Node&) = 0;
    virtual void parse_data(const uint8_t *data, size_t len);

    const std::string &name() const
    {
        return m_name;
    }
    Manager &mgr() const
    {
        return m_mgr;
    }
    virtual uint32_t refresh_restart()
    {
        return 0;
    }
private:

    Manager &m_mgr;
    const std::string m_name;
};

}

#endif
