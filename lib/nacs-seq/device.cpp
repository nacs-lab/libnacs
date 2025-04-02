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

#include "device.h"

#include <map>

namespace NaCs::Seq {

namespace {

using factory_map_t = std::map<std::string,
                               std::function<std::unique_ptr<Device>(Manager&,
                                                                     const std::string&)>>;
static factory_map_t &get_factories()
{
    // Use a local static variable instead of a global one since
    // this has a guaranteed initialization order based on first use.
    // We need to make sure this gets initialized before any backend registration happens.
    static factory_map_t factories;
    return factories;
}

}

void Device::register_factory(
    std::string backend_name,
    std::function<std::unique_ptr<Device>(Manager&,
                                          const std::string&)> factory)
{
    get_factories().emplace(std::move(backend_name), factory);
}

std::unique_ptr<Device> Device::create(Manager &mgr, const std::string &backend_name,
                                       const std::string &dev_name)
{
    auto &factories = get_factories();
    auto it = factories.find(backend_name);
    if (it == factories.end())
        return {};
    return it->second(mgr, dev_name);
}

Device::Device(Manager &mgr, std::string name)
    : m_mgr(mgr),
      m_name(name)
{
}

Device::~Device()
{
}

bool Device::check_noramp(uint32_t, const std::string&)
{
    return false;
}

void Device::prepare(Manager::ExpSeq&, Compiler&)
{
}

void Device::init_run(HostSeq&)
{
}

void Device::prepare_run(HostSeq&)
{
}

void Device::post_run(HostSeq&)
{
}

void Device::parse_data(const uint8_t*, size_t)
{
    throw std::runtime_error("Unknown backend data");
}

namespace {

// A dummy backend to help testing the frontend.
class Dummy : public Device {
public:
    using Device::Device;

private:
    void add_channel(uint32_t chn_id, const std::string &chn_name) override {}
    // bool check_noramp(uint32_t chn_id, const std::string &chn_name) override;
    // void prepare(Manager::ExpSeq &expseq, Compiler &compiler) override;
    void generate(Manager::ExpSeq &expseq, Compiler &compiler) override {}

    // void init_run(HostSeq &host_seq) override;
    // void prepare_run(HostSeq &host_seq) override;
    void pre_run(HostSeq &host_seq) override {}
    void start(HostSeq &host_seq) override {}
    void cancel(HostSeq &host_seq) override {}
    void wait(HostSeq &host_seq) override {}
    // void post_run(HostSeq &host_seq) override;

    void config(const YAML::Node&) override {}
    // void parse_data(const uint8_t *data, size_t len) override;
};

static Device::Register<Dummy> register_backend("dummy");

}

}
