/*************************************************************************
 *   Copyright (c) 2017 - 2019 Yichao Yu <yyc1992@gmail.com>             *
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

namespace NaCs {

class UnknownCPUInfo : public CPUInfo {
public:
    UnknownCPUInfo(std::string arch, std::string name, std::string ext_features)
        : CPUInfo(std::move(name), std::move(ext_features)),
          m_arch(arch)
    {}

private:
    const char *get_arch() const override
    {
        return m_arch.c_str();
    }

    std::string m_arch;
};

NACS_EXPORT() const CPUInfo &CPUInfo::get_host()
{
    static const UnknownCPUInfo host_info(LLVM::get_cpu_arch(),
                                          LLVM::get_cpu_name(),
                                          LLVM::get_cpu_features());
    return host_info;
}

} // Nacs
