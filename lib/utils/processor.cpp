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

#include "processor.h"
#include "llvm/cpu_p.h"
#include "streams.h"

#include <array>
#include <iostream>

namespace NaCs {

namespace {

static inline std::vector<std::string>&
append_features(std::vector<std::string> &features, const std::string &ext_features)
{
    if (ext_features.empty())
        return features;
    auto add_feature =
        [&] (const char *p, size_t n) {
            if (*p == '+' || *p == '-') {
                features.emplace_back(p, n);
            }
            else {
                std::string s("+");
                s.append(p, n);
                features.push_back(std::move(s));
            }
        };
    const char *start = ext_features.c_str();
    const char *p = start;
    for (; *p; p++) {
        if (*p == ',') {
            add_feature(start, p - start);
            start = p + 1;
        }
    }
    if (p > start)
        add_feature(start, p - start);
    return features;
}

static inline std::vector<std::string>&
append_features(std::vector<std::string> &&features, const std::string &ext_features)
{
    return append_features(features, ext_features);
}

} // (anonymous)

bool CPUInfo::test_feature(int) const
{
    return false;
}

int CPUInfo::get_vector_size() const
{
    return 8;
}

const char *CPUInfo::get_name() const
{
    return name.c_str();
}

void CPUInfo::dump(std::ostream &stm) const
{
    stm << name;
    for (auto &feature: append_features({}, ext_features)) {
        stm << "," << feature;
    }
}

NACS_EXPORT() void CPUInfo::dump() const
{
    dump(std::cerr);
}

NACS_EXPORT() CPUInfo::operator std::string() const
{
    string_ostream stm;
    dump(stm);
    return stm.get_buf();
}

NACS_EXPORT() void CPUInfo::dump_llvm(std::ostream &stm) const
{
    auto target = get_llvm_target(UINT32_MAX);
    stm << "Arch: " << get_arch() << std::endl;
    stm << "CPU: " << target.first << std::endl;
    stm << "Features:";
    bool first = true;
    for (auto &feature: target.second)
        stm << (first ? " " : ", ") << feature;
    stm << std::endl;
}

NACS_EXPORT() void CPUInfo::dump_llvm() const
{
    dump_llvm(std::cerr);
}

std::pair<std::string,std::vector<std::string>> CPUInfo::get_llvm_target(uint32_t) const
{
    return {name, append_features({}, ext_features)};
}

CPUInfo::CPUInfo(std::string name, std::string ext_features)
    : name(std::move(name)),
      ext_features(std::move(ext_features))
{
}

CPUInfo::~CPUInfo()
{
}

} // NaCs
