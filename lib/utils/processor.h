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

// This processor detection code is a port of the version I wrote for julia.

#include "utils.h"

#ifndef __NACS_UTILS_PROCESSOR_H__
#define __NACS_UTILS_PROCESSOR_H__

#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

namespace NaCs {

#define NACS_FEATURE_DEF_NAME(name, bit, llvmver, str) NACS_FEATURE_DEF(name, bit, llvmver)
#define NACS_FEATURE_DEF(name, bit, llvmver) constexpr int name = bit;
namespace X86 {
namespace Feature {
#include "features_x86.h"
}
}
namespace AArch32 {
namespace Feature {
#include "features_aarch32.h"
}
}
namespace AArch64 {
namespace Feature {
#include "features_aarch64.h"
}
}
#undef NACS_FEATURE_DEF
#undef NACS_FEATURE_DEF_NAME

class CPUInfo {
public:
    virtual bool test_feature(int bit) const;
    // Largest (vector) register size in bytes.
    virtual int get_vector_size() const;
    virtual const std::string &get_arch() const = 0;
    virtual std::pair<std::string,std::vector<std::string>>
    get_llvm_target(uint32_t llvmver) const;
    virtual ~CPUInfo();
    virtual void dump(std::ostream&) const;

    NACS_EXPORT(utils) const std::string &get_name() const;
    NACS_EXPORT(utils) void dump() const;
    NACS_EXPORT(utils) operator std::string() const;
    NACS_EXPORT(utils) void dump_llvm(std::ostream&) const;
    NACS_EXPORT(utils) void dump_llvm() const;

    NACS_EXPORT(utils) static const CPUInfo &get_host();
    NACS_EXPORT(utils) static std::unique_ptr<CPUInfo> create(const char *arch,
                                                              const char *str);
protected:
    CPUInfo() = default;
    CPUInfo(std::string name, std::string ext_features);

    std::string name;
    std::string ext_features;
};

}

#endif
