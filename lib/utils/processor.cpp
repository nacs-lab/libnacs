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
#include <cstring>
#include <iostream>

namespace NaCs {

namespace {

template<typename T1, typename T2, typename T3>
static inline bool test_bits(T1 v, T2 mask, T3 test)
{
    return T3(v & mask) == test;
}

template<typename T1, typename T2>
static inline bool test_all_bits(T1 v, T2 mask)
{
    return test_bits(v, mask, mask);
}

template<typename T1, typename T2>
static inline bool test_nbit(const T1 &bits, T2 _bitidx)
{
    auto bitidx = static_cast<uint32_t>(_bitidx);
    auto u32idx = bitidx / 32;
    auto bit = bitidx % 32;
    return (bits[u32idx] & (1 << bit)) != 0;
}

template<typename T>
static inline void unset_bits(T &bits)
{
    (void)bits;
}

template<typename T, typename T1, typename... Rest>
static inline void unset_bits(T &bits, T1 _bitidx, Rest... rest)
{
    auto bitidx = static_cast<uint32_t>(_bitidx);
    auto u32idx = bitidx / 32;
    auto bit = bitidx % 32;
    bits[u32idx] = bits[u32idx] & ~uint32_t(1 << bit);
    unset_bits(bits, rest...);
}

template<typename T, typename T1>
static inline void set_bit(T &bits, T1 _bitidx, bool val)
{
    auto bitidx = static_cast<uint32_t>(_bitidx);
    auto u32idx = bitidx / 32;
    auto bit = bitidx % 32;
    if (val) {
        bits[u32idx] = bits[u32idx] | uint32_t(1 << bit);
    }
    else {
        bits[u32idx] = bits[u32idx] & ~uint32_t(1 << bit);
    }
}

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

template<size_t n>
struct FeatureList : std::array<uint32_t,n> {
private:
    static inline constexpr uint32_t _add_u32(uint32_t mask, uint32_t)
    {
        return mask;
    }

    template<typename T, typename... Rest>
    static inline constexpr uint32_t _add_u32(uint32_t mask, uint32_t u32idx,
                                              T bit, Rest... args)
    {
        return _add_u32(mask | ((int(bit) >= 0 && int(bit) / 32 == (int)u32idx) ?
                                (1 << (int(bit) % 32)) : 0),
                        u32idx, args...);
    }

    template<typename... Args>
    static inline constexpr uint32_t _get_u32(uint32_t u32idx, Args... args)
    {
        return _add_u32(uint32_t(0), u32idx, args...);
    }

    template<size_t... I, typename... Args>
    static inline constexpr FeatureList<n>
    _get(std::index_sequence<I...>, Args... args)
    {
        return FeatureList<n>(_get_u32(I, args...)...);
    }

    template<size_t... I>
    inline constexpr FeatureList<n>
    _or(std::index_sequence<I...>, const FeatureList<n> &other) const
    {
        return FeatureList<n>(((*this)[I] | other[I])...);
    }

    template<size_t... I>
    inline constexpr FeatureList<n>
    _and(std::index_sequence<I...>, const FeatureList<n> &other) const
    {
        return FeatureList<n>(((*this)[I] & other[I])...);
    }

    template<size_t... I>
    inline constexpr FeatureList<n>
    _not(std::index_sequence<I...>) const
    {
        return FeatureList<n>((~(*this)[I])...);
    }

public:
    using std::array<uint32_t,n>::array;
    template<typename... Args>
    constexpr FeatureList(Args... args)
        : std::array<uint32_t,n>{args...}
    {}

    template<typename... Args>
    static inline constexpr FeatureList<n> get(Args... args)
    {
        return _get(std::make_index_sequence<n>(), args...);
    }

    inline constexpr FeatureList<n> operator|(const FeatureList<n> &other) const
    {
        return _or(std::make_index_sequence<n>(), other);
    }

    inline constexpr FeatureList<n> operator&(const FeatureList<n> &other) const
    {
        return _and(std::make_index_sequence<n>(), other);
    }

    inline constexpr FeatureList<n> operator~() const
    {
        return _not(std::make_index_sequence<n>());
    }
};

struct FeatureName {
    const char *name;
    uint32_t bit; // bit index into a `uint32_t` array;
    uint32_t llvmver; // 0 if it is available on the oldest LLVM version we support
};

template<typename CPU, size_t n>
struct CPUSpec {
    const char *name;
    CPU cpu;
    CPU fallback;
    uint32_t llvmver;
    FeatureList<n> features;
};

struct FeatureDep {
    uint32_t feature;
    uint32_t dep;
};

// Recursively enable all features that the current feature set depends on.
template<typename FeatureList, typename Deps>
static inline void enable_depends(FeatureList &features, const Deps &deps)
{
    bool changed = true;
    while (changed) {
        changed = false;
        for (auto &dep: deps) {
            if (!test_nbit(features, dep.feature) || test_nbit(features, dep.dep))
                continue;
            set_bit(features, dep.dep, true);
            changed = true;
        }
    }
}

// Recursively disable all features that the current feature set does not provide.
template<typename FeatureList, typename Deps>
static inline void disable_depends(FeatureList &features, const Deps &deps)
{
    bool changed = true;
    while (changed) {
        changed = false;
        for (auto &dep: deps) {
            if (!test_nbit(features, dep.feature) || test_nbit(features, dep.dep))
                continue;
            unset_bits(features, dep.feature);
            changed = true;
        }
    }
}

template<typename Specs>
static auto find_cpu(uint32_t cpu, const Specs &specs) -> decltype(&*std::begin(specs))
{
    for (auto &spec: specs) {
        if (cpu == uint32_t(spec.cpu)) {
            return &spec;
        }
    }
    return nullptr;
}

template<typename Specs>
static auto find_cpu(const char *name, const Specs &specs) -> decltype(&*std::begin(specs))
{
    for (auto &spec: specs) {
        if (strcmp(name, spec.name) == 0) {
            return &spec;
        }
    }
    return nullptr;
}

template<typename Specs>
static const char *find_cpu_name(uint32_t cpu, const Specs &specs)
{
    if (auto *spec = find_cpu(cpu, specs))
        return spec->name;
    return "generic";
}

template<typename Features>
static uint32_t find_feature_bit(const Features &features, const char *str, size_t len)
{
    for (auto &feature: features) {
        if (strncmp(feature.name, str, len) == 0 && feature.name[len] == 0) {
            return feature.bit;
        }
    }
    return (uint32_t)-1;
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

#include "processor_fallback.cpp"

#include "processor_x86.cpp"

namespace NaCs {

namespace {

struct InfoBuilder {
    std::unique_ptr<CPUInfo> parse(const char *str);

protected:
    virtual void add_feature(bool enable, const char *feature, size_t len)
    {
        if (!features.empty())
            features.push_back(',');
        features.push_back(enable ? '+' : '-');
        features.append(feature, len);
    }
    virtual CPUInfo *create() = 0;

    std::string name;
    std::string features;
};

struct UnknownInfoBuilder : InfoBuilder {
    UnknownInfoBuilder(std::string arch)
        : m_arch(std::move(arch))
    {
    }

private:
    CPUInfo *create() override
    {
        return new UnknownCPUInfo(std::move(m_arch), std::move(name), std::move(features));
    }

    std::string m_arch;
};

template<typename FeatureList>
struct KnownInfoBuilder : InfoBuilder {
protected:
    template<typename Features>
    void _add_feature(const Features &feature_names, bool enable,
                      const char *feature, size_t len)
    {
        auto fbit = find_feature_bit(feature_names, feature, len);
        if (fbit != (uint32_t)-1) {
            set_bit(enable ? en : dis, fbit, true);
            return;
        }
        InfoBuilder::add_feature(enable, feature, len);
    }

    FeatureList en{};
    FeatureList dis{};
};

struct X86InfoBuilder : KnownInfoBuilder<X86::FeatureList> {
    X86InfoBuilder(bool is_x64)
        : m_is_x64(is_x64)
    {
    }

private:
    void add_feature(bool enable, const char *feature, size_t len) override
    {
        _add_feature(X86::Feature::names, enable, feature, len);
    }
    CPUInfo *create() override
    {
        return new X86::Info(m_is_x64, std::move(name), std::move(features), en, dis);
    }

    bool m_is_x64;
};

std::unique_ptr<CPUInfo> InfoBuilder::parse(const char *str)
{
    bool name_set = false;
    bool done = false;
    auto start = str;
    for (auto p = str; !done; p++) {
        auto c = *p;
        done = !c;
        if (!done && c != ',')
            continue;
        if (!name_set) {
            name.append(start, p - start);
            name_set = true;
        }
        else {
            bool enable = true;
            if (*start == '+') {
                start++;
            }
            else if (*start == '-') {
                enable = false;
                start++;
            }
            add_feature(enable, start, p - start);
        }
        start = p + 1;
    }
    return std::unique_ptr<CPUInfo>(create());
}

} // (anonymous)

NACS_EXPORT() std::unique_ptr<CPUInfo> CPUInfo::create(const char *arch, const char *str)
{
    if (strcmp(arch, "x86_64") == 0 || strcmp(arch, "x86-64") == 0)
        return X86InfoBuilder(true).parse(str);
    if (strcmp(arch, "i386") == 0 || strcmp(arch, "i486") == 0 ||
        strcmp(arch, "i586") == 0 || strcmp(arch, "i686") == 0)
        return X86InfoBuilder(false).parse(str);
    return UnknownInfoBuilder(arch).parse(str);
}

} // NaCs
