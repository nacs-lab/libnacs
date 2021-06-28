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
namespace X86 {

enum class CPU : uint32_t {
    generic = 0,
    intel_nocona,
    intel_prescott,
    intel_atom_bonnell,
    intel_atom_silvermont,
    intel_atom_goldmont,
    intel_atom_goldmont_plus,
    intel_atom_tremont,
    intel_core2,
    intel_core2_penryn,
    intel_yonah,
    intel_corei7_nehalem,
    intel_corei7_westmere,
    intel_corei7_sandybridge,
    intel_corei7_ivybridge,
    intel_corei7_haswell,
    intel_corei7_broadwell,
    intel_corei7_skylake,
    intel_corei7_skylake_avx512,
    intel_corei7_cascadelake,
    intel_corei7_cooperlake,
    intel_corei7_cannonlake,
    intel_corei7_icelake_client,
    intel_corei7_icelake_server,
    intel_corei7_tigerlake,
    intel_knights_landing,
    intel_knights_mill,

    amd_fam10h,
    amd_athlon_fx,
    amd_athlon_64,
    amd_athlon_64_sse3,
    amd_bdver1,
    amd_bdver2,
    amd_bdver3,
    amd_bdver4,
    amd_btver1,
    amd_btver2,
    amd_k8,
    amd_k8_sse3,
    amd_opteron,
    amd_opteron_sse3,
    amd_barcelona,
    amd_znver1,
    amd_znver2,
};

static constexpr size_t feature_sz = 11;
using FeatureList = NaCs::FeatureList<feature_sz>;
using CPUSpec = NaCs::CPUSpec<CPU,feature_sz>;

namespace Feature {
static constexpr auto mask = FeatureList::get(
#define NACS_FEATURE_DEF_NAME(name, bit, llvmver, str) NACS_FEATURE_DEF(name, bit, llvmver)
#define NACS_FEATURE_DEF(name, bit, llvmver) bit,
#include "features_x86.h"
#undef NACS_FEATURE_DEF
#undef NACS_FEATURE_DEF_NAME
    -1);

static constexpr FeatureName names[] = {
#define NACS_FEATURE_DEF(name, bit, llvmver) {#name, bit, llvmver},
#define NACS_FEATURE_DEF_NAME(name, bit, llvmver, str) {str, bit, llvmver},
#include "features_x86.h"
#undef NACS_FEATURE_DEF
#undef NACS_FEATURE_DEF_NAME
};

static constexpr FeatureDep deps[] = {
    {ssse3, sse3},
    {fma, avx},
    {sse41, ssse3},
    {sse42, sse41},
    {avx, sse42},
    {f16c, avx},
    {avx2, avx},
    {vaes, avx},
    {vaes, aes},
    {vpclmulqdq, avx},
    {vpclmulqdq, pclmul},
    {avx512f, avx2},
    {avx512dq, avx512f},
    {avx512ifma, avx512f},
    {avx512pf, avx512f},
    {avx512er, avx512f},
    {avx512cd, avx512f},
    {avx512bw, avx512f},
    {avx512bf16, avx512bw},
    {avx512bitalg, avx512bw},
    {avx512vl, avx512f},
    {avx512vbmi, avx512bw},
    {avx512vbmi2, avx512bw},
    {avx512vnni, avx512f},
    {avx512vp2intersect, avx512f},
    {avx512vpopcntdq, avx512f},
    {sse4a, sse3},
    {xop, fma4},
    {fma4, avx},
    {fma4, sse4a},
};
static constexpr int x64_only[] = {cx16, sahf};

constexpr auto generic = FeatureList::get();
constexpr auto bonnell = FeatureList::get(sse3, ssse3, cx16, movbe, sahf);
constexpr auto silvermont = bonnell | FeatureList::get(sse41, sse42, popcnt,
                                                       pclmul, prfchw, rdrnd);
constexpr auto goldmont = silvermont | FeatureList::get(aes, sha, rdseed, xsave, xsaveopt,
                                                        xsavec, xsaves, clflushopt, fsgsbase);
constexpr auto goldmont_plus = goldmont | FeatureList::get(ptwrite, rdpid); // sgx
constexpr auto tremont = goldmont_plus | FeatureList::get(clwb, gfni);
constexpr auto knl = FeatureList::get(sse3, ssse3, sse41, sse42, cx16, sahf, popcnt,
                                      aes, pclmul, avx, xsave, xsaveopt, rdrnd, f16c, fsgsbase,
                                      avx2, bmi, bmi2, fma, lzcnt, movbe, adx, rdseed, prfchw,
                                      avx512f, avx512er, avx512cd, avx512pf, prefetchwt1);
constexpr auto knm = knl | FeatureList::get(avx512vpopcntdq);
constexpr auto yonah = FeatureList::get(sse3);
constexpr auto prescott = yonah;
constexpr auto core2 = FeatureList::get(sse3, ssse3, cx16, sahf);
constexpr auto nocona = FeatureList::get(sse3, cx16);
constexpr auto penryn = nocona | FeatureList::get(ssse3, sse41, sahf);
constexpr auto nehalem = penryn | FeatureList::get(sse42, popcnt);
constexpr auto westmere = nehalem | FeatureList::get(pclmul);
constexpr auto sandybridge = westmere | FeatureList::get(avx, xsave, xsaveopt);
constexpr auto ivybridge = sandybridge | FeatureList::get(rdrnd, f16c, fsgsbase);
constexpr auto haswell = ivybridge | FeatureList::get(avx2, bmi, bmi2, fma, lzcnt, movbe);
constexpr auto broadwell = haswell | FeatureList::get(adx, rdseed, prfchw);
constexpr auto skylake = broadwell | FeatureList::get(aes, xsavec, xsaves, clflushopt); // sgx
constexpr auto skx = skylake | FeatureList::get(avx512f, avx512cd, avx512dq, avx512bw, avx512vl,
                                                pku, clwb);
constexpr auto cascadelake = skx | FeatureList::get(avx512vnni);
constexpr auto cooperlake = cascadelake | FeatureList::get(avx512bf16);
constexpr auto cannonlake = skylake | FeatureList::get(avx512f, avx512cd, avx512dq, avx512bw,
                                                       avx512vl, pku, avx512vbmi, avx512ifma,
                                                       sha); // sgx
constexpr auto icelake = cannonlake | FeatureList::get(avx512bitalg, vaes, avx512vbmi2,
                                                       vpclmulqdq, avx512vpopcntdq,
                                                       gfni, clwb, rdpid);
constexpr auto icelake_server = icelake | FeatureList::get(pconfig, wbnoinvd);
constexpr auto tigerlake = icelake | FeatureList::get(avx512vp2intersect, movdiri,
                                                      movdir64b, shstk);

constexpr auto k8_sse3 = FeatureList::get(sse3, cx16);
constexpr auto amdfam10 = k8_sse3 | FeatureList::get(sse4a, lzcnt, popcnt, sahf);

constexpr auto btver1 = amdfam10 | FeatureList::get(ssse3, prfchw);
constexpr auto btver2 = btver1 | FeatureList::get(sse41, sse42, avx, aes, pclmul, bmi, f16c,
                                                  movbe, xsave, xsaveopt);

constexpr auto bdver1 = amdfam10 | FeatureList::get(xop, fma4, avx, ssse3, sse41, sse42, aes,
                                                    prfchw, pclmul, xsave, lwp);
constexpr auto bdver2 = bdver1 | FeatureList::get(f16c, bmi, tbm, fma);
constexpr auto bdver3 = bdver2 | FeatureList::get(xsaveopt, fsgsbase);
constexpr auto bdver4 = bdver3 | FeatureList::get(avx2, bmi2, mwaitx);

constexpr auto znver1 = haswell | FeatureList::get(adx, aes, clflushopt, clzero, mwaitx, prfchw,
                                                   rdseed, sha, sse4a, xsavec, xsaves);
constexpr auto znver2 = znver1 | FeatureList::get(clwb, rdpid, wbnoinvd);

} // Feature

namespace {

static constexpr CPUSpec cpus[] = {
    {"generic", CPU::generic, CPU::generic, 0, Feature::generic},
    {"bonnell", CPU::intel_atom_bonnell, CPU::generic, 0, Feature::bonnell},
    {"silvermont", CPU::intel_atom_silvermont, CPU::generic, 0, Feature::silvermont},
    {"goldmont", CPU::intel_atom_goldmont, CPU::generic, 0, Feature::goldmont},
    {"goldmont-plus", CPU::intel_atom_goldmont_plus, CPU::generic, 0, Feature::goldmont_plus},
    {"tremont", CPU::intel_atom_tremont, CPU::generic, 0, Feature::tremont},
    {"core2", CPU::intel_core2, CPU::generic, 0, Feature::core2},
    {"yonah", CPU::intel_yonah, CPU::generic, 0, Feature::yonah},
    {"prescott", CPU::intel_prescott, CPU::generic, 0, Feature::prescott},
    {"nocona", CPU::intel_nocona, CPU::generic, 0, Feature::nocona},
    {"penryn", CPU::intel_core2_penryn, CPU::generic, 0, Feature::penryn},
    {"nehalem", CPU::intel_corei7_nehalem, CPU::generic, 0, Feature::nehalem},
    {"westmere", CPU::intel_corei7_westmere, CPU::generic, 0, Feature::westmere},
    {"sandybridge", CPU::intel_corei7_sandybridge, CPU::generic, 0, Feature::sandybridge},
    {"ivybridge", CPU::intel_corei7_ivybridge, CPU::generic, 0, Feature::ivybridge},
    {"haswell", CPU::intel_corei7_haswell, CPU::generic, 0, Feature::haswell},
    {"broadwell", CPU::intel_corei7_broadwell, CPU::generic, 0, Feature::broadwell},
    {"skylake", CPU::intel_corei7_skylake, CPU::generic, 0, Feature::skylake},
    {"knl", CPU::intel_knights_landing, CPU::generic, 0, Feature::knl},
    {"knm", CPU::intel_knights_mill, CPU::generic, 0, Feature::knm},
    {"skylake-avx512", CPU::intel_corei7_skylake_avx512, CPU::generic, 0, Feature::skx},
    {"cascadelake", CPU::intel_corei7_cascadelake, CPU::generic, 0, Feature::cascadelake},
    {"cooperlake", CPU::intel_corei7_cooperlake, CPU::intel_corei7_cascadelake,
     90000, Feature::cooperlake},
    {"cannonlake", CPU::intel_corei7_cannonlake, CPU::generic, 0, Feature::cannonlake},
    {"icelake-client", CPU::intel_corei7_icelake_client, CPU::generic, 0, Feature::icelake},
    {"icelake-server", CPU::intel_corei7_icelake_server, CPU::generic, 0,
     Feature::icelake_server},
    {"tigerlake", CPU::intel_corei7_tigerlake, CPU::intel_corei7_icelake_client, 100000,
     Feature::tigerlake},

    {"athlon64", CPU::amd_athlon_64, CPU::generic, 0, Feature::generic},
    {"athlon-fx", CPU::amd_athlon_fx, CPU::generic, 0, Feature::generic},
    {"k8", CPU::amd_k8, CPU::generic, 0, Feature::generic},
    {"opteron", CPU::amd_opteron, CPU::generic, 0, Feature::generic},

    {"athlon64-sse3", CPU::amd_athlon_64_sse3, CPU::generic, 0, Feature::k8_sse3},
    {"k8-sse3", CPU::amd_k8_sse3, CPU::generic, 0, Feature::k8_sse3},
    {"opteron-sse3", CPU::amd_opteron_sse3, CPU::generic, 0, Feature::k8_sse3},

    {"amdfam10", CPU::amd_fam10h, CPU::generic, 0, Feature::amdfam10},
    {"barcelona", CPU::amd_barcelona, CPU::generic, 0, Feature::amdfam10},

    {"btver1", CPU::amd_btver1, CPU::generic, 0, Feature::btver1},
    {"btver2", CPU::amd_btver2, CPU::generic, 0, Feature::btver2},

    {"bdver1", CPU::amd_bdver1, CPU::generic, 0, Feature::bdver1},
    {"bdver2", CPU::amd_bdver2, CPU::generic, 0, Feature::bdver2},
    {"bdver3", CPU::amd_bdver3, CPU::generic, 0, Feature::bdver3},
    {"bdver4", CPU::amd_bdver4, CPU::generic, 0, Feature::bdver4},

    {"znver1", CPU::amd_znver1, CPU::generic, 0, Feature::znver1},
    {"znver2", CPU::amd_znver2, CPU::amd_znver1, 90000, Feature::znver2},
};

struct Info : CPUInfo {
    Info(bool is_x64, std::string name, std::string ext_features,
         FeatureList en, FeatureList dis)
        : CPUInfo(normalize_name(std::move(name), is_x64), std::move(ext_features)),
          m_is_x64(is_x64),
          m_en(en),
          m_dis(dis)
    {
        // Use a separate function to avoid accidentally using the moved arguments...
        init();
    }
private:
    void init()
    {
        auto spec = find_cpu(name.c_str(), cpus);
        if (spec)
            m_en = m_en | spec->features;
        if (!m_is_x64)
            for (auto bit: Feature::x64_only)
                set_bit(m_dis, bit, true);
        enable_depends(m_en, Feature::deps);
        m_en = m_en & ~m_dis;
        disable_depends(m_en, Feature::deps);
        if (spec) {
            // If the base feature is known, disable everything that we didn't
            // enable explicitly.
            m_dis = Feature::mask & ~m_en;
        }
    }
    static inline std::string normalize_name(std::string name, bool is_x64)
    {
        if (name == "atom")
            return "bonnell";
        if (name == "slm")
            return "silvermont";
        if (name == "glm")
            return "goldmont";
        if (name == "corei7")
            return "nehalem";
        if (name == "corei7-avx")
            return "sandybridge";
        if (name == "core-avx-i")
            return "ivybridge";
        if (name == "core-avx2")
            return "haswell";
        if (name == "skx")
            return "skylake-avx512";
        if (is_x64 ? (name == "x86-64" || name == "x86_64") :
            (name == "pentium4" || name == "i686"))
            return "generic";
        return name;
    }
    bool test_feature(int bit) const final override
    {
        return test_nbit(m_en, bit);
    }
    int get_vector_size() const final override
    {
        if (test_feature(Feature::avx512f))
            return 64;
        if (test_feature(Feature::avx))
            return 32;
        // Require sse2
        return 16;
    }
    const std::string &get_arch() const final override
    {
        static const std::string name_64 = "x86-64";
        static const std::string name_32 = "i386";
        return m_is_x64 ? name_64 : name_32;
    }
    void dump(std::ostream &stm) const override
    {
        CPUInfo::dump(stm);
        for (auto &fename: Feature::names) {
            if (test_nbit(m_en, fename.bit)) {
                stm << ",+" << fename.name;
            }
            else if (test_nbit(m_dis, fename.bit)) {
                stm << ",-" << fename.name;
            }
        }
    }
    std::pair<std::string,std::vector<std::string>>
    get_llvm_target(uint32_t llvmver) const final override;

    bool m_is_x64;
    FeatureList m_en;
    FeatureList m_dis;
};

std::pair<std::string,std::vector<std::string>> Info::get_llvm_target(uint32_t llvmver) const
{
    std::string llvm_name = name;
    auto *spec = find_cpu(name.c_str(), cpus);
    while (spec) {
        if (spec->llvmver <= llvmver)
            break;
        spec = find_cpu((uint32_t)spec->fallback, cpus);
        llvm_name = spec->name;
    }
    // Use translate `generic` into what we actually require
    if (llvm_name == "generic")
        llvm_name = m_is_x64 ? "x86-64" : "pentium4";
    std::vector<std::string> features;
    for (auto &fename: Feature::names) {
        if (fename.llvmver > llvmver)
            continue;
        if (test_nbit(m_en, fename.bit)) {
            features.insert(features.begin(), std::string("+") + fename.name);
        }
        else if (test_nbit(m_dis, fename.bit)) {
            features.push_back(std::string("-") + fename.name);
        }
    }
    features.push_back("+sse2");
    features.push_back("+mmx");
    features.push_back("+fxsr");
    // This is required to make LLVM happy if LLVM's feature based CPU arch guess
    // returns a value that may not have 64bit support.
    // This can happen with virtualization.
    if (m_is_x64)
        features.push_back("+64bit");
    if (llvmver >= 90000)
        features.push_back("+cx8");
    append_features(features, ext_features);
    return std::make_pair(std::move(llvm_name), std::move(features));
}

#if NACS_CPU_X86 || NACS_CPU_X86_64

static inline void cpuid(int32_t info[4], int32_t typ)
{
#if defined _MSC_VER
    __cpuid(info, typ);
#else
    asm volatile (
#if NACS_CPU_X86 && defined(__PIC__)
        "xchg %%ebx, %%esi;"
        "cpuid;"
        "xchg %%esi, %%ebx;" :
        "=S" (info[1]),
#else
        "cpuid" :
        "=b" (info[1]),
#endif
        "=a" (info[0]),
        "=c" (info[2]),
        "=d" (info[3]) :
        "a" (typ)
        );
#endif
}

static inline void cpuidex(int32_t info[4], int32_t typ, int32_t subtyp)
{
#if defined _MSC_VER
    __cpuidex(info, typ, subtyp);
#else
    asm volatile (
#if NACS_CPU_X86 && defined(__PIC__)
        "xchg %%ebx, %%esi;"
        "cpuid;"
        "xchg %%esi, %%ebx;" :
        "=S" (info[1]),
#else
        "cpuid" :
        "=b" (info[1]),
#endif
        "=a" (info[0]),
        "=c" (info[2]),
        "=d" (info[3]) :
        "a" (typ),
        "c" (subtyp)
        );
#endif
}

static inline uint64_t get_xcr0()
{
#if defined _MSC_VER
    return _xgetbv(_XCR_XFEATURE_ENABLED_MASK);
#else
    uint32_t eax, edx;
    asm volatile ("xgetbv" : "=a" (eax), "=d" (edx) : "c" (0));
    return (uint64_t(edx) << 32) | eax;
#endif
}

// For CPU model and feature detection on X86

const int SIG_INTEL = 0x756e6547; // Genu
const int SIG_AMD = 0x68747541; // Auth

static CPU get_intel_processor_name(uint32_t family, uint32_t model, uint32_t brand_id,
                                    const FeatureList &features)
{
    if (brand_id != 0)
        return CPU::generic;
    switch (family) {
    case 3:
    case 4:
    case 5:
        return CPU::generic;
    case 6:
        switch (model) {
        case 0x01: // Pentium Pro processor
        case 0x03: // Intel Pentium II OverDrive processor, Pentium II processor, model 03
        case 0x05: // Pentium II processor, model 05, Pentium II Xeon processor,
            // model 05, and Intel Celeron processor, model 05
        case 0x06: // Celeron processor, model 06
        case 0x07: // Pentium III processor, model 07, and Pentium III Xeon processor, model 07
        case 0x08: // Pentium III processor, model 08, Pentium III Xeon processor,
            // model 08, and Celeron processor, model 08
        case 0x0a: // Pentium III Xeon processor, model 0Ah
        case 0x0b: // Pentium III processor, model 0Bh
        case 0x09: // Intel Pentium M processor, Intel Celeron M processor model 09.
        case 0x0d: // Intel Pentium M processor, Intel Celeron M processor, model
            // 0Dh. All processors are manufactured using the 90 nm process.
        case 0x15: // Intel EP80579 Integrated Processor and Intel EP80579
            // Integrated Processor with Intel QuickAssist Technology
            return CPU::generic;
        case 0x0e: // Intel Core Duo processor, Intel Core Solo processor, model
            // 0Eh. All processors are manufactured using the 65 nm process.
            return CPU::intel_yonah;
        case 0x0f: // Intel Core 2 Duo processor, Intel Core 2 Duo mobile
            // processor, Intel Core 2 Quad processor, Intel Core 2 Quad
            // mobile processor, Intel Core 2 Extreme processor, Intel
            // Pentium Dual-Core processor, Intel Xeon processor, model
            // 0Fh. All processors are manufactured using the 65 nm process.
        case 0x16: // Intel Celeron processor model 16h. All processors are
            // manufactured using the 65 nm process
            return CPU::intel_core2;
        case 0x17: // Intel Core 2 Extreme processor, Intel Xeon processor, model
            // 17h. All processors are manufactured using the 45 nm process.
            //
            // 45nm: Penryn , Wolfdale, Yorkfield (XE)
        case 0x1d: // Intel Xeon processor MP. All processors are manufactured using
            // the 45 nm process.
            return CPU::intel_core2_penryn;
        case 0x1a: // Intel Core i7 processor and Intel Xeon processor. All
            // processors are manufactured using the 45 nm process.
        case 0x1e: // Intel(R) Core(TM) i7 CPU         870  @ 2.93GHz.
            // As found in a Summer 2010 model iMac.
        case 0x1f:
        case 0x2e: // Nehalem EX
            return CPU::intel_corei7_nehalem;
        case 0x25: // Intel Core i7, laptop version.
        case 0x2c: // Intel Core i7 processor and Intel Xeon processor. All
            // processors are manufactured using the 32 nm process.
        case 0x2f: // Westmere EX
            return CPU::intel_corei7_westmere;
        case 0x2a: // Intel Core i7 processor. All processors are manufactured
            // using the 32 nm process.
        case 0x2d:
            return CPU::intel_corei7_sandybridge;
        case 0x3a:
        case 0x3e: // Ivy Bridge EP
            return CPU::intel_corei7_ivybridge;

            // Haswell:
        case 0x3c:
        case 0x3f:
        case 0x45:
        case 0x46:
            return CPU::intel_corei7_haswell;

            // Broadwell:
        case 0x3d:
        case 0x47:
        case 0x4f:
        case 0x56:
            return CPU::intel_corei7_broadwell;

            // Skylake:
        case 0x4e: // Skylake mobile
        case 0x5e: // Skylake desktop
        case 0x8e: // Kaby Lake mobile
        case 0x9e: // Kaby Lake desktop
        case 0xa5: // Comet Lake-H/S
        case 0xa6: // Comet Lake-U
            return CPU::intel_corei7_skylake;

            // Skylake Xeon:
        case 0x55:
            if (test_nbit(features, Feature::avx512bf16))
                return CPU::intel_corei7_cooperlake;
            if (test_nbit(features, Feature::avx512vnni))
                return CPU::intel_corei7_cascadelake;
            return CPU::intel_corei7_skylake_avx512;

            // Cannonlake:
        case 0x66:
            return CPU::intel_corei7_cannonlake;

            // Icelake:
        case 0x7d:
        case 0x7e:
        case 0x9d:
            return CPU::intel_corei7_icelake_client;

            // Icelake Xeon:
        case 0x6a:
        case 0x6c:
            return CPU::intel_corei7_icelake_server;

            // Tiger Lake
        case 0x8c:
        case 0x8d:
            return CPU::intel_corei7_tigerlake;

        case 0x1c: // Most 45 nm Intel Atom processors
        case 0x26: // 45 nm Atom Lincroft
        case 0x27: // 32 nm Atom Medfield
        case 0x35: // 32 nm Atom Midview
        case 0x36: // 32 nm Atom Midview
            return CPU::intel_atom_bonnell;

            // Atom Silvermont codes from the Intel software optimization guide.
        case 0x37:
        case 0x4a:
        case 0x4d:
        case 0x5d:
            // Airmont
        case 0x4c:
        case 0x5a:
        case 0x75:
            return CPU::intel_atom_silvermont;

            // Goldmont:
        case 0x5c:
        case 0x5f:
            return CPU::intel_atom_goldmont;
        case 0x7a:
            return CPU::intel_atom_goldmont_plus;
        case 0x86:
        case 0x96:
        case 0x9c:
            return CPU::intel_atom_tremont;

        case 0x57:
            return CPU::intel_knights_landing;

        case 0x85:
            return CPU::intel_knights_mill;

        default:
            return CPU::generic;
        }
        break;
    case 15: {
        switch (model) {
        case 0: // Pentium 4 processor, Intel Xeon processor. All processors are
            // model 00h and manufactured using the 0.18 micron process.
        case 1: // Pentium 4 processor, Intel Xeon processor, Intel Xeon
            // processor MP, and Intel Celeron processor. All processors are
            // model 01h and manufactured using the 0.18 micron process.
        case 2: // Pentium 4 processor, Mobile Intel Pentium 4 processor - M,
            // Intel Xeon processor, Intel Xeon processor MP, Intel Celeron
            // processor, and Mobile Intel Celeron processor. All processors
            // are model 02h and manufactured using the 0.13 micron process.
        default:
            return CPU::generic;

        case 3: // Pentium 4 processor, Intel Xeon processor, Intel Celeron D
            // processor. All processors are model 03h and manufactured using
            // the 90 nm process.
        case 4: // Pentium 4 processor, Pentium 4 processor Extreme Edition,
            // Pentium D processor, Intel Xeon processor, Intel Xeon
            // processor MP, Intel Celeron D processor. All processors are
            // model 04h and manufactured using the 90 nm process.
        case 6: // Pentium 4 processor, Pentium D processor, Pentium processor
            // Extreme Edition, Intel Xeon processor, Intel Xeon processor
            // MP, Intel Celeron D processor. All processors are model 06h
            // and manufactured using the 65 nm process.
#if NACS_CPU_X86_64
            return CPU::intel_nocona;
#else
            return CPU::intel_prescott;
#endif
        }
    }
    default:
        break; /*"generic"*/
    }
    return CPU::generic;
}

static CPU get_amd_processor_name(uint32_t family, uint32_t model, const FeatureList &features)
{
    switch (family) {
    case 4:
    case 5:
    case 6:
    default:
        return CPU::generic;
    case 15:
        if (test_nbit(features, Feature::sse3))
            return CPU::amd_k8_sse3;
        switch (model) {
        case 1:
            return CPU::amd_opteron;
        case 5:
            return CPU::amd_athlon_fx;
        default:
            return CPU::amd_athlon_64;
        }
    case 16:
        switch (model) {
        case 2:
            return CPU::amd_barcelona;
        case 4:
        case 8:
        default:
            return CPU::amd_fam10h;
        }
    case 20:
        return CPU::amd_btver1;
    case 21:
        if (model >= 0x50 && model <= 0x6f)
            return CPU::amd_bdver4;
        if (model >= 0x30 && model <= 0x3f)
            return CPU::amd_bdver3;
        if (model >= 0x10 && model <= 0x1f)
            return CPU::amd_bdver2;
        if (model <= 0x0f)
            return CPU::amd_bdver1;
        return CPU::amd_btver1; // fallback
    case 22:
        return CPU::amd_btver2;
    case 23:
        if ((model >= 0x30 && model <= 0x3f) || model == 0x71)
            return CPU::amd_znver2;
        if (model <= 0x0f)
            return CPU::amd_znver1;
        return CPU::amd_btver1;
    }
}

static inline void features_disable_avx512(FeatureList &features)
{
    using namespace Feature;
    unset_bits(features, avx512f, avx512dq, avx512ifma, avx512pf, avx512er, avx512cd,
               avx512bw, avx512vl, avx512vbmi, avx512vpopcntdq, avx512vbmi2, avx512vnni,
               avx512bitalg, avx512vp2intersect, avx512bf16);
}

static inline void features_disable_avx(FeatureList &features)
{
    using namespace Feature;
    unset_bits(features, avx, Feature::fma, f16c, xsave, avx2, xop, fma4,
               xsaveopt, xsavec, xsaves, vaes, vpclmulqdq);
}

static std::pair<CPU,FeatureList> _get_host_cpu()
{
    FeatureList features = {};

    int32_t info0[4];
    cpuid(info0, 0);
    uint32_t maxleaf = info0[0];
    if (maxleaf < 1)
        return std::make_pair(CPU::generic, features);
    int32_t info1[4];
    cpuid(info1, 1);

    auto vendor = info0[1];
    auto brand_id = info1[1] & 0xff;

    auto family = (info1[0] >> 8) & 0xf; // Bits 8 - 11
    auto model = (info1[0] >> 4) & 0xf;  // Bits 4 - 7
    if (family == 6 || family == 0xf) {
        if (family == 0xf)
            // Examine extended family ID if family ID is F.
            family += (info1[0] >> 20) & 0xff; // Bits 20 - 27
        // Examine extended model ID if family ID is 6 or F.
        model += ((info1[0] >> 16) & 0xf) << 4; // Bits 16 - 19
    }

    // Fill in the features
    features[0] = info1[2];
    features[1] = info1[3];
    if (maxleaf >= 7) {
        int32_t info7[4];
        cpuidex(info7, 7, 0);
        features[2] = info7[1];
        features[3] = info7[2];
        features[4] = info7[3];
    }
    int32_t infoex0[4];
    cpuid(infoex0, 0x80000000);
    uint32_t maxexleaf = infoex0[0];
    if (maxexleaf >= 0x80000001) {
        int32_t infoex1[4];
        cpuid(infoex1, 0x80000001);
        features[5] = infoex1[2];
        features[6] = infoex1[3];
    }
    if (maxleaf >= 0xd) {
        int32_t infod[4];
        cpuidex(infod, 0xd, 0x1);
        features[7] = infod[0];
    }
    if (maxexleaf >= 0x80000008) {
        int32_t infoex8[4];
        cpuidex(infoex8, 0x80000008, 0);
        features[8] = infoex8[1];
    }
    if (maxleaf >= 7) {
        int32_t info7[4];
        cpuidex(info7, 7, 1);
        features[9] = info7[0];
    }
    if (maxleaf >= 0x14) {
        int32_t info14[4];
        cpuidex(info14, 0x14, 0);
        features[10] = info14[1];
    }

    // Fix up AVX bits to account for OS support and match LLVM model
    uint64_t xcr0 = 0;
    const uint32_t avx_mask = (1 << 27) | (1 << 28);
    bool hasavx = test_all_bits(features[0], avx_mask);
    if (hasavx) {
        xcr0 = get_xcr0();
        hasavx = test_all_bits(xcr0, 0x6);
    }
    unset_bits(features, 32 + 27);
    if (!hasavx)
        features_disable_avx(features);
    bool hasavx512save = hasavx && test_all_bits(xcr0, 0xe0);
    if (!hasavx512save)
        features_disable_avx512(features);
    // Ignore feature bits that we are not interested in.
    features = features & Feature::mask;

    CPU cpu;
    if (vendor == SIG_INTEL) {
        cpu = get_intel_processor_name(family, model, brand_id, features);
    }
    else if (vendor == SIG_AMD) {
        cpu = get_amd_processor_name(family, model, features);
    }
    else {
        cpu = CPU::generic;
    }

    return std::make_pair(cpu, features);
}

static Info create_host_info()
{
    auto host = _get_host_cpu();
    std::string name = host.first != CPU::generic ?
        std::string(find_cpu_name(uint32_t(host.first), cpus)) :
        LLVM::get_cpu_name();
    return Info(NACS_CPU_X86_64, std::move(name), "", host.second,
                Feature::mask & ~host.second);
}

#endif // x86 || x64

} // (anonymous)

} // X86

#if NACS_CPU_X86 || NACS_CPU_X86_64
NACS_EXPORT() const CPUInfo &CPUInfo::get_host()
{
    static const auto host_info = X86::create_host_info();
    return host_info;
}
#endif // x86 || x64

} // NaCs
