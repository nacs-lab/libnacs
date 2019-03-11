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
namespace ARM {

enum class CPU : uint32_t {
    generic = 0,

    // Architecture targets
    armv7_a,
    armv7_m,
    armv7e_m,
    armv7_r,
    armv8_a,
    armv8_m_base,
    armv8_m_main,
    armv8_r,
    armv8_1_a,
    armv8_2_a,
    armv8_3_a,
    // armv8_4_a,

    // ARM
    // armv6l
    arm_mpcore,
    arm_1136jf_s,
    arm_1156t2f_s,
    arm_1176jzf_s,
    arm_cortex_m0,
    arm_cortex_m1,
    // armv7ml
    arm_cortex_m3,
    arm_cortex_m4,
    arm_cortex_m7,
    // armv7l
    arm_cortex_a5,
    arm_cortex_a7,
    arm_cortex_a8,
    arm_cortex_a9,
    arm_cortex_a12,
    arm_cortex_a15,
    arm_cortex_a17,
    arm_cortex_r4,
    arm_cortex_r5,
    arm_cortex_r7,
    arm_cortex_r8,
    // armv8ml
    arm_cortex_m23,
    arm_cortex_m33,
    // armv8l
    arm_cortex_a32,
    arm_cortex_r52,
    // aarch64
    arm_cortex_a35,
    arm_cortex_a53,
    arm_cortex_a55,
    arm_cortex_a57,
    arm_cortex_a72,
    arm_cortex_a73,
    arm_cortex_a75,

    // Cavium
    // aarch64
    cavium_thunderx,
    cavium_thunderx88,
    cavium_thunderx88p1,
    cavium_thunderx81,
    cavium_thunderx83,
    cavium_thunderx2t99,
    cavium_thunderx2t99p1,

    // NVIDIA
    // aarch64
    nvidia_denver1,
    nvidia_denver2,

    // AppliedMicro
    // aarch64
    apm_xgene1,
    apm_xgene2,
    apm_xgene3,

    // Qualcomm
    // armv7l
    qualcomm_scorpion,
    qualcomm_krait,
    // aarch64
    qualcomm_kyro,
    qualcomm_falkor,
    qualcomm_saphira,

    // Samsung
    // aarch64
    samsung_exynos_m1,
    samsung_exynos_m2,
    samsung_exynos_m3,

    // Apple
    // armv7l
    apple_swift,
    // aarch64
    apple_cyclone,
    apple_typhoon,
    apple_twister,
    apple_hurricane,

    // Marvell
    // armv7l
    marvell_pj4,

    // Intel
    // armv7l
    intel_3735d,
};

static bool is_generic_cpu_name(CPU cpu)
{
    switch (cpu) {
    case CPU::generic:
    case CPU::armv7_a:
    case CPU::armv7_m:
    case CPU::armv7e_m:
    case CPU::armv7_r:
    case CPU::armv8_a:
    case CPU::armv8_m_base:
    case CPU::armv8_m_main:
    case CPU::armv8_r:
    case CPU::armv8_1_a:
    case CPU::armv8_2_a:
    case CPU::armv8_3_a:
        return true;
    default:
        return false;
    }
}

} // ARM

namespace AArch32 {
using CPU = ARM::CPU;
static constexpr size_t feature_sz = 3;
using FeatureList = NaCs::FeatureList<feature_sz>;
using CPUSpec = NaCs::CPUSpec<CPU,feature_sz>;

namespace Feature {
static constexpr auto mask = FeatureList::get(
#define NACS_FEATURE_DEF_NAME(name, bit, llvmver, str) NACS_FEATURE_DEF(name, bit, llvmver)
#define NACS_FEATURE_DEF(name, bit, llvmver) bit,
#include "features_aarch32.h"
#undef NACS_FEATURE_DEF
#undef NACS_FEATURE_DEF_NAME
    -1);

static constexpr FeatureName names[] = {
#define NACS_FEATURE_DEF(name, bit, llvmver) {#name, bit, llvmver},
#define NACS_FEATURE_DEF_NAME(name, bit, llvmver, str) {str, bit, llvmver},
#include "features_aarch32.h"
#undef NACS_FEATURE_DEF
#undef NACS_FEATURE_DEF_NAME
};

// This does not cover all dependencies (e.g. the ones that depends on arm versions)
static constexpr FeatureDep deps[] = {
    {neon, vfp3},
    {vfp4, vfp3},
    {crypto, neon},
};

constexpr auto generic = FeatureList::get();

constexpr auto armv7m = FeatureList::get(v7, mclass, hwdiv);
constexpr auto armv7a = FeatureList::get(v7, aclass);
constexpr auto armv7r = FeatureList::get(v7, rclass);
constexpr auto armv8m = FeatureList::get(v7, v8, mclass, hwdiv);
constexpr auto armv8a = FeatureList::get(v7, v8, aclass, neon, vfp3, vfp4, d32,
                                         hwdiv, hwdiv_arm);
constexpr auto armv8r = FeatureList::get(v7, v8, rclass, neon, vfp3, vfp4, d32,
                                         hwdiv, hwdiv_arm);

// armv7l
constexpr auto arm_cortex_a5 = armv7a;
constexpr auto arm_cortex_a7 = armv7a | FeatureList::get(vfp3, vfp4, neon);
constexpr auto arm_cortex_a8 = armv7a | FeatureList::get(d32, vfp3, neon);
constexpr auto arm_cortex_a9 = armv7a;
constexpr auto arm_cortex_a12 = armv7a | FeatureList::get(d32, vfp3, vfp4, neon);
constexpr auto arm_cortex_a15 = armv7a | FeatureList::get(d32, vfp3, vfp4, neon);
constexpr auto arm_cortex_a17 = armv7a | FeatureList::get(d32, vfp3, vfp4, neon);
constexpr auto arm_cortex_r4 = armv7r | FeatureList::get(vfp3, hwdiv);
constexpr auto arm_cortex_r5 = armv7r | FeatureList::get(vfp3, hwdiv, hwdiv_arm);
constexpr auto arm_cortex_r7 = armv7r | FeatureList::get(vfp3, hwdiv, hwdiv_arm);
constexpr auto arm_cortex_r8 = armv7r | FeatureList::get(vfp3, hwdiv, hwdiv_arm);
constexpr auto qualcomm_scorpion = armv7a | FeatureList::get(v7, aclass, vfp3, neon);
constexpr auto qualcomm_krait = armv7a | FeatureList::get(vfp3, vfp4, neon, hwdiv, hwdiv_arm);
constexpr auto apple_swift = armv7a | FeatureList::get(d32, vfp3, vfp4, neon, hwdiv, hwdiv_arm);
constexpr auto marvell_pj4 = armv7a | FeatureList::get(vfp3);
constexpr auto intel_3735d = armv7a | FeatureList::get(vfp3, neon);
// armv8ml
constexpr auto arm_cortex_m23 = armv8m; // unsupported
constexpr auto arm_cortex_m33 = armv8m | FeatureList::get(v8_m_main); // unsupported
// armv8l
constexpr auto armv8a_crc = armv8a | FeatureList::get(crc);
constexpr auto armv8_1a = armv8a_crc | FeatureList::get(v8_1a);
constexpr auto armv8_2a = armv8_1a | FeatureList::get(v8_2a);
constexpr auto armv8a_crc_crypto = armv8a_crc | FeatureList::get(crypto);
constexpr auto armv8_2a_crypto = armv8_2a | FeatureList::get(crypto);
constexpr auto armv8_3a = armv8_2a | FeatureList::get(v8_3a);
constexpr auto armv8_3a_crypto = armv8_3a | FeatureList::get(crypto);

constexpr auto arm_cortex_a32 = armv8a; // TODO? (crc, crypto)
constexpr auto arm_cortex_r52 = armv8r; // TODO? (crc, crypto)
constexpr auto arm_cortex_a35 = armv8a; // TODO? (crc, crypto)
constexpr auto arm_cortex_a53 = armv8a_crc;
constexpr auto arm_cortex_a55 = armv8_2a_crypto;
constexpr auto arm_cortex_a57 = armv8a_crc;
constexpr auto arm_cortex_a72 = armv8a_crc;
constexpr auto arm_cortex_a73 = armv8a_crc;
constexpr auto arm_cortex_a75 = armv8_2a_crypto;
constexpr auto cavium_thunderx = armv8a_crc_crypto;
constexpr auto cavium_thunderx88 = armv8a_crc_crypto;
constexpr auto cavium_thunderx88p1 = armv8a_crc_crypto;
constexpr auto cavium_thunderx81 = armv8a_crc_crypto;
constexpr auto cavium_thunderx83 = armv8a_crc_crypto;
constexpr auto cavium_thunderx2t99 = armv8a_crc_crypto | FeatureList::get(v8_1a);
constexpr auto cavium_thunderx2t99p1 = armv8a_crc_crypto | FeatureList::get(v8_1a);
constexpr auto nvidia_denver1 = armv8a; // TODO? (crc, crypto)
constexpr auto nvidia_denver2 = armv8a_crc_crypto;
constexpr auto apm_xgene1 = armv8a;
constexpr auto apm_xgene2 = armv8a; // TODO?
constexpr auto apm_xgene3 = armv8a; // TODO?
constexpr auto qualcomm_kyro = armv8a_crc_crypto;
constexpr auto qualcomm_falkor = armv8a_crc_crypto;
constexpr auto qualcomm_saphira = armv8_3a_crypto;
constexpr auto samsung_exynos_m1 = armv8a_crc_crypto;
constexpr auto samsung_exynos_m2 = armv8a_crc_crypto;
constexpr auto samsung_exynos_m3 = armv8a_crc_crypto;
constexpr auto apple_cyclone = armv8a_crc_crypto;
constexpr auto apple_typhoon = armv8a_crc_crypto;
constexpr auto apple_twister = armv8a_crc_crypto;
constexpr auto apple_hurricane = armv8a_crc_crypto;

static void enable_depends(FeatureList &features)
{
    if (test_nbit(features, v8_3a))
        set_bit(features, v8_2a, true);
    if (test_nbit(features, v8_2a))
        set_bit(features, v8_1a, true);
    if (test_nbit(features, v8_1a))
        set_bit(features, crc, true);
    if (test_nbit(features, v8_1a)) {
        set_bit(features, v8, true);
        set_bit(features, aclass, true);
    }
    if (test_nbit(features, v8_m_main)) {
        set_bit(features, v8, true);
        set_bit(features, mclass, true);
    }
    if (test_nbit(features, v8)) {
        set_bit(features, v7, true);
        if (test_nbit(features, aclass)) {
            set_bit(features, neon, true);
            set_bit(features, vfp3, true);
            set_bit(features, vfp4, true);
            set_bit(features, hwdiv_arm, true);
            set_bit(features, hwdiv, true);
            set_bit(features, d32, true);
        }
    }
    NaCs::enable_depends(features, deps);
}

static void disable_depends(FeatureList &features)
{
    NaCs::disable_depends(features, deps);
}

} // Feature

namespace {

static constexpr CPUSpec cpus[] = {
    {"generic", CPU::generic, CPU::generic, 0, Feature::generic},
    // armv6
    {"mpcore", CPU::arm_mpcore, CPU::generic, 0, Feature::generic},
    {"arm1136jf-s", CPU::arm_1136jf_s, CPU::generic, 0, Feature::generic},
    {"arm1156t2f-s", CPU::arm_1156t2f_s, CPU::generic, 0, Feature::generic},
    {"arm1176jzf-s", CPU::arm_1176jzf_s, CPU::generic, 0, Feature::generic},
    {"cortex-m0", CPU::arm_cortex_m0, CPU::generic, 0, Feature::generic},
    {"cortex-m1", CPU::arm_cortex_m1, CPU::generic, 0, Feature::generic},
    // armv7ml
    {"armv7-m", CPU::armv7_m, CPU::generic, 0, Feature::armv7m},
    {"armv7e-m", CPU::armv7e_m, CPU::generic, 0, Feature::armv7m},
    {"cortex-m3", CPU::arm_cortex_m3, CPU::generic, 0, Feature::armv7m},
    {"cortex-m4", CPU::arm_cortex_m4, CPU::generic, 0, Feature::armv7m},
    {"cortex-m7", CPU::arm_cortex_m7, CPU::generic, 0, Feature::armv7m},
    // armv7l
    {"armv7-a", CPU::armv7_a, CPU::generic, 0, Feature::armv7a},
    {"armv7-r", CPU::armv7_r, CPU::generic, 0, Feature::armv7r},
    {"cortex-a5", CPU::arm_cortex_a5, CPU::generic, 0, Feature::arm_cortex_a5},
    {"cortex-a7", CPU::arm_cortex_a7, CPU::generic, 0, Feature::arm_cortex_a7},
    {"cortex-a8", CPU::arm_cortex_a8, CPU::generic, 0, Feature::arm_cortex_a8},
    {"cortex-a9", CPU::arm_cortex_a9, CPU::generic, 0, Feature::arm_cortex_a9},
    {"cortex-a12", CPU::arm_cortex_a12, CPU::generic, 0, Feature::arm_cortex_a12},
    {"cortex-a15", CPU::arm_cortex_a15, CPU::generic, 0, Feature::arm_cortex_a15},
    {"cortex-a17", CPU::arm_cortex_a17, CPU::generic, 0, Feature::arm_cortex_a17},
    {"cortex-r4", CPU::arm_cortex_r4, CPU::generic, 0, Feature::arm_cortex_r4},
    {"cortex-r5", CPU::arm_cortex_r5, CPU::generic, 0, Feature::arm_cortex_r5},
    {"cortex-r7", CPU::arm_cortex_r7, CPU::generic, 0, Feature::arm_cortex_r7},
    {"cortex-r8", CPU::arm_cortex_r8, CPU::generic, 0, Feature::arm_cortex_r8},
    {"scorpion", CPU::qualcomm_scorpion, CPU::armv7_a, UINT32_MAX, Feature::qualcomm_scorpion},
    {"krait", CPU::qualcomm_krait, CPU::generic, 0, Feature::qualcomm_krait},
    {"swift", CPU::apple_swift, CPU::generic, 0, Feature::apple_swift},
    {"pj4", CPU::marvell_pj4, CPU::armv7_a, UINT32_MAX, Feature::marvell_pj4},
    {"3735d", CPU::intel_3735d, CPU::armv7_a, UINT32_MAX, Feature::intel_3735d},

    // armv8ml
    {"armv8-m.base", CPU::armv8_m_base, CPU::generic, 0, Feature::armv8m},
    {"armv8-m.main", CPU::armv8_m_main, CPU::generic, 0, Feature::armv8m},
    {"cortex-m23", CPU::arm_cortex_m23, CPU::generic, 0, Feature::arm_cortex_m23},
    {"cortex-m33", CPU::arm_cortex_m33, CPU::generic, 0, Feature::arm_cortex_m33},

    // armv8l
    {"armv8-a", CPU::armv8_a, CPU::generic, 0, Feature::armv8a},
    {"armv8-r", CPU::armv8_r, CPU::generic, 0, Feature::armv8r},
    {"armv8.1-a", CPU::armv8_1_a, CPU::generic, 0, Feature::armv8_1a},
    {"armv8.2-a", CPU::armv8_2_a, CPU::generic, 0, Feature::armv8_2a},
    {"armv8.3-a", CPU::armv8_3_a, CPU::generic, 0, Feature::armv8_3a},
    {"cortex-a32", CPU::arm_cortex_a32, CPU::generic, 0, Feature::arm_cortex_a32},
    {"cortex-r52", CPU::arm_cortex_r52, CPU::generic, 0, Feature::arm_cortex_r52},
    {"cortex-a35", CPU::arm_cortex_a35, CPU::generic, 0, Feature::arm_cortex_a35},
    {"cortex-a53", CPU::arm_cortex_a53, CPU::generic, 0, Feature::arm_cortex_a53},
    {"cortex-a55", CPU::arm_cortex_a55, CPU::arm_cortex_a53, 60000, Feature::arm_cortex_a55},
    {"cortex-a57", CPU::arm_cortex_a57, CPU::generic, 0, Feature::arm_cortex_a57},
    {"cortex-a72", CPU::arm_cortex_a72, CPU::generic, 0, Feature::arm_cortex_a72},
    {"cortex-a73", CPU::arm_cortex_a73, CPU::generic, 0, Feature::arm_cortex_a73},
    {"cortex-a75", CPU::arm_cortex_a75, CPU::arm_cortex_a73, 60000, Feature::arm_cortex_a75},
    {"thunderx", CPU::cavium_thunderx, CPU::armv8_a, UINT32_MAX, Feature::cavium_thunderx},
    {"thunderx88", CPU::cavium_thunderx88, CPU::armv8_a, UINT32_MAX, Feature::cavium_thunderx88},
    {"thunderx88p1", CPU::cavium_thunderx88p1, CPU::armv8_a, UINT32_MAX,
     Feature::cavium_thunderx88p1},
    {"thunderx81", CPU::cavium_thunderx81, CPU::armv8_a, UINT32_MAX,
     Feature::cavium_thunderx81},
    {"thunderx83", CPU::cavium_thunderx83, CPU::armv8_a, UINT32_MAX,
     Feature::cavium_thunderx83},
    {"thunderx2t99", CPU::cavium_thunderx2t99, CPU::armv8_a, UINT32_MAX,
     Feature::cavium_thunderx2t99},
    {"thunderx2t99p1", CPU::cavium_thunderx2t99p1, CPU::armv8_a, UINT32_MAX,
     Feature::cavium_thunderx2t99p1},
    {"denver1", CPU::nvidia_denver1, CPU::arm_cortex_a53, UINT32_MAX, Feature::nvidia_denver1},
    {"denver2", CPU::nvidia_denver2, CPU::arm_cortex_a57, UINT32_MAX, Feature::nvidia_denver2},
    {"xgene1", CPU::apm_xgene1, CPU::armv8_a, UINT32_MAX, Feature::apm_xgene1},
    {"xgene2", CPU::apm_xgene2, CPU::armv8_a, UINT32_MAX, Feature::apm_xgene2},
    {"xgene3", CPU::apm_xgene3, CPU::armv8_a, UINT32_MAX, Feature::apm_xgene3},
    {"kyro", CPU::qualcomm_kyro, CPU::armv8_a, UINT32_MAX, Feature::qualcomm_kyro},
    {"falkor", CPU::qualcomm_falkor, CPU::armv8_a, UINT32_MAX, Feature::qualcomm_falkor},
    {"saphira", CPU::qualcomm_saphira, CPU::armv8_a, UINT32_MAX, Feature::qualcomm_saphira},
    {"exynos-m1", CPU::samsung_exynos_m1, CPU::generic, 0, Feature::samsung_exynos_m1},
    {"exynos-m2", CPU::samsung_exynos_m2, CPU::generic, 0, Feature::samsung_exynos_m2},
    {"exynos-m3", CPU::samsung_exynos_m3, CPU::generic, 0, Feature::samsung_exynos_m3},
    {"cyclone", CPU::apple_cyclone, CPU::generic, 0, Feature::apple_cyclone},
    {"typhoon", CPU::apple_typhoon, CPU::apple_cyclone, UINT32_MAX, Feature::apple_typhoon},
    {"twister", CPU::apple_twister, CPU::apple_typhoon, UINT32_MAX, Feature::apple_twister},
    {"hurricane", CPU::apple_hurricane, CPU::apple_twister, UINT32_MAX, Feature::apple_hurricane},
};

struct Info : CPUInfo {
    Info(std::string name, std::string ext_features, FeatureList en, FeatureList dis)
        : CPUInfo(std::move(name), std::move(ext_features)),
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
        Feature::enable_depends(m_en);
        m_en = m_en & ~m_dis;
        Feature::disable_depends(m_en);
        if (spec) {
            // If the base feature is known, disable everything that we didn't
            // enable explicitly.
            m_dis = Feature::mask & ~m_en;
        }
    }
    bool test_feature(int bit) const final override
    {
        return test_nbit(m_en, bit);
    }
    int get_vector_size() const final override
    {
        if (test_feature(Feature::neon))
            return 16;
        return 8;
    }
    const std::string &get_arch() const final override
    {
        static const std::string name = "arm";
        return name;
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
    auto en_features = m_en;
    if (spec && is_generic_cpu_name(spec->cpu)) {
        llvm_name = "generic";
        en_features = en_features | (spec->features & ~m_dis);
    }
    std::vector<std::string> features;
    for (auto &fename: Feature::names) {
        if (fename.llvmver > llvmver)
            continue;
        bool enable = test_nbit(en_features, fename.bit);
        bool disable = test_nbit(m_dis, fename.bit);
        if (fename.bit == Feature::d32) {
            if (enable) {
                features.push_back("-d16");
            }
            else if (disable) {
                features.push_back("+d16");
            }
            continue;
        }
        if (enable) {
            features.insert(features.begin(), std::string("+") + fename.name);
        }
        else if (disable) {
            features.push_back(std::string("-") + fename.name);
        }
    }
    if (test_nbit(en_features, Feature::v8_2a))
        features.push_back("+v8.2a");
    if (test_nbit(en_features, Feature::v8_1a))
        features.push_back("+v8.1a");
    if (test_nbit(en_features, Feature::v8_m_main)) {
        features.push_back("+v8m.main");
        features.push_back("+armv8-m.main");
    }
    if (test_nbit(en_features, Feature::aclass))
        features.push_back("+aclass");
    if (test_nbit(en_features, Feature::rclass))
        features.push_back("+rclass");
    if (test_nbit(en_features, Feature::mclass))
        features.push_back("+mclass");
    if (test_nbit(en_features, Feature::v8)) {
        features.push_back("+v8");
        if (test_nbit(en_features, Feature::aclass))
            features.push_back("+armv8-a");
        if (test_nbit(en_features, Feature::rclass))
            features.push_back("+armv8-r");
        if (test_nbit(en_features, Feature::mclass)) {
            features.push_back("+v8m");
            features.push_back("+armv8-m.base");
        }
    }
    if (test_nbit(en_features, Feature::v7)) {
        features.push_back("+v7");
        if (test_nbit(en_features, Feature::aclass))
            features.push_back("+armv7-a");
        if (test_nbit(en_features, Feature::rclass))
            features.push_back("+armv7-r");
        if (test_nbit(en_features, Feature::mclass))
            features.push_back("+armv7-m");
    }
    features.push_back("+v6");
    features.push_back("+vfp2");
    append_features(features, ext_features);
    return std::make_pair(std::move(llvm_name), std::move(features));
}

} // (anonymous)

} // AArch32

namespace AArch64 {
using CPU = ARM::CPU;
static constexpr size_t feature_sz = 3;
using FeatureList = NaCs::FeatureList<feature_sz>;
using CPUSpec = NaCs::CPUSpec<CPU,feature_sz>;

namespace Feature {
static constexpr auto mask = FeatureList::get(
#define NACS_FEATURE_DEF_NAME(name, bit, llvmver, str) NACS_FEATURE_DEF(name, bit, llvmver)
#define NACS_FEATURE_DEF(name, bit, llvmver) bit,
#include "features_aarch64.h"
#undef NACS_FEATURE_DEF
#undef NACS_FEATURE_DEF_NAME
    -1);

static constexpr FeatureName names[] = {
#define NACS_FEATURE_DEF(name, bit, llvmver) {#name, bit, llvmver},
#define NACS_FEATURE_DEF_NAME(name, bit, llvmver, str) {str, bit, llvmver},
#include "features_aarch64.h"
#undef NACS_FEATURE_DEF
#undef NACS_FEATURE_DEF_NAME
};

// This does not cover all dependencies (e.g. the ones that depends on arm versions)
static constexpr FeatureDep deps[] = {};

constexpr auto generic = FeatureList::get();
constexpr auto armv8a_crc = FeatureList::get(crc);
constexpr auto armv8a_crc_crypto = armv8a_crc | FeatureList::get(crypto);
constexpr auto armv8_1a = armv8a_crc | FeatureList::get(v8_1a, lse, rdm); // lor, hpd
constexpr auto armv8_2a = armv8_1a | FeatureList::get(v8_2a); // ras
constexpr auto armv8_2a_crypto = armv8_2a | FeatureList::get(crypto);
constexpr auto armv8_3a = armv8_2a | FeatureList::get(v8_3a, rcpc);
constexpr auto armv8_3a_crypto = armv8_3a | FeatureList::get(crypto);

constexpr auto arm_cortex_a32 = generic; // TODO? (crc, crypto)
constexpr auto arm_cortex_a35 = generic; // TODO? (crc, crypto)
constexpr auto arm_cortex_a53 = armv8a_crc;
constexpr auto arm_cortex_a55 = armv8_2a_crypto | FeatureList::get(rcpc); // dotprod;
constexpr auto arm_cortex_a57 = armv8a_crc;
constexpr auto arm_cortex_a72 = armv8a_crc;
constexpr auto arm_cortex_a73 = armv8a_crc;
constexpr auto arm_cortex_a75 = armv8_2a_crypto | FeatureList::get(rcpc); // dotprod;
constexpr auto cavium_thunderx = armv8a_crc_crypto;
constexpr auto cavium_thunderx88 = armv8a_crc_crypto;
constexpr auto cavium_thunderx88p1 = armv8a_crc_crypto;
constexpr auto cavium_thunderx81 = armv8a_crc_crypto;
constexpr auto cavium_thunderx83 = armv8a_crc_crypto;
constexpr auto cavium_thunderx2t99 = armv8a_crc_crypto | FeatureList::get(v8_1a);
constexpr auto cavium_thunderx2t99p1 = armv8a_crc_crypto | FeatureList::get(v8_1a);
constexpr auto nvidia_denver1 = generic; // TODO? (crc, crypto)
constexpr auto nvidia_denver2 = armv8a_crc_crypto;
constexpr auto apm_xgene1 = generic;
constexpr auto apm_xgene2 = generic; // TODO?
constexpr auto apm_xgene3 = generic; // TODO?
constexpr auto qualcomm_kyro = armv8a_crc_crypto;
constexpr auto qualcomm_falkor = armv8a_crc_crypto;
constexpr auto qualcomm_saphira = armv8_3a_crypto;
constexpr auto samsung_exynos_m1 = armv8a_crc_crypto;
constexpr auto samsung_exynos_m2 = armv8a_crc_crypto;
constexpr auto samsung_exynos_m3 = armv8a_crc_crypto;
constexpr auto apple_cyclone = armv8a_crc_crypto;
constexpr auto apple_typhoon = armv8a_crc_crypto;
constexpr auto apple_twister = armv8a_crc_crypto;
constexpr auto apple_hurricane = armv8a_crc_crypto;

static void enable_depends(FeatureList &features)
{
    if (test_nbit(features, v8_3a))
        set_bit(features, v8_2a, true);
    if (test_nbit(features, v8_2a))
        set_bit(features, v8_1a, true);
    if (test_nbit(features, v8_1a))
        set_bit(features, crc, true);
    if (test_nbit(features, v8_1a)) {
        set_bit(features, lse, true);
        set_bit(features, rdm, true);
    }
    NaCs::enable_depends(features, deps);
}

static void disable_depends(FeatureList &features)
{
    NaCs::disable_depends(features, deps);
}

} // Feature

namespace {

static constexpr CPUSpec cpus[] = {
    {"generic", CPU::generic, CPU::generic, 0, Feature::generic},
    {"armv8.1-a", CPU::armv8_1_a, CPU::generic, 0, Feature::armv8_1a},
    {"armv8.2-a", CPU::armv8_2_a, CPU::generic, 0, Feature::armv8_2a},
    {"armv8.3_a", CPU::armv8_3_a, CPU::generic, 0, Feature::armv8_3a},
    {"cortex-a35", CPU::arm_cortex_a35, CPU::generic, 0, Feature::arm_cortex_a35},
    {"cortex-a53", CPU::arm_cortex_a53, CPU::generic, 0, Feature::arm_cortex_a53},
    {"cortex-a55", CPU::arm_cortex_a55, CPU::arm_cortex_a53, UINT32_MAX, Feature::arm_cortex_a55},
    {"cortex-a57", CPU::arm_cortex_a57, CPU::generic, 0, Feature::arm_cortex_a57},
    {"cortex-a72", CPU::arm_cortex_a72, CPU::generic, 0, Feature::arm_cortex_a72},
    {"cortex-a73", CPU::arm_cortex_a73, CPU::generic, 0, Feature::arm_cortex_a73},
    {"cortex-a75", CPU::arm_cortex_a75, CPU::arm_cortex_a73, UINT32_MAX, Feature::arm_cortex_a75},
    {"thunderx", CPU::cavium_thunderx, CPU::generic, 0, Feature::cavium_thunderx},
    {"thunderxt88", CPU::cavium_thunderx88, CPU::generic, 0, Feature::cavium_thunderx88},
    {"thunderxt88p1", CPU::cavium_thunderx88p1, CPU::cavium_thunderx88, UINT32_MAX,
     Feature::cavium_thunderx88p1},
    {"thunderxt81", CPU::cavium_thunderx81, CPU::generic, 0, Feature::cavium_thunderx81},
    {"thunderxt83", CPU::cavium_thunderx83, CPU::generic, 0, Feature::cavium_thunderx83},
    {"thunderx2t99", CPU::cavium_thunderx2t99, CPU::generic, 0, Feature::cavium_thunderx2t99},
    {"thunderx2t99p1", CPU::cavium_thunderx2t99p1, CPU::cavium_thunderx2t99, UINT32_MAX,
     Feature::cavium_thunderx2t99p1},
    {"denver1", CPU::nvidia_denver1, CPU::generic, UINT32_MAX, Feature::nvidia_denver1},
    {"denver2", CPU::nvidia_denver2, CPU::generic, UINT32_MAX, Feature::nvidia_denver2},
    {"xgene1", CPU::apm_xgene1, CPU::generic, UINT32_MAX, Feature::apm_xgene1},
    {"xgene2", CPU::apm_xgene2, CPU::generic, UINT32_MAX, Feature::apm_xgene2},
    {"xgene3", CPU::apm_xgene3, CPU::generic, UINT32_MAX, Feature::apm_xgene3},
    {"kyro", CPU::qualcomm_kyro, CPU::generic, 0, Feature::qualcomm_kyro},
    {"falkor", CPU::qualcomm_falkor, CPU::generic, 0, Feature::qualcomm_falkor},
    {"saphira", CPU::qualcomm_saphira, CPU::qualcomm_falkor, 60000, Feature::qualcomm_saphira},
    {"exynos-m1", CPU::samsung_exynos_m1, CPU::generic, 0, Feature::samsung_exynos_m1},
    {"exynos-m2", CPU::samsung_exynos_m2, CPU::generic, 0, Feature::samsung_exynos_m2},
    {"exynos-m3", CPU::samsung_exynos_m3, CPU::generic, 0, Feature::samsung_exynos_m3},
    {"cyclone", CPU::apple_cyclone, CPU::generic, 0, Feature::apple_cyclone},
    {"typhoon", CPU::apple_typhoon, CPU::apple_cyclone, UINT32_MAX, Feature::apple_typhoon},
    {"twister", CPU::apple_twister, CPU::apple_typhoon, UINT32_MAX, Feature::apple_twister},
    {"hurricane", CPU::apple_hurricane, CPU::apple_twister, UINT32_MAX, Feature::apple_hurricane},
};

struct Info : CPUInfo {
    Info(std::string name, std::string ext_features, FeatureList en, FeatureList dis)
        : CPUInfo(std::move(name), std::move(ext_features)),
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
        Feature::enable_depends(m_en);
        m_en = m_en & ~m_dis;
        Feature::disable_depends(m_en);
        if (spec) {
            // If the base feature is known, disable everything that we didn't
            // enable explicitly.
            m_dis = Feature::mask & ~m_en;
        }
    }
    bool test_feature(int bit) const final override
    {
        return test_nbit(m_en, bit);
    }
    int get_vector_size() const final override
    {
        return 16;
    }
    const std::string &get_arch() const final override
    {
        static const std::string name = "aarch64";
        return name;
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
    auto en_features = m_en;
    if (spec && is_generic_cpu_name(spec->cpu)) {
        llvm_name = "generic";
        en_features = en_features | (spec->features & ~m_dis);
    }
    std::vector<std::string> features;
    for (auto &fename: Feature::names) {
        if (fename.llvmver > llvmver)
            continue;
        if (test_nbit(en_features, fename.bit)) {
            features.insert(features.begin(), std::string("+") + fename.name);
        }
        else if (test_nbit(m_dis, fename.bit)) {
            features.push_back(std::string("-") + fename.name);
        }
    }
    if (test_nbit(en_features, Feature::v8_2a))
        features.push_back("+v8.2a");
    if (test_nbit(en_features, Feature::v8_1a))
        features.push_back("+v8.1a");
    features.push_back("+neon");
    features.push_back("+fp-armv8");
    append_features(features, ext_features);
    return std::make_pair(std::move(llvm_name), std::move(features));
}

} // (anonymous)

} // AArch64

} // NaCs

#if NACS_CPU_AARCH32 || NACS_CPU_AARCH64

#include <dirent.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#if NACS_CPU_AARCH64 || __GLIBC_PREREQ(2, 16)
#  include <sys/auxv.h>
#else
#  define DYN_GETAUXVAL
#  include "dlload.h"
#endif

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <set>

namespace NaCs {

namespace ARM {

namespace {

namespace AArch {
#if NACS_CPU_AARCH64
using namespace AArch64;
#else
using namespace AArch32;
#endif
}

// auxval reader

#ifndef AT_HWCAP
#  define AT_HWCAP 16
#endif
#ifndef AT_HWCAP2
#  define AT_HWCAP2 26
#endif

#if defined(DYN_GETAUXVAL)
static bool getauxval_dlsym(unsigned long type, unsigned long *val)
{
    static auto getauxval_p = (unsigned long (*)(unsigned long))DL::sym(nullptr, "getauxval");
    if (getauxval_p) {
        *val = getauxval_p(type);
        return true;
    }
    return false;
}

static unsigned long getauxval_procfs(unsigned long type)
{
    int fd = open("/proc/self/auxv", O_RDONLY);
    if (fd == -1)
        return 0;
    unsigned long val = 0;
    unsigned long buff[2];
    while (read(fd, buff, sizeof(buff)) == sizeof(buff)) {
        if (buff[0] == 0)
            break;
        if (buff[0] == type) {
            val = buff[1];
            break;
        }
    }
    close(fd);
    return val;
}

static inline unsigned long getauxval(unsigned long type)
{
    unsigned long val;
    if (getauxval_dlsym(type, &val))
        return val;
    return getauxval_procfs(type);
}
#endif

struct CPUID {
    uint8_t implementer;
    uint8_t variant;
    uint16_t part;
    bool operator<(const CPUID &right) const
    {
        if (implementer < right.implementer)
            return true;
        if (implementer > right.implementer)
            return false;
        if (part < right.part)
            return true;
        if (part > right.part)
            return false;
        return variant < right.variant;
    }
};

// /sys/devices/system/cpu/cpu<n>/regs/identification/midr_el1 reader
static inline void get_cpuinfo_sysfs(std::set<CPUID> &res)
{
    // This only works on a 64bit 4.7+ kernel
    auto dir = opendir("/sys/devices/system/cpu");
    if (!dir)
        return;
    while (auto entry = readdir(dir)) {
        if (entry->d_type != DT_DIR)
            continue;
        if (strncmp(entry->d_name, "cpu", 3) != 0)
            continue;
        string_ostream stm;
        stm << "/sys/devices/system/cpu/" << entry->d_name << "/regs/identification/midr_el1";
        std::ifstream file(stm.get_buf());
        if (!file)
            continue;
        uint64_t val = 0;
        file >> std::hex >> val;
        if (!file)
            continue;
        CPUID cpuid = {
            uint8_t(val >> 24),
            uint8_t((val >> 20) & 0xf),
            uint16_t((val >> 4) & 0xfff)
        };
        res.insert(cpuid);
    }
    closedir(dir);
}

// Use an external template since lambda's can't be templated in C++11
template<typename T, typename F>
static inline bool try_read_procfs_line(const std::string &line,
                                        const char *prefix, T &out,
                                        bool &flag, F &&reset)
{
    size_t prefix_len = strlen(prefix);
    if (line.size() < prefix_len)
        return false;
    if (memcmp(&line[0], prefix, prefix_len) != 0)
        return false;
    if (flag)
        reset();
    const char *p = &line[prefix_len];
    while (*p == '\t' || *p == ' ' || *p == ':')
        p++;
    char *str_end;
    auto num = std::strtoull(p, &str_end, 0);
    out = (T)num;
    if (str_end == p) {
        flag = false;
    }
    else if (num > (unsigned long long)std::numeric_limits<T>::max()) {
        flag = false;
    }
    else {
        flag = true;
    }
    return true;
}

// /proc/cpuinfo reader
static inline void get_cpuinfo_procfs(std::set<CPUID> &res)
{
    std::ifstream file("/proc/cpuinfo");
    CPUID cpuid = {0, 0, 0};
    bool impl = false;
    bool part = false;
    bool var = false;
    auto reset = [&] () {
        if (impl && part)
            res.insert(cpuid);
        impl = false;
        part = false;
        var = false;
        memset(&cpuid, 0, sizeof(cpuid));
    };
    for (std::string line; std::getline(file, line);) {
        if (line.empty()) {
            reset();
            continue;
        }
        try_read_procfs_line(line, "CPU implementer", cpuid.implementer, impl, reset) ||
            try_read_procfs_line(line, "CPU variant", cpuid.variant, var, reset) ||
            try_read_procfs_line(line, "CPU part", cpuid.part, part, reset);
    }
    reset();
}

static std::set<CPUID> get_cpuinfo(void)
{
    std::set<CPUID> res;
    get_cpuinfo_sysfs(res);
    if (res.empty())
        get_cpuinfo_procfs(res);
    return res;
}

static CPU get_cpu_name(CPUID cpuid)
{
    switch (cpuid.implementer) {
    case 0x41: // ARM
        switch (cpuid.part) {
        case 0xb02: return CPU::arm_mpcore;
        case 0xb36: return CPU::arm_1136jf_s;
        case 0xb56: return CPU::arm_1156t2f_s;
        case 0xb76: return CPU::arm_1176jzf_s;
        case 0xc20: return CPU::arm_cortex_m0;
        case 0xc21: return CPU::arm_cortex_m1;
        case 0xc23: return CPU::arm_cortex_m3;
        case 0xc24: return CPU::arm_cortex_m4;
        case 0xc27: return CPU::arm_cortex_m7;
        case 0xd20: return CPU::arm_cortex_m23;
        case 0xd21: return CPU::arm_cortex_m33;
        case 0xc05: return CPU::arm_cortex_a5;
        case 0xc07: return CPU::arm_cortex_a7;
        case 0xc08: return CPU::arm_cortex_a8;
        case 0xc09: return CPU::arm_cortex_a9;
        case 0xc0d: return CPU::arm_cortex_a12;
        case 0xc0f: return CPU::arm_cortex_a15;
        case 0xc0e: return CPU::arm_cortex_a17;
        case 0xc14: return CPU::arm_cortex_r4;
        case 0xc15: return CPU::arm_cortex_r5;
        case 0xc17: return CPU::arm_cortex_r7;
        case 0xc18: return CPU::arm_cortex_r8;
        case 0xd13: return CPU::arm_cortex_r52;
        case 0xd01: return CPU::arm_cortex_a32;
        case 0xd04: return CPU::arm_cortex_a35;
        case 0xd03: return CPU::arm_cortex_a53;
        case 0xd05: return CPU::arm_cortex_a55;
        case 0xd07: return CPU::arm_cortex_a57;
        case 0xd08: return CPU::arm_cortex_a72;
        case 0xd09: return CPU::arm_cortex_a73;
        case 0xd0a: return CPU::arm_cortex_a75;
        default: return CPU::generic;
        }
    case 0x42: // Broadcom (Cavium)
        switch (cpuid.part) {
        case 0x516: return CPU::cavium_thunderx2t99p1;
        default: return CPU::generic;
        }
    case 0x43: // Cavium
        switch (cpuid.part) {
        case 0xa0: return CPU::cavium_thunderx;
        case 0xa1:
            if (cpuid.variant == 0)
                return CPU::cavium_thunderx88p1;
            return CPU::cavium_thunderx88;
        case 0xa2: return CPU::cavium_thunderx81;
        case 0xa3: return CPU::cavium_thunderx83;
        case 0xaf: return CPU::cavium_thunderx2t99;
        default: return CPU::generic;
        }
    case 0x4e: // NVIDIA
        switch (cpuid.part) {
        case 0x000: return CPU::nvidia_denver1;
        case 0x003: return CPU::nvidia_denver2;
        default: return CPU::generic;
        }
    case 0x50: // AppliedMicro
        // x-gene 2
        // x-gene 3
        switch (cpuid.part) {
        case 0x000: return CPU::apm_xgene1;
        default: return CPU::generic;
        }
    case 0x51: // Qualcomm
        switch (cpuid.part) {
        case 0x00f:
        case 0x02d:
            return CPU::qualcomm_scorpion;
        case 0x04d:
        case 0x06f:
            return CPU::qualcomm_krait;
        case 0x201:
        case 0x205:
        case 0x211:
            return CPU::qualcomm_kyro;
        case 0x800:
        case 0x801:
            return CPU::arm_cortex_a73; // second-generation Kryo
        case 0xc00:
            return CPU::qualcomm_falkor;
        case 0xc01:
            return CPU::qualcomm_saphira;
        default: return CPU::generic;
        }
    case 0x53: // Samsung
        // exynos-m2
        // exynos-m3
        switch (cpuid.part) {
        case 0x001: return CPU::samsung_exynos_m1;
        default: return CPU::generic;
        }
    case 0x56: // Marvell
        switch (cpuid.part) {
        case 0x581:
        case 0x584:
            return CPU::marvell_pj4;
        default: return CPU::generic;
        }
    case 0x67: // Apple
        // swift
        // cyclone
        // twister
        // hurricane
        switch (cpuid.part) {
        case 0x072: return CPU::apple_typhoon;
        default: return CPU::generic;
        }
    case 0x69: // Intel
        switch (cpuid.part) {
        case 0x001: return CPU::intel_3735d;
        default: return CPU::generic;
        }
    default:
        return CPU::generic;
    }
}

static std::pair<int,char> get_elf_arch(void)
{
#if NACS_CPU_AARCH64
    return std::make_pair(8, 'A');
#else
    int ver = 0;
    char profile = 0;
    struct utsname name;
    if (uname(&name) >= 0) {
        // name.machine is the elf_platform in the kernel.
        if (strcmp(name.machine, "armv6l") == 0) {
            ver = 6;
        }
        else if (strcmp(name.machine, "armv7l") == 0) {
            ver = 7;
        }
        else if (strcmp(name.machine, "armv7ml") == 0) {
            ver = 7;
            profile = 'M';
        }
        else if (strcmp(name.machine, "armv8l") == 0 || strcmp(name.machine, "aarch64") == 0) {
            ver = 8;
        }
    }
    if (__ARM_ARCH > ver)
        ver = __ARM_ARCH;
#  if __ARM_ARCH > 6 && defined(__ARM_ARCH_PROFILE)
    profile = __ARM_ARCH_PROFILE;
#  endif
    return std::make_pair(ver, profile);
#endif
}

static std::pair<int,bool> feature_arch_version(const AArch::FeatureList &feature)
{
#if NACS_CPU_AARCH64
    (void)feature;
    return std::make_pair(8, false);
#else
    if (test_nbit(feature, AArch::Feature::v8))
        return std::make_pair(8, test_nbit(feature, AArch::Feature::mclass));
    if (test_nbit(feature, AArch::Feature::v7))
        return std::make_pair(7, test_nbit(feature, AArch::Feature::mclass));
    return std::make_pair(6, false);
#endif
}

static CPU generic_for_arch(std::pair<int,bool> arch)
{
#if NACS_CPU_AARCH64
    (void)arch;
    return CPU::generic;
#else
#  if defined(__ARM_ARCH_PROFILE)
    char klass = __ARM_ARCH_PROFILE;
#  else
    char klass = arch.second ? 'M' : 'A';
#  endif
    if (arch.first >= 8) {
        if (klass == 'M') {
            return CPU::armv8_m_base;
        }
        else if (klass == 'R') {
            return CPU::armv8_r;
        }
        else {
            return CPU::armv8_a;
        }
    }
    else if (arch.first == 7) {
        if (klass == 'M') {
            return CPU::armv7_m;
        }
        else if (klass == 'R') {
            return CPU::armv7_r;
        }
        else {
            return CPU::armv7_a;
        }
    }
    return CPU::generic;
#endif
}

static bool check_cpu_arch_ver(uint32_t cpu, std::pair<int,bool> arch)
{
    auto spec = find_cpu(cpu, AArch::cpus);
    // This happens on AArch64 and indicates that the cpu name isn't a valid aarch64 CPU
    if (!spec)
        return false;
    auto cpu_arch = feature_arch_version(spec->features);
    if (arch.second != cpu_arch.second)
        return false;
    if (arch.first > cpu_arch.first)
        return false;
    return true;
}

template<size_t ncpu>
static void shrink_big_little(std::vector<std::pair<uint32_t,CPUID>> &list,
                              const CPU (&cpus)[ncpu])
{
    auto find = [&] (uint32_t name) {
        for (uint32_t i = 0; i < ncpu; i++) {
            if (cpus[i] == CPU(name)) {
                return (int)i;
            }
        }
        return -1;
    };
    int maxidx = -1;
    for (auto &ele: list) {
        int idx = find(ele.first);
        if (idx > maxidx) {
            maxidx = idx;
        }
    }
    if (maxidx >= 0) {
        list.erase(std::remove_if(list.begin(), list.end(), [&] (std::pair<uint32_t,CPUID> &ele) {
                    int idx = find(ele.first);
                    return idx != -1 && idx < maxidx;
                }), list.end());
    }
}

static std::pair<CPU,AArch::FeatureList> _get_host_cpu()
{
    AArch::FeatureList features = {};
    // Here we assume that only the lower 32bit are used on aarch64
    // Change the cast here when that's not the case anymore (and when there's features in the
    // high bits that we want to detect).
    features[0] = (uint32_t)getauxval(AT_HWCAP);
    features[1] = (uint32_t)getauxval(AT_HWCAP2);
    auto cpuinfo = get_cpuinfo();
    auto arch = get_elf_arch();
#if NACS_CPU_AARCH32
    if (arch.first >= 7) {
        if (arch.second == 'M') {
            set_bit(features, AArch::Feature::mclass, true);
        }
        else if (arch.second == 'R') {
            set_bit(features, AArch::Feature::rclass, true);
        }
        else if (arch.second == 'A') {
            set_bit(features, AArch::Feature::aclass, true);
        }
    }
    switch (arch.first) {
    case 8:
        set_bit(features, AArch::Feature::v8, true);
        /* FALLTHROUGH */
    case 7:
        set_bit(features, AArch::Feature::v7, true);
        break;
    default:
        break;
    }
#endif

    std::set<uint32_t> cpus;
    std::vector<std::pair<uint32_t,CPUID>> list;
    for (auto info: cpuinfo) {
        auto name = (uint32_t)get_cpu_name(info);
        if (name == 0)
            continue;
        if (!check_cpu_arch_ver(name, arch))
            continue;
        if (cpus.insert(name).second) {
            features = features | find_cpu(name, AArch::cpus)->features;
            list.emplace_back(name, info);
        }
    }
    // Not all elements/pairs are valid
    static constexpr CPU v8order[] = {
        CPU::arm_cortex_a32,
        CPU::arm_cortex_a35,
        CPU::arm_cortex_a53,
        CPU::arm_cortex_a55,
        CPU::arm_cortex_a57,
        CPU::arm_cortex_a72,
        CPU::arm_cortex_a73,
        CPU::arm_cortex_a75,
        CPU::nvidia_denver2,
        CPU::samsung_exynos_m1
    };
    shrink_big_little(list, v8order);
#if NACS_CPU_AARCH32
    // Not all elements/pairs are valid
    static constexpr CPU v7order[] = {
        CPU::arm_cortex_a5,
        CPU::arm_cortex_a7,
        CPU::arm_cortex_a8,
        CPU::arm_cortex_a9,
        CPU::arm_cortex_a12,
        CPU::arm_cortex_a15,
        CPU::arm_cortex_a17
    };
    shrink_big_little(list, v7order);
#endif
    CPU cpu;
    if (list.empty()) {
        cpu = generic_for_arch(arch);
    }
    else {
        // This also covers `list.size() > 1` case which means there's a unknown combination
        // consists of CPU's we know. Unclear what else we could try so just randomly return
        // one...
        cpu = (CPU)list[0].first;
    }
    // Ignore feature bits that we are not interested in.
    features = features & AArch::Feature::mask;

    return std::make_pair(cpu, features);
}

static AArch::Info create_host_info()
{
    auto host = _get_host_cpu();
    std::string name = is_generic_cpu_name(host.first) ?
        std::string(find_cpu_name(uint32_t(host.first), AArch::cpus)) :
        LLVM::get_cpu_name();
    return AArch::Info(std::move(name), "", host.second,
                       AArch::Feature::mask & ~host.second);
}

} // (anonymous)

} // ARM

NACS_EXPORT() const CPUInfo &CPUInfo::get_host()
{
    static const auto host_info = ARM::create_host_info();
    return host_info;
}

} // NaCs
#endif // aarch32 || aarch64
