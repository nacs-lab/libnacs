/*************************************************************************
 *   Copyright (c) 2018 - 2019 Yichao Yu <yyc1992@gmail.com>             *
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

#include "cpu_p.h"

namespace NaCs {
namespace LLVM {

const std::string &get_cpu_arch()
{
#if NACS_CPU_X86
    static const std::string arch = "i386";
#elif NACS_CPU_X86_64
    static const std::string arch = "x86_64";
#elif NACS_CPU_AARCH32
#  ifdef __BIG_ENDIAN__
    static const std::string arch = "armeb";
#  else
    static const std::string arch = "arm";
#  endif
#elif NACS_CPU_AARCH64
#  ifdef __BIG_ENDIAN__
    static const std::string arch = "aarch64_be";
#  else
    static const std::string arch = "aarch64";
#  endif
#elif NACS_CPU_PPC32
    static const std::string arch = "powerpc";
#elif NACS_CPU_PPC64
#  ifdef __BIG_ENDIAN__
    static const std::string arch = "powerpc64";
#  else
    static const std::string arch = "powerpc64le";
#  endif
#else
    static const std::string arch = "unknown";
#endif
    return arch;
}

const std::string &get_cpu_name()
{
    static const std::string name("generic");
    return name;
}

const std::string &get_cpu_features()
{
    static const std::string features("");
    return features;
}

}
}
