/*************************************************************************
 *   Copyright (c) 2019 - 2019 Yichao Yu <yyc1992@gmail.com>             *
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

#include "macros.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if !NACS_OS_WINDOWS
#  include <dlfcn.h>
#  include <link.h>

static const ElfW(Dyn) *find_section(const ElfW(Dyn) *dyn, ElfW(Sxword) tag)
{
    for (; dyn->d_tag != DT_NULL; ++dyn) {
        if (dyn->d_tag == tag) {
            return dyn;
        }
    }
    return nullptr;
}

static char *get_soname(const char *name)
{
    auto handle = dlopen(name, RTLD_LAZY | RTLD_LOCAL);
    if (!handle)
        return nullptr;
    // Get the link map
    const link_map *link_map;
    if (dlinfo(handle, RTLD_DI_LINKMAP, &link_map) != 0) {
        dlclose(handle);
        return nullptr;
    }
    auto strtab_sec = find_section(link_map->l_ld, DT_STRTAB);
    auto soname_sec = find_section(link_map->l_ld, DT_SONAME);
    if (!strtab_sec || !soname_sec) {
        dlclose(handle);
        return nullptr;
    }
    auto strtab = (const char*)strtab_sec->d_un.d_ptr;
    auto res = strdup(&strtab[soname_sec->d_un.d_val]);
    dlclose(handle);
    return res;
}
#endif

int main(int argc, const char **argv)
{
    if (argc < 2) {
        fputs("Too few arguments.", stderr);
        return 1;
    }
#if !NACS_OS_WINDOWS
    if (auto soname = get_soname(argv[1])) {
        puts(soname);
        free(soname);
        return 0;
    }
#endif
    puts(argv[1]);
    return 0;
}
