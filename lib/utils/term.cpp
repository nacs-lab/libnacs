/*************************************************************************
 *   Copyright (c) 2017 - 2018 Yichao Yu <yyc1992@gmail.com>             *
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

#include "utils.h"
#include "term.h"

#if !NACS_OS_WINDOWS
#  include <spawn.h>
#  include <sys/stat.h>
#  include <fcntl.h>
#  include <unistd.h>
#else
#  include <io.h>
#endif

#include <stdio.h>

#include <iostream>

namespace NaCs {
namespace Term {

static bool has_color()
{
    static bool _has_color = [] {
        bool res = false;
#if !NACS_OS_WINDOWS
        pid_t pid;
        char *argv[] = {(char*)"tput", (char*)"colors", NULL};
        posix_spawn_file_actions_t faction;
        if (posix_spawn_file_actions_init(&faction))
            goto ret;
        if (posix_spawn_file_actions_addopen(&faction, 0, "/dev/null",
                                             O_RDONLY, 0))
            goto close_faction;
        if (posix_spawn_file_actions_addopen(&faction, 2, "/dev/null",
                                             O_WRONLY, 0))
            goto close_faction;
        int fds[2];
        if (pipe(fds))
            goto close_faction;
        if (posix_spawn_file_actions_adddup2(&faction, fds[1], 1))
            goto close_pipe;
        if (posix_spawnp(&pid, "tput", &faction, nullptr, argv, environ))
            goto close_pipe;
        close(fds[1]);
        fds[1] = -1;
        {
            char buff[32];
            auto r = read(fds[0], buff, sizeof(buff));
            if (r <= 0)
                goto close_pipe;
            buff[r] = 0;
            res = atoi(buff) >= 8;
        }
    close_pipe:
        close(fds[0]);
        close(fds[1]);
    close_faction:
        posix_spawn_file_actions_destroy(&faction);
    ret:
#endif
        if (!res) {
            if (auto term = getenv("TERM")) {
                res = (strcmp(term, "xterm") == 0 ||
                       strcmp(term, "xterm-256color") == 0);
            }
        }
        return res;
    }();
    return _has_color;
}

static bool isattystm(std::ostream &stm)
{
    // This is the most reliable way I can find to get the corresponding fd from a std::ostream
    auto buf = stm.rdbuf();
#if !NACS_OS_WINDOWS
    if (buf == std::cout.rdbuf())
        return isatty(1);
    if (buf == std::cerr.rdbuf())
        return isatty(2);
#else
    if (buf == std::cout.rdbuf())
        return _isatty(_fileno(stdout));
    if (buf == std::cerr.rdbuf())
        return _isatty(_fileno(stderr));
#endif
    return false;
}

NACS_EXPORT() std::ostream &operator<<(std::ostream &stm, _TermState state)
{
    if (!has_color() || !isattystm(stm))
        return stm;

    char buff[128];
    auto p = buff;
    p[0] = '\033';
    p[1] = '[';
    p += 2;
    if (state.bold) {
        p[0] = '1';
        p[1] = ';';
        p += 2;
    }
    if (state.has_fg) {
        p[0] = '3';
        p[1] = char('0' + state.fg);
        p[2] = ';';
        p += 3;
    }
    if (state.has_bg) {
        p[0] = '4';
        p[1] = char('0' + state.bg);
        p[2] = ';';
        p += 3;
    }
    if (p[-1] == ';')
        p--;
    p[0] = 'm';
    p[1] = 0;

    stm << buff;

    return stm;
}

NACS_EXPORT() std::ostream &reset(std::ostream &stm)
{
    if (!has_color() || !isattystm(stm))
        return stm;
    stm << "\033[0m";
    return stm;
}

}
}
