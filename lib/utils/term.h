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

#ifndef __NACS_UTILS_TERM_H__
#define __NACS_UTILS_TERM_H__

#include <ostream>

namespace NaCs {
namespace Term {

struct _TermState {
    static constexpr _TermState get_fg(int fg)
    {
        return _TermState(fg, 0, false, true, false);
    }
    static constexpr _TermState get_bg(int bg)
    {
        return _TermState(0, bg, false, false, true);
    }
    static constexpr _TermState get_bold()
    {
        return _TermState(0, 0, true, false, false);
    }
    constexpr _TermState operator()(bool bold) const
    {
        return _TermState(fg, bg, bold, has_fg, has_bg);
    }
    constexpr _TermState operator|(_TermState other) const
    {
        _TermState res = *this;
        if (other.has_fg) {
            res.has_fg = true;
            res.fg = other.fg;
        }
        if (other.has_bg) {
            res.has_bg = true;
            res.bg = other.bg;
        }
        if (other.bold)
            res.bold = true;
        return res;
    }
    constexpr _TermState(const _TermState&) = default;
    constexpr _TermState &operator=(const _TermState&) = default;
    constexpr _TermState(int fg, int bg, bool bold, bool has_fg, bool has_bg)
    : fg(fg), bg(bg), bold(bold), has_fg(has_fg), has_bg(has_bg)
    {
    }
    int fg;
    int bg;
    bool bold;
    bool has_fg;
    bool has_bg;
};

constexpr auto black = _TermState::get_fg(0);
constexpr auto red = _TermState::get_fg(1);
constexpr auto green = _TermState::get_fg(2);
constexpr auto yellow = _TermState::get_fg(3);
constexpr auto blue = _TermState::get_fg(4);
constexpr auto magenta = _TermState::get_fg(5);
constexpr auto cyan = _TermState::get_fg(6);
constexpr auto white = _TermState::get_fg(7);

constexpr auto black_bg = _TermState::get_bg(0);
constexpr auto red_bg = _TermState::get_bg(1);
constexpr auto green_bg = _TermState::get_bg(2);
constexpr auto yellow_bg = _TermState::get_bg(3);
constexpr auto blue_bg = _TermState::get_bg(4);
constexpr auto magenta_bg = _TermState::get_bg(5);
constexpr auto cyan_bg = _TermState::get_bg(6);
constexpr auto white_bg = _TermState::get_bg(7);

constexpr auto bold = _TermState::get_bold();

NACS_EXPORT(utils) std::ostream &operator<<(std::ostream &stm, _TermState state);
NACS_EXPORT(utils) std::ostream &reset(std::ostream &stm);

}
}

#endif
