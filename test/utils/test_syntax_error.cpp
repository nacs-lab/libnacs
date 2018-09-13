/*************************************************************************
 *   Copyright (c) 2018 - 2018 Yichao Yu <yyc1992@gmail.com>             *
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

#include "../../lib/utils/errors.h"

#include <assert.h>

#include <iostream>
#include <sstream>

using namespace NaCs;

static std::pair<std::string,std::string> get_strings(const SyntaxError &err)
{
    std::stringstream stm;
    stm << err;
    return {err.what(), stm.str()};
}

static void check_err_strings(const SyntaxError &err, const std::string &what,
                              const std::string &printmsg)
{
    auto res = get_strings(err);
    assert(res.first == what);
    assert(res.second == printmsg);
}

int main()
{
    check_err_strings(SyntaxError("Missing `)`", "a = f(1, 2) + g(12];", 10, 19, 16, 20),
                      "SyntaxError: Missing `)`\nL10:19 a = f(1, 2) + g(12];",
                      "SyntaxError: Missing `)`\n"
                      "Line: 10, Column: 19\n"
                      "a = f(1, 2) + g(12];\n"
                      "               ~~~^~\n");
    check_err_strings(SyntaxError("Literal too long", "a = 12345678901234567890;", 3, -1, 5, 24),
                      "SyntaxError: Literal too long\nL3 a = 12345678901234567890;",
                      "SyntaxError: Literal too long\n"
                      "Line: 3\n"
                      "a = 12345678901234567890;\n"
                      "    ~~~~~~~~~~~~~~~~~~~~\n");
    return 0;
}
