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

#include "execute.h"

namespace NaCs {
namespace LLVM {
namespace Exe {

/**
 */

MemMgr::MemMgr()
{
}

MemMgr::~MemMgr()
{
}

uint8_t *MemMgr::allocateCodeSection(uintptr_t sz, unsigned align,
                                     unsigned sid, StringRef sname)
{
    return tmp.allocateCodeSection(sz, align, sid, sname);
}

uint8_t *MemMgr::allocateDataSection(uintptr_t sz, unsigned align,
                                     unsigned sid, StringRef sname, bool ro)
{
    return tmp.allocateDataSection(sz, align, sid, sname, ro);
}

void MemMgr::registerEHFrames(uint8_t*, uint64_t, size_t)
{
}

void MemMgr::deregisterEHFrames()
{
}

bool MemMgr::finalizeMemory(std::string *ErrMsg)
{
    return tmp.finalizeMemory(ErrMsg);
}

}
}
}
