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

#include <ostream>

// This file contains a few classes/functions that are useful for both bytecode
// and command list.

namespace {

struct Printer {
    void ttl(uint32_t ttl, uint64_t t)
    {
        stm << "ttl=0x" << std::hex << ttl << std::dec << " t=" << t << std::endl;
    }
    void ttl1(uint8_t chn, bool val, uint64_t t)
    {
        stm << "ttl(" << int(chn) << ")=" << int(val) << " t=" << t << std::endl;
    }
    void dds_freq(uint8_t chn, uint32_t freq)
    {
        stm << "freq(" << int(chn) << ")=0x"
            << std::hex << freq << std::dec << std::endl;
    }
    void dds_amp(uint8_t chn, uint16_t amp)
    {
        stm << "amp(" << int(chn) << ")=0x"
            << std::hex << amp << std::dec << std::endl;
    }
    void dds_phase(uint8_t chn, uint16_t phase)
    {
        stm << "phase(" << int(chn) << ")=0x"
            << std::hex << phase << std::dec << std::endl;
    }
    void dds_detphase(uint8_t chn, uint16_t detphase)
    {
        stm << "phase(" << int(chn) << ")+=0x"
            << std::hex << detphase << std::dec << std::endl;
    }
    void dds_reset(uint8_t chn)
    {
        stm << "reset(" << int(chn) << ")" << std::endl;
    }
    void dac(uint8_t chn, uint16_t V)
    {
        stm << "dac(" << int(chn) << ")=" << std::hex << V << std::dec << std::endl;
    }
    void wait(uint64_t t)
    {
        stm << "wait(" << t << ")" << std::endl;
    }
    void clock(uint8_t period)
    {
        stm << "clock(" << int(period) << ")" << std::endl;
    }
    std::ostream &stm;
};

struct TimeKeeper {
    void ttl(uint32_t, uint64_t t)
    {
        total_t += t;
    }
    void ttl1(uint8_t, bool, uint64_t t)
    {
        total_t += t;
    }
    void dds_freq(uint8_t, uint32_t)
    {
        total_t += 50;
    }
    void dds_amp(uint8_t, uint16_t)
    {
        total_t += 50;
    }
    void dds_phase(uint8_t, uint16_t)
    {
        total_t += 50;
    }
    void dds_detphase(uint8_t, uint16_t)
    {
        total_t += 50;
    }
    void dds_reset(uint8_t)
    {
        total_t += 50;
    }
    void dac(uint8_t, uint16_t)
    {
        total_t += 45;
    }
    void wait(uint64_t t)
    {
        total_t += t;
    }
    void clock(uint8_t)
    {
        total_t += 5;
    }
    uint64_t total_t = 0;
};

// Simple wrapper around malloc for C interop.
// Implement an API close enough to `std::vector` for our purpose.
// Note that the destructor of this type does **NOT** free the memory.
template<typename ET>
struct MallocVector {
    MallocVector()
        : m_data(nullptr),
          m_len(0)
    {}
    static_assert(std::is_trivial<ET>::value, "");
    ET &operator[](size_t i)
    {
        return m_data[i];
    }
    const ET &operator[](size_t i) const
    {
        return m_data[i];
    }
    size_t size() const
    {
        return m_len;
    }
    void resize(size_t sz)
    {
        m_data = (ET*)realloc(m_data, sizeof(ET) * sz);
        m_len = sz;
    }
private:
    ET *m_data;
    size_t m_len;
};

}
