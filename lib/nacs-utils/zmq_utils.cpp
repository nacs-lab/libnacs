/*************************************************************************
 *   Copyright (c) 2021 - 2021 Yichao Yu <yyc1992@gmail.com>             *
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

#include "zmq_utils.h"

#include <inttypes.h>
#include <stdio.h>

namespace NaCs::ZMQ {

NACS_EXPORT() zmq::context_t &global_context()
{
    static zmq::context_t context;
    return context;
}

NACS_EXPORT() std::pair<zmq::socket_t,zmq::socket_t> inproc_socketpair(zmq::context_t &ctx)
{
    // This should be unique among concurrent callers
    void *ptr = __builtin_frame_address(0);
    constexpr char fmt[] = "inproc://nacs/zmq/pair/%" PRIxPTR;
    char addr[sizeof(fmt) + sizeof(ptr) * 2];
    snprintf(addr, sizeof(addr), fmt, (uintptr_t)ptr);
    std::pair<zmq::socket_t,zmq::socket_t> sockets{zmq::socket_t{ctx, ZMQ_PAIR},
        zmq::socket_t{ctx, ZMQ_PAIR}};
    sockets.first.bind(addr);
    sockets.second.connect(addr);
    sockets.first.unbind(addr);
    return sockets;
}

NACS_EXPORT() MultiClient &MultiClient::global()
{
    static MultiClient multi_client(global_context());
    return multi_client;
}

NACS_EXPORT_ MultiClient::MultiClient(zmq::context_t &context)
    : m_context(context),
      m_cmd_sockets(inproc_socketpair(context))
{
}

NACS_EXPORT() MultiClient::~MultiClient()
{
    if (m_worker.joinable()) {
        send(m_cmd_sockets.second, bits_msg((void*)nullptr));
        m_worker.join();
    }
    m_sockets.clear();
}

NACS_EXPORT() uint64_t MultiClient::SockRef::send_addr()
{
    // Assume the client lock is held.
    uint64_t id = ++m_info.msg_counter;
    auto &sock = m_info.client.m_cmd_sockets.second;
    send_more(sock, bits_msg(&m_info));
    send_more(sock, bits_msg(id));
    send_more(sock, m_info.empty);
    return id;
}

NACS_EXPORT() void MultiClient::SockRef::add_wait(uint64_t id, finish_cb_t finish_cb)
{
    // Assume the client lock is held.
    m_info.callbacks.emplace(id, finish_cb);
}

NACS_EXPORT() MultiClient::SockRef MultiClient::get_socket(const std::string &addr)
{
    std::lock_guard<std::mutex>locker(m_lock);
    auto it = m_sockets.find(addr);
    if (it == m_sockets.end()) {
        it = m_sockets.try_emplace(addr, *this).first;
        auto &info = it->second;
        set_linger(info.sock, 0);
        info.sock.connect(addr);
        ensure_worker();
    }
    return SockRef(it->second);
}

void MultiClient::ensure_worker()
{
    // Assume the client lock is held.
    if (m_running)
        return;
    if (m_worker.joinable())
        m_worker.join();
    m_running = true;
    m_worker = std::thread(&MultiClient::worker_func, this);
}

void MultiClient::worker_func()
{
    // TODO: use poller when it's finalized.
    std::vector<zmq::pollitem_t> poll_items;
    poll_items.push_back({(void*)m_cmd_sockets.first, 0, ZMQ_POLLIN, 0});
    std::vector<decltype(m_sockets.begin())> sockets;
    zmq::message_t msg;
    std::vector<zmq::message_t> msgs;
    while (true) {
        poll_items.resize(1);
        // Cache the socket since other thread may have added new ones when we are polling.
        sockets.clear();
        std::unique_lock<std::mutex> locker(m_lock);
        for (auto it = m_sockets.begin(), end = m_sockets.end(); it != end;) {
            auto &info = it->second;
            if (info.ref_count == 0 && info.callbacks.empty()) {
                it = m_sockets.erase(it);
            }
            else {
                poll_items.push_back({(void*)info.sock, 0, ZMQ_POLLIN, 0});
                sockets.push_back(it);
                ++it;
            }
        }
        locker.unlock();
        if (poll_items.size() == 1)
            break;
        // `zmq::poll(poll_items)` causes deprecation warning.
        // Ref https://github.com/zeromq/cppzmq/issues/494
        zmq::poll(poll_items.data(), poll_items.size());
        if (poll_items[0].revents) {
            assert(poll_items[0].revents == ZMQ_POLLIN);
            recv(m_cmd_sockets.first, msg);
            assert(msg.size() == sizeof(void*));
            SocketInfo *info;
            memcpy(&info, msg.data(), sizeof(void*));
            if (!info)
                break;
            while (true) {
                recv(m_cmd_sockets.first, msg);
                if (has_more(m_cmd_sockets.first)) {
                    send_more(info->sock, msg);
                }
                else {
                    send(info->sock, msg);
                    break;
                }
            }
        }
        auto nsockets = sockets.size();
        for (size_t i = 0; i < nsockets; i++) {
            auto &item = poll_items[i + 1];
            if (!item.revents)
                continue;
            assert(item.revents == ZMQ_POLLIN);
            auto &info = sockets[i]->second;
            recv(info.sock, msg);
            if (msg.size() != 8) {
                readall(info.sock);
                continue;
            }
            uint64_t id;
            memcpy(&id, msg.data(), 8);
            if (!recv_more(info.sock, msg))
                continue;
            // Make sure we get the delimiter frame.
            if (msg.size() != 0) {
                readall(info.sock);
                continue;
            }
            locker.lock();
            auto cb_it = info.callbacks.find(id);
            if (cb_it == info.callbacks.end()) {
                locker.unlock();
                readall(info.sock);
                continue;
            }
            auto cb = cb_it->second;
            info.callbacks.erase(cb_it);
            locker.unlock();
            while (recv_more(info.sock, msg))
                msgs.push_back(std::move(msg));
            try {
                cb(nullptr, std::move(msgs));
            }
            catch (...) {
                // Ignore
            }
            msgs.clear();
        }
    }
    m_running = false;
}

}
