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

#define CATCH_CONFIG_MAIN

#include "../../lib/nacs-utils/zmq_utils.h"

#include <atomic>
#include <chrono>
#include <thread>

#include <catch2/catch.hpp>

using namespace NaCs;

struct Server {
    Server()
        : sock(ctx, ZMQ_ROUTER)
    {
    }
    int start()
    {
        for (int i = 4000; i < 8000; i++) {
            try {
                sock.bind("tcp://127.0.0.1:" + std::to_string(i));
                thread = std::thread([&] {
                    zmq::message_t empty;
                    std::vector<std::vector<zmq::message_t>> waiters;
                    auto send_reply = [&] (auto &&addr, auto &&msg) {
                        ZMQ::send_addr(sock, addr, empty);
                        ZMQ::send(sock, msg);
                    };
                    while (true) {
                        auto addr = ZMQ::recv_addr(sock);
                        zmq::message_t msg;
                        if (!ZMQ::recv_more(sock, msg))
                            goto err;
                        if (ZMQ::match(msg, "ping")) {
                            send_reply(addr, ZMQ::str_msg("pong"));
                            goto out;
                        }
                        else if (ZMQ::match(msg, "exit")) {
                            send_reply(addr, ZMQ::str_msg("exit"));
                            ZMQ::readall(sock);
                            break;
                        }
                        else if (ZMQ::match(msg, "wait")) {
                            waiters.push_back(std::move(addr));
                            goto out;
                        }
                        else if (ZMQ::match(msg, "release")) {
                            send_reply(addr, ZMQ::str_msg("releasing"));
                            for (auto &addr: waiters)
                                send_reply(addr, ZMQ::str_msg("finished"));
                            waiters.clear();
                            goto out;
                        }
                    err:
                        send_reply(addr, ZMQ::bits_msg<uint8_t>(1));
                    out:
                        ZMQ::readall(sock);
                    }
                    for (auto &addr: waiters)
                        send_reply(addr, ZMQ::str_msg("cancalled"));
                    waiters.clear();
                });
                printf("Server started on port %d\n", i);
                return i;
            }
            catch (...) {
                continue;
            }
        }
        fprintf(stderr, "Cannot find free port to bind\n");
        abort();
    }

    std::thread thread;
    zmq::context_t ctx;
    zmq::socket_t sock;
};

static void test_str_msg(zmq::message_t &msg, const char *str)
{
    REQUIRE(msg.size() == strlen(str));
    REQUIRE(memcmp(msg.data(), str, strlen(str)) == 0);
}

static void ping_server(ZMQ::MultiClient::SockRef &sock)
{
    auto reply = sock.send_msg([&] (auto &sock) {
        ZMQ::send(sock, ZMQ::str_msg("ping"));
    });
    auto msgs = reply.get();
    REQUIRE(msgs.size() == 1);
    auto &msg = msgs[0];
    test_str_msg(msg, "pong");
}

static void test_multi_ping(ZMQ::MultiClient &client, const char *addr, int nthreads)
{
    std::vector<std::thread> threads(nthreads);
    std::atomic<bool> start{false};
    std::atomic<int> finished{0};
    for (int i = 0; i < nthreads; i++) {
        threads[i] = std::thread([&] {
            while (!start.load(std::memory_order_relaxed))
                CPU::pause();
            auto sock = client.get_socket(addr);
            for (int i = 0; i < 1000; i++)
                ping_server(sock);
            finished.fetch_add(1, std::memory_order_relaxed);
            while (finished.load(std::memory_order_relaxed) < nthreads) {
                ping_server(sock);
            }
        });
    }
    // This of course doesn't guarantee anything
    // but it should increase the contention on the client
    // by letting the threads start closer to each other
    start.store(true, std::memory_order_relaxed);
    CPU::wake();
    for (auto &t: threads) {
        t.join();
    }
}

static void test_interleave(ZMQ::MultiClient::SockRef &sock)
{
    auto reply = sock.send_msg([&] (auto &sock) {
        ZMQ::send(sock, ZMQ::str_msg("wait"));
    });
    for (int i = 0; i < 10; i++)
        ping_server(sock);
    auto status = reply.wait_for(std::chrono::milliseconds(1));
    REQUIRE(status == std::future_status::timeout);
    auto rel_reply = sock.send_msg([&] (auto &sock) {
        ZMQ::send(sock, ZMQ::str_msg("release"));
    });
    auto rel_msgs = rel_reply.get();
    REQUIRE(rel_msgs.size() == 1);
    auto &rel_msg = rel_msgs[0];
    test_str_msg(rel_msg, "releasing");

    auto msgs = reply.get();
    REQUIRE(msgs.size() == 1);
    auto &msg = msgs[0];
    test_str_msg(msg, "finished");
}

static void test_interleave2(ZMQ::MultiClient::SockRef &sock, ZMQ::MultiClient::SockRef &sock2)
{
    auto reply = sock.send_msg([&] (auto &sock) {
        ZMQ::send(sock, ZMQ::str_msg("wait"));
    });
    auto reply2 = sock2.send_msg([&] (auto &sock) {
        ZMQ::send(sock, ZMQ::str_msg("wait"));
    });
    for (int i = 0; i < 10; i++)
        ping_server(sock);
    for (int i = 0; i < 10; i++)
        ping_server(sock2);
    auto status = reply.wait_for(std::chrono::milliseconds(1));
    REQUIRE(status == std::future_status::timeout);
    auto status2 = reply2.wait_for(std::chrono::milliseconds(1));
    REQUIRE(status2 == std::future_status::timeout);

    {
        auto rel_reply = sock.send_msg([&] (auto &sock) {
            ZMQ::send(sock, ZMQ::str_msg("release"));
        });
        auto rel_msgs = rel_reply.get();
        REQUIRE(rel_msgs.size() == 1);
        auto &rel_msg = rel_msgs[0];
        test_str_msg(rel_msg, "releasing");

        auto msgs = reply.get();
        REQUIRE(msgs.size() == 1);
        auto &msg = msgs[0];
        test_str_msg(msg, "finished");
    }

    status2 = reply2.wait_for(std::chrono::milliseconds(1));
    REQUIRE(status2 == std::future_status::timeout);

    {
        auto rel_reply2 = sock2.send_msg([&] (auto &sock) {
            ZMQ::send(sock, ZMQ::str_msg("release"));
        });
        auto rel_msgs2 = rel_reply2.get();
        REQUIRE(rel_msgs2.size() == 1);
        auto &rel_msg2 = rel_msgs2[0];
        test_str_msg(rel_msg2, "releasing");

        auto msgs2 = reply2.get();
        REQUIRE(msgs2.size() == 1);
        auto &msg2 = msgs2[0];
        test_str_msg(msg2, "finished");
    }
}

static void exit_server(ZMQ::MultiClient::SockRef &sock)
{
    auto reply = sock.send_msg([&] (auto &sock) {
        ZMQ::send(sock, ZMQ::str_msg("exit"));
    });
    auto msgs = reply.get();
    REQUIRE(msgs.size() == 1);
    auto &msg = msgs[0];
    test_str_msg(msg, "exit");
}

ANON_TEST_CASE() {
    Server server;
    int port = server.start();
    std::string addr("tcp://127.0.0.1:" + std::to_string(port));
    Server server2;
    int port2 = server2.start();
    std::string addr2("tcp://127.0.0.1:" + std::to_string(port2));

    auto &client = ZMQ::MultiClient::global();
    auto sock = client.get_socket(addr.c_str());
    auto sock2 = client.get_socket(addr2.c_str());

    // Single thread ping.
    ping_server(sock);
    // Multi thread ping.
    for (int i = 0; i < 100; i++) {
        test_multi_ping(client, addr.c_str(), 2);
        test_multi_ping(client, addr.c_str(), 4);
        test_multi_ping(client, addr.c_str(), 8);
    }
    // Make sure we can recieve things out-of-order
    test_interleave(sock);
    test_interleave2(sock, sock2);

    exit_server(sock);
    exit_server(sock2);
    server.thread.join();
    server2.thread.join();
}
