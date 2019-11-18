#!/usr/bin/python

# Copyright (c) 2019 - 2019 Yichao Yu <yyc1992@gmail.com>
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3.0 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library. If not,
# see <http://www.gnu.org/licenses/>.

import asyncio
import concurrent.futures
from enum import Enum
import json
import threading
import time

class JSONStream:
    def __init__(self):
        self.__buff = ''
        self.__decoder = json.JSONDecoder()

    def __iter__(self):
        return self

    def __next__(self):
        try:
            obj, pos = self.__decoder.raw_decode(self.__buff)
        except json.JSONDecodeError as e:
            # If the error happens past the last character, we have an incomplete object.
            # Stash it and continue. Otherwise:
            if e.pos != len(self.__buff):
                # Shouldn't really happen but if it does, throw away invalid data.
                self.__buff = ''
            raise StopIteration
        self.__buff = self.__buff[pos:].lstrip()
        return obj

    def add_bytes(self, data):
        if self.__buff:
            self.__buff += data.decode()
        else:
            self.__buff = data.decode().lstrip()

class MSquaredProtocol(asyncio.Protocol):
    def __init__(self):
        loop = asyncio.get_running_loop()
        self.__done = loop.create_future()
        self.__transport = None
        self.__id = 1
        self.__json = JSONStream()
        self.__handlers = {}
        self.__background_task = None

    async def __background(self):
        while True:
            await asyncio.sleep(1)
            self.send_raw('ping', text_in="Ping")
            t = time.time()
            for key in list(self.__handlers):
                val = self.__handlers[key]
                if not val:
                    del self.__handlers[key]
                    continue
                i = 0
                while i < len(val):
                    if val[i][1] < t:
                        try:
                            val.pop(i)[0].set_result(None)
                        except:
                            pass
                    i += 1

    def connection_made(self, transport):
        self.__transport = transport
        self.__background_task = asyncio.shield(self.__background())

    def data_received(self, data):
        self.__json.add_bytes(data)
        for obj in self.__json:
            try:
                msg = obj['message']
                param = msg['parameters'] if 'parameters' in msg else {}
                self.message_recieved(msg['op'], param)
            except:
                pass

    def message_recieved(self, op, msg):
        if not op in self.__handlers:
            return
        handlers = self.__handlers[op]
        if not handlers:
            return
        hdl = handlers.pop(0)
        try:
            hdl[0].set_result(msg)
        except:
            pass

    def add_handler(self, op, future, timeout=30):
        if not op in self.__handlers:
            handlers = []
            self.__handlers[op] = handlers
        else:
            handlers = self.__handlers[op]
        handlers.append((future, time.time() + timeout))

    def connection_lost(self, exc):
        self.__transport.close()
        self.__background_task.cancel()
        self.__done.set_result(True)

    def send_raw(self, op, **params):
        if op == 'ping':
            # Don't do useless counting with `ping`.
            tid = 0
        else:
            tid = self.__id
            self.__id += 1
        if params:
            msg = (b'{"message":{"transmission_id":[' + str(tid).encode() + b'],"op":"' +
                       op.encode() + b'","parameters":' +
                       json.dumps(params, separators=(',', ':')).encode() + b'}}')
        else:
            msg = (b'{"message":{"transmission_id":[' + str(tid).encode() + b'],"op":"' +
                       op.encode() + b'"}}')
        self.__transport.write(msg)

    def done(self):
        return self.__done

    def stop(self):
        self.__transport.close()

class MSquared:
    class State(Enum):
        Stopped = 0
        Starting = 1
        Started = 2
        Stopping = 3

    def run_in_loop(func):
        def f(self, *args, **kwargs):
            async def func1():
                await self.__initialized
                return func(self, *args, **kwargs)
            future = asyncio.run_coroutine_threadsafe(func1(), self.__loop)
            try:
                return future.result()
            except concurrent.futures.CancelledError:
                return
        return f

    def with_reply(*names):
        def real_dec(func):
            def f(self, *args, **kwargs):
                def new_future(name):
                    future = concurrent.futures.Future()
                    self.__protocol.add_handler(name, future)
                    future.set_running_or_notify_cancel()
                    return future
                res = tuple(new_future(name) for name in names)
                func(self, *args, **kwargs)
                return res
            return f
        return real_dec

    def __init__(self, addr, host_addr):
        self.__addr = addr
        self.__host_addr = host_addr
        self.__loop = None
        self.__initialized = None
        self.__protocol = None
        self.__thread = None
        self.__state = self.State.Stopped
        self.__init_cv = threading.Condition()
        self.start()

    @run_in_loop
    def send_raw(self, *args, **kwargs):
        self.__protocol.send_raw(*args, **kwargs)

    @run_in_loop
    @with_reply('move_wave_t_reply', 'move_wave_t_f_r')
    def move_wave(self, wl):
        self.__protocol.send_raw('move_wave_t', wavelength=[wl], report="finished")

    @run_in_loop
    @with_reply('poll_move_wave_t_reply')
    def poll_move_wave(self):
        self.__protocol.send_raw('poll_move_wave_t')

    @run_in_loop
    @with_reply('tune_etalon_reply', 'tune_etalon_f_r')
    def tune_etalon(self, percent):
        self.__protocol.send_raw('tune_etalon', setting=[percent], report="finished")

    @run_in_loop
    @with_reply('tune_resonator_reply', 'tune_resonator_f_r')
    def tune_resonator(self, percent):
        self.__protocol.send_raw('tune_resonator', setting=[percent], report="finished")

    @run_in_loop
    @with_reply('fine_tune_resonator_reply', 'fine_tune_resonator_f_r')
    def fine_tune_resonator(self, percent):
        self.__protocol.send_raw('fine_tune_resonator', setting=[percent], report="finished")

    @run_in_loop
    @with_reply('etalon_lock_reply', 'etalon_lock_f_r')
    def etalon_lock(self, on):
        self.__protocol.send_raw('etalon_lock', operation="on" if on else "off",
                                     report="finished")

    @run_in_loop
    @with_reply('etalon_lock_status_reply')
    def etalon_lock_status(self):
        self.__protocol.send_raw('etalon_lock_status')

    @run_in_loop
    @with_reply('select_profile_reply')
    def select_etalon_profile(self, profile):
        self.__protocol.send_raw('select_profile', profile=[profile])

    @run_in_loop
    @with_reply('get_status_reply')
    def system_status(self):
        self.__protocol.send_raw('get_status')

    @run_in_loop
    @with_reply('get_alignment_status_reply')
    def alignment_status(self):
        self.__protocol.send_raw('get_alignment_status')

    @run_in_loop
    @with_reply('beam_alignment_reply', 'beam_alignment_f_r')
    def set_alignment_mode(self, mode):
        self.__protocol.send_raw('beam_alignment', mode=[mode], report="finished")

    @run_in_loop
    @with_reply('beam_adjust_x_reply', 'beam_adjust_x_f_r')
    def alignment_adjust_x(self, val):
        self.__protocol.send_raw('beam_adjust_x', x_value=[val], report="finished")

    @run_in_loop
    @with_reply('beam_adjust_y_reply', 'beam_adjust_y_f_r')
    def alignment_adjust_y(self, val):
        self.__protocol.send_raw('beam_adjust_y', y_value=[val], report="finished")

    @run_in_loop
    @with_reply('get_wavelength_range_reply')
    def wavelength_range(self):
        self.__protocol.send_raw('get_wavelength_range')

    async def __run(self):
        self.__loop = asyncio.get_running_loop()
        self.__initialized = self.__loop.create_future()
        with self.__init_cv:
            if self.__state != self.State.Starting:
                return
            self.__state = self.State.Started
            self.__init_cv.notify_all()

        try:
            transport, self.__protocol = await self.__loop.create_connection(
                MSquaredProtocol, *self.__addr)
        except:
            await asyncio.sleep(1)
            return

        self.__protocol.send_raw('start_link', ip_address=self.__host_addr)
        self.__initialized.set_result(True)

        try:
            await self.__protocol.done()
        finally:
            pass

    def __thread_fun(self):
        while True:
            asyncio.run(self.__run())
            with self.__init_cv:
                if self.__state == self.State.Stopping:
                    return
                self.__state = self.State.Starting

    def start(self):
        # `__thread` is only written to by the controller thread. no need for synchronization
        if self.__thread is not None:
            return
        with self.__init_cv:
            self.__state = self.State.Starting
            self.__thread = threading.Thread(target = self.__thread_fun)
            self.__thread.start()
            self.__init_cv.wait_for(lambda: self.__state == self.State.Started)

    @run_in_loop
    def __req_stop(self):
        self.__protocol.stop()

    def stop(self):
        if self.__thread is None:
            return
        # `__state` can be `Starting` or `Started`
        # The worker thread checks the state before changing it
        # and the `__req_stop` below forces a state change so we should stop it soon enough
        with self.__init_cv:
            self.__state = self.State.Stopping
        self.__req_stop()
        self.__thread.join()
        self.__thread = None
        self.__state = self.State.Stopped

    # This can throw `CancelledError` if the connection is restarted
    @staticmethod
    def wait(futures, timeout=1):
        try:
            return tuple(f.result(timeout) for f in futures)
        except concurrent.futures.TimeoutError:
            return

    del run_in_loop
    del with_reply
