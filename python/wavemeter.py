#

from . import handle
import ctypes
import array

nacs_seq_new_wavemeter_parser = handle.nacs_seq.nacs_seq_new_wavemeter_parser
nacs_seq_new_wavemeter_parser.restype = ctypes.c_void_p
nacs_seq_new_wavemeter_parser.argtypes = []

nacs_seq_free_wavemeter_parser = handle.nacs_seq.nacs_seq_free_wavemeter_parser
nacs_seq_free_wavemeter_parser.restype = None
nacs_seq_free_wavemeter_parser.argtypes = [ctypes.c_void_p]

nacs_seq_wavemeter_parse = handle.nacs_seq.nacs_seq_wavemeter_parse
nacs_seq_wavemeter_parse.restype = ctypes.c_size_t
nacs_seq_wavemeter_parse.argtypes = [ctypes.c_void_p, ctypes.c_char_p,
                                     ctypes.POINTER(ctypes.c_void_p),
                                     ctypes.POINTER(ctypes.c_void_p), ctypes.c_int]

PyBytes_FromStringAndSize = ctypes.pythonapi['PyBytes_FromStringAndSize']
PyBytes_FromStringAndSize.restype = ctypes.py_object
PyBytes_FromStringAndSize.argtypes = [ctypes.c_void_p, ctypes.c_ssize_t]

class WavemeterParser:
    def __init__(self):
        self.hdl = nacs_seq_new_wavemeter_parser()

    def parse(self, name, use_cache):
        ptime = ctypes.c_void_p()
        pdata = ctypes.c_void_p()
        sz = nacs_seq_wavemeter_parse(self.hdl, name.encode(),
                                      ctypes.byref(ptime), ctypes.byref(pdata),
                                      use_cache)
        if sz == 0:
            return array.array('d'), array.array('d')
        btime = PyBytes_FromStringAndSize(ptime, sz * 8)
        bdata = PyBytes_FromStringAndSize(pdata, sz * 8)
        return array.array('d', btime), array.array('d', bdata)

    def __del__(self):
        nacs_seq_free_wavemeter_parser(self.hdl)
