#

from . import handle
import ctypes
import array

nacs_utils_new_wavemeter = handle.nacs_utils.nacs_utils_new_wavemeter
nacs_utils_new_wavemeter.restype = ctypes.c_void_p
nacs_utils_new_wavemeter.argtypes = [ctypes.c_double, ctypes.c_double]

nacs_utils_free_wavemeter = handle.nacs_utils.nacs_utils_free_wavemeter
nacs_utils_free_wavemeter.restype = None
nacs_utils_free_wavemeter.argtypes = [ctypes.c_void_p]

nacs_utils_wavemeter_parse = handle.nacs_utils.nacs_utils_wavemeter_parse
nacs_utils_wavemeter_parse.restype = ctypes.c_size_t
nacs_utils_wavemeter_parse.argtypes = [ctypes.c_void_p, ctypes.c_char_p,
                                       ctypes.POINTER(ctypes.c_void_p),
                                       ctypes.POINTER(ctypes.c_void_p),
                                       ctypes.c_double, ctypes.c_double]

nacs_utils_wavemeter_clear = handle.nacs_utils.nacs_utils_wavemeter_clear
nacs_utils_wavemeter_clear.restype = None
nacs_utils_wavemeter_clear.argtypes = [ctypes.c_void_p]

PyBytes_FromStringAndSize = ctypes.pythonapi['PyBytes_FromStringAndSize']
PyBytes_FromStringAndSize.restype = ctypes.py_object
PyBytes_FromStringAndSize.argtypes = [ctypes.c_void_p, ctypes.c_ssize_t]

class WavemeterParser:
    def __init__(self, lo=0, hi=1.7976931348623157e308):
        self.hdl = nacs_utils_new_wavemeter(lo, hi)

    def parse(self, name, tstart=0, tend=1.7976931348623157e308):
        ptime = ctypes.c_void_p()
        pdata = ctypes.c_void_p()
        sz = nacs_utils_wavemeter_parse(self.hdl, name.encode(),
                                        ctypes.byref(ptime), ctypes.byref(pdata),
                                        tstart, tend)
        if sz == 0:
            return array.array('d'), array.array('d')
        btime = PyBytes_FromStringAndSize(ptime, sz * 8)
        bdata = PyBytes_FromStringAndSize(pdata, sz * 8)
        return array.array('d', btime), array.array('d', bdata)

    def clear(self):
        nacs_utils_wavemeter_clear(self.hdl)

    def __del__(self):
        nacs_utils_free_wavemeter(self.hdl)
