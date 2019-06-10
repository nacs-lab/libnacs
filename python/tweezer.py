#

from . import handle
import ctypes
import array

nacs_tweezer_max_sideband_ratios = handle.nacs_utils.nacs_tweezer_max_sideband_ratios
nacs_tweezer_max_sideband_ratios.restype = None
nacs_tweezer_max_sideband_ratios.argtypes = [ctypes.c_double, ctypes.c_int,
                                             ctypes.POINTER(ctypes.c_double)]

nacs_tweezer_sideband = handle.nacs_utils.nacs_tweezer_sideband
nacs_tweezer_sideband.restype = ctypes.c_double
nacs_tweezer_sideband.argtypes = [ctypes.c_double, ctypes.c_int, ctypes.c_int]

def max_sideband_ratios(eta, orders):
    ptr = (ctypes.c_double * orders)()
    nacs_tweezer_max_sideband_ratios(eta, orders, ptr)
    return array.array('d', ptr)

def sideband(eta, n1, n2):
    return nacs_tweezer_sideband(eta, n1, n2)
