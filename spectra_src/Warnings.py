#-------------------------------------------------------------------------------
# warning function for normal/numba mode
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

from .Types import *
from . import Configurations as CFG

import warnings as WARNINGS
WARNINGS.simplefilter('once', UserWarning)

from numba import objmode # type: ignore

from numba import njit as _nb_njit # type: ignore
_NB_NJIT_KWGS : T_DICT[T_STR, T_BOOL] = {
    "cache" : CFG._IS_CACHE, 
    "nogil" : CFG._IS_NOGIL,
    "fastmath" : CFG._IS_FASTMATH,
    "parallel" : False,
}


def WARN_(text : T_STR):

    if CFG._IS_JIT:
        with objmode():
            WARNINGS.warn( text )
    else:
        WARNINGS.warn( text )


if CFG._IS_JIT:
    WARN_ = _nb_njit(**_NB_NJIT_KWGS) (WARN_)
