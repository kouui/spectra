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

from numba import objmode


def WARN_(text : T_STR):

    if CFG._IS_JIT:
        with objmode():
            WARNINGS.warn( text )
    else:
        WARNINGS.warn( text )
