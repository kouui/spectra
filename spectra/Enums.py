
#-------------------------------------------------------------------------------
# Enum definition in spectra
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------



from enum import Enum

class E_DATA(Enum):
    INTERPOLATE = 1
    CALCULATE   = 2

class E_ATOM(Enum):
    HYDROGEN      = 1
    HYDROGEN_LIKE = 2
    NORMAL        = 3






#import numba as _numba
#T_DATA_MEMTYPE = _numba.typeof( T_DATA.CALCULATE )
#T_ATOM_MEMTYPE = _numba.typeof( T_ATOM.HYDROGEN )