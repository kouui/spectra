
#-------------------------------------------------------------------------------
# Enum definition in spectra
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------



from enum import IntEnum as _IntEnum

class E_ATOMIC_DATA_SOURCE(_IntEnum):
    EXPERIMENT  : int = 1
    CALCULATE   : int = 2

class E_ATOM(_IntEnum):
    HYDROGEN          : int = 1
    HYDROGEN_LIKE     : int = 2
    NORMAL            : int = 3

class E_COLLISIONAL_TRANSITION(_IntEnum):
    EXCITATION     : int = 1
    IONIZATION     : int = 2

class E_COLLISIONAL_TRANSITION_SOURCE(_IntEnum):
    ELECTRON     : int = 1
    PROTON       : int = 2
    CHARGE_TRANSFER : int = 3


class E_COLLISIONAL_TRANSITION_FORMULA(_IntEnum):
    OMEGA          : int = 1

class E_ABSORPTION_PROFILE_TYPE(_IntEnum):
    VOIGT          : int = 0
    GAUSSIAN       : int = 1

class E_ATMOSPHERE_COORDINATE_TYPE(_IntEnum):
    POINT          : int = 0
    CARTESIAN      : int = 1



#import numba as _numba
#T_DATA_MEMTYPE = _numba.typeof( T_DATA.CALCULATE )
#T_ATOM_MEMTYPE = _numba.typeof( T_ATOM.HYDROGEN )