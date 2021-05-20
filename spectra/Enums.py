
#-------------------------------------------------------------------------------
# Enum definition in spectra
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------



from enum import Enum

class E_ATOMIC_DATA_SOURCE(Enum):
    INTERPOLATE : int = 1
    CALCULATE   : int = 2

class E_ATOM(Enum):
    HYDROGEN          : int = 1
    HYDROGEN_LIKE     : int = 2
    NORMAL            : int = 3

class E_COLLISIONAL_TRANSITION(Enum):
    EXCITATION     : int = 1
    IONIZATION     : int = 2

class E_COLLISIONAL_TRANSITION_SOURCE(Enum):
    EXCITATION     : int = 1
    IONIZATION     : int = 2


class E_COLLISIONAL_TRANSITION_FORMULA(Enum):
    OMEGA          : int = 1
    UNDEFINED      : int = 2

class E_ABSORPTION_PROFILE_TYPE(Enum):
    VOIGT          : int = 0
    GAUSSIAN       : int = 1



#import numba as _numba
#T_DATA_MEMTYPE = _numba.typeof( T_DATA.CALCULATE )
#T_ATOM_MEMTYPE = _numba.typeof( T_ATOM.HYDROGEN )