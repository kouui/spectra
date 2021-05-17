
#-------------------------------------------------------------------------------
# Global imports for the sake of convenience in spectra
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

#from .ImportExternalModule import *

from .Types import *

from .Enums import E_DATA, E_ATOM
from .Elements import ELEMENT_DICT

from . import Constants as _CST
from . import Configurations as _CFG

#-------------------------------------------------------------------------------
# numba related functions/class
#-------------------------------------------------------------------------------

## : numba TypedList slower than built in list
# https://github.com/numba/numba/issues/4584
#from numba.typed import List as nb_List # type: ignore
nb_List = list

## : numba TypedDict much slower than numpy struct array
# comment : currently we will not use dictionary as data struct in spectra
# https://github.com/numba/numba/issues/4364
#from numba.typed import Dict as nb_Dict # type: ignore

from numba import njit as nb_njit       # type: ignore
from numba import vectorize as nb_vec   # type: ignore
from numba.experimental import jitclass as nb_jitclass # type: ignore
