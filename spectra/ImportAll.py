
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

from . import Constants as CST
from . import Configurations as CFG

#-------------------------------------------------------------------------------
# numba related functions/class
#-------------------------------------------------------------------------------

## : numba TypedList slower than built in list
# https://github.com/numba/numba/issues/4584
#from numba.typed import List as nb_List # type: ignore
from numba.typed import List # type: ignore
nb_List : T_TYPE[ T_UNION[List,T_LIST] ]
if CFG._IS_JIT:
    #from numba.typed import List # type: ignore
    nb_List = List
else:
    nb_List = list
del List

## : numba TypedDict much slower than numpy struct array
# comment : currently we will not use dictionary as a data struct in spectra
# https://github.com/numba/numba/issues/4364
#from numba.typed import Dict as nb_Dict # type: ignore

from numba import njit as nb_njit # type: ignore
from numba import vectorize as nb_vec # type: ignore

#from numba.experimental import jitclass as nb_jitclass # type: ignore
# comment : currently we will not use jitclass as a data struct in spectra

NB_VEC_KWGS : T_DICT[T_STR, T_UNION[T_BOOL, T_STR]] = {
    "cache"  : CFG._IS_CACHE, 
    "target" : CFG._VEC_TARGET, 
    "nogil"  : CFG._IS_NOGIL, 
    "nopython" : True,
    "fastmath" : CFG._IS_FASTMATH,
}

NB_NJIT_KWGS : T_DICT[T_STR, T_BOOL] = {
    "cache" : CFG._IS_CACHE, 
    "nogil" : CFG._IS_NOGIL,
    "fastmath" : CFG._IS_FASTMATH,
    "parallel" : False,
}

NB_NJIT_KWGS_PARALLEL : T_DICT[T_STR, T_BOOL] = {
    "cache" : CFG._IS_CACHE, 
    "nogil" : CFG._IS_NOGIL,
    "fastmath" : CFG._IS_FASTMATH,
    "parallel" : CFG._IS_PARALLEL,
}


#-------------------------------------------------------------------------------
# numpy related functions/class
#-------------------------------------------------------------------------------

from numpy import vectorize as np_vec

NP_VEC_KWGS : T_DICT[T_STR, T_ANY] = {
    "otypes" : [T_FLOAT],
    "cache" : True,
}

#-------------------------------------------------------------------------------
# logging and warning
#-------------------------------------------------------------------------------

import warnings
warnings.simplefilter('always')