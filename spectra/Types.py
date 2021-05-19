
#-------------------------------------------------------------------------------
# type definition in spectra
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

from typing import Dict as T_DICT, Type
from typing import List as T_LIST
from typing import Any as T_ANY
from typing import Union as T_UNION 
from typing import Tuple as T_TUPLE
from typing import Callable as T_CALLABLE
from typing import NoReturn as T_NORETURN
from typing import TypedDict as T_TYPEDICT


from typing import Type as T_TYPE



#-------------------------------------------------------------------------------
# fundamental types
#-------------------------------------------------------------------------------

T_FLOAT       = float
T_INT         = int      # _numpy.int64
T_STR         = str
T_BOOL        = bool

T_SLICE       = slice

# reference : https://numpy.org/devdocs/reference/typing.html
from numpy import ndarray as _ndarray
T_ARRAY       = _ndarray
del _ndarray


#-------------------------------------------------------------------------------
# advanced types
#-------------------------------------------------------------------------------

T_VEC_IFA   = T_UNION[ T_INT, T_FLOAT, T_ARRAY ]
T_VEC_FA    = T_UNION[ T_FLOAT, T_ARRAY ]
T_VEC_IA    = T_UNION[ T_INT, T_ARRAY ]

T_UNION_IS = T_UNION[ T_INT, T_SLICE ]

#-------------------------------------------------------------------------------
# numba types
#-------------------------------------------------------------------------------

#import numba.types as nb_types     # type: ignore
#
#T_NB_FLOAT1D = nb_types.float64[:]         # type: ignore
#T_NB_FLOAT2D = nb_types.float64[:,:]       # type: ignore
#T_NB_FLOAT3D = nb_types.float64[:,:,:]     # type: ignore
#
#T_NB_INT1D = nb_types.int64[:]             # type: ignore
#T_NB_INT2D = nb_types.int64[:,:]           # type: ignore
#T_NB_INT3D = nb_types.int64[:,:,:]         # type: ignore
#
#T_NB_STR   = nb_types.unicode_type         # type: ignore


## comment : using numba.typed.Dict will dramatically slow down 
##           import time
## comment : incase of List:
##             - preprocessing : list --> tuple, these are homogeneous tuples
##             - njit          : list (--> tuple), list initialized in njit function will
##                                be converted to numba.typed.List automatically
##                                might be heterogeneous tuples, so do not convert the list to tuple

## : example of using numba types without conflicting with mypy
#
# from numba.typed import List as nb_List # type: ignore
# from numba.typed import Dict as nb_Dict # type: ignore
#nb_List = list

# a : T_LIST[T_ARRAY] = nb_List()
# a.append( numpy.arange(6).reshape(2,3) )
# 
# b : T_DICT[T_STR,T_ARRAY] = nb_Dict.empty(
#     key_type=T_NB_STR,
#     value_type=T_NB_FLOAT1D,
# )
# 
# b['posx'] = numpy.asarray([1, 0.5, 2], dtype=T_FLOAT)
# b['posy'] = numpy.asarray([1.5, 3.5, 2], dtype=T_FLOAT)
# b['velx'] = numpy.asarray([0.5, 0, 0.7], dtype=T_FLOAT)
# b['vely'] = numpy.asarray([0.2, -0.2, 0.1], dtype=T_FLOAT)
