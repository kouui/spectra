

#-------------------------------------------------------------------------------
# function definition of basic math
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

from ..ImportAll import *

def is_odd_(n : T_INT) -> T_INT:
    """A fast method to check whether a number is odd.

    Parameters
    ----------
    n : T_INT
        input integer number

    Returns
    -------
    T_INT
        1 if input number is an odd number, otherwise 0
    """
    return n & 0x1

#-----------------------------------------------------------------------------
# numba optimization
#-----------------------------------------------------------------------------

if CFG._IS_JIT:

    is_odd_ = nb_vec( **NB_VEC_KWGS )( is_odd_ )



