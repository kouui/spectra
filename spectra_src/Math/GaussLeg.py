

#-------------------------------------------------------------------------------
# Compute abscissas and weights for N-point Gauss Legendre quadrature formula
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

from ..ImportAll import *

_ERROR = 3.0E-14

import numpy as _numpy

def gauss_quad_coe_( a : T_FLOAT, b : T_FLOAT, n : T_INT ) -> T_TUPLE[T_ARRAY,T_ARRAY]:
    """Compute abscissas and weights for N-point Gauss Legendre quadrature formula

    Parameters
    ----------
    a : T_FLOAT
        lower boundary of x
    b : T_FLOAT
        upper boundary of x
    n : T_INT
        number of point

    Returns
    -------
    T_TUPLE[T_ARRAY,T_ARRAY]
        
    x_arr : T_ARRAY, 1d
        abscissas array
    w_arr : T_ARRAY, 1d
        weight array
    """
    m  = (n + 1) // 2
    xm = 0.5 * (a + b)
    xl = 0.5 * (b - a)

    x_arr = _numpy.empty( n, dtype=DT_NB_FLOAT )
    w_arr = _numpy.empty( n, dtype=DT_NB_FLOAT )

    for i in range( m ):

        zz = _numpy.cos( CST.pi_ * (i + 0.75) / (n + 0.5) )
        
        while True:
            p1 = 1.0
            p2 = 0.0
            for k in range(1, n+1):
                p3 = p2
                p2 = p1
                p1 = (2.0*(k - 0.5)*zz*p2 - (k - 1.0)*p3) / k
            
            pp = n * (zz*p1 - p2) / (zz*zz - 1.0)
            z1 = zz
            zz = z1 - p1/pp
            dz = abs(p1/pp)
            
            if dz < _ERROR:
                break

        x_arr[i] = (xm - xl*zz)
        w_arr[i] = (2.0 * xl / ((1.0 - zz*zz)*pp*pp))
        x_arr[n-1 - i] = (xm + xl*zz)
        w_arr[n-1 - i] = w_arr[i]

    return x_arr, w_arr

#-----------------------------------------------------------------------------
# numba optimization
#-----------------------------------------------------------------------------

if CFG._IS_JIT:

    gauss_quad_coe_ = nb_njit( **NB_NJIT_KWGS )( gauss_quad_coe_ )