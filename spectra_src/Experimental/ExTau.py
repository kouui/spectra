#-------------------------------------------------------------------------------
# experimental function manipulating optical depth tau
#-------------------------------------------------------------------------------
# VERSION
# 0.1.1
#    2021/05/29   k.i., u.k. 
#        - 
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

from ..ImportAll import *

import numpy as _numpy

## : ichimoto's solution for integrating xc, forward difference
def tau_v1_(Z : T_ARRAY, xi : T_ARRAY) -> T_ARRAY:
    r"""given height `Z` and opacity `xi`, calculate optical depth tau

    Parameters
    ----------
    Z : T_ARRAY
        height with direction : [interior, surface]
        [:math:`cm`]
    xi : T_ARRAY
        opacity
        [:math:`cm^{-1}`]

    Returns
    -------
    T_ARRAY
        optical depth
        [-]
    """

    nZ = Z.shape[0]
    tau = _numpy.zeros( nZ, dtype=DT_NB_FLOAT )
    
    tau[-1] = xi[-1] * ( Z[-1] - Z[-2] ) * 0.5
    for i in range(nZ-2,-1,-1):
        tau[i] = tau[i+1] + ( Z[i+1] - Z[i] ) * ( xi[i] + xi[i+1] ) * 0.5

    return tau

## : kouui's solution for integrating xc, central difference
def tau_v2_(Z : T_ARRAY, xi : T_ARRAY) -> T_ARRAY:
    r"""given height `Z` and opacity `xi`, calculate optical depth tau

    Parameters
    ----------
    Z : T_ARRAY
        height with direction : [interior, surface]
        [:math:`cm`]
    xi : T_ARRAY
        opacity
        [:math:`cm^{-1}`]

    Returns
    -------
    T_ARRAY
        optical depth
        [-]
    """

    nZ = Z.shape[0]
    tau = _numpy.zeros( nZ, dtype=DT_NB_FLOAT )


    tau[-1] = xi[-1] * ( Z[-1] - Z[-2] ) * 0.5
    for i in range(nZ-2,0,-1):
        tau[i] = tau[i+1] + xi[i] * ( Z[i+1] - Z[i-1] ) * 0.5
    tau[0] = tau[1] + xi[0] * ( Z[1] - Z[0] ) * 0.5
    
    return tau


tau_ = tau_v1_
