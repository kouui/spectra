#-------------------------------------------------------------------------------
# function/struct for computing optical depth
#-------------------------------------------------------------------------------
# VERSION
#
# 0.1.0 
#    2021/07/06   u.k.   ...
#-------------------------------------------------------------------------------

from ..ImportAll import *

import numpy as _numpy

def make_tau_(ND  : T_INT, e_max : T_FLOAT, relax : T_FLOAT = 0.):
    """
    make a 1D mesh of log scale from 0 to 1Ee_max with ND grid points.

    Input:
        ND: (,), number of grid points
        e_max: (,), with minvalue of 1Ee_max
        relax: (,), 0 -> 1  logspace -> linspace

    Output:
        Tau: (ND,), optical depth,
    """
    if relax < 0.:
        relax = 0
    if relax > 1.:
        relax = 1.

    e_min = -4.
    tau = _numpy.empty(ND,dtype=DT_NB_FLOAT)
    mid = e_max - _numpy.log10(2.)
    tau[:ND//2+1] = _numpy.logspace(e_min,mid,ND//2+1)*(1-relax) + _numpy.linspace(10.**e_min,10.**mid,ND//2+1)*relax
        
    tau[0] = 0
    tau[ND//2+1:] = (-tau[:ND//2+1][::-1]+2*tau[ND//2])[1:]

    return tau

def z_to_dtau_(z : T_ARRAY, alpha : T_ARRAY):
    r"""[summary]

    Parameters
    ----------
    z : T_ARRAY, [cm]
        1D depth with z[0] as the surface and z[i] < z[i+1]
    alpha : T_ARRAY, [cm^{-1}]
        extinction coefficient

    Returns
    --------
    dtau : T_ARRAY, [-]
        difference of optical depth
    """
    nZ = z.shape[0]
    dtau = _numpy.zeros(nZ, dtype=DT_NB_FLOAT)
    