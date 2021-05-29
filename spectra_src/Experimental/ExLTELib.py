#-------------------------------------------------------------------------------
# experimental function/struct for manipulating line 
#-------------------------------------------------------------------------------
# VERSION
# 0.1.1
#    2021/05/29   k.i., u.k. 
#        - 
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

from ..ImportAll import *

from ..Atomic import LTELib as _LTELib
from ..Math import Integrate as _Integrate

from . import ExTau as _Tau
from ..Util import ElementUtil as _ElementUtil
from ..Util import RomanUtil as _RomanUtil

import numpy as _numpy

def lte_integ_(Z : T_ARRAY, T : T_ARRAY, x : T_ARRAY, wl : T_FLOAT, 
               um : T_FLOAT) -> T_TUPLE[T_FLOAT, T_ARRAY, T_ARRAY]:
    r"""calculate lte output intensity

    Parameters
    ----------
    Z : T_ARRAY
        height
        [:math:`cm`]
    T : T_ARRAY
        temperature
        [:math:`K`]
    x : T_ARRAY
        continuum + line opacity
        [:math:`cm^{-1}`]
    wl : T_FLOAT
        wavelength
        [:math:`cm`]
    um : T_FLOAT
        direction of line of sight = cos(th)
        [-]

    Returns
    -------
    T_TUPLE[T_FLOAT, T_ARRAY, T_ARRAY]

    intensity : T_FLOAT
        intensity
        [:math:`erg \cdot cm^{-2} \cdot Sr^{-1} \cdot cm^{-1} \cdot s^{-1}`]
    
    contrib : T_ARRAY
        contribution function
        [?]
    
    tau : T_ARRAY
        optical depth
        [-]

    Notes
    ------
        
        k.i.  '94/12/05   add contrib keyword
        k.i.  '95/10/26   add tau keyword
        2021.05.16     ki   from idl

    """

    um_m1 = 1. / um

    tau = _Tau.tau_( Z[:], x[:] )
    tau[:] *= um_m1
    contrib : T_ARRAY = _LTELib.planck_cm_( wl, T[:] ) * _numpy.exp( -tau[:] ) * x[:] * um_m1
    intensity = _Integrate.trapze_(contrib[:], Z[:])

    return intensity, contrib, tau

#@OVERLOAD
#def log_saha_(ion : T_STR, T : T_FLOAT) -> T_FLOAT: ...
#@OVERLOAD
#def log_saha_(ion : T_STR, T : T_ARRAY) -> T_ARRAY: ...
#
#def log_saha_(ion : T_STR, T : T_VEC_FA) -> T_VEC_FA:
def log_saha_(ion : T_STR, T : T_FLOAT) -> T_FLOAT :
    ci = 2.07e-16
    ion1  = _ElementUtil.format_ion_( ion )
    ion1p = _ElementUtil.shfit_ion_( ion, 1 )

    ionization_potential : T_FLOAT = _ElementUtil.ion_to_ioniz_potential_( ion )
    lphais = _numpy.log( _LTELib.Ufunc_(ion1,T) / _LTELib.Ufunc_(ion1p,T) ) + \
             11604.52 * ionization_potential / T + \
             _numpy.log((T**(-1.5))*ci)
    
    return lphais


def population_to_H_(ion : T_STR, ep : T_FLOAT, T : T_FLOAT, Ne : T_FLOAT) -> T_FLOAT:
    r"""abandance of an element in some ionize. & exitation state
    relative to hydrogen devided by g (statistical weight)

    Parameters
    ----------
    ion : T_STR
        element & ionization stage as Ca_I, Fe_II, S_I, etc.
    ep : T_FLOAT
        excitation potential of the level
        [:math:`erg`]
    T : T_FLOAT
        temperature
        [:math:`K`]
    Ne : T_FLOAT
        electron number density
        [:math:`cm^{-3}`]

    Returns
    -------
    T_FLOAT
        population ratio
        [-]
    """

    sym, stage = _ElementUtil.ion_to_sym_and_stage_( ion )
    stage_int = _RomanUtil.roman_to_index_( stage ) # ex. I -> 1, II -> 2

    if stage_int > 3:
        raise ValueError("Only ionization stage I, II, III are available")
    si = stage_int - 1

    si_upper_limit = 2
    if sym == "H":
        si_upper_limit = 1

    ## : log (n(j)/n(j+1)/ne)   by lte saha's eq.
    lpjk = 0.
    for l in range(si, si_upper_limit):
        ion1 = _ElementUtil.sym_and_stage_to_ion_( sym, _RomanUtil.index_to_roman_( l+1 ) )
        lpjk += _numpy.log( Ne ) + log_saha_( ion1, T )
    
    ss = 0.
    for m in range(0, si_upper_limit+1):
        lsm = 0.
        for l in range(m, si_upper_limit):
            ion1 = _ElementUtil.sym_and_stage_to_ion_( sym, _RomanUtil.index_to_roman_( l+1 ) )
            lsm += _numpy.log( Ne ) + log_saha_( ion1, T )
            diff = lsm - lpjk
            if diff > 70.:
                diff = 70.
            ss += _numpy.exp( diff )
    
    rijk = _ElementUtil.sym_to_abun_( sym ) * \
           _numpy.exp( -11_604.52 * ep / T ) / \
           _LTELib.Ufunc_(ion, T) / ss

    return rijk



