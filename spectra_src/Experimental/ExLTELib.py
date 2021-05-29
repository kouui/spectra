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

from ..Struct import Atmosphere as _Atmosphere
from . import ExLine as _Line
from ..Atomic import ContinuumOpacity as _ContinuumOpacity
from ..Atomic import BasicP as _BasicP
from ..RadiativeTransfer import Profile as _Profile

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

@OVERLOAD
def log_saha_(ion : T_STR, T : T_FLOAT) -> T_FLOAT: ...
@OVERLOAD
def log_saha_(ion : T_STR, T : T_ARRAY) -> T_ARRAY: ...

def log_saha_(ion : T_STR, T : T_VEC_FA) -> T_VEC_FA:
#def log_saha_(ion : T_STR, T : T_FLOAT) -> T_FLOAT :
    ci = 2.07e-16
    ion1  = _ElementUtil.format_ion_( ion )
    ion1p = _ElementUtil.shfit_ion_( ion, 1 )

    ionization_potential : T_FLOAT = _ElementUtil.ion_to_ioniz_potential_( ion )
    lphais = _numpy.log( _LTELib.Ufunc_(ion1,T) / _LTELib.Ufunc_(ion1p,T) ) + \
             CST.eV2erg_ * ionization_potential / (CST.k_ * T) + \
             _numpy.log((T**(-1.5))*ci)
    
    return lphais

#@OVERLOAD
#def population_to_H_(ion : T_STR, ep : T_FLOAT, T : T_FLOAT, Ne : T_FLOAT) -> T_FLOAT: ...
#@OVERLOAD
#def population_to_H_(ion : T_STR, ep : T_FLOAT, T : T_ARRAY, Ne : T_ARRAY) -> T_ARRAY: ...

def population_to_H_(ion : T_STR, ep : T_FLOAT, T : T_VEC_FA, Ne : T_VEC_FA) -> T_VEC_FA:
    r"""abandance of an element in some ionize. & exitation state
    relative to hydrogen devided by g (statistical weight)

    Parameters
    ----------
    ion : T_STR
        element & ionization stage as Ca_I, Fe_II, S_I, etc.
    ep : T_FLOAT
        excitation potential of the level
        [:math:`erg`]
    T : T_VEC_FA
        temperature
        [:math:`K`]
    Ne : T_VEC_FA
        electron number density
        [:math:`cm^{-3}`]

    Returns
    -------
    T_VEC_FA
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
    
    ## : ?
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
    
    ep_eV = ep / CST.eV2erg_
    rijk = _ElementUtil.sym_to_abun_( sym ) * \
           _numpy.exp( -11_604.52 * ep_eV / T ) / \
           _LTELib.Ufunc_(ion, T) / ss

    return rijk

@OVERLOAD
def xl_lte_(ion : T_STR, wl0 : T_FLOAT, ep : T_FLOAT, gf : T_FLOAT, T : T_FLOAT, Nh : T_FLOAT, Ne : T_FLOAT, dld : T_FLOAT) -> T_FLOAT: ...
@OVERLOAD
def xl_lte_(ion : T_STR, wl0 : T_FLOAT, ep : T_FLOAT, gf : T_FLOAT, T : T_ARRAY, Nh : T_ARRAY, Ne : T_ARRAY, dld : T_ARRAY) -> T_ARRAY: ...

def xl_lte_(ion : T_STR, wl0 : T_FLOAT, ep : T_FLOAT, gf : T_FLOAT, T : T_VEC_FA, Nh : T_VEC_FA, Ne : T_VEC_FA, dld : T_VEC_FA) -> T_VEC_FA:
    """line absorption coeff. at the line center (cm-1) at LTE
    stimulated emission is included

    Parameters
    ----------
    ion : T_STR
        element & ionization stage (Ca_I, Fe_II etc.)
    wl0 : T_FLOAT
        central wavelength of the line
        [:math:`cm`]
    ep : T_FLOAT
        excitation potential
        [:math:`erg`]
    gf : T_FLOAT
        gf-value
        [-]
    T : T_VEC_FA
        temperature
        [:math:`K`]
    Nh : T_VEC_FA
        hydrogen density
        [:math:`cm^{-3}`]
    Ne : T_VEC_FA
        electron density
        [:math:`cm^{-3}`]
    dld : T_VEC_FA
        doppler width
        [:math:`cm`] ??

    Returns
    -------
    T_VEC_FA
        [description]
    """
    wl0_AA = wl0 * 1.E+8 # [cm] -> [AA]
    dld_AA = dld * 1.E+8 # [cm] -> [AA] ?
    s0 = 4.99468E-21 * wl0_AA*wl0_AA / dld_AA
    stim = 1. - _numpy.exp( -1.4388E+8 / T / wl0_AA ) # correction of stimilated ...
    rijk = population_to_H_( ion, ep, T, Ne )

    return Nh * rijk * gf * s0 * stim

def line_prof_lte_(atmos : _Atmosphere.AtmosphereC1D, line : _Line.Line, 
                   dw : T_ARRAY, um : T_FLOAT = 1.) -> T_TUPLE[T_ARRAY, T_FLOAT, T_ARRAY] :
    """calculate LTE line profile

    Parameters
    ----------
    atmos : _Atmosphere.AtmosphereC1D
        1D Atmosphere struct
    line : _Line.Line
        Line struct
    dw : T_ARRAY
        wavelength array, distance from line center in cm
        [:math:`cm`]
    um : T_FLOAT, optional
        [description], by default 1.
        [-]

    Returns
    -------
    T_TUPLE[T_ARRAY, T_FLOAT, T_ARRAY]
        [description]
    """
    ion = line.ion
    sym, stage = _ElementUtil.ion_to_sym_and_stage_( ion )
    

    Z  : T_ARRAY = atmos.Z[:]
    Te : T_ARRAY = atmos.Te[:]
    Ne : T_ARRAY = atmos.Ne[:]
    Nh : T_ARRAY = atmos.Nh[:]
    Vt : T_ARRAY = atmos.Vt[:]
    Vd : T_ARRAY = atmos.Vd[:]

    wl0 = line.wl0
    ep  = line.ep
    gf  = line.gf

    nZ = Z.shape[0]
    nw = dw.shape[0]
    # find center wavelength
    j_min = dw.argmin()



    xc : T_ARRAY = _ContinuumOpacity.H_LTE_continuum_opacity_( Te[:], Ne[:], Nh[:], wl0 )
    # ic : [erg/cm^2/Sr/cm/s]
    ic, cntrbc, tauc = lte_integ_( Z[:], Te[:], xc[:], wl0, um )

    # am : atomic mass
    am = _ElementUtil.sym_to_mass_( sym )
    # dld : doppler width in [cm]
    dld : T_ARRAY = _BasicP.doppler_width_( wl0, Te[:], Vt[:], am )
    #xl0 :  opacity at line center
    xl0 : T_ARRAY = xl_lte_(ion , wl0, ep, gf, Te[:], Nh[:], Ne[:], dld[:])  
    # vv : doppler shift in unit of doppler width
    vv : T_ARRAY = _BasicP.dop_vel_to_shift_(wl0, Vd[:]) / dld[:]

    prof = _numpy.empty( nw, dtype=DT_NB_FLOAT )
    contrib : T_ARRAY
    for j in range(0,nw):
        v = dw[j] / dld   # normalize 
        #h, f = _Profile.hf(l.gamma, v-vv)
        h : T_FLOAT
        # hf,l.a,v-vv,h,f
        # xl : line opacity at a specific wavelength along line-of-sight
        xl = xl0[:] * h
        x = xl[:] + xc[:]
        intens1, cntrb, tau = lte_integ_(Z[:], Te[:], x[:], wl0, um)
        prof[j] = intens1
        # prof[j] = lteinteg(a.h,a.t,x,l.wl0,um,contrib=cntr,tau=tau);
        if j == j_min:
            contrib = cntrb

    
    return prof, ic, contrib


    
    




