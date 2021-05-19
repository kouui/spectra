
#-------------------------------------------------------------------------------
# function definition of opacity calculation
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------


from ..ImportAll import *
from . import Hydrogen

import numpy as _numpy
from numpy import sqrt as _sqrt
from numpy import exp as _exp

def thomson_scattering_(n_e : T_VEC_IFA) -> T_VEC_FA:
    r"""
    Thomson scattering of free electrons (non-relativistic).
    Returns absorption coefficient instead of absorption cross section.

    Parameters
    ----------

    n_e : T_VEC_IFA
        electron density, 
        [:math:`cm^{-3}`]

    Returns
    -------

    T_VEC_FA
        absorption coefficient, 
        [:math:`cm^{-1}`]

    Notes
    -----
    from [1]_.

    .. math:: \kappa = \sigma_{T} n_{e}

    where :math:`\sigma_{T}` is the absorption cross section for Thomson scattering

    .. math::
        \sigma_{T} = \frac{8 \pi e^{4}}{3 m_{e}^{2} c^{4}} = 6.6524 \times 10^{-25} \quad cm^{2}

    References
    ----------

    .. [1] Ivan Hubeny, Dimitri Mihalas, "Theory of Stellar Atmosphere:
        An Introduction to Astrophysical Non-equilibrium
        Quantitative Spectroscopic Analysis",
        Princeton University Press, pp. 149, 2015.
    """

    return 6.6524E-25 * n_e

#-------------------------------------------------------------------------------
# HI atom bound-free cross section from lower level ni
#-------------------------------------------------------------------------------


def hydrogenic_bf_cross_sec_n_(ni : T_VEC_IA, w : T_VEC_IFA, Z : T_INT) -> T_VEC_FA:
    r"""
    return bound-free cross section of a hydrogenic atom
    in unit of e-18 cm**2 ? cm**2

    Parameters
    ----------
    ni : T_VEC_IA, 
        principle quantum number of HI atom
        [:math:`-`]

    w  : T_VEC_IFA, 
        wavelength
        [:math:`cm`]

    Z : T_INT, 
        nucleus charge, +Ze
        [:math:`-`]

    Returns
    ----------
    alpha : T_VEC_FA
        absorption cross section
        [:math:`1E-18 cm^{2}`] ? [:math:`cm^{2}`]

    Notes
    -----
    absorption cross section formula from [3]_ page 188, Eq(7.92)

    from [1]_, [2]_, to calculate absorption coefficient given the cross section,

    - non-LTE ok
    - stimulated emission is not corrected
    - absorption coeff.:

            kbf = nk * alpha * (1.-exp(-hv/kT)) (/cm)

    Modification history:

    - k.ichimoto 15 jun.1987,	6 Jan.1992
    - k.ichimoto 19 Feb.1994
    - 2019.9.15   k.ichimoto from IDL ahic.pro

    References
    ----------
    .. [1] varnazza et al. (1976) ap.j.suppl. vol.30,1.
    .. [2] polynomial fitting for b-f gaunt factor is taken from Carbon & Gingerich (1969) "theory and observation of normal stellar atmosphere"  mit press
    .. [3] Ivan Hubeny, Dimitri Mihalas, "Theory of Stellar Atmosphere:
        An Introduction to Astrophysical Non-equilibrium
        Quantitative Spectroscopic Analysis",
        Princeton University Press, 2015.

    """

    ## ionization wavelength limit
    w_ni = CST.ni2cm_ * ni*ni

    if w > w_ni:
        alpha = 0.

    else:
        # ratio of the transition energy to ionization energy
        ratio = w_ni / w
        ratio_m1 = 1. / ratio

        # bound-free gaunt factor
        g_ni = Hydrogen.gaunt_factor_(ni, ratio)

        alpha = 7.904E-18 * (ni/Z/Z) * ratio_m1**3 * g_ni

    return alpha

#-----------------------------------------------------------------------------
# HI bound-free Cross section in LTE per 1 HI atom in unit of cm^2
#-----------------------------------------------------------------------------

def HI_bf_cross_sec_(Te : T_VEC_IFA, w : T_VEC_IFA) -> T_VEC_FA:
    r"""
    HI bound-free Cross section in LTE per 1 HI atom in unit of cm^2

    I recommend to use the
        sum_of_ni_from_1_to_5 {nHi * hydrogenic_bf_cross_sec_n_(ni, w, 1) * correction_factor}
    where
        correction_factor is ( 1. - exp(-1. * Cst.h_* Cst.c_/_w /_kT ) )

    Parameters
    -----------
    Te : float, 
        electron temperature
        [:math:`K`]

    w : float, 
        wavelength
        [:math:`cm`]

    Returns
    ----------
    alpha : float, 
        total absorption cross section
        [:math:`cm^{2}`]

    Notes
    -----
    from [1]_.

    HI bound-free opacity in LTE is

          k(T,wl) = n_HI * HI_bf_cross_sec_(Te,w) (/cm)

    - n_HI: hydrogen number density in cm^(-3)
    - stimulated emission is corrected
    - partition function of HI is assumed to be 2.0, which is valid for T < 2.e4 K

    Modification history:

    - 2019.9.15  K.Ichimoto  from IDL ahic.pro


    References
    ----------

    .. [1] varnazza et al. (1976) ap.j.suppl. vol.30, 1.

    .. [2] Robert J. Rutten, "Radiative Transfer in Stellar Atmosphere", 2003.
    """
    if Te > 2E4:
        
        raise ValueError("partition function of HI is assumed to be 2.0, which is only valid for T < 2E4 [K]")

    n_limit = int( _sqrt( w / CST.ni2cm_ ) ) + 1

    if n_limit > 7:
        alpha = 0.

    else:
        alpha = 0.
        kT = CST.k_ * Te
        fac = CST.E_Rydberg_ * ( 1. - n_limit**(-2) ) / kT
        #for n in range(_n_limit, min([_n_limit+3, 7])+1 ):
        for n in range(n_limit, 8 ):
            ## what is this `n*n` here ?
            alpha += n*n * _exp( -1. * fac ) * hydrogenic_bf_cross_sec_n_(n, w, 1)

        ## Cst.h_ * Cst.c_ / Cst.k ~ 1.438787
        ## correction factor
        alpha *= ( 1. - _exp(-1. * CST.h_* CST.c_/ w / kT ) )

    return alpha
#-----------------------------------------------------------------------------
# HI free-free CrossSection in LTE per 1 HI atom
#-----------------------------------------------------------------------------

def HI_ff_cross_sec_(Te : T_VEC_IFA, w : T_VEC_IFA) -> T_VEC_FA:
    r"""
    HI bound-free Cross section in LTE per 1 HI atom in unit of cm^2

    Parameters
    -----------
    Te : T_VEC_IFA, 
        electron temperature
        [:math:`K`]

    w : T_VEC_IFA, 
        wavelength
        [:math:`cm`]

    Returns
    ---------
    alpha : T_VEC_FA, 
        absorption cross section
        [:math:`cm^{2}`]

    Notes
    -----
    from [1]_, [2]_.

    HI free-free opacity in LTE is

          kff = n_HI * HI_ff_cross_sec_(Te,w) (/cm)

    - n_HI: hydrogen number density in cm^(-3)
    - stimulated emission is corrected
    - partition function of HI is assumed to be 2.0, which is valid for T < 2.e4 K

    Modification history:

    - 2019.9.15  K.Ichimoto  from IDL ahic.pro

    Original doc-string

    ;  HI free-free absorption coefficient per 1 proton and unit pe (dyne/cm**2)
        in unit of e-26
    ;  stimulated emission is not corrected
    ;  absorption coeff.:   kff = np*pe*ahiff*(1.-exp(-hv/kt))*1.e-26 (/cm)
    ;  ref.:  varnazza et al. (1976) ap.j.suppl. vol.30,1.
    ;         porinomial fitting for f-f gaunt factor is taken from
    ;         Carbon & Gingerich (1969) "theory and observation of
    ;           normal stellar atmosphere"  MIT press
    ;                           k.ichimoto 15 jun.1987,	6 Jan.1992
    ;                           k.ichimoto 19 Feb.1994

    References
    ----------
    .. [1] varnazza et al. (1976) ap.j.suppl. vol.30, 1.
    .. [2] porinomial fitting for f-f gaunt factor is taken from Carbon & Gingerich (1969) "theory and observation of
           normal stellar atmosphere"  MIT press

    """
    if Te > 2E4:
        
        raise ValueError("partition function of HI is assumed to be 2.0, which is only valid for T < 2E4 [K]")

    kT = CST.k_ * Te

    w_AA = w * 1.E8

    g0 = 1.0828 + 3.865E-6 * Te
    g1 = 7.564E-7 + (4.920E-10 - 2.482E-15 * Te) * Te
    g2 = 5.326E-12 + (-3.904E-15 + 1.8790E-20 * Te) * Te
    gff = g0 + ( g1 + g2 * w_AA ) * w_AA
    
    ahiff = 9.9264E-6 * Te**(-1.5) * w_AA**3 * gff

    ## exp( -157779./_Te ) ? exp( -157779./_w/_Te ) ?
    alpha = 0.66699 * Te**2.5 * _exp( -157779./Te ) * ahiff

    alpha *= 0.5 * ( 1. - _exp(-1. * CST.h_* CST.c_ / w  / kT ) )

    return alpha * 1.E-26

#-----------------------------------------------------------------------------
# H-minus (negative hidrogen) Cross Section per HI atom in unit of cm^2 (*LTE*)
#-----------------------------------------------------------------------------

def Hminus_cross_sec_(Te : T_VEC_IFA, w : T_VEC_IFA, Ne : T_VEC_IFA) -> T_VEC_FA:
    r"""
    H-minus (negative hidrogen) Cross Section per HI atom in unit of cm^2 (*LTE*)

    Parameters
    -----------
    Te  : float, 
        electron temperature
        [:math:`K`]

    w   : float, 
        wavelength
        [:math:`cm`]

    Ne  : float, 
        electron temperature
        [:math:`cm^{-3}`]

    Returns
    ------------
    alpha : float, 
        absorption cross section
        [:math:`cm^2`]

    Notes
    -----
    from [1]_, [2]_.

    H- opacity in LTE is

          khm = n_HI * Hminus_cross_sec_(T,wl,n_e) (/cm)

    - n_HI: hydrogen number density in cm^(-3)
    - stimulated emission is corrected
    - partition function of HI is assumed to be 2.0, which is valid for T < 2.e4 K

    Modification history:

    - k.ichimoto 18 jun. 1987, 	6 Jan. 1992
    - k.ichimoto  19 Feb.1994
    - 2019.9.15  K.Ichimoto  from IDL ahic.pro

    References
    ----------
    .. [1] varnazza et al. (1976) ap.j.suppl. vol.30, 1.
    .. [2] porinomial fitting for f-f gaunt factor is taken from Carbon & Gingerich (1969) "theory and observation of
           normal stellar atmosphere"  MIT press
    """
    w_AA = w * 1.E8
    theta = 5039.778 / Te
    w_AA3 = 1.E-3 * w_AA

    #------ bound-free
    if w_AA3 > 16.419:
        sigma = 0.
    elif w_AA3 < 16.419 and w_AA3 > 14.2:
        xl = 16.419 - w_AA3
        sigma = (0.269818+(0.220190+(-0.0411288+0.00273236)*xl)*xl)*xl
    else:# _w_AA3 <= 14.2:
        sigma = 0.00680133 + (0.178708+(0.16479 + (-0.0204842 + 5.95244e-4 * w_AA3) * w_AA3) *w_AA3)*w_AA3

    kbf = 0.41590 * theta**2.5 * _exp(1.738*theta) * (1.- _exp(-28.5486*theta/w_AA3)) * sigma

    #------ free-free
    g0 = 0.005366+(-0.011493+0.027039*theta)*theta
    g1 = -3.2062+(11.924-5.939*theta)*theta
    g2 = -0.40192+(7.0355-0.34592*theta)*theta
    kff = g0 + (g1+g2*w_AA3)*w_AA3 * 1.E-3

    ahm = kbf + kff
    pe = 1.38066E-16 * Ne * Te
    alpha = ahm * pe

    return alpha * 1.E-26

#-----------------------------------------------------------------------------
# HI Rayleigh scattering CrossSection per 1 HI atom in unit cm^2
#-----------------------------------------------------------------------------

def HI_rayleigh_cross_sec_(w : T_VEC_IFA) -> T_VEC_FA:
    r"""
    return Rayleigh scattering cross section of a HI atom in unit of cm^2

    Parameters
    ----------
    w : T_VEC_IFA, 
        wavelength
        [:math:`cm`]

    Returns
    ----------
    alpha : T_VEC_FA, 
        absorption cross section
        [:math:`cm^2`]

    Notes
    -----
    from [1]_, [2]_.

    HI Rayleigh opacity in LTE is

         khray = n_HI * HI_rayleigh_cross_sec_(wl) (/cm)

    - n_HI: hydrogen number density in cm^(-3)
    - stimulated emission is corrected
    - partition function of HI is assumed to be 2.0, which is valid for T < 2.e4 K

    modification history:

    - k.ichimoto 18 jun. 1987, 	6 Jan. 1992
    - k.ichimoto  19 Feb.1994
    - 2019.9.10  K.Ichimoto  from IDL avray.pro

    References
    ----------
    .. [1] Gingerich SAO special report 167,21,1964.
    .. [2] porinomial fitting for f-f gaunt factor is taken from Carbon & Gingerich (1969) "theory and observation of
           normal stellar atmosphere"  MIT press
    """
    w_AA = w * 1.E8
    if w_AA < 1026.:
       alpha = 0.
    else:
       w_AA2_m1 = 1. / (w_AA * w_AA)
       alpha = 5.799E13 * w_AA2_m1**2 + 1.422E20 * w_AA2_m1**3 + 2.784 * w_AA2_m1**4

    return alpha * 1.E-26

#-----------------------------------------------------------------------------
# H2+ CrossSection per 1 HI atom and per 1 H+ in unit of cm^5?
#-----------------------------------------------------------------------------

#;/* -----  data from Carbon & Gingerich (1969)  ----- */
_AVH2P_E = _numpy.array(
        [ 3.0,    2.852,  2.58,   2.294,  2.023,  1.774,
        1.547,  1.344,  1.165,  1.007,  0.8702, 0.7516,
        0.6493, 0.5610, 0.4848, 0.4198, 0.3620, 0.3128,
        0.2702, 0.2332, 0.2011, 0.1732, 0.1491, 0.1281,
        0.1100, 0.09426,0.0809, 0.06906,0.05874,0.04994,
        0.04265,0.03635,0.0308, 0.026,  0.02195,0.01864,
        0.01581,0.01332,0.01118,0.00938,0.00793,0.00669,
        0.00561,0.00469,0.00392,0.0033 ], dtype=T_FLOAT)
_AVH2P_UP = _numpy.array(
        [  85.,       9.99465,   4.97842,   3.28472,   2.41452,
        1.87038,   1.48945,   1.20442,   0.98279,   0.80665,
        0.66493,   0.54997,   0.45618,   0.37932,   0.31606,
        0.26382,   0.22057,   0.18446,   0.15473,   0.12977,
        1.08890e-1,9.14000e-2,7.67600e-2,6.44500e-2,5.41200e-2,
        4.54000e-2,3.81000e-2,3.19500e-2,2.67600e-2,2.23700e-2,
        1.86900e-2,1.56100e-2,1.30200e-2,1.08300e-2,8.99000e-3,
        7.45000e-3,6.15000e-3,5.08000e-3,4.16000e-3,3.42000e-3,
        2.77000e-3,2.21000e-3,1.78000e-3,1.45000e-3,1.24000e-3,
        1.14000e-3 ], dtype=T_FLOAT)
_AVH2P_US = _numpy.array(
        [ -85.,      -7.1426,   -2.3984,   -0.99032,  -0.39105,
        -0.09644,  0.05794,   0.13996,   0.18186,   0.20052,
        0.20525,   0.20167,   0.19309,   0.18167,   0.16871,
        0.15511,   0.14147,   0.12815,   0.11542,   0.10340,
        0.09216,   8.18000e-2,7.22900e-2,6.36700e-2,5.58400e-2,
        4.88400e-2,4.25700e-2,3.69900e-2,3.20700e-2,2.77500e-2,
        2.39400e-2,2.06100e-2,1.77200e-2,1.52200e-2,1.30500e-2,
        1.11900e-2,9.58000e-3,8.21000e-3,7.01000e-3,6.00000e-3,
        5.11000e-3,4.35000e-3,3.72000e-3,3.22000e-3,2.86000e-3,
        2.63000e-3 ], dtype=T_FLOAT)
_AVH2P_FR = _numpy.array(
        [ 0.,      4.30272e-18,1.51111e-17,4.02893e-17,8.89643e-17,
        1.70250e-16,2.94529e-16,4.77443e-16,7.25449e-16,1.06238e-15,
        1.50501e-15,2.08046e-15,2.82259e-15,3.76256e-15,4.93692e-15,
        6.38227e-15,8.17038e-15,1.02794e-14,1.28018e-14,1.57371e-14,
        1.91217e-14,2.30875e-14,2.75329e-14,3.27526e-14,3.85481e-14,
        4.52968e-14,5.18592e-14,5.99825e-14,6.92092e-14,7.94023e-14,
        9.01000e-14,1.01710e-13,1.14868e-13,1.29969e-13,1.46437e-13,
        1.63042e-13,1.81440e-13,2.02169e-13,2.25126e-13,2.49637e-13,
        2.73970e-13,3.00895e-13,3.30827e-13,3.64140e-13,3.99503e-13,
        4.34206e-13 ], dtype=T_FLOAT)

def _avH2p_(Te : T_VEC_IFA, w : T_VEC_IFA) -> T_VEC_FA:
    r"""

    H2+ opacity per 1 HI atom and per 1 H+

    Parameters
    ----------
    T  : T_VEC_IFA, 
        temperature
        [:math:`K`]

    wl : T_VEC_IFA, 
        wavelength
        [:math:`cm`]

    Returns
    ---------
    T_VEC_FA
        ?
        [:math:`cm^{5}`]

    Notes
    -----
    from [1]_, [2]_.

    - stimulated emission is corrected.

    absorption coeff,
        kh2 = nhi * np * _avh2p_, [:math:`cm^{-1}`]

    Modification history:

    - k.ichimoto 18 jun. 1987, 	6 Jan. 1992
    - k.ichimoto  19 Feb.1994
    - 2019.9.11  K.Ichimoto  from IDL avh2p.pro


    References
    ----------
    .. [1] Gingerich SAO special report 167,21,1964.
    .. [2] Carbon & Gingerich (1969) "theory and observation of normal stellar atmosphere"  MIT press
    """
    w_AA = w * 1.E8


    ev = 911.3047 / w_AA
    Tk = 6.3348e-6 * Te

    for i in range(_AVH2P_E.shape[0]):
        if _AVH2P_E[i] <= ev:
            break

    if i == _AVH2P_E.shape[0]-1:
        n = 45
    else:
        n = i

    d   = (ev-_AVH2P_E[n])/(_AVH2P_E[n+1]-_AVH2P_E[n])
    usq = _AVH2P_US[n] + (_AVH2P_US[n+1]-_AVH2P_US[n])*d
    upq = _AVH2P_UP[n] + (_AVH2P_UP[n+1]-_AVH2P_UP[n])*d
    frq = _AVH2P_FR[n] + (_AVH2P_FR[n+1]-_AVH2P_FR[n])*d

    alpha = _numpy.abs( frq * ( _exp(usq/Tk) - _exp(-upq/Tk) ) )

    return alpha * 1.E-26

#-----------------------------------------------------------------------------
# H2+ CrossSection per 1 HI atom n unit  cm^2 in LTE
#-----------------------------------------------------------------------------
def H2p_cross_sec_(Te : T_VEC_IFA, w : T_VEC_IFA, Ne : T_VEC_IFA, 
           Nh : T_VEC_IFA) -> T_VEC_FA:
    r"""
    return H2+ CrossSection per 1 HI atom n unit of cm^2 in LTE

    Parameters
    ----------
    Te  : T_VEC_IFA, 
        electron temperature
        [:math:`K`]

    w : T_VEC_IFA, 
        wavelength
        [:math:`cm`]

    Ne : T_VEC_IFA, 
        electron number density
        [:math:`cm^{-3}`]

    Nh : T_VEC_IFA, 
        hydrogen number density
        [:math:`cm^{-3}``]

    Returns
    --------
    alpha : T_VEC_FA, 
        absorption cross srctions
        [:math:`cm^{2}`]

    Notes
    -----
    H2+ opacity in LTE is

         kH2p = n_HI * H2p_cross_sec_(T,wl,n_e,n_H) (/cm)

    use `_avH2p_(T,wl)`

    Modification history:

    - 2019.9.15  K.Ichimoto
    """
    Pe = 1.38066e-16 * Ne * Te

    #/*  -----   solving LTE Saha's eq. for hydrogen   -----	*/
    q = 2.0*2.07e-16 * Te**(-1.5) * 10.**(5040.*13.6/Te) * Ne
    Np = 1./(q+1.) * Nh

    alpha = Np * _avH2p_(Te, w)
    return alpha

#-----------------------------------------------------------------------------
# hydrogen LTE continuum opacity in cm^{-1}
#-----------------------------------------------------------------------------
def H_LTE_opacity(Te : T_VEC_IFA, Ne : T_VEC_IFA, Nh : T_VEC_IFA, 
                  w : T_VEC_IFA) -> T_VEC_FA:
    r"""
    LTE continuum opacity in cm-1, incruding HI, H-, H2, rayleigh scat.

    Parameters
    ----------
    Te  : T_VEC_IFA, 
        electron temperature
        [:math:`K`]

    Ne  : T_VEC_IFA, 
        electron number density
        [:math:`cm^{-3}`]

    Nh  : T_VEC_IFA, 
        hydrogen number density
        [:math:`cm^{-3}`]

    w   : T_VEC_IFA, 
        wavelength
        [:math:`cm`]

    Returns
    --------
    kappa : T_VEC_FA, 
        opacity
        [:math:`cm^{-1}`]

    Notes
    -----
    hydrogen total opacity in LTE is [1]_.

         kappa = H_LTE_opacity_(_Te, _Ne, _Nh, _w) (/cm)

    - n_HI: hydrogen number density in cm^(-3)
    - stimulated emission is corrected
    - partition function of HI is assumed to be 2.0, which is valid for T < 2.e4 K

    Modification history:

    - k.ichimoto 18 jun. 1987, 	6 Jan. 1992
    - k.ichimoto  19 Feb.1994
    - 2019.9.15  K.Ichimoto  from IDL avray.pro

    References
    ----------
    .. [1] Gingerich SAO special report 167,21,1964. porinomial fitting for f-f gaunt factor is taken from Carbon & Gingerich (1969) "theory and observation of normal stellar atmosphere"  MIT press
    """
    #/*  -----   solving LTE Saha's eq. for hydrogen   -----	*/
    q = 2. * 2.07E-16 * Te**(-1.5) * 10.**(5040.*13.6/Te) * Ne
    N_HI =  q /(q+1.0) * Nh
    kappa = N_HI * (
             HI_bf_cross_sec_(Te, w) +
             HI_ff_cross_sec_(Te, w) +
             Hminus_cross_sec_(Te, w, Ne) +
             H2p_cross_sec_(Te, w, Ne, Nh) +
             HI_rayleigh_cross_sec_(w)
             )

    return kappa

#-------------------------------------------------------------------------------
# numba optimization
#-------------------------------------------------------------------------------

if CFG._IS_JIT:

    thomson_scattering_ = nb_vec(**NB_VEC_KWGS) (thomson_scattering_)

    hydrogenic_bf_cross_sec_n_ = nb_vec(**NB_VEC_KWGS) (hydrogenic_bf_cross_sec_n_)
    HI_bf_cross_sec_ = nb_vec(**NB_VEC_KWGS) (HI_bf_cross_sec_)
    HI_ff_cross_sec_ = nb_vec(**NB_VEC_KWGS) (HI_ff_cross_sec_)
    Hminus_cross_sec_= nb_vec(**NB_VEC_KWGS) (Hminus_cross_sec_)
    _avH2p_ = nb_vec(**NB_VEC_KWGS) (_avH2p_)
    H2p_cross_sec_ = nb_vec(**NB_VEC_KWGS) (H2p_cross_sec_)
    H_LTE_opacity = nb_vec(**NB_VEC_KWGS) (H_LTE_opacity)

else:

    hydrogenic_bf_cross_sec_n_ = np_vec(hydrogenic_bf_cross_sec_n_, **NP_VEC_KWGS)
    HI_bf_cross_sec_ = np_vec(HI_bf_cross_sec_, **NP_VEC_KWGS)
    HI_ff_cross_sec_ = np_vec(HI_ff_cross_sec_, **NP_VEC_KWGS)
    Hminus_cross_sec_= np_vec(Hminus_cross_sec_, **NP_VEC_KWGS)
    _avH2p_ = np_vec(_avH2p_, **NP_VEC_KWGS)
    H2p_cross_sec_ = np_vec(H2p_cross_sec_, **NP_VEC_KWGS)
    H_LTE_opacity = np_vec(H_LTE_opacity, **NP_VEC_KWGS)
