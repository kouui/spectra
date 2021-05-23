#-------------------------------------------------------------------------------
# function definition of process of hydrogen atom
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#        - `Gaunt_factor_Gingerich_cm_()`, `w_um = w * 1E5` -> `w_um = w * 1.E4`
# 0.0.1
#    2020/11/10   u.k.
#        - `PI_cross_section_cm()` and `PI_cross_section()`, if `x<1.0` then cross section `alpha=0.`
#-------------------------------------------------------------------------------

from ..ImportAll import *
from ..Math import Special

import numpy as _numpy
from numpy import sqrt as _sqrt
from numpy import exp as _exp
from numpy import log as _log


from debtcollector import removals as _removals  # type: ignore

@OVERLOAD
def ratio_Etran_to_Eionize_(ni : T_INT, w : T_FLOAT) -> T_FLOAT : ...
@OVERLOAD
def ratio_Etran_to_Eionize_(ni : T_ARRAY, w : T_FLOAT) -> T_ARRAY : ...
@OVERLOAD
def ratio_Etran_to_Eionize_(ni : T_ARRAY, w : T_ARRAY) -> T_ARRAY : ...
@OVERLOAD
def ratio_Etran_to_Eionize_(ni : T_INT, w : T_ARRAY) -> T_ARRAY : ...

def ratio_Etran_to_Eionize_(ni : T_VEC_IA, w : T_VEC_FA) -> T_VEC_FA:
    """
    Compute the "ratio of transition energy to ionization energy"

    Parameters
    -----------
    ni : T_VEC_IA
        principal quantum number of lower level
        [:math:`-`]

    w : T_VEC_FA
        wavelength
        [:math:`cm`]

    Returns
    --------
    ratio : T_VEC_FA
        ratio of transition energy to ionization energy
        [:math:`-`]
    """
    # ionization energy
    #Eik = E_Rydberg_ * (1./ni**2)
    # transition energy
    #E_tran = h_ * c_ / w

    #ratio = E_tran / Eik

    w_limit = CST.ni2cm_ * ni * ni
    ratio = w_limit / w
    return ratio

#-----------------------------------------------------------------------------
# Gaunt factor
#-----------------------------------------------------------------------------
@_removals.remove
def gaunt_factor_gingerich_cm_(ni : T_VEC_IA, w : T_VEC_FA) -> T_VEC_FA:
    """
    Gaunt factor

    Parameters
    -----------
    ni : T_VEC_IA, 
        principal quantum number of lower level
        [:math:`-`]

    w : T_VEC_FA, 
        wavelength
        [:math:`cm`]

    Returns
    --------
    g : T_VEC_FA
        Gaunt factor
        [:math:`-`]

    References
    ------------
    .. [1] Gingerich, March 1964
    """
    # wavelength in ?
    #w_um = w * 1E5
    w_um = w * 1.E4

    if ni == 1:
        C1, C2, C3 = 0.9916, 9.068E-3, -0.2524
    elif ni == 2:
        C1, C2, C3 = 1.105, -7.922E-2, 4.536E-3
    elif ni == 3:
        C1, C2, C3 = 1.101, -3.290E-2, 1.152E-3
    elif ni == 4:
        C1, C2, C3 = 1.101, -1.923E-2, 5.110E-4
    elif ni == 5:
        C1, C2, C3 = 1.102, -0.01304, 2.638E-4
    elif ni == 6:
        C1, C2, C3 = 1.0986, -0.00902, 1.367E-4
    else:
        C1, C2, C3 = 1., 0., 0.

    g = C1 + ( C2 + C3 * w_um ) * w_um
    return g

@_removals.remove
def gaunt_factor_gingerich_(ni : T_VEC_IA, x : T_VEC_IFA) -> T_VEC_FA:
    r"""
    Gaunt factor

    Parameters
    -----------
    ni : T_VEC_IA
        principal quantum number of lower level
        [:math:`-`]

    x : T_VEC_IFA
        ratio of the transition energy to the ionization energy of the lower level.
        [:math:`-`]

    Returns
    --------
    g : T_VEC_FA
        Gaunt factor
        [:math:`-`]

    Notes
    ------
    refer to [1]_.

    `x` is the ratio of the transition energy to the ionization energy.

    ionization energy is `Cst.E_Rydberg_ * (1/ni**2)`

    transition and `hc/w` for bound-free transition where `w` wavelength.

    References
    ------------
    .. [1] Gingerich, March 1964
    """
    # ionization energy
    Eik = CST.E_Rydberg_ * (1./ni**2)
    # transition energy
    E_tran = x * Eik
    # wavelength
    w = CST.h_ * CST.c_ / E_tran

    g = gaunt_factor_gingerich_cm_(ni, w)
    return g

def gaunt_factor_coe_(i : T_VEC_IA, ni : T_VEC_IA) -> T_VEC_FA:
    """
    coefficients to calculate Gaunt factor

    Parameters
    -----------
    i  : T_VEC_IA, 
        order of the coefficient, gi
        [:math:`-`]

    ni : T_VEC_IA, 
        principal quantum number of lower level
        [:math:`-`]

    Returns
    --------
    gi : T_VEC_FA, 
        Gaunt factor coefficient
        [:math:`-`]

    Notes
    ------
    refer to [1]_ Table 1.


    References
    ------------
    .. [1] L. C. Johnson, "Approximations for collisional and
           radiative transition rates in atomic hydrogen",
           Astrophysical Journal, vol. 174, p.227, May 1972.
           1972ApJ...174..227J
    """
    gi : T_FLOAT
    if ni == 1:
        gi = (1.1330, -0.4059, 0.07014)[i] # type: ignore
    elif ni == 2:
        gi = (1.0785, -0.2319, 0.02947)[i] # type: ignore
    else:
        gi = (
            0.9935 + ( 0.2328 - 0.1296 / ni ) / ni,
            - ( 0.6282 - ( 0.5598 - 0.5299 / ni ) / ni ) / ni,
            ( 0.3887 - ( 1.181 - 1.470 / ni ) / ni ) / ni / ni,
        )[i] # type: ignore

    return gi

def gaunt_factor_(ni : T_VEC_IA, x : T_VEC_IFA) -> T_VEC_FA:
    """
    Gaunt factor

    Parameters
    -----------
    ni : T_VEC_IA, 
        principal quantum number of lower level
        [:math:`-`]

    x : T_VEC_IFA, 
        ratio of the transition energy to the ionization energy of the lower level.
        [:math:`-`]


    Returns
    --------
    g : T_VEC_FA, 
        Gaunt factor
        [:math:`-`]

    Notes
    ------
    refer to [1]_ Eq(4).

    `x` is the ratio of the transition energy to the ionization energy.

    ionization energy is `Cst.E_Rydberg_ * (1/ni**2)`

    transition energy is `Cst.E_Rydberg_ * (1/ni**2 - 1/nj**2)` for bound-bound
    transition and `hc/w` for bound-free transition where `w` wavelength




    References
    ------------
    .. [1] L. C. Johnson, "Approximations for collisional and
           radiative transition rates in atomic hydrogen",
           Astrophysical Journal, vol. 174, p.227, May 1972.
           1972ApJ...174..227J
    """

    g0 = gaunt_factor_coe_(0, ni)
    g1 = gaunt_factor_coe_(1, ni)
    g2 = gaunt_factor_coe_(2, ni)
    # Gaunt factor
    g = g0 + ( g1 + g2 / x ) / x

    return g

#-----------------------------------------------------------------------------
# Einstein Aji coefficient
#-----------------------------------------------------------------------------
def absorption_oscillator_strength_(ni : T_VEC_IA, nj : T_VEC_IA) -> T_VEC_FA:
    r"""
    absorption oscillator strength given by (Bethe and Salpeter 1957)
    without the correction of Gaunt factor.

    Parameters
    -----------
    ni : T_VEC_IA, 
        principal quantum number of lower level
        [:math:`-`]

    nj : T_VEC_IA, 
        principal quantum number of upper level
        [:math:`-`]

    Returns
    --------
    fij : T_VEC_FA 
        absorption oscillator strength
        [:math:`-`]

    Notes
    ------
    refer to [1]_ Eq(1).


    References
    ------------
    .. [1] L. C. Johnson, "Approximations for collisional and
           radiative transition rates in atomic hydrogen",
           Astrophysical Journal, vol. 174, p.227, May 1972.
           1972ApJ...174..227J
    """
    coef = 32. / (3. * CST.sqrt3_ * CST.pi_ )
    x = 1. - (ni / nj)**2

    fij = coef * ni / ( x * nj )**3

    return fij

def einstein_A_coefficient_(ni : T_VEC_IA, nj : T_VEC_IA) -> T_VEC_FA:
    r"""
    Einstein coefficient for spontaneous emission in the hydrogen atom

    Parameters
    ------------
    ni : T_VEC_IA
        principal quantum number of lower level
        [:math:`-`]

    nj : T_VEC_IA,  
        principal quantum number of upper level
        [:math:`-`]

    Returns
    ---------
    Aji : T_VEC_FA, 
        Einstein Aji coefficient
        [:math:`s^{-1}`]

    Notes
    ------
    refer to [1]_ page 117 Eq(5.7) and Eq(5.8), page 135 Eq(5.105)

    .. math:: A_{ji} = C_{1} \cdot \frac{g_{i}}{g_{j}} \lambda^{-2} f_{ij}

    where

    .. math:: C_{1} = \frac{8 \pi^{2} e^{2} }{m_{e} c} \sim 0.6670


    References
    -----------

    .. [1] Ivan Hubeny, Dimitri Mihalas, "Theory of Stellar Atmosphere:
        An Introduction to Astrophysical Non-equilibrium
        Quantitative Spectroscopic Analysis",
        Princeton University Press, 2015.

    """
    # Gaunt factor
    x = 1. - (ni / nj)**2
    g = gaunt_factor_(ni, x)

    # corrected absorption oscillator strength
    fij = absorption_oscillator_strength_(ni, nj) * g

    # excitation energy
    Eij = CST.E_Rydberg_ * (1./ni**2 - 1./nj**2)

    # wavelength
    w = CST.h_ * CST.c_ / Eij

    # constant factor
    C1 = 0.667025 # 8 * Cst.pi_**2 * Cst.e_**2 / Cst.me_ / Cst.c_

    Aji = (ni/nj)**2 * C1 / (w*w) * fij

    return Aji

#-----------------------------------------------------------------------------
# Collisional Excitation rate coefficient
#-----------------------------------------------------------------------------
@OVERLOAD
def CE_rate_coe_(ni : T_INT, nj : T_INT, Te : T_FLOAT) -> T_FLOAT: ...
@OVERLOAD
def CE_rate_coe_(ni : T_INT, nj : T_INT, Te : T_ARRAY) -> T_ARRAY: ...
@OVERLOAD
def CE_rate_coe_(ni : T_ARRAY, nj : T_ARRAY, Te : T_FLOAT) -> T_ARRAY: ...

def CE_rate_coe_(ni : T_VEC_IA, nj : T_VEC_IA, Te : T_VEC_IFA) -> T_VEC_FA:
    r"""
    Collisional Excitation rate coefficient qij for the hydrogen atom

    Parameters
    ------------
    ni : T_VEC_IA,
        principal quantum number of lower level
        [:math:`-`]

    nj : T_VEC_IA,  
        principal quantum number of upper level
        [:math:`-`]

    Te : T_VEC_IFA, 
        electron temperature
        [:math:`K`]

    Returns
    ---------
    qij : T_VEC_FA, 
        Colisional excitation rate coefficient
        [:math:`cm^{3}s^{-1}`]

    Notes
    ------
    refer to [1]_ Eq(36) and [2] page 275 Eq(9.51), Eq(9.54).

    the :math:`S_{e}(n,n')` is the Collisional Excitation rate coefficient :math:`q_{ij}`,

    with

    .. math:: n_{i} C_{ij} = n_{i} n_{e} q_{ij}


    References
    -----------
    .. [1] L. C. Johnson, "Approximations for collisional and
           radiative transition rates in atomic hydrogen",
           Astrophysical Journal, vol. 174, p.227, May 1972.
           1972ApJ...174..227J

    .. [2] Ivan Hubeny, Dimitri Mihalas, "Theory of Stellar Atmosphere:
        An Introduction to Astrophysical Non-equilibrium
        Quantitative Spectroscopic Analysis",
        Princeton University Press, 2015.

    """

    kT = CST.k_ * Te

    # Gaunt factor
    x = 1. - (ni / nj)**2
    x_m = 1. / x
    g = gaunt_factor_(ni, x)

    # corrected absorption oscillator strength
    fij = absorption_oscillator_strength_(ni, nj) * g

    # Eq(11)
    Aij = 2. * ni**2 * x_m * fij

    # Eq(25, 26, 31, 32)
    if ni >= 2:
        ni_m_ = 1 / ni
        bi = ni_m_ * ( 4.0 + ni_m_ * ( -18.63 + ni_m_ * ( 36.24 - ni_m_ * 28.09 ) ) )
        ri = 1.94 * ni_m_**(1.57)
    else:
        bi = -0.603
        ri = 0.45

    # Eq(23)
    Bij = 4. * ni * (ni/nj)**3 * x_m**2 * ( 1. + x_m * ( 4/3 + x_m * bi ) )

    # excitation energy
    Eij = CST.E_Rydberg_ * (1./ni**2 - 1./nj**2)

    # Eq(37)
    y = Eij / kT
    # Eq(30)
    rij = ri * x
    # Eq(38)
    z = rij + y

    term1 = Aij * ( ( 1./y + 0.5 ) * Special.E1_(y) - ( 1./z + 0.5 ) * Special.E1_(z) )
    term2 = ( Bij - Aij *  _log( 2. * ni * ni * x_m ) ) * ( Special.E2_(y) / y - Special.E2_(z) / z )

    # Eq(36)
    Sij = CST.C0_ * Te**(0.5) * 2 * ni * ni * x_m * y * y * (term1 + term2)
    qij = Sij

    return qij

#-----------------------------------------------------------------------------
# Collisional Ionization rate coefficient
#-----------------------------------------------------------------------------
@OVERLOAD
def CI_rate_coe_(ni : T_INT, Te : T_FLOAT) -> T_FLOAT: ...
@OVERLOAD
def CI_rate_coe_(ni : T_ARRAY, Te : T_FLOAT) -> T_ARRAY: ...
@OVERLOAD
def CI_rate_coe_(ni : T_INT, Te : T_ARRAY) -> T_ARRAY: ...

def CI_rate_coe_(ni : T_VEC_IA, Te : T_VEC_IFA) -> T_VEC_FA:
    r"""
    Collisional ionization rate coefficient qik for the hydrogen atom

    Parameters
    ------------
    ni : T_VEC_IA,  
        principal quantum number of lower level
        [:math:`-`]

    Te : T_VEC_IFA, 
        electron temperature
        [:math:`K`]

    Returns
    ---------
    qik : T_VEC_FA, 
        Colisional ionization rate coefficient
        [:math:`cm^{3}s^{-1}`]

    Notes
    ------
    refer to [1]_ Eq(39) and [2] page 275 Eq(9.51), Eq(9.54).

    the :math:`S_{e}(n,n')` is the Collisional ionization rate coefficient :math:`q_{ik}`,

    with

    .. math:: n_{i} C_{ik} = n_{i} n_{e} q_{ik}


    References
    -----------
    .. [1] L. C. Johnson, "Approximations for collisional and
           radiative transition rates in atomic hydrogen",
           Astrophysical Journal, vol. 174, p.227, May 1972.
           1972ApJ...174..227J

    .. [2] Ivan Hubeny, Dimitri Mihalas, "Theory of Stellar Atmosphere:
        An Introduction to Astrophysical Non-equilibrium
        Quantitative Spectroscopic Analysis",
        Princeton University Press, 2015.

    """

    kT = CST.k_ * Te

    # Eq(20)
    i_arr = _numpy.array([0,1,2],dtype=DT_NB_INT)
    gi_arr = gaunt_factor_coe_(i_arr, ni)
    Ai = 32. / (3. * CST.sqrt3_ * CST.pi_) * ni * ( gi_arr / (i_arr[:]+3) ).sum()

    # Eq(25, 26, 31, 32)
    if ni >= 2:
        ni_m_ = 1 / ni
        bi = ni_m_ * ( 4.0 + ni_m_ * ( -18.63 + ni_m_ * ( 36.24 - ni_m_ * 28.09 ) ) )
        ri = 1.94 * ni_m_**(1.57)
    else:
        bi = -0.603
        ri = 0.45

    # Eq(24)
    Bi = 2./3. * ni * ni * (5. + bi)

    # ionization energy
    Eik = CST.E_Rydberg_ * (1./ni**2)

    # Eq(40)
    y = Eik / kT
    # Eq(41)
    z = ri + y

    E1_y = Special.E1_(y)
    E1_z = Special.E1_(z)
    E2_y = Special.E2_(y)
    E2_z = Special.E2_(z)

    # Eq(42)
    xi_y = Special.E0_(y) - 2.*E1_y + E2_y
    xi_z = Special.E0_(z) - 2.*E1_z + E2_z

    term1 = Ai * ( 1./y * E1_y - 1./z * E1_z )
    term2 = ( Bi - Ai*_log(2*ni*ni) ) * ( xi_y - xi_z )
    # Eq(39)
    Sik = CST.C0_ * Te**(0.5) * 2 * ni * ni * y * y * (term1 + term2)
    qik = Sik

    return qik

@_removals.remove
def CI_rate_coe_clark_(ni : T_VEC_IA, Te : T_VEC_IFA) -> T_VEC_FA:
    r"""
    Collisional ionization rate coefficient qik for the hydrogen atom

    Parameters
    ------------
    ni : T_VEC_IA,  
        principal quantum number of lower level
        [:math:`-`]

    Te : T_VEC_IFA, 
        electron temperature
        [:math:`K`]

    Returns
    ---------
    qik : T_VEC_FA, 
        Colisional ionization rate coefficient
        [:math:`cm^{3}s^{-1}`]

    Notes
    ------
    refer to [1]_ Eq(8) and [2] page 275 Eq(9.51), Eq(9.54).

    the :math:`C(y)` is the Collisional ionization rate coefficient :math:`q_{ik}`,

    with

    .. math:: n_{i} C_{ik} = n_{i} n_{e} q_{ik}


    References
    -----------
    .. [1] Clark, Abdallah & Mann, "Integral and differential
           cross sections for electron impact ionization",
           Astrophysical Journal, vol. 381, 597-600, Nov 1972.
           1991ApJ...381...597C

    .. [2] Ivan Hubeny, Dimitri Mihalas, "Theory of Stellar Atmosphere:
        An Introduction to Astrophysical Non-equilibrium
        Quantitative Spectroscopic Analysis",
        Princeton University Press, 2015.

    """

    kT = CST.k_ * Te

    C1 =  1.53690
    C2 =  0.99656
    C3 = -0.61916
    C4 =  2.44630
    C5 = -2.47730
    C6 =  3.21510
    C7 = -1.45120
    C8 =  1.72300
    C9 = -0.47075

    G = ( ni - 1 ) * ( 4 * ni + 1 ) / ( 6 * ni )

    # ionization energy
    Eik = CST.E_Rydberg_ * (1./ni**2)

    y = Eik / kT
    TeV = Te * CST.K2eV_
    F = 5.89E-9 * _sqrt( TeV ) * y * ni**4# / Eik**2 # mistake in the paper?

    E_y  = _exp(-y)
    E1_y = Special.E1_(y)
    E2_y = Special.E2_(y)

    # Eq(8)
    term1 = ( C1 + (C2 + C3 * G) / ni ) * E1_y
    term2 = ( C4 + (C5 + C6 * G) / ni ) * (E_y - y * E1_y)
    term3 = ( C7 + (C8 + C9 * G) / ni ) * (E_y - 2 * y * E1_y + y * E2_y)
    C_y = F * (term1 + term2 + term3)

    qik = C_y
    return qik

#-----------------------------------------------------------------------------
# Photoionization cross section
#-----------------------------------------------------------------------------
@OVERLOAD
def PI_cross_section_cm_(ni : T_INT, w : T_FLOAT, Z : T_INT) -> T_FLOAT: ...
@OVERLOAD
def PI_cross_section_cm_(ni : T_ARRAY, w : T_FLOAT, Z : T_INT) -> T_ARRAY: ...
@OVERLOAD
def PI_cross_section_cm_(ni : T_ARRAY, w : T_ARRAY, Z : T_INT) -> T_ARRAY: ...
@OVERLOAD
def PI_cross_section_cm_(ni : T_INT, w : T_ARRAY, Z : T_INT) -> T_ARRAY: ...

def PI_cross_section_cm_(ni : T_VEC_IA, w : T_VEC_IFA, Z : T_VEC_IA) -> T_VEC_FA:
    r"""
    Photoionization cross-section for hydrogen from lower level ni at wavelength w.

    Parameters
    ------------
    ni : T_VEC_IA,
        principal quantum number of lower level
        [:math:`-`]

    w : T_VEC_IFA,
        wavelength
        [:math:`cm`]

    Z : T_VEC_IA, 
        net charge
        [:math:`-`]
    Returns
    ---------
    alpha : T_VEC_FA,
        photoionization cross section
        [:math:`cm^{2}`]

    Notes
    -------

    refer to [1]_ page 187, Eq(7.84)


    References
    -----------

    .. [1] Ivan Hubeny, Dimitri Mihalas, "Theory of Stellar Atmosphere:
        An Introduction to Astrophysical Non-equilibrium
        Quantitative Spectroscopic Analysis",
        Princeton University Press, 2015.
    """
    # ionization energy
    Eik = CST.E_Rydberg_ * (1./ni**2)

    # frequency
    v = CST.c_ / w
    # ratio of transition energy to ionization energy
    x = CST.h_ * v / Eik

    if x < 1.0:
        alpha = 0.
    else:
        alpha = 2.815E29 * Z**4  / (v**3 * ni**5) * gaunt_factor_(ni, x)
    return alpha

@OVERLOAD
def PI_cross_section_(ni : T_INT, x : T_FLOAT, Z : T_INT) -> T_FLOAT: ...
@OVERLOAD
def PI_cross_section_(ni : T_ARRAY, x : T_ARRAY, Z : T_INT) -> T_ARRAY: ...

def PI_cross_section_(ni : T_VEC_IA, x : T_VEC_IFA, Z : T_VEC_IA) -> T_VEC_FA:
    r"""
    Photoionization cross-section for hydrogen from lower level ni at wavelength w.

    Parameters
    ------------
    ni : T_VEC_IA,  
        principal quantum number of lower level
        [:math:`-`]

    x : T_VEC_IFA, 
        ratio of the transition energy to ioniozation energy
        [:math:`-`]

    Z : T_VEC_IA, 
        net charge
        [:math:`-`]
    Returns
    ---------
    alpha : T_VEC_FA, 
        photoionization cross section
        [:math:`cm^{2}`]

    Notes
    -------

    refer to [1]_ page 187, Eq(7.84)


    References
    -----------

    .. [1] Ivan Hubeny, Dimitri Mihalas, "Theory of Stellar Atmosphere:
        An Introduction to Astrophysical Non-equilibrium
        Quantitative Spectroscopic Analysis",
        Princeton University Press, 2015.
    """
    # ionization energy
    Eik = CST.E_Rydberg_ * (1./ni**2)

    # frequency
    v = x * Eik / CST.h_

    # cross section
    if x < 1.0:
        alpha = 0.
    else:
        alpha = 2.815E29 * Z**4  / (v**3 * ni**5) * gaunt_factor_(ni, x)
    
    return alpha

#-----------------------------------------------------------------------------
# spontaneous radiative recombination
#-----------------------------------------------------------------------------

def Rki_spon_rate_coe_(ni : T_VEC_IA, Te : T_VEC_IFA) -> T_VEC_FA:
    r"""
    Parameters
    ------------
    ni : T_VEC_IA,  
        principal quantum number of lower level
        [:math:`-`]

    Te : T_VEC_IFA,  
        electron temperature
        [:math:`K`]

    Returns
    ---------
    RCki : T_VEC_FA, 
        spontaneous radiative recombination rate coefficient
        [:math:`cm^{3}s^{-1}`]

    Notes
    ------
    refer to [1]_ Eq(7)

    with

    .. math:: n_{k} R_{ki}^{spon} = n_{k} n_{e} {RC}_{ki}

    References
    -----------
    .. [1] L. C. Johnson, "Approximations for collisional and
           radiative transition rates in atomic hydrogen",
           Astrophysical Journal, vol. 174, p.227, May 1972.
           1972ApJ...174..227J
    """
    kT = CST.k_ * Te
    # ionization energy
    Eik = CST.E_Rydberg_ * (1./ni**2)
    #
    r = Eik / kT

    summation  = gaunt_factor_coe_(0,ni) * Special.E1_(r)
    summation += gaunt_factor_coe_(1,ni) * Special.E2_(r)
    summation += gaunt_factor_coe_(2,ni) * Special.E3_(r)

    Ski = 5.197E-14 * r**(1.5) * _exp(r) * summation

    RCki = Ski
    return RCki

#-----------------------------------------------------------------------------
# collisional line broadening
# for hydrogen, they are
#   1. Resonance broadening
#   2. Van der Waals broadening
#   3. Linear Stark broadening
#-----------------------------------------------------------------------------
@OVERLOAD
def collisional_broadening_Res_and_Van_( ni : T_INT, nj : T_INT, nH_I_ground : T_FLOAT, Te : T_FLOAT) -> T_FLOAT: ...
@OVERLOAD
def collisional_broadening_Res_and_Van_( ni : T_INT, nj : T_INT, nH_I_ground : T_ARRAY, Te : T_ARRAY) -> T_ARRAY: ...
@OVERLOAD
def collisional_broadening_Res_and_Van_( ni : T_ARRAY, nj : T_ARRAY, nH_I_ground : T_FLOAT, Te : T_FLOAT) -> T_ARRAY: ...

def collisional_broadening_Res_and_Van_( ni : T_VEC_IA, nj : T_VEC_IA, 
           nH_I_ground : T_VEC_FA, Te : T_VEC_IFA) -> T_VEC_FA:
    r"""

    collisional broadening caused by
      1. Resonance broadening
      2. Van der Waals broadening

    Parameters
    ------------
    ni : T_VEC_IA, 
        principal quantum number of lower level
        [:math:`-`]

    nj : T_VEC_IA,  
        principal quantum number of lower level
        [:math:`-`]

    nH_I_ground : T_VEC_FA,  
        population of neutral hydrogem ground level
        [:math:`cm^{-3}`]

    Te : T_VEC_IFA, 
        electron temperature
         [:math:`K`]

    Returns
    ---------
    gamma : T_VEC_FA,
        half life time of the damping (line width of the Lorentz profile)
         [:math:`s^{-1}`]

    Notes
    ------
    refer to [1]_

    References
    -----------
    .. [1] M. C. Lortet & E. Roueff, "Broadening of Hydrogen Lines
           in a Neutral Medium", A&A, 3, 462-467, 1969
           1969A&A......3...462L
    """

    psr = _numpy.array([ 0.0, 4.94E-8, 7.93E-9, 2.75E-9, 1.29E-9, 7.14E-10 ], dtype=DT_NB_FLOAT)

    if ni == 1:
        n = nj
    else:
        n = ni

    cvdw = 1.61E-33 * ( nj**4 - ni**4 )

    V_HH = 20596. * _sqrt( Te )
    psi_w = 17. * cvdw**0.4 * V_HH**0.6

    if n <= 6 :
        psi_r = psr[ n-1 ]
        psi = ( psi_r**2.65 + psi_w**2.65 ) ** (1./2.65)
    else:
        psi = psi_w

    gamma = nH_I_ground * psi / ( 4.0 * CST.pi_ )

    return gamma

@OVERLOAD
def collisional_broadening_LinearStark_( ni : T_INT, nj : T_INT, Ne : T_FLOAT) -> T_FLOAT: ...
@OVERLOAD
def collisional_broadening_LinearStark_( ni : T_INT, nj : T_INT, Ne : T_ARRAY) -> T_ARRAY: ...
@OVERLOAD
def collisional_broadening_LinearStark_( ni : T_ARRAY, nj : T_ARRAY, Ne : T_FLOAT) -> T_ARRAY: ...

def collisional_broadening_LinearStark_( ni : T_VEC_IA, nj : T_VEC_IA, Ne : T_VEC_IFA) -> T_VEC_FA:
    r"""

    collisional broadening caused by
      1. Linear Stark broadening

    Parameters
    ------------
    ni : T_VEC_IA,  
        principal quantum number of lower level
        [:math:`-`]

    nj : T_VEC_IA,  
        principal quantum number of lower level
        [:math:`-`]

    Ne : T_VEC_IFA,  
        electron density
        [:math:`cm^{-3}`]

    Returns
    ---------
    gamma : T_VEC_FA, 
        half life time of the damping (line width of the Lorentz profile)
        [:math:`s^{-1}`]

    Notes
    ------
    refer to [1]_

    References
    -----------
    .. [1] K. Sutton, "Approximate line shapes for hydrogen",
           J. Quant. Spectrosc. Radiat. Transfer,
           Vol. 20, pp. 333-343, 1978JQSRT...20...33
    """

    if nj - ni == 1:
        a1 = 0.642
    else:
        a1 = 1.

    gamma = 0.255 * a1 * (nj*nj - ni*ni) * Ne**(2./3.)
    return gamma


#-------------------------------------------------------------------------------
# numba optimization
#-------------------------------------------------------------------------------

if CFG._IS_JIT:

    ratio_Etran_to_Eionize_ = nb_vec(**NB_VEC_KWGS) (ratio_Etran_to_Eionize_)
    
    #gaunt_factor_gingerich_cm_ = nb_vec(**NB_VEC_KWGS) (gaunt_factor_gingerich_cm_)
    #gaunt_factor_gingerich_ = nb_vec(**NB_VEC_KWGS) (gaunt_factor_gingerich_)

    gaunt_factor_coe_ = nb_vec(**NB_VEC_KWGS) (gaunt_factor_coe_)
    gaunt_factor_ = nb_vec(**NB_VEC_KWGS) (gaunt_factor_)

    absorption_oscillator_strength_ = nb_vec(**NB_VEC_KWGS) (absorption_oscillator_strength_)
    einstein_A_coefficient_ = nb_vec(**NB_VEC_KWGS) (einstein_A_coefficient_)

    CE_rate_coe_ = nb_vec(**NB_VEC_KWGS) (CE_rate_coe_)
    CI_rate_coe_ = nb_vec(**NB_VEC_KWGS) (CI_rate_coe_)

    PI_cross_section_cm_ = nb_vec(**NB_VEC_KWGS) (PI_cross_section_cm_)
    PI_cross_section_ = nb_vec(**NB_VEC_KWGS) (PI_cross_section_)
    Rki_spon_rate_coe_ = nb_vec(**NB_VEC_KWGS) (Rki_spon_rate_coe_)

    collisional_broadening_Res_and_Van_ = nb_vec(**NB_VEC_KWGS) (collisional_broadening_Res_and_Van_)
    collisional_broadening_LinearStark_ = nb_vec(**NB_VEC_KWGS) (collisional_broadening_LinearStark_)

else:
    
    #gaunt_factor_gingerich_cm_ = np_vec( gaunt_factor_gingerich_cm_, **NP_VEC_KWGS )
    gaunt_factor_coe_ = np_vec( gaunt_factor_coe_, **NP_VEC_KWGS )
    CE_rate_coe_ = np_vec( CE_rate_coe_, **NP_VEC_KWGS )
    CI_rate_coe_ = np_vec( CI_rate_coe_, **NP_VEC_KWGS )
    PI_cross_section_cm_ = np_vec( PI_cross_section_cm_, **NP_VEC_KWGS )
    PI_cross_section_ = np_vec( PI_cross_section_, **NP_VEC_KWGS )
    collisional_broadening_Res_and_Van_ = np_vec(collisional_broadening_Res_and_Van_,**NP_VEC_KWGS )
    collisional_broadening_LinearStark_ = np_vec(collisional_broadening_LinearStark_,**NP_VEC_KWGS)