
#-------------------------------------------------------------------------------
# function definition of LTE process
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------
from ..ImportAll import *

import numpy as _numpy

from debtcollector import removals as _removals  # type: ignore


def boltzmann_distribution_(gi : T_VEC_IA, gj : T_VEC_IA, 
                            Eji : T_VEC_FA, Te : T_VEC_IFA) -> T_VEC_FA:
    """calculate the population ratio between upper level j and lower level i under LTE.

    Parameters
    ----------
    gi : T_VEC_IA
        statistical weight of lower level i, [-]
    gj : T_VEC_IA
        statistical weight of upper level j, [-]
    Eji : T_VEC_FA
        the gap of level energy between upper level j and lower level i
        :math:`E_{ji}=E_j-E_i, \quad [erg]`
    Te : T_VEC_IFA
        electron temperature
        [:math:`K`]

    Returns
    -------
    T_VEC_FA
        population ratio of upper level j and lower level i.
        :math:`rt = n_j / n_i`, [-]
    """
    kT = CST.k_ * Te
    ratio = gj / gi * _numpy.exp( -Eji / kT )

    return ratio

def saha_distribution_(gi : T_VEC_IA, gk : T_VEC_IA, 
                       chi : T_VEC_FA, Ne : T_VEC_IFA, 
                       Te : T_VEC_IFA) -> T_VEC_FA:
    """calculate the population ratio between the ground states of two subsequent ionization stage under LTE.

    Parameters
    ----------
    gi : T_VEC_IA
        statistical weight of the ground state in ionization stage I, [-]
        [-]
    gk : T_VEC_IA
        statistical weight of the ground state in ionization stage I+1 
        [-]
    chi : T_VEC_FA
        ionization energy from ground state in ionization stage I
        to ground state in ionization stage I, 
        [:math:`erg`]
    Ne : T_VEC_IFA
        electron density, 
        [:math:`cm^{-3}`]
    Te : T_VEC_IFA
        electron temperature, 
        [:math:`K`]

    Returns
    -------
    T_VEC_FA
        :math:`n_k / n_i`, :math:`n_k` the population of the ground state of ionization stage I+1
        and :math:`n_i` the population of the ground state of ionization stage I. 
        [-]
    """
    kT = CST.k_ * Te
    ratio = 1. / Ne * CST.saha_ * Te**(1.5) * gk / gi * _numpy.exp( -chi / kT)

    return ratio

_removals.remove
def LTE_ratio_Line_(g : T_ARRAY, idxI : T_ARRAY, idxJ : T_ARRAY, 
                    w0 : T_ARRAY, Te : T_UNION[T_FLOAT,T_INT]) -> T_ARRAY:
    """Compute LTE population ratio nj/ni for each line transition

    Parameters
    ----------
    g : T_ARRAY
        statistical weight, 
        [-]
    idxI : T_ARRAY
        lower level index for each line transition
        [-]
    idxJ : T_ARRAY
        upper level index for each line transition
        [-]
    w0 : T_ARRAY
        wavelength of each line transition, 
        [:math:`cm`]
    Te : T_UNION[T_FLOAT,T_INT]
        electron temperature, 
        [:math:`K`]

    Returns
    -------
    T_ARRAY
        population ratio. 
        [-]
    """
    nLine = w0.size
    hc = CST.h_ * CST.c_
    nRatio = _numpy.ones(nLine, dtype=T_FLOAT)
    for k in range(nLine):
        gi = g[ idxI[k] ]
        gj = g[ idxJ[k] ]
        chi = hc / w0[k]
        nRatio[k] = boltzmann_distribution_(gi, gj, chi, Te)

    return nRatio

_removals.remove
def LTE_ratio_Cont_(g : T_ARRAY, idxI : T_ARRAY, idxJ : T_ARRAY,
                     w0 : T_ARRAY, Te : T_UNION[T_FLOAT,T_INT], 
                     Ne : T_UNION[T_FLOAT,T_INT]) -> T_ARRAY:
    """Compute LTE population ratio nj/ni for each continuum transition

    Parameters
    ----------
    g : T_ARRAY
        statistical weight, 
        [-]
    idxI : T_ARRAY
        lower level index for each line transition
        [-]
    idxJ : T_ARRAY
        upper level index for each line transition
        [-]
    w0 : T_ARRAY
        wavelength of each line transition, 
        [:math:`cm`]
    Te : T_UNION[T_FLOAT,T_INT]
        electron temperature, 
        [:math:`K`]
    Ne : T_UNION[T_FLOAT,T_INT]
        electron density, [:math:`cm^{-3}`]

    Returns
    -------
    T_ARRAY
        population ratio. 
        [-]
    """
    nLine = w0.size
    hc = CST.h_ * CST.c_
    nRatio = _numpy.ones(nLine, dtype=T_FLOAT)
    for k in range(nLine):
        gi = g[ idxI[k] ]
        gk = g[ idxJ[k] ]
        chi = hc / w0[k]
        nRatio[k] = saha_distribution_(gi, gk, chi, Ne, Te)

    return nRatio

_removals.remove
def LTE_ratio_(erg : T_ARRAY, g : T_ARRAY, stage : T_ARRAY, 
               Te : T_UNION[T_FLOAT,T_INT], Ne : T_UNION[T_FLOAT,T_INT]) -> T_ARRAY:
    """Compute normalized LTE population

    Parameters
    ----------
    erg : T_ARRAY
        level energy relative to 1st level, 
        [:math:`erg`]
    g : T_ARRAY
        statistical weight, 
        [-]
    stage : T_ARRAY
        ionization stage, 
        [-]
    Te : T_UNION[T_FLOAT,T_INT]
        electron temperature, 
        [:math:`K`]
    Ne : T_UNION[T_FLOAT,T_INT]
        electron density, 
        [:math:`cm^{-3}`]

    Returns
    -------
    T_ARRAY
        normalized population ratio. [-]
    """
    nLevel = erg.size
    nRatio = _numpy.empty(nLevel, T_FLOAT)
    nRatio[0] = 1.
    for i in range(1, nLevel):
        gj, gi = g[i], g[i-1]
        if stage[i] - stage[i-1] == 0:
            nRatio[i] = nRatio[i-1] * boltzmann_distribution_(gi, gj, erg[i]-erg[i-1], Te)
        elif stage[i] - stage[i-1] == 1:
            nRatio[i] = nRatio[i-1] * saha_distribution_(gi, gj, erg[i]-erg[i-1], Ne, Te)

        nRatio[:] /= nRatio[:].sum()

    return nRatio

_removals.remove
def einsteinA_to_einsteinBs_hz_(Aji : T_UNION[T_ARRAY, T_FLOAT], 
    f0 : T_UNION[T_ARRAY, T_FLOAT], gi : T_UNION[T_ARRAY, T_INT], 
    gj : T_UNION[T_ARRAY, T_INT]) ->  T_UNION[ T_TUPLE[T_ARRAY,T_ARRAY], T_TUPLE[T_FLOAT, T_FLOAT] ]:
    
    """

    given Einstein A coefficient Aij,
    calculate Einstein B coefficients Bij and Bji.

    Parameters
    ----------

    Aji : T_UNION[T_ARRAY, T_FLOAT]
        Einstein A coefficient Aji, 
        [:math:`s^{-1}`]
    f0 : T_UNION[T_ARRAY, T_FLOAT]
        central frequency of corresponding line transition, 
        [:math:`Hz`]
    gi : T_UNION[T_ARRAY, T_INT]
        statistical weight of lower level, 
        [-]
    gj : T_UNION[T_ARRAY, T_INT]
        statistical weight of upper level, 
        [-]

    Returns
    -------

    Bji : T_UNION[T_ARRAY, T_FLOAT]
        Einstein B coefficient Bji, 
        [:math:`s^{-1}/(erg/cm^{2}/Sr/Hz/s)`]
    Bij : T_UNION[T_ARRAY, T_FLOAT]
        Einstein B coefficient Bji, 
        [:math:`s^{-1}/(erg/cm^{2}/Sr/Hz/s)`]

    Notes
    -----

    Einstein Relations for bound-bound transitions [1]_.

    .. math:: B_{ji} = A_{ji} / \frac{2h\nu^{3}}{c^{2}}
    .. math:: B_{ij} = B_{ji} \frac{g_{j}}{g_{i}}

    References
    ----------

    .. [1] Ivan Hubeny, Dimitri Mihalas, "Theory of Stellar Atmosphere:
        An Introduction to Astrophysical Non-equilibrium
        Quantitative Spectroscopic Analysis",
        Princeton University Press, pp. 117, 2015.
    """

    factor = 2. * CST.h_ * f0*f0*f0 / CST.c_*CST.c_
    Bji = Aji / factor
    Bij = Bji * T_FLOAT(gj) / T_FLOAT(gi)

    return Bji, Bij


def einsteinA_to_einsteinBs_cm_(Aji : T_UNION[T_ARRAY, T_FLOAT], 
       w0 : T_UNION[T_ARRAY, T_FLOAT], gi : T_UNION[T_ARRAY, T_INT],
       gj : T_UNION[T_ARRAY, T_INT])->  T_UNION[ T_TUPLE[T_ARRAY,T_ARRAY], T_TUPLE[T_FLOAT, T_FLOAT] ]:
    """

    given Einstein A coefficient Aij,
    calculate Einstein B coefficients Bij and Bji.

    Parameters
    ----------

    Aji : T_UNION[T_ARRAY, T_FLOAT]
        Einstein A coefficient Aji, 
        [:math:`s^{-1}`]
    w0 : T_UNION[T_ARRAY, T_FLOAT]
        central wavelength of corresponding line transition, 
        [:math:`cm`]
    gi : T_UNION[T_ARRAY, T_INT]
        statistical weight of lower level, 
        [-]
    gj : T_UNION[T_ARRAY, T_INT]
        statistical weight of upper level, 
        [-]

    Returns
    -------

    Bji : T_UNION[T_ARRAY, T_FLOAT]
        Einstein B coefficient Bji, 
        [:math:`s^{-1}/(erg/cm^{2}/Sr/Hz/s)`]
    Bij : T_UNION[T_ARRAY, T_FLOAT]
        Einstein B coefficient Bji, 
        [:math:`s^{-1}/(erg/cm^{2}/Sr/Hz/s)`]

    Notes
    -----

    Einstein Relations for bound-bound transitions.
    Since we are in wavelength unit, the relation between
    [:math:`A_{ji}`] and [:math:`B_{ji}`] should be evaluated using
    Planck function also in wavelength unit

    .. math:: B_{ji} = A_{ji} / (2 h c^{2} \lambda^{5} )
    .. math:: B_{ij} = B_{ji} \frac{g_{j}}{g_{i}}

    """
    factor = 2. * CST.h_ * CST.c_*CST.c_ / w0**5
    Bji = Aji / factor
    Bij = Bji * gj / gi

    return Bji, Bij

_removals.remove
def Aji_to_Bji_cm_(Aji : T_UNION[T_ARRAY, T_FLOAT], 
       w0 : T_UNION[T_ARRAY, T_FLOAT]) -> T_UNION[T_ARRAY, T_FLOAT]:
    
    """Given Einstein Aji coefficient, compute Einstein Bji coefficient 

    Parameters
    ----------
    Aji : T_UNION[T_ARRAY, T_FLOAT]
        Einstein A coefficient Aji,
        [:math:`s^{-1}`]
    
    w0 : T_UNION[T_ARRAY, T_FLOAT]
        wavelength
        [:math:`cm`]

    Returns
    -------
    T_UNION[T_ARRAY, T_FLOAT]
        Einstein B coefficient Bji,
        [?]
    """
    factor = 2. * CST.h_ * CST.c_*CST.c_ / w0**5
    Bji = Aji / factor

    return Bji


_removals.remove
def Bji_to_Bij_(Bji : T_UNION[T_ARRAY, T_FLOAT], 
                gi : T_UNION[T_ARRAY, T_INT], 
                gj : T_UNION[T_ARRAY, T_INT]) -> T_UNION[T_ARRAY, T_FLOAT]:
    """[summary]

    Parameters
    ----------
    Bji : T_UNION[T_ARRAY, T_FLOAT]
        Einstein B coefficient Bji,
        [?]
    gi : T_UNION[T_ARRAY, T_INT]
        statistical weight of the lower level
    gj : T_UNION[T_ARRAY, T_INT]
        statistical weight of the upper level

    Returns
    -------
    T_UNION[T_ARRAY, T_FLOAT]
        Einstein B coefficient Bij,
        [?]
    """
    Bij = Bji * gj / gi

    return Bij

_removals.remove
def planck_hz_(F : T_UNION[T_ARRAY, T_FLOAT, T_INT],
               T : T_UNION[T_ARRAY, T_FLOAT, T_INT]) -> T_UNION[T_ARRAY, T_FLOAT, T_INT]:
    """
    given frequency and temperature,
    calculate the frequency based planck function.


    Parameters
    ----------

    F : T_UNION[T_ARRAY, T_FLOAT, T_INT]
        frequency, 
        [:math:`Hz`]
    T : T_UNION[T_ARRAY, T_FLOAT, T_INT]
        temperature, 
        [:math:`K`]

    Returns
    -------

    intensity : T_UNION[T_ARRAY, T_FLOAT, T_INT]
        frequency based intensity, 
        [:math:`erg/cm^2/Sr/Hz/s`]

    Notes
    -----

    The Planck function [1]_.

    .. math:: B_{\nu}(T) = \frac{2h\nu^3}{c^2} \frac{1}{e^{h\nu/kT}-1} \quad [erg/cm^2/Sr/Hz/s]

    References
    ----------

    .. [1] Ivan Hubeny, Dimitri Mihalas, "Theory of Stellar Atmosphere:
        An Introduction to Astrophysical Non-equilibrium
        Quantitative Spectroscopic Analysis",
        Princeton University Press, pp. 102, 2015.
    """
    F_T = F / T
    if F_T > 1.04183e+13:                                        # ignore exponential part, exponential > 500
        intensity = 0
    #elif F_T > 2.08366e+12:                                      # Wien approximation, exponential > 100
    #    intensity = 2.0*Cst.h_*F*F*F/Cst.c_/Cst.c_ * np.exp( -Cst.h_*F_T/(Cst.k_) )
    #elif F_T < 2.08366e+08:                                      # Rayleighâ€“Jeans approximation, exponential < 0.01
    #    intensity = 2.0*F*F/Cst.c_/Cst.c_ * Cst.k_ * T
    else:                                                        # normal formula
        intensity = 2.0*CST.h_*F*F*F/CST.c_/CST.c_ / ( _numpy.exp( CST.h_*F_T/(CST.k_) ) - 1. )

    return intensity

def planck_cm_(W : T_UNION[T_ARRAY, T_FLOAT, T_INT],
               T : T_UNION[T_ARRAY, T_FLOAT, T_INT]) -> T_UNION[T_ARRAY, T_FLOAT, T_INT]:
    """
    given wavelength and temperature,
    calculate the wavelength based planck function.


    Parameters
    ----------

    W : T_UNION[T_ARRAY, T_FLOAT, T_INT]
        wavelength, 
        [:math:`cm`]
    T : T_UNION[T_ARRAY, T_FLOAT, T_INT]
        temperature, 
        [:math:`K`]

    Returns
    -------

    intensity : T_UNION[T_ARRAY, T_FLOAT, T_INT]
        wavelength based intensity, 
        [:math:`erg/cm^2/Sr/cm/s`]

    Notes
    -----

    The Planck function [1]_.

    .. math:: B_{\lambda}(T) = \frac{2hc^2}{\lambda^5} \frac{1}{e^{hc/\lambda kT}-1} \quad [erg/cm^2/Sr/cm/s]

    References
    ----------

    .. [1] Robert J. Rutten, "Introduction to Astrophysical Radiative Transfer", pp. 43, 2015.
    """
    WT = W*T
    if WT < 2.87755e-03:                                         # ignore exponential part, exponential > 500
        intensity = 0
    #elif WT < 1.43878e-02:                                       # Wien approximation, exponential > 100
    #    intensity = 2.0*Cst.h_*Cst.c_*Cst.c_ /(W)**5 * np.exp( -Cst.h_*Cst.c_/(Cst.k_*WT) )
    #elif WT > 1.43878e+02:                                       # Rayleighâ€“Jeans approximation, exponential < 0.01
    #    intensity = 2.0*Cst.c_*Cst.k_*T / (W)**4
    else:                                                        # normal formula
        intensity = 2.0*CST.h_*CST.c_*CST.c_ /(W)**5 /  ( _numpy.exp( CST.h_*CST.c_/(CST.k_*WT) ) - 1.0 )

    return intensity


#-----------------------------------------------------------------------------
# partition function
#-----------------------------------------------------------------------------

_UFUNC_COEFFICIENT_TABLE : T_DICT[ T_STR, T_TUPLE[T_FLOAT, ...] ] = {
    'H_II'   : (0., 0., 0., 0., 0.),
    'H_I'    : (0.30103, 0., 0., 0., 0.),
    'He_I'   : (0., 0., 0., 0., 0.),
    'He_II'  : (0.30103, 0., 0., 0., 0.),
    'Li_I'   : (0.31804, -0.20616, 0.91456, -1.66121, 1.04195),
    'Be_I'   : (0.00801, -0.17135, 0.62921, -0.58945, 0.),
    'Be_II'  : (0.30389, -0.00819, 0., 0., 0.),
    'B_i'    : (0.78028, -0.01622, 0., 0., 0.),
    'C_I'    : (0.96752, -0.09452, 0.08055, 0., 0.),
    'C_II'   : (0.77239, -0.02540, 0., 0., 0.),
    'N_I'    : (0.60683, -0.08674, 0.30565, -0.28114, 0.),
    'N_II'   : (0.94968, -0.06463, -0.1291, 0., 0.),
    'O_I'    : (0.05033, -0.05703, 0., 0., 0.),
    'O_II'   : (0.60405, -0.03025, 0.04525, 0., 0.),
    'F_I'    : (0.76284, -0.03582, -0.05619, 0., 0.),
    'Ne_I'   : (0., 0., 0., 0., 0.),  # ufunc = 1,
    'Ne_II'  : (0.74847, -0.06562, -0.07088, 0., 0.),
    'Na_I'   : (0.30955, -0.17778,  1.10594, -2.42847, 1.70721),
    'Na_II'  : (0., 0., 0., 0., 0.),  # ufunc = 1,
    'Na_III' : (_numpy.log10(6), 0., 0., 0., 0.),  # ufunc1 = 6.0   # Allen, 197,
    'Mg_I'   : (0.00556, -0.12840,  0.81506, -1.79635, 1.26292),
    'Mg_II'  : (0.30257, -0.00451,  0.,  0., 0.),
    'Mg_III' : (0., 0., 0., 0., 0.),  # ufunc = 1,
    'Al_I'   : (0.76786, -0.05207,  0.14713,  -0.21376, 0.),
    'Al_II'  : (0.00334, -0.00995,  0.,  0., 0.),
    'Si_I'   : (0.97896, -0.19208,  0.04753,  0., 0.),
    'Si_II'  : (0.75647, -0.05490,  -0.10126,  0., 0.),
    'P_I'    : (0.64618, -0.31132,  0.68633,  -0.47505, 0.),
    'P_II'   : (0.93588, -0.18848,  0.08921,  -0.22447, 0.),
    'S_I'    : (0.95254, -0.15166,  0.02340,  0., 0.),
    'S_II'   : (0.61971, -0.17465,  0.48283,  -0.39157, 0.),
    'Cl_I'   : (0.74465, -0.07389,  -0.06965,  0., 0.),
    'Cl_II'  : (0.92728, -0.15913,  -0.01983,  0., 0.),
    'K_I'    : (0.34419, -0.48157,  1.92563,  -3.17826, 1.83211),
    'Ca_I'   : (0.07460, -0.75759, 2.58494, -3.53170, -1.65240),
    'Ca_II'  : (0.34383, -0.41472, 1.01550, 0.31930, 0.),
    'Ca_III' : (0., 0., 0., 0., 0.),  # ufunc = 1,
    'Sc_I'   : (1.08209, -0.77814, 1.78504, -1.39179, 0.),
    'Sc_II'  : (1.35894, -0.51812, 0.15634, 0., 0.),
    'Ti_I'   : (1.47343, -0.97220, 1.47986, -0.93275, 0.),
    'Ti_II'  : (1.74561, -0.51230, 0.27621, 0., 0.),
    'Ti_III' : (0., 0., 0., 0., 0.),  # ufunc = 1,
    'V_I'    : (1.68359, -0.82055, 0.92361, -0.78342, 0.),
    'V_II'   : (1.64112, -0.74045, 0.49148, 0., 0.),
    'V_III'  : (_numpy.log10(28), 0., 0., 0., 0.),  # ufunc = 28,
    'Cr_I'   : (1.02332, -1.02540, 2.02181, -1.32723, 0.),
    'Cr_II'  : (0.85381, -0.71166, 2.18621, -0.97590, -2.72893),
    'Cr_III' : (_numpy.log10(25), 0., 0., 0., 0.),  # ufunc = 25,
    'Mn_I'   : (0.80810, -0.39108, 1.74756, -3.13517, 1.93514),
    'Mn_II'  : (0.88861, -0.36398, 1.39674, -1.86424, -2.32389),
    'Mn_III' : (_numpy.log10(6), 0., 0., 0., 0.),  # ufunc = 6,
    'Fe_I'   : (1.44701, -0.67040, 1.01267, -0.81428, 0.),
    'Fe_II'  : (1.63506, -0.47118, 0.57918, -0.12293, 0.),
    'Fe_III' : (_numpy.log10(25), 0., 0., 0., 0.),  # ufunc = 25,
    'Co_I'   : (1.52929, -0.71430, 0.37210, -0.23278, 0.),
    'Ni_I'   : (1.49063, -0.33662, 0.08553, -0.19277, 0.),
    'Ni_II'  : (1.03800, -0.69572, 0.53893, 0.28861, 0.),
    'Ni_III' : (_numpy.log10(21), 0., 0., 0., 0.),  # ufunc = 21,
    'Ba_I'   : (4.83034, -10.8244, 10.0364, -4.34979, 0.719568),
    'Ba_II'  : (2.54797, -6.38871, 7.98813, -4.38907, 0.862666),
    'Ba_iii' : (0., 0., 0., 0., 0.)  # ufunc = 1,
}

def Ufunc_(elm : T_STR, T : T_UNION[T_ARRAY,T_FLOAT,T_INT]) -> T_UNION[T_ARRAY,T_FLOAT]:
    """
    partition function of neutral hydrogen


    Parameters
    ----------
    elm : T_STR
        element & ionization stage as ca_i, fe_ii, etc.
    T   : T_UNION[T_ARRAY,T_FLOAT,T_INT]
        [:math:`K`]

    Returns
    --------
    T_UNION[T_ARRAY, T_FLOAT]
        partition function of neutral hydrogen
        [-]

    Notes
    -----
    Wanring : this polynomial fitting formula(table) is only avaible for low temperature (:math:`<10000 K`)

    from [1]_, [2]_.
    .. math:: log(u) = c_0 + c_1 \cdot log(\theta) + c_2 \cdot log(\theta)^2 + c_3 \cdot log(\theta)^3 + c_4 \cdot log(\theta)^4
    .. math:: log = log_{10}
    .. math:: \theta = 5040 / T

    modification history:
    - 2006.5.23     k.i.
    - 2015.7.5      k.i.	'h_ii'
    - 2019.11.26    k.i.	'Ba' from Gary 2009 (use poly_ufunc.pro)
    - 2020.2.13     k.i.	from IDL ufunc_gray.pro

    References
    ----------

    .. [1] Gray 1992, "The observation and analysis of stellar photospheres", app.D
    .. [2] 2009  table
    """

    c : T_TUPLE[T_FLOAT, ...] = _UFUNC_COEFFICIENT_TABLE[elm]
    
    th = 5040. / T
    s : T_UNION[T_ARRAY, T_FLOAT] = c[0]
    for i in range(1,5):
        s += c[i] * (_numpy.log10(th))**i
    ufunc1 = 10.**s

    return ufunc1


#-----------------------------------------------------------------------------
# numba optimization
#-----------------------------------------------------------------------------

if CFG._IS_JIT:

    boltzmann_distribution_ = nb_vec( **NB_VEC_KWGS ) ( boltzmann_distribution_ )
    saha_distribution_      = nb_vec( **NB_VEC_KWGS ) ( saha_distribution_ )
    planck_cm_              = nb_vec( **NB_VEC_KWGS ) ( planck_cm_ )