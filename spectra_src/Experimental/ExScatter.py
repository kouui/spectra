
#-------------------------------------------------------------------------------
# compute 
# - continuum absorption coefficient
# - scattering coefficient 
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/08/10   u.k.   
#-------------------------------------------------------------------------------

from ..ImportAll import *
from ..Atomic import ContinuumOpacity as _ContinuumOpacity
from ..Atomic import Hydrogen as _Hydrogen
from ..Atomic import LTELib as _LTELib

from numpy import exp as _exp

def HI_bf_cross_sec_(ns : T_VEC_IA, nHI_pop : T_VEC_FA, Te : T_FLOAT):
    r"""
    compute H bound-free total cross section in NLTE case
    """
    pass



#def func_(ns : T_VEC_IA, nH_pop : T_VEC_FA, nH : T_FLOAT, ne : T_FLOAT, Te : T_FLOAT):

def HI_bf_emissivity_LTE_(wave : T_FLOAT, nHI_pop_LTE: T_VEC_FA, Te: T_FLOAT):

    n_HI_Level = nHI_pop_LTE.size
    
    eta = 0.0
    kT = CST.k_ * Te
    factor = (1-_exp(-CST.h_*CST.c_/(wave*kT)))
    for k in range(n_HI_Level):
        ni = k
        alpha = _Hydrogen.PI_cross_section_cm_(ni, wave, 1)
        eta += nHI_pop_LTE[k] * alpha * factor * _LTELib.planck_cm_(wave, Te)

    return eta

def background_opacity_(wave : T_FLOAT, nH_pop_LTE : T_VEC_FA, nH_pop : T_ARRAY, nH : T_FLOAT, ne : T_FLOAT, Te : T_FLOAT):

    n_Hp = nH * nH_pop[-1]
    n_HI = nH - n_Hp
    # alpha : line extinction coefficient per particle [cm^{-1}]
    # sigma : 

    # Thomson scattering
    alpha_thomson = _ContinuumOpacity.thomson_scattering_(ne)
    # Rayleigh scattering
    alpha_rayleigh = _ContinuumOpacity.HI_rayleigh_cross_sec_(wave) * n_HI
    # H- LTE assumption
    alpha_Hminus = _ContinuumOpacity.Hminus_cross_sec_(Te, wave, ne) * n_HI
    # H2+ LTE assumption
    alpha_H2plus = _ContinuumOpacity.H2p_cross_sec_(Te, wave, nH_pop[-1]) * n_HI
    # HI bf LTE assumption
    #alpha_HIbf = _ContinuumOpacity.HI_bf_LTE_cross_sec_(Te, wave) * nH
    # HII ff LTE assumption)
    alpha_HIff = _ContinuumOpacity.Hp_ff_cross_sec_(Te, wave, ne, 1) * n_Hp

    # eta :  emissivity
    # HI bf LTE assumption
    #eta_HI_bf_LTE =   HI_bf_emissivity_LTE_(wave, nH_pop_LTE[:-1], Te) * nH
    # planck function
    plk = _LTELib.planck_cm_(wave, Te)


    # coherent scattering absorption coefficient
    alpha_coherent = alpha_thomson + alpha_rayleigh
    # local process absorption coefficient
    alpha_local = alpha_H2plus + alpha_HIff + alpha_Hminus
    # total background absorption coefficient
    alpha_background = alpha_coherent + alpha_local

    # background source function
    #src_background = () / sigma_local
    s_background = plk

    # total emissivity (without scattering)
    eta_local = s_background * alpha_local

    epsilon_background = alpha_coherent / alpha_background # [-]

    return alpha_background, eta_local, epsilon_background
    







