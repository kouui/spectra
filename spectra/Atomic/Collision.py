
#-------------------------------------------------------------------------------
# function definition of collisional process
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

from ..ImportAll import *

import numpy as _numpy

@OVERLOAD
def interp_omega_(table : T_ARRAY, Te : T_FLOAT, Te_table: T_ARRAY ,f1 : T_FLOAT, f2 : T_FLOAT) -> T_FLOAT: ...
@OVERLOAD
def interp_omega_(table : T_ARRAY, Te : T_ARRAY, Te_table: T_ARRAY ,f1 : T_FLOAT, f2 : T_FLOAT) -> T_ARRAY: ...

def interp_omega_(table : T_ARRAY, Te : T_UNION[T_FLOAT,T_INT,T_ARRAY], Te_table: T_ARRAY ,
                 f1 : T_UNION[T_FLOAT,T_INT], f2 : T_UNION[T_FLOAT,T_INT]) -> T_UNION[T_FLOAT,T_ARRAY]:
    """given electron temperature, interpolate CE/CI coefficient

    Parameters
    ----------
    table : T_ARRAY
        a table of CE coefficient as a function of temperature for interpolation
    Te : T_UNION[T_FLOAT,T_INT,T_ARRAY]
        electron termperature, [:math:`K`]
    Te_table : T_ARRAY
        corresponding electron termperature points in _table, [:math:`K`]
    f1 : T_UNION[T_FLOAT,T_INT]
        multiply factor needed to compute CE/CI rate coefficient
    f2 : T_UNION[T_FLOAT,T_INT]
        division factor needed to compute CE/CI rate coefficient

    Returns
    -------
    T_UNION[T_FLOAT,T_ARRAY]
        CE/CI coefficient
    """
    omega = _numpy.interp( Te, Te_table[:], table[:] )  * f1 / f2
    return omega

def Cij_to_Cji_(Cij : T_VEC_FA,  nj_by_ni_LTE : T_VEC_FA) -> T_VEC_FA:
    """Given the LTE population ration and collisional coefficient Cij, 
       calculate collisional coefficient Cji

    Parameters
    ----------
    Cij : T_VEC_FA
        collisional coefficient Cij, [:math:`cm^{-3} s^{-1}`]
    nj_by_ni_LTE : T_VEC_FA
        population ratio in LTE, nj/ni

    Returns
    -------
    T_VEC_FA
        collisional coefficient Cji, [:math:`cm^{-3} s^{-1}`]
    """
    Cji = Cij / nj_by_ni_LTE

    return Cji

@OVERLOAD
def CE_rate_coe_(omega : T_FLOAT, Te : T_FLOAT, gi : T_INT, dEij : T_FLOAT) -> T_FLOAT: ...
@OVERLOAD
def CE_rate_coe_(omega : T_FLOAT, Te : T_ARRAY, gi : T_INT, dEij : T_FLOAT) -> T_ARRAY: ...
@OVERLOAD
def CE_rate_coe_(omega : T_ARRAY, Te : T_FLOAT, gi : T_ARRAY, dEij : T_ARRAY) -> T_ARRAY: ...

def CE_rate_coe_(omega : T_VEC_FA, Te : T_VEC_IFA, gi : T_VEC_IA, dEij : T_VEC_FA) -> T_VEC_FA:
    """[summary]

    Parameters
    ----------
    omega : T_VEC_FA
        the collisional strength we interpolate from data
    Te : T_VEC_IFA
        electron temperature, [:math:`K`]
    gi : T_VEC_IA
        statistical weight of lower level of the transition, [:math:`-`]
    dEij : T_VEC_FA
        excitation energy, [:math:`erg`]

    Returns
    -------
    T_VEC_FA
        collisional excitation rate coefficient, [:math:`cm^{3} s^{-1}`]
    """
    #Cst.pi_ * Cst.a0_**2 * (8*Cst.k_/Cst.pi_/Cst.me_)**0.5 * Cst.E_Rydberg_ / Cst.k_
    #=8.629140599482176e-06 [in CGS unit]

    kT = CST.k_ * Te
    CEij = (8.63E-06 * omega) / (gi * Te**0.5) * _numpy.exp( - dEij / kT )

    return CEij

@OVERLOAD
def CI_rate_coe_(omega : T_FLOAT, Te : T_FLOAT, dEik : T_FLOAT) -> T_FLOAT: ...
@OVERLOAD
def CI_rate_coe_(omega : T_FLOAT, Te : T_ARRAY, dEik : T_FLOAT) -> T_ARRAY: ...
@OVERLOAD
def CI_rate_coe_(omega : T_ARRAY, Te : T_FLOAT, dEik : T_ARRAY) -> T_ARRAY: ...

def CI_rate_coe_(omega : T_VEC_FA, Te : T_VEC_FA, dEik : T_VEC_FA) -> T_VEC_FA:
    """[summary]

    Parameters
    ----------
    omega : T_VEC_FA
        the collisional strength (?) we interpolate from data
    Te : T_VEC_FA
        electron temperature, [:math:`K`]
    dEik : T_VEC_FA
        ionization energy, [:math:`erg`]

    Returns
    -------
    T_VEC_FA
        collisional ionization rate coefficient, [:math:`cm^{3} s^{-1}`]
    """
    kT = CST.k_ * Te
    CIij = omega * Te**0.5 * _numpy.exp( - dEik / kT )

    return CIij

#-----------------------------------------------------------------------------
# numba optimization
#-----------------------------------------------------------------------------

if CFG._IS_JIT:
    interp_omega_     = nb_njit( **NB_NJIT_KWGS ) ( interp_omega_ )
    Cij_to_Cji_       = nb_vec( **NB_VEC_KWGS ) ( Cij_to_Cji_ )
    CE_rate_coe_      = nb_vec( **NB_VEC_KWGS ) ( CE_rate_coe_ )
    CI_rate_coe_      = nb_vec( **NB_VEC_KWGS ) ( CI_rate_coe_ )
    