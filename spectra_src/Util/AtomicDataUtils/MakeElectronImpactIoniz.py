#-------------------------------------------------------------------------------
# definition of functions for make electron impact ionoization atomic data
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/07/13   u.k.   
#-------------------------------------------------------------------------------

from spectra_src.ImportAll import *
import os

import numpy as _numpy
from numpy import sqrt as _sqrt
from numpy import exp as _exp

import warnings as _warnings

def _is_skip_(line: T_STR) -> T_BOOL:

    # comment lines
    if line[0] == '#':
        return True

    # empty lines
    if len( line.strip() ) == 0:
        return True

    return False

#-------------------------------------------------------------------------------
# A practical fit formula for ionization rate coefficients of atoms and ions 
# by electron impact: Z = 1-28.
# G. S. Voronov (1997, Atomic Data and Nuclear Data Tables, 65, 1). 
# ground.Z_1_28.txt
#-------------------------------------------------------------------------------

def _v1_read_data_() -> T_DICT[T_INT, T_DICT[T_INT,T_ARRAY] ]:

    file = os.path.join( CFG._ROOT_DIR, "data/atomic_raw/ionization/electron_impact_ionization/ground.Z_1_28.txt" )

    data : T_DICT[T_INT, T_DICT[T_INT,T_ARRAY] ] = {}

    with open(file, 'r') as f:

        for line in f:
            #print( line )

            if _is_skip_(line):
                continue

            words = [w.strip() for w in line[5:].strip().split()]
            #print(words)
            nZ = int( words[0] )
            nE = int( words[1] )
            ioniz_stage = nZ - nE + 1 # neutral --> 1
            params = _numpy.array( [float(v) for v in words[2:9]], dtype=DT_NB_FLOAT )

            try:
                _ = data[nZ]
            except KeyError:
                data[nZ] = {}

            data[nZ][ioniz_stage] = params

    return data

_V1_DATA = _v1_read_data_()

#def _v1_fit_func_(Te: T_VEC_FA, params: T_ARRAY) -> T_VEC_FA:
#
#    p_dE    = params[0]
#    p_P     = params[1]
#    p_A     = params[2]
#    p_X     = params[3]
#    p_K     = params[4]
#    p_T_min = params[5]
#    p_T_max = params[6]
def _v1_fit_func_(Te: T_VEC_FA, 
        p_dE : T_FLOAT, p_P : T_FLOAT, p_A : T_FLOAT, p_X : T_FLOAT,
        p_K : T_FLOAT, p_T_min : T_FLOAT, p_T_max : T_FLOAT) -> T_VEC_FA:


    Te_eV = (Te * CST.k_) / CST.eV2erg_
#    if Te_eV < p_T_min:
#        _warnings.warn(f"{Te}[K] < {p_T_min}[eV] (min)")
#    if Te_eV > p_T_max * 1.E3:
#        _warnings.warn(f"{Te}[K] > {p_T_max}[keV] (max)")

    u = p_dE / Te_eV # U is dimemsionless

    # [cm^{3} s*{-1}]
    rate_coe = p_A * ( 1. + p_P * _sqrt(u) ) / (p_X + u) * u**p_K * _exp( -u )

    return rate_coe

_v1_fit_func_ = np_vec( _v1_fit_func_, **NP_VEC_KWGS)

def v1_rate_coe_(Te : T_VEC_IFA, nZ : T_INT, ioniz_stage: T_INT) -> T_VEC_FA:

    params = _V1_DATA[nZ][ioniz_stage]

#    rate_coe = _numpy.empty_like( Te )
#
#    for k in range(rate_coe.size):
#        rate_coe[k] = _v1_fit_func_(Te[k], params)
    rate_coe = _v1_fit_func_(Te, *params)

    return rate_coe


#-------------------------------------------------------------------------------
# quantities conversion
# collision rate coefficient --> collisional strength 
#-------------------------------------------------------------------------------

def rate_coe_to_collision_strength_(dEik : T_FLOAT, Te : T_VEC_FA, rate_coe : T_VEC_FA) -> T_VEC_FA:

    kT = CST.k_ * Te
    omega = rate_coe * Te**(-0.5) * _numpy.exp( dEik / kT )
    
    return omega







    
            
                
            