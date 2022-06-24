#-------------------------------------------------------------------------------
# definition of functions for make photoionization atomic data
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

def _is_skip_(line: T_STR) -> T_BOOL:

    # comment lines
    if line[0] == '#':
        return True

    # empty lines
    if len( line.strip() ) == 0:
        return True

    return False

#-------------------------------------------------------------------------------
# Atomic data for astrophysics. II. New analytic fits for photoionization
# cross sections of atoms and ions.
# Verner, D. A., Ferland, G. J., Korista, K. T., & Yakovlev, D. G. 1996, ApJ, in press.
# ground.Z_1_26.txt
#-------------------------------------------------------------------------------

def _v1_read_data_() -> T_DICT[T_INT, T_DICT[T_INT,T_ARRAY] ]:

    file = os.path.join( CFG._ROOT_DIR, "data/atomic_raw/ionization/photoionization/ground.Z_1_26.txt" )

    data : T_DICT[T_INT, T_DICT[T_INT,T_ARRAY] ] = {}

    with open(file, 'r') as f:

        for line in f:

            if _is_skip_(line):
                continue

            words = [w.strip() for w in line.split()]
            nZ = int( words[0] )
            nE = int( words[1] )
            ioniz_stage = nZ - nE + 1 # neutral --> 1
            params = _numpy.array( [float(v) for v in words[2:11]], dtype=DT_NB_FLOAT )

            try:
                _ = data[nZ]
            except KeyError:
                data[nZ] = {}

            data[nZ][ioniz_stage] = params

    return data

_V1_DATA = _v1_read_data_()

#def _v1_fit_func_(wave: T_FLOAT, params: T_ARRAY) -> T_VEC_FA:
#
#    p_E_th    = params[0] # eV
#    p_E_max   = params[1] # eV
#    p_E_0     = params[2] # eV
#    p_sigma0  = params[3] # Mb
#    p_y_a     = params[4]
#    p_P       = params[5]
#    p_y_w     = params[6]
#    p_y_0     = params[7]
#    p_y_1     = params[8]

def _v1_fit_func_(wave: T_FLOAT, 
        p_E_th : T_FLOAT, p_E_max : T_FLOAT, p_E_0 : T_FLOAT, p_sigma0 : T_FLOAT,
        p_y_a : T_FLOAT, p_P : T_FLOAT, p_y_w : T_FLOAT, p_y_0 : T_FLOAT, p_y_1 : T_FLOAT) -> T_VEC_FA:

    E_eV = CST.h_ * CST.c_ / wave / CST.eV2erg_# - p_E_th

    if E_eV < p_E_th:
        return 0.
    
#    if E_eV > p_E_max:
#        WARN_(f"E_eV = {E_eV}[eV] > p_E_max = {p_E_max}[eV] (max), set to p_E_max")

    x = E_eV / p_E_0 - p_y_0
    y = _sqrt( x*x + p_y_1*p_y_1 )
    F = ( (x-1)*(x-1) + p_y_w*p_y_w ) * y**(0.5*p_P - 5.5) * ( 1. + _sqrt(y/p_y_a) )**(-p_P)
    alpha = p_sigma0 * F
    alpha *= 1.E-18  # [Mb] to [cm^2]

    return alpha

_v1_fit_func_ = np_vec( _v1_fit_func_, **NP_VEC_KWGS)

def v1_photoioniz_cross_section_( wave : T_VEC_FA, nZ : T_INT, ioniz_stage: T_INT) -> T_VEC_FA:

    params = _V1_DATA[nZ][ioniz_stage]

    rate_coe = _v1_fit_func_(wave, *params)

    return rate_coe


