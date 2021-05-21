
#-------------------------------------------------------------------------------
# definition of functions for hydrogen atom model with degeneracy of N
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

from ...ImportAll import *
from ...Atomic import Hydrogen as _Hydrogen

import numpy as _numpy


_Hydrogen.ratio_Etran_to_Eionize_

def ratio_Etran_to_Eionize_( ni_arr : T_ARRAY, wave_arr : T_ARRAY ) -> T_ARRAY:
    
    nCont = wave_arr.shape[0]
    ratio = _numpy.empty(wave_arr.shape, T_FLOAT)

    for k in range(nCont):
        #E_ionize = CST.E_Rydberg_ / ni_arr[k]**2
        #E_tran =  CST.h_ * CST.c_ / wave_arr[k,:]
        #E_tran[:] +=  E_ionize - E_tran[0]
        #ratio[k,:] = E_tran[:] / E_ionize
        ratio[k,:] = _Hydrogen.ratio_Etran_to_Eionize_(ni_arr[k], wave_arr[k,:])

    return ratio

def compute_PI_cross_section(ni : T_ARRAY, meshCont : T_ARRAY) -> T_ARRAY:

    ## compute ratio of the transition energy to ionization energy
    Eratio = ratio_Etran_to_Eionize_(ni[:], meshCont[::])

    PI_alpha = _numpy.empty(Eratio.shape, dtype=T_FLOAT)
    for k in range(Eratio.shape[0]):
        PI_alpha[k,:] = _Hydrogen.PI_cross_section_(ni[k], Eratio[k,:], 1)

    return PI_alpha




