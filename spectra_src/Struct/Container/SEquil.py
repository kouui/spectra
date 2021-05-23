
#-------------------------------------------------------------------------------
# definition of struct for storing Statistical Equilibrium result
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

from ...ImportAll import *
from dataclasses import dataclass as _dataclass


@_dataclass(**STRUCT_KWGS_UNFROZEN)
class SE_Container:
    """Statistical Equilibrium result container for
        - single spatial point
    """

    n_SE                 : T_ARRAY # 1d (nLevel,), [/cm^3]
    n_LTE                : T_ARRAY # 1d (nLevel,), [/cm^3]

    nj_by_ni             : T_ARRAY # 1d (nLine+nCont,), [-]

    ## wavelength mesh of line transition
    wave_mesh_shifted_1d : T_ARRAY # 1d (sum_of_line_wavelength_mesh,), [cm]
    ## absorption profile of line transition
    absorb_prof_1d       : T_ARRAY # 1d (sum_of_line_wavelength_mesh,), [/cm]
    ## index array of line transition
    Line_mesh_idxs       : T_ARRAY # 1d (sum_of_line_wavelength_mesh,), [-]


    Jbar                 : T_ARRAY # 1d (nLine,), [erg/cm^2/Sr/s]

@_dataclass(**STRUCT_KWGS_UNFROZEN)
class TranRates_Container:
    """Statistical Equilibrium (Transition Rates) result container for
        - single spatial point
    """
    Rji_spon             : T_ARRAY # 1d (nLine+nCont), [/s]
    Rji_stim             : T_ARRAY # 1d (nLine+nCont), [/s]
    Rij                  : T_ARRAY # 1d (nLine+nCont), [/s]

    Cji_Ne               : T_ARRAY # 1d (nLine+nCont), [/s]
    Cij_Ne               : T_ARRAY # 1d (nLine+nCont), [/s]

