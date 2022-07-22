
#-------------------------------------------------------------------------------
# definition of struct for storing cloud model result
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#    2022/07/20   k.i.   add Src, tau_1D
#-------------------------------------------------------------------------------

from ...ImportAll import *
from dataclasses import dataclass as _dataclass

@_dataclass(**STRUCT_KWGS_UNFROZEN)
class CloudModel_Container:
    """Cloud Model result container for
        - multiple transition
        - single spatial point
    """
    w0                 : T_ARRAY     # 1d, (nLine,) wavelength, [cm]
    tau_max            : T_ARRAY     # 1d, (nLine,), maximum optical depth, [-]
    Ibar               : T_ARRAY     # 1d, (nLine,), intensity profile integrated over wavelength, [erg/cm^2/Sr/s]
    Src                : T_ARRAY     # 1d, (nLine,), source function, [erg/cm^2/Sr/s]
    tau_1D             : T_ARRAY     # 1d, (sum_of_line_wavelength_mesh,), optical thickness profile, [-]
    prof_1D            : T_ARRAY     # 1d, (sum_of_line_wavelength_mesh,), out intensity profile, [erg/cm^2/Sr/cm/s]
    wl_1D              : T_ARRAY     # 1d, (sum_of_line_wavelength_mesh,), doppler shifted wavelength mesh , [cm]
    Line_mesh_idxs     : T_ARRAY     # 1d, (sum_of_line_wavelength_mesh,), [-]

