
#-------------------------------------------------------------------------------
# definition of struct for Atmosphere
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------


from ..ImportAll import *

from dataclasses import dataclass as _dataclass

#-------------------------------------------------------------------------------
# struct
#-------------------------------------------------------------------------------

@_dataclass(**STRUCT_KWGS_UNFROZEN)
class Atmosphere0D:

    Nh : T_FLOAT
    Ne : T_FLOAT
    Te : T_FLOAT
    Vd : T_FLOAT
    Vt : T_FLOAT

    ndim : T_INT     = 0
    is_gray : T_BOOL = True

    Tr : T_FLOAT     = 6.E3
    use_Tr : T_BOOL  = False

    doppler_shift_continuum : T_BOOL = False

    _coord_type : T_E_ATMOSPHERE_COORDINATE_TYPE = E_ATMOSPHERE_COORDINATE_TYPE.POINT

@_dataclass(**STRUCT_KWGS_UNFROZEN)
class AtmosphereC1D:

    Nh : T_ARRAY
    Ne : T_ARRAY
    Te : T_ARRAY
    Vd : T_ARRAY
    Vt : T_ARRAY

    is_gray : T_BOOL
    ndim : T_INT     = 1
    
    Tr : T_FLOAT     = 6.E3
    use_Tr : T_BOOL  = False

    doppler_shift_continuum : T_BOOL = False

    _coord_type : T_E_ATMOSPHERE_COORDINATE_TYPE = E_ATMOSPHERE_COORDINATE_TYPE.CARTESIAN
