
#-------------------------------------------------------------------------------
# definition of struct for Wavelength Mesh
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------


from ..ImportAll import *
from dataclasses import dataclass as _dataclass

@_dataclass(**STRUCT_KWGS)
class Radiative_Line:

    nRadiativeLine   : T_INT
    Coe              : T_ARRAY # struct array

@_dataclass(**STRUCT_KWGS_UNFROZEN)
class Wavelength_Mesh:

    ##: initialized given the Atom.Cont
    #   could be modified due to doppler shift
    #   currently, we assume that Cont_mesh will not be affected by doppler shift
    Cont_mesh       : T_ARRAY # 2d

    ##: initialized given the atmosphere model
    Line_mesh       : T_ARRAY # 1d, doppler width unit, without doppler shift
    Line_mesh_idxs  : T_ARRAY # 2d,   (nLine, 2)
    
    ##: we need Line_mesh_share to deal with non-gray atmoshere
    Line_mesh_share      : T_ARRAY # 1d,
    Line_mesh_share_idxs : T_ARRAY # 2d,   (nLine, 2)
