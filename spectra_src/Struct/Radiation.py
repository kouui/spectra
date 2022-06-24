

#-------------------------------------------------------------------------------
# definition of struct for storing radiation
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
# 0.1.1
#    2022/01/07   u.k.   modified atlas(backRad)
#-------------------------------------------------------------------------------


from ..ImportAll import *
from dataclasses import dataclass as _dataclass


@_dataclass(**STRUCT_KWGS_UNFROZEN)
class Radiation:

    backRad : T_ARRAY  # 2d, (2, n_wavelength)

    PI_intensity : T_ARRAY  # 2d, (nCont, _N_CONT_MESH )


from . import Atom as _Atom
from . import WavelengthMesh as _WavelengthMesh
from . import Atmosphere as _Atmosphere

from ..Atomic import PhotoIonize as _PhotoIonize
from ..Atomic import LTELib as _LTELib

import numpy as _numpy
import os

def init_Radiation_(atmos : T_UNION[_Atmosphere.Atmosphere0D,_Atmosphere.AtmosphereC1D],
                    wMesh : _WavelengthMesh.Wavelength_Mesh) -> Radiation:

    root=CFG._ROOT_DIR

    backRad = _numpy.load(  os.path.join(root,"data/intensity/atlas/QS/atlas_QS/20220107.npy") )
    #backRad[0,:] *= 1E-8
    #backRad[1,:] *= 2.5*intensity_fac
    
    Cont_mesh = wMesh.Cont_mesh
    if atmos.use_Tr:
        Tr = atmos.Tr
        PI_intensity = _LTELib.planck_cm_(Cont_mesh[:,:], Tr)
    else:
        #  TODO : average backRad for PI_intensity ?
        PI_intensity = _PhotoIonize.interpolate_PI_intensity_(backRad[:,:], Cont_mesh[:,:])

    radiation = Radiation(
        backRad      = backRad,
        PI_intensity = PI_intensity,
    )

    return radiation

