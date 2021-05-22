

#-------------------------------------------------------------------------------
# definition of struct for storing radiation
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------


from ..ImportAll import *
from dataclasses import dataclass as _dataclass


@_dataclass(**STRUCT_KWGS_UNFROZEN)
class Radiation:

    backRad : T_ARRAY  # 2d, (2, n_wavelength)

    PI_intensity : T_ARRAY  # 2d, (nCont, _N_CONT_MESH )


from . import Atom as _Atom
from . import WavelengthMesh as _WavelengthMesh
from ..Atomic import PhotoIonize as _PhotoIonize

import numpy as _numpy
import os

def init_Radiation_(atom : _Atom.Atom, wMesh : _WavelengthMesh.Wavelength_Mesh) -> Radiation:

    root=CFG._ROOT_DIR

    backRad = _numpy.load(  os.path.join(root,"data/intensity/atlas/QS/atlas_QS.npy") )
    backRad[0,:] *= 1E-8
    
    Cont_mesh = wMesh.Cont_mesh
    PI_intensity = _PhotoIonize.interpolate_PI_intensity_(backRad[:,:], Cont_mesh[:,:])

    radiation = Radiation(
        backRad      = backRad,
        PI_intensity = PI_intensity,
    )

    return radiation

