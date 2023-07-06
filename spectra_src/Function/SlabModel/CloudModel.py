
#-------------------------------------------------------------------------------
# definition of functions to perform statistical equilibrium
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#    2022/07/20   k.i.   add Src, tau_1D
# 0.1.2
#    2023/07/04   u.k.   when Aji equals to 0, set Src to 0 to avoid the zero division warning
#-------------------------------------------------------------------------------

from typing import Container
from ...ImportAll import *

#from ..SEquil import SELib as _SELib
from ...Struct import Atom as _Atom
from ...Struct import Atmosphere as _Atmosphere
from ...Struct import Container as _Container
from ...Math import Integrate as _Integrate

import numpy as _numpy


def SE_to_slab_0D_(atom : _Atom.Atom, atmos : _Atmosphere.Atmosphere0D,
                   SE_con : _Container.SE_Container, depth : T_FLOAT) -> _Container.CloudModel_Container:

    nLine = atom.nLine
    Line = atom.Line

    N_ele = atmos.Nh * atom.Abun

    Line_mesh_idxs       = SE_con.Line_mesh_idxs
    wave_mesh_shifted_1d = SE_con.wave_mesh_shifted_1d
    absorb_prof_1d       = SE_con.absorb_prof_1d

    ## 1. obtain the upper/lower level population for line transitions
    nj : T_ARRAY = SE_con.n_SE[Line["idxJ"][:]]
    ni : T_ARRAY = SE_con.n_SE[Line["idxI"][:]]

    ## 2. compute extinction coefficient alpha
    hv    : T_ARRAY = CST.h_ * Line["f0"][:]
    Bij   : T_ARRAY = Line["BIJ"][:]
    Bji   : T_ARRAY = Line["BJI"][:]
    alp0  : T_ARRAY = hv / (4. * CST.pi_) * (Bij * ni - Bji * nj ) * N_ele

    ## 3. compute line source function
    Aji   : T_ARRAY = Line["AJI"][:]
    #Src   : T_ARRAY = ( Aji * nj ) / ( Bij * ni - Bji * nj )
    Src   : T_ARRAY = _numpy.zeros_like(Aji)
    for k in range(0,nLine):
        if Aji[k] <= 0. :
            Src[k] = 0.
        else:
            Src[k] = ( Aji[k] * nj[k] ) / ( Bij[k] * ni[k] - Bji[k] * nj[k] )

    ## 4. compute optical depth given the thichness of the slab
    ## 5. compute the line profile
    arr_w0        = _numpy.empty(nLine, dtype=DT_NB_FLOAT)
    arr_tau_max   = _numpy.empty(nLine, dtype=DT_NB_FLOAT)
    arr_Ibar      = _numpy.empty(nLine, dtype=DT_NB_FLOAT)
    arr_prof_1D   = _numpy.empty_like( absorb_prof_1d )
    arr_tau_1D    = _numpy.empty_like( absorb_prof_1d )
    #arr_wl_1D     = _numpy.empty_like( absorb_prof_1d )
    for k in range(0,nLine):
        i1 = Line_mesh_idxs[k,0]
        i2 = Line_mesh_idxs[k,1]

        tau = depth * alp0[k] * absorb_prof_1d[i1:i2]
        wl  = wave_mesh_shifted_1d[i1:i2]
        prof = Src[k] * (1. - _numpy.exp(-tau[:]))
        #l.tau0[i] = np.max(tau)
        #l.prof[i][:] = S[i] * (1. - np.exp(-tau[:]))
        Ibar = _Integrate.trapze_(prof[:], wl[:])

        # store value
        arr_w0[k] = Line["w0"][k]
        arr_tau_max[k] = tau[:].max()
        arr_Ibar[k] = Ibar
        arr_prof_1D[i1:i2] = prof[:]
        arr_tau_1D[i1:i2] = tau[:]


    cloud_con = _Container.CloudModel_Container(
        w0                 = arr_w0, 
        tau_max            = arr_tau_max, 
        Ibar               = arr_Ibar, 
        Src                = Src, 
        tau_1D             = arr_tau_1D,
        prof_1D            = arr_prof_1D,
        wl_1D              = wave_mesh_shifted_1d.copy(),
        Line_mesh_idxs     = Line_mesh_idxs.copy(),
    )

    return cloud_con














