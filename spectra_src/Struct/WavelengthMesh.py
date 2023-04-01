
#-------------------------------------------------------------------------------
# definition of struct for Wavelength Mesh
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------


from ..ImportAll import *
from dataclasses import dataclass as _dataclass


@_dataclass(**STRUCT_KWGS_UNFROZEN)
class Wavelength_Mesh:

    ##: initialized given the Atom.Cont
    #   could be modified due to doppler shift
    #   currently, we assume that Cont_mesh will not be affected by doppler shift
    Cont_mesh       : T_ARRAY # 2d
    Cont_Coe        : T_ARRAY # struct

    ##: initialized given the atmosphere model
    Line_mesh       : T_ARRAY # 1d, doppler width unit, without doppler shift
    Line_mesh_idxs  : T_ARRAY # 2d,   (nLine, 2)
    Line_absorb_prof: T_ARRAY # 1d, un-normalized. `/= dopWidth_cm` --> normalized

    Line_Coe        : T_ARRAY # struct
    
    ##: we need Line_mesh_share to deal with non-uniform atmoshere
    Line_mesh_share      : T_ARRAY # 1d,
    Line_mesh_share_idxs : T_ARRAY # 2d,   (nLine, 2)


import numpy as _numpy
from ..Util import MeshUtil as _MeshUtil

_N_CONT_MESH : T_INT = 41
_N_LINE_MESH : T_INT = 41
_LINE_MESH_QCORE  : T_FLOAT = 2.5
_LINE_MESH_QWING  : T_FLOAT = 10.
_LINE_MESH_TYPE   : T_E_ABSORPTION_PROFILE_TYPE = \
    E_ABSORPTION_PROFILE_TYPE.GAUSSIAN


def init_Wave_Mesh_(Cont : T_ARRAY, Line : T_ARRAY, RL_Coe : T_ARRAY) -> Wavelength_Mesh:
    

    ## : Cont_mesh
    dtype  = _numpy.dtype([
                          ('idxI',T_INT),         #: level index, the Level index of lower level
                          ('idxJ',T_INT),         #: level index, the Level index of upper level
                          ('w0',T_FLOAT),         #: central wavelength, [cm]
                          ('nLambda', T_INT),     #: number of meaningful wavelength mesh point
                          ])
    nCont = Cont.shape[0]
    Cont_mesh = _numpy.zeros((nCont, _N_CONT_MESH), dtype=DT_NB_FLOAT)
    Cont_Coe  = _numpy.empty(nCont, dtype=dtype)
    if nCont > 0:  # has_continuum == True

        for k in range(nCont):
            #mesh = _MeshUtil.make_continuum_mesh_( _N_CONT_MESH ) # in limit wavelength unit
            #w0 = Cont["w0"][k]
            Cont_mesh[k,:] = Cont["w0"][k] * _MeshUtil.make_continuum_mesh_( _N_CONT_MESH )

            Cont_Coe["idxI"][k] = Cont["idxI"][k]
            Cont_Coe["idxJ"][k] = Cont["idxJ"][k]
            Cont_Coe["w0"][k]   = Cont["w0"][k]
            Cont_Coe["nLambda"][k] = _N_CONT_MESH

    ## : Line_mesh, only transition defined in RL will be considered
    nLine = Line.shape[0]
    n_total_line_mesh = 0
    i_RL : T_INT
    for k in range(nLine):
        res = _numpy.where( RL_Coe["lineIndex"][:] == k )[0]
        try:
            i_RL = res[0]
        except IndexError:
            i_RL = -1
        
        ##: Line_mesh contains mesh for all lines
        if i_RL != -1:
            n_total_line_mesh += RL_Coe["nLambda"][i_RL]
        else:
            n_total_line_mesh += _N_LINE_MESH

    dtype  = _numpy.dtype([
                          ('idxI',T_INT),         #: level index, the Level index of lower level
                          ('idxJ',T_INT),         #: level index, the Level index of upper level
                          ('w0',T_FLOAT),         #: central wavelength, [cm]
                          ('ProfileType',T_INT),  #: 0: Voigt; 1: Gaussian
                          ('qcore', T_FLOAT),
                          ('qwing', T_FLOAT),
                          ('nLambda', T_INT),     #: number of meaningful wavelength mesh point
                          ])
    Line_mesh        = _numpy.empty(n_total_line_mesh, dtype=DT_NB_FLOAT)
    Line_absorb_prof = _numpy.empty(n_total_line_mesh, dtype=DT_NB_FLOAT)  # not defined
    Line_mesh_idxs   = _numpy.empty((nLine, 2), dtype=DT_NB_INT)
    Line_Coe         = _numpy.empty(nLine, dtype=dtype)
    nLambda : T_INT
    qcore : T_FLOAT
    qwing : T_FLOAT
    proftype : T_INT
    bias : T_INT = 0
    for k in range(nLine):
        res = _numpy.where( RL_Coe["lineIndex"][:] == k )[0]
        try:
            i_RL = res[0]
        except IndexError:
            i_RL = -1

        if i_RL != -1:
            nLambda = RL_Coe["nLambda"][i_RL]
            qcore   = RL_Coe["qcore"][i_RL]
            qwing   = RL_Coe["qwing"][i_RL]
            proftype= RL_Coe["ProfileType"][i_RL]
        else:
            nLambda = _N_LINE_MESH
            qcore   = _LINE_MESH_QCORE
            qwing   = _LINE_MESH_QWING
            proftype= int( _LINE_MESH_TYPE )
        
        wave_mesh = _MeshUtil.make_full_line_mesh_( nLambda, qcore, qwing )
        Line_mesh_idxs[k,0] = bias
        Line_mesh_idxs[k,1] = Line_mesh_idxs[k,0] + nLambda
        bias = Line_mesh_idxs[k,1]
        Line_mesh[Line_mesh_idxs[k,0] : Line_mesh_idxs[k,1]] = wave_mesh[:]
        
        Line_Coe["idxI"][k] = Line["idxI"][k]
        Line_Coe["idxJ"][k] = Line["idxJ"][k]
        Line_Coe["w0"][k] = Line["w0"][k]
        Line_Coe["ProfileType"][k] = proftype
        Line_Coe["qcore"][k] = qcore
        Line_Coe["qwing"][k] = qwing
        Line_Coe["nLambda"][k] = nLambda

    wavelength_mesh = Wavelength_Mesh(
        Cont_mesh = Cont_mesh,
        Cont_Coe  = Cont_Coe,
        Line_mesh = Line_mesh,
        Line_Coe  = Line_Coe,
        Line_absorb_prof = Line_absorb_prof,
        Line_mesh_idxs = Line_mesh_idxs,
        Line_mesh_share = _numpy.empty(0, dtype=DT_NB_FLOAT),
        Line_mesh_share_idxs = _numpy.empty((nLine,2), dtype=DT_NB_INT)
    )

    return wavelength_mesh



    


