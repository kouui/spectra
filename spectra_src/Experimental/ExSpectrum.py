
#-------------------------------------------------------------------------------
# spectrum struct
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/08/10   u.k.   
#-------------------------------------------------------------------------------

from ..ImportAll import *

from ..Struct import Atom as _Atom
from ..Struct import WavelengthMesh as _WavelengthMesh

from ..Atomic import BasicP as _BasicP

from dataclasses import dataclass as _dataclass

import numpy as _numpy

import matplotlib.pyplot as plt 

_MAX_OVERLAP_LINE = 10
_MAX_OVERLAP_CONT = 20

@_dataclass(**STRUCT_KWGS_UNFROZEN)
class Spectrum:                                                         ##: struct contains line & cont mesh

    nSpectrum : T_INT                                                   ##: number of wavelength mesh
    spectrum : T_ARRAY                                                  ##: (nSpectrum,), wavelength mesh in cm unit
    line_mesh_dop : T_ARRAY                                                ##: (nSpectrum-#cont_mesh.size,), wavelength mesh in doppler width unit

    belonging_line : T_ARRAY                                            ##: (nSpectrum,_MAX_OVERLAP_LINE), the wavelength point belonging to which b-b transition
    belonging_cont : T_ARRAY                                            ##: (nSpectrum,_MAX_OVERLAP_CONT), the wavelength point belonging to which b-f transition
    
##    index_line : T_ARRAY                                                ##: forgot what this is used for
##    index_cont : T_ARRAY                                                ##: forgot what this is used for

    RL_lineindex : T_ARRAY

    dop_width_cm : T_ARRAY


def init_spectrum_(atom :_Atom.Atom, wMesh: _WavelengthMesh.Wavelength_Mesh, Te : T_FLOAT, Vt : T_FLOAT, Vd: T_FLOAT ):
    
    Line_mesh_idxs = wMesh.Line_mesh_idxs
    Line_mesh = wMesh.Line_mesh
    Cont_mesh = wMesh.Cont_mesh
    nSpec = wMesh.Line_mesh.size + wMesh.Cont_mesh.reshape(-1).size
    spectrum = _numpy.empty(nSpec, dtype=DT_NB_FLOAT)

    mass = atom.Mass
    nCont = atom.nCont

    ##: append line mesh (cm)
    nRL = atom.nRL
    RL_coe = atom.RL.Coe
    RL_lineindex  = RL_coe["lineIndex"].copy()
    
    bias = 0
    line_mesh_dop= _numpy.empty((nSpec,),dtype=DT_NB_FLOAT)
    line_mesh_cm = _numpy.empty_like(Line_mesh)
    dop_width_cm = _numpy.empty((atom.nLine,),dtype=DT_NB_FLOAT)

    for k in RL_coe["lineIndex"][:]:                                              ##: loop over b-b lines defined in RL
        w0 = atom.Line["w0"][k]
        dopWidth_cm = _BasicP.doppler_width_(w0, Te, Vt, mass)                    ##: doppler width [cm] with typical Te and Vt
        i1, i2 = Line_mesh_idxs[k,:]
        nwave = i2 - i1
        spectrum[bias:bias+nwave] = Line_mesh[i1 : i2] * dopWidth_cm + w0         ##: wavelength mesh [cm]
        spectrum[bias:bias+nwave] += ( w0 * Vd / CST.c_ )
        line_mesh_dop[bias:bias+nwave] = Line_mesh[i1 : i2]
        #line_mesh_cm[bias:bias+nwave] = spectrum[bias:bias+nwave]
        line_mesh_cm[i1:i2] = spectrum[bias:bias+nwave]
        dop_width_cm[k] = dopWidth_cm

        bias += nwave

    ##: append cont mesh (cm)
    for k in range(atom.nCont):
        nwave = Cont_mesh[k,:].size
        spectrum[bias:bias+nwave] = Cont_mesh[k,:]                                ##: since we will sort spectrum later, the order of mesh does not matter
        bias += nwave
    spectrum = spectrum[:bias]                                                    ##: cut part of n_nonRL_line_mesh
    line_mesh_dop = line_mesh_dop[:bias]
    
    
    # sort and remove duplicates (almost no duplicates since we are working with float64)
    ids = _numpy.argsort(spectrum)                                                ##: sort spectrum and line_mesh_dop with spectrum ascending order   
    spectrum = spectrum[ids].copy()
    line_mesh_dop = line_mesh_dop[ids].copy()
    # idxs = _numpy.where(spectrum[1:] == spectrum[:-1])[0]
    # idxs[:] += 1
    # mask = _numpy.ones_like(spectrum, dtype=T_BOOL)
    # mask[idxs] = False
    # spectrum = spectrum[mask]

    nSpectrum = spectrum.size

    # decide belongings line (lines of single atom)
    belonging_line = -1 * _numpy.ones((nSpectrum,_MAX_OVERLAP_LINE), dtype=DT_NB_INT)    ## value > -1 : overlapping line index
##    index_line = -1 * _numpy.ones((nSpectrum,nRL,2), dtype=DT_NB_INT)
    for k in RL_coe["lineIndex"][:]:
        i1, i2 = Line_mesh_idxs[k,:]
        blue = line_mesh_cm[i1]
        red = line_mesh_cm[i2-1]
        i_spec_blue = _numpy.where(spectrum[:]>=blue)[0][0]
        i_spec_red  = _numpy.where(spectrum[:]<=red)[0][-1] + 1
        for ispec in range(i_spec_blue, i_spec_red):
            for iline in range(_MAX_OVERLAP_LINE):
                
                if belonging_line[ispec,iline] == -1:
                    belonging_line[ispec,iline] = k
                    break
                if iline == _MAX_OVERLAP_LINE-1:
                    raise ValueError("_MAX_OVERLAP_LINE not enough")
        
##        index_line[i_spec_blue:i_spec_red,k,:] = i_spec_blue, i_spec_red

    # decide belongings cont
    belonging_cont = -1 * _numpy.ones((nSpectrum,_MAX_OVERLAP_LINE), dtype=DT_NB_INT)
##    index_cont = _numpy.empty((nSpectrum,nRL,2), dtype=DT_NB_INT)
    for k in range(nCont):
        blue = Cont_mesh[k,-1]
        red  = Cont_mesh[k,0]
        i_spec_blue = _numpy.where(spectrum[:]>=blue)[0][0]
        i_spec_red  = _numpy.where(spectrum[:]<=red)[0][-1] + 1
        for ispec in range(i_spec_blue, i_spec_red):
            for icont in range(_MAX_OVERLAP_CONT):

                if belonging_cont[ispec,icont] == -1:
                    belonging_cont[ispec,icont] = k
                    break
                if icont == _MAX_OVERLAP_CONT-1:
                    raise ValueError("_MAX_OVERLAP_CONT not enough")

##        index_cont[i_spec_blue:i_spec_red,k,:] = i_spec_blue, i_spec_red
    Spec = Spectrum(
        nSpectrum=nSpectrum,
        spectrum = spectrum,
        belonging_line = belonging_line,
        belonging_cont = belonging_cont,
##        index_line = index_line,
##        index_cont = index_cont,
        RL_lineindex = RL_lineindex,
        dop_width_cm = dop_width_cm,
        line_mesh_dop=line_mesh_dop,
    )

    return Spec
