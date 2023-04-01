

from ..ImportAll import *
#from ..Elements import ELEMENT_DICT as _ELEMENT_DICT
#from ..Struct import Atmosphere as _Atmosphere
#from dataclasses import dataclass as _dataclass
from ..Struct.Container import CloudModel as _CloudModel
from ..Struct import Atom as _Atom

import numpy as _numpy

#-------------------------------------------------------------
def each_prof(l:_CloudModel.CloudModel_Container, il:int):
#-------------------------------------------------------------
# extract one line profile from l.wl_1D, l.prof_1D

    idx = l.Line_mesh_idxs[il,:]
    wl = l.wl_1D[idx[0]:idx[1]]
    prof = l.prof_1D[idx[0]:idx[1]]  #  Src[k] * (1. - _numpy.exp(-tau[:])) in CloudModel.SE_to_slab_0D_()
    return wl,prof

#-------------------------------------------------------------
def extract_lprof(l:_CloudModel.CloudModel_Container,wmin:float,wmax:float,dw:float, Ic:float=0.):
#-------------------------------------------------------------
# extract line profile in a wl interval from l-structure
# wmin,wmax,dw  -- in AA 
# Ic  -- backgroung intensity

    nw = int((wmax-wmin)/dw) + 1
    wl = _numpy.linspace(wmin,wmax,nw)
    prof = _numpy.zeros(nw)

    w0 = l.w0 *1e8    # [cm] -> [AA]
    jj = _numpy.where((w0 > wmin) & (w0 < wmax))[0]

    for j in jj:
        idx1 = l.Line_mesh_idxs[j,0] 
        idx2 = l.Line_mesh_idxs[j,1]
        tau1 = _numpy.interp(wl, l.wl_1D[idx1:idx2]*1e8, l.tau_1D[idx1:idx2])
        # prof1 = np.interp(wl, l.wl_1D[idx1:idx2]*1e8, l.prof_1D[idx1:idx2])
        prof1 = Ic * _numpy.exp(-tau1) + l.Src[j] * (1.-_numpy.exp(-tau1))
        print(f"l.Src[{j}]={l.Src[j]:+.2E}")
        prof[:] += prof1
    return wl,prof

#-------------------------------------------------------------
def level_info(atom:_Atom.Atom,wl0:float):  # get level info
#-------------------------------------------------------------
#  return lower/upper level index and term-J of a line at wl0 [A]

    w0_AA = atom.Line["w0_AA"]/1.00027    # wl in atmosphere
    line_idx = _numpy.argmin(abs(wl0-w0_AA))
    print("Line# =",line_idx)
    idxI = atom.Line[line_idx]["idxI"]
    idxJ = atom.Line[line_idx]["idxJ"]
    ltj = atom._ctj_table.Level[idxI]
    utj = atom._ctj_table.Level[idxJ]
    lower = ltj[1]+ltj[2]
    upper = utj[1]+utj[2]
    
    return idxI,lower,idxJ,upper

        