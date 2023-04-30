
#-------------------------------------------------------------------------------
# definition of functions for Level/Cont/Line indexing/table 
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#    2022/07/24   k.i.,u.k.   each_prof, extract_lprof, level_info
#-------------------------------------------------------------------------------

from ...ImportAll import *
from ...Struct import Atom as _Atom
from ...Struct.Container import CloudModel as _CloudModel
import numpy as _numpy



def Level_ctj_table_to_Level_info_( Level_ctj_table : T_LIST[T_TUPLE[T_STR,T_STR,T_STR]] ) -> T_DICT[T_STR, T_LIST[T_STR]]:

    Level_info : T_DICT[T_STR, T_LIST[T_STR]] = {
        "configuration" : [],
        "term" : [],
        "J" : [],
    }

    for i, ctj in enumerate( Level_ctj_table ):
        Level_info["configuration"].append( ctj[0] )
        Level_info["term"].append( ctj[1] )
        Level_info["J"].append( ctj[2] )

    return Level_info

#-------------------------------------------------------------
def each_prof(l:_CloudModel.CloudModel_Container, il:T_INT):
#-------------------------------------------------------------
# extract one line profile from l.wl_1D, l.prof_1D

    idx = l.Line_mesh_idxs[il,:]
    wl = l.wl_1D[idx[0]:idx[1]]
    prof = l.prof_1D[idx[0]:idx[1]]  #  Src[k] * (1. - _numpy.exp(-tau[:])) in CloudModel.SE_to_slab_0D_()
    return wl,prof

#-------------------------------------------------------------
def extract_lprof(l:_CloudModel.CloudModel_Container, wmin:T_FLOAT, wmax:T_FLOAT, dw:T_FLOAT, Ic:T_FLOAT=0.):
#-------------------------------------------------------------
# extract line profile in a wl interval from l-structure
# wmin,wmax,dw  -- in A 
# Ic  -- background intensity
# tau1[:] for return

    nw = int((wmax-wmin)/dw) + 1
    wl = _numpy.linspace(wmin,wmax,nw)
    profe = _numpy.zeros(nw)
    tau = _numpy.zeros(nw)

    w0 = l.w0 *1e8
    jj = _numpy.where((w0 > wmin) & (w0 < wmax))[0]

    for j in jj:
        idx1 = l.Line_mesh_idxs[j,0] 
        idx2 = l.Line_mesh_idxs[j,1]
        tau1 = _numpy.interp(wl, l.wl_1D[idx1:idx2]*1e8, l.tau_1D[idx1:idx2])
        # prof1 = np.interp(wl, l.wl_1D[idx1:idx2]*1e8, l.prof_1D[idx1:idx2])
        prof1 = l.Src[j] * (1.-_numpy.exp(-tau1))
        print('Ic, Src=',Ic,l.Src[j])
        tau += tau1
        profe += prof1
        
    prof = Ic*_numpy.exp(-tau) + profe
    
    return wl,prof,tau

#-------------------------------------------------------------
def level_info(atom:_Atom.Atom,wl0:T_FLOAT):  # get level info
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

def line_list(atom:_Atom.Atom, Aji_min:T_FLOAT = 1E4):
    selects : T_LIST[T_DICT[T_STR,T_INT | T_FLOAT | T_TUPLE[T_STR,T_STR,T_STR]]] = []
    for i in range(atom.nLine):
        wv0 = atom.Line['w0_AA'][i]
        Aji = atom.Line['AJI'][i]
        if Aji < Aji_min: continue
        ctj1, ctj2 = atom._ctj_table.Line[i]
        #print(f"{i:3d}  wavelength={wv0:10.2F}[A]  Aji={Aji:.1E}    {ctj1}<--{ctj2}")
        selects.append( {'id':i, 'w0_AA':wv0, 'Aji':Aji, 'ctj1':ctj1, 'ctj2':ctj2} )
    selects = sorted(selects, key=lambda x: x['w0_AA'])
    print(f"id   wavelength[A]  Aji[s-1]   lower level ctj    -     upper level ctj ")
    for sel in selects:
        id = sel['id']
        wv0 = sel['w0_AA']
        Aji = sel['Aji']
        ctj1 = sel['ctj1']
        ctj2 = sel['ctj2']
        #print(f"wavelength={wv0:10.2F}[A]  Aji={Aji:.1E}    {ctj1}<--{ctj2}")
        print(f"{id:3d}  {wv0:10.2F}  {Aji:.1E}    {ctj1}  -  {ctj2}")
    print(f"# selected = {len(selects)}")
    return(selects)
