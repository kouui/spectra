
#-------------------------------------------------------------------------------
# definition of struct for Atom
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

@_dataclass(**STRUCT_KWGS)
class Collisional_Transition:

    _transition_type    : T_INT
    _transition_source  : T_INT
    _transition_formula : T_INT

    Te_table            : T_ARRAY # 1d
    Omega_table         : T_ARRAY # 2d
    Coe                 : T_ARRAY # struct array

@_dataclass(**STRUCT_KWGS)
class Photo_Ionization:

    alpha_table      : T_ARRAY # 1d
    alpha_table_idxs : T_ARRAY # 2d,  (nCont, 2)

    #intensity   : T_ARRAY # 2d   (nCont, nContMesh)

    alpha_interp     : T_ARRAY # 2d   (nCont, nContMesh)
    # interpolated photoionization cross section
    # no matter there is doppler shift, alpha_interp always 
    # starts from the cross section at ionization limit



@_dataclass(**STRUCT_KWGS)
class ATOMIC_DATA_SOURCE:

    AJI : T_INT
    CE  : T_INT
    CI  : T_INT
    PI  : T_INT


@_dataclass(**STRUCT_KWGS)
class CTJ_Table:

    Level : None
    Line  : None
    Cont  : None

@_dataclass(**STRUCT_KWGS)
class Index_Table:

    Line  : None
    Cont  : None



@_dataclass(**STRUCT_KWGS)
class Atom:

    #Title : T_STR
    #Element : T_STR

    Z      : T_INT
    Mass   : T_FLOAT
    Abun   : T_FLOAT

    nLevel : T_INT
    nLine  : T_INT
    nCont  : T_INT
    nTran  : T_INT
    

    Level  : T_ARRAY    # struct array
    Line   : T_ARRAY    # struct array
    Cont   : T_ARRAY    # struct array

    _has_continuum : T_BOOL

    _atomic_data_source : ATOMIC_DATA_SOURCE
    _atom_type          : T_INT

    _ctj_table : CTJ_Table
    _idx_table : Index_Table

    CE         : Collisional_Transition
    CI         : Collisional_Transition
    PI         : Photo_Ionization

#-------------------------------------------------------------------------------
# init function
#-------------------------------------------------------------------------------

from ..Util.AtomUtils import AtomIO as _AtomIO
import numpy as _numpy
from collections import OrderedDict as _OrderedDict

def init_Atom_(conf_path : T_STR, is_hydrogen : T_BOOL = False ) -> T_TUPLE[Atom, T_DICT[T_STR, T_UNION[None,T_STR]]] : 
    """given the *.conf file, create the Atom struct

    Parameters
    ----------
    conf_path : T_STR
        path to the *.conf file
    is_hydrogen : T_BOOL, optional
        whether is Hydrogen atom, by default False

    Returns
    -------
    T_TUPLE[Atom, T_DICT[T_STR, T_UNION[None,T_STR]]]
    
    atom : Atom
        the Atom struct
    path_dict : T_DICT[T_STR, T_UNION[None,T_STR]]
        a dictionary of the path of data files
    """
    # path dict
    #--------------------
    path_dict = _AtomIO.read_conf_(conf_path)
    
    if is_hydrogen:
        atom_type = E_ATOM.HYDROGEN
    else:
        atom_type = E_ATOM.NORMAL

    # read Level
    #--------------------
    if path_dict["Level"] is None:
        raise ValueError("Lack of .Level file")
        
    Z, Mass, Abun, nLevel, Level, Level_info_table = _AtomIO.read_Atom_Level_(path_dict["Level"])

    # nTran nLine nCont
    #--------------------
    nLine, nCont, nTran, _has_continuum = _AtomIO.nLine_nCont_nTran( Level["stage"] )
    if not _has_continuum:
        raise ValueError("Currently we don't support Atomic Model without comtinuum.")
    


    atom = Atom()
    return atom , path_dict



