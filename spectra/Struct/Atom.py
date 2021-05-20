
#-------------------------------------------------------------------------------
# definition of struct for Atom
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

from ..ImportAll import *
from ..Util.AtomUtils import AtomIO as _AtomIO

from dataclasses import dataclass as _dataclass

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

    Level : None
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
    if not conf_path.endswith(".conf"):
        raise ValueError("`conf_path` should be a string ends with '.conf'")

    return atom , path_dict



