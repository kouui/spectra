
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
class Radiative_Line:

    nRadiativeLine   : T_INT
    Coe              : T_ARRAY # struct array


@_dataclass(**STRUCT_KWGS)
class Collisional_Transition:

    _transition_type    : T_E_COLLISIONAL_TRANSITION
    _transition_source  : T_E_COLLISIONAL_TRANSITION_SOURCE
    _transition_formula : T_E_COLLISIONAL_TRANSITION_FORMULA

    Te_table            : T_ARRAY # 1d
    Omega_table         : T_ARRAY # 2d
    Coe                 : T_ARRAY # struct array

@_dataclass(**STRUCT_KWGS)
class Photo_Ionization:

    alpha_table      : T_ARRAY # 2d,  (2, ?)
    alpha_table_idxs : T_ARRAY # 2d,  (nCont, 2)
    Coe              : T_ARRAY # struct

    #intensity   : T_ARRAY # 2d   (nCont, nContMesh)

    alpha_interp     : T_ARRAY # 2d   (nCont, nContMesh)
    # interpolated photoionization cross section
    # no matter there is doppler shift, alpha_interp always 
    # starts from the cross section at ionization limit



@_dataclass(**STRUCT_KWGS)
class ATOMIC_DATA_SOURCE:

    AJI : T_E_ATOMIC_DATA_SOURCE
    CE  : T_E_ATOMIC_DATA_SOURCE
    CI  : T_E_ATOMIC_DATA_SOURCE
    PI  : T_E_ATOMIC_DATA_SOURCE


@_dataclass(**STRUCT_KWGS)
class CTJ_Table:

    Level : T_CTJ_TABLE
    Line  : T_CTJ_PAIR_TABLE
    Cont  : T_CTJ_PAIR_TABLE

@_dataclass(**STRUCT_KWGS)
class Index_Table:

    Line  : T_IDX_PAIR_TABLE
    Cont  : T_IDX_PAIR_TABLE



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
    nRL    : T_INT 
    

    Level  : T_ARRAY    # struct array
    Line   : T_ARRAY    # struct array
    Cont   : T_ARRAY    # struct array

    _has_continuum      : T_BOOL

    _atomic_data_source : ATOMIC_DATA_SOURCE
    _atom_type          : T_E_ATOM

    _ctj_table : CTJ_Table
    _idx_table : Index_Table

    CE         : Collisional_Transition
    CI         : Collisional_Transition
    PI         : Photo_Ionization
    RL         : Radiative_Line

#-------------------------------------------------------------------------------
# init function
#-------------------------------------------------------------------------------

from ..Util.AtomUtils import AtomIO as _AtomIO
from . import WavelengthMesh as _WavelengthMesh

def init_Atom_(conf_path : T_STR, is_hydrogen : T_BOOL = False 
              ) -> T_TUPLE[Atom, _WavelengthMesh.Wavelength_Mesh, T_DICT[T_STR, T_UNION[None,T_STR]]] : 
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

    _atom_type : T_E_ATOM
    if is_hydrogen:
        _atom_type = E_ATOM.HYDROGEN
    else:
        _atom_type = E_ATOM.NORMAL

    # read Level
    #--------------------
    if path_dict["Level"] is None:
        raise ValueError("Lack of .Level file")
        
    Z, Mass, Abun, nLevel, Level, Level_info_table = _AtomIO.make_Atom_Level_(path_dict["Level"])

    # nTran nLine nCont
    #--------------------
    nLine, nCont, nTran, _has_continuum = _AtomIO.nLine_nCont_nTran_( Level["stage"] )
    if not _has_continuum:
        raise ValueError("Currently we don't support Atomic Model without comtinuum.")
    
    # ctj and idx table
    #--------------------
    Line_idx_table, Line_ctj_table, Cont_idx_table, Cont_ctj_table = \
        _AtomIO.prepare_idx_ctj_mapping_(Level_info_table, Level["stage"], Level["isGround"], nLine, nCont)

    _ctj_table = CTJ_Table(Level=Level_info_table, Line=Line_ctj_table, Cont=Cont_ctj_table)
    _idx_table = Index_Table(Line=Line_idx_table, Cont=Cont_idx_table)
    
    # make Cont
    #--------------------
    Cont = _AtomIO.make_Atom_Cont_(nCont, Cont_idx_table, Level)
    
    # read Aji
    #--------------------
    Line, data_source_Aji = \
        _AtomIO.make_Atom_Line_(path_dict["Aji"], Level,Line_idx_table,Line_ctj_table,_atom_type)
    
    # read CE
    #--------------------
    Te_table, Omega_table, Coe, _transition_type, _transition_source, _transition_formula = \
        _AtomIO.make_Atom_CECI_(path_dict["CEe"],"CE",nLine,Line,Level,Level_info_table,Line_ctj_table)
    CE = Collisional_Transition(
        _transition_type = _transition_type,
        _transition_source = _transition_source,
        _transition_formula = _transition_formula,
        Te_table = Te_table,
        Omega_table = Omega_table,
        Coe = Coe
    )
    data_source_CE : T_E_ATOMIC_DATA_SOURCE
    if Te_table.size == 0:
        data_source_CE = E_ATOMIC_DATA_SOURCE.CALCULATE
    else:
        data_source_CE = E_ATOMIC_DATA_SOURCE.EXPERIMENT
    
    del Te_table, Omega_table, Coe, _transition_type, _transition_source, _transition_formula
    # read CI
    #--------------------
    Te_table, Omega_table, Coe, _transition_type, _transition_source, _transition_formula = \
        _AtomIO.make_Atom_CECI_(path_dict["CIe"],"CI",nCont,Cont,Level,Level_info_table,Cont_ctj_table)
    CI = Collisional_Transition(
        _transition_type = _transition_type,
        _transition_source = _transition_source,
        _transition_formula = _transition_formula,
        Te_table = Te_table,
        Omega_table = Omega_table,
        Coe = Coe
    )
    data_source_CI : T_E_ATOMIC_DATA_SOURCE
    if Te_table.size == 0:
        data_source_CI = E_ATOMIC_DATA_SOURCE.CALCULATE
    else:
        data_source_CI = E_ATOMIC_DATA_SOURCE.EXPERIMENT
    
    del Te_table, Omega_table, Coe, _transition_type, _transition_source, _transition_formula
    
    # read radiative line
    #--------------------
    Coe, nRadiativeLine = \
        _AtomIO.make_Atom_RL_(path_dict["RadiativeLine"],Level_info_table,Line_ctj_table)
    RL = Radiative_Line(nRadiativeLine=nRadiativeLine, Coe=Coe)
    nRL = nRadiativeLine
    del Coe, nRadiativeLine
    
    # make mesh
    #--------------------
    waveMesh = _WavelengthMesh.init_Wave_Mesh_(Cont, Line, RL.Coe)

    # read PI
    #--------------------
    Cont_mesh : T_ARRAY = waveMesh.Cont_mesh
    alpha_table, alpha_table_idxs, Coe, alpha_interp, data_source_PI = \
        _AtomIO.make_Atom_PI_(path_dict["PI"],Level,Cont,Cont_mesh,_atom_type,Level_info_table,Cont_ctj_table)
    PI = Photo_Ionization(
        alpha_table = alpha_table,
        alpha_table_idxs = alpha_table_idxs,
        Coe = Coe,
        alpha_interp = alpha_interp
    )
    del alpha_table, alpha_table_idxs, Coe, alpha_interp

    # make ATOMIC_DATA_SOURCE
    _atomic_data_source = ATOMIC_DATA_SOURCE(
        AJI = data_source_Aji,
        CE  = data_source_CE,
        CI  = data_source_CI,
        PI  = data_source_PI
    )

    atom = Atom(
        Z  = Z,
        Mass = Mass,
        Abun = Abun,
        nLevel = nLevel,
        nLine  = nLine,
        nCont  = nCont,
        nTran  = nTran,
        nRL    = nRL,
        Level  = Level,
        Line   = Line,
        Cont   = Cont,
        _has_continuum = _has_continuum,
        _atomic_data_source = _atomic_data_source,
        _atom_type = _atom_type,
        _ctj_table = _ctj_table,
        _idx_table = _idx_table,
        CE = CE,
        CI = CI,
        PI = PI,
        RL = RL,
    )
    return atom, waveMesh , path_dict



