
#-------------------------------------------------------------------------------
# definition of functions for creating Atom struct
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

from ...ImportAll import *

import numpy as _numpy
import os

from collections import OrderedDict as _OrderedDict

#-------------------------------------------------------------------------------
# basic function to manipulate data files
#-------------------------------------------------------------------------------

def skip_line_(ln : T_STR) -> T_BOOL:
    """
    skip (1). comment line. start with '#' 
         (2). empty line. ln.strip() is ''
    """
    if ln[0] == "#" or ln.strip() == '':
        return True
    return False

def check_end_(ln : T_STR) -> T_BOOL:
    """
    end of data reading
    (1). line starts with "END" end data reading
    """
    if ln[:3].upper() == "END":
        return True
    return False

def read_general_info_(rs : T_INT, 
    lns : T_LIST[T_STR]) -> T_TUPLE[T_INT, T_STR, T_INT, T_STR, T_INT]:
    """read general information in *.Level config file
        1. Title
        2. Z
        3. Element
        4. nLevel
    """
    for i, ln in enumerate(lns[rs:]):

        if skip_line_(ln):
            continue
        if check_end_(ln):
            break

        if ln.split()[0].strip().lower() == "title:":
            Title = ' '.join( ln.split()[1:] )

        if ln.split()[0].strip().lower() == "z":
            Z = int( ln.split()[-1].strip() )

        if ln.split()[0].strip().lower() == "element":
            Element = ln.split()[-1].strip()

        if ln.split()[0].strip().lower() == "nlevel":
            nLevel = int( ln.split()[-1].strip() )

    re = rs + i + 1

    return re, Title, Z, Element, nLevel

def read_Level_info_(rs : T_INT, lns : T_LIST[T_STR], 
                     Level_info : T_DICT[T_STR,T_LIST], 
                     erg : T_ARRAY, g : T_ARRAY, stage : T_ARRAY, 
                     n : T_ARRAY) -> T_INT:
    """read level information to
        1. Level_info
            - ["configuration"]
            - ["term"]
            - ["J"]
            - ["2S+1"]
        2. erg
        3. g
        4. stage
        5. n
    """
    idx = 0
    #_prefix = ''
    for i, ln in enumerate(lns[rs:]):

        if skip_line_(ln):
            continue
        elif check_end_(ln):
            break

        words = ln.split()
        words = [v.strip() for v in words]

        if words[0] == "prefix":
            prefix = words[1] if words[1] != '-' else ''
            continue

        if words[0] == '-' :
            if prefix == '':
                Level_info["configuration"].append( '-' )
            elif prefix[-1] == '.':
                Level_info["configuration"].append( prefix[:-1] )
            else:
                assert False
        else:
            Level_info["configuration"].append( prefix+words[0] )

        Level_info["term"].append( words[1] )
        Level_info["J"].append( words[2] )
        Level_info["2S+1"].append( words[5] )

        n[idx] = int(words[3]) if words[3] != '-' else 0
        erg[idx] = float( words[8] )
        g[idx] = int( words[6] )
        stage[idx] = int( words[7] )

        idx += 1

    re = rs + i + 1

    return re

def read_Line_info_(lns : T_LIST[T_STR], Aji : T_ARRAY, 
                    line_ctj_table : T_CTJ_PAIR_TABLE):
    """read Line information
    """
    #_count = 0
    #_prefix = ''
    for i, ln in enumerate(lns[:]):

        if skip_line_(ln):
            continue
        elif check_end_(ln):
            break

        words = ln.split()
        words = [v.strip() for v in words]

        if words[0] == "prefix":
            prefix = words[1] if words[1] != '-' else ''
            continue

        # get ctj pair
        ctj_ij = ( (prefix+words[0],words[1],words[2]), (prefix+words[3],words[4],words[5]) )

        if ctj_ij in line_ctj_table:
            line_index = line_ctj_table.index( ctj_ij )
            Aji[line_index] += float( words[6] )
            #_count += 1

def read_CE_temperature_(lns : T_LIST[T_STR]) -> T_TUPLE[T_INT,T_INT,T_LIST[T_FLOAT],T_STR]:
    """read Temperature grid for interpolation
    """
    for i, ln in enumerate(lns[:]):

        if skip_line_(ln):
            continue
        elif check_end_(ln):
            break

        words = ln.split()
        words = [v.strip() for v in words]

        if words[0].lower() == "type":
            CE_type = words[1]

        if words[0].lower() == "temperature":
            Te = [float(v) for v in words[1:]]

    nTe = len( Te )
    re = i + 1

    return re, nTe, Te, CE_type

def read_CE_table_(rs : T_INT, lns : T_LIST[T_STR], CE_table : T_ARRAY, 
                   idxI : T_ARRAY, idxJ : T_ARRAY, f1 : T_ARRAY, f2 : T_ARRAY, 
                   level_info_table : T_CTJ_TABLE, line_ctj_table : T_CTJ_PAIR_TABLE):
    """read CE table for interpolation
    """
    #_count = 0
    #_prefix = ''
    for i, ln in enumerate(lns[rs:]):

        if skip_line_(ln):
            continue
        elif check_end_(ln):
            break

        words = ln.split()
        words = [v.strip() for v in words]

        if words[0] == "prefix":
            prefix = words[1] if words[1] != '-' else ''
            continue

        # get ctj pair
        ctj_ij = ( (prefix+words[0],words[1],words[2]), (prefix+words[3],words[4],words[5]) )

        if ctj_ij in line_ctj_table:
            line_index = line_ctj_table.index( ctj_ij )
            idxI[line_index] = level_info_table.index( ctj_ij[0] )
            idxJ[line_index] = level_info_table.index( ctj_ij[1] )
            CE_table[line_index,:] += [float(v) for v in words[6:-2]]
            f1[line_index] = float(words[-2])
            f2[line_index] = float(words[-1])

            #_count += 1

def read_CI_temperature_(lns : T_LIST[T_STR]) -> T_TUPLE[T_INT,T_INT,T_LIST[T_FLOAT]]:
    r"""
    read Temperature grid for interpolation
    """
    for i, ln in enumerate(lns[:]):

        if skip_line_(ln):
            continue
        elif check_end_(ln):
            break

        words = ln.split()
        words = [v.strip() for v in words]

        #if _words[0].lower() == "ncont":
        #    _nCont = int( _words[1] )

        if words[0].lower() == "temperature":
            Te = [float(_v) for _v in words[1:]]

    nTe = len( Te )
    re = i + 1

    return re, nTe, Te

def read_CI_table_(rs : T_INT, lns : T_LIST[T_STR], CI_table : T_ARRAY, 
                   f2 : T_ARRAY, idxI : T_ARRAY, idxJ : T_ARRAY, 
                   level_info_table : T_CTJ_TABLE, cont_ctj_table : T_CTJ_PAIR_TABLE):
    """read CI table for interpolation
    """
    #_count = 0
    #_prefix = ''
    for i, ln in enumerate(lns[rs:]):

        if skip_line_(ln):
            continue
        elif check_end_(ln):
            break

        words = ln.split()
        words = [v.strip() for v in words]

        if words[0] == "prefix":
            prefix = words[1] if words[1] != '-' else ''
            continue

        # get ctj pair
        if words[3] == '-':
            if prefix == '':
                ctj_ij = ( (prefix+words[0],words[1],words[2]), (words[3],words[4],words[5]) )
            elif prefix[-1] == '.':
                ctj_ij = ( (prefix+words[0],words[1],words[2]), (prefix[:-1],words[4],words[5]) )
            else:
                assert False
        else:
            ctj_ij = ( (prefix+words[0],words[1],words[2]), (prefix+words[3],words[4],words[5]) )

        if ctj_ij in cont_ctj_table:
            contIndex = cont_ctj_table.index( ctj_ij )
            idxI[contIndex] = level_info_table.index( ctj_ij[0] )
            idxJ[contIndex] = level_info_table.index( ctj_ij[1] )
            CI_table[contIndex,:] += [float(v) for v in words[6:-1]]
            f2[contIndex] = float(words[-1])

        #_count += 1

def read_PI_info_(lns : T_LIST[T_STR]) -> T_TUPLE[T_INT,T_INT]:
    """read nCont and nMesh for Photoionization
    """
    for i, ln in enumerate(lns[:]):

        if skip_line_(ln):
            continue
        elif check_end_(ln):
            break

        words = ln.split()
        words = [v.strip() for v in words]

        if words[0].lower() == "ncont":
            nCont = int( words[1] )

        #if _words[0].lower() == "nmesh":
        #    _nMesh = int( _words[1] )

    re = i + 1

    return nCont, re#, _nMesh

def read_PI_table_(rs : T_INT, lns : T_LIST[T_STR], PI_table_dict : T_DICT[T_INT,T_ARRAY], 
                   PI_coe : T_ARRAY, level_info_table : T_CTJ_TABLE, cont_ctj_table : T_CTJ_PAIR_TABLE):
    """read PI table for interpolation
    """
    countMesh = 0
    readMesh = False

    #_prefix = ''
    for i, ln in enumerate(lns[rs:]):

        if skip_line_(ln):
            continue
        elif check_end_(ln):
            break

        words = ln.split()
        words = [v.strip() for v in words]

        if words[0] == "prefix":
            prefix = words[1] if words[1] != '-' else ''
            continue

        # get ctj pair
        if len(words) > 2:

            if words[3] == '-':
                if prefix == '':
                    ctj_ij = ( (prefix+words[0],words[1],words[2]), (words[3],words[4],words[5]) )
                elif prefix[-1] == '.':
                    ctj_ij = ( (prefix+words[0],words[1],words[2]), (prefix[:-1],words[4],words[5]) )
                else:
                    assert False
            else:
                ctj_ij = ( (prefix+words[0],words[1],words[2]), (prefix+words[3],words[4],words[5]) )

            if ctj_ij in cont_ctj_table:
                contIndex = cont_ctj_table.index( ctj_ij )
                PI_coe["idxI"][contIndex] = level_info_table.index( ctj_ij[0] )
                PI_coe["idxJ"][contIndex] = level_info_table.index( ctj_ij[1] )

                nLambda = int(words[6])
                PI_coe["nLambda"][contIndex] = nLambda
                PI_coe["alpha0"][contIndex] = float(words[8])

                readMesh = True
                countMesh = 0
                mesh_array = _numpy.zeros((2,nLambda),T_FLOAT)
            else:
                readMesh = False

        else:
            if not readMesh:
                continue
            #_PI_table[_countLine,_countMesh,:] += [float(v) for v in _words]
            mesh_array[:,countMesh] = [float(v) for v in words]

            if countMesh == (nLambda-1):
                PI_table_dict[contIndex] = mesh_array

            countMesh += 1

def read_Radiative_Line_number_(lns : T_LIST[T_STR]) -> T_TUPLE[T_INT,T_INT]:
    """read nRadiativeLine
    """
    for i, ln in enumerate(lns[:]):

        if skip_line_(ln):
            continue
        elif check_end_(ln):
            break

        words = ln.split()
        words = [v.strip() for v in words]

        if words[0].lower() == "nradiative":
            nRadiativeLine = int( words[1] )

    re = i + 1

    return nRadiativeLine, re

def read_Mesh_info_(rs : T_INT, lns : T_LIST[T_STR], Mesh_coe : T_ARRAY, 
                    filename : T_LIST[T_STR], level_info_table : T_CTJ_TABLE, line_ctj_table : T_CTJ_PAIR_TABLE):
    r"""read CI table for interpolation
    """
    count = 0
    #_prefix = ''
    for i, ln in enumerate(lns[rs:]):

        if skip_line_(ln):
            continue
        elif check_end_(ln):
            break

        words = ln.split()
        words = [v.strip() for v in words]

        if words[0] == "prefix":
            prefix = words[1] if words[1] != '-' else ''
            continue

        # get ctj pair
        ctj_ij = ( (prefix+words[0],words[1],words[2]), (prefix+words[3],words[4],words[5]) )

        if ctj_ij not in line_ctj_table:
            raise ValueError(f"cannot find {ctj_ij} in line_ctj_table")
        #assert ctj_ij in line_ctj_table

        Mesh_coe["lineIndex"][count] = line_ctj_table.index( ctj_ij )
        Mesh_coe["idxI"][count] = level_info_table.index( ctj_ij[0] )
        Mesh_coe["idxJ"][count] = level_info_table.index( ctj_ij[1] )

        if words[6] == "Voigt":
            Mesh_coe["ProfileType"][count] = int( E_ABSORPTION_PROFILE_TYPE.VOIGT )
        elif words[6] == "Gaussian":
            Mesh_coe["ProfileType"][count] = int( E_ABSORPTION_PROFILE_TYPE.GAUSSIAN )
        else:
            raise ValueError("Profile type should be either 'Voigt' or 'Gaussian'")

        Mesh_coe["nLambda"][count] = int(words[7])
        Mesh_coe["qcore"][count] = float(words[10])
        Mesh_coe["qwing"][count] = float(words[11])

        filename.append( words[12] )

        count += 1

def read_conf_(conf_path : T_STR) -> T_DICT[T_STR, T_UNION[None,T_STR]]:

    if not conf_path.endswith(".conf"):
        raise ValueError("`conf_path` should be a string ends with '.conf'")

    path_dict_keys = ("folder","conf","Level","Aji","CEe","CIe","PI","RadiativeLine","Grotrian")
    path_dict : T_DICT[T_STR, T_UNION[None,T_STR]] = {
        key : None for key in path_dict_keys
    }

    i = conf_path.rfind('/')
    folder = conf_path[:i+1]
    path_dict["conf"] = os.path.abspath( conf_path )

    with open(conf_path, "r") as f:
        lines = f.readlines()

    for ln in lines:
        words =ln.split()
        words = [v.strip() for v in words]

        try:
            _ = path_dict[words[0]]
        except KeyError:
            raise ValueError(f"{words[0]} is not a valid key")

        if words[0] == "folder":
            path_dict[words[0]] = os.path.abspath( os.path.join( folder, words[1] ) )
        else:
            if path_dict["folder"] is None:
                raise ValueError("`folder` configuration should be read first.")
            path_dict[words[0]] = os.path.abspath( os.path.join( path_dict["folder"], words[1] ) )

    return path_dict


#-------------------------------------------------------------------------------
# functions for constructing structs
#-------------------------------------------------------------------------------

def make_Atom_Level_(path : T_STR) -> T_TUPLE[T_INT,T_FLOAT,T_FLOAT,T_INT,T_ARRAY,T_TUPLE[T_TUPLE[T_STR,T_STR,T_STR],...]]:

    with open(path, 'r') as f:
        fLines = f.readlines()
    #--- read general info
    rs, title, Z, element, nLevel = read_general_info_(rs=0, lns=fLines)
    Mass : T_FLOAT = ELEMENT_DICT[element]["Mass"]
    Abun : T_FLOAT = 10**(ELEMENT_DICT[element]["Abundance"]-12.0)
    #--- read Level info
    dtype  = _numpy.dtype([
                          ('erg',T_FLOAT),            #: level energy, erg
                          ('g',T_INT),                #: g=2J+1, statistical weight
                          ('stage',T_INT),            #: ionization stage
                          ('gamma',T_FLOAT),          #: radiative damping constant of Level
                          ("isGround",T_BOOL),        #: whether a level is ground level
                          ('n',T_INT),                #: quantum number n
                          ])
    Level = _numpy.zeros(nLevel, dtype=dtype)
    Level_info : T_DICT[T_STR, T_LIST[T_STR]] = {"configuration" : [], "term" : [], "J": [], "2S+1": []}
    rs = read_Level_info_(rs, lns=fLines, Level_info=Level_info,
                                  erg=Level["erg"][:], g=Level["g"][:],
                                  stage=Level["stage"][:], n=Level["n"][:])
    Level["erg"][:] *= CST.eV2erg_

    Level["isGround"][:] = 1
    for k in range(1,nLevel):
        if Level["stage"][k] == Level["stage"][k-1]:
            Level["isGround"][k] = 0

    Level["gamma"][:] = 0

    #--- make tuple of tuple (configuration, term, J)
    _Level_info_table : T_LIST[T_TUPLE[T_STR,T_STR,T_STR]] = []
    for k in range(nLevel):
        _Level_info_table.append((Level_info["configuration"][k],
                                  Level_info["term"][k],
                                  Level_info["J"][k]))
    Level_info_table : T_TUPLE[T_TUPLE[T_STR,T_STR,T_STR],...] = tuple(_Level_info_table)

    return Z, Mass, Abun, nLevel, Level, Level_info_table

def nLine_nCont_nTran_(stage : T_ARRAY) -> T_TUPLE[T_INT,T_INT,T_INT,T_BOOL] :

    #stage = Level["stage"]
    nLevel = stage.size
    count = 0
    current_stage = stage[0]
    
    counter = _OrderedDict()
    for k,s in enumerate(stage):
        if s == current_stage:
            count += 1
            if k == nLevel-1:
                counter[current_stage] = count
                break
        elif s > current_stage:
            counter[current_stage] = count
            count = 1
            current_stage = s
            counter[current_stage] = count

    nLine, nCont = 0, 0
    for k, v in counter.items():
        nLine += v * (v-1) // 2
        if k+1 in counter.keys():
            nCont += v
    nTran = nLine + nCont

    has_continuum = True if nCont > 0 else False

    return nLine, nCont, nTran, has_continuum

def prepare_idx_ctj_mapping_(Level_info_table : T_CTJ_TABLE, stage : T_ARRAY, isGround : T_ARRAY, 
    nLine : T_INT, nCont : T_INT) -> T_TUPLE[ T_IDX_PAIR_TABLE,T_CTJ_PAIR_TABLE,T_IDX_PAIR_TABLE,T_CTJ_PAIR_TABLE ]:
    r"""make tuples for mapping
    - lineIndex <--> (ctj_i, ctj_j)
    - lineIndex <--> (idxI, idxJ)
    - contIndex <--> (ctj_i, ctj_j)
    - contIndex <--> (idxI, idxJ)
    """

    _T_INFO_DICT = T_DICT[ T_STR, T_LIST[T_ANY] ]

    def _idx_ctj_into_dict(_dict : _T_INFO_DICT, _i : T_INT, _j : T_INT):
        _dict["idxI"].append( _i )
        _dict["idxJ"].append( _j )
        _dict["ctj_i"].append( Level_info_table[_i] )
        _dict["ctj_j"].append( Level_info_table[_j] )

    nLevel = stage.size

    Tran_dict : _T_INFO_DICT = {}
    Line_dict : _T_INFO_DICT = {}
    Cont_dict : _T_INFO_DICT = {}
    for key in ("idxI", "idxJ", "ctj_i", "ctj_j"):
        Tran_dict[key] = []
        Line_dict[key] = []
        Cont_dict[key] = []

    for i in range(0, nLevel):
        for j in range(i+1, nLevel):
            # i : lower level
            # j : upper level
            if stage[i] == stage[j]:
                _idx_ctj_into_dict(Tran_dict, i, j)
                _idx_ctj_into_dict(Line_dict, i, j)
            elif stage[i] == stage[j]-1 and isGround[j]:
                _idx_ctj_into_dict(Tran_dict, i, j)
                _idx_ctj_into_dict(Cont_dict, i, j)

    if nLine != len(Line_dict["idxI"]):
        raise ValueError("incompatible size of Line_dict['idxI']")
    if nCont != len(Cont_dict["idxI"]):
        raise ValueError("incompatible size of Cont_dict['idxI']")

    Line_idx_table1 : T_LIST[T_TUPLE[T_INT,T_INT]] = []
    Line_ctj_table1 : T_LIST[T_CTJ_PAIR] = []
    for k in range(nLine):
        Line_idx_table1.append( ( Line_dict["idxI"][k] , Line_dict["idxJ"][k]  ) )
        Line_ctj_table1.append( ( Line_dict["ctj_i"][k], Line_dict["ctj_j"][k] ) )
    Line_idx_table = tuple( Line_idx_table1 )
    Line_ctj_table = tuple( Line_ctj_table1 )

    Cont_idx_table1 : T_LIST[T_TUPLE[T_INT,T_INT]] = []
    Cont_ctj_table1 : T_LIST[T_CTJ_PAIR] = []
    for k in range(nCont):
        Cont_idx_table1.append( ( Cont_dict["idxI"][k] , Cont_dict["idxJ"][k]  ) )
        Cont_ctj_table1.append( ( Cont_dict["ctj_i"][k], Cont_dict["ctj_j"][k] ) )
    Cont_idx_table = tuple( Cont_idx_table1 )
    Cont_ctj_table = tuple( Cont_ctj_table1 )

    return Line_idx_table, Line_ctj_table, Cont_idx_table, Cont_ctj_table

def make_Atom_Cont_(nCont : T_INT, Cont_idx_table : T_IDX_PAIR_TABLE, Level : T_ARRAY) -> T_ARRAY:

    dtype = _numpy.dtype([('idxI',T_INT),            #: level index, the Level index of lower level
                          ('idxJ',T_INT),            #: level index, the Level index of lower level
                          ('f0',T_FLOAT),            #: central frequency
                          ('w0',T_FLOAT),            #: central wavelength in cm
                          ('w0_AA',T_FLOAT),         #: central wavelength in Angstrom
                          ('gi',T_INT),              #: statistical weight of lower level
                          ('gj',T_INT),              #: statistical weight of upper level
                          ('ni',T_INT),              #: quantum number n of lower level
                          ('nj',T_INT),              #: quantum number n of upper level
                          ])
    Cont = _numpy.zeros(nCont, dtype=dtype)
    
    if nCont > 0:
        for k in range(nCont):
            i, j = Cont_idx_table[k]
            Cont["idxI"][k], Cont["idxJ"][k] = i, j
            Cont["f0"][k] = (Level["erg"][j]-Level["erg"][i]) / CST.h_
        Cont["w0"][:] = CST.c_ / Cont["f0"][:]
        Cont["w0_AA"][:] = Cont["w0"][:] * 1E+8

        # read gi,gj,ni,nj
        Cont["gi"][:] = Level["g"][ Cont["idxI"][:] ]
        Cont["gj"][:] = Level["g"][ Cont["idxJ"][:] ]
        Cont["ni"][:] = Level["n"][ Cont["idxI"][:] ]
        Cont["nj"][:] = Level["n"][ Cont["idxJ"][:] ]

    return Cont

def make_Atom_Line_(path : T_UNION[T_STR,None], Level : T_ARRAY, 
                    Line_idx_table : T_IDX_PAIR_TABLE, Line_ctj_table : T_CTJ_PAIR_TABLE,
                    atom_type : T_E_ATOM) -> T_TUPLE[T_ARRAY,T_E_ATOMIC_DATA_SOURCE]:

    dtype = _numpy.dtype([('idxI',T_INT),            #: level index, the Level index of lower level
                          ('idxJ',T_INT),            #: level index, the Level index of lower level
                          ('AJI',T_FLOAT),           #: Einstein Aji coefficient
                          ('f0',T_FLOAT),            #: central frequency
                          ('w0',T_FLOAT),            #: central wavelength in cm
                          ('w0_AA',T_FLOAT),         #: central wavelength in Angstrom
                          #("isContinuum",T_INT),    #: continuum tansition identifier, 0: same stage, 1: continuum transition, 2: others
                          ('Gamma',T_FLOAT),         #: radiative damping constant of Line
                          ('gi',T_INT),              #: statistical weight of lower level
                          ('gj',T_INT),              #: statistical weight of upper level
                          ('ni',T_INT),              #: quantum number n of lower level
                          ('nj',T_INT),              #: quantum number n of upper level
                          ('BJI',T_FLOAT),           #: Einstein Bji coefficient
                          ('BIJ',T_FLOAT),           #: Einstein BIJ coefficient
                          ])
    nLine = len( Line_idx_table )
    Line = _numpy.recarray(nLine, dtype=dtype)

    # idxI and idxJ
    for k in range(nLine):
        Line["idxI"][k], Line["idxJ"][k] = Line_idx_table[k]

    Line["AJI"][:] = 0

    data_source_Aji : T_E_ATOMIC_DATA_SOURCE
    if path is not None: # normal case
        with open(path, 'r') as f:
            fLines = f.readlines()
        read_Line_info_(fLines, Line["AJI"][:], Line_ctj_table)
        data_source_Aji = E_ATOMIC_DATA_SOURCE.EXPERIMENT

    else: # in case of we want to calculate Aji numerically
        data_source_Aji = E_ATOMIC_DATA_SOURCE.CALCULATE
        ## hydrogen has function prepared
        if atom_type == E_ATOM.HYDROGEN:
            from ...Atomic import Hydrogen as _Hydrogen
            Line["AJI"][:] = _Hydrogen.einstein_A_coefficient_(Line['ni'][:], Line['nj'][:])
        else:
            raise ValueError("We don't have function to calculate Aji for non-hydrogen atoms")

    # calculate f0, w0, w0_AA
    for k in range(nLine):
        i : int = Line["idxI"][k]
        j : int = Line["idxJ"][k]
        Line["f0"][k] = (Level["erg"][j]-Level["erg"][i]) / CST.h_
    Line["w0"][:] = CST.c_ / Line["f0"][:]
    Line["w0_AA"][:] = Line["w0"][:] * 1.E+8

    # read gi,gj,ni,nj
    Line["gi"][:] = Level["g"][ Line["idxI"][:] ]
    Line["gj"][:] = Level["g"][ Line["idxJ"][:] ]
    Line["ni"][:] = Level["n"][ Line["idxI"][:] ]
    Line["nj"][:] = Level["n"][ Line["idxJ"][:] ]


    # compute Bji, Bij
    from ...Atomic import LTELib as _LTELib
    Line["BJI"][:], Line["BIJ"][:] = _LTELib.einsteinA_to_einsteinBs_cm_(Line["AJI"][:], Line["w0"][:], Line["gi"][:], Line["gj"][:])

    # compute Level gamma and Line Gamma
    Line["Gamma"][:] = 0.
    from ...Atomic import BasicP as _BasicP
    _BasicP.update_level_gamma_(Line["AJI"][:],Line["idxJ"][:],Level["gamma"][:])
    _BasicP.update_line_gamma_(Line["idxI"][:],Line["idxJ"][:],Level["gamma"][:],Line["Gamma"][:])


    return Line, data_source_Aji

def make_Atom_CECI_(path_electron : T_UNION[T_STR, None], tran_type : T_STR, 
               n_transition : T_INT, Tran : T_ARRAY, Level : T_ARRAY,
               Level_info_table : T_CTJ_TABLE, Tran_ctj_table : T_CTJ_PAIR_TABLE
              ) -> T_TUPLE[T_ARRAY,T_ARRAY,T_ARRAY,T_E_COLLISIONAL_TRANSITION,T_E_COLLISIONAL_TRANSITION_SOURCE,T_E_COLLISIONAL_TRANSITION_FORMULA]:

    ## : transition type
    _transition_type : T_E_COLLISIONAL_TRANSITION
    if tran_type == "CE":
        _transition_type     = E_COLLISIONAL_TRANSITION.EXCITATION
    elif tran_type == "CI":
        _transition_type     = E_COLLISIONAL_TRANSITION.IONIZATION
    else:
        raise ValueError("only 'CE' and 'CI' are supported for `transition_type` argument.")

    ## : transition source
    ## currently only electron collisional transition is available
    _transition_source  : T_E_COLLISIONAL_TRANSITION_SOURCE = E_COLLISIONAL_TRANSITION_SOURCE.ELECTRON

    ## : transition formula
    ## currently only formula OMEGA is available
    _transition_formula : T_E_COLLISIONAL_TRANSITION_FORMULA = E_COLLISIONAL_TRANSITION_FORMULA.OMEGA
    
    nTe : T_INT
    dtype  = _numpy.dtype([
                          ('idxI',T_INT),      #: level index, the Level index of lower level
                          ('idxJ',T_INT),      #: level index, the Level index of upper level
                          ('f1',T_INT),        #: a factor for ESC calculation due to fine structure, \Omega * f1 / f2
                          ('f2',T_INT),        #: a factor for ESC calculation due to fine structure, \Omega * f1 / f2
                          ('gi',T_INT),        #: statistical weight of lower level
                          ('gj',T_INT),        #: statistical weight of upper level
                          ('dEij',T_FLOAT)     #: excitation energy, [:math:`erg`]
                          ])
    
    if path_electron is None:
        
        nTe = 0        
        Te_table    = _numpy.zeros(nTe, dtype=DT_NB_FLOAT)
        Omega_table = _numpy.zeros((n_transition, nTe), dtype=DT_NB_FLOAT)
        Coe         = _numpy.zeros(n_transition, dtype=dtype)
        Coe["f1"][:] = 1
        Coe["f2"][:] = 1
        Coe["idxI"][:] = Tran["idxI"][:]
        Coe["idxJ"][:] = Tran["idxJ"][:]
    
    else:
        with open(path_electron, 'r') as file:
            fLines = file.readlines()

        # read Temperature grid for interpolation
        if tran_type == "CE":
            rs, nTe, Te, CE_type = read_CE_temperature_(fLines)
        elif tran_type == "CI":
            rs, nTe, Te = read_CI_temperature_(fLines)

        Te_table    = _numpy.array(Te, dtype=DT_NB_FLOAT)
        Omega_table = _numpy.zeros((n_transition, nTe), dtype=DT_NB_FLOAT)
        Coe         = _numpy.zeros(n_transition, dtype=dtype)

        if tran_type == "CE":
            read_CE_table_(rs=rs, lns=fLines, CE_table=Omega_table,
                    idxI=Coe["idxI"][:],idxJ=Coe["idxJ"][:],
                    f1=Coe["f1"][:], f2=Coe["f2"][:],
                    level_info_table=Level_info_table,
                    line_ctj_table=Tran_ctj_table)
        
        elif tran_type == "CI":
            read_CI_table_(rs=rs, lns=fLines, CI_table=Omega_table,
                    f2=Coe["f2"][:],
                    idxI=Coe["idxI"][:],idxJ=Coe["idxJ"][:],
                    level_info_table=Level_info_table,
                    cont_ctj_table=Tran_ctj_table)


    for k in range(n_transition):
        Coe["gi"][k] = Level["g"][Coe["idxI"][k]]
        Coe["gj"][k] = Level["g"][Coe["idxJ"][k]]
        Coe["dEij"][k] = Level["erg"][Coe["idxJ"][k]] - Level["erg"][Coe["idxI"][k]]

    return Te_table, Omega_table, Coe, _transition_type, _transition_source, _transition_formula
        
def make_Atom_RL_(path : T_UNION[T_STR, None], 
                  Level_info_table : T_CTJ_TABLE, 
                  Line_ctj_table : T_CTJ_PAIR_TABLE) -> T_TUPLE[T_ARRAY, T_INT]:

    dtype  = _numpy.dtype([
                          ('idxI',T_INT),         #: level index, the Level index of lower level
                          ('idxJ',T_INT),         #: level index, the Level index of upper level
                          ('lineIndex',T_INT),    #: line index
                          ('ProfileType',T_INT),  #: 0: Voigt; 1: Gaussian
                          ('qcore', T_FLOAT),
                          ('qwing', T_FLOAT),
                          ('nLambda', T_INT),     #: number of meaningful wavelength mesh point
                          ])
    
    nRadiativeLine : T_INT
    if path is None:
        nRadiativeLine = 0
        Coe = _numpy.zeros(nRadiativeLine, dtype=dtype)
    else:
        with open(path, 'r') as f:
            fLines = f.readlines()

        nRadiativeLine, rs = read_Radiative_Line_number_(fLines)

        Coe = _numpy.zeros(nRadiativeLine, dtype=dtype)

        RadLine_filenames : T_LIST[T_STR] = []
        read_Mesh_info_(rs=rs, lns=fLines,
                        Mesh_coe = Coe,
                        filename = RadLine_filenames,
                        level_info_table=Level_info_table,
                        line_ctj_table=Line_ctj_table)
    
    return Coe, nRadiativeLine

def make_Atom_PI_(path : T_UNION[T_STR, None], Level : T_ARRAY, Cont : T_ARRAY, Cont_mesh : T_ARRAY, 
        atom_type : T_E_ATOM, Level_info_table : T_CTJ_TABLE, Cont_ctj_table : T_CTJ_PAIR_TABLE
        ) -> T_TUPLE[T_ARRAY,T_ARRAY,T_ARRAY,T_ARRAY, T_E_ATOMIC_DATA_SOURCE]:

    nCont = Cont.shape[0]
    dtype  = _numpy.dtype([
                            ('idxI',T_INT),        #: level index, the Level index of lower level
                            ('idxJ',T_INT),        #: level index, the Level index of upper level
                            #('lineIndex',T_INT),  #: line index
                            ('nLambda', T_INT),    #: number of meaningful wavelength mesh point
                            ('alpha0', T_FLOAT),   #: Photoionization cross section at frequency edge
                            ('gi',T_INT),          #: statistical weight of lower level
                            ('gj',T_INT),          #: statistical weight of upper level
                            ('dEij',T_FLOAT)       #: ionization limit, [:math:`erg`]
                            ])
    Coe = _numpy.empty(nCont, dtype=dtype) # in an order of Cont
    data_source_PI : T_E_ATOMIC_DATA_SOURCE
    alpha_table : T_ARRAY
    alpha_table_idxs : T_ARRAY
    alpha_interp : T_ARRAY

    ## no continuum
    if nCont == 0:
        data_source_PI = E_ATOMIC_DATA_SOURCE.CALCULATE
        alpha_table = _numpy.empty((2,0), dtype=DT_NB_FLOAT)
        alpha_table_idxs = _numpy.empty( (nCont,2), dtype=DT_NB_INT )
        alpha_interp = _numpy.empty( (nCont,0), dtype=DT_NB_INT )
        
        return alpha_table, alpha_table_idxs, Coe, alpha_interp, data_source_PI

    
    if path is None:
        data_source_PI = E_ATOMIC_DATA_SOURCE.CALCULATE
        alpha_table = _numpy.empty((2,0), dtype=DT_NB_FLOAT)
        alpha_table_idxs = _numpy.empty( (nCont,2), dtype=DT_NB_INT )

        Coe["idxI"][:] = Cont["idxI"][:]
        Coe["idxJ"][:] = Cont["idxJ"][:]
        Coe["nLambda"][:] = Cont_mesh.shape[1]

    else:
        data_source_PI = E_ATOMIC_DATA_SOURCE.EXPERIMENT

        with open(path, 'r') as f:
            fLines = f.readlines()

        nCont = Cont.shape[0]
        nCont1, rs = read_PI_info_(fLines)

        if nCont != nCont1:
            raise ValueError( "incompatible nCont from `Cont` and `read_PI_info_`" )

        alpha_table_dict : T_DICT[T_INT,T_ARRAY] = _OrderedDict()
        
    

        read_PI_table_(rs=rs, lns=fLines,
                       PI_table_dict = alpha_table_dict,
                       PI_coe = Coe,
                       level_info_table=Level_info_table,
                       cont_ctj_table=Cont_ctj_table)

        
        sorted_keys = sorted( list(alpha_table_dict.keys()) )
        if len(sorted_keys) != nCont:
            raise ValueError("number of read alpha table != nCont")

        n_total_mesh_length = 0
        for i, key in enumerate(sorted_keys):
            n_total_mesh_length += alpha_table_dict[key].shape[1]

        alpha_table = _numpy.empty((2,n_total_mesh_length), dtype=DT_NB_FLOAT)
        alpha_table_idxs = _numpy.empty( (nCont,2), dtype=DT_NB_INT )
        bias = 0
        for i, key in enumerate(sorted_keys):
            arr2d : T_ARRAY = alpha_table_dict[key]
            alpha_table_idxs[i,0] = bias
            alpha_table_idxs[i,1] = alpha_table_idxs[i,0] + arr2d.shape[1]
            bias = alpha_table_idxs[i,1]

            arr2d[0,:] *= 1.E-7 # unit conversion of wavelength, nm --> cm   (cm, cm^2)
            arr2d[0,:] += Cont["w0"][i] - arr2d[0,0] # shift edge wavelength to the computed wavelength w0

            alpha_table[:, alpha_table_idxs[i,0] : alpha_table_idxs[i,1] ] = arr2d[:,:]

        for k in range(nCont):
            Coe["gi"][k]   = Level["g"][Coe["idxI"][k]]
            Coe["gj"][k]   = Level["g"][Coe["idxJ"][k]]
            Coe["dEij"][k] = Level["erg"][Coe["idxJ"][k]] - Level["erg"][Coe["idxI"][k]]

    
    if data_source_PI == E_ATOMIC_DATA_SOURCE.CALCULATE:
        
        if atom_type != E_ATOM.HYDROGEN:
            raise ValueError("We don't have function to calculate Photoionization cross section for non-hydrogen atom")
        
        from ...Function.Hydrogen import DegenerateN as _DegenerateN
        alpha_interp = _DegenerateN.compute_PI_cross_section_(Cont["ni"][:], Cont_mesh[:,:])

    else:
        from ...Atomic import PhotoIonize as _PhotoIonize
        alpha_interp = _PhotoIonize.interpolate_PI_alpha_(alpha_table[:,:], alpha_table_idxs[:,:], Cont_mesh[:,:])

    return alpha_table, alpha_table_idxs, Coe, alpha_interp, data_source_PI


