
#-------------------------------------------------------------------------------
# definition of functions for creating Atom struct
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

from ...ImportAll import *

import numpy as _numpy

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
                     Level_info : T_DICT[T_STR,T_STR], 
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
                    line_ctj_table : T_CTJ_TABLE):
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
                   level_info_table : T_CTJ_TABLE, line_ctj_table : T_CTJ_TABLE):
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
                   level_info_table : T_CTJ_TABLE, cont_ctj_table : T_CTJ_TABLE):
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

def read_PI_nfo_(lns : T_LIST[T_STR]) -> T_TUPLE[T_INT,T_INT]:
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
                   PI_coe : T_ARRAY, level_info_table : T_CTJ_TABLE, cont_ctj_table : T_CTJ_TABLE):
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
                PI_coe.idxI[contIndex] = level_info_table.index( ctj_ij[0] )
                PI_coe.idxJ[contIndex] = level_info_table.index( ctj_ij[1] )

                nLambda = int(words[6])
                PI_coe.nLambda[contIndex] = nLambda
                PI_coe.alpha0[contIndex] = float(words[8])

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
                    filename : T_LIST[T_STR], level_info_table : T_CTJ_TABLE, line_ctj_table : T_CTJ_TABLE):
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

        Mesh_coe.lineIndex[count] = line_ctj_table.index( ctj_ij )
        Mesh_coe.idxI[count] = level_info_table.index( ctj_ij[0] )
        Mesh_coe.idxJ[count] = level_info_table.index( ctj_ij[1] )

        if words[6] == "Voigt":
            Mesh_coe.ProfileType[count] = E_ABSORPTION_PROFILE_TYPE.VOIGT
        elif words[6] == "Gaussian":
            Mesh_coe.ProfileType[count] = E_ABSORPTION_PROFILE_TYPE.GAUSSIAN
        else:
            raise ValueError("Profile type should be either 'Voigt' or 'Gaussian'")

        Mesh_coe.nLambda[count] = int(words[7])
        Mesh_coe.qcore[count] = float(words[10])
        Mesh_coe.qwing[count] = float(words[11])

        filename.append( words[12] )

        count += 1











