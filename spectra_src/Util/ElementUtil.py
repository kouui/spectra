
#-------------------------------------------------------------------------------
# function definition of Element Util
#-------------------------------------------------------------------------------
# VERSION
#
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

from ..ImportAll import *

from . import RomanUtil as _RomanUtil

def _sym_to_idx_( ele : T_STR ) -> T_INT:

    idx = ELEMENT_SYMBOL.index( ele )
    if idx < 0:
        raise ValueError(f"Can't find element : {ele}")

    return idx

def sym_to_z_( ele : T_STR ) -> T_INT:

    idx = _sym_to_idx_( ele )
    return ELEMENT_Z[idx]

def sym_to_mass_( ele : T_STR ) -> T_FLOAT:

    idx = _sym_to_idx_( ele )
    return ELEMENT_MASS[idx]

def sym_to_abun_( ele : T_STR ) -> T_FLOAT:

    idx = _sym_to_idx_( ele )
    return ELEMENT_ABUN[idx]

#def sym_to_ioniz_potential_( ele : T_STR ):
#
#    idx = _sym_to_idx_( ele )
#    return ELEMENT_IONIZPOTENTIAL[idx]


def format_sym_(sym : T_STR) -> T_STR:
    return sym.lower().capitalize()

def format_stage_(stage : T_STR) -> T_STR:
    return stage.upper()

def format_ion_( ion : T_STR ) -> T_STR:

    sym, stage = ion.split('_')
    #return format_sym_(sym) + '_' + format_stage_(stage)
    return sym_and_stage_to_ion_( sym, stage )

def ion_to_sym_and_stage_( ion : T_STR ) -> T_TUPLE[T_STR, T_STR]:

    sym, stage = ion.split('_')
    sym = format_sym_( sym )
    stage = format_stage_( stage )

    return sym, stage

def sym_and_stage_to_ion_( sym : T_STR, stage : T_STR ) -> T_STR:

    return format_sym_(sym) + '_' + format_stage_(stage)

def ion_to_ioniz_potential_( ion : T_STR ) -> T_FLOAT:

    sym, stage = ion_to_sym_and_stage_(ion)
    ip_arr = ELEMENT_IONIZPOTENTIAL[ _sym_to_idx_(sym) ]
    stage_int = _RomanUtil.roman_to_index_( stage )
    
    return ip_arr[ stage_int - 1 ] * CST.eV2erg_

def shfit_ion_(ion : T_STR, k : T_INT ) -> T_STR:

    sym, stage = ion_to_sym_and_stage_(ion)
    stagep1 = _RomanUtil.shift_roman_( stage, k )
    
    return sym_and_stage_to_ion_( sym, stagep1 )


