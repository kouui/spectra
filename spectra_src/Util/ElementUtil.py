
#-------------------------------------------------------------------------------
# function definition of Element Util
#-------------------------------------------------------------------------------
# VERSION
#
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

from ..ImportAll import *

def _sym_to_idx_( ele : T_STR ) -> T_INT:

    idx = ELEMENT_SYMBOL.index( ele )
    if idx < 0:
        raise ValueError(f"Can't find element : {ele}")

    return idx

def sym_to_z_( ele : T_STR ):

    idx = _sym_to_idx_( ele )
    return ELEMENT_Z[idx]

def sym_to_mass_( ele : T_STR ):

    idx = _sym_to_idx_( ele )
    return ELEMENT_MASS[idx]

def sym_to_abun_( ele : T_STR ):

    idx = _sym_to_idx_( ele )
    return ELEMENT_ABUN[idx]

def sym_to_ioniz_potential_( ele : T_STR ):

    idx = _sym_to_idx_( ele )
    return ELEMENT_IONIZPOTENTIAL[idx]