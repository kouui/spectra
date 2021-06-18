
#-------------------------------------------------------------------------------
# definition of functions for Level/Cont/Line indexing/table 
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

from ...ImportAll import *

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