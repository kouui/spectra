
#-------------------------------------------------------------------------------
# experimental function/struct for line profiles
#-------------------------------------------------------------------------------
# VERSION
# 0.1.1
#    2021/05/29   k.i., u.k. 
#        - class Line
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

from ..ImportAll import *

from dataclasses import dataclass as _dataclass

@_dataclass(**STRUCT_KWGS_UNFROZEN)
class Level:

    _2Splus1     : T_INT
    _L           : T_INT
    _2Jplus1     : T_INT

@_dataclass(**STRUCT_KWGS_UNFROZEN)
class Line:

    ion        : T_STR          #  element, ex. 'Fe_I', 'H_I'


    wl0        : T_FLOAT        #  central wavelength  [cm]
    ep         : T_FLOAT        #  excitation potential, [erg]
    gf         : T_FLOAT        #  gf value
    gamma      : T_FLOAT        #  dumping width in unit of. Dopp./w
    
    geff       : T_FLOAT = 0.        #  effective Lande-g
    G2m        : T_FLOAT = 0.        #  second order moment of Zeeman comp.
    ew         : T_FLOAT = 0.        #  equivalent width (mA)


    upper      : T_STR = ''          #  upper state, ex. '3L0.5'
    lower      : T_STR = ''          #  lower state, ex. '3L0.5'

    name       : T_STR = ''          #  name (Ha) etc.
    mode       : T_STR = ''          #  profile mode ('HF','v1', 'Gaus' etc)


