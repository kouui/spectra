
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

    ion        : T_STR          #  element, ex. 'FeI', 'H_I'
    name       : T_STR          #  name (Ha) etc.
    mode       : T_STR          #  profile mode ('HF','v1', 'Gaus' etc)

    upper      : T_STR          #  upper state, ex. '3L0.5'
    lower      : T_STR          #  lower state, ex. '3L0.5'


    wl0        : T_FLOAT        #  central wavelength  [cm]
    ep         : T_FLOAT        #  excitation potential, [erg]
    gf         : T_FLOAT        #  gf value
    gamma      : T_FLOAT        #  dumping width in unit of. Dopp./w
    geff       : T_FLOAT        #  effective Lande-g
    G2m        : T_FLOAT        #  second order moment of Zeeman comp.
    ew         : T_FLOAT        #  equivalent width (mA)


