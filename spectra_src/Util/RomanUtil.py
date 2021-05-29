
#-------------------------------------------------------------------------------
# function definition of Roman characters for ionization stage 
#-------------------------------------------------------------------------------
# VERSION
# 0.1.1
#    2021/05/29   u.k.
#        - created this module for ionization stage
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

from ..ImportAll import *

_ROMAN_CHARACTERS = (
    "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
    "XI", "XII", "XIII", "XIV", "XV", "XVI", "XVII", "XVIII", "XIX", "XX",
    "XXI", "XXII", "XXIII", "XXIV", "XXV", "XXVI", "XXVII", "XXVIII", "XXIX", "XXX",
    "XXXI", "XXXII", "XXXIII", "XXXIV", "XXXV", "XXXVI", "XXXVII", "XXXVIII", "XXXIX", "XL",
    "XLI", "XLII", "XLIII", "XLIV", "XLV", "XLVI", "XLVII", "XLVIII", "XLIX", "L",
)

def index_to_roman_(index : T_INT) -> T_STR:

    if index < 1 :
        raise ValueError("index < 1 unavailable")

    return _ROMAN_CHARACTERS[ index - 1 ]

def roman_to_index_(roman : T_STR) -> T_INT:

    index = _ROMAN_CHARACTERS.index( roman )
    
    if index < 0:
        raise ValueError("Cannot find your roman character in our table")

    return index + 1

def shift_roman_(roman : T_STR, k : T_INT) -> T_STR:

    index = roman_to_index_( roman )
    index_new = index + k

    return index_to_roman_( index_new )



    






