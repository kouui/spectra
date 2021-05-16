
#-------------------------------------------------------------------------------
# Elements definition in spectra
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

from .Types import *

# Z : atomic number
# Mass : relative atom mass
# Abundance : relative abundance, n_ele = 10**(Abundance-12) * n_H
ELEMENT_DICT : T_DICT[T_STR, T_DICT] = \
{
    "He" : {"Z":2, "Mass":4.003, "Abundance":11.00},
    "C"  : {"Z":6, "Mass":12.01, "Abundance":8.54},
    "N"  : {"Z":7, "Mass":14.01, "Abundance":8.06},
    "O"  : {"Z":8, "Mass":16.00, "Abundance":8.83},
    "Ne" : {"Z":10, "Mass":20.18, "Abundance":7.55},
    "Na" : {"Z":11, "Mass":23.00, "Abundance":6.45},
    "Mg" : {"Z":12, "Mass":24.32, "Abundance":7.54},
    "Al" : {"Z":13, "Mass":26.97, "Abundance":6.45},
    "Si" : {"Z":14, "Mass":28.09, "Abundance":7.65},
    "P"  : {"Z":15, "Mass":30.97, "Abundance":5.45},
    "S"  : {"Z":16, "Mass":32.07, "Abundance":7.21},
    "Ar" : {"Z":18, "Mass":39.94, "Abundance":6.75},
    "K"  : {"Z":19, "Mass":39.10, "Abundance":5.70},
    "Ca" : {"Z":20, "Mass":40.08, "Abundance":6.40},
    "Cr" : {"Z":24, "Mass":52.01, "Abundance":6.00},
    "Mn" : {"Z":25, "Mass":54.93, "Abundance":5.55},
    "Fe" : {"Z":26, "Mass":55.85, "Abundance":7.72},
    "Co" : {"Z":27, "Mass":58.94, "Abundance":5.35},
    "Ni" : {"Z":28, "Mass":58.69, "Abundance":6.40},
    "H"  : {"Z":1 , "Mass":1.008, "Abundance":12.00},
}


#-- element table
_ELEMENT   : T_TUPLE[T_STR, ...] = \
    ('He','C','N','O','Ne','Na','Mg','Al','Si','P', 'S','Ar','K','Ca','Cr','Mn','Fe','Co','Ni','H')

#-- atomic number
_NZ        : T_TUPLE[T_INT, ...] = \
    (2, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20, 24, 25, 26, 27, 28, 1)

#-- relative atom mass, where 'C' has value of 12
_AM        : T_TUPLE[T_FLOAT, ...] = \
    (4.003, 12.01, 14.01, 16.00, 20.18, 23.00, 24.32, 26.97, 28.09, 30.97,
     32.07, 39.94, 39.10, 40.08, 52.01, 54.93, 55.85, 58.94, 58.69, 1.008)

#--  relative abundance
_ABUNDANCE : T_TUPLE[T_FLOAT, ...] = \
    (11.00, 8.54, 8.06, 8.83, 7.55, 6.45, 7.54, 6.45, 7.65, 5.45, 
    7.21, 6.75, 5.70, 6.40, 6.00, 5.55, 7.72, 5.35, 6.40, 12.00)