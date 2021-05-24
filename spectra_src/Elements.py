
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
_ELEMENT_DICT_OLD : T_DICT[T_STR, T_DICT[T_STR,T_UNION[T_INT,T_FLOAT]]] = \
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

#-------------------------------------------------------------------------------
# more Element data
#-------------------------------------------------------------------------------

from collections import OrderedDict as _OrderedDict
import numpy as _numpy

ELEMENT_DICT : T_DICT = _OrderedDict()
ELEMENT_DICT["H"]    = {"Z":1 , "Mass":1.0080, "Abundance":12.00, "IonizeEV" : [13.508]}
ELEMENT_DICT["He"]   = {"Z":2 , "Mass":4.0026, "Abundance":10.93, "IonizeEV" : [24.587, 54.416]}
ELEMENT_DICT["Li"]   = {"Z":3 , "Mass":6.9410, "Abundance":0.7  , "IonizeEV" : []}
ELEMENT_DICT["Be"]   = {"Z":4 , "Mass":9.0112, "Abundance":1.1  , "IonizeEV" : []}
#ELEMENT_DICT["B"]    = {"Z":5 , "Mass":10.811, "Abundance":?}
ELEMENT_DICT["C"]    = {"Z":6 , "Mass":12.0111, "Abundance":8.52, "IonizeEV" : []}
ELEMENT_DICT["N"]    = {"Z":7 , "Mass":14.0067, "Abundance":7.96, "IonizeEV" : []}
ELEMENT_DICT["O"]    = {"Z":8 , "Mass":15.9994, "Abundance":8.82, "IonizeEV" : []}
ELEMENT_DICT["F"]    = {"Z":9 , "Mass":18.9984, "Abundance":4.6 , "IonizeEV" : []}
ELEMENT_DICT["Ne"]   = {"Z":10, "Mass":20.1790, "Abundance":7.92, "IonizeEV" : []}
ELEMENT_DICT["Na"]   = {"Z":11, "Mass":22.9898, "Abundance":6.25, "IonizeEV" : [5.139, 47.269]}
ELEMENT_DICT["Mg"]   = {"Z":12, "Mass":24.3050, "Abundance":7.42, "IonizeEV" : [7.646, 15.035]}
ELEMENT_DICT["Al"]   = {"Z":13, "Mass":26.9815, "Abundance":6.39, "IonizeEV" : [5.986, 18.826]}
ELEMENT_DICT["Si"]   = {"Z":14, "Mass":28.0860, "Abundance":7.52, "IonizeEV" : [8.151, 16.345]}
ELEMENT_DICT["P"]    = {"Z":15, "Mass":30.9738, "Abundance":5.52, "IonizeEV" : []}
ELEMENT_DICT["S"]    = {"Z":16, "Mass":32.0600, "Abundance":7.20, "IonizeEV" : []}
ELEMENT_DICT["Cl"]   = {"Z":17, "Mass":35.4530, "Abundance":5.6 , "IonizeEV" : []}
ELEMENT_DICT["Ar"]   = {"Z":18, "Mass":39.9480, "Abundance":6.8 , "IonizeEV" : []}
ELEMENT_DICT["K"]    = {"Z":19, "Mass":39.1020, "Abundance":4.95, "IonizeEV" : []}
ELEMENT_DICT["Ca"]   = {"Z":20, "Mass":40.0800, "Abundance":6.30, "IonizeEV" : [6.113, 11.871]}
ELEMENT_DICT["Sc"]   = {"Z":21, "Mass":44.9560, "Abundance":3.22, "IonizeEV" : []}
ELEMENT_DICT["Ti"]   = {"Z":22, "Mass":47.9000, "Abundance":5.13, "IonizeEV" : [6.82, 13.58]}
ELEMENT_DICT["V"]    = {"Z":23, "Mass":50.9414, "Abundance":4.40, "IonizeEV" : []}
ELEMENT_DICT["Cr"]   = {"Z":24, "Mass":51.9960, "Abundance":5.85, "IonizeEV" : [6.766, 16.50]}
ELEMENT_DICT["Mn"]   = {"Z":25, "Mass":54.9380, "Abundance":5.40, "IonizeEV" : [7.435, 15.640]}
ELEMENT_DICT["Fe"]   = {"Z":26, "Mass":55.8470, "Abundance":7.60, "IonizeEV" : [7.870, 16.16]}
ELEMENT_DICT["Co"]   = {"Z":27, "Mass":58.9332, "Abundance":5.1 , "IonizeEV" : []}
ELEMENT_DICT["Ni"]   = {"Z":28, "Mass":58.7100, "Abundance":6.30, "IonizeEV" : [7.635, 18.168]}
ELEMENT_DICT["Cu"]   = {"Z":29, "Mass":63.5460, "Abundance":4.5 , "IonizeEV" : []}
ELEMENT_DICT["Zn"]   = {"Z":30, "Mass":65.3700, "Abundance":4.2 , "IonizeEV" : []}
ELEMENT_DICT["Ga"]   = {"Z":31, "Mass":69.7200, "Abundance":2.4 , "IonizeEV" : []}
ELEMENT_DICT["Ge"]   = {"Z":32, "Mass":72.5900, "Abundance":2.9 , "IonizeEV" : []}
ELEMENT_DICT["As"]   = {"Z":33, "Mass":74.9216, "Abundance":2.3 , "IonizeEV" : []}
ELEMENT_DICT["Se"]   = {"Z":34, "Mass":78.9600, "Abundance":3.2 , "IonizeEV" : []}
ELEMENT_DICT["Br"]   = {"Z":35, "Mass":79.9040, "Abundance":2.6 , "IonizeEV" : []}
ELEMENT_DICT["Kr"]   = {"Z":36, "Mass":83.8000, "Abundance":3.2 , "IonizeEV" : []}
ELEMENT_DICT["Rb"]   = {"Z":37, "Mass":85.4678, "Abundance":2.4 , "IonizeEV" : []}
ELEMENT_DICT["Sr"]   = {"Z":38, "Mass":87.6200, "Abundance":2.85, "IonizeEV" : []}
ELEMENT_DICT["Y"]    = {"Z":39, "Mass":88.9059, "Abundance":1.8 , "IonizeEV" : []}
ELEMENT_DICT["Zr"]   = {"Z":40, "Mass":91.2200, "Abundance":2.5 , "IonizeEV" : []}
ELEMENT_DICT["Nb"]   = {"Z":41, "Mass":92.9060, "Abundance":2.0 , "IonizeEV" : []}
ELEMENT_DICT["Mo"]   = {"Z":42, "Mass":95.9400, "Abundance":1.92, "IonizeEV" : []}
#ELEMENT_DICT["Tc"]   = {"Z":43, "Mass":98.9060, "Abundance":?}
ELEMENT_DICT["Ru"]   = {"Z":44, "Mass":101.0700, "Abundance":1.60, "IonizeEV" : []}
ELEMENT_DICT["Rh"]   = {"Z":45, "Mass":102.9050, "Abundance":1.2 , "IonizeEV" : []}
## -- start : raplace Mass higher accuracy
ELEMENT_DICT["Pd"]   = {"Z":46, "Mass":106.4000, "Abundance":1.45, "IonizeEV" : []}
ELEMENT_DICT["Ag"]   = {"Z":47, "Mass":107.9000, "Abundance":0.80, "IonizeEV" : []}
ELEMENT_DICT["Cd"]   = {"Z":48, "Mass":112.4000, "Abundance":1.8 , "IonizeEV" : []}
ELEMENT_DICT["In"]   = {"Z":49, "Mass":114.8000, "Abundance":1.4 , "IonizeEV" : []}
ELEMENT_DICT["Sn"]   = {"Z":50, "Mass":118.7000, "Abundance":1.5 , "IonizeEV" : []}
ELEMENT_DICT["Sb"]   = {"Z":51, "Mass":121.8000, "Abundance":1.0 , "IonizeEV" : []}
ELEMENT_DICT["Te"]   = {"Z":52, "Mass":127.6000, "Abundance":2.0 , "IonizeEV" : []}
ELEMENT_DICT["I"]    = {"Z":53, "Mass":126.9000, "Abundance":1.4 , "IonizeEV" : []}
ELEMENT_DICT["Xe"]   = {"Z":54, "Mass":131.3000, "Abundance":2.0 , "IonizeEV" : []}
ELEMENT_DICT["Cs"]   = {"Z":55, "Mass":132.9000, "Abundance":1.1 , "IonizeEV" : []}
## -- end : raplace Mass higher accuracy
ELEMENT_DICT["Ba"]   = {"Z":56, "Mass":137.3400, "Abundance":1.95, "IonizeEV" : [5.212, 10.004]}

## make tuples
ELEMENT_SYMBOL  = tuple( ( key for key, val in ELEMENT_DICT.items() ) )
ELEMENT_Z       = tuple( ( val["Mass"] for key, val in ELEMENT_DICT.items() ) )
ELEMENT_MASS    = tuple( ( val["Mass"] for key, val in ELEMENT_DICT.items() ) )
ELEMENT_ABUN    = tuple( ( val["Abundance"] for key, val in ELEMENT_DICT.items() ) )
# in unit of [eV]
ELEMENT_IONIZPOTENTIAL    = tuple( ( _numpy.array(val["IonizeEV"], dtype=DT_NB_FLOAT) for key, val in ELEMENT_DICT.items() ) )

