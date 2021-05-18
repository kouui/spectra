

#-------------------------------------------------------------------------------
# function definition of Photoionization process
#-------------------------------------------------------------------------------
# VERSION
#
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
# 0.0.1
#    2021/05/08   u.k.
#        - func interpolate_PI_alpha : flip boundary value `_fill_value`
#-------------------------------------------------------------------------------

from ..ImportAll import *
from ..Math import Integrate

import numpy as _numpy
from scipy.interpolate import splrep as _splrep
from scipy.interpolate import splev as _splev



