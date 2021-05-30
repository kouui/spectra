
from numpy import array as _array
from numpy import isclose as _ISCLOSE
from numpy import allclose as _ALLCLOSE
_KWGS_CLOSE = { "rtol" : 1.E-05,  "atol" : 1.E-20}

import os, unittest

from spectra_src.ImportAll import *

from spectra_src.RadiativeTransfer import Profile as _Profile

import numpy as _numpy


class Test_Profile(unittest.TestCase):

    def test_voigt_gaussian(self):

        x = _numpy.linspace(-3, 3, 61)
        a = 0
        gaussian = _Profile.gaussian_(x)
        voigt = _Profile.voigt_(a, x)
        self.assertTrue( _ALLCLOSE( voigt, gaussian , **_KWGS_CLOSE ) )

    def test_voigt_hf(self):

        x = _numpy.linspace(-3, 3, 61)
        a = 0.1
        h, f = _Profile.hf_(a, x)
        voigt = _Profile.voigt_(a, x)
        self.assertTrue( _ALLCLOSE( voigt, h , **{ "rtol" : 1.E-03,  "atol" : 1.E-20} ) )

if __name__ == "__main__":

    unittest.main()