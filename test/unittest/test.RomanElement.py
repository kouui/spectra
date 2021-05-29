
from numpy import array as _array
from numpy import isclose as _ISCLOSE
from numpy import allclose as _ALLCLOSE
_KWGS_CLOSE = { "rtol" : 1.E-05,  "atol" : 1.E-20}

import os, unittest

from spectra_src.ImportAll import *

from spectra_src.Util import RomanUtil as _RomanUtil
from spectra_src.Util import ElementUtil as _ElementUtil


class Test_Roman(unittest.TestCase):

    def test_index_to_roman(self):

        self.assertTrue( _RomanUtil.index_to_roman_(10) == "X"  )

    def test_roman_to_index(self):

        self.assertTrue( _RomanUtil.roman_to_index_("IV") == 4  )

    def test_shift_roman(self):

        self.assertTrue( _RomanUtil.shift_roman_("IV", 10) == "XIV"  )

class Test_Element(unittest.TestCase):

    def test_sym_to_z(self):
        
        self.assertTrue( _ElementUtil.sym_to_z_("He") == 2  )

    def test_sym_to_mass(self):

        self.assertTrue( _ISCLOSE( _ElementUtil.sym_to_mass_("He"), 4.0026 ) )

    def test_sym_to_Abun(self):

        self.assertTrue( _ISCLOSE( _ElementUtil.sym_to_abun_("He"), 10.93 ) )

    def test_format_ion(self):

        self.assertTrue( _ElementUtil.format_ion_("fe_iii"), "Fe_III" )

    def test_ion_to_sym_and_stage(self):

        self.assertTrue( _ElementUtil.ion_to_sym_and_stage_("fe_iii"), ("Fe", "III") )

    def test_sym_and_stage_to_ion(self):

        self.assertTrue( _ElementUtil.sym_and_stage_to_ion_("He", "i"), "He_I" )

    def test_ion_to_ioniz_potential(self):

        self.assertTrue( _ISCLOSE( _ElementUtil.ion_to_ioniz_potential_("H_I"), 13.508*CST.eV2erg_ ) )

    def test_shift_ion(self):

        self.assertTrue( _ElementUtil.shfit_ion_("fe_iii", -2) == "Fe_I" )
        self.assertTrue( _ElementUtil.shfit_ion_("He_I", 1) == "He_II" )
    

if __name__ == "__main__":

    unittest.main()