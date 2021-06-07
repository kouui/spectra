
from numpy import array as _array
from numpy import isclose as _ISCLOSE
from numpy import allclose as _ALLCLOSE
_KWGS_CLOSE = { "rtol" : 1.E-05,  "atol" : 1.E-20}

import os, unittest
from spectra_src.ImportAll import *

from spectra_src.Math import GaussLeg


class Test_Gauss_Leg(unittest.TestCase):

    def test_gauss_quad_coe_odd_(self):

        a, b, n = -1, 1., 3
        xs, ws = GaussLeg.gauss_quad_coe_(a, b, n)
        xs_correct = (-0.7745966692414834, 0., +0.7745966692414834)
        ws_correct = (0.5555555555555556, 0.8888888888888888, 0.5555555555555556)
        
        self.assertTrue( _ALLCLOSE( xs, xs_correct , **_KWGS_CLOSE ) )
        self.assertTrue( _ALLCLOSE( ws, ws_correct , **_KWGS_CLOSE ) )
    
    def test_gauss_quad_coe_even_(self):
        
        a, b, n = -1, 1., 4
        xs, ws = GaussLeg.gauss_quad_coe_(a, b, n)
        xs_correct = (-0.8611363115940526, -0.3399810435848563, +0.3399810435848563, +0.8611363115940526)
        ws_correct = (0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538)
        
        self.assertTrue( _ALLCLOSE( xs, xs_correct , **_KWGS_CLOSE ) )
        self.assertTrue( _ALLCLOSE( ws, ws_correct , **_KWGS_CLOSE ) )




if __name__ == "__main__":

    unittest.main()

