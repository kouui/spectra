
#-------------------------------------------------------------------------------
# function definition of absorption profile related process
#-------------------------------------------------------------------------------
# VERSION
#
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

from ..ImportAll import *
from ..Math import BasicM

import numpy as _numpy


#-------------------------------------------------------------------------------
# absorption profile
#   - voigt
#   - gaussian
#-------------------------------------------------------------------------------


#@nb_vec( **NB_VEC_KWGS )
def voigt_(a : T_VEC_IFA, x : T_VEC_IFA) -> T_VEC_IFA:
    r"""
    Calculate Doppler width normalized voigt function using polynomial fitting formula.

    `voigt(a,x)` function itself is not a vectorized function, only available to scalar operation.
    so we applied

    `nb.vectorize`.

    Parameters
    ----------

    a : T_VEC_IFA
        damping constant normalized by Doppler width, [-]
    x : T_VEC_IFA
        Doppler width normalized mesh, 
        [-]

    Returns
    -------

    res : float
        voigt function, normalized to 1, 
        [-]

    Notes
    -----

    The Voigt function is not normalized but has area :math:`\sqrt{\pi}` in `x` unit

    This is a combination of Humlicek(1982) [1]_ and Hui et al.(1978) [2]_ methods.

    When :math:`a > 0.1`, one can not ignore the wing component of Voigt function.
    That is, to guarantee its normalization, one has to take care of
    whether the mesh points are wide enough.

    References
    ----------

    .. [1] J.Humlicek, 'Optimized computation of the voigt and complex probability functions',
        Journal of Quantitative Spectroscopy and Radiative Transfer (JQSRT),
        Volume 27, Issue 4, April 1982, Pages 437-444.

    .. [2] A.K.Hui, B.H.Armstrong, A.A.Wray, 'Rapid computation of the Voigt and complex error functions',
        Journal of Quantitative Spectroscopy and Radiative Transfer (JQSRT),
        Volume 19, Issue 5, May 1978, Pages 509-516.
    """

    if a < 0.01:
        Z = a - 1j*x
        U = Z * Z
        A0 = 36183.31; B0 = 32066.6
        A1 = 3321.9905;B1 = 24322.84
        A2 = 1540.787; B2 = 9022.228
        A3 = 219.0313; B3 = 2186.181
        A4 = 35.76683; B4 = 364.2191
        A5 = 1.320522; B5 = 61.57037
        A6 = .56419  ; B6 = 1.841439
        F = (_numpy.exp(U) - Z * ( A0 - U * ( A1 - U * ( A2 - U * ( A3 - U * ( A4 - U * ( A5 - U * A6 )))))) /
         ( B0 - U * ( B1 - U * ( B2 - U * ( B3 - U * ( B4 - U * ( B5 - U * ( B6 - U ))))))))
    else:
        A0 = 122.607931777; B0 = 122.607931774
        A1 = 214.382388695; B1 = 352.730625111
        A2 = 181.928533092; B2 = 457.334478784
        A3 = 93.1555804581; B3 = 348.703917719
        A4 = 30.1801421962; B4 = 170.354001821
        A5 = 5.91262620977; B5 = 53.9929069129
        A6 = .564189583563; B6 = 10.4798571143

        Z = a - 1j*_numpy.abs(x)# * C.k[0]
        F = ( ( A0 + Z * ( A1 + Z * ( A2 + Z * ( A3 + Z * ( A4 + Z * ( A5 + Z * A6 ) ) ) ) ) ) /
             ( B0 + Z * ( B1 + Z * ( B2 + Z * ( B3 + Z * ( B4 + Z * ( B5 + Z * ( B6 + Z ) ) ) ) ) ) ) )

    res = F.real / CST.sqrtPi_
    return res

def gaussian_(x : T_VEC_IFA) -> T_VEC_IFA:
    r"""
    Calculate Doppler width normalized gaussian profile.

    Parameters
    ----------

    x : T_VEC_IFA
        Doppler width normalized mesh, 
        [-]

    Returns
    -------

    T_VEC_IFA
        gaussian profile, normalized to 1, 
        [-]

    Notes
    -----

    This formula refers to [1]_.

    References
    ----------

    .. [1] Ivan Hubeny, Dimitri Mihalas, "Theory of Stellar Atmosphere:
        An Introduction to Astrophysical Non-equilibrium
        Quantitative Spectroscopic Analysis",
        Princeton University Press, pp. 203, Eq(8.24), 2015.
    """
    return _numpy.exp(-x*x) / CST.sqrtPi_



#-----------------------------------------------------------------------------
# numba optimization
#-----------------------------------------------------------------------------

if CFG._IS_JIT:
    voigt_    = nb_vec(**NB_VEC_KWGS) (voigt_)
    gaussian_ = nb_vec(**NB_VEC_KWGS) (gaussian_)

else:
    voigt_    = _numpy.vectorize(voigt_, otypes=[T_FLOAT], cache=True)