
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

from numpy import abs as _abs
from numpy import exp as _exp


#-------------------------------------------------------------------------------
# absorption profile
#   - voigt
#   - gaussian
#-------------------------------------------------------------------------------

@OVERLOAD
def voigt_(a : T_FLOAT, x : T_FLOAT) -> T_FLOAT : ...
@OVERLOAD
def voigt_(a : T_ARRAY, x : T_ARRAY) -> T_ARRAY : ...
@OVERLOAD
def voigt_(a : T_ARRAY, x : T_FLOAT) -> T_ARRAY : ...
@OVERLOAD
def voigt_(a : T_FLOAT, x : T_ARRAY) -> T_ARRAY : ...

def voigt_(a : T_VEC_FA, x : T_VEC_FA) -> T_VEC_FA:
    r"""
    Calculate Doppler width normalized voigt function using polynomial fitting formula.

    `voigt(a,x)` function itself is not a vectorized function, only available to scalar operation.
    so we applied

    `nb.vectorize`.

    Parameters
    ----------

    a : T_VEC_FA
        damping constant normalized by Doppler width, [-]
    x : T_VEC_FA
        Doppler width normalized mesh, 
        [-]

    Returns
    -------

    res : T_VEC_FA
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

        Z = a - 1j*_abs(x)# * C.k[0]
        F = ( ( A0 + Z * ( A1 + Z * ( A2 + Z * ( A3 + Z * ( A4 + Z * ( A5 + Z * A6 ) ) ) ) ) ) /
             ( B0 + Z * ( B1 + Z * ( B2 + Z * ( B3 + Z * ( B4 + Z * ( B5 + Z * ( B6 + Z ) ) ) ) ) ) ) )

    res = F.real / CST.sqrtPi_
    return res

@OVERLOAD
def gaussian_(x : T_FLOAT) -> T_FLOAT : ...
@OVERLOAD
def gaussian_(x : T_ARRAY) -> T_ARRAY : ...

def gaussian_(x : T_VEC_FA) -> T_VEC_FA:
    r"""
    Calculate Doppler width normalized gaussian profile.

    Parameters
    ----------

    x : T_VEC_FA
        Doppler width normalized mesh, 
        [-]

    Returns
    -------

    T_VEC_FA
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
    return _exp(-x*x) / CST.sqrtPi_


def hf_(a : T_FLOAT, x : T_ARRAY) -> T_TUPLE[T_ARRAY,T_ARRAY]:
    """compute Voigt and Faraday-Voigt function
    by comprex probability function w(z)=exp(-z^2)*erfc(-i*z)

    z=x+i*y  for y>=0
                     |oo    exp(-t^2)
    H(a,v) = a/!pi * |  ---------------- dt    = u(x,y)  {a=y,v=x}
                     |-oo  a^2+(v-t)^2

                     |oo (v-t)exp(-t^2)
    F(a,v) = 1/!pi * |  ---------------- dt    = v(x,y)  {a=y,v=x}
                     |-oo  a^2+(v-t)^2
    
    return H(a,v)/sqrt(!pi),  F(a,v)/sqrt(!pi)

    Parameters
    ----------
    a : T_FLOAT
        dampling constant a normalized doppler width (Hz)
    x : T_ARRAY, 1d
        wavelength mesh normalized by its doppler width

    Returns
    -------
    T_TUPLE[T_ARRAY,T_ARRAY]
        
    h : T_ARRAY, 1d
        normalized voigt function
    f : T_ARRAY, 1d
        normalized Faraday-Voigt function

    Notes
    ------
    	k.i. '94/02/27
    	k.i. '04/11/03	double precesion
    	k.i. '05/01/05	bug fixed for s>15.
    	k.i. '20/12/11	to python
        k.i. '20/12/26  return H(a,v)/sqrt(pi), F(a,v)/sqrt(pi), x2 for F

    
    References
    ------------
    .. [1] Humlicek,J. 1982, J.Quant.Spectrosc.Radiat.Transfer, vol.27, No.4, pp437-444.
    """
    
    x_1d = x.reshape(-1)

    w4_1d = _numpy.empty( x.size, dtype = DT_NB_COMPLEX )

    for i in range(x.size):
        v = x_1d[i]
        
        t = a - 1j * v
        s = _abs(v) + a
        
        if s >= 15.:
            val = t * 0.5641896 / ( 0.5 + t*t )
        elif 5.5 <= s < 15.:
            u = t*t
            val = t * ( 1.410474 + u*0.5641896 ) / ( 0.75 + u * ( 3.+u ) ) 
        else:
            if a >= 0.195 * _abs(v) - 0.176 :
                val = (16.4955+t*(20.20933+t*(11.96482+t*(3.778987+t*0.5642236)))) /  \
                      (16.4955+t*(38.82363+t*(39.27121+t*(21.69274 +t*(6.699398+t)))))
            else:
                u = t*t
                val = _exp(u)-t * \
                    (36183.31-u*(3321.9905-u*(1540.787-u*(219.0313-u*(35.76683-u*(1.320522-u*0.56419)))))) / \
                    (32066.6-u*(24322.84-u*(9022.228-u*(2186.181-u*(364.2191-u*(61.57037-u*(1.841439-u)))))))

        w4_1d[i] = val
    
    w4 = w4_1d.reshape(x.shape)
    h = w4.real / CST.sqrtPi_
    f = w4.imag / CST.sqrtPi_

    return h, f

#-----------------------------------------------------------------------------
# numba optimization
#-----------------------------------------------------------------------------

if CFG._IS_JIT:
    voigt_    = nb_vec(**NB_VEC_KWGS) (voigt_)
    gaussian_ = nb_vec(**NB_VEC_KWGS) (gaussian_)
    hf_       = nb_njit(**NB_NJIT_KWGS) (hf_)

else:
    voigt_    = np_vec(voigt_, **NP_VEC_KWGS)