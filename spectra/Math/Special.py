
#-------------------------------------------------------------------------------
# function definition of Special functions
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

from ..ImportAll import *
import numpy as _numpy


#-------------------------------------------------------------------------------
# Exponential integral E1(x) and E2(x)
#
#       En(x) = \int_{1}^{\infty} t^{1/n}e^{-xt}dt, x>0, n=0,1,...
#
#       with relation: En+1 = 1/n (exp(-x) - xEn(x))
#-------------------------------------------------------------------------------

_A53 = _numpy.array([-0.57721566,  0.99999193, -0.24991055,
                0.05519968, -0.00976004,  0.00107857], dtype=T_FLOAT)
_A56 = _numpy.array([8.5733287401, 18.0590169730, 8.6347608925,  0.2677737343],dtype=T_FLOAT)
_B56 = _numpy.array([9.5733223454, 25.6329561486,21.0996530827,  3.9584969228],dtype=T_FLOAT)

def E0_(x : T_VEC_IFA) -> T_VEC_FA:
    """E_0(x)

    Parameters
    ----------
    x : T_VEC_IFA
        x

    Returns
    -------
    T_VEC_FA
        E_0(x)
    """
    return _numpy.exp(-x) / x


def E1_(x : T_VEC_IFA) -> T_VEC_FA:
    """
    Approximated formula for Exponential integral :math:`E_1(x)`.

    Parameters
    ----------

    x : T_VEC_IFA
        independent variable x

    Returns
    -------

    T_VEC_FA
         1st order Exponential integral of x

    Notes
    -----
    Refer to [1]_.

    References
    ----------

    .. [1] D.ABarry, J.-Y.Parlange, "Approximation for the exponential integral (Theis well function)",
           Journal of Hydrology, Volume 227, Issues 1â€“4, 31 January 2000, Pages 287-291

    """
    if x <= 0. or x> 80.0:
        print(x)
        raise ValueError("argument x should be a positive number smaller than 80.0")
    

    if x <= 1.0:
        E1 = -_numpy.log(x) + _A53[0] + x*(_A53[1] + x*(_A53[2] + x*(_A53[3] + x*(_A53[4] + x*_A53[5]))))
    else:
        E1  = _A56[3]/x +  _A56[2] + x*(_A56[1] + x*(_A56[0] + x))
        E1 /= _B56[3] + x*(_B56[2] + x*(_B56[1] + x*(_B56[0] + x)))
        E1 *= _numpy.exp(-x)

    return E1

def E2_(x : T_VEC_IFA) -> T_VEC_FA:
    """
    Calculate Exponential integral :math:`E_2(x)` from :math:`E_1(x)`.

    Parameters
    ----------

    x : T_VEC_IFA
        independent variable x

    Returns
    -------

    T_VEC_FA
         2nd order Exponential integral of x

    Notes
    ------
    Using relation [1]_ Equation(11.127):

    .. math:: E_{n+1}(x) = [e^{-x} - xE_{n}(x)] / n

    References
    ----------
    .. [1] Ivan Hubeny, Dimitri Mihalas, "Theory of Stellar Atmosphere:
        An Introduction to Astrophysical Non-equilibrium
        Quantitative Spectroscopic Analysis",
        Princeton University Press, pp. 365, 2015.
    """
    return _numpy.exp(-x) - x * E1_(x)

def E3_(x : T_VEC_IFA) -> T_VEC_FA:
    """Calculate Exponential integral :math:`E_3(x)` from :math:`E_2(x)`.

    Parameters
    ----------
    x : T_VEC_IFA
        independent variable x

    Returns
    -------
    T_VEC_FA
        3rd order Exponential integral of x
    """
    return 0.5 * ( _numpy.exp(-x) - x * E2_(x) )

#-----------------------------------------------------------------------------
# numba optimization
#-----------------------------------------------------------------------------

if CFG._IS_JIT:

    E0_ = nb_vec( **NB_VEC_KWGS ) ( E0_ )
    E1_ = nb_vec( **NB_VEC_KWGS ) ( E1_ )
    E2_ = nb_vec( **NB_VEC_KWGS ) ( E2_ )
    E3_ = nb_vec( **NB_VEC_KWGS ) ( E3_ )

else:
    E1_ = _numpy.vectorize(E1_, otypes=[T_FLOAT], cache=True)
