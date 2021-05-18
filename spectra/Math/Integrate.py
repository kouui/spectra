
#-------------------------------------------------------------------------------
# function definition of integration method
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

from ..ImportAll import *

import numpy as _numpy

#-----------------------------------------------------------------------------
# Trapzoidal integral
#-----------------------------------------------------------------------------

def trapze_(y : T_ARRAY, x : T_ARRAY) -> T_FLOAT:
    """
    Integration using Trapzoidal rule

    Parameters
    ----------
    y : T_ARRAY
        integrand as a function of variable x

    x : T_ARRAY
        independent variable x

    Returns
    -------

    sum : T_FLOAT
         result of integration

    Notes
    ------
    Refer to [1]_.

    References
    ----------
    .. [1] Trapzoidal rule, wikipedia, https://en.wikipedia.org/wiki/Trapezoidal_rule
    """

    n = x.size
    dx = _numpy.empty(n, dtype=T_FLOAT)
    for i in range(1, n-1):
        dx[i] = 0.5 * (x[i+1]-x[i-1])
    dx[0] = 0.5 * (x[1]-x[0])
    dx[n-1] = 0.5 * (x[n-1] - x[n-2])

    sum = 0.0
    for i in range(n):
        sum += dx[i] * y[i]
    
    return sum

#-----------------------------------------------------------------------------
# simpson integral
#-----------------------------------------------------------------------------

def _simps_odd_evenspaced_(y : T_ARRAY, dx : T_UNION[T_INT,T_FLOAT]) -> T_FLOAT:
    """ simpson integral for evenly spaced odd size array

    Parameters
    ----------
    y : T_ARRAY
        integrand
    dx : T_UNION[T_INT,T_FLOAT]
        interval of independent variable

    Returns
    -------
    T_FLOAT
        integrated result
    """
    N = y.size

    result = dx / 3. * ( y[0:N-2:2] + 4*y[1:N-1:2] + y[2:N:2] )
    result = _numpy.sum( result )

    return result

def _simps_odd_(y : T_ARRAY, x : T_ARRAY) -> T_FLOAT:
    """simpson integral for odd size array

    Parameters
    ----------
    y : T_ARRAY
        integrand
    x : T_ARRAY
        independent variable

    Returns
    -------
    T_FLOAT
        integrated result
    """
    N = y.size

    h = x[1:] - x[:-1]
    h0 = h[0:N-2:2]
    h1 = h[1:N-1:2]
    hsum = h0 + h1
    hprod = h0 * h1
    h0divh1 = h0 / h1
    result = hsum / 6. * ( y[0:N-2:2] * (2. - 1./h0divh1) +
                           y[1:N-1:2] * hsum * hsum / hprod +
                           y[2:N:2] * (2. - h0divh1) )
    result = _numpy.sum( result )

    return result

def _simps_even_evenspaced_(y : T_ARRAY, dx : T_UNION[T_INT,T_FLOAT]) -> T_FLOAT:
    """simpson integral for evenly spaced even size array

    Parameters
    ----------
    y : T_ARRAY
        integrand
    dx : T_UNION[T_INT,T_FLOAT]
        interval of independent variable

    Returns
    -------
    T_FLOAT
        integrated variable
    """

    N = y.size

    last_dx = dx
    first_dx = dx
    val = 0.
    result : T_FLOAT = 0.

    val += 0.5 * last_dx * (y[-2] + y[-1])
    result += _simps_odd_evenspaced_(y[:N-1], dx)
    val += 0.5 * first_dx * (y[0] + y[1])
    result += _simps_odd_evenspaced_(y[1:], dx)

    result = 0.5 * (result + val)
    
    return result

def _simps_even_(y : T_ARRAY, x : T_ARRAY) -> T_FLOAT:
    """simpson integral for even size array

    Parameters
    ----------
    y : T_ARRAY
        integrand
    dx : T_ARRAY
        independent variable

    Returns
    -------
    T_FLOAT
        integrated variable
    """
    N = y.size

    last_dx = x[-1] - x[-2]
    first_dx = x[1] - x[0]
    val = 0.
    result = 0.

    val += 0.5 * last_dx * (y[-2] + y[-1])
    result += _simps_odd_(y[:N-1], x[:N-1])
    val += 0.5 * first_dx * (y[0] + y[1])
    result += _simps_odd_(y[1:], x[1:])

    result = 0.5 * (result + val)
    return result


from . import BasicM

def simpson_(y : T_ARRAY, 
             x : T_UNION[T_ARRAY, None] = None, 
             dx : T_UNION[T_FLOAT, T_INT, None] = None) -> T_FLOAT:
    """
    Integration using Simpson's rule

    Parameters
    ----------
    y : float, 1darray
        integrand as a function of variable x

    x : float, 1darray
        independent variable x

    dx : float,
        inteval of x if x is evenly spaced

    Returns
    -------

    _sum : float
         result of integration
    """
    N = y.size

    if BasicM.is_odd_(N):

        if dx is not None:
            return _simps_odd_evenspaced_(y, dx)
        elif x is not None:
            return _simps_odd_(y, x)
        else:
            dx = 1
            return _simps_odd_evenspaced_(y, dx)

    else:
        if dx is not None:
            return _simps_even_evenspaced_(y, dx)
        elif x is not None:
            return _simps_even_(y, x)
        else:
            dx = 1
            return _simps_even_evenspaced_(y, dx)

#-----------------------------------------------------------------------------
# numba optimization
#-----------------------------------------------------------------------------

if CFG._IS_JIT:

    trapze_                   = nb_njit( **NB_NJIT_KWGS ) ( trapze_ )
    _simps_odd_evenspaced_    = nb_njit( **NB_NJIT_KWGS ) ( _simps_odd_evenspaced_ )
    _simps_odd_               = nb_njit( **NB_NJIT_KWGS ) ( _simps_odd_ )
    _simps_even_evenspaced_   = nb_njit( **NB_NJIT_KWGS ) ( _simps_even_evenspaced_ )
    _simps_even_              = nb_njit( **NB_NJIT_KWGS ) ( _simps_even_ )
    simpson_                  = nb_njit( **NB_NJIT_KWGS ) ( simpson_ )