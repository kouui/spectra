

#-------------------------------------------------------------------------------
# function definition of optically thin process
#-------------------------------------------------------------------------------
# VERSION
#
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

from ..ImportAll import *
import numpy as _numpy

def relative_flux_(AJI : T_UNION[T_FLOAT,T_ARRAY], 
                   f0  : T_UNION[T_FLOAT,T_ARRAY], 
                   nj : T_UNION[T_INT,T_ARRAY]) -> T_UNION[T_FLOAT,T_ARRAY]:
    """
    calculate optically thin relative flux (some constants are removed)

    Parameters
    ----------

    AJI : T_UNION[T_FLOAT,T_ARRAY]
        Einstein A coefficient, 
        [:math:`s^{-1}`]

    f0 : T_UNION[T_FLOAT,T_ARRAY]
        Transition line frequency, 
        [:math:`hz`]

    nj : T_UNION[T_INT,T_ARRAY]
        Population of the upper level of line transition, 
        [:math:`cm^{-3}`]

    Returns
    -------

    T_UNION[T_FLOAT,T_ARRAY]
        relative flux under the assumption of optically thin. 
        [:math:`erg \; cm^{-3} \; s^{-1}`]

    Notes
    -----

    The absolute flux under the assumption of optically thin [1]_.

    .. math:: F_{ji} = \frac{1}{4 \pi R^{2}} \int_{\Delta V} \epsilon_{ji} dV \quad [erg \; cm^{-2} \; s^{-1}]

    where the emissivity :math:`\epsilon_{ji}` is

    .. math:: \epsilon_{ji} = h \nu n_{j} A_{ji} \quad [erg \; cm^{-3} \; s^{-1}]

    and

    .. math:: \epsilon_{ji} = \int_{\nu} \epsilon_{\nu} d \nu = \int_{\nu} h \nu n_{j} A_{ji} \psi d \nu

    So our relative flux is given by

    .. math: h \nu n_{j} A_{ji} \quad [erg \; cm^{-3} \; s^{-1}]

    References
    ----------

    .. [1] John T. Mariska, "The Solar Transition Region",
        Cambridge University Press, pp. 19, 1992

    """

    return CST.h_ * f0 * nj * AJI


#-----------------------------------------------------------------------------
# numba optimization
#-----------------------------------------------------------------------------

if CFG._IS_JIT:
    relative_flux_ = nb_vec(**NB_VEC_KWGS) ( relative_flux_ )
