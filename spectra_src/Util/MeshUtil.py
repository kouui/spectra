

#-------------------------------------------------------------------------------
# function definition of wavelength mesh utility
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
# construct wavelength mesh
#-------------------------------------------------------------------------------

def make_half_line_mesh_(nLambda : T_INT, qcore : T_FLOAT, qwing : T_FLOAT, 
                         q : T_ARRAY):
    """
    Construct half line mesh. Following RH's `getlambda.c`.
    called by it's wrapper function `makeLineMesh_Full`

    Parameters
    ----------

    nLambda : T_INT
        Number of points in full line mesh. Must be an odd number.

    qcore : T_FLOAT
        a parameter to control mesh density.
        there are always half of all points inside [-qcore,qcore].

    qwing : T_FLOAT
        a parameter to control how far away your mesh points reach.

    q : T_ARRAY
        a array to store output half line mesh distribution.

    Notes
    -------
    with q[0] = 0 and q[-1] = qwing, nLambda//2 points belonging to -qcore< q < +qcore. So,

    - change qwing, you change how far your Doppler width mesh reach
    - change qcore, you change how dense in line core.
    """
    if not BasicM.is_odd_(nLambda) : 
        raise ValueError( "`nLambda` should be an odd number." )

    nLhalf = nLambda//2 + 1

    if qwing <= 2*qcore:
        beta = 1.0
    else:
        beta = 0.5 * qwing / qcore

    y = beta + (beta*beta + (beta-1.)*nLhalf + 2. - 3.*beta)**(0.5)
    b = 2.0 * _numpy.log(y) / (nLhalf - 1)
    a = qwing / (nLhalf - 2. + y*y)

    for i in range(0, nLhalf):
        q[i] = a * (i + (_numpy.exp(b*i)-1.))


def make_full_line_mesh_(nLambda : T_INT, qcore : T_FLOAT = 2.5, qwing : T_FLOAT = 10.) -> T_ARRAY:
    """
    Construct Full line mesh.
    Calls inner function `makeLineMesh_Half` to construct half part
    and then make the use of symmetry.

    Parameters
    ----------

    nLambda : T_INT
        Number of points in full line mesh. Must be an odd number.

    qcore : T_FLOAT
        a parameter to control mesh density.
        there are always half of all points inside [-qcore,qcore].

    qwing : T_FLOAT
        a parameter to control how far away your mesh points reach.

    Returns
    --------
    mesh : T_ARRAY
        Full line mesh distribution.

    Notes
    -------
    with q[0] = 0 and q[-1] = qwing, nLambda//2 points belonging to -qcore< q < +qcore. So,

    - change qwing, you change how far your Doppler width mesh reach
    - change qcore, you change how dense in line core.
    """

    mesh = _numpy.empty(nLambda, dtype=DT_NB_FLOAT)
    nLmid = nLambda // 2
    make_half_line_mesh_(nLambda, qcore, qwing, mesh[nLmid:])
    mesh[:nLmid] = -1. * mesh[nLmid+1:][::-1]

    return mesh


def make_continuum_mesh_(nLambda : T_INT) -> T_ARRAY:
    """
    Given number of mesh points,
    compute a mesh distribution sampling most aroung wavelength edge.

    Parameters
    ----------

    nLambda : T_INT
        Number of points in Continuum mesh.

    Returns
    --------

    Mesh : T_ARRAY
        Output mesh distribution.

    """
    mesh = _numpy.empty(nLambda, dtype=DT_NB_FLOAT)
    for j in range(1, nLambda+1):
        qj = ( nLambda + 1.- j ) / nLambda
        mesh[j-1] = qj**(0.5)

    return mesh

def half_to_full_(arr_half : T_ARRAY, isMinus : T_BOOL = False) -> T_ARRAY:
    r"""
    create full (anti-)symmetric full array according to half array

    Parameters
    ----------

    arr_half : T_ARRAY, 1d
        half array.

    isMinus : T_BOOL
        True : anti-symmetric full array; False : symmetric full array

    Returns
    --------

    arr_full : T_ARRAY, 1d,
        full array.

    """
    _nLmid = arr_half.shape[0]
    _nLfull = (_nLmid-1) * 2 + 1
    _arr_full = _numpy.zeros(_nLfull, dtype=arr_half.dtype)

    fac = -1 if isMinus else 1

    _arr_full[_nLmid:] = arr_half[1:]
    _arr_full[:_nLmid] = fac * arr_half[::-1]

    return _arr_full

def array_from_1D_(arr_1D : T_ARRAY, mesh_idxs : T_ARRAY, k : T_INT) -> T_ARRAY :
    """extract target sub-array from a 1D array containing sub-arrays with different size

    Parameters
    ----------
    arr_1D : T_ARRAY
        the 1D array containing sub-arrays with different size
    mesh_idxs : T_ARRAY, 2d, (nLine, 2)
        start and end index of each sub-array
    k : T_INT
        the index of line you want to extract

    Returns
    -------
    T_ARRAY
        target sub-array
    """
    i1 : T_INT = mesh_idxs[k,:]
    i2 : T_INT = mesh_idxs[k,:]

    return arr_1D[i1:i2].copy()


#-----------------------------------------------------------------------------
# numba optimization
#-----------------------------------------------------------------------------

if CFG._IS_JIT:

    make_half_line_mesh_      = nb_njit( **NB_NJIT_KWGS ) (make_half_line_mesh_)
    make_full_line_mesh_      = nb_njit( **NB_NJIT_KWGS ) (make_full_line_mesh_)