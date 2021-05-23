

#-------------------------------------------------------------------------------
# function definition of Solving Statistial Equilibrium equations
#-------------------------------------------------------------------------------
# VERSION
#
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

from ..ImportAll import *

import numpy as _numpy

#-------------------------------------------------------------------------------
# fill radiative/collisional transition matrix 
#-------------------------------------------------------------------------------

def set_matrixC_(
         Cmat : T_ARRAY, 
         Cji  : T_ARRAY, 
         Cij  : T_ARRAY, 
         idxI : T_ARRAY, 
         idxJ : T_ARRAY, 
         Ne   : T_UNION[T_FLOAT, T_INT]):
    r"""
    Compute the collisional rate matrix.

    Parameters
    ----------

    Cmat : T_ARRAY, 2d
        collisional rate matrix, a 2D Array to store computed results, 
        [:math:`s^{-1}`]

    Cji : T_ARRAY, 1d
        downward collisional transition rate, 
        [:math:`s^{-1} \cdot cm^{3}`]

    Cij : T_ARRAY, 1d
        upward collisional transition rate, 
        [:math:`s^{-1} \cdot cm^{3}`]

    idxI : T_ARRAY, 1d
        level index of lower level i, [-]

    idxJ : T_ARRAY, 1d
        level index of upper level j, [-]

    Ne: T_UNION[T_FLOAT, T_INT]
        electron density, 
        [:math:`cm^{-3}`]

    Notes
    -----
    Refer to [1]_ Equation(9.80).

    .. math:: \sum_{j \neq i} n_j (R_{ji}+n_{e} C_{ji}) - n_i \sum_{j \neq i}(R_{ij}+n_{e} C_{ij}) = 0


    References
    ----------

    .. [1] Ivan Hubeny, Dimitri Mihalas, "Theory of Stellar Atmosphere:
        An Introduction to Astrophysical Non-equilibrium
        Quantitative Spectroscopic Analysis",
        Princeton University Press, pp. 282, 2015.

    """

    n_row, n_col = Cmat.shape
    if n_row != n_col:
        raise ValueError( "Cmat should be a squared matrix" )

    nTran = Cji.size
    for k in range(nTran):
        i, j = idxI[k], idxJ[k]
        Cmat[i,j] += Ne * Cji[k]
        Cmat[j,i] += Ne * Cij[k]

def set_matrixR_(
         Rmat      : T_ARRAY, 
         Rji_spon  : T_ARRAY,
         Rji_stim  : T_ARRAY, 
         Rij       : T_ARRAY, 
         idxI      : T_ARRAY, 
         idxJ      : T_ARRAY):
    r"""
    Compute the radiative rate matrix.

    Parameters
    ----------

    Rmat : T_ARRAY, 2d
        radiative rate matrix, a 2D Array to store computed results, [:math:`s^{-1}`]

    Rji_spon : T_ARRAY, 1d
        spontaneous radiative transition rate, [:math:`s^{-1}`]

    Rji_stim : T_ARRAY, 1d
        spontaneous radiative transition rate, [:math:`s^{-1}`]

    Rij : T_ARRAY, 1d
        upward radiative transition rate, [:math:`s^{-1}`]

    idxI : T_ARRAY, 1d
        level index of lower level i, [-]

    idxJ : T_ARRAY, 1d
        level index of upper level j, [-]

    Notes
    -----
    Refer to [1]_ Equation(9.80).

    .. math:: \sum_{j \neq i} n_j (R_{ji}+C_{ji}) - n_i \sum_{j \neq i}(R_{ij}+C_{ij}) = 0


    References
    ----------

    .. [1] Ivan Hubeny, Dimitri Mihalas, "Theory of Stellar Atmosphere:
        An Introduction to Astrophysical Non-equilibrium
        Quantitative Spectroscopic Analysis",
        Princeton University Press, pp. 282, 2015.

    """
    n_row, n_col = Rmat.shape
    if n_row != n_col:
        raise ValueError( "Rmat should be a squared matrix" )

    nTran = Rji_spon.size
    for k in range(nTran):
        i, j = idxI[k], idxJ[k]
        Rmat[i,j] += Rji_spon[k] + Rji_stim[k]
        Rmat[j,i] += Rij[k]

#-------------------------------------------------------------------------------
# solving system of equations
#-------------------------------------------------------------------------------

def solve_SE_(Rmat : T_ARRAY, Cmat : T_ARRAY) -> T_ARRAY:
    r"""Solve the linear equation system of statistical equilibrium.

    Parameters
    ----------

    Rmat : T_ARRAY, (nLevel,nLevel)
        radiative transition rate matrix, 
        [:math:`s^{-1}`]

    Cmat : T_ARRAY, (nLevel,nLevel)
        collisional transition rate matrix,
        [:math:`s^{-1}`]

    Returns
    -------

    nArr : T_ARRAY, (nLevel,)
        normalized level population. 
        [:math:`cm^{-3}`]

    """

    nLevel = Rmat.shape[0]
    A = Cmat[:,:] + Rmat[:,:]
    b = _numpy.zeros(nLevel, dtype=DT_NB_FLOAT)

    #-------------------------------------------------------------
    # diagnal components
    #-------------------------------------------------------------
    for k in range(nLevel):
        A[k,k] = -A[:,k].sum()

    #-------------------------------------------------------------
    # abundance definition equation
    #-------------------------------------------------------------
    A[-1,:] = 1.
    b[-1] = 1.

    nArr = _numpy.linalg.solve(A, b)

    return nArr

#-------------------------------------------------------------------------------
# numba optimization
#-------------------------------------------------------------------------------

if CFG._IS_JIT:

    set_matrixC_ = nb_njit(**NB_NJIT_KWGS) (set_matrixC_)
    set_matrixR_ = nb_njit(**NB_NJIT_KWGS) (set_matrixR_)
    solve_SE_    = nb_njit(**NB_NJIT_KWGS) (solve_SE_)