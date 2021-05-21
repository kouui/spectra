

#-------------------------------------------------------------------------------
# function definition of Photoionization process
#-------------------------------------------------------------------------------
# VERSION
#
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
# 0.0.1
#    2021/05/08   u.k.
#        - func interpolate_PI_alpha : flip boundary value `_fill_value`
#-------------------------------------------------------------------------------

from ..ImportAll import *
from ..Math import Integrate

import numpy as _numpy
#from scipy.interpolate import splrep as _splrep
#from scipy.interpolate import splev as _splev
from scipy.interpolate import interp1d as _interp1d # type: ignore


#-------------------------------------------------------------------------------
# interpolate data
#-------------------------------------------------------------------------------
def interpolate_PI_intensity_(backRad : T_ARRAY, continuum_mesh : T_ARRAY) -> T_ARRAY:
    """Given continuum mesh, interpolate background intensity profile

    Parameters
    ----------

    backRad : T_ARRAY
        a table of background intensity, wavelength_cm vs intensity_cm

    continuum_mesh : T_ARRAY
        continuum mesh

    Returns
    -------

    T_ARRAY
        interpolated continuum intensity profile.
    """

    ## scipy interpolation
    #_fill_value = (_backRad[1,0],_backRad[1,-1])
    #_bsp_obj = interp1d(x=_backRad[0,:], y=_backRad[1,:], bounds_error=False, fill_value=_fill_value)

    intensity_mesh = _numpy.empty(continuum_mesh.shape, dtype=T_FLOAT)
    intensity_mesh_1d = intensity_mesh.reshape(-1)
    intensity_mesh_1d[:] = _numpy.interp(continuum_mesh[:,:].reshape(-1), backRad[0,:], backRad[1,:])
    #for k in range(continuum_mesh.shape[0]):
        ## scipy interpolation
        #_intensity_mesh[k,:] = splev(_continuum_mesh[k,:], _bsp_obj, ext=3)
        #_intensity_mesh[k,:] = _bsp_obj(_continuum_mesh[k,:])

        #intensity_mesh[k,:] = _numpy.interp(continuum_mesh[k,:], backRad[0,:], backRad[1,:])

    return intensity_mesh



def interpolate_PI_alpha_(alpha_table : T_ARRAY, alpha_table_idxs : T_ARRAY,
                          continuum_mesh : T_ARRAY) -> T_ARRAY:
    """Given continuum mesh, interpolate photoionization cross section

    only for the photoionization cross section, we use scipy cubic interpolation
    for the sake of \nu^{3} dependency

    Parameters
    ----------

    alpha_table : T_ARRAY, 2d
        concanated array of table of photoionization cross section, wavelength_cm vs alpha_cm^2
    
    alpha_table_idxs : T_ARRAY, 2d
        index array to find the array of each continuum transition in alpha_table

    continuum_mesh : T_ARRAY
        continuum mesh

    Returns
    -------

    T_ARRAY
        interpolated photoionization cross section.
    """

    alpha_mesh = _numpy.empty(continuum_mesh.shape, dtype=T_FLOAT)
    for k in range(continuum_mesh.shape[0]):
        
        i, j = alpha_table_idxs[k,:]
        alpha_table_sub = alpha_table[:, i:j]

        ## scipy cubic interpolation
        #fill_value = alpha_table[1,0], alpha_table[1,-1]
        fill_value = alpha_table_sub[1,-1], alpha_table_sub[1,0] # could be `alpha_table[1,-1],0`
        bsp_obj : _interp1d = _interp1d(x=alpha_table_sub[0,:], y=alpha_table_sub[1,:], kind="cubic",
                            bounds_error=False, fill_value=fill_value)  # no extrapolate
        alpha_mesh[k,:] = bsp_obj(continuum_mesh[k,:])

        ## numpy linear interpolation
        #_alpha_mesh[k,:] = np.interp(_continuum_mesh[k,:], _alpha_table[0,:], _alpha_table[1,:])


    return alpha_mesh


#-------------------------------------------------------------------------------
# bound free radiative trasition coefficient
#-------------------------------------------------------------------------------

def bound_free_radiative_transition_coefficient_(
                wave : T_ARRAY, 
                J : T_ARRAY, 
                alpha : T_ARRAY, 
                Te : T_UNION[T_FLOAT,T_INT], 
                nk_by_ni_LTE : T_FLOAT
                ) -> T_TUPLE[T_FLOAT, T_FLOAT, T_FLOAT]:
    """Given wavelength mesh, mean intensity (as function of wavelength),
    photoionization cross section, compute

    - radiative ionization rate,
    - stimulated radiative recombination rate,
    - spontaneous radiative recombination rate.

    Parameters
    ----------

    wave :  T_ARRAY,
        wavelength mesh, 
        [:math:`cm`]

    J :  T_ARRAY,
        mean intensity as function of wavelength, 
        [:math:`erg/cm^2/Sr/cm/s`]

    alpha :  T_ARRAY,
        photoionization cross section as function of wavelength, 
        [:math:`cm^{2}`]

    Te : T_UNION[T_FLOAT,T_INT]
        Temperature, 
        [:math:`K`]

    nk_by_ni_LTE : T_FLOAT
        population ratio in LTE, nk/ni
        [-]

    Returns
    -------
    T_TUPLE[T_FLOAT, T_FLOAT, T_FLOAT]

    Rik : T_FLOAT
        Radiative ionization rate, 
        [:math:`s^{-1} \cdot cm^{-3}`] ? [:math:`s^{-1}`]

    Rki_stim : T_FLOAT
        Stimulated Radiative ionization rate,
        [:math:`s^{-1} \cdot cm^{-3}`] ? [:math:`s^{-1}`]

    Rki_spon : T_FLOAT
        Spontaneous Radiative ionization rate,
        [:math:`s^{-1} \cdot cm^{-3}`] ? [:math:`s^{-1}`]


    Notes
    -----

    Radiative ionization rate [1]_ Equation(9.43),

        .. math:: R_{ik} = 4\pi \int_{\nu_0}^{\infty} \alpha_{ik}(\nu) (h\nu)^{-1} J(\nu) d\nu

    Stimulated radiative recombination rate [1]_ Equation(9.47),

        .. math:: R_{ki}^{stim} = \frac{n_i^{LTE}}{n_k^{LTE}} 4\pi \int_{\nu_0}^{\infty} \alpha_{ik}(\nu) (h\nu)^{-1} J(\nu) e^{-h\nu/kT} d\nu

    Spontaneous radiative recombination rate [1]_ Equation(9.45)

        .. math:: R_{ki}^{stim} = \frac{n_i^{LTE}}{n_k^{LTE}} 4\pi \int_{\nu_0}^{\infty} \alpha_{ik}(\nu) (h\nu)^{-1} (2h\nu^{3}/c^2) e^{-h\nu/kT} d\nu

    Warnings
    ---------

    the temperature in exponential term comes from the "stimulated correction", therefore it should be electron temperature.

    In case of using planck function with radiation temperature for intensity replacement, radiation temperature should only go into the planck function, not the exponential term.

    References
    -----------

    .. [1] Ivan Hubeny, Dimitri Mihalas, "Theory of Stellar Atmosphere:
        An Introduction to Astrophysical Non-equilibrium
        Quantitative Spectroscopic Analysis",
        Princeton University Press, pp. 273, 2015.
    """
    # factor : h\nu
    hv = CST.h_ * CST.c_ / wave

    expo = _numpy.exp( -hv / ( CST.k_*Te ) )

    integrand_ik = alpha * J / hv
    Rik = 4. * CST.pi_ * Integrate.trapze_(integrand_ik, wave)

    # factor : [ni/nk]_{LTE}
    #factor_ = ne * (gi/gk) / expo[-1] * Te**(-1.5) * Cst.saha_**(-1)
    factor = 1 / nk_by_ni_LTE

    integrand_ki_stim = integrand_ik * expo
    Rki_stim = factor * 4. * CST.pi_ * Integrate.trapze_(integrand_ki_stim, wave)

    #integrand_ki_spon = alpha/hv * (2.*Cst.h_*Cst.c_/wave**3) * expo
    integrand_ki_spon = alpha/hv * (2.*CST.h_*CST.c_*CST.c_/wave**5) * expo
    Rki_spon = factor * 4. * CST.pi_ * Integrate.trapze_(integrand_ki_spon, wave)

    return Rik, Rki_stim, Rki_spon


#-----------------------------------------------------------------------------
# numba optimization
#-----------------------------------------------------------------------------

if CFG._IS_JIT:

    bound_free_radiative_transition_coefficient_ = nb_njit(**NB_NJIT_KWGS) (bound_free_radiative_transition_coefficient_)