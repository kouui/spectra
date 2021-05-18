
#-------------------------------------------------------------------------------
# function definition of naive/basic physics process
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

from ..ImportAll import *


def wave_to_freq_(wave : T_VEC_IFA) -> T_VEC_IFA:
    """convert Wavelength to Frequency.

    Parameters
    ----------
    wave : T_VEC_IFA
        Wavelength, [:math:`cm`]

    Returns
    -------
    T_VEC_IFA
        Frequency, [:math:`Hz`]
    """

    return CST.c_ / wave


def freq_to_wave_(freq : T_VEC_IFA) -> T_VEC_IFA:
    """convert Frequency to Wavelength.

    Parameters
    ----------
    freq : T_VEC_IFA
        Frequency, [:math:`Hz`]

    Returns
    -------
    T_VEC_IFA
        Wavelength, [:math:`cm`]
    """
    return CST.c_ / freq


def dop_vel_to_shift_(p0 : T_VEC_IFA, v : T_VEC_IFA) -> T_VEC_IFA:
    """given Doppler velocity and the line central wavelength/frequency,
    compute Doppler shift in wavelength/frequency.

    Parameters
    ----------
    p0 : T_VEC_IFA
        Frequency in any frequency unit or Wavelength in any length unit
    v : T_VEC_IFA
        line of sight velocity, [:math:`cm/s`]

    Returns
    -------
    T_VEC_IFA
        Frequency or Wavelength, same unit with input `p0`
    """

    return p0 * v / CST.c_


def doppler_width_(p0 : T_VEC_IFA, Te : T_VEC_IFA, 
                       Vt : T_VEC_IFA, am : T_FLOAT) -> T_VEC_IFA:
    """Given central wavelength/frequency, relative atomic mass of a line,
    and the temperature, turbulent velocity, compute the corresponding
    Doppler Width.

    Parameters
    ----------
    p0 : T_VEC_IFA
        Frequency in any frequency unit or Wavelength in any length unit
    Te : T_VEC_IFA
        Temperature, [:math:`K`]
    Vt : T_VEC_IFA
        Turbulent velocity, [:math:`cm/s`]
    am : T_FLOAT
        atomic mass relative to hydrogen atom. [-]

    Returns
    -------
    T_VEC_IFA
        [description]
    """
    eta0 = (2. * CST.k_ * Te / ( CST.mH_ * am ) + Vt * Vt )**(0.5)

    return p0 * eta0 / CST.c_


def update_level_gamma_(Aji : T_ARRAY, idxJ : T_ARRAY, 
                        gamma : T_ARRAY):
    """Given Einstein A coefficient of Levels, upper index of Lines,
    compute radiative damping constant for each Level.

    Parameters
    ----------
    Aji : T_ARRAY
        Einstein A coefficient of Lines, [:math:`s^{-1}`]
    idxJ : T_ARRAY
        index of upper level of Lines.
    gamma : T_ARRAY
        radiative damping constant of each Level, [:math:`s^{-1}`].
        An array to store output result.
    """
    gamma[:] = 0
    for i in range(Aji.shape[0]):
        gamma[idxJ[i]] += Aji[i]


def update_line_gamma_(idxI : T_ARRAY, idxJ : T_ARRAY, gamma_level : T_ARRAY, 
                       gamma_line : T_ARRAY):
    """Given Einstein radiative damping constant of Levels,
    compute radiative damping constant of Lines.

    Parameters
    ----------
    idxI : T_ARRAY
        index of lower level of Lines.
    idxJ : T_ARRAY
        index of upper level of Lines.
    gamma_level : T_ARRAY
        radiative damping constant of each Level, [:math:`s^{-1}`].
    gamma_line : T_ARRAY
        radiative damping constant of each Lines, [:math:`s^{-1}`].
        An array to store output result.
    """
    gamma_line[:] = 0
    for i in range(gamma_line.shape[0]):
        gamma_line[i] = gamma_level[idxI[i]] + gamma_level[idxJ[i]]


def damping_const_a_(gamma_line : T_VEC_FA, dop_width_hz : T_VEC_FA) -> T_VEC_FA:
    """Given the radiative damping constant and
    the Doppler Width (in frequency unit) of the line,
    compute damping constant a.

    Parameters
    ----------
    gamma_line : T_VEC_FA
        radiative damping constant of each Lines, [:math:`s^{-1}`].
    dop_width_hz : T_VEC_FA
        Doppler Width in frequency unit, [:math:`s^{-1}`].

    Returns
    -------
    T_VEC_FA
        damping constant, [-]
    """
    return gamma_line / ( 4 * CST.pi_ * dop_width_hz )


def refractive_index_in_air_(wave : T_UNION[T_FLOAT, T_INT, T_ARRAY],
                             unit : T_STR) -> T_UNION[T_FLOAT, T_INT, T_ARRAY]:
    """Given wavelength, calculate the refraction index in air

    J. Opt. Soc. Am. 62, 958 (1972)

    Parameters
    ----------
    wave : T_UNION[T_FLOAT, T_INT, T_ARRAY]
        wavelength
    unit : T_STR
        wavelength unit, "cm" or "um" or "nm" or "AA"

    Returns
    -------
    T_UNION[T_FLOAT, T_ARRAY]
        refraction index in air
    """
    fac : T_FLOAT = {
        "cm" : 1.E4,
        "um" : 1.,
        "nm" : 1E-3,
        "AA" : 1E-4,
    }[unit]

    sigma = wave * fac # unit --> um

    out = 8342.13 + 2406030 / (130 - sigma * sigma ) + 15997. / (38.9 - sigma * sigma)
    
    return  out * 1E-8 + 1.


def air_to_vacuum_(wave : T_UNION[T_FLOAT, T_INT, T_ARRAY],
                   unit : T_STR) -> T_UNION[T_FLOAT, T_ARRAY]:
    """Given the wavelength in air, compute the wavelength in vacuum

    https://physics.nist.gov/PhysRefData/ASD/Html/lineshelp.html#AIR

    Parameters
    ----------
    wave : T_UNION[T_FLOAT, T_INT, T_ARRAY]
        wavelength
    unit : T_STR
        wavelength unit, "cm" or "um" or "nm" or "AA"

    Returns
    -------
    T_UNION[T_FLOAT, T_INT, T_ARRAY]
        refraction index in air
    """
    
    return wave * refractive_index_in_air_(wave, unit)

def vacuum_to_air_(wave : T_UNION[T_FLOAT, T_INT, T_ARRAY],
                   unit : T_STR) -> T_UNION[T_FLOAT, T_ARRAY]:
    """Given the wavelength in vacuum, compute the wavelength in air

    https://physics.nist.gov/PhysRefData/ASD/Html/lineshelp.html#AIR

    Parameters
    ----------
    wave : T_UNION[T_FLOAT, T_INT, T_ARRAY]
        wavelength
    unit : T_STR
        wavelength unit, "cm" or "um" or "nm" or "AA"

    Returns
    -------
    T_UNION[T_FLOAT, T_INT, T_ARRAY]
        refraction index in air
    """
    
    return wave / refractive_index_in_air_(wave, unit)


#-----------------------------------------------------------------------------
# numba optimization
#-----------------------------------------------------------------------------

if CFG._IS_JIT:
    doppler_width_        = nb_vec( **NB_VEC_KWGS ) ( doppler_width_ )
    damping_const_a_      = nb_vec( **NB_VEC_KWGS ) ( damping_const_a_ )
    # update_level_gamma_
    # update_line_gamma_

