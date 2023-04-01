
#-------------------------------------------------------------------------------
# initialization of FAL model by reading data file
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/08/10   u.k.   
#-------------------------------------------------------------------------------

from ..ImportAll import *
from ..Elements import ELEMENT_DICT as _ELEMENT_DICT
from ..Struct import Atmosphere as _Atmosphere

from dataclasses import dataclass as _dataclass

import numpy as _numpy

@_dataclass(**STRUCT_KWGS_UNFROZEN)
class Population_Container:
    """
    """

    model : T_STR
    nLevel : T_INT
    n_space_dim : T_INT
    space_dim : T_TUPLE[T_INT,...]
    n_population : T_ARRAY
    n_total : T_ARRAY

def weight_per_H_():

    avg_weight = 0.0
    for name, val in _ELEMENT_DICT.items():
        avg_weight += val["Mass"] * 10**( val["Abundance"] - 12.0 )

    return avg_weight


def init_FAL_(file : T_STR):

    def _is_skip_(_s: T_STR):
        return _s[0]=="*" or len(_s) <= 1

    with open(file, "r") as f:

        # read model id
        for line in f:
            if _is_skip_(line):
                continue
            model = line.strip()
            break
        
        # read scale type
        for line in f:
            if _is_skip_(line):
                continue
            scale_str = line.strip()
            break

        # read gravity
        for line in f:
            if _is_skip_(line):
                continue
            gravity = 10**( float( line.strip() ) )
            break

        # read number of depth mesh
        for line in f:
            if _is_skip_(line):
                continue
            nDep = int( line.strip() )
            break

        # read atmosphere parameters
        dscale      = _numpy.empty(nDep, dtype=DT_NB_FLOAT)
        temperature = _numpy.empty(nDep, dtype=DT_NB_FLOAT)
        ne          =  _numpy.empty(nDep, dtype=DT_NB_FLOAT)
        vel         = _numpy.empty(nDep, dtype=DT_NB_FLOAT)
        vturb       = _numpy.empty(nDep, dtype=DT_NB_FLOAT)

        count = 0
        for line in f:
            if _is_skip_(line):
                continue
            
            dscale[count], temperature[count], ne[count], vel[count], vturb[count] = \
                [float(w) for w in line.strip().split()]
            
            count += 1
            if count == nDep:
                break
        if count < nDep:
            raise ValueError(f"only {count} point read for nDep = {nDep}")
        vel[:]   *= 1E5 # km -> cm
        vturb[:] *= 1E5 # km -> cm

        # read hydrogen population
        count = 0
        for line in f:
            if _is_skip_(line):
                continue

            if count == 0:
                words = line.strip().split()
                nLevel = len(words)
                n_pop = _numpy.zeros((nDep, nLevel), dtype=DT_NB_FLOAT)
            
            n_pop[count, :] = [float(w) for w in line.strip().split()]
            if count == nDep:
                break
            count += 1

        if count < nDep:
            raise ValueError(f"only {count} point read for nDep = {nDep}")

    # convert scale
    if scale_str == "Mass scale":
        column_mass = 10**(dscale[:])
    else:
        raise ValueError("unknown scale_str={scale_str}")

    # calculate hydrogen total number density
    n_total = n_pop.sum(axis=1)
    
    # calculate mass density and height mesh with column mass
    rho = CST.mH_ * weight_per_H_() * n_total[:]
    height = _numpy.zeros(nDep, dtype=DT_NB_FLOAT)
    for k in range(1, nDep):
        height[k] = height[k-1] - 2.0 * (column_mass[k]-column_mass[k-1]) / (rho[k-1]+rho[k])

    # construct struct
    pop_con = Population_Container(
        model = model,
        nLevel = nLevel,
        n_space_dim = 1,
        space_dim = (nDep,),#_numpy.array([nDep,]),
        n_population = n_pop[:,:] / n_total[:].reshape(-1,1),
        n_total = n_total,
    )
    atmos = _Atmosphere.AtmosphereC1D(
        model = model,
        Nh = n_total,
        Ne = ne,
        Te = temperature,
        Vd = vel,
        Vt = vturb,
        Z  = height,
        tau5 = _numpy.zeros(nDep, dtype=DT_NB_FLOAT),
        column_mass = column_mass,
        is_uniform = False,
        ndim = 1,
    )

    return pop_con, atmos



    