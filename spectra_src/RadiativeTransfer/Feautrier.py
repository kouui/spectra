#-------------------------------------------------------------------------------
# function/struct for feautrier method
#-------------------------------------------------------------------------------
# VERSION
#
# 0.1.0 
#    2021/07/06   u.k.   ...
#-------------------------------------------------------------------------------

from ..ImportAll import *

import numpy as _numpy



def formal_improved_RH_(tau : T_ARRAY, S: T_ARRAY, mu : T_FLOAT, r0 : T_FLOAT, h0: T_FLOAT, rn: T_FLOAT, hn: T_FLOAT):
    """
    Purpose :
        Evaluate monochromatic intensities j = 0.5*(I(+)+I(-))
        along a ray with given optical depth scale tau and source function.
        [G.B. Rybicki & D.G Hummer, A&A 245, 171-181]
        with:
            finite slab :
                I0(-) = 0             --> r0 = h0 = 0
                In(+) = I_incident(+) --> rn = 0, hn = I_incident(+)
            semi-finite :
                In(+) = In(-)         --> hn = 0, rn = 1
                or
                In(+)=B+mu*dB/dtau    --> rn = 0, hn = B+mu*dB/dtau (diffusion limit)
            illuminated medium :
                I0(-) = I_ill         --> r0 = 0, h0 = I_ill
            symmetrical slab :
                I(+) = I(-)           --> hn = 0, rn = 1, last grid locates at the center of the slab

    Input :
        tau: (ND,), optical depth scale
        source: (ND,), monochromatic source function
        mu: (,), mu=cos(theta)
        r0: (,), 0 or 1
        h0: (,), incident intensity at upper boundary
        <upper boundary condition: I_upper(-)= r0*I_upper(+)+ h0,   at tau[0]>
        rn: (,), 0 or 1
        hn: (,), incident intensity at lower boundary
        <lower boundary condition: I_lower(+)= rn*I_lower(-)+ hn,   at tau[nd-1]>

    Output :
        j: (ND,), 0.5*(I(+) + I(-)) at each optical depth, mean intensity like
    """

    ND = tau.shape[0]
    #-- dtau
    dtau_m = (tau[1:] - tau[:-1])/mu          # 0 -> ND-2 : 1/2 -> ND-3/2
    dtau = 0.5 * (dtau_m[:-1]+dtau_m[1:])     # 0 -> ND-3 : 1 -> ND-2

    E = _numpy.empty(ND, dtype=DT_NB_FLOAT)
    F = _numpy.empty(ND, dtype=DT_NB_FLOAT)
    #-- forward-elimination
    C = 2.0/dtau_m[0]/dtau_m[0]
    H = 1.0 + ( 2.0 / dtau_m[0] ) * (1.0 -r0) / (1.0 + r0)
    R = S[0] + 2.0*h0 / ( (1.0+r0)*dtau_m[0] )
    F[0] = H / C
    E[0] = R / (H+C)
    for d in range(1,ND-1):
        A = 1.0/dtau_m[d-1]/dtau[d-1]
        C = 1.0/dtau_m[d]/dtau[d-1]
        H = 1.0
        R = S[d]
        F[d] = ( H + A*F[d-1]/(1.0 + F[d-1]) ) / C
        E[d] = ( R + A*E[d-1] ) / ( C * (1.0 + F[d]) )
    A = 2.0/dtau_m[ND-2]/dtau_m[ND-2]
    H = 1.0 + ( 2.0 / dtau_m[ND-2] ) * (1.0 - rn) / (1.0 + rn)
    R = S[ND-1] + 2.0*hn / ((1.0 + rn)*dtau_m[ND-2])
    E[ND-1] = ( R + A*E[ND-2] ) / ( H + A*(F[ND-2]/(1.0+F[ND-2])) )

    #-- backward-substitution
    j = _numpy.empty(ND, dtype=DT_NB_FLOAT)
    j[ND-1] = E[ND-1]
    for d in range(ND-2,-1,-1):
        j[d] = (1.0+F[d])**(-1) * j[d+1] + E[d]


    return j

#-----------------------------------------------------------------------------
# numba optimization
#-----------------------------------------------------------------------------

if CFG._IS_JIT:
    formal_improved_RH_ = nb_njit(**NB_NJIT_KWGS) (formal_improved_RH_)

else:
    pass