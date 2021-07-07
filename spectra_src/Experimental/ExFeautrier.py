
#-------------------------------------------------------------------------------
# experimental function/struct for feautrier method
#-------------------------------------------------------------------------------
# VERSION
#
# 0.1.0 
#    2021/07/06   u.k.   ...
#-------------------------------------------------------------------------------

from ..ImportAll import *

import numpy as np

#-------------------------------------------------------------------------------
# create tau profile along depth
#-------------------------------------------------------------------------------
def make_tau_(ND, e_max, relax=0):
    """
    make a 1D mesh of log scale from 0 to 1Ee_max with ND grid points.

    Input:
        ND: (,), number of grid points
        e_max: (,), with minvalue of 1Ee_max
        relax: (,), 0 -> 1  logspace -> linspace

    Output:
        Tau: (ND,), optical depth,
    """
    assert 0<=relax<=1, "let 0 <= relax <= 1 ."

    e_min = -4
    Tau = np.empty(ND,dtype=np.double)
    mid = e_max - np.log10(2)
    if relax ==0:
        Tau[:ND//2+1] = np.logspace(e_min,mid,ND//2+1)
    else:
        Tau[:ND//2+1] = np.logspace(e_min,mid,ND//2+1)*(1-relax) + np.linspace(10**e_min,10**mid,ND//2+1)*relax
    Tau[0] = 0
    Tau[ND//2+1:] = (-Tau[:ND//2+1][::-1]+2*Tau[ND//2])[1:]

    return Tau



#-------------------------------------------------------------------------------
# direct solution
#-------------------------------------------------------------------------------

def direct_feautrier_(tau, I_upper, I_lower, planckB, eps):
    """
    direct feautrier method which applies to only a single frequency point.
    Angle quadrature has 4 angle position.

    Input:
        tau: (ND,), optical depth, [0,ND]~[upper-lower]
        I_upper: (NA,), incident intensity at upper surface
        I_lower: (NA,), incident intensity at lower surface
        planckB: (ND,), local planck function, [0,ND]~[upper-lower]
        eps: (ND,), destruction coefficient, [0,ND]~[upper-lower]

    Output:
        S: (ND,), source function, [0,ND]~[upper-lower]
    """
    assert I_upper.size==4 and I_lower.size==4, "Angle quadrature has 4 angle position."
    ND = tau.shape[0]

    #-- angle quadrature mus and weights
    mus = np.array([0.06943184,  0.33000948,  0.66999052,  0.93056816],dtype=np.double)
    ws = np.array([0.17392742,  0.32607258,  0.32607258,  0.17392742],dtype=np.double)

    #-- array initialization
    D = np.zeros((ND,4,4),dtype=np.double)
    E = np.zeros((ND,4),dtype=np.double)
    j = np.zeros((ND,4),dtype=np.double)

    #-- dtau
    dtau_m = tau[1:] - tau[:-1]               # 0 -> ND-2 : 1/2 -> ND-3/2
    dtau = 0.5 * (dtau_m[:-1]+dtau_m[1:])     # 0 -> ND-3 : 1 -> ND-2

    A = np.zeros((4,4),dtype=np.double)
    B = np.zeros((4,4),dtype=np.double)
    C = np.zeros((4,4),dtype=np.double)
    R = np.zeros(4,dtype=np.double)

    #-- forward-elimination
    for d in range(ND):
        if d==0:
            for i in range(4):
                C[i,i] = 2*mus[i]*mus[i]/dtau_m[d]/dtau_m[d]
                B[i,:] = -(1-eps[d]) * ws[:]
                B[i,i] += 1 + 2*mus[i]/dtau_m[d] + 2*mus[i]*mus[i]/dtau_m[d]/dtau_m[d]
                R[i] = eps[d]*planckB[d] + 2*mus[i]/dtau_m[d]*I_upper[i]
            Binv = np.linalg.inv(B)
            D[d,:,:] = Binv @ C
            E[d,:] = Binv @ R

        elif d==(ND-1):
            for i in range(4):
                A[i,i] = 2*mus[i]*mus[i]/dtau_m[d-1]/dtau_m[d-1]
                B[i,:] = -(1-eps[d]) * ws[:]
                B[i,i] += 1 + 2*mus[i]/dtau_m[d-1] + 2*mus[i]*mus[i]/dtau_m[d-1]/dtau_m[d-1]
                R[i] = eps[d]*planckB[d] + 2*mus[i]/dtau_m[d-1]*I_lower[i]
            inv = np.linalg.inv( B - A @ D[d-1,:,:] )
            D[d,:,:] = 0
            E[d,:] = inv @ (R + A @ E[d-1,:] )

        else:
            for i in range(4):
                A[i,i] = mus[i]*mus[i]/dtau_m[d-1]/dtau[d-1]
                C[i,i] = mus[i]*mus[i]/dtau_m[d]/dtau[d-1]
                B[i,:] = -(1-eps[d]) * ws[:]
                B[i,i] += 1 + A[i,i] + C[i,i]
                R[i] = eps[d]*planckB[d]
            inv = np.linalg.inv( B - A @ D[d-1,:,:] )
            D[d,:,:] = inv @ C
            E[d,:] = inv @ (R + A @ E[d-1,:] )

    #-- backward-substitution
    d = ND-1; j[d,:] = E[d,:]
    for d in range(ND-2,-1,-1):
        j[d,:] = D[d,:,:] @ j[d+1,:] + E[d,:]

    #-- compute sourece function
    S = np.zeros(ND, dtype=np.double)
    for d in range(ND):
        S[d] = (1-eps[d]) * (j[d,:]*ws[:]).sum() + eps[d]*planckB[d]

    return S

#-------------------------------------------------------------------------------
# formal solution, operation
#-------------------------------------------------------------------------------

def formal_featrier_2nd_(tau, S, mu, I_upper, I_lower, method="elimination"):
    """
    Formal 2nd order feautrier method which applies to
    only single frequency point and single angle.

    Input:
        tau: (ND,), optical depth, [0,ND]~[upper-lower]
        S: (ND,), source function, [0,ND]~[upper-lower]
        mu: (1,), angle mu=cos(theta)
        I_upper: (1,), incident intensity at upper surface
        I_lower: (1,), incident intensity at lower surface
        method: "matrix", "elimination", how to solve the linear system

    Output:
        j: (ND,), 0.5*(I(+) + I(-)) at each optical depth, mean intensity like
    """

    ND = tau.shape[0]
    #-- dtau
    dtau_m = (tau[1:] - tau[:-1])/mu          # 0 -> ND-2 : 1/2 -> ND-3/2
    dtau = 0.5 * (dtau_m[:-1]+dtau_m[1:])     # 0 -> ND-3 : 1 -> ND-2

    if method=="matrix":

        T = np.zeros((ND,ND), dtype=np.double)
        R = np.empty(ND, dtype=np.double)
        #E = np.empty(ND, dtype=np.double)
        for d in range(ND):
            if d==0:
                C = 2/dtau_m[d]/dtau_m[d]
                B = 1 + 2/dtau_m[d] + C
                R[d] = S[d] + 2/dtau_m[d] * I_upper
                T[d,d] = B; T[d,d+1] = -C
            elif d==ND-1:
                A = 2/dtau_m[d-1]/dtau_m[d-1]
                B = 1 + 2/dtau_m[d-1] + A
                R[d] = S[d] + 2/dtau_m[d-1] * I_lower
                T[d,d-1] = -A; T[d,d] = B
            else:
                A = 1/dtau_m[d-1]/dtau[d-1]
                C = 1/dtau_m[d]/dtau[d-1]
                B = 1 + A + C
                R[d] = S[d]
                T[d,d-1] = -A; T[d,d] = B; T[d,d+1] = -C

        j = np.linalg.solve(T,R)

    elif method=="elimination":
        D = np.empty(ND, dtype=np.double)
        E = np.empty(ND, dtype=np.double)
        j = np.empty(ND, dtype=np.double)

        #-- forward-elimination
        for d in range(ND):
            if d==0:
                C = 2/dtau_m[d]/dtau_m[d]
                B = 1 + 2/dtau_m[d] + C
                R = S[d] + 2/dtau_m[d] * I_upper
                inv = 1/B
                D[d] = inv * C
                E[d] = inv * R
            elif d==ND-1:
                A = 2/dtau_m[d-1]/dtau_m[d-1]
                B = 1 + 2/dtau_m[d-1] + A
                R = S[d] + 2/dtau_m[d-1] * I_lower
                inv = 1/( B - A*D[d-1] )
                D[d] = 0
                E[d] = inv * (R + A * E[d-1] )
            else:
                A = 1/dtau_m[d-1]/dtau[d-1]
                C = 1/dtau_m[d]/dtau[d-1]
                B = 1 + A + C
                R = S[d]
                inv = 1/( B - A * D[d-1] )
                D[d] = inv * C
                E[d] = inv * (R + A * E[d-1] )

        #-- backward-substitution

        d = ND-1; j[d] = E[d]
        for d in range(ND-2,-1,-1):
            j[d] = D[d] * j[d+1] + E[d]

    else:
        assert False, "method keyword surports 'matrix' or 'elimination'. "

    return j

def formal_improved_RH_(tau, S, mu, r0=0, h0=0, rn=0, hn=0):
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
        <upper boundary condition: I_upper(-)= r0*I_lower(+)+ h0,   at tau[0]>
        rn: (,), 0 or 1
        hn: (,), incident intensity at lower boundary
        <lower boundary condition: I_lower(+)= rn*I_upper(-)+ hn,   at tau[nd-1]>

    Output :
        j: (ND,), 0.5*(I(+) + I(-)) at each optical depth, mean intensity like
    """
    assert r0 in (0,1), "r0 = 0 or 1"
    assert rn in (0,1), "rn = 0 or 1"

    ND = tau.shape[0]
    #-- dtau
    dtau_m = (tau[1:] - tau[:-1])/mu          # 0 -> ND-2 : 1/2 -> ND-3/2
    dtau = 0.5 * (dtau_m[:-1]+dtau_m[1:])     # 0 -> ND-3 : 1 -> ND-2

    E = np.empty(ND, dtype=np.double)
    F = np.empty(ND, dtype=np.double)
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
    j = np.empty(ND, dtype=np.double)
    j[ND-1] = E[ND-1]
    for d in range(ND-2,-1,-1):
        j[d] = (1.0+F[d])**(-1) * j[d+1] + E[d]


    return j

#-------------------------------------------------------------------------------
# formal solution, T matrix
#-------------------------------------------------------------------------------

def tmat_improved_RH_(tau, mu, r0=0, rn=0):
    """
    Purpose :
        Compute monochromatic T matrix of lambda transformation based on FormalimprovedRH method
        with:
            finite slab:     j(tau) = T(tau,mu,r0=0,rn=0) dot S + 0.5*(I_lower*exp(-(Tv-tau)/mu)+I_upper*exp(-tau/mu))

            simi-infinite:   j(tau) = T(tau,mu,r0=0,rn=1) dot S + 0.5*I_upper*exp(-tau/mu)
                        or   j(tau) = T(tau,mu,r0=0,rn=0) dot S + 0.5*(I_lower*exp(-(Tv-tau)/mu)+I_upper*exp(-tau/mu))
                             where I_lower = B+mu*dB/dtau (difussion limit)
            symmetry slab:   j(tau) = T(tau,mu,r0=0,rn=1) dot S + 0.5*I_upper*exp(-tau/mu)
                             where the deepest grid n is a the center of the slab

    Input :
        tau: (ND,), optical depth scale
        mu: (,), mu=cos(theta)
        r0: (,), 0 or 1
        rn: (,), 0 or 1
        finite: 0; infinite or midpoint of symmetry slab: 1

    Output :
        T: (ND,ND), T matrix
    """
    assert r0 in (0,1), "r0 = 0 or 1"
    assert rn in (0,1), "rn = 0 or 1"

    ND = tau.shape[0]
    #-- dtau
    dtau_m = (tau[1:] - tau[:-1])/mu          # 0 -> ND-2 : 1/2 -> ND-3/2
    dtau = 0.5 * (dtau_m[:-1]+dtau_m[1:])     # 0 -> ND-3 : 1 -> ND-2

    F = np.empty(ND, dtype=np.double)
    X_matrix = np.zeros((ND,ND), dtype=np.double)

    #-- forward-elimination
    C = 2.0/dtau_m[0]/dtau_m[0]
    H = 1.0 + ( 2.0 / dtau_m[0] ) * (1.0 -r0) / (1.0 + r0)
    F[0] = H / C
    X_matrix[0,0] = 1.0 / (H+C)
    for d in range(1,ND-1):
        A = 1.0/dtau_m[d-1]/dtau[d-1]
        C = 1.0/dtau_m[d]/dtau[d-1]
        H = 1.0
        F[d] = ( H + A*F[d-1]/(1.0 + F[d-1]) ) / C
        X_matrix[d,0:d] =  1.0 / ( C*(1.0+F[d]) ) * A * X_matrix[d-1,0:d]
        X_matrix[d,d] = 1.0 / ( C*(1.0+F[d]) )
    A = 2.0/dtau_m[ND-2]/dtau_m[ND-2]
    H = 1.0 + ( 2.0 / dtau_m[ND-2] ) * (1.0 - rn) / (1.0 + rn)
    X_matrix[ND-1,:ND-1] = 1.0 / ( H + A*F[ND-2]/(1.0+F[ND-2]) ) * A * X_matrix[ND-2,0:ND-1]
    X_matrix[ND-1,ND-1] = 1.0 / ( H + A*F[ND-2]/(1.0+F[ND-2]) )

    #-- backward-substitution
    T_matrix = np.zeros((ND,ND), dtype=np.double)
    T_matrix[ND-1,:] = X_matrix[ND-1,:]
    for d in range(ND-2,-1,-1):
        T_matrix[d,:] = (1.0+F[d])**(-1.0) * T_matrix[d+1,:] + X_matrix[d,:]

    return T_matrix

#-------------------------------------------------------------------------------
# formal solution, iteration
#-------------------------------------------------------------------------------

def angle_averaged_tmat_(TmatirxFunc, tau, I_upper=np.zeros(4,dtype=np.double), I_lower=np.zeros(4,dtype=np.double), r0=0, rn=0):
    """
    Purpose :
        Given a function to compute monochromatic T matrix of lambda transformation,
        calculate its 4 point Angle Quadrature
        with:
            finite slab:     j(tau) = T(tau,mu,r0=0,rn=0) dot S + 0.5*(I_lower*exp(-(Tv-tau)/mu)+I_upper*exp(-tau/mu))

            simi-infinite:   j(tau) = T(tau,mu,r0=0,rn=1) dot S + 0.5*I_upper*exp(-tau/mu)
                        or   j(tau) = T(tau,mu,r0=0,rn=0) dot S + 0.5*(I_lower*exp(-(Tv-tau)/mu)+I_upper*exp(-tau/mu))
                             where I_lower = B+mu*dB/dtau (difussion limit)
            symmetry slab:   j(tau) = T(tau,mu,r0=0,rn=1) dot S + 0.5*I_upper*exp(-tau/mu)
                             where the deepest grid n is a the center of the slab

    Input :
        TmatirxFunc: function to compute Tmatrix
        tau: (ND,), optical depth scale
        mu: (,), mu=cos(theta)
        r0: (,), 0 or 1
        rn: (,), 0 or 1
        I_upper: (4,), I0(-) at upper boundary.
        I_lower: (4,), In(+) at lower boundary.

    Output :
        T_average: (ND,ND), Angle averaged T matrix
        Tmats: (4,ND,ND), T matrix for 4 mu's
    """
    assert I_upper.size==4 and I_lower.size==4, "bad boundary condition."
    ND = tau.shape[0]
    mus = np.array([0.06943184,  0.33000948,  0.66999052,  0.93056816],dtype=np.double)
    ws = np.array([0.17392742,  0.32607258,  0.32607258,  0.17392742],dtype=np.double)

    T_average = np.zeros((ND,ND), dtype=np.double)

    BC = np.zeros(ND, dtype=np.double)
    for i in range(4):
        Tmat = TmatirxFunc(tau, mus[i], r0=r0, rn=rn)
        BC[0] += ws[i] * 2*mus[i]/(tau[1]-tau[0])*I_upper[i] * Tmat[:,0].sum() * (1-r0)
        BC[ND-1] += ws[i] * 2*mus[i]/(tau[ND-1]-tau[ND-2])*I_lower[i] * Tmat[:,ND-1].sum() * (1-rn)
        T_average[:,:] += ws[i] * Tmat

    return T_average, BC
