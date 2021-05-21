

#-------------------------------------------------------------------------------
# definition of functions to perform statistical equilibrium
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
# 0.0.3
#    2021/05/08   u.k.
#        - func B_Jbar : _meshInfo assigment
#        - func B_Jbar : _Gamma = numpy.atleast_1d( _Gamma )
# 0.0.2 
#    2020/11/16   u.k.
#        - Te and Tr identification in bf_R_rate() and bf_R_rate_loop()
#        - 0.5*intensity --> 1.0*intensity in the integration in B_Jbar_CRD()
# 0.0.1
#    2020/11/16   u.k.
#        - by defalut the line absorption profile is Gaussian shape so
#          no damping effect in line wings
#          ```B_Jbar()
#          _meshInfo[3] = 2
#          ```
#        - in LevelN.collisional_broadening(), ground hydrogen population is set(fixed) to 1E10
#-------------------------------------------------------------------------------

from ...ImportAll import *
from ...Atomic import LTELib as _LTELib
from ...Atomic import PhotoIonize as _PhotoIonize
from ...Atomic import Hydrogen as _Hydrogen
from ...Atomic import BasicP as _BasicP
from ...RadiativeTransfer import Profile as _Profile
from ...Util import MeshUtil as _MeshUtil
from ...Math import Integrate as _Integrate

import numpy as _numpy

#-----------------------------------------------------------------------------
# mid level functions with array as function argument
#-----------------------------------------------------------------------------

def _ni_nj_LTE_(Level : T_ARRAY, Line : T_ARRAY, Cont : T_ARRAY, 
       Te : T_FLOAT, Ne : T_FLOAT) -> T_TUPLE[T_ARRAY,T_ARRAY,T_ARRAY]:

    nLevel = Level.shape[0]
    nLine  = Line.shape[0]
    nCont  = Cont.shape[0]
    nTran  = nLine + nCont

    ## : initilize _nj_by_ni
    nj_by_ni = _numpy.empty( nTran, dtype=T_FLOAT )
    idxI     = _numpy.empty( nTran, dtype=T_INT )
    idxJ     = _numpy.empty( nTran, dtype=T_INT )

    ## : for line transitions
    gi = Line['gi'][:]
    gj = Line['gj'][:]
    Eji = CST.h_ * Line['f0'][:]

    nj_by_ni[:nLine] = _LTELib.boltzmann_distribution_(gi[:], gj[:], Eji[:], Te)
    idxI[:nLine]     = Line['idxI'][:]
    idxJ[:nLine]     = Line['idxJ'][:]

    ## : if there is continuum transition
    if nCont > 0:

        gi = Cont['gi'][:]
        gj = Cont['gj'][:]
        chi = CST.h_ * Cont['f0'][:]

        nj_by_ni[nLine:] = _LTELib.saha_distribution_ (gi[:], gj[:], chi[:], Ne, Te)
        idxI[nLine:]     = Cont['idxI'][:]
        idxJ[nLine:]     = Cont['idxJ'][:]

    isGround = Level['isGround'][:]

    ni = _nj_by_ni_To_ni_(nj_by_ni[:], idxI[:], idxJ[:], isGround[:], nLine)

    return ni, nj_by_ni[:nLine], nj_by_ni[nLine:]


def _nj_by_ni_To_ni_(nj_by_ni : T_ARRAY, idxI : T_ARRAY, idxJ : T_ARRAY, isGround : T_ARRAY, nLine : T_INT) -> T_ARRAY:

    nLevel = isGround.shape[0]
    nTran  = idxI.shape[0]

    ni = _numpy.empty(nLevel, dtype=T_FLOAT)
    ni[0] = 1.

    # loop through continuum transition first, 
    # must be sorted from lower ionization stage to higher ionization stage
    # if so, we have all ground level population relative to the very first ground level
    for k in range(nLine, nTran):
        i = idxI[k]
        j = idxJ[k]
        if isGround[ i ]:
            ni[ j ] = nj_by_ni[ k ] * ni[ i ]
    
    # loop through line transition
    for k in range(0, nLine):
        i = idxI[k]
        j = idxJ[k]
        if isGround[ i ]:
            ni[ j ] = nj_by_ni[ k ] * ni[ i ]

    return ni[:] / ni.sum(axis=0)

def _bf_R_rate_(Cont : T_ARRAY, Cont_mesh : T_ARRAY, Te : T_FLOAT, 
                nj_by_ni_Cont : T_ARRAY, alpha_interp : T_ARRAY, PI_I : T_ARRAY, 
                backRad : T_ARRAY, Tr : T_FLOAT, 
                use_Tr : T_BOOL, update_intensity : T_BOOL,
                ) -> T_TUPLE[T_ARRAY,T_ARRAY,T_ARRAY]:

    if update_intensity:
        if use_Tr:
            PI_I[:,:] = _LTELib.planck_cm_(Cont_mesh[:,:], Tr)
        else:
            PI_I[:,:] = _PhotoIonize.interpolate_PI_intensity_(backRad[:,:], Cont_mesh[:,:])
    
    #-------------------------------------------------------------------------
    # we compute/interpolate photoionizatoin cross section only once
    # and assume that while suffering Doppler shift
    #    - continuum wavelength mesh might shift
    #        (but for the sake of simplicity, we assume they do not shift)
    #    - photoionizatoin cross section keep constant
    #-------------------------------------------------------------------------
    nCont = Cont.shape[0]

    Rik      = _numpy.empty( nCont, dtype=T_FLOAT )
    Rki_stim = _numpy.empty( nCont, dtype=T_FLOAT )
    Rki_spon = _numpy.empty( nCont, dtype=T_FLOAT )
    ## loop over continuum transition
    for kL in range(nCont):
        res = _PhotoIonize.bound_free_radiative_transition_coefficient_(
                                wave = Cont_mesh[kL,::-1],
                                J = PI_I[kL,::-1],
                                alpha = alpha_interp[kL,::-1],
                                Te = Te,
                                nk_by_ni_LTE=nj_by_ni_Cont[kL])
        Rik[kL]      = res[0]
        Rki_stim[kL] = res[1]
        Rki_spon[kL] = res[2]

    return Rik, Rki_stim, Rki_spon

def _B_Jbar_with_(Line : T_ARRAY, Line_mesh_Coe : T_ARRAY,
                  Line_mesh : T_ARRAY, Line_mesh_idxs : T_ARRAY, 
                  Te : T_FLOAT, Vt : T_FLOAT, Vd : T_FLOAT, 
                  Ne : T_FLOAT, Nh_I_ground : T_FLOAT, 
                  Mass : T_FLOAT, atom_type : T_E_ATOM,
                  backRad : T_ARRAY, Tr : T_FLOAT, use_Tr : T_BOOL,
                  ) -> T_TUPLE[T_ARRAY,T_ARRAY,T_ARRAY,T_ARRAY,T_ARRAY]:

    nLine = Line.shape[0]

    wave_mesh_cm_shifted_all : T_ARRAY = _numpy.empty_like(Line_mesh)
    absorb_prof_cm_all       : T_ARRAY = _numpy.empty_like(Line_mesh)
    Jbar_all                 : T_ARRAY = _numpy.empty(nLine, dtype=T_FLOAT)

    Bji_Jbar = _numpy.empty( nLine, dtype=T_FLOAT )
    Bij_Jbar = _numpy.empty( nLine, dtype=T_FLOAT )

    absorb_prof_cm : T_ARRAY
    for k in range(nLine):

        ## collisional broadening line width
        gamma = Line["Gamma"][k]
        if atom_type == E_ATOM.HYDROGEN:
            gamma += _Hydrogen.collisional_broadening_Res_and_Van_(Line["ni"], Line["nj"], Nh_I_ground, Te)
            gamma += _Hydrogen.collisional_broadening_LinearStark_(Line["ni"], Line["nj"], Ne)
        
        Bij = Line['BIJ'][k]
        Bji = Line['BJI'][k]
        w0  = Line['w0'][k]
        f0  = Line['f0'][k]
        dopWidth_cm = _BasicP.doppler_width_(w0, Te, Vt, Mass)

        i_start, i_end = Line_mesh_idxs[k,:]
        Line_mesh[i_start:i_end]
        proftype = Line_mesh_Coe["ProfileType"][k]
        nLambda  = Line_mesh_Coe["nLambda"][k]
        qcore    = Line_mesh_Coe["qcore"][k]
        qwing    = Line_mesh_Coe["qwing"][k]
        # wm : wave mesh
        wm = _MeshUtil.make_full_line_mesh_( nLambda, qcore, qwing )
        if proftype == E_ABSORPTION_PROFILE_TYPE.VOIGT:
            dopWidth_hz = dopWidth_cm * f0 / w0
            a = gamma / ( 4. * CST.pi_ * dopWidth_hz )
            absorb_prof_cm = _Profile.voigt_(a, wm[:]) 

        elif proftype == E_ABSORPTION_PROFILE_TYPE.GAUSSIAN:
            absorb_prof_cm = _Profile.gaussian_(wm[:]) 
        
        absorb_prof_cm[:] = absorb_prof_cm[:] / dopWidth_cm # normalization
        ## -> could save to `absorb_prof_cm_all` here
        absorb_prof_cm_all[i_start:i_end] = absorb_prof_cm[:]

        wm_cm =  wm[:] * dopWidth_cm + w0
        wm_cm_shifted = wm_cm[:] + ( w0 * Vd / CST.c_ )
        ## -> could save to `wave_mesh_cm_shifted_all` here
        wave_mesh_cm_shifted_all[i_start:i_end] = wm_cm_shifted[:]

        ## : interpolate background intensity
        if use_Tr:
            # make sure this integral eqauls to 1.
            Jbar0 = _Integrate.trapze_(absorb_prof_cm[:], wm_cm_shifted[:])
            Jbar0 *= _LTELib.planck_cm_(w0, Tr)
        else:
            I_cm_interp = _numpy.interp( wm_cm_shifted[:], backRad[0,:], backRad[1,:] )
            integrand = I_cm_interp[:] * absorb_prof_cm[:]
            Jbar0 = _Integrate.trapze_(integrand[:], wm_cm_shifted[:])
        
        Jbar_all[k] = Jbar0

        Bij_Jbar[k] = Bij * Jbar0
        Bji_Jbar[k] = Bji * Jbar0

    return Bij_Jbar, Bji_Jbar, wave_mesh_cm_shifted_all, absorb_prof_cm_all, Jbar_all 


        




