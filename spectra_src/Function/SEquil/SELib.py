

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
from ...Atomic import Collision as _Collision
from ...Atomic import SEsolver as _SEsolver
from ...RadiativeTransfer import Profile as _Profile
from ...Util import MeshUtil as _MeshUtil
from ...Math import Integrate as _Integrate

import numpy as _numpy

#-----------------------------------------------------------------------------
# high level functions with struct as function argument
#-----------------------------------------------------------------------------

from ...Struct import Atom as _Atom
from ...Struct import Atmosphere as _Atmosphere
from ...Struct import WavelengthMesh as _WavelengthMesh
from ...Struct import Radiation as _Radiation
from ...Struct import Container as _Container

def cal_SE_with_Nh_Te_(atom : _Atom.Atom, atmos : _Atmosphere.Atmosphere0D, 
                       wMesh : _WavelengthMesh.Wavelength_Mesh, 
                       radiation : _Radiation.Radiation,
                       Nh_SE : T_UNION[T_ARRAY, None],
                       ) -> T_TUPLE[_Container.SE_Container,_Container.TranRates_Container] :
    
    Nh       = atmos.Nh                  # [/cm^{3}] 
    Ne0      = 1.E-4 * Nh                # [/cm^{3}]
    if Nh_SE is None:
        atmos.Ne = 0.5 * Nh              # [/cm^{3}]
    else:
        atmos.Ne = Ne0 + Nh * Nh_SE[0]   # [/cm^{3}]
    
    is_hydrogen = ( atom._atom_type ==  E_ATOM.HYDROGEN )
    
    while True:
        
        SE_con, tran_rate_con = cal_SE_(atom, atmos, wMesh, radiation, Nh_SE)
        
        n_SE = SE_con.n_SE
        
        if is_hydrogen:
            Ne_SE  = Ne0 + Nh * n_SE[-1]
            Ne_new = 0.5 * ( Ne_SE + atmos.Ne )
            
            if ( abs( Ne_new - atmos.Ne ) / atmos.Ne ) < 0.01:
                atmos.Ne = Ne_new
                break
            else:
                atmos.Ne = Ne_new
        else:
            break

    return SE_con, tran_rate_con

def cal_SE_with_Ne_Te_(atom : _Atom.Atom, atmos : _Atmosphere.Atmosphere0D, 
                       wMesh : _WavelengthMesh.Wavelength_Mesh, 
                       radiation : _Radiation.Radiation,
                       Nh_SE : T_UNION[T_ARRAY, None],
                       ) -> T_TUPLE[_Container.SE_Container,_Container.TranRates_Container] :
    
##    is_hydrogen = ( atom._atom_type ==  E_ATOM.HYDROGEN )

## : this comment out block tries to compute Nh with iteration
##   currently, we assume Nh does not change much in the iteration 
##   (in collisional brodenning functions), so we fix Nh
##   maybe we need this iteration when we include collisional with proton and H I
##
##    if (Nh_SE is None) & (~is_hydrogen) :
##        raise ValueError("could not perform SE for non-hydrogen atom without given `Nh_SE`.")
##        
##    if (Nh_SE is None) & (is_hydrogen):
##            atmos.Nh = 2 * atmos.Ne
##    else:
##        atmos.Nh  = atmos.Ne / ( 1.E-4 + Nh_SE[-1] )          # [cm^{-3}]

    if (Nh_SE is None):
        atmos.Nh = 2 * atmos.Ne
    else:
        atmos.Nh = atmos.Ne / ( 1.E-4 + Nh_SE[-1] )

    SE_con, tran_rate_con = cal_SE_(atom, atmos, wMesh, radiation, Nh_SE)

    if is_hydrogen := ( atom._atom_type ==  E_ATOM.HYDROGEN ):
        atmos.Nh = atmos.Ne / ( 1.E-4 + SE_con.n_SE[-1] )

    return SE_con, tran_rate_con

    
def cal_SE_(atom : _Atom.Atom, atmos : _Atmosphere.Atmosphere0D, 
            wMesh : _WavelengthMesh.Wavelength_Mesh, 
            radiation : _Radiation.Radiation,
            Nh_SE : T_UNION[T_ARRAY, None],
            ) -> T_TUPLE[_Container.SE_Container,_Container.TranRates_Container] :
    
    ## : extract variable from structs

    Mass            = atom.Mass

    atom_type       = atom._atom_type

    Level           = atom.Level
    Line            = atom.Line
    Cont            = atom.Cont

    nLevel          = atom.nLevel
    nLine           = atom.nLine
    nCont           = atom.nCont

    data_src_CE     = atom._atomic_data_source.CE
    data_src_CI     = atom._atomic_data_source.CI

    CE_Omega_table  = atom.CE.Omega_table
    CE_Te_table     = atom.CE.Te_table
    CE_Coe          = atom.CE.Coe

    CI_Omega_table  = atom.CI.Omega_table
    CI_Te_table     = atom.CI.Te_table
    CI_Coe          = atom.CI.Coe

    Cont_mesh       = wMesh.Cont_mesh

    alpha_interp    = atom.PI.alpha_interp

    PI_intensity    = radiation.PI_intensity
    backRad         = radiation.backRad

    Line_mesh_Coe   = wMesh.Line_Coe
    Line_mesh       = wMesh.Line_mesh
    Line_mesh_idxs  = wMesh.Line_mesh_idxs

    Te              = atmos.Te
    Ne              = atmos.Ne
    Vt              = atmos.Vt
    Vd              = atmos.Vd
    Tr              = atmos.Tr

    use_Tr          = atmos.use_Tr

    doppler_shift_continuum = atmos.doppler_shift_continuum
    if doppler_shift_continuum :
        raise NotImplementedError("Doppler shift of continuum wavelength mesh not yet implemented.")

    Nh_I_ground : T_FLOAT
    if Nh_SE is None:
        Nh_I_ground = 0.5 * atmos.Nh # half hydrogen atoms are in its H I ground Level
    else:
        Nh_I_ground = atmos.Nh * Nh_SE[0] / Nh_SE.sum()

    Aji             = Line["AJI"][:]

    ## : append idxI, idxJ
    idxI = _numpy.empty(nLine+nCont, dtype=DT_NB_INT)
    idxJ = _numpy.empty(nLine+nCont, dtype=DT_NB_INT)
    idxI[:nLine] = Line["idxI"][:]
    idxJ[:nLine] = Line["idxJ"][:]
    idxI[nLine:] = Cont["idxI"][:]
    idxJ[nLine:] = Cont["idxJ"][:]

    ## : Given ..., perform SE to calculate n_SE

    n_LTE , nj_by_ni = _ni_nj_LTE_(Level, Line, Cont, Te, Ne)
    #nj_by_ni_Line = nj_by_ni[:nLine]
    nj_by_ni_Cont = nj_by_ni[nLine:]

    Rik, Rki_stim, Rki_spon = _bf_R_rate_(
        Cont, Cont_mesh[:,:], Te, nj_by_ni_Cont[:], alpha_interp[:,:],
        PI_intensity[:,:], backRad[:,:], Tr, use_Tr, doppler_shift_continuum
        )

    Bij_Jbar, Bji_Jbar, wave_mesh_cm_shifted_all, absorb_prof_cm_all, Jbar_all  = \
        _B_Jbar_(
            Line, Line_mesh_Coe, Line_mesh[:], Line_mesh_idxs[:,:],
            Te, Vt, Vd, Ne, Nh_I_ground, Mass, atom_type, backRad[:,:], Tr, use_Tr
        )

    Cij = _get_Cij_(
        Line, Cont, Te, atom_type, 
        CE_Omega_table, CE_Te_table, CE_Coe, data_src_CE,
        CI_Omega_table, CI_Te_table, CI_Coe, data_src_CI
    )
    Cji = _Collision.Cij_to_Cji_(Cij[:], nj_by_ni[:])

    Rij, Rji_stim, Rji_spon = _make_Rji_Rij_(
        Aji[:], Bji_Jbar[:], Bij_Jbar[:], 
        Rki_spon[:], Rki_stim[:], Rik[:]
        )
    n_SE = _solve_SE_( 
        nLevel, idxI[:], idxJ[:], 
        Rji_spon[:], Rji_stim[:], Rij[:], 
        Cji[:], Cij[:], Ne
        )

    SE_con = _Container.SE_Container(
        n_SE = n_SE,
        n_LTE = n_LTE,
        nj_by_ni = nj_by_ni,
        wave_mesh_shifted_1d = wave_mesh_cm_shifted_all,
        absorb_prof_1d = absorb_prof_cm_all,
        Line_mesh_idxs = Line_mesh_idxs,
        Jbar = Jbar_all,
    )

    tran_rate_con = _Container.TranRates_Container(
        Rji_spon = Rji_spon[:],
        Rji_stim=Rji_stim[:],
        Rij=Rij[:],
        Cji_Ne = Cji[:] * Ne,
        Cij_Ne = Cij[:] * Ne,
    )

    return SE_con, tran_rate_con
    
    

    


#-----------------------------------------------------------------------------
# mid level functions with array as function argument
#-----------------------------------------------------------------------------

def _ni_nj_LTE_(Level : T_ARRAY, Line : T_ARRAY, Cont : T_ARRAY, 
       Te : T_FLOAT, Ne : T_FLOAT) -> T_TUPLE[T_ARRAY,T_ARRAY]:

    nLevel = Level.shape[0]
    nLine  = Line.shape[0]
    nCont  = Cont.shape[0]
    nTran  = nLine + nCont

    ## : initilize _nj_by_ni
    nj_by_ni = _numpy.empty( nTran, dtype=DT_NB_FLOAT )
    idxI     = _numpy.empty( nTran, dtype=DT_NB_INT )
    idxJ     = _numpy.empty( nTran, dtype=DT_NB_INT )

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

    return ni, nj_by_ni


def _nj_by_ni_To_ni_(nj_by_ni : T_ARRAY, idxI : T_ARRAY, idxJ : T_ARRAY, isGround : T_ARRAY, nLine : T_INT) -> T_ARRAY:

    nLevel = isGround.shape[0]
    nTran  = idxI.shape[0]

    ni = _numpy.empty(nLevel, dtype=DT_NB_FLOAT)
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
                use_Tr : T_BOOL, doppler_shift_continuum : T_BOOL,
                ) -> T_TUPLE[T_ARRAY,T_ARRAY,T_ARRAY]:

    if doppler_shift_continuum:
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

    Rik      = _numpy.empty( nCont, dtype=DT_NB_FLOAT )
    Rki_stim = _numpy.empty( nCont, dtype=DT_NB_FLOAT )
    Rki_spon = _numpy.empty( nCont, dtype=DT_NB_FLOAT )
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

def _B_Jbar_(Line : T_ARRAY, Line_mesh_Coe : T_ARRAY,
             Line_mesh : T_ARRAY, Line_mesh_idxs : T_ARRAY, 
             Te : T_FLOAT, Vt : T_FLOAT, Vd : T_FLOAT, 
             Ne : T_FLOAT, Nh_I_ground : T_FLOAT, 
             Mass : T_FLOAT, atom_type : T_E_ATOM,
             backRad : T_ARRAY, Tr : T_FLOAT, use_Tr : T_BOOL,
             ) -> T_TUPLE[T_ARRAY,T_ARRAY,T_ARRAY,T_ARRAY,T_ARRAY]:

    nLine = Line.shape[0]

    wave_mesh_cm_shifted_all : T_ARRAY = _numpy.empty_like(Line_mesh)
    absorb_prof_cm_all       : T_ARRAY = _numpy.empty_like(Line_mesh)
    Jbar_all                 : T_ARRAY = _numpy.empty(nLine, dtype=DT_NB_FLOAT)

    Bji_Jbar = _numpy.empty( nLine, dtype=DT_NB_FLOAT )
    Bij_Jbar = _numpy.empty( nLine, dtype=DT_NB_FLOAT )

    absorb_prof_cm : T_ARRAY
    for k in range(nLine):

        ## collisional broadening line width
        gamma : T_FLOAT = Line["Gamma"][k]
        if atom_type == E_ATOM.HYDROGEN:
            gamma += _Hydrogen.collisional_broadening_Res_and_Van_(Line["ni"][k], Line["nj"][k], Nh_I_ground, Te)
            gamma += _Hydrogen.collisional_broadening_LinearStark_(Line["ni"][k], Line["nj"][k], Ne)
        
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
        
        else:
            raise ValueError("Only 'VOIGT' and 'GAUSSIAN' are valid E_ABSORPTION_PROFILE_TYPE")
        
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

def _get_Cij_(Line : T_ARRAY, Cont : T_ARRAY, Te : T_FLOAT, atom_type : T_E_ATOM,
              CE_Omega_table : T_ARRAY, CE_Te_table : T_ARRAY, CE_Coe : T_ARRAY, data_src_CE : T_E_ATOMIC_DATA_SOURCE,
              CI_Omega_table : T_ARRAY, CI_Te_table : T_ARRAY, CI_Coe : T_ARRAY, data_src_CI : T_E_ATOMIC_DATA_SOURCE):
        

    nLine = Line.shape[0]
    nCont = Cont.shape[0]
    Cij = _numpy.empty(nLine+nCont, dtype=DT_NB_FLOAT)

    ## : for line transition
    if data_src_CE == E_ATOMIC_DATA_SOURCE.EXPERIMENT:

        for k in range(nLine):
            omega = _Collision.interp_omega_(CE_Omega_table[k,:],Te,CE_Te_table[:],CE_Coe["f1"][k],CE_Coe["f2"][k])
            Cij[k] = _Collision.CE_rate_coe_(omega,Te,CE_Coe["gi"][k],CE_Coe["dEij"][k])

    elif data_src_CE == E_ATOMIC_DATA_SOURCE.CALCULATE:

        if atom_type != E_ATOM.HYDROGEN:
            raise ValueError("we don't have function to calculate collisional rate coefficient for non-hydrogen atom.")
        
        for k in range(nLine):
            Cij[k] = _Hydrogen.CE_rate_coe_(Line["ni"][k],Line["nj"][k],Te)

    else:

        raise ValueError("only 'CALCULATE' and 'EXPERIMENT' are valid E_ATOMIC_DATA_SOURCE.")

    if nCont > 0:

        ## : for line transition
        if data_src_CI == E_ATOMIC_DATA_SOURCE.EXPERIMENT:

            for k in range(nCont):
                omega = _Collision.interp_omega_(CI_Omega_table[k,:],Te,CI_Te_table[:],1.,CI_Coe["f2"][k])
                Cij[k+nLine] = _Collision.CI_rate_coe_(omega,Te,CI_Coe["dEij"][k])

        elif data_src_CI == E_ATOMIC_DATA_SOURCE.CALCULATE:

            if atom_type != E_ATOM.HYDROGEN:
                raise ValueError("we don't have function to calculate collisional rate coefficient for non-hydrogen atom.")

            for k in range(nCont):
                Cij[k+nLine] = _Hydrogen.CI_rate_coe_(Cont["ni"][k],Te)

        else:

            raise ValueError("only 'CALCULATE' and 'EXPERIMENT' are valid E_ATOMIC_DATA_SOURCE.")

    else:

        raise ValueError("currently, atomic model without continuum is not yet supported.")

    return Cij

def _make_Rji_Rij_(Aji : T_ARRAY, Bji_Jbar : T_ARRAY, Bij_Jbar : T_ARRAY,
                   Rki_spon : T_ARRAY, Rki_stim : T_ARRAY, Rik : T_ARRAY
                  ) -> T_TUPLE[T_ARRAY,T_ARRAY,T_ARRAY]:
    nLine = Aji.shape[0]
    nCont = Rik.shape[0]
    nTran = nLine + nCont
    Rji_spon = _numpy.empty(nTran, dtype=DT_NB_FLOAT)
    Rji_stim = _numpy.empty(nTran, dtype=DT_NB_FLOAT)
    Rij      = _numpy.empty(nTran, dtype=DT_NB_FLOAT)

    Rji_spon[:nLine] = Aji[:]
    Rji_spon[nLine:] = Rki_spon[:]
    Rji_stim[:nLine] = Bji_Jbar[:]
    Rji_stim[nLine:] = Rki_stim[:]
    Rij[:nLine]      = Bij_Jbar[:]
    Rij[nLine:]      = Rik[:]
    
    return Rij, Rji_stim, Rji_spon

def _solve_SE_(nLevel : T_INT, idxI : T_ARRAY, idxJ : T_ARRAY, 
               Rji_spon : T_ARRAY, Rji_stim : T_ARRAY, Rij : T_ARRAY,
               Cji : T_ARRAY, Cij : T_ARRAY, Ne : T_FLOAT) -> T_ARRAY :

    Cmat = _numpy.zeros((nLevel,nLevel), dtype=DT_NB_FLOAT)
    _SEsolver.set_matrixC_(Cmat[:,:],Cji[:],Cij[:],idxI[:],idxJ,Ne)

    Rmat = _numpy.zeros((nLevel,nLevel), dtype=DT_NB_FLOAT)
    _SEsolver.set_matrixR_(Rmat[:,:], Rji_spon[:], Rji_stim[:], Rij[:], idxI[:], idxJ[:])

    n_SE = _SEsolver.solve_SE_(Rmat, Cmat)

    return n_SE



#-------------------------------------------------------------------------------
# numba optimization
#-------------------------------------------------------------------------------

if CFG._IS_JIT:

    _nj_by_ni_To_ni_ = nb_njit(**NB_NJIT_KWGS) (_nj_by_ni_To_ni_)
    _ni_nj_LTE_      = nb_njit(**NB_NJIT_KWGS) (_ni_nj_LTE_)
    _bf_R_rate_      = nb_njit(**NB_NJIT_KWGS) (_bf_R_rate_)
    _B_Jbar_         = nb_njit(**NB_NJIT_KWGS) (_B_Jbar_)
    _get_Cij_        = nb_njit(**NB_NJIT_KWGS) (_get_Cij_)
    _solve_SE_       = nb_njit(**NB_NJIT_KWGS) (_solve_SE_)