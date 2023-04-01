
import numpy as np

from spectra_src.ImportAll import *
import warnings
warnings.simplefilter('ignore', UserWarning)
warnings.simplefilter('ignore', DeprecationWarning)

from spectra_src.Util import HelpUtil

from spectra_src.Struct import Atom as Atom
from spectra_src.Experimental import ExFAL as _ExFAL
from spectra_src.Experimental import ExSpectrum as _ExSpectrum
from spectra_src.Experimental import ExScatter as _ExScatter
from spectra_src.Atomic import BasicP as _BasicP
from spectra_src.RadiativeTransfer import Profile as _Profile
from spectra_src.Atomic import Hydrogen as _Hydrogen
from spectra_src.Atomic import LTELib as _LTELib

from spectra_src.Function.SEquil import SELib as _SELib
from spectra_src.Math import GaussLeg 
from spectra_src.Experimental import ExFeautrier2 as _ExFeautrier

import os


def main():
#
    conf_path = os.path.join( CFG._ROOT_DIR, "data/conf/H6.conf" )
    atom, wMesh, path_dict = Atom.init_Atom_(conf_path , is_hydrogen=True)

    FAL_path = os.path.join( CFG._ROOT_DIR, "data/atmos/FAL/FALC_82.atmos" )
    pop_con, atmos = _ExFAL.init_FAL_(FAL_path)
    ##: TODO : interpolate FAL atmosphere parameters to an symmetric log z mesh
    spectrum = _ExSpectrum.init_spectrum_(atom=atom, wMesh=wMesh, Te=atmos.Te.min(), Vt=atmos.Vt.min(), Vd=atmos.Vd[0] )
    nspec = spectrum.nSpectrum
    spec1d= spectrum.spectrum

    nZ = atmos.Z.size

    # pre-calculate LTE population
    nH_pop_LTE = np.empty((nZ, atom.nLevel), dtype=DT_NB_FLOAT)
    for kz in range(nZ):
        nH_pop_LTE[kz,:], _ = _SELib._ni_nj_LTE_(atom.Level, atom.Line, atom.Cont, atmos.Te[kz], atmos.Ne[kz])

    # background opacity, emissivity and scattering coefficient
    alp_c = np.empty((nspec,nZ), dtype=DT_NB_FLOAT)
    eta_c = np.empty((nspec,nZ), dtype=DT_NB_FLOAT)
    eps_c = np.empty((nspec,nZ), dtype=DT_NB_FLOAT)
    for ks in range(nspec):
        #print(f"\r{ks}/{spectrum.nSpectrum}", end="")
        wave = spec1d[ks]
        
        # pre-calculate background opacity
        for kz in range(nZ):

            nH_pop = pop_con.n_population[kz,:]
            nH = atmos.Nh[kz]
            Te = atmos.Te[kz]
            ne = atmos.Ne[kz]
            alp_c[ks,kz], eta_c[ks,kz], eps_c[ks,kz] = _ExScatter.background_opacity_(wave, nH_pop_LTE[kz,:], nH_pop, nH, ne, Te)


    # (base opacity)  background : thomson, rayleigh, H-, H2+, H+-ff, (LTE)... 
    # (active opacity)radiative transfer : (atomic-model) bound-bound (line), bound-free (continuum)

    # calculate transition (radiative transfer) opacity and emissivity (bound-bound + bound-free)
    Line = atom.Line
    Level = atom.Level
    Cont = atom.Cont


    # calculate bound-free cross section
    ##: TODO: if not hydrogen, then interpolate from data?
    cross_sec = np.zeros((atom.nCont, nspec))
    for ks in range(nspec):
        wave = spectrum.spectrum[ks]
        cross_sec[:,ks] = _Hydrogen.PI_cross_section_cm_(Cont["ni"][:], wave, 1)

    #idx = np.where(np.logical_and((spectrum.spectrum<6565*1E-8),(spectrum.spectrum>6560*1E-8)))
    #print(idx)
    #print(spectrum.spectrum[idx])
    #assert False, ""

    alp_rt = np.zeros((nspec,nZ), dtype=DT_NB_FLOAT)
    eta_rt = np.zeros((nspec,nZ), dtype=DT_NB_FLOAT)

    fourpi = 4. * CST.pi_

    #spec_dop_width_cm = spectrum.dop_width_cm   ##: not used
    spec_line_mesh_dop = spectrum.line_mesh_dop

    for ks in range(nspec):
    #    print(f"\r{ks}/{spectrum.nSpectrum}", end="")
        wave = spec1d[ks] # [cm]

        ##: line alp and eta
        for kl in spectrum.belonging_line[ks,:]:
            if kl == -1:
                break
            Aji = Line["AJI"][kl]
            Bji = Line["BJI"][kl]
            Bij = Line["BIJ"][kl]
            f0 = Line["f0"][kl]
            hv = CST.h_ * f0
            w0 = Line["w0"][kl]
            
            dopWidth_cm = _BasicP.doppler_width_(w0, atmos.Te[:], atmos.Vt[:], atom.Mass)
            dopWidth_hz = dopWidth_cm * f0 / w0
            ni = pop_con.n_population[ :,Line["idxI"][kl] ].copy() # [-]
            nj = pop_con.n_population[ :,Line["idxJ"][kl] ].copy() # [-]
            alp0 = hv / fourpi * (Bij * ni[:] - Bji * nj[:] ) * pop_con.n_total[:]
            a = Line["Gamma"][kl] / ( fourpi * dopWidth_hz[:] )
            ##: TODO: add collisional broadening
            
            #wm =  (wave - w0) / spectrum.dop_width_cm[kl]
            wm  = spec_line_mesh_dop[ks]
            phi = _Profile.voigt_(a, wm) / dopWidth_cm
            
            alp_rt[ks,:] += alp0 * phi
            eta_rt[ks,:] += hv / fourpi * ( Aji * nj ) * phi * pop_con.n_total[:]
        del kl
        ## cont alp and eta
        for kc in spectrum.belonging_cont[ks,:]:
            if kc == -1:
                break
            f0 = CST.c_ / wave
            ni = pop_con.n_population[ :,Cont["idxI"][kc] ]
            nj = pop_con.n_population[ :,Cont["idxJ"][kc] ]
            ni_LTE = nH_pop_LTE[:,Cont["idxI"][kc]]
            nj_LTE = nH_pop_LTE[:,Cont["idxJ"][kc]]
            bi = ni / ni_LTE
            bj = nj / nj_LTE
            exp_hnu_kT = np.exp(CST.h_*f0/CST.k_/atmos.Te[:])
            
            src_bf = 2*CST.h_*(CST.c_*CST.c_) / (wave)**5 / (bi/bj*exp_hnu_kT-1.)
            alp_bf = ni * cross_sec[kc,ks] * (1. - bj/bi/exp_hnu_kT)
            eta_bf = src_bf * alp_bf
            alp_rt[ks,:] += alp_bf
            eta_rt[ks,:] += eta_bf
            
        del kc

    ##: alp_bf ~ 1E-23
    ##: alp0   ~ 1E-9
    ##: is it correct to have such large difference

    ##: debug
    # import matplotlib.pyplot as plt
    # fig, ax = plt.subplots(1,1)
    # ax.plot(alp_rt[551:700,41],'.',markersize=1)
    # plt.show(block=False)
    # #ax.plot(spec1d[551:700]*1e8, alp_rt[551:700,41],'.',markersize=1)
    # #plt.show(block=False)
    # breakpoint()


    # angles

    nray = 4

    mus, weights = GaussLeg.gauss_quad_coe_(0.,1.,nray)

    # calculate total extinction and emmisivity
    alp_total = alp_rt + alp_c
    #src_total = (eta_rt + eta_c + eps_c * Jbar) / (alp_total) # compute in iteration

    # optical depth
    dtau = np.zeros((nspec,nZ), dtype=DT_NB_FLOAT)

    for ks in range(nspec):
        wave = spectrum.spectrum[ks]
        dtau[ks,0]  = 0.5 * alp_total[ks,0] * (atmos.Z[0]-atmos.Z[1])
        dtau[ks,1:] = 0.5 * (alp_total[ks,:-1]+alp_total[ks,1:]) * (atmos.Z[:-1]-atmos.Z[1:])

    # compute tau
    tau = dtau.cumsum(axis=1)

    Jbar = np.zeros((nspec, nZ), dtype=DT_NB_FLOAT)

    # solve radiative transfer for (ray, wavelength)
    mask = np.ones(nspec, dtype="uint8")
    niter = 0

    ##-------------------------------------------------------------------------------------------
    ## inconsistency : population has to be computed from (fit to) some mean intensity J
    ## since we already know the source function, perhaps we don't need the iteration to update source function
    ##-------------------------------------------------------------------------------------------

    while mask.sum()>0:
        print(f"niter = {niter}")
        for ks in range(spectrum.nSpectrum):
            if ks != 900:
                continue
            
            if (mask[ks]==0): 
                continue

            #print(f"updating ks = {ks}")

            src_total_ks = (eta_rt[ks,:] + eta_c[ks,:] + eps_c[ks,:] * Jbar[ks,:]) / (alp_total[ks,:])
            breakpoint()
            
            Pnew = np.zeros(nZ, dtype=DT_NB_FLOAT)
            for kr in range(nray):
                mu = mus[kr]
                weight = weights[kr]
                #I_l = I_lower[kr, ks]
                #wave = spectrum.spectrum[ks]
                Pnew[:] += weight * _ExFeautrier.formal_improved_RH_(tau[ks,:], src_total_ks, mu, 0, 0, 1, 0) ##: r0, h0, rn, hn
            print(Pnew)
            dJmax_ratio = np.abs((Jbar[ks,:] - Pnew[:]) / Pnew[:]).max()
            print(dJmax_ratio)
            Jbar[ks,:] = Pnew[:]
            if dJmax_ratio < 1E-2:
                mask[ks] = 0
                print(f"convergence of frequency no {ks}")

        niter += 1

    return 0

if __name__ == '__main__':
    main()


