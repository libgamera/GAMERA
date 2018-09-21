#!/usr/local/bin/python
import sys
import os
sys.path.append(os.path.abspath('../lib'))
import gappa as gp
import numpy as np
import copy
import matplotlib.pyplot as plt


def t_esc(e,t):
    '''
        arbitrary function of the escape time
    '''
    e_dash = 5e-3
    t_esc_0 = 200 * gp.yr_to_sec
    t_dash = t_esc_0 / 1 / gp.yr_to_sec
    t_esc = t_esc_0*(1+e/e_dash)**(-0.5*(1-0.2*t / t_dash)) * (np.cos(t / t_dash)**2+0.1)
    return t_esc

if __name__ == "__main__":
    fu = gp.Utils(); fu.DrawGamera();
    mpik_green = '#057775'

    '''
        define injection spectrum
    '''
    e = np.logspace(-5,3,100) * gp.TeV_to_erg
    power_law = 1e40 * ((e/gp.TeV_to_erg)**-2)
    power_law = np.array(list(zip(e,power_law)))
    t_max = 1700
    t = np.linspace(1,t_max,1000)
    t2 = np.linspace(2,t_max,10)
    '''
        define and set up Particles object
    '''
    pa = gp.Particles()
    pa.ToggleQuietMode()
    pa.SetCustomInjectionSpectrum(power_law)


    '''
        create and fill the lookup holding the spectral and temporal evolution 
        of the escape time.
        order is important! Outer loop over energy, inner loop over time! 
    '''
    lookup = []
    for ee in e:
        for tt in t:
            lookup.append([tt,ee,t_esc(ee,tt)])
    pa.SetTimeAndEnergyDependentEscapeTime(lookup)
    

    '''
        another possibility, especially useful for numpy users, is to 
        give a meshgrid plus axis vectors to SetTimeAndEnergyDependentEscapeTime
    '''
    t_m, e_m = np.meshgrid(t, e)
    pa.SetTimeAndEnergyDependentEscapeTime(t, e, t_esc(e_m, t_m))

    '''
        plot the escape time vs age and energy
    '''
    f = plt.figure(figsize=(6,5))
    plt.pcolormesh(t_m,e_m,np.log10(t_esc(e_m, t_m)))
    b = plt.colorbar()
    CS = plt.contour(t_m,e_m,np.log10(t_esc(e_m, t_m)),10,cmap="plasma")
    plt.yscale("log")
    plt.ylabel("energy (erg)")
    plt.xlabel("time (yrs)")
    b.set_label("log10(escape time scale (s))")
    f.savefig("tesc_t_e.png",bbox_inches="tight")


    '''
        calculate proton spectra at different time steps in presence of 
        energy dependent particle escape
    '''
    p_sed = []
    for tt in t2:
            pa.SetAge(tt)
            pa.CalculateProtonSpectrum()
            p_sed.append(np.array(pa.GetParticleSED()))

    '''
        make a plot comparing numerical to analytical solution
    '''
    f = plt.figure(figsize=(4,4))
    for p,tt,a in list(zip(p_sed,t2,np.linspace(0.2,1,len(t2)))):
        plt.loglog(p[:,0],p[:,1],c=mpik_green,alpha=a,label=str(int(tt)))
                
    plt.xlabel("energy (TeV)")
    plt.ylabel("E"+r"$^2$""dN/dE (erg)")
    plt.legend(title="age (yrs)",ncol=2,fontsize=9,loc="lower left")
    plt.grid()
    f.savefig("esc_spectra_edep_tdep.png",bbox_inches="tight")

    '''
        calculate proton spectra at different time steps in presence of 
        energy dependent particle escape, now finer t-bins for the mesh plot
    '''
    t2 = np.linspace(2,t_max,100)
    p_sed = []
    for tt in t2:
            pa.SetAge(tt)
            pa.CalculateProtonSpectrum()
            p_sed.append(np.array(pa.GetParticleSED()))

    mesh_val = []
    for p in p_sed:
        mesh_val.append([])
        for pp in p:
            mesh_val[len(mesh_val)-1].append(pp[1])
    mesh_val = np.array(mesh_val)
    t_m, e_m = np.meshgrid(t2, p_sed[0][:,0])
    
    '''
        make the mesh plot
    '''
    f = plt.figure(figsize=(6,5))
    plt.pcolormesh(t_m,e_m,np.log10(mesh_val.T))
    b = plt.colorbar()
    CS = plt.contour(t_m,e_m,np.log10(mesh_val.T),10,cmap="plasma")
    plt.yscale("log")
    plt.ylabel("energy (TeV)")
    plt.xlabel("time (yrs)")
    b.set_label("log10(E"+r"$^2$""dN/dE (erg))")
    f.savefig("spec_t_e.png",bbox_inches="tight")
        

    
