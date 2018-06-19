#!/usr/local/bin/python
import sys
import os
sys.path.append(os.path.abspath('../lib'))
import gappa as gp
import numpy as np
import copy
import matplotlib.pyplot as plt


def analytical(t,pl):
    '''
        analytical solution (see e.g. Atoyan&Aharonian, 
        Mon. Not. R. Astron. Soc. 302, 253-276 (1999))
    '''
    spec = pl
    t_e = np.interp(spec[:,0],t_esc_lookup[:,0],t_esc_lookup[:,1])
    spec[:,1] *=  t_e* (1 - np.exp(-t*gp.yr_to_sec/t_e))
    spec[:,1] *= spec[:,0]**2 
    spec[:,0] /= gp.TeV_to_erg

    return np.array(spec)

if __name__ == "__main__":
    fu = gp.Utils(); fu.DrawGamera();
    mpik_green = '#057775'

    '''
        define injection spectrum
    '''
    e = np.logspace(-8,3,100) * gp.TeV_to_erg
    power_law = 1e40 * ((e/gp.TeV_to_erg)**-2)
    power_law = np.array(zip(e,power_law))
 
    '''
        define escape time lookup (here as function of energy)
    '''
    e_dash = 5e-3
    t_esc = 200 * gp.yr_to_sec
    tmax = 0.8 * t_esc / gp.yr_to_sec
    t = np.linspace(5,tmax,8)
    global t_esc_lookup
    t_esc = t_esc * (1+e/e_dash)**-0.5
    t_esc_lookup = np.array(zip(e,t_esc))
    
    '''
        plot the escape time scale spectrum
    '''

    f = plt.figure(figsize=(3,3))
    plt.loglog(e,t_esc,c=mpik_green)
    plt.ylabel("t"+r"$_{esc}$"+"(s)")
    plt.xlabel("time (yrs)")
    plt.grid(alpha=0.5)
    f.savefig("tesc_e.png",bbox_inches="tight")

    '''
        define and set up Particles object
    '''
    pa1 = gp.Particles()
    pa1.ToggleQuietMode()
    pa1.SetCustomInjectionSpectrum(power_law)
    pa1.SetEnergyDependentEscapeTime(t_esc_lookup)

    '''
        calculate proton spectra at different time steps in presence of 
        energy dependent particle escape
    '''
    p_sed = []
    for tt in t:
            pa1.SetAge(tt)
            pa1.CalculateProtonSpectrum()
            p_sed.append(np.array(pa1.GetParticleSED()))

    '''
        make a plot comparing numerical to analytical solution
    '''
    f = plt.figure(figsize=(4,4))
    for p,tt in zip(p_sed,t):
        a = analytical(tt,copy.copy(power_law))
        if tt == t[0]:
            plt.loglog(a[:,0],a[:,1],c=mpik_green,lw=3,alpha=0.6,label="analytical")
            plt.loglog(p[:,0],p[:,1],c="black",label="GAMERA")
        else:
            plt.loglog(a[:,0],a[:,1],c=mpik_green,lw=3,alpha=0.6)
            plt.loglog(p[:,0],p[:,1],c="black")
                
    plt.xlabel("energy (TeV)")
    plt.ylabel("E"+r"$^2$""dN/dE (erg)")
    plt.legend(title="age (yrs)")
    plt.legend()
    plt.grid()
    f.savefig("esc_spectra_edep.png",bbox_inches="tight")
    
