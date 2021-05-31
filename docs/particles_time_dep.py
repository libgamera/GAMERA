#!/usr/local/bin/python

import sys
import os
sys.path.append(os.path.abspath('../lib'))
import gappa as gp
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

    fu = gp.Utils(); fu.DrawGamera();
    mpik_green = '#057775'

    fp = gp.Particles()
    '''
        example: define a power-law spectrum and set it as injection spectrum. 
        Units: E(erg) - dN/dE (1/erg/s)
    '''
    # define reference energy as 1 TeV
    e_ref = gp.TeV_to_erg
    energy_in_erg_pl = np.logspace(-5,3,100) * gp.TeV_to_erg
    alpha_pl = 2; e_total_pl = 1e37 # erg
    power_law = (energy_in_erg_pl/e_ref)**-alpha_pl
    # renormalise to e_total_pl (integrate E*dN/dE over E)
    power_law *= e_total_pl / fu.Integrate(list(zip(energy_in_erg_pl,power_law * energy_in_erg_pl)))
    # cast into a 2D array
    power_law_spectrum = np.array(list(zip(energy_in_erg_pl,power_law)))
    # set it up
    fp.SetCustomInjectionSpectrum(power_law_spectrum)

    '''
        define and setup time evolution lookups of several environmental
        observables
    '''
    age = 1e5 # yrs
    age_ref = 1e3 # reference age for time evolution

    # create an array of time steps
    t_steps = np.logspace(0,np.log10(age),1000)

    # create exmaplatory time lookups 

    # -- for b-field
    b_field = 1e-4 * (t_steps / age_ref)**-0.5
    # -- for ambient density
    density = 1e4 * (t_steps / age_ref)**-0.8
    # -- for expansion velocity
    vel = 1e9 * (t_steps / age_ref)**-1
    # -- for source radius
    radius = np.cumsum(vel[:-1] * np.diff(t_steps) * gp.yr_to_sec) / gp.pc_to_cm
    radius = np.concatenate((radius,np.array([radius[-1]]))) 
    # -- for luminosity   
    lum = e_total_pl * (t_steps / age_ref)**1 * np.exp(-(t_steps / age_ref))

    '''
        set up these lookups
    '''
    fp.SetLuminosityLookup(list(zip(t_steps,lum)))
    fp.SetBField(list(zip(t_steps,b_field)))
    fp.SetAmbientDensity(list(zip(t_steps,density)))
    fp.SetExpansionVelocity(list(zip(t_steps,vel)))
    fp.SetRadius(list(zip(t_steps,radius)))

    # set CMB as radiation field
    fp.AddThermalTargetPhotons(2.7,0.255*gp.eV_to_erg)

    ''' 
        Calculate the particle spectrum (in this case electrons. for protons, type
        CalculateProtonSpectrum instead, in which case the only loss mechanism 
        applied is adiabatic expansion)
    '''
    p_seds = []
    for t in np.logspace(1,np.log10(age),9):
        fp.SetAge(t)
        fp.CalculateElectronSpectrum()
        p_seds.append(np.array(fp.GetParticleSED()))

    '''
        make plots
    '''
    alphas = np.linspace(0.2,1,9)
    f = plt.figure(figsize=(5,5))
    for s,t,a in list(zip(p_seds,np.logspace(1,np.log10(age),9),alphas)):
        plt.loglog(s[:,0],s[:,1],c=mpik_green,alpha=a,label=str(int(t)))
    plt.xlabel("E (TeV)")
    plt.ylabel("E"+r"$^2$"+"dN/dE (erg)")
    plt.grid()
    plt.xlim(xmin=1e-6,xmax=2e3)
    plt.ylim(ymin=1e38,ymax=1e49)
    plt.legend(title="Age(yrs)",ncol=3,fontsize=9)
    plt.title("Particle SED")
    f.savefig("particle_time_evolution_basic.png",bbox_inches="tight")

    f = plt.figure(figsize=(5,5))
    ls = ["--","-.",":","-"]
    plt.semilogx(t_steps,b_field/np.max(b_field),c=mpik_green,alpha=.3,ls=ls[0],label="B-field")
    plt.semilogx(t_steps,density/np.max(density),c=mpik_green,alpha=.3,ls=ls[1],label="ambient\ndensity")
    plt.semilogx(t_steps,vel/np.max(vel),c=mpik_green,alpha=.3,ls=ls[2],label="velocity")
    plt.semilogx(t_steps,radius/np.max(radius),c=mpik_green,alpha=.3,ls=ls[3],label="radius")
    plt.semilogx(t_steps,lum/np.max(lum),c=mpik_green,label="luminosity")
    plt.xlabel("time (years)")
    plt.ylabel("arb. units")
    plt.grid()
    plt.legend(fontsize=9)
    f.savefig("time_evolution_basic.png",bbox_inches="tight")
    
