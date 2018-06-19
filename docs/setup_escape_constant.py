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

    '''
        define injection spectrum
    '''
    e = np.logspace(-5,3,100) * gp.TeV_to_erg
    power_law = 1e40 * ((e/gp.TeV_to_erg)**-2)
    power_law = np.array(zip(e,power_law))
 
    '''
        define escape time and some time steps
    '''
    t_esc = 200 * gp.yr_to_sec
    tmax = 3 * t_esc / gp.yr_to_sec
    t = np.linspace(2,tmax,30)
    t2 = np.linspace(1,tmax,1000)
    
    '''
        analytical solution
    '''
    N0 = np.interp(1,power_law[:,0],power_law[:,1])
    analytical = N0 * t_esc * ( 1 - np.exp(-t2*gp.yr_to_sec/t_esc) )    
    analytical_no_escape = N0 * t2 * gp.yr_to_sec

    '''
        define and set up Particles object
    '''
    pa = gp.Particles()
    pa.SetCustomInjectionSpectrum(power_law)
    pa.SetConstantEscapeTime(t_esc)

    '''
        calculate proton spectra at different time steps with 
        constant particle escape time scale
    '''
    p_spec = []
    for tt in t:
            pa.SetAge(tt)
            pa.CalculateProtonSpectrum()
            p_spec.append(np.array(pa.GetParticleSpectrum()))

    '''
        make a plot comparing numerical to analytical solution
    '''
    f = plt.figure(figsize=(4,4))
    plt.plot(t2,analytical_no_escape,c="gray",label="no escape",zorder=0,alpha=0.5)
    plt.plot(t2,analytical,c=mpik_green,label="analytical",zorder=0)
    for tt,p in zip(t,p_spec):
        N0 = np.interp(1,p[:,0],p[:,1])
        if tt == t[0]:
            plt.scatter(tt,N0,c="black",s=10,label="GAMERA",zorder=1)
        else:
            plt.scatter(tt,N0,c="black",s=10,zorder=1)
    plt.legend(loc="lower right")
    plt.grid(alpha=0.5)
    plt.ylabel("dN/dE[@1erg] (1/erg)")
    plt.xlabel("time (yrs)")
    plt.yscale("log")
    f.savefig("esc_amp_const.png",bbox_inches="tight")

