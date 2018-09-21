#!/usr/local/bin/python
import sys
import os
sys.path.append(os.path.abspath('../lib'))
import gappa as gp
import numpy as np
import matplotlib.pyplot as plt

def int_fac(x):
    '''
        integral needed in the analytical solution
    '''
    if x >= 1:
        xx = np.linspace(1,x,100)
    else:
        xx = np.linspace(x,1,100)
    yy = 1 / np.interp(xx/gp.yr_to_sec,t_esc_lookup[:,0],t_esc_lookup[:,1]) 
    return fu.Integrate(list(zip(xx,yy)))

def analytical(t):
    '''
        analytical solution for y'(t) = Q - y(t) / t_esc(t)
    '''
    t *= gp.yr_to_sec
    N0 = np.interp(1,power_law[:,0],power_law[:,1])

    xx = np.linspace(1.001,t,100)
    int_f = []
    for x in xx:
        int_f.append(int_fac(x))
    int_f = np.array(int_f)
    yy = N0 * np.exp(int_f)
    f = fu.Integrate(list(zip(xx,yy)))

    xx = np.linspace(0,0.999,100)
    int_f2 = []
    for x in xx:
        int_f2.append(-int_fac(x))
    int_f2 = np.array(int_f2)
    yy = N0 * np.exp(int_f2)
    f2 = fu.Integrate(list(zip(xx,yy)))


    return  np.exp(-int_fac(t))*(f-f2)
    


if __name__ == "__main__":
    fu = gp.Utils(); fu.DrawGamera();
    mpik_green = '#057775'

    '''
        define injection spectrum
    '''
    e = np.logspace(-5,3,100) * gp.TeV_to_erg
    power_law = 1e40 * ((e/gp.TeV_to_erg)**-2)
    power_law = np.array(list(zip(e,power_law)))
 
    '''
        define escape time
    '''
    global tmin
    tmin = 1
    t_esc = 200 * gp.yr_to_sec
    tmax = 3 * t_esc / gp.yr_to_sec
    t_dash = t_esc / 2 / gp.yr_to_sec
    t = np.linspace(2*tmin,tmax,50)
    t2 = np.linspace(tmin,tmax,1000)
    t_esc = t_esc * (np.cos(t2 / t_dash)**2+0.1)

    N0 = np.interp(1,power_law[:,0],power_law[:,1])
    analytical_no_escape = N0 * t2 * gp.yr_to_sec

    '''
        plot escape time vs. time
    '''
    f = plt.figure(figsize=(3,3))
    plt.plot(t2,t_esc,c=mpik_green)
    plt.ylabel("t"+r"$_{esc}$"+"(s)")
    plt.xlabel("time (yrs)")
    plt.grid(alpha=0.5)
    f.savefig("tst.png",bbox_inches="tight")
    global t_esc_lookup
    t_esc_lookup = np.array(list(zip(t2,t_esc)))


    '''
        define and set up Particles object
    '''
    pa1 = gp.Particles()
    pa1.ToggleQuietMode()
    pa1.SetCustomInjectionSpectrum(power_law)
    pa1.SetTimeDependentEscapeTime(t_esc_lookup)

    '''
        calculate proton spectra at different time steps in presence of 
        time dependent particle escape
    '''
    p_spec = []
    for tt in t:
            print("time "+str(round(tt,0))+"yrs")
            pa1.SetAge(tt)
            pa1.CalculateProtonSpectrum()
            p_spec.append(np.array(pa1.GetParticleSpectrum()))

    '''
        make a plot comparing numerical to analytical solution
    '''
    f = plt.figure(figsize=(4,4))
    plt.grid(alpha=0.5)
    a = []
    for tt,p in list(zip(t,p_spec)):
        N0 = np.interp(1,p[:,0],p[:,1])
        a.append(analytical(tt))
        if tt == t[0]:
            plt.scatter(tt,N0,c="black",s=10,zorder=0,label="GAMERA")
        else:
            plt.scatter(tt,N0,c="black",s=10,zorder=1)
    plt.plot(t2,analytical_no_escape,zorder=2,c="gray",label="no escape",alpha=0.5)
    plt.plot(t,a,zorder=1,c=mpik_green,label="analytical")

    plt.yscale("log")
    plt.ylabel("dN/dE[@1erg] (1/erg)")
    plt.xlabel("time (yrs)")
    plt.legend(loc="lower right")
    f.savefig("esc_amp_tdep.png",bbox_inches="tight")

