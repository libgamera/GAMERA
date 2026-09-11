#!/usr/local/bin/python

import sys
import os
sys.path.append(os.path.abspath('../lib'))
import gappa as gp
import numpy as np
import math
import matplotlib.pyplot as plt
#import ConfigParser
import configparser


global lum0, age, tc, mej, e0, etab, eps, timesteps

def CalculateTimeDependentStuffNebula():
    '''
        Model parameter evolution, following Torres et al 2014, Journal of High Energy Astrophysics, Volume 1, p. 31-62.
    '''
    t = np.logspace(0,math.log10(1.e1*age),1000)
    gammap = 1.3333
    vej = math.sqrt(10.*e0/(3.*mej))
    c = math.pow((6./(15.*(gammap-1.)))+289./240.,-0.2);

    lum = (1.-etab)*lum0*(1.+t/tz)**(-1.*(brind+1.)/(brind-1.))
    emax = 3.*eps*gp.el_charge*np.sqrt(etab*lum/((1.-etab)*gp.c_speed))
    r = c*(lum0*t*gp.yr_to_sec/e0)**0.2 * vej*t*gp.yr_to_sec
    v = 1.2*r/(gp.yr_to_sec*t)
    b = np.sqrt(gp.yr_to_sec*etab*6./r**4 * np.concatenate(([0], ((lum * r)[1:] * np.diff(t)).cumsum())))
    r /= gp.pc_to_cm
    lum = np.vstack((t, lum)).T
    b = np.vstack((t, b)).T
    emax = np.vstack((t, emax)).T
    r = np.vstack((t, r)).T 
    v = np.vstack((t, v)).T

    return lum, b, emax, r, v

def broken_powerlaw(ebreak,index_low,index_high, emaxt, bins):
    '''
        broken power law used as injection spectrum
    '''
    e = np.logspace(np.log10(gp.m_e),np.log10(3*np.max(emaxt[:,1])),bins)
    n = np.zeros(len(e))
    e_low = [e<ebreak]
    e_high = [e>=ebreak] 
    n[e_low] += (e[e_low]/ebreak)**-spindlow
    n[e_high] += (e[e_high]/ebreak)**-spindhigh

    return np.array(list(zip(e,n)))

if __name__ == "__main__":

    fu = gp.Utils()
    fu.DrawGamera()

    # Read in parameter file
    configFile = os.path.abspath(sys.argv[1])
    #configParser = ConfigParser.RawConfigParser()
    configParser = configparser.RawConfigParser()
    configParser.read(configFile)

    ## read in paramters
    lum0 = configParser.getfloat('Parameters','InitialLuminosity')
    age = configParser.getfloat('Parameters','Age')
    tz = configParser.getfloat('Parameters','TauZero')
    dist = configParser.getfloat('Parameters','Distance')
    dens = configParser.getfloat('Parameters','AmbientDensity')
    tNIR = configParser.getfloat('Parameters','tNIR')
    eNIR = gp.eV_to_erg*configParser.getfloat('Parameters','edensNIR')
    tFIR = configParser.getfloat('Parameters','tFIR')
    eFIR = gp.eV_to_erg*configParser.getfloat('Parameters','edensFIR')
    ebreak = gp.m_e*configParser.getfloat('Parameters','gammabreak')
    spindlow = configParser.getfloat('Parameters','SpectralIndexLow')
    spindhigh = configParser.getfloat('Parameters','SpectralIndexHigh')
    mej = gp.mSol*configParser.getfloat('Parameters','Mej')
    e0 = configParser.getfloat('Parameters','E0')
    etab = configParser.getfloat('Parameters','etaB')
    eps = configParser.getfloat('Parameters','epsilon')
    brind = configParser.getfloat('Parameters','BrakingIndex')
    mwl = np.loadtxt(configParser.get('Files','mwl'))

    # get time-dependent model parameter evolutions
    lumt,bt,emaxt,r,v = CalculateTimeDependentStuffNebula()

    # setup particle object
    fp = gp.Particles()
    fp.SetCustomInjectionSpectrum(broken_powerlaw(ebreak,spindlow,spindhigh,emaxt,200))
    fp.SetLuminosity(lumt)
    fp.SetBField(bt)
    fp.SetEmax(emaxt)
    fp.SetRadius(r)
    fp.SetExpansionVelocity(v)
    fp.SetAmbientDensity(dens)
    fp.SetAge(age)
    fp.AddThermalTargetPhotons(2.7,0.25*gp.eV_to_erg) # CMB
    fp.AddThermalTargetPhotons(tFIR,eFIR) # far-IR
    fp.AddThermalTargetPhotons(tNIR,eNIR) # near-IR
    fp.SetTmin(2e-1 * age) # time at which iteration begins (this line is only necessary for super high B-field sources)


    # set up radiation object
    bins = 100 # number of bins of radiation spectrum
    erad = np.logspace(-20,3.5,bins) * gp.TeV_to_erg # energies(in ergs) where radiation will be calculated 
    fr = gp.Radiation()
    fr.SetDistance(dist)
    fr.AddThermalTargetPhotons(2.7,0.25*gp.eV_to_erg)
    fr.AddThermalTargetPhotons(tFIR,eFIR)
    fr.AddThermalTargetPhotons(tNIR,eNIR)
    fr.SetAmbientDensity(dens)
    fr.SetSynchrotronEmissionModel(1) # 90degree angle btw electrons and B-field

    # calculate spectra at the different time steps
    rad_spectra = []
    part_spectra = []
    bins = 100 # number of bins of electron spectrum
    fp.SetAge(age) # set age of system
        
    # calculate electron spectrum 
    fp.CalculateElectronSpectrum(bins)
    el_spec = np.array(fp.GetParticleSpectrum())
    el_sed = np.array(fp.GetParticleSED())

    fr.SetElectrons(el_spec)
    fr.SetBField(fp.GetBField())
    fr.AddSSCTargetPhotons(fp.GetRadius())

    # set up parameters required to calculate radiation spectrum

    # whole electron energy range
    fr.CalculateDifferentialPhotonSpectrum(erad)
    tot = np.array(fr.GetTotalSED())
    ic = np.array(fr.GetICSED())
    ic_cmb = np.array(fr.GetICSED(0))
    ic_fir = np.array(fr.GetICSED(1))
    ic_nir = np.array(fr.GetICSED(2))
    ic_ssc = np.array(fr.GetICSED(3))
    synch = np.array(fr.GetSynchrotronSED())
    brems = np.array(fr.GetBremsstrahlungSED())


    f = plt.figure(figsize=(7,4))
    plt.loglog(tot[:,0],tot[:,1],c="black",lw=3,alpha=0.3)
    plt.loglog(ic[:,0],ic[:,1],c="red",lw=1,alpha=0.8,label="IC")
    plt.loglog(ic_cmb[:,0],ic_cmb[:,1],c="gray",lw=1,alpha=0.8,label="IC-CMB")
    plt.loglog(ic_fir[:,0],ic_fir[:,1],c="gray",lw=1,alpha=0.8,ls=":",label="IC-FIR")
    plt.loglog(ic_nir[:,0],ic_nir[:,1],c="gray",lw=1,alpha=0.8,ls="-.",label="IC-NIR")
    plt.loglog(ic_ssc[:,0],ic_ssc[:,1],c="gray",lw=1,alpha=0.8,ls="--",label="IC-SSC")
    plt.loglog(brems[:,0],brems[:,1],c="cyan",lw=1,alpha=0.8,label="Brems.")
    plt.loglog(synch[:,0],synch[:,1],c="blue",lw=1,alpha=0.8,label="Synch.")

    plt.legend(ncol=2)
    plt.scatter(mwl[:,0]/1e12,mwl[:,3],c="black",alpha=1)

    plt.title("current time spectrum")
    plt.ylim(ymin=1e-5*np.max(tot[:,1]),ymax=10*np.max(tot[:,1]))
    plt.xlabel("E (TeV)")
    plt.ylabel("N"+r"$^2$"+"dN/dE (erg/cm"+r"$^2$"+"/s)")

    plt.grid()
    f.savefig("sed.png",bbox_inches="tight")

    f = plt.figure(figsize=(5,5))
    plt.loglog(el_sed[:,0],el_sed[:,1],c="black",lw=2,alpha=1)
    plt.title("particle spectrum")
    plt.xlabel("E (TeV)")
    plt.ylabel("N"+r"$^2$"+"dN/dE (erg)")
    plt.grid()
    f.savefig("el_sed.png",bbox_inches="tight")



