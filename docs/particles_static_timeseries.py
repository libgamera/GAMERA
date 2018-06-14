#!/usr/local/bin/python

import sys
import os
sys.path.append(os.path.abspath('path/to/your/lib/directory'))
import gappa as gp
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

    fu = gp.Utils()
    fu.DrawGamera()

    bins = 500

    # define reference energy as 1 TeV
    e_ref = gp.TeV_to_erg

    '''
        example: define a power-law spectrum. Units:
        E(erg) - dN/dE (1/erg/s)
    '''
    logemin = -5; logemax = 3
    # it is good to add a margin at the upper end of the energy range. That is,
    # due to the algorithm, at the highest bins a downwards bias can occur. 
    energy_in_erg_pl = np.logspace(logemin,logemax,bins) * gp.TeV_to_erg 
    alpha_pl = 2; e_total_pl = 1e37 # erg
    power_law = (energy_in_erg_pl/e_ref)**-alpha_pl
    # renormalise to e_total_pl (integrate E*dN/dE over E)
    power_law *= e_total_pl / fu.Integrate(zip(energy_in_erg_pl,power_law * energy_in_erg_pl))
    # cast into a 2D array
    power_law_spectrum = np.array(zip(energy_in_erg_pl,power_law))


    fp = gp.Particles()

    '''
        define loss terms:
        - Synchrotron losses     <-> B-Field strength
        - Inverse-Compton losses <-> Radiation fields
        - Bremsstrahlung losses
    '''
    b_field = 1e-5 # Gauss
    density = 5 # 1/cm^3
    age = 1e4 # yrs
    distance = 1e3 # pc

    fp.SetBField(b_field)
    fp.SetAmbientDensity(density)

    ''' 
        the following are the parameters of the CMB
        See also the tutorial on how to set target fields in the GAMERA docu!
    '''
    t_cmb = 2.7; edens_cmb = 0.4*gp.eV_to_erg

    # set the fields in the Radiation object
    fp.AddThermalTargetPhotons(t_cmb,edens_cmb,bins)
    fp.SetCustomInjectionSpectrum(power_law_spectrum)

    # set up radiation object
    fr = gp.Radiation()
    fr.AddArbitraryTargetPhotons(fp.GetTargetPhotons()) # output from 'Particles' is in the right format to be used in 'Radiation'
    fr.SetBField(b_field)
    fr.SetAmbientDensity(density)
    fr.SetDistance(distance) # optional, in parsec. If not set or equals zero, luminosity instead of flux will be calculated


    time_steps = np.logspace(1,6,10)

    ''' 
        Calculate the electron and radiation spectrum at different time steps
    '''
    rad = []
    part = []
    for t in time_steps:
        fp.SetAge(t)
        fp.CalculateElectronSpectrum()
        
        sp = np.array(fp.GetParticleSpectrum())
        part.append(np.array(fp.GetParticleSED()))

        fr.SetElectrons(sp)
        fr.CalculateDifferentialPhotonSpectrum(np.logspace(-19,3,bins/2) * gp.TeV_to_erg)
        rad.append(np.array(fr.GetTotalSED()))
        
    


    ''' 
        Calculate the radiation spectrum from the previously calculated electron distribution
    ''' 

    #### make a plot #####
    f, (ax1, ax2) = plt.subplots(1, 2,figsize=(12,6))
    ax1.set_prop_cycle('color',plt.get_cmap('plasma_r')(np.linspace(0.2, .8, len(time_steps))))  #
    for p,t in zip(part,time_steps):
        ax1.loglog(p[:,0],p[:,1],label=str(int(t)))
    ax1.set_xlabel("E (TeV)")
    ax1.set_ylabel("E"+r"$^2$"+"dN/dE (erg)")
    ax1.grid()
    ax1.legend(ncol=2,prop={'size':10},title="age(yrs):")
    ax1.set_xlim(xmin=1e-6,xmax=2e3)
    ax1.set_ylim(ymin=1e40,ymax=1e52)
    ax1.set_title("Particle SEDs")

    ax2.set_prop_cycle('color',plt.get_cmap('plasma_r')(np.linspace(0.2, .8, len(time_steps))))  #
    for r,t in zip(rad,time_steps):
        ax2.loglog(r[:,0],r[:,1],label=str(int(t)))
    ax2.set_xlabel("E (TeV)")
    ax2.set_ylabel("E"+r"$^2$"+"dN/dE (erg/cm^2/s)")
    ax2.legend(ncol=2,prop={'size':10},title="age(yrs):")
    ax2.grid()
    ax2.set_xlim(xmin=1e-19,xmax=1e4)
    ax2.set_ylim(ymin=1e-15,ymax=1e-6)
    ax2.set_title("Radiation SEDs")
    f.savefig("particles_static_timeseries.png",bbox_inches='tight')
