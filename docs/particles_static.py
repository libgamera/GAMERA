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
    energy_in_erg_pl = np.logspace(-5,3,bins) * gp.TeV_to_erg 
    alpha_pl = 2; e_total_pl = 1e37 # erg
    power_law = (energy_in_erg_pl/e_ref)**-alpha_pl
    # renormalise to e_total_pl (integrate E*dN/dE over E)
    power_law *= e_total_pl / fu.Integrate(list(zip(energy_in_erg_pl,power_law * energy_in_erg_pl)))
    # cast into a 2D array
    power_law_spectrum = np.array(list(zip(energy_in_erg_pl,power_law)))


    fp = gp.Particles()

    '''
        define loss terms:
        - Synchrotron losses     <-> B-Field strength
        - Inverse-Compton losses <-> Radiation fields
        - Bremsstrahlung losses
    '''
    b_field = 1e-5 # Gauss
    density = 5. # 1/cm^3
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
    fp.SetAge(age)

    ''' 
        Calculate the particle spectrum (in this case electrons. for protons, type
        CalculateProtonSpectrum instead, in which case the only loss mechanism 
        applied is adiabatic expansion)
    '''
    fp.CalculateElectronSpectrum()
    #fp.CalculateProtonSpectrum() # uncomment, if you wish to calculate protons instead
    sp = np.array(fp.GetParticleSpectrum())
    sed = np.array(fp.GetParticleSED())


    '''
        In this example, energy losses are constant and no particle escape has 
        been set. You can therefore also use the semi-analytic approach, which
        is much faster when loss rates are hight (e.g. in stron radiation fields
        or B-fields)
    '''
    fp.SetSolverMethod(1)
    fp.CalculateElectronSpectrum()
    sp_s = np.array(fp.GetParticleSpectrum())
    sed_s = np.array(fp.GetParticleSED())

    ''' 
        Calculate the radiation spectrum from the previously calculated electron distribution
    ''' 

    fr = gp.Radiation()
    fr.AddArbitraryTargetPhotons(fp.GetTargetPhotons()) # output from 'Particles' is in the right format to be used in 'Radiation'
    fr.SetBField(b_field)
    fr.SetAmbientDensity(density)
    fr.SetDistance(distance) # optional, in parsec. If not set or equals zero, luminosity instead of flux will be calculated
    fr.SetElectrons(sp)
    fr.CalculateDifferentialPhotonSpectrum(np.logspace(-19,3,bins/2) * gp.TeV_to_erg)

    tot = np.array(fr.GetTotalSED())
    ic = np.array(fr.GetICSED())
    brems = np.array(fr.GetBremsstrahlungSED())
    synch = np.array(fr.GetSynchrotronSED())

    #### make a plot #####
    f, (ax1, ax2) = plt.subplots(1, 2,figsize=(12,6))
    ax1.set_prop_cycle('color',plt.get_cmap('plasma_r')(np.linspace(0.2, .8, 2)))  #
    ax1.loglog(sed[:,0],sed[:,1],ls="-",label="numerical")
    ax1.loglog(sed_s[:,0],sed_s[:,1],ls=":",label="semianalytic")
    ax1.set_xlabel("E (TeV)")
    ax1.set_ylabel("E"+r"$^2$"+"dN/dE (erg)")
    ax1.grid()
    ax1.legend()
    ax1.set_xlim(xmin=1e-6,xmax=2e3)
    ax1.set_ylim(ymin=1e44,ymax=1e48)
    ax1.set_title("Particle SED ("+ str(round(age,0))+"yrs)")

    ax2.set_prop_cycle('color',plt.get_cmap('plasma')(np.linspace(0., .8, 4)))  #
    ax2.loglog(tot[:,0],tot[:,1],ls="-",label="total",lw=3,alpha=0.3)
    ax2.loglog(synch[:,0],synch[:,1],ls=":",label="synch")
    ax2.loglog(brems[:,0],brems[:,1],ls="-.",label="brems")
    ax2.loglog(ic[:,0],ic[:,1],linestyle="--",label="IC")
    ax2.set_xlabel("E (TeV)")
    ax2.set_ylabel("E"+r"$^2$"+"dN/dE (erg/cm^2/s)")
    ax2.legend(ncol=2,prop={'size':12})
    ax2.grid()
    ax2.set_xlim(xmin=1e-19,xmax=1e4)
    ax2.set_ylim(ymin=1e-15,ymax=1e-8)
    ax2.set_title("Radiation SED ("+ str(round(age,0))+"yrs)")
    f.savefig("particles_static.png",bbox_inches='tight')
