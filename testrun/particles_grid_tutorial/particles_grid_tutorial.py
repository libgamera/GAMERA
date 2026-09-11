#!/usr/local/bin/python

import sys
import os
sys.path.append(os.path.abspath('../../lib'))
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
    energy_in_erg_pl = np.logspace(0,2,bins) * gp.TeV_to_erg 
    alpha_pl = 2; e_total_pl = 1e37 # erg
    power_law = (energy_in_erg_pl/e_ref)**-alpha_pl
    # renormalise to e_total_pl (integrate E*dN/dE over E)
    power_law *= e_total_pl / fu.Integrate(list(zip(energy_in_erg_pl,power_law * energy_in_erg_pl)))
    # cast into a 2D array
    power_law_spectrum = np.array(list(zip(energy_in_erg_pl,power_law)))

    age = 1e4 # yrs
    cut = 1e2
    t = np.logspace(0,np.log10(age)+4,1e4)
    lum = np.array(list(zip(t,1e38*np.exp(-t**2/cut**2))))

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
    fp.SetLuminosity(lum)
    fp.SetAge(age)
    

    ''' 
        Calculate the particle spectrum (in this case electrons. for protons, type
        CalculateProtonSpectrum instead, in which case the only loss mechanism 
        applied is adiabatic expansion)
    '''
    fp.CalculateElectronSpectrum()
    sp = np.array(fp.GetParticleSpectrum())
    sed = np.array(fp.GetParticleSED())
    print(sed)


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
    print(sed_s)
  

    #### make a plot #####
    f = plt.figure(figsize=(5,5))
    plt.semilogx(sed[:,0],sed[:,1],ls="-",label="numerical")
    plt.semilogx(sed_s[:,0],sed_s[:,1],ls=":",label="semianalytic")
    plt.xlabel("E (TeV)")
    plt.ylabel("E"+r"$^2$"+"dN/dE (erg)")
    plt.grid()
    plt.legend()
    plt.xlim(xmin=1e-0,xmax=1e2)
    #plt.ylim(ymin=1e44,ymax=1e48)
    plt.ylim(ymin=0, ymax=np.max(sed[:,1]))
    plt.title("Particle SED ("+ str(round(age,0))+"yrs)")
    f.savefig("particles_grid_tutorial.png",bbox_inches='tight')

