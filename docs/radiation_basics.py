#!/usr/local/bin/python

import sys
sys.path.append('path/to/your/lib/directory')
import gappa as gp
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":


    '''
        define arbitrary particle spectra. they need to be 
        2-dimensional arrays (or vectors if using C++) with each
        entry holding a tuple: E(erg) - dN/dE (1/erg). Let's define a 
        logparabola and a power-law as example spectra. Of course the 
        output of the Particle class is designed to use as an input for
        the Radiation class. This is explained in another tutorial.
    '''

    # define reference energy as 1 TeV
    e_ref = gp.TeV_to_erg

    '''
        example: define a power-law spectrum
    '''
    # define the energy range of the particles (total energy), here between 1 GeV and 1 PeV
    energy_in_erg_pl = np.logspace(-3,3,200) * gp.TeV_to_erg

    # spectral parameters for the power-law: index and total energy in particles
    alpha_pl = 2; e_total_pl = 1e49 # erg

    # make the power-law
    power_law = energy_in_erg_pl**-alpha_pl

    # renormalise it to e_total_protons using the Integrator (GSL) 
    # in the Utils class
    fu = gp.Utils()
    power_law *= e_total_pl / fu.Integrate(zip(energy_in_erg_pl,power_law * energy_in_erg_pl))

    # zip into a 2D-list
    power_law_spectrum = np.array(zip(energy_in_erg_pl,power_law))


    '''
        example: define a logparabola
    '''
    energy_in_erg_logpar = np.logspace(-5,3,200) * gp.TeV_to_erg

    # spectral parameters: index(alpha_logpar), curvature(beta_logpar), total energy in particles
    alpha_logpar = 2.3; beta_logpar = 0.3; e_total_logpar = 1e46

    # calculate logparabola
    logparabola = (energy_in_erg_logpar/e_ref)**-(alpha_logpar + beta_logpar * np.log10((energy_in_erg_logpar/e_ref)))

    logparabola *= e_total_logpar / fu.Integrate(zip(energy_in_erg_logpar,energy_in_erg_logpar * logparabola))

    logparabola_spectrum = np.array(zip(energy_in_erg_logpar,logparabola))
    

    '''
        now set the parameters needed for the calculation 
        of the photon emission
    '''
    b_field = 1e-5 # Gauss
    ambient_density = 1 # 1/cm^3
    # the following are the parameters of the CMB
    t_cmb = 2.7; edens_cmb = 0.25 * gp.eV_to_erg #erg
    distance = 1e3 # optional, in parsec. If not set or equals zero, differential 
                   # photon production rate instead of flux will be calculated

    '''
        create a Radiation object and set it up
    '''
    fr = gp.Radiation()
    fr.SetAmbientDensity(ambient_density)
    fr.SetBField(b_field)
    fr.AddThermalTargetPhotons(t_cmb,edens_cmb)
    fr.SetDistance(distance)

    '''
        define the particle type. this will determine which radiation processes
        will be calculated. you can even set both!
    '''
    fr.SetProtons(power_law_spectrum)
    fr.SetElectrons(logparabola_spectrum)

    '''
        calculate the flux at an arbitrary range of gamma-ray energies (in erg)
    '''
    # define energies at which gamma-ray emission should be calculated 
    e = np.logspace(-6,15,200) * gp.eV_to_erg

    # do the calculation
    fr.CalculateDifferentialPhotonSpectrum(e)


    '''
        extract the different SEDs
    '''
    total_sed = fr.GetTotalSED() # SED, E^2dNdE (erg/s/cm^2) vs E (TeV)
    pp_sed = fr.GetPPSED() 
    synch_sed = fr.GetSynchrotronSED() 
    brems_sed = fr.GetBremsstrahlungSED()
    ic_sed = fr.GetICSED()

    # in order to conveniently use them with numpy, transpose to numpy arrays
    total_sed = np.array(total_sed)
    pp_sed = np.array(pp_sed)
    synch_sed = np.array(synch_sed)
    brems_sed = np.array(brems_sed)
    ic_sed = np.array(ic_sed)


    '''
        ########## PLOTS ################################
        plot particle spectra, radiation spectra and SEDs
        make a plot showing the particle spectra
    '''
    # particle spectra plot
    f,ax = plt.subplots(figsize=(5,5))
    ax.set_prop_cycle('color',plt.get_cmap('plasma')(np.linspace(0., .8, 2))) #
    plt.loglog(power_law_spectrum[:,0],power_law_spectrum[:,1],label="protons: power-law")
    plt.loglog(logparabola_spectrum[:,0],logparabola_spectrum[:,1],label="electrons: log-parabola")
    plt.xlabel("Energy (erg)")
    plt.ylabel("dN/dE (1/erg)")
    plt.ylim(ymin=1e31,ymax=1e55)
    plt.legend()
    plt.title("particle spectra")
    plt.grid()
    f.savefig("particle_spectra.png",bbox_inches='tight')

    # SED plot
    f,ax = plt.subplots(figsize=(5,5))
    ax.set_prop_cycle('color',plt.get_cmap('plasma')(np.linspace(0., .8, 5)))  #
    plt.loglog(total_sed[:,0],total_sed[:,1],lw=4,alpha=0.2,label="sum")
    plt.loglog(pp_sed[:,0],pp_sed[:,1],label="p-p")
    plt.loglog(synch_sed[:,0],synch_sed[:,1],label="synchrotron")
    plt.loglog(brems_sed[:,0],brems_sed[:,1],label="bremsstrahlung")
    plt.loglog(ic_sed[:,0],ic_sed[:,1],label="inverse-compton")
    plt.xlabel("Energy (TeV)")
    plt.ylabel("E^2dN/dE (erg/cm^2/s)")
    plt.ylim(ymin=1e-18,ymax=1e-8)
    plt.legend()
    plt.grid()
    plt.title("Radiation SEDs")
    f.savefig("radiation_SEDs.png",bbox_inches='tight')

