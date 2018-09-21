#!/usr/local/bin/python

import sys
sys.path.append('lib')
import gappa as gp
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    e_ref = gp.TeV_to_erg
    energy_in_erg_pl = np.logspace(-3,3,200) * gp.TeV_to_erg
    alpha_pl = 2; 
    power_law = energy_in_erg_pl**-alpha_pl
    fu = gp.Utils()
    power_law /= fu.Integrate(zip(energy_in_erg_pl,power_law * energy_in_erg_pl))
    b_field = 1e-5
    ambient_density = 1
    t_cmb = 2.7; edens_cmb = 0.25 * gp.eV_to_erg
    distance = 1e3 

    fr = gp.Radiation()
    fr.ToggleQuietMode()
    fr.SetAmbientDensity(ambient_density)
    fr.SetBField(b_field)
    fr.AddThermalTargetPhotons(t_cmb,edens_cmb)
    fr.SetDistance(distance)

    fr.SetProtons(list(zip(energy_in_erg_pl,1e49*power_law)))
    fr.SetElectrons(list(zip(energy_in_erg_pl,1e47*power_law)))

    fr.CalculateDifferentialPhotonSpectrum(np.logspace(-6,15,50) * gp.eV_to_erg)

    total_sed = np.array(fr.GetTotalSED())
    pp_sed = np.array(fr.GetPPSED())
    synch_sed = np.array(fr.GetSynchrotronSED())
    brems_sed = np.array(fr.GetBremsstrahlungSED())
    ic_sed = np.array(fr.GetICSED())


    print "total_sed:"
    print total_sed
    print "-------------"
    print "pp_sed:"
    print pp_sed
    print "-------------"
    print "synch_sed:"
    print synch_sed
    print "-------------"
    print "brems_sed:"
    print brems_sed
    print "-------------"
    print "ic_sed:"
    print ic_sed
    print "-------------"


    f,ax = plt.subplots(figsize=(5,5))
    ax.set_prop_cycle('color',plt.get_cmap('plasma')(np.linspace(0., .8, 2))) #
    plt.loglog(energy_in_erg_pl,1e49*power_law,label="protons")
    plt.loglog(energy_in_erg_pl,1e47*power_law,label="electrons")
    plt.xlabel("Energy (erg)")
    plt.ylabel("dN/dE (1/erg)")
    plt.ylim(ymin=1e31,ymax=1e55)
    plt.legend()
    plt.title("particle spectra")
    plt.grid()
    f.savefig("autotest/particle_spectra.png",bbox_inches='tight')

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
    f.savefig("autotest/radiation_SEDs.png",bbox_inches='tight')

