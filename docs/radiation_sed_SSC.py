# get_ipython().magic(u'pylab')
import sys
sys.path.append('<your gamera path>/lib/')
import gappa as gp
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    green = '#057775'
    fu = gp.Utils()
    
    # Set up the parameters for the electron spectrum
    e_ref = gp.TeV_to_erg # reference 1 TeV
    el_energy_in_erg_pl = np.logspace(-5,2,200) * gp.TeV_to_erg
    alpha_pl_el = 2.3
    beta_el = 0.3
    e_total_pl = 1e49
    
    # Set the electron spectrum
    elLogP = el_energy_in_erg_pl**-(alpha_pl_el+beta_el*np.log10((el_energy_in_erg_pl/e_ref)))
    elLogP *= e_total_pl/fu.Integrate(list(zip(el_energy_in_erg_pl,elLogP*el_energy_in_erg_pl)))
    elLogPsp = np.array(list(zip(el_energy_in_erg_pl,elLogP)))
    
    # Set other ambient parameters
    b_field = 1e-6  #Gauss
    t_1 = 30000   # Kelvin
    edens_1 = 2. * gp.eV_to_erg  # final value in erg/cm3
    distance = 1e4  # in pc
    
    # Initialize the radiation class
    fr = gp.Radiation() 
    fr.SetBField(b_field)
    fr.AddThermalTargetPhotons(t_1,edens_1)  # add a thermal field
    fr.SetDistance(distance)
    
    fr.SetElectrons(elLogPsp)  # set the electrons
    
    e = np.logspace(-6,15,200) * gp.eV_to_erg  # energy range of the gamma-ray emission
    fr.CalculateDifferentialPhotonSpectrum(e)  # do a first calculation of the radiation
    
    # Add as target photon the Synchrotron radiation previously computed
    # assuming a size of the region of 1 pc and then recompute the
    # radiation output
    fr.AddSSCTargetPhotons(1.)
    fr.CalculateDifferentialPhotonSpectrum(e)
    
    total_sed = np.array(fr.GetTotalSED()) # SED, E^2dNdE (erg/s/cm^2) vs E (TeV)
    synch_sed = np.array(fr.GetSynchrotronSED())
    ic_sed = np.array(fr.GetICSED(0))
    ssc_sed = np.array(fr.GetICSED(1))
    
    # plot the particle spectrum
    plt.figure()
    ax1=plt.gca()
    ax1.loglog(elLogPsp[:,0],elLogPsp[:,1],c='C2',label="electrons",linestyle="--")
    ax1.set_xlabel("Energy (erg)")
    ax1.set_ylabel("dN/dE (1/erg)")
    ax1.set_ylim(ymin=1e38,ymax=1e55)
    ax1.legend()
    ax1.grid()
    ax1.set_title("Particle Spectra")
    
    # SED plot
    plt.figure()
    ax3=plt.gca()
    ax3.loglog(total_sed[:,0],total_sed[:,1],alpha=0.2,c='k',label="sum",lw=3)
    ax3.loglog(synch_sed[:,0],synch_sed[:,1],c='C2',label="Synchrotron",linestyle="--")
    ax3.loglog(ic_sed[:,0],ic_sed[:,1],c='C2',label="Inverse Compton (thermal)",linestyle=":")
    ax3.loglog(ssc_sed[:,0],ssc_sed[:,1],c='C2',label="SSC",linestyle="-.")
    
    ax3.set_xlabel("Energy (TeV)")
    ax3.set_ylabel("E"+r"$^2$"+"dN/dE (erg/cm"+r"$^2$"+"/s)")
    ax3.set_ylim(ymin=1e-18,ymax=1e-10)
    ax3.legend()
    ax3.grid()
    ax3.set_title("Radiation SEDs")
    
    plt.show()
