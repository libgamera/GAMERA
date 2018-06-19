#!/usr/local/bin/python

import sys
sys.path.append('PATH/TO/YOUR/GAMERA_LIB/DIRECTORY')
import gappa as gp
import numpy as np
import matplotlib.pyplot as plt


''' 
    this little script will showcase how to use the different hadronic interaction
    models for the calculation of gamma-ray spectra. The method is taken from
    Kafexhiu et al., Physical Review D, Volume 90, Issue 12, id.123014, 2014 and includes
    prediction models GEANT 4.10.0, PYTHIA 8.18, SIBYLL 2.1 and QGSJET-I.
'''
if __name__ == "__main__":

    green = '#057775'
    fu = gp.Utils()
    fu.DrawGamera()

    '''
        define an arbitrary particle spectrum. It needsto be 
        2-dimensional arrays (or vectors if using C++) with each
        entry holding a tuple: E(erg) - dN/dE (1/erg). Let's define a 
        simple power-law
    '''

    # define reference energy as 1 TeV
    e_ref = gp.TeV_to_erg

    '''
        example: define a power-law spectrum
    '''
    # define the energy range of the particles (total energy), here between 1 GeV and 1 PeV
    energy_in_erg_pl = np.logspace(-3,3,100) * gp.TeV_to_erg

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
        now set the parameters needed for the calculation 
        of the photon emission
    '''
    ambient_density = 1 # 1/cm^3
    distance = 1e3 # optional, in parsec. If not set or equals zero, differential photon production rate
                   # instead of flux will be calculated

    '''
        create a Radiation object and set it up
    '''
    fr = gp.Radiation()
    fr.SetAmbientDensity(ambient_density)
    fr.SetDistance(distance)

    fr.SetProtons(power_law_spectrum)

    

    '''
        Set the hadronic emission model and
        calculate the flux at an arbitrary range of gamma-ray energies (in erg)
    '''

    # define energies at which gamma-ray emission should be calculated 
    e = np.logspace(-5,3,100) * gp.TeV_to_erg


    # You can also pick the desired emission models via the 'PiModel' variable
    # in the Radiation class. PiModel is an integer, the values mean:
    # 0 - GEANT 4.10.0
    # 1 - PYTHIA 8.18
    # 2 - SIBYLL 2.1
    # 3 - QGSJET-I
    # The default and pre-set value is PiModel = 0 (GEANT 4)

    f = plt.figure(figsize=(6,6))
    models = [0,1,2,3] 
    model_names = ["GEANT 4.10.0","PYTHIA 8.18","SIBYLL 2.1","3 - QGSJET-I"]
    line_styles = ["-","--",":","-."]
    for m,n,ls in zip(models,model_names,line_styles):
        print("Using "+n+"\n")
        fr.SetPPEmissionModel(m) # SET THE EMISSION MODEL
        fr.CalculateDifferentialPhotonSpectrum(e)
        pp_sed = fr.GetPPSED() # extract the radiation SED
        pp_sed = np.array(pp_sed) # convert to numpy array for convenience
        plt.loglog(pp_sed[:,0],pp_sed[:,1],c=green,ls=ls,label=n)
   


    plt.xlabel("Energy (TeV)")
    plt.ylabel("E^2dN/dE (erg/cm^2/s)")
    plt.ylim(ymin=1e-13,ymax=2e-11)
    plt.legend()
    plt.grid()
    plt.title("Different Hadronic Emission Models")
    f.savefig("HadronicEmissionModelsSEDs.png",bbox_inches='tight')

