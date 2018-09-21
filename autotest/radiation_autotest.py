#!/usr/local/bin/python

import sys
sys.path.append('lib')
import gappa as gp
import numpy as np

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

