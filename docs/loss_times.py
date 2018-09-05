#!/usr/local/bin/python
import sys
sys.path.append(YOUR_GAMERA_LIB)
import gappa as gp
import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":
    
    # set up particle spectrum with random environmental parameters
    fp = gp.Particles()
    fp.AddThermalTargetPhotons(10,10*gp.eV_to_erg)
    fp.SetAmbientDensity(1)
    fp.SetRadius(1)
    fp.SetExpansionVelocity(1e9)
    fp.SetBField(5e-6)
    fp.SetAge(1e5)

    # extract the cooling time scales at energy points 'e'
    e = np.logspace(-5,5,150)
    total = np.array(fp.GetCoolingTimeScale(e,"sum"))
    ic = np.array(fp.GetCoolingTimeScale(e,"inverse_compton"))
    brems = np.array(fp.GetCoolingTimeScale(e,"bremsstrahlung"))
    ad = np.array(fp.GetCoolingTimeScale(e,"adiabatic_losses"))
    synch = np.array(fp.GetCoolingTimeScale(e,"synchrotron"))

    # make a plot
    f = plt.figure(figsize=(5,5))
    plt.loglog(brems[:,0],brems[:,1],c="red",label="brems.")
    plt.loglog(ic[:,0],ic[:,1],c="orange",label="IC")
    plt.loglog(ad[:,0],ad[:,1],c="green",label="adiab.")
    plt.loglog(synch[:,0],synch[:,1],c="blue",label="synch")
    plt.loglog(total[:,0],total[:,1],c="black",lw=1,ls="--",label="sum")
    plt.xlabel("energy (erg)")
    plt.ylabel("cooling time scale (yrs)")
    plt.grid()
    plt.legend()
    f.savefig("loss_times.png",bbox_inches="tight")
