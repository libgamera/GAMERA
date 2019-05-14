import sys
sys.path.append("../GAMERA/lib/")
import gappa as gp
import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, acos


# Functions being used later to define the anisotropy field
def angular_difference(phi_0, th_0,phi,th):
    ang_dif = acos(sin(th)*cos(phi)*sin(th_0)*cos(phi_0) + sin(th)*sin(phi)*
            sin(th_0)*sin(phi_0) + cos(th)*cos(th_0))
    return ang_dif

def is_star(phi, th, phi_0, th_0, width):
    if (angular_difference(phi_0, th_0, phi,th) > width):
        return 0
    else: 
        return 1
    
    

fr = gp.Radiation() # Defining a new Radiation class object

energy = np.logspace(-3,4,100) * gp.TeV_to_erg
powerlaw = 1.e39*energy**-2
electrons = np.array(list(zip(energy, powerlaw)))
fr.SetElectrons(electrons)
fr.SetDistance(5)

temp  = 500  # Temperature in K
edens = 40 * gp.eV_to_erg    # Photon energy density in erg/cm^3

fr.AddThermalTargetPhotons(temp, edens)
fr.AddThermalTargetPhotons(temp, edens)

# Defining a circular top-hat anisotropy function located around ph_0 and th_0
ph_0 = 0.25*gp.pi
th_0 = 0.5*gp.pi
angular_size = 0.5    # Size of the tophead function in radian
phis = np.linspace(ph_0 - angular_size*2, ph_0 + angular_size*2, 20)
thetas = np.linspace(th_0 - angular_size*2, th_0 + angular_size*2, 20)

ph_m, th_m = np.meshgrid(phis, thetas)

is_star_vec = np.vectorize(is_star)
flat = is_star_vec(ph_m, th_m, ph_0, th_0, angular_size)
dphi =  phis[1] - phis[0]
dtheta = thetas[1] - thetas[0]
norm = np.sum(flat * dphi * dtheta* np.sin(th_m))   # Normalising the anisotropy distribution
flat = flat/norm

plt.figure(figsize=(7,7))
plt.pcolormesh(ph_m/gp.pi,th_m/gp.pi,flat)
plt.ylabel("theta [rad/pi]")
plt.xlabel("phi [rad/pi]")
plt.tight_layout()
plt.savefig("anisotropy_distribution.png")

# Set the anisotropy map for the target field 1
vis = np.array([0*gp.pi,0.5*gp.pi])
fr.SetTargetPhotonAnisotropy(1, vis, phis, thetas, flat)

## if you need isotropic electrons.
## WARNING at the moment this is very slow
# fr.SetElectronsIsotropic()

e = np.logspace(5,15, 100) * gp.eV_to_erg
# Unset the fast mode flag for computation of single components
fr.UnsetICFastMode()
fr.CalculateDifferentialPhotonSpectrum(e)

SEDtot = fr.GetTotalSED() 
SED1 =fr.GetICSED(0)      # Get the single components of the IC spectrum
SED2 =fr.GetICSED(1)

plt.figure(figsize=(7,7))
plt.xlabel("Energy [TeV]")
plt.ylabel("Flux [erg cm$^{-2}$ s$^{-1}$]")
plt.loglog(np.array(SEDtot)[:,0], np.array(SEDtot)[:,1],'C0',label="TOT case")
plt.loglog(np.array(SED1)[:,0], np.array(SED1)[:,1],'C1',label="Field 0 isotropic")
plt.loglog(np.array(SED2)[:,0], np.array(SED2)[:,1],'C2',label="Field 0 anisotropic")
plt.xlim(0.9e-5,1.1e3)
plt.ylim(4.9e-15,3.1e-12)
plt.legend()
plt.tight_layout()
plt.grid()
plt.savefig("SED_iso_aniso.png")

