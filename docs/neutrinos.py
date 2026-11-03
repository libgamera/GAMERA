##########################################################################
# Tutorial explaining the calculation of neutrinos with GAMERA.
# The implementation is using the parametrisation from Kelner et al. 2006
# https://ui.adsabs.harvard.edu/abs/2006PhRvD..74c4018K/abstract
# author: Mischa Breuhaus
##########################################################################

import gappa as gp
import numpy as np
import matplotlib.pyplot as plt


fr = gp.Radiation()

#--------------------- Add protons ----------------------------------
evalues = np.logspace(-4,15, 500)*gp.TeV_to_erg
e_cut = 60.*gp.TeV_to_erg
dn_dE_proton = evalues**(-2)*np.exp(-evalues/e_cut)

fu = gp.Utils()
norm = fu.EnergyContent(list(zip(evalues, dn_dE_proton)))
dn_dE_proton = dn_dE_proton / norm * 1e50       # Normalised to e50 erg

protons = list(zip(evalues, dn_dE_proton))

fr.SetProtons(protons)

#-------------------- Add other particles ---------------------------
e_cut_HE = e_cut * 2
dn_dE_HE = evalues**(-2)*np.exp(-evalues/e_cut_HE)
norm_HE = fu.EnergyContent(list(zip(evalues, dn_dE_HE)))
dn_dE_HE = dn_dE_HE / norm_HE * 0.5e50       # Normalised to 0.5e50 erg


fr.AddHadrons(list(zip(evalues, dn_dE_proton)),  1.0)   # A = 1 (protons)
fr.AddHadrons(list(zip(evalues, dn_dE_HE)), 4.0)   # A = 4 (helium)




#------------- Set the ambient medium density--------------------------------
fr.SetAmbientDensity(1.)
fr.SetDistance(1.0e3)

#---------------- Calculate neutrinos ---------------------------------------
E = np.logspace(-1, 2, 100) * gp.TeV_to_erg
fr.CalculateNeutrinoSpectrum(E)


# Total neutrino SED from each source
sed_protons = np.array(fr.GetProtonNeutrinoSED())       # from SetProtons()
sed_H       = np.array(fr.GetHadronNeutrinoSED(0))      # AddHadrons A=1
sed_He      = np.array(fr.GetHadronNeutrinoSED(1))      # AddHadrons A=4
sed_total   = np.array(fr.GetNeutrinoSED())             # sum of all




# ---------- Plot ----------
fig, ax = plt.subplots(figsize=(8, 5))
 
ax.loglog(sed_total[:, 0],   sed_total[:, 1],
          color="black",      lw=2,   label="Total")
ax.loglog(sed_protons[:, 0], sed_protons[:, 1],
          color="tab:blue",   lw=1.5, ls="--", label="Protons (SetProtons)")
ax.loglog(sed_H[:, 0],       sed_H[:, 1],
          color="tab:red",  lw=1.5, ls=":", label="Protons (AddHadrons, A=1)")
ax.loglog(sed_He[:, 0],      sed_He[:, 1],
          color="tab:green",    lw=1.5, ls="-.",  label="Helium (AddHadrons, A=4)")
ax.set_ylim(1.0e-13, 1.0e-10) 
ax.set_xlabel(r"$E\;[\mathrm{TeV}]$", fontsize=13)
ax.set_ylabel(r"$E^2\,\mathrm{d}N/\mathrm{d}E\;[\mathrm{erg\,cm^{-2}\,s^{-1}}]$", fontsize=13)
ax.set_title("Total neutrino SED from pp interactions", fontsize=13)
ax.legend(fontsize=11)
ax.grid(True, which="both", alpha=0.3)
plt.tight_layout()
plt.savefig("neutrino_sed_total.png", dpi=150)
plt.show()