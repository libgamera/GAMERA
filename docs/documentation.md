[(back to main page)](main_page.md)

Units
=====
UNDER CONSTRUCTION
------------------

Units in GAMERA are mostly in cgs, with several exceptions that are meaningful
in the astrophysical context.

- source ages are in years
- distances and extensions are in parsecs
- when retrieving SEDs, for convenience the units are
  - erg / s / cm^2 vs TeV for radiation flux SEDs
  - erg / s vs TeV for luminosity SEDs (i.e. when the distance to the source is not set)
  - erg vs TeV for particle SEDs

  (Differential spectra are again in cgs units)
-


The units of the spectrum will impact the unit of the calculated radiation spectra.
For instance, if the unit of the particle spectrum only differential in energy, 
i.e. 1/erg, the output radiation spectra will have the unit of a flux if a 
source distance is specified (see below) or differential photon count per energy and time if not. 

On the other hand, if the input unit is differential also in volume, i.e. 1/erg/cm^3, 
then also the output radiation spectrum will be. Therefore, if no distance is 
specified in this case, a volume photon emissivity will be calculated, i.e. d³N/dEdVdt. 

Radiation models
================

In GAMERA, radiation mechanisms are implemented for
- Electrons (including Positrons)
- Protons (actually assuming a generic hadron mix, see ...)
- Arbitrary mix of hadronic species for both projectiles and target nuclei. 

Supported gamma-ray production mechanisms:

__for Electrons__

- Synchrotron Emission
  - isotropised pitch angle distribution of electrons [(following Ghisellini et al. 1988)](http://adsabs.harvard.edu/abs/1988ApJ...334L...5G) 
  - custom fixed pitch angles [(see e.g. Blumenthal & Gould 1970)](http://adsabs.harvard.edu/abs/1970RvMP...42..237B)

- Bremsstrahlung
  - Both electron-electron and electron-ion Bremsstrahlung [(following Baring et al. 1999)](http://adsabs.harvard.edu/abs/1999ApJ...513..311B)

- Inverse-Compton Emission 
  - using the full Klein-Nishina cross-section [(see e.g. Blumenthal & Gould 1970)](http://adsabs.harvard.edu/abs/1970RvMP...42..237B)
  - allowing for arbitrary target field spectra
  - supporting synchrotron-self-Compton(SSC) emission [(using Atoyan & Aharonian 1996)](http://adsabs.harvard.edu/abs/1996MNRAS.278..525A)
  - anisotropic radiation fields

__for Hadrons__

- π⁰- and η-Decay
  - Following the parameterisation of [Kafexhiu et al. 2014](http://adsabs.harvard.edu/abs/2014PhRvD..90l3014K)
 

Particle evolution
==================

[(back to main page)](main_page.md)


