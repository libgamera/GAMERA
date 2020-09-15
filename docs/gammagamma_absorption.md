[(back to main page)](main_page.md)

[(one page up)](tutorials_main.md)

Gamma-gamma absorption Tutorial
===============================

PRELIMINARY
-----------
This section is still under construction and the documentation will be improved with time.



In this tutorial we explain how to use the latest functionality added to GAMERA: the iclusion of gamma-gamma absorption effects
on the photon spectra.

All the new functionalities have been implemented in the `Radiation` class.

When setting up the target photon fields for the IC scattering, it is possible to compute also the effect of the gamma-gamma
absorption on the resulting radiation.

If the photon fiels is not isotropic but has an angular dependency, defined in the same way as done for the IC scattering
(see the tutorial on [Inverse Compton scattering](inverse_compton.md)), GAMERA can compute the absorbed flux using the full
angular dependent cross section, while it will use the integrated version otherwise.

To make the calculation possible, the photon field must be given a linear size to perform the spatial integration of
the absoprtion coefficient so as to obtain the actual value of the optical depth.

The size of the photon field is given with the following function:
```python
fr.SetSizePhotonField(0,sizeph)
```
where the first argument is the photon field counter and the second is the size of the photon field given in parsecs.

Once this is done, you can obtain the absorbed differential spectrum or the absorbed SED through the following functions:
```python
abs_sed = fr.GetTotalAbsorbedSED([0,1,..])
abs_spectrum = fr.GetTotalAbsorbedSpectrum([0,1,..])
```
where the argument is the array with the index of the photon fields we are taking into account for the absorption
This can be useful in case a user would like to consider a certain photon field for the scattering,
but not for the absoption and viceversa.

It is also possible to define a certain spatial dependency of the photon field through the function
```python
fr.SetTargetFieldSpatialDep(field, dependency)
```
where `field` is the number of the photon field and `dependency` is a vector of tuples 
in which the first element is the distance from the source (given in pc) and the second element is 
the fraction of the energy density that was given when setting up the photon field.

If you have more complicated cases, you can always do the spatial integration directlz in your code by calculating the absorption coefficient
for the field via `fr.ComputeAbsCoeff(energy, field)` which takes the energy of the gamma-ray and the field index as arguments, returning the 
absorption coefficient in units of `1/cm`.

<!-- **CAVEAT 1:** At the moment it is not implement a spatial dependency of the photon fields.
To take this into account, you could implement it in the scripting step, through the creation multiple fields with intensities
following a certain spatial function. -->

**CAVEAT:** No effect of the secondaries on the radiation spectrum is implemented.


