Effect of heavier nuclei for gamma-ray emission
===============================================


**Preliminary**


In this page we describe the new functionalities introduced for the calculation of the pi0 component
of a gamma ray spectrum coming in case we are in presence of different hadronic species or with a 
different ambient medium composition.

In the previus tutorials we have seen that to set up the contribution of pi0 emission, we need to
define a proton spectrum and set the ambient density of the material. This is done via 
(using the gappa python wrapper):

```python
import numpy as np
import gappa as gp
fr = gp.Radiation()

# Setting up the proton spectrum
e = np.logspace(-2, 5, 100) * gp.TeV_to_erg
dnde = e**(-2)
protons = np.array(list(zip(e, dnde)))

fr.SetProtons(protons)
fr.SetAmbientDensity(1.) # number density of the protons in cm^-3
```

In this case, we are defining a proton spectrum with energies between 0.01 and 1e5 TeV with a -2 spectral
index and an ambient density for protons of 1 cm-3.
The computations done by GAMERA at this point assume though an Interstellar Medium (ISM) that where there
is also 10% of Helium, so when we write
```python
fr.SetAmbientDensity(1.)
```
we are assiming 1 proton per cubic centimeter plus 0.1 Helium nuclei per cubic centimeter. These are also the values
used when computing `Bremsstrahlung` emission and `Ionization` losses.

New functionalities
-------------------

The new functionalities allow the user to specify additional hadronic species as relativistic particles
and arbitrary composition of the ISM.

**Add more hadronic species as relativistic particles**

This is done with the following function.

```python
fr.AddHadrons(np.array(list(zip(e, dnde))), 4) # to add a hadronic species with mass number 4
```
so the function takes 2 arguments, the differential spectrum of the particles and the mass number of them.
This function can be called multiple times and will keep add species to the flux.
When called with mass number one, it has exactly the same effect of `SetProtons`, so be careful not to double
accidentally the proton flux when using both functions.
It is possible to delete all previously added hadronic species with the function below:

```python
fr.ClearHadrons()
```

**Set arbitrary composition of the ISM**

This is done via

```python
composition = [[1, 1.], [4, 0.1], ...]
fr.SetAmbientMediumComposition(composition)
```
where `composition` is a list of tuples of the form (mass number, density in cm-3).

Be aware that calling `SetAmbientMediumComposition` will overwrite `SetAmbientDensity` and viceversa.
So these two functions should not be used together.

**Get back the previously defined hadronic spectra, mass numbers and the composition of the ISM**

The extract the previously defined hadronic spectra, the mass numbers and the composition of the ambient medium the following functions can be used:

```python
i = 0
hadron_spectrum_i = fr.GetHadrons(i)    # This returns the spectrum of the first hadron species defined with the AddHadrons-function.
                                        # i=1 will return the spectrum of the second species and so on.
hadron_masses = fr.GetHadronMasses()    # This returns a list with the masses of all previously defined hadron species
composition = fr.GetAmbientMediumComposition()    # This returns the list defined with the function SetAmbientMediumComposition(composition)
```

**Get the Gamma-ray spectra**

One can retrieve the Gamma-ray spectra and SEDs from each hadronic species defined with the function `AddHadrons`. But first, the emission has to be calculated, which is done with the same function as usually:
```python
e_values = np.logspace(-1,4,50)*gp.TeV_to_erg
fr.CalculateDifferentialPhotonSpectrum(e_values)
```
The total spectrum and the total SED, which can be retrieved as in other tutorials with
```python
total_spectrum = fr.GetTotalSpectrum()  # Returns the total spectrum
total_SED = fr.GetTotalSED()            # Returns the total SED
```
It already contains the contributions from all hadron species. To get the contribution from individual species, one can use the following function:

```python
i = 0
total_hadron_spectrum = fr.GetHadronSpectrum()  # To get the total spectrum from all hadron species
total_hadron_SED = fr.GetHadronSED()            # To get the total SED from all hadron species

Spectrum = fr.GetHadronSpectrum(i)    # To get the spectrum from hadron species number i
SED = fr.GetHadronSED(i)              # To get the SED from hadron species number i
```



**Caveat**

The densities that are being set up with `SetAmbientMediumComposition` will affect only the pi0 emission.
So far the `Bremsstrahlung` emission has not been yet modified accordingly. This will happen in the future.

Same is true if you want to compute the effects on the cooling times of the particles. We will be working on
it in the future.




