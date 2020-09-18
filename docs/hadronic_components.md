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

**Set arbitrary composition of the ISM**

This is done via

```python
composition = [[1, 1.], [4, 0.1], ...]
fr.SetAmbientMediumComposition(composition)
```
where `composition` is a list of tuples of the form (mass number, density in cm-3).

Be aware that calling `SetAmbientMediumComposition` will overwrite `SetAmbientDensity` and viceversa.
So these two functions should not be used together.

**Caveat**

The densities that are being set up with `SetAmbientMediumComposition` will affect only the pi0 emission.
So far the `Bremsstrahlung` emission has not been yet modified accordingly. This will happen in the future.

Same is true if you want to compute the effects on the cooling times of the particles. We will be working on
it in the future.




