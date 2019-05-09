[(back to main page)](main_page.md)

[(one page up)](tutorials_main.md)

Inverse Compton Tutorial
========================

In this tutorial we will see in some more detail how to set up the radiation fields for the Inverse Compton scattering
for both isotropic and anisotropic case. In the following we assume that `fr` is an instance of the class `Radiation`.

GAMERA allows the user to add directly a thermal photon field knowing the temperature and the energy density through the
function:
```
fr.AddThermalTargetPhotons(temperature,energy_density)
```

For a more particular usage, it is also possible to add arbitrary photon fields via
```
fr.AddArbitraryTargetPhotons(photon_array)
```
where `photon_array` is an array of tuples \(E,photon_density\) in units of `erg` vs `erg^-1 cm^-3`.

A further possibility is to import a photon field from file with:
```
fr.ImportTargetPhotonsFromFile(filename)
```
where the structure of the file must be `energy [eV]` vs `number density [cm^-3]`.


Case 1: Isotropic photon fields
-------------------------------

Let's start with the simple case of a thermal photon field. In principle we can add multiple photon fields. In the default
case, when asking for the Inverse Compton component of the gamma ray sectrum, it will be computed on the sum of the
photon fields.

<!-- Add example! -->

If we explicitily want the component associated to the i-th field, the code will recompute the component associated only
to that field

<b>to be continued</b>


[(one page up)](tutorials_main.md)

[(back to main page)](main_page.md)
