Step-by-step: Time independent modeling of Particle Spectra
===========================================================


This tutorial will cover the scenario where a injected population of electron evolves in a static 
system. By that it is meant that magnetic field and other environmental parameters 
as well as the injection spectrum don't change over time. 
Please note, however, that in the following for explanatory reasons also adiabatic expansion 
is assumed, albeit with constant velocity - so not quite static. This strictly leads to 
time-dependent losses, decreasing as $$\dot{E}_{ad}\sim t^{-1}$$.

 
The first thing one has to do is to create a `Particles`-object 
```
import gappa as gp
fp = gp.Particles()
```

--- Setting Up The Source Term / Injection Spectrum
An injection spectrum has to be specified. The most convenient and generic way is to define 
a 2d-numpy array and set it like for example here:

```
# powerlaw with spectral index = 2, norm = 1e37 erg/s, 
# e_ref = 1TeV, range between 10MeV and 1PeV
e = np.logspace(-5,3,bins) * gp.TeV_to_erg 
power_law = (e/gp.TeV_to_erg)**-2

# renormalise to 1e37 erg/s
power_law *= 1e37/fu.Integrate(zip(e,e*power_law))

# zip and set it up in the Particles-object
fp.SetCustomInjectionSpectrum(zip(e,power_law))
```

--- Setting Up Environmental Observables / The Energy Loss Rate
Environmental observables like B-Field and ambient density are then set via `Setter`-functions: 
```
fp.SetBField(b_field) # in Gauss
fp.SetAmbientDensity(density) # in particles/cm^3
fp.SetExpansionVelocity(vel) # in cm/s
fp.SetRadius(start_radius) # in pc
fp.SetAge(age) # in yrs
fp.AddThermalTargetPhotons(temperature,energy_density) #in K, erg/cm^3. 
```
(Note that there is a [https://www.mpi-hd.mpg.de/personalhomes/jhahn/Tutorials/Calculate%20Radiation%20From%20A%20Particle%20Population/Setting%20Up%20IC-Target%20Photons/index.html](separate tutorial) on how to set radiation fields for the IC process) 
 
The above commands determine the energy loss rate of the injected particles. You access 
this information with 

```
e_loss_rate  = fp.GetEnergyLossRate(e,age) # (erg,yrs)
cooling_time = fp.GetCoolingTimeScale(e,age)
```
which will return 2D-arrays containing `[energy(erg)~energy_loss_rate(erg/s)]` and 
`[energy(erg)~cooling time(yrs)]`. The cooling time is defined as `t_cool(E) = |E / (dE/dt)|`. 
The input energy bins should match the energy range of the input source spectrum, otherwise 
you might get a crash! Also, per default these functions will show losses that electrons, would 
suffer in the specified environment. You can toggle that by 
```
fp.SetType("electrons")
fp.SetType("protons")
```
In the case of protons, only possible adiabatic losses will be accounted for. 

The result for electrons looks like this: 

@LossRates.png[width:550px] 
 
Note that the energy loss rate at high energy transitions from $$\dot{E}\sim E$$ at 
early ages, where adiabatic expansion is the dominating cooling contribution, to 
the synchrotron and IC expectation $$\dot{E}\sim E^2$$. 

Unfortunately, it is currently not possible to show the cooling times and rates of 
the different cooling processes independently with one single `Particles`-object. However, 
you can define several objects where you set only B-field, expansion velocity etc. and 
extract the corresponding cooling plots - they will then correspond only to Synchrotron 
losses, adiabatic losses etc. 

Setting Up The Particle Escape Term
-----------------------------------
If you want to consider particle escape, you can use the following methods: 

**1: Apply an energy-independent escape time** 
```
fp.SetConstantEscapeTime(t_esc) # in seconds
```
**2: Apply an energy-dependent escape time** 
This is possible by specifying a 2D-vector holding the escape time as a function 
of energy and setting it by specifying
```
t_esc_0 = 1e3 * gp.yr_to_sec; t_dash = 1e3;
t_esc = zip(e,t_esc_0 * (t / t_dash)^0.5)
fp.SetEnergyDependentEscapeTime(t_esc)
```
There is a [https://www.mpi-hd.mpg.de/personalhomes/jhahn/Tutorials/Modeling%20Of%20Particle%20Spectra/Particle%20Escape/index.html](dedicated tutorial) on particle escape.


--- Compute The Time Evolution

Finally, the spectrum can be calculated via 
```
fp.CalculateElectronSpectrum()
```

This will calculate the spectrum at a time `age`, assuming the particles in questions are electrons, 
which means that Synchrotron, Bremsstrahlung, IC and adiabatic losses will be accounted for. 

If you want to calculate a hadron spectrum, you should instead specify 
```
fp.CalculateProtonSpectrum()
```
Please note that in case the losses are time-dependent a starting age for the particle spectrum 
evolution has to be specified. Per default, this is set to 1 year. However, in certain cases it 
might be necessary to change this number, which can be done by stating 
```
fp.SetTmin(tmin) # yrs
```


You can access the result via two options:
```
sp  = fp.GetParticleSpectrum() # returns diff. spectrum: E(erg) vs dN/dE (1/erg)
sed = fp.GetParticleSED()  # returns SED: E(TeV) vs E**2*dN/dE (erg)
```



[ParticlesBasic.py](Here) is a working script which will produce the following 
output (assuming electrons). Note that particle escape was not set in this example.

@ParticlesBasic.png[width:550px] 



