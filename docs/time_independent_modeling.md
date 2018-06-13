Step-by-step: Time independent modeling of Particle Spectra
===========================================================


This tutorial will cover the scenario where a injected population of electron evolves in a static 
system. By that it is meant that magnetic field and other environmental parameters 
as well as the injection spectrum don't change over time. 

Step 1: create a `Particles`-object
----------------------------------

```
import gappa as gp
fp = gp.Particles()
```

Step 2: Set up the injection spectrum
-------------------------------------

Here we will set a power law injection spectrum
```
# powerlaw with spectral index = 2, norm = 1e37 erg/s, 
# e_ref = 1TeV, range between 10MeV and 1PeV
e = np.logspace(-5,3,bins) * gp.TeV_to_erg 
power_law = (e/gp.TeV_to_erg)**-2

# renormalise to 1e37 erg/s
fu = gp.Utils()
power_law *= 1e37/fu.Integrate(zip(e,e*power_law))

# zip and set it up in the Particles-object
fp.SetCustomInjectionSpectrum(zip(e,power_law))
```

Step 3: Set up environmental parameters
---------------------------------------

environmental parameters like B-Field and ambient density are set via `Setter`-functions: 
```
fp.SetBField(b_field) # in Gauss
fp.SetAmbientDensity(density) # in particles/cm^3
fp.AddThermalTargetPhotons(temperature,energy_density) #in K, erg/cm^3. 
```
>Note: 
>There is a [dedicated tutorial](inverse_compton.md) on how to set radiation fields for the IC process 
 
The above commands determine the energy loss rate of the injected particles. You access this information with 

```
e_loss_rate  = fp.GetEnergyLossRate(e,age) # (erg,yrs)
cooling_time = fp.GetCoolingTimeScale(e,age)
```
which will return 2D-arrays containing `[energy(erg)~energy_loss_rate(erg/s)]` and 
`[energy(erg)~cooling time(yrs)]`. The cooling time is defined as `t_cool(E) = |E / (dE/dt)|`. 

@LossRates.png[width:550px] 
 
Note that the energy loss rate at high energy transitions from $$\dot{E}\sim E$$ at 
early ages, where adiabatic expansion is the dominating cooling contribution, to 
the synchrotron and IC expectation $$\dot{E}\sim E^2$$. 


Step 4 (Optional): Set up the particle escape term
--------------------------------------------------

there are several options for this. Either a constant, an energy-dependent, a 
time-dependent or a both energy- and time-dependent escape time can be applied.
[Here](particle_escape.md) you can learn how to do this. In the following we
just assume a constant escape time value:

```
fp.SetConstantEscapeTime(t_esc) # in seconds
```

Step 5: Set starting point of iteration
---------------------------------------

The program needs to know at which point in time to start the iteration.
It will then set the initial condition of the system accordingly. 

```
fp.SetTmin(tmin) # yrs
```
The default is set to one year. 
The initial condition is assumed to be the result of loss-free injection of
particles until t = tmin. The spectral shape of the injection up until this point
is fixed to the value at t = tmin, i.e. Q(t \le tmin) =  Q(t=tmin). 
Only at t \ge tmin will cooling, escape and the actual evolution of Q be taken
into account.


Step 5: Set the source age
--------------------------

the specified system will be evolved to the point in time t = age. You can set
the age via

```
fp.SetAge(age) # in yrs
```




Step 5: Compute The Time Evolution
----------------------------------

Finally, the spectrum can be calculated:

__For electrons:__
```
fp.CalculateElectronSpectrum()
```

__For protons:__
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



