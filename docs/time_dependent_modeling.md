# Time-Dependent Modeling

Time-dependent modeling, i.e. modeling the spectral evolution of sources while 
physical parameters like magnetic fields are changing over time, can be done 
in `GAMERA` in much the same way as static modeling described [in a separate tutorial](time_independent_modeling.md).

## Time-Dependent Environmental Observables And Energy Losses
If you want to model time-dependent observables such as the ambient magnetic field 
and therefore also time-dependent energy losse, the only difference is instead 
of giving skalars to the `Setter`-functions you have to provide 2D-vectors. 
These vectors need to constist of two columns, with the first column holding the 
time steps and the second the value of the environmental parameter at that time. 
For instance, if you wanted to define a magnetic field declining with time as 
 
![bfield_tdep](bfield_tdep.png)
 
you would have to implement this like that:

```
fp = Particles()
t_ref = 100
time_steps = np.logspace(1,4,100) #in years
b = 1e-4 * (time_steps / t_ref) ** -0.5 #in Gauss
b_lookup = zip(time_steps,b)
fp.SetBField(b_lookup)
```

At the moment, time-dependent 
- ambient density 
- source luminosity 
- expansion velocity and radius
can be implemented in this fashion. 

CurrentlyTime-dependent radiation fields are only possible via an iterative
approach. There is an [extra tutorial on iterating](iteration.md).
 
## Time-Dependent But Energy-Independent Particle Escape Time-Scales
Very similar is the treatment for energy-independent but time-dependent 
particle escape time, here an example:
```
t_esc_0 = 1000 * gp.yr_to_sec
t_esc = t_esc_0 * (time_steps / t_ref) ** 0.5
t_esc_lookup = zip(time_steps)
fp.SetTimeDependentEscapeTime(t_esc_lookup)
```

## Time-Dependent Spectral Shapes Of Injection Spectra and Particle Escape Time-Scales
`GAMERA` can handle a time-dependent spectral shapes of the injection spectrum 
as well as the particle ecape. 
The former point is described [in this tutorial](tt.md). 
A dynamically changing spectral shape of the particle escape time-scale is covered 
in this [tutorial](particle_escape.md).
 

That's the whole difference to static modeling as far as the `GAMERA` interface 
is concerned. A working script can be found [ParticlesTimeDep.py](here), which will output the following plots: 

Temporal evolution of the particle spectrum assuming the following time dependences 
(escape time not set in this example): 
 
![time_evolution_basic](time_evolution_basic.png)

![particle_time_evolution_basic](particle_time_evolution_basic.png)
