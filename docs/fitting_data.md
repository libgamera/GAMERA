[(one page up)](tutorials_main.md)

Fitting data with GAMERA models
===============================

Here we document a quick example that shows how the pytohn framework can be used to fit data.

The example script is [here](gappa_fitter_example.py) and can be used as guide for further user's applications.

In the script we have provided 2 possibilities to obtain the final results:
 - fitting with scipy.optimize.curve_fit (see doc [here](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html))
 - fitting with the MCMC implementation given by the package [emcee](https://emcee.readthedocs.io/en/stable/)

In order to fully run the script, you need the following additional packages:
 - scipy (scripts works for version 1.6.0)
 - emcee (script works for version 3.0.2)
 - corner
 - multiprocessing

Description of the script
-------------------------

If you open the script, you can see that everything boils down to writing a nice model function that will accept
external parameters and will return as a result an array of spectral values:

```python
def model_for_fitting(ene, a, b, c, d, e):
    """
    Build model for data fitting.
    The parameters need to be separated

    The function is based on the previous one and can be modified by the user depending on the parameters that are
    needed to be fitted. The first argument is the photon energy of the points to be fitted (given in TeV)
    
    Arguments:
        - ene : array of energies of photons (in TeV)
        - a, b, c, d, e : parameters to be fitted
    """

    # print("Function call")  # DEBUG only
    fu = gp.Utils()
    fp = gp.Particles()
    fr = gp.Radiation()
    fp.ToggleQuietMode()
    fr.ToggleQuietMode()
    fp.SetSolverMethod(1)  # Use the semi-analytical model (which is faster) given that we don't have time evolution
                           # in the environment

    e_total_pl = 10 ** a  # total energy of the particles (logspace) in erg
    alpha_pl = b  # index of the particle PL
    bfield = 10 ** c  # B field of the environment (logspace) in Gauss
    temp = 10 ** d  # temperature of the photon field in K (linspace)
    edens = 10 ** e  # energy density of the photon field in erg/cm3 (linspace)
    age = 5  # age of the system in years (logspace) - FIXED
    distance = 1e3  # in pc - FIXED

    energy_in_erg_pl = np.logspace(-3, 3, 75) * gp.TeV_to_erg  # emin and emax of the PL are fixed
    power_law = energy_in_erg_pl ** -alpha_pl
    power_law *= e_total_pl / fu.Integrate(list(zip(energy_in_erg_pl, power_law * energy_in_erg_pl)))
    power_law_spectrum = np.array(list(zip(energy_in_erg_pl, power_law)))

    fp.SetBField(bfield)
    fp.AddThermalTargetPhotons(temp, edens, 150)
    fp.SetCustomInjectionSpectrum(
        power_law_spectrum)  # the injection rate is constant and given by the normalization of injected PL
    fr.AddArbitraryTargetPhotons(fp.GetTargetPhotons())
    fr.SetBField(bfield)
    fr.SetDistance(distance)

    time_steps = np.logspace(1, age, 2)
    # Need to have at least a timestep otherwise the evolution
    # is not computed
    for time in time_steps:
        fp.SetAge(time)
        fp.CalculateElectronSpectrum()
        sp = np.array(fp.GetParticleSpectrum())
        fr.SetElectrons(sp)
        # compute the spectrum on the points given by the data
        fr.CalculateDifferentialPhotonSpectrum(ene * gp.TeV_to_erg)
        rad = np.array(fr.GetTotalSED())
    return rad[:, 1]  # returns the y value of the function
```

Once this function is defined, the fits can be performed with the external packages:

**Using curve_fit**

This fit is a simple least square fit using a well tested routine available in scipy

```python
 print("Start the curve_fit part")
 # For testing, start from a point a bit further from true position
 p1 = [35.1, 2.3, -5, 1.1, np.log10(1 * gp.eV_to_erg)]
 res = opt.curve_fit(model_for_fitting, x, y, p0=p1, sigma=yerr, absolute_sigma=True)
 print(res)  # print the results
```

**Using emcee**

This method is more involved and will run a markov-chain Montecarlo (MCMC) on the dataset
in order not only to fit the data, but also to derive the posterior probability distribution of the
fitted paramerers in a Bayesian approach (see the package documentation for additional details).

Without going into many details, we can just say that in this approach we need some probability 
priors on the parameters (a good way to do it is to have uninformative flat priors in the interval
we think the parameter lives) and a likelihood function to be maximised by the fit.

To notice that this method can be quite slow, depending on how complex is your model, because it requires running
the model multiples times in order to get a good estimate of the final probability distribution.

Depending on the model, the prior and the likelihood can be written like:
```python
## Auxiliary functions for the MCMC
def log_prior(theta):
    """
    Uninformative flat prior.
    Needs to be adjusted depending on the chosen parameters
    
    Arguments:
        - theta : list of parameters
    """

    a, b, c, d, e = theta  # extract the parameters
    if 30 < a < 40 and 1.5 < b < 3 and -8 < c < -3 and 0 < d < 1.5 and -14 < e < -10:
        return 0.0
    return -np.inf


def log_likelihood(theta, x, y, yerr):
    """
    Compute the chi2.
    Needs to be adjusted depending on the chosen parameters
    
    Arguments:
        - theta : list of parameters
        - x, y, yerr : data points with errorbar
    """

    a, b, c, d, e = theta  # adjust here in case you have different parameters
    model = model_for_fitting(x, a, b, c, d, e)
    sigma2 = yerr ** 2 + model ** 2
    return -0.5 * np.sum((y - model) ** 2 / sigma2)


def log_prob(theta, x, y, yerr):
    """
    Compute the total log probability.
    
    Arguments:
        - theta : list of parameters
        - x, y, yerr : data points with errorbar
    """

    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y, yerr)
```

at this point we can run the MCMC chain:

```python
print("Start the MCMC part")
import os
os.environ["OMP_NUM_THREADS"] = "1"  # to avoid problems with numpy and Pool. Follows emcee guidelines.

pos = np.array(res[0]) + 1e-3 * np.random.randn(32, 5)  # initial position of the walkers
nwalkers, ndim = pos.shape

####
## Note:
## In this example, due to the long time to run the model, the burn-in phase is quite short and
## would be better to run everything for a larger number of steps
## (like 100 burn-in and 1000 of the chain could be a good option)
burn_in_steps = 10
chain_steps = 100

with Pool() as pool:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, pool=pool, args=(x, y, yerr))
    state = sampler.run_mcmc(pos, burn_in_steps,
                             progress=True)  # saves the position of the walkers in the state variable
    print("End burn-in phase")
    sampler.reset()  # resets the chain
    sampler.run_mcmc(state, chain_steps, progress=True)  # start again the chain form after the burn-in
    print("End MCMC")
```
because this MCMC implementation uses multiple *walkers* to test the method and compute likelihoods,
it can be parallelised. Hence the call to the Pool() method. 

If you want to run on a single CPU, remove the with statement and the `os.environ` line.
