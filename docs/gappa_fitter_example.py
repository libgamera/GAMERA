"""
Script to illustrate how the python wrapper
of GAMERA can be used to fit data.

BETA version

The script is structured as follows:
 - produce a synthetic dataset
 - Fit the data with the routines coded in scipy.optimize.curve_fit
   (tested with scipy version: 1.6.0. 
   See: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html)
 - Fit the data using a MCMC implementation from the package emcee
   (tested with version 3.0.2. See: https://emcee.readthedocs.io/en/stable/)

The parameters that are left to vary in the fit are:
   - total energy in electrons
   - spectral index of the electrons
   - magnetic field
   - temperature of the external photon field
   - energy density of the external photon field
Other fixed parameters:
   - age of the system (10^5 years)
   - distance 10^3 pc

Author: Carlo Romoli - MPIK
"""

import sys

sys.path.append('<your-GAMERA-path>/lib/')
import gappa as gp
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import emcee
from multiprocessing import Pool
import corner

# np.random.seed(1)  # Fixed seed for testing

def sample_final_sed(emin, emax, npoints):
    """
    Function to create a synthetic dataset starting from the model
    Useful for testing

    Arguments:
        - emin : log10 of the minimum energy needed in the photon spectrum (TeV)
        - emax : log10 of the maximum energy needed in the photon spectrum (TeV)
        - npoints : number of synthetic points randomly chose in the range
    """

    logene = emin * np.random.rand(npoints) + emax  # random number for the energy of the photons (in log TeV)
    logene.sort()
    energy = 10 ** logene
    return energy


def model(par, synthetic = False):
    """
    Produce the full model.

    Arguments:
        - par : array of the model parameters.
                Except the spectral index, all the parameters are given in log10 of the quantity.
                Parameters are:
                - total electron energy
                - spectral index of the electron PL
                - magnetic field
                - temperature of the target photon field
                - age of the system
                - distance
        - synthetic : bool argument to have a full evenly sampled spectrum or a randomly sampled distribution of points
    """

    fu = gp.Utils()
    fp = gp.Particles()
    fr = gp.Radiation()
    fp.ToggleQuietMode()  # do no print the output
    fr.ToggleQuietMode()

    e_total_pl = 10 ** par[0]  # erg
    alpha_pl = par[1]
    bfield = 10 ** par[2]  # gauss
    temp = 10 ** par[3]  # Kelvin
    edens = 10 ** par[4]  # erg/cm3
    age = par[5]  # in log10(1/yr)
    distance = 10 ** par[6]  # pc

    energy_in_erg_pl = np.logspace(-3, 3, 200) * gp.TeV_to_erg  # emin and emax of the PL are fixed
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

    time_steps = np.logspace(1, age, 10)
    rad = []  # auxiliary lists to store the distributions for every step
    part = []
    for time in time_steps:
        fp.SetAge(time)
        fp.CalculateElectronSpectrum()
        sp = np.array(fp.GetParticleSpectrum())
        part.append(np.array(fp.GetParticleSED()))
        fr.SetElectrons(sp)
        if synthetic:
            fr.CalculateDifferentialPhotonSpectrum(sample_final_sed(-19, 3, 35) * gp.TeV_to_erg)
        else:
            fr.CalculateDifferentialPhotonSpectrum(np.logspace(-19, 3, 100) * gp.TeV_to_erg)
        rad.append(np.array(fr.GetTotalSED()))
    return part, rad  # returns the SEDs (TeV vs erg for particles and TeV vs erg/cm2/s for photons)


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


if __name__ == "__main__":

    ### Produce the full model on evenly sampled range to show original
    p0 = [35, 2.2, -5, 1, np.log10(1 * gp.eV_to_erg), 5, 3]
    part, rad = model(p0, synthetic=False)

    ## plot the model with the evolution of the particles. To show full result
    time_steps = np.logspace(1, 5, 10)
    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    ax1.set_prop_cycle('color', plt.get_cmap('plasma_r')(np.linspace(0.2, .8, len(time_steps))))  #
    for p, t in list(zip(part, time_steps)):
        ax1.loglog(p[:, 0], p[:, 1], label=str(int(t)))
    ax1.set_xlabel("E (TeV)")
    ax1.set_ylabel("E" + r"$^2$" + "dN/dE (erg)")
    ax1.grid()
    ax1.legend(ncol=2, prop={'size': 10}, title="age(yrs):")
    ax1.set_xlim(xmin=1e-6, xmax=2e3)
    ax1.set_ylim(ymin=1e40, ymax=1e48)
    ax1.set_title("Particle SEDs")

    ax2.set_prop_cycle('color', plt.get_cmap('plasma_r')(np.linspace(0.2, .8, len(time_steps))))  #
    for r, t in list(zip(rad, time_steps)):
        ax2.loglog(r[:, 0], r[:, 1], label=str(int(t)))
    ax2.set_xlabel("E (TeV)")
    ax2.set_ylabel("E" + r"$^2$" + "dN/dE (erg/cm^2/s)")
    ax2.legend(ncol=2, prop={'size': 10}, title="age(yrs):")
    ax2.grid()
    ax2.set_xlim(xmin=1e-19, xmax=1e4)
    ax2.set_ylim(ymin=1e-16, ymax=1e-10)
    ax2.set_title("Radiation SEDs")

    ## load the synthetic dataset and plot it
    part, sed = model(p0, synthetic=True)
    x = sed[-1][:, 0]  # energy in TeV
    y = sed[-1][:, 1]  # y value un erg/cm2/s
    yerr = sed[-1][:, 1] * 0.1  # assume an arbitrary 10% error in the y axis

    f2 = plt.figure()
    plt.loglog(x, y, 'ro-', label='synthetic dataset')

    #####
    # Fit the synthetic data using
    # scipy.optimize.curve_fit
    print("Start the curve_fit part")
    # For testing, start from a point a bit further from true position
    p1 = [35.1, 2.3, -5, 1.1, np.log10(1 * gp.eV_to_erg)]
    res = opt.curve_fit(model_for_fitting, x, y, p0=p1, sigma=yerr, absolute_sigma=True)
    print(res)  # print the results

    # Plot the best fit results for some reasons the number of points is not exactly the same, so there is the need
    # for small adjustment
    plt.loglog(np.logspace(-19, 3, 100)[:-1],
               model_for_fitting(np.logspace(-19, 3, 100) * gp.TeV_to_erg, res[0][0], res[0][1], res[0][2], res[0][3],
                                 res[0][4]), 'b-', label='best fit')
    plt.xlabel("Energy [TeV]")
    plt.ylabel("$E^2dN/dE$ [erg/cm2/s]")

    ######
    # Test using the emcee package and the MCMC chain
    print("Start the MCMC part")
    import os
    os.environ["OMP_NUM_THREADS"] = "1"  # to avoid problems with numpy and Pool.

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

    flat_samples = sampler.get_chain(
        flat=True)  # the burn in could also be set here, with the discard argument. thin option?
    print(flat_samples.shape)

    labels = ["log10(tot energy)", "sp. index", "log10(Bfield/1G)", "log10(Temp/1K)", "log10(e. density/(1erg/cm3))"]
    truth = [35, 2.2, -5, 1, np.log10(1 * gp.eV_to_erg)]

    ## This shows the correlation plot between the parameters
    ## The lines are the original true values that were used to
    ## obtain the models
    fig = corner.corner(flat_samples, labels=labels, truths=truth)

    plt.show()
