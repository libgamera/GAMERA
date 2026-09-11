#!/usr/local/bin/python
import sys
import os
sys.path.append(os.path.abspath('../lib'))
import gappa as gp
import numpy as np
import copy
import matplotlib.pyplot as plt
import matplotlib as mpl

CMAP = mpl.cm.get_cmap('plasma_r')




'''
The function 'gaus' calculates gaussian distributed
values for the angular values ph_m and th_m with means
ph_0 and th_0 and width s
'''
def gaus(ph_0,th_0,s,ph_m,th_m):

    dphi = phi[1]-phi[0]
    dtheta = theta[1]-theta[0]

    g = np.exp(-((ph_m-ph_0)**2+(th_m-th_0)**2)/(2*gp.pi*s**2))
    g_n = np.exp(-((ph_m-ph_0)**2+(th_m-th_0)**2)/(2*gp.pi*s**2))*np.sin(th_m)

    norm = np.sum(g_n * dphi * dtheta)
    return g / norm




def flat_map(ph_0,th_0,s,ph_m,th_m):

    dphi = phi[1]-phi[0]
    dtheta = theta[1]-theta[0]
    flat = np.sqrt((ph_m-ph_0)**2+(th_m-th_0)**2)**0
    norm = np.sum(flat * np.sin(th_m) * dphi * dtheta)
    print(4*gp.pi/norm)
    return flat / norm 


if __name__ == "__main__":
    fu = gp.Utils(); fu.DrawGamera();
    mpik_green = '#057775'



    '''
        another possibility, especially useful for numpy users, is to 
        give a meshgrid plus axis vectors to SetTimeAndEnergyDependentEscapeTime
    '''

    t = 3e4
    edens = 1#gp.eV_to_erg
    


    phis = np.linspace(0,gp.pi,10)
    th_0 = 0.5*gp.pi
    width =0.004*gp.pi # rad
    ph_0 = 1*gp.pi

    vis = np.array([0,0.5*gp.pi+4./3.])#*gp.pi # Observation angle
    widths = np.linspace(0.01,1,10)*gp.pi
    bins = 20
    ebins = 10
    '''
        define injection spectrum
    '''
    e = np.logspace(-5,3,ebins) * gp.TeV_to_erg
    power_law = 1e40 * ((e/gp.TeV_to_erg)**-3)
    power_law = np.array(list(zip(e,power_law)))


    phi = np.linspace(ph_0-2*width,ph_0+2*width,bins)
    theta = np.linspace(th_0-2*width,th_0+2*width,bins)

    ph_m, th_m = np.meshgrid(phi, theta)
    flat = gaus(ph_0,th_0,width,ph_m,th_m)#**0 - 1
#    flat[len(flat)/2][len(flat)/2] = 1

    f = plt.figure(figsize=(5,5))
    plt.pcolormesh(ph_m,th_m,flat)
    plt.ylabel("theta [pi]")
    plt.xlabel("phi [pi]")
    f.savefig("flat_map.png")


    # Set distributions for the anisotropic case
    ra_aniso = gp.Radiation()
    ra_aniso.SetElectrons(power_law)
    ra_aniso.AddThermalTargetPhotons(t,edens)
    ra_aniso.SetTargetPhotonAnisotropy(0,vis,phi,theta,flat)

    # Set distributions for the isotropic case
    ra_iso = gp.Radiation()
    ra_iso.SetElectrons(power_law)
    ra_iso.AddThermalTargetPhotons(t,edens)


    # Calculate the resulting spectra
    e = np.logspace(-9,3,ebins) * gp.TeV_to_erg
    ra_iso.CalculateDifferentialPhotonSpectrum(e)
    ra_aniso.CalculateDifferentialPhotonSpectrum(e)

    # Get the SEDs for plotting
    icsed_iso = np.array(ra_iso.GetICSED())
    icsed_aniso = np.array(ra_aniso.GetICSED())

    f = plt.figure(figsize=(5,5))
    plt.loglog(icsed_iso[:,0],icsed_iso[:,1],label="iso")
    plt.loglog(icsed_aniso[:,0],icsed_aniso[:,1],label="aniso")
    plt.legend()
    plt.ylabel("E^2 dN/dE (arb)")
    plt.xlabel("E (TeV)")
    plt.savefig("spectra1_4.png")

    f = plt.figure(figsize=(5,5))
    plt.semilogx(icsed_iso[:,0],icsed_aniso[:,1]/icsed_iso[:,1],label="ratio aniso / iso")
    plt.legend()
    plt.xlabel("E (TeV)")
    plt.savefig("spectra2_4.png")


