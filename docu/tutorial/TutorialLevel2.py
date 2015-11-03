#!/usr/local/bin/python

import sys
import os
sys.path.append(os.path.abspath('./lib'))
import gamerapy
import numpy as np
import math
import matplotlib.pyplot as plt
import ConfigParser


## This is a small script that performs the evolution of constantly injected
## electrons with luminosity 'lum' between emin and emax and shaped in a
## broken power-law with 'low_index', 'break_energy' and 'high_index' over a
## time 'age'.
## Of course, a single power law is a special case and can either be achieved
## by setting both indices to the same value or by setting low_index = 0.
## electrons are subjected to Synchr., IC and Bremsstrahlung losses, determined
## by density dens, B-Field bfield, radiation field e-density edensir and 
## temperature tir. The CMB is hardcoded, just like in real life.
## Finally, the radiation spectrum is calculated given a certain distance 
## and plotted into outfile.
## All variables are detemined in the configFile and have the following units:
## Units are: 
## tir - K
## edensir - eV/cm^3
## dens - cm^-3
## bfield - G
## dist - pc
## emax - TeV
## ebreak - TeV
## age - yrs

if __name__ == "__main__":

  # Read in parameter file
  configFile = os.path.abspath(sys.argv[1])
  configParser = ConfigParser.RawConfigParser()
  configParser.read(configFile)
  lum = float(configParser.get('Parameters','Luminosity'))
  age = float(configParser.get('Parameters','Age'))
  dist = gamerapy.pc_to_cm*float(configParser.get('Parameters','Distance'))
  dens = float(configParser.get('Parameters','AmbientDensity'))
  bfield = float(configParser.get('Parameters','BField'))
  t = float(configParser.get('Parameters','tRAD'))
  e = gamerapy.TeV_to_erg*1.e-12*float(configParser.get('Parameters','edensRAD'))
  ebins = float(configParser.get('Parameters','Ebins'))
  emax = gamerapy.TeV_to_erg*float(configParser.get('Parameters','Emax'))
  emin = gamerapy.TeV_to_erg*float(configParser.get('Parameters','Emin'))
  spind = float(configParser.get('Parameters','SpectralIndex'))
  outfile = configParser.get('Files','outfile')

  fr = gamerapy.Radiation()
  fp = gamerapy.Particles()
  fu = gamerapy.Utils()
  fu.DrawGamera()
  # set particle stuff
  fp.SetLuminosity(lum)
  fp.SetBField(bfield)
  fp.SetEmax(emax)
  fp.SetEmin(emin)
  fp.SetSpectralIndex(spind)
  fp.SetEnergyBinsForNumericalSolver(ebins)
  fp.SetAmbientDensity(dens)
  fp.SetAge(age) 

  # set radiation stuff
  fr.SetDistance(dist)
  fr.SetAmbientDensity(fp.GetAmbientDensity())
  fr.SetBField(fp.GetBField())
  fr.AddThermalTargetPhotons(2.7,0.25*1.602*1.e-12)
  fr.AddThermalTargetPhotons(t,e)
  fr.CreateICLossLookup()
  fp.SetICLossLookup(fr.GetICLossLookup())

  # calculate stuff 
  fp.CalculateParticleSpectrum("electrons")
  fr.SetElectrons(fp.GetParticleSpectrum())
  fp.SetLuminosity(1.e2*lum)
  fp.SetEmax(gamerapy.TeV_to_erg*1.e3)
  fp.SetSpectralIndex(spind-0.1)
  fp.CalculateParticleSpectrum("protons")
  fr.SetProtons(fp.GetParticleSpectrum())
  fr.CalculateDifferentialPhotonSpectrum()

  ElectronSED = np.array(fr.GetElectronSED())
  ProtonSED = np.array(fr.GetProtonSED())
  TotalSED = np.array(fr.GetTotalSED())
  ICSED = np.array(fr.GetICSED())
  BremsSED = np.array(fr.GetBremsstrahlungSED())
  SynchSED = np.array(fr.GetSynchrotronSED())
  PPSED = np.array(fr.GetPPSED())

  ## plot stuff ##
  f, (ax1, ax2) = plt.subplots(1, 2, sharey=False,figsize=(15, 6))  
  ax1.set_yscale("log")
  ax1.set_xscale("log")
  ax2.set_yscale("log")
  ax2.set_xscale("log")
  ax1.set_title('particles')
  ax2.set_title('gammas')
  ax1.set_ylabel(r'$E^{2} \mathrm{d}N/\mathrm{d}E(\mathrm{erg})$', fontsize=13)
  ax1.set_xlabel(r'$E(\mathrm{TeV})$', fontsize=13)
  ax2.set_ylabel(r'$\nu F_{\nu}(\mathrm{erg} \cdot \mathrm{cm}^{-2} \cdot \mathrm{s}^{-1})$', fontsize=13)
  ax2.set_xlabel(r'$E(\mathrm{TeV})$', fontsize=13)
  ax1.set_xlim([0.7*min(min(ElectronSED[:,0]),min(ProtonSED[:,0])),
                1.5*max(max(ElectronSED[:,0]),max(ProtonSED[:,0]))])
  ax1.set_ylim([0.7*min(min(ElectronSED[:,1]),min(ProtonSED[:,1])),
                1.5*max(max(ElectronSED[:,1]),max(ProtonSED[:,1]))])
  ax2.set_xlim([0.7*min(TotalSED[:,0]),1.5*max(TotalSED[:,0])])
  ax2.set_ylim([0.7*min(TotalSED[:,1]),1.5*max(TotalSED[:,1])])
  ax1.plot(ElectronSED[:,0],ElectronSED[:,1],color='blue',alpha=1.,label="electrons")
  ax1.plot(ProtonSED[:,0],ProtonSED[:,1],color='red',alpha=1.,label="protons")
  ax1.legend(prop={'size':12},loc="lower left")
  ax2.plot(TotalSED[:,0],TotalSED[:,1],color='black',alpha=1.)
  ax2.plot(PPSED[:,0],PPSED[:,1],color='green',alpha=.6,label="pion-decay")
  ax2.plot(SynchSED[:,0],SynchSED[:,1],color='blue',alpha=.6,label="Synch")
  ax2.plot(BremsSED[:,0],BremsSED[:,1],color='orange',alpha=.6,label="Brems")
  ax2.plot(ICSED[:,0],ICSED[:,1],color='red',alpha=.6,label="IC")
  ax2.legend(prop={'size':12},loc="upper right")
  f.savefig(outfile)

