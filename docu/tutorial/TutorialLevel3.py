#!/usr/local/bin/python

import sys
import os
sys.path.append(os.path.abspath('./lib'))
import gappa as gp
import numpy as np
import math
import matplotlib.pyplot as plt
import ConfigParser

def GetPulsarSpindown(tc, age, l0):

  t = np.logspace(0,math.log10(1.e6*age),300)
  puls = l0/(1.+t/tc)**2

  puls = np.vstack((t, puls)).T
  return puls

if __name__ == "__main__":

  # Read in parameter file
  configFile = os.path.abspath(sys.argv[1])
  configParser = ConfigParser.RawConfigParser()
  configParser.read(configFile)
  lum = float(configParser.get('Parameters','Luminosity'))
  age = float(configParser.get('Parameters','Age'))
  tc = float(configParser.get('Parameters','CharAge'))
  dist = float(configParser.get('Parameters','Distance'))
  dens = float(configParser.get('Parameters','AmbientDensity'))
  bfield = float(configParser.get('Parameters','BField'))
  t = float(configParser.get('Parameters','tRAD'))
  e = gp.eV_to_erg*float(configParser.get('Parameters','edensRAD'))
  ebins = float(configParser.get('Parameters','Ebins'))
  emax = gp.TeV_to_erg*float(configParser.get('Parameters','Emax'))
  emin = gp.TeV_to_erg*float(configParser.get('Parameters','Emin'))
  spind = float(configParser.get('Parameters','SpectralIndex'))
  outfile = configParser.get('Files','outfile')

  puls = np.array(GetPulsarSpindown(tc,age,lum))

  fr = gp.Radiation()
  fp = gp.Particles()
  fu = gp.Utils()
  fu.DrawGamera()
  # set particle stuff
  fp.SetLuminosityLookup(puls)
  fp.SetBField(bfield)
  fp.SetEmax(emax)
  fp.SetEmin(emin)
  fp.SetSpectralIndex(spind)
  fp.SetEnergyBins(ebins)
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
  fr.CalculateDifferentialPhotonSpectrum()

  ElectronSED = np.array(fr.GetElectronSED())
  TotalSED = np.array(fr.GetTotalSED())
  ICSED = np.array(fr.GetICSED())
  BremsSED = np.array(fr.GetBremsstrahlungSED())
  SynchSED = np.array(fr.GetSynchrotronSED())

  ## plot stuff ##
  f, (ax1, ax2, ax3) = plt.subplots(3, 1, sharey=False,figsize=(7, 15))
  ax1.set_yscale("log")
  ax1.set_xscale("log")
  ax2.set_yscale("log")
  ax2.set_xscale("log")
  ax3.set_yscale("log")
  ax3.set_xscale("log")
  ax1.set_ylabel(r'$E^{2} \mathrm{d}N/\mathrm{d}E(\mathrm{erg})$', fontsize=13)
  ax1.set_xlabel(r'$E (\mathrm{TeV})$', fontsize=13)
  ax2.set_ylabel(r'$\nu F_{\nu}(\mathrm{erg} \cdot \mathrm{cm}^{-2} \cdot \mathrm{s}^{-1})$', fontsize=13)
  ax2.set_xlabel(r'$E (\mathrm{TeV})$', fontsize=13)
  ax3.set_ylabel(r'$L (\mathrm{erg} \cdot \mathrm{s}^{-1})$', fontsize=13)
  ax3.set_xlabel(r'$t (\mathrm{yrs})$', fontsize=13)
  ax1.set_xlim([0.7*min(ElectronSED[:,0]),1.5*max(ElectronSED[:,0])])
  ax1.set_ylim([0.7*min(ElectronSED[:,1]),1.5*max(ElectronSED[:,1])])
  ax2.set_xlim([0.7*min(TotalSED[:,0]),1.5*max(TotalSED[:,0])])
  ax2.set_ylim([0.7*min(TotalSED[:,1]),1.5*max(TotalSED[:,1])])
  ax3.set_xlim([0.7*min(puls[:,0]),1.5*max(puls[:,0])])
  ax3.set_ylim([0.7*min(puls[:,1]),1.5*max(puls[:,1])])
  ax1.plot(ElectronSED[:,0],ElectronSED[:,1],color='black',alpha=1.,label="electrons")
  ax1.legend(prop={'size':12},loc="lower left")
  ax2.plot(TotalSED[:,0],TotalSED[:,1],color='black',alpha=1.)
  ax2.plot(SynchSED[:,0],SynchSED[:,1],color='blue',alpha=.6,label="Synch")
  ax2.plot(BremsSED[:,0],BremsSED[:,1],color='orange',alpha=.6,label="Brems")
  ax2.plot(ICSED[:,0],ICSED[:,1],color='red',alpha=.6,label="IC")
  ax2.legend(title="gammas",prop={'size':12},loc="lower left")
  ax3.plot(puls[:,0],puls[:,1],color='black',alpha=1.,label="luminosity")
  ax3.legend(prop={'size':12},loc="lower left")
  f.savefig(outfile)
