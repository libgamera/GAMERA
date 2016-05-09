#!/usr/local/bin/python

import sys
import os
sys.path.append(os.path.abspath('./lib'))
import gappa as gp
import numpy as np
import math
import matplotlib.pyplot as plt
import ConfigParser

def GetSpindownStuff(tc, age, l0, b0, emax0, n0):

  t = np.logspace(0,math.log10(1000.*age),300)
  lumt = []
  emaxt = []
  bt = []
  nt = []
  n = 0
  for i in t:
    l = l0/math.pow(1.+i/tc,2.)
    btt = b0*math.sqrt(l/l0)*(1.+0.5*math.sin(0.1*n*3.14))
    lumtt = l*(1.+0.5*math.cos(0.1*n*3.14))
    emaxtt = emax0*math.pow(l/l0,0.25)*(1.+0.5*math.sin(0.05*n*3.14))
    ntt = n0*math.pow(l/l0,0.25)*(1.+0.5*math.cos(0.05*n*3.14))

    bt.append([])
    bt[n].append(i)
    bt[n].append(btt)
    lumt.append([])
    lumt[n].append(i)
    lumt[n].append(lumtt)
    emaxt.append([])
    emaxt[n].append(i)
    emaxt[n].append(emaxtt)
    nt.append([])
    nt[n].append(i)
    nt[n].append(ntt)

    n = n+1
  return lumt,bt,emaxt,nt

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
  ebins = int(configParser.get('Parameters','Ebins'))
  emax = gp.TeV_to_erg*float(configParser.get('Parameters','Emax'))
  emin = gp.TeV_to_erg*float(configParser.get('Parameters','Emin'))
  ebreak = gp.TeV_to_erg*float(configParser.get('Parameters','Ebreak'))
  spindlow = float(configParser.get('Parameters','SpectralIndexLow'))
  spindhigh = float(configParser.get('Parameters','SpectralIndexHigh'))
  outfile = configParser.get('Files','outfile')

  lumt,bt,emaxt,denst = GetSpindownStuff(tc, age, lum, bfield, emax, dens)
  lumt = np.array(lumt)
  bt = np.array(bt)
  emaxt = np.array(emaxt)
  denst = np.array(denst)

  fr = gp.Radiation()
  fp = gp.Particles()
  fu = gp.Utils()
  fu.DrawGamera()
  # set particle stuff
  fp.SetLuminosityLookup(lumt)
  fp.SetBFieldLookup(bt)
  fp.SetEmaxLookup(emaxt)
  fp.SetAmbientDensityLookup(denst)
  fp.SetEmin(emin)
  fp.SetBreakEnergy(ebreak)
  fp.SetLowSpectralIndex(spindlow)
  fp.SetSpectralIndex(spindhigh)
  fp.SetAge(age)

  # set radiation stuff
  fr.SetDistance(dist)
  fr.AddThermalTargetPhotons(2.7,0.25*1.602*1.e-12)
  fr.AddThermalTargetPhotons(t,e)
  fr.CreateICLossLookup()
  fp.SetICLossLookup(fr.GetICLossLookup())

  fp.ToggleQuietMode()
  fr.ToggleQuietMode()

  f, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6, 1, sharey=False,figsize=(7, 30))
  ax1.set_yscale("log")
  ax1.set_xscale("log")
  ax2.set_yscale("log")
  ax2.set_xscale("log")

  # calculate stuff
  t = np.logspace(1,math.log10(age),20)
  n = 0
  for i in t:
    print i,n
    fp.SetAge(i)
    fr.SetAmbientDensity(fp.GetAmbientDensity())
    fr.SetBField(fp.GetBField())
    fp.CalculateElectronSpectrum(ebins)
    fr.SetElectrons(fp.GetParticleSpectrum())
    fr.CalculateDifferentialPhotonSpectrum()


    ElectronSED = np.array(fr.GetElectronSED())
    TotalSED = np.array(fr.GetTotalSED())
    ax1.plot(ElectronSED[:,0],ElectronSED[:,1],color='black',alpha=1.-n,label=str(round(i,0)))
    ax2.plot(TotalSED[:,0],TotalSED[:,1],color='black',alpha=1.-n,label=str(round(i,0)))
    n = n+0.04

  ax1.legend(title="age",prop={'size':6},loc="lower left")
  ax2.legend(title="age",prop={'size':6},loc="lower left")
  ## plot stuff ##
  ax3.set_yscale("log")
  ax3.set_xscale("log")
  ax4.set_yscale("log")
  ax4.set_xscale("log")
  ax5.set_yscale("log")
  ax5.set_xscale("log")
  ax6.set_yscale("log")
  ax6.set_xscale("log")
  ax1.set_ylabel(r'$E^{2} \mathrm{d}N/\mathrm{d}E(\mathrm{erg})$', fontsize=13)
  ax1.set_xlabel(r'$E (\mathrm{TeV})$', fontsize=13)
  ax2.set_ylabel(r'$\nu F_{\nu}(\mathrm{erg} \cdot \mathrm{cm}^{-2} \cdot \mathrm{s}^{-1})$', fontsize=13)
  ax2.set_xlabel(r'$E (\mathrm{TeV})$', fontsize=13)
  ax3.set_ylabel(r'$L (\mathrm{erg} \cdot \mathrm{s}^{-1})$', fontsize=13)
  ax3.set_xlabel(r'$t (\mathrm{yrs})$', fontsize=13)
  ax4.set_ylabel(r'$B (\mu \mathrm{G})$', fontsize=13)
  ax4.set_xlabel(r'$t (\mathrm{yrs})$', fontsize=13)
  ax5.set_ylabel(r'$E_{max} (\mathrm{erg})$', fontsize=13)
  ax5.set_xlabel(r'$t (\mathrm{yrs})$', fontsize=13)
  ax6.set_ylabel(r'$n (\mathrm{cm}^{-3})$', fontsize=13)
  ax6.set_xlabel(r'$t (\mathrm{yrs})$', fontsize=13)
  ax1.set_xlim([0.7*min(ElectronSED[:,0]),1000.])
  ax1.set_ylim([0.01*min(ElectronSED[:,1]),1.5*max(ElectronSED[:,1])])
  ax2.set_xlim([0.7*min(TotalSED[:,0]),1000.])
  ax2.set_ylim([1.e-20,5.e-10])
  ax3.set_xlim([0.7*min(lumt[:,0]),1.5*max(lumt[:,0])])
  ax3.set_ylim([0.7*min(lumt[:,1]),1.5*max(lumt[:,1])])
  ax4.set_xlim([0.7*min(bt[:,0]),1.5*max(bt[:,0])])
  ax4.set_ylim([0.7*min(bt[:,1]),1.5*max(bt[:,1])])
  ax5.set_xlim([0.7*min(emaxt[:,0]),1.5*max(emaxt[:,0])])
  ax5.set_ylim([0.7*min(emaxt[:,1]),1.5*max(emaxt[:,1])])
  ax6.set_xlim([0.7*min(denst[:,0]),1.5*max(denst[:,0])])
  ax6.set_ylim([0.7*min(denst[:,1]),1.5*max(denst[:,1])])
  ax3.plot(lumt[:,0],lumt[:,1],color='black',alpha=1.,label="luminosity")
  ax3.legend(prop={'size':12},loc="lower left")
  ax4.plot(bt[:,0],bt[:,1],color='black',alpha=1.,label="B-field")
  ax4.legend(prop={'size':12},loc="lower left")
  ax5.plot(emaxt[:,0],emaxt[:,1],color='black',alpha=1.,label="E_{max}")
  ax5.legend(prop={'size':12},loc="lower left")
  ax6.plot(denst[:,0],denst[:,1],color='black',alpha=1.,label="density")
  ax6.legend(prop={'size':12},loc="lower left")
  f.savefig(outfile)
