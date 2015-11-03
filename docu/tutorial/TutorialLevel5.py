#!/usr/local/bin/python

import sys
import os
sys.path.append(os.path.abspath('./lib'))
import gamerapy
import numpy as np
import math
import matplotlib.pyplot as plt
import ConfigParser

global lum0, age, tc, mej, e0, etab, eps

def CalculateTimeDependentStuff():
  t = np.logspace(0,math.log10(1.e6*age),80)
  gammap = 1.3333
  vej = math.sqrt(10.*e0/(3.*mej))
  c = math.pow((6./(15.*(gammap-1.)))+289./240.,-0.2);

  lum = (1.-etab)*lum0*(1.+t/tc)**(-1.*(brind+1.)/(brind-1.))
  emax = 3.*eps*gamerapy.el_charge*np.sqrt(etab*lum/((1.-etab)*gamerapy.c_speed))
  r = c*(lum0*t*gamerapy.yr_to_sec/e0)**0.2 * vej*t*gamerapy.yr_to_sec
  v = 1.2*r/(gamerapy.yr_to_sec*t)
  b = np.sqrt(gamerapy.yr_to_sec*etab*6./r**4 * np.concatenate(([0], ((lum * r)[1:] * np.diff(t)).cumsum())))

  lum = np.vstack((t, lum)).T
  b = np.vstack((t, b)).T
  emax = np.vstack((t, emax)).T
  r = np.vstack((t, r)).T
  v = np.vstack((t, v)).T

  return lum, b, emax, r, v


if __name__ == "__main__":
  # Read in parameter file
  configFile = os.path.abspath(sys.argv[1])
  configParser = ConfigParser.RawConfigParser()
  configParser.read(configFile)

  f, (ax1, ax2, ax3, ax4, ax5, ax6, ax7) = plt.subplots(7, 1, sharey=False,figsize=(7, 35))
  lum0 = configParser.getfloat('Parameters','InitialLuminosity')
  age = configParser.getfloat('Parameters','Age')
  tc = configParser.getfloat('Parameters','CharAge')
  dist = gamerapy.pc_to_cm*configParser.getfloat('Parameters','Distance')
  dens = configParser.getfloat('Parameters','AmbientDensity')
  tNIR = configParser.getfloat('Parameters','tNIR')
  eNIR = gamerapy.TeV_to_erg*1.e-12*configParser.getfloat('Parameters','edensNIR')
  tFIR = configParser.getfloat('Parameters','tFIR')
  eFIR = gamerapy.TeV_to_erg*1.e-12*configParser.getfloat('Parameters','edensFIR')
  ebins = configParser.getfloat('Parameters','Ebins')
  emin = gamerapy.TeV_to_erg*configParser.getfloat('Parameters','Emin')
  ebreak = gamerapy.TeV_to_erg*configParser.getfloat('Parameters','Ebreak')
  spindlow = configParser.getfloat('Parameters','SpectralIndexLow')
  spindhigh = configParser.getfloat('Parameters','SpectralIndexHigh')
  mej = gamerapy.mSol*configParser.getfloat('Parameters','Mej')
  e0 = configParser.getfloat('Parameters','E0')
  etab = configParser.getfloat('Parameters','etaB')
  eps = configParser.getfloat('Parameters','epsilon')
  brind = configParser.getfloat('Parameters','BrakingIndex')
  outfile = configParser.get('Files','outfile')

  lumt,bt,emaxt,r,v = CalculateTimeDependentStuff()

  fu = gamerapy.Utils()
  fu.DrawGamera()

  # set particle stuff
  fp = gamerapy.Particles()
  fp.SetLuminosityLookup(lumt)
  fp.SetBFieldLookup(bt)
  fp.SetEmaxLookup(emaxt)
  fp.SetRadiusLookup(r)
  fp.SetVelocityLookup(v)
  fp.SetAmbientDensity(dens)
  fp.SetEmin(emin)
  fp.SetBreakEnergy(ebreak)
  fp.SetLowSpectralIndex(spindlow)
  fp.SetSpectralIndex(spindhigh)
  fp.SetEnergyBinsForNumericalSolver(ebins)
  fp.SetAge(age)

  # set radiation stuff
  fr = gamerapy.Radiation()
  fr.SetDistance(dist)
  fr.AddThermalTargetPhotons(2.7,0.25*1.602*1.e-12)
  fr.AddThermalTargetPhotons(tFIR,eFIR)
  fr.AddThermalTargetPhotons(tNIR,eNIR)
  fr.SetAmbientDensity(fp.GetAmbientDensity())
  fr.SetSynchrotronEmissionModel(1)
  fr.CreateICLossLookup()
  fp.SetICLossLookup(fr.GetICLossLookup())

  # calculate stuff
  fp.SetAge(age)
  fr.SetBField(fp.GetBField())
  fp.CalculateParticleSpectrum("electrons")
  fr.SetElectrons(fp.GetParticleSpectrum())
  fr.AddSSCTargetPhotons(fp.GetRadius())
  fr.CalculateDifferentialPhotonSpectrum()
  ElectronSED = np.array(fr.GetElectronSED())
  TotalSED = np.array(fr.GetTotalSED())
  ICSED = np.array(fr.GetICSED())
  BremsSED = np.array(fr.GetBremsstrahlungSED())
  SynchSED = np.array(fr.GetSynchrotronSED())

  # plot stuff
  ax1.loglog(ElectronSED[:,0],ElectronSED[:,1],color='black',alpha=1.,label=str(age))
  ax1.set_xlim([0.7*min(ElectronSED[:,0]),2.e4])
  ax1.set_ylim([1.e42,1.e49])
  ax1.legend(title="age",prop={'size':6},loc="lower left")
  ax1.set_ylabel(r'$E^{2} \mathrm{d}N/\mathrm{d}E(\mathrm{erg})$', fontsize=13)
  ax1.set_xlabel(r'$E (\mathrm{TeV})$', fontsize=13)
  ax1.set_xlim([0.7*min(ElectronSED[:,0]),2.e4])

  ax2.loglog(TotalSED[:,0],TotalSED[:,1],color='black',alpha=1.,label=str(age))
  ax2.set_xlim([0.7*min(TotalSED[:,0]),2.e3])
  ax2.set_ylim([1.e-16,1.e-6])
  ax2.legend(title="age",prop={'size':6},loc="lower left")
  ax2.set_ylabel(r'$\nu F_{\nu}(\mathrm{erg} \cdot \mathrm{cm}^{-2} \cdot \mathrm{s}^{-1})$', fontsize=13)
  ax2.set_xlabel(r'$E (\mathrm{TeV})$', fontsize=13)

  ax3.loglog(lumt[:,0],lumt[:,1],color='black',alpha=1.,label="luminosity")
  ax3.set_xlim([2.,1.5*age])
  ax3.set_ylim([1.e36,1.1*max(lumt[:,1])])
  ax3.legend(prop={'size':12},loc="lower left")
  ax3.set_ylabel(r'$L (\mathrm{erg} \cdot \mathrm{s}^{-1})$', fontsize=13)
  ax3.set_xlabel(r'$t (\mathrm{yrs})$', fontsize=13)

  ax4.loglog(bt[:,0],bt[:,1],color='black',alpha=1.,label="B-field")
  ax4.set_xlim([2.,1.5*age])
  ax4.set_ylim([1.e-6,1.1*max(bt[:,1])])
  ax4.legend(prop={'size':12},loc="lower left")
  ax4.set_ylabel(r'$B (\mu \mathrm{G})$', fontsize=13)
  ax4.set_xlabel(r'$t (\mathrm{yrs})$', fontsize=13)

  ax5.semilogx(emaxt[:,0],emaxt[:,1],color='black',alpha=1.,label=r'$E_{max}$')
  ax5.legend(prop={'size':12},loc="lower left")
  ax5.set_xlim([2.,1.5*age])
  ax5.set_ylim([0.7*min(emaxt[:,1]),1.1*max(emaxt[:,1])])
  ax5.set_ylabel(r'$E_{max} (\mathrm{erg})$', fontsize=13)
  ax5.set_xlabel(r'$t (\mathrm{yrs})$', fontsize=13)

  ax6.semilogx(r[:,0],r[:,1]/gamerapy.pc_to_cm,color='black',alpha=1.,label="Radius")
  ax6.set_xlim([2.,1.5*age])
  ax6.set_ylim([0.,100])
  ax6.legend(prop={'size':12},loc="lower left")
  ax6.set_ylabel(r'$R (\mathrm{pc})$', fontsize=13)
  ax6.set_xlabel(r'$t (\mathrm{yrs})$', fontsize=13)

  ax7.semilogx(v[:,0],v[:,1],color='black',alpha=1.,label="Velocity")
  ax7.set_xlim([2.,1.5*age])
  ax7.set_ylim([0.,1.e9])
  ax7.legend(prop={'size':12},loc="lower right")
  ax7.set_ylabel(r'$V (\mathrm{cm/s})$', fontsize=13)
  ax7.set_xlabel(r'$t (\mathrm{yrs})$', fontsize=13)

  f.savefig(outfile)
