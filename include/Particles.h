#ifndef _PARTICLES_
#define _PARTICLES_

#include "Utils.h"
#include "Radiation.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <time.h>


using namespace std;

/**  @file Particles.C
 *   @author Joachim Hahn
 *   @short Class that calculates particle spectra
 *   This class allows for the time-dependent model of particle spectra,
 *   using as a center piece a numerical algorithm (a piece-wise linear
 *   advection scheme in energy space. Thus, in principle any time evolution of
 *   energy losses and particle injection can be assumed and the resulting
 * particle
 *   distribution can be calculated.
 *   What this code cannot do is the simultaneous spatial evolution of
 * particles, e.g.
 *   energy dependent diffusion. Here, no assumptions on the spatial
 * distribution of
 *   particles is made whatsoever.
 */

class Particles {
  typedef double (Particles::*fPointer)(double, void *);

  static double evaluate(double x, void* params);
  static fPointer _funcPtr;
  static Particles *_partPtr;


struct timespec time0, time1, time2, time3;
 private:
  Utils *fUtils;
  Radiation *fRadiation;
  int Type;  ///< integer indicating particle type. supported: 0-electrons and
             ///1-protons
  int METHOD; /// calculation method. This is determined automatically, but can
              /// be force set via SetSolverMethod(). However, if an inappropriate
              /// method is chosen manually, the result might be wrong or
              /// the program might crash.
              /// 0: grid solver,
              /// 1: semi-analytical method in the presence of constant losses
              /// 2: no losses (simply adds up according to luminosity)
  bool FASTMODE; /// grid solver: sacrifice accuracy for speed. -> true: disables slope limiters. default: false

  bool ESCTIME;/// boolean that indicates whether an escape time has been set. This will be automatically set to 'true' if any method for escape is used. default: false
  double theta;           ///< CR acceleration efficiency
  double SpectralIndex;   ///< spectral index of injected particles
  double SpectralIndex2;  ///< low-energy spectral index for broken power law
                          ///injection spectrum
  double SpectralIndexConstant;   ///< spectral index of injected particles
                                  ///manually set to a constant value
  double SpectralIndex2Constant;  ///< low-energy spectral index for broken
                                  ///power law injection spectrum manually set
                                  ///to a constant value
  double Kep;                     ///< electron-to-proton _number_ ratio
  double Tmin;  ///< optimised minimal time from which to start particle
                ///injection (given by the 'maximal' starting time that doesn't
                ///change results)
  double TminInternal;  ///< internal minimal time for which all lookups
                        ///(B,N,R,V etc) are filled. If no lookups are provided
                        ///(constant losses etc), it is set to 0.
  double TmaxInternal;  ///< Same as TminInternal but here it is the maximal
                        ///time
  double TActual;  ///< Actual time of iteration time step
  double TIterationStart; ///< start time when using grid solver
  double EminInternal;  ///< internal parameter that is relevant to the grid-
                        /// solver. Essentially, it is used to determine the
                        /// starting time of the iteration. The starting time
                        /// is then determined in a way that all particles that
                        /// are injected before that time are cooled to an
                        /// energy of maximal EminInternal after a time t=age.
  double MinELookup;    ///< Lower energy bound of custom injection lookup.
                        ///  If time evolution of custom spectrum is given, this
                        ///  is reset at every T = Age step.
  double MaxELookup;    ///< Upper energy bound of custom injection lookup.
                        ///  If time evolution of custom spectrum is given, this
                        ///  is reset at every T = Age step.
  double TminLookup;    ///< Lower time bound of custom injection lookup.
  double TmaxLookup;    ///< Higher time bound of custom injection lookup.
  double BField;        ///< amplified B-field immediately dowmstream
  double N;             ///< ambient density
  double R;             ///< source extension (cm)
  double V;             ///< source expansion velocity (cm/s)
  double Lum;           ///< energy injection rate at a given time
  double Emin;  ///< internally and dynamically determined lower energy boundary
                ///within which the particle spectrum is calculated
  bool EminConstantForNormalisationOnly;  ///< flag that, if set to 'true'
                                          ///indicates that emin will only be
                                          ///used in the calculation of the
                                          ///particle norm, but not the spectral
                                          ///range.
  double Emax;  ///< internally and dynamically determined upper energy boundary
                ///within which the particle spectrum is calculated
  double EminConstant;  ///< Constantly set lower energy boundary within which
                        ///the particle spectrum is calculated (optional)
  double EmaxConstant;  ///< Constantly set upper energy boundary within which
                        ///the particle spectrum is calculated (optional)
  double TminExternal;  ///< Constantly set lower time boundary at which the
                        ///particle iteration is started
  double BConstant;     ///< Constantly set constant B-Field
  double NConstant;     ///< Constantly set constant ambient density
  int iclosslookupbins; ///< bins of the IC-Loss lookup
  double LumConstant;   ///< Constantly set source luminosity
  double RConstant;     ///< Constantly set source extension
  double VConstant;     ///< Constantly set source expansion speed
  double eMax;          ///< maximum energy of particles the shock can contain
  double eBreak;  ///< break energy for a broken power-law particle injection
                  ///spectrum
  double CustomInjectionNorm;///< energy content in injection spectrum before normalisation
  double eBreakConstant;  ///< break energy for a broken power-law particle
                          ///injection spectrum manually set to a constant value
  double energyMarginFactor;  ///< energy safety margin above the spectral
                              ///cut-off. E.g. energyMarginFactor = 1.e-3 means
                              ///that the upper boundary of the energy spectrum
                              ///corresponds to the point where the has dropped
                              ///to 1/1000 exponential
  double energyMargin;  ///< this is the dynamically determined margin that
                        ///results from energyMarginFactor and CutOffFactor. It
                        ///determines the upper energy boundary of the particle
                        ///spectrum as upper_boundary = emax*energyMargin.
  double CutOffFactor;  ///< factor determining the strength of the exponential
                        ///cut-off in the particle injection spectrum:
                        ///f~exp[-(e/ecut)^CutOffFactor]
  bool sharpEnergyCut;  ///< toggle a sharp cut-off at emax
  double Age;           ///< source age
  double eBreakS2;      ///< constant used in 'SourceSpectrum'. Used for speed
                        ///optimisation
  double eBreak2mS2;    ///< constant used in 'SourceSpectrum'. Used for speed
                        ///optimisation
  double eBreakS;       ///< constant used in 'SourceSpectrum'. Used for speed
                        ///optimisation
  double eBreak2mS;     ///< constant used in 'SourceSpectrum'. Used for speed
                        ///optimisation
  double emin2mS2;      ///< constant used in 'SourceSpectrum'. Used for speed
                        ///optimisation
  double emin2mS;       ///< constant used in 'SourceSpectrum'. Used for speed
                        ///optimisation
  double emineBreak2mS2;  ///< constant used in 'SourceSpectrum'. Used for speed
                          ///optimisation
  double eBreak2mSInd2;   ///< constant used in 'SourceSpectrum'. Used for speed
                          ///optimisation
  double emin2mSInd2;     ///< constant used in 'SourceSpectrum'. Used for speed
                          ///optimisation
  double emineBreak2mSInd2;  ///< constant used in 'SourceSpectrum'. Used for
                             ///speed optimisation
  double fS;           ///< constant used in 'SourceSpectrum'. Used for speed
                       ///optimisation
  double fS2;          ///< constant used in 'SourceSpectrum'. Used for speed
                       ///optimisation
  double bremsl_epf;   ///< constant used in 'EnergyLossRate'. Used for speed
                       ///optimisation
  double bremsl_eef;   ///< constant used in 'EnergyLossRate'. Used for speed
                       ///optimisation

  double synchl,icl,adl,bremsl,ionization,ppcol; ///< loss rates for synchrotron, IC, adiabatic exp.,
                                                 /// bremsstrahlug and ionizarion and pp collision for protons
  bool IONIZATION;  ///< Boolean to decide if to compute or not the ionization losses.

  bool logarithmicCRLumLookupTimeBins;  ///< boolean indicating that time steps
                                        ///in CRLumLookup are logarithmic. In
                                        ///this case, the precise entry in the
                                        ///Lookup up for a given time can be
                                        ///calculated rather than iterated
  bool linearCRLumLookupTimeBins;  ///< boolean indicating that time steps in
                                   ///CRLumLookup are linear. In this case, the
                                   ///precise entry in the Lookup up for a given
                                   ///time can be calculated rather than
                                   ///iterated
  double CRLumLookupTimeMin;  ///< starting time (earliest entry) of CRLumLookup
  double CRLumLookupTimeMax;  ///< ending time (last entry) of CRLumLookup
  double logCRLumLookupTimeMin;         ///< log10 of CRLumLookupTimeMin (speed
                                        ///optimisation)
  double logCRLumLookupTimeMax;         ///< log10 of CRLumLookupTimeMax (speed
                                        ///optimisation)
  unsigned int CRLumLookupSize;         ///< entries of CRLumLookup
  double CRLumLookupdeltat;             ///< time bin size in CRLumLookup
  vector<vector<double> > CRLumLookup;  ///< 2D-vector containing quantities
                                        ///needed to calculate particle spectra.
                                        ///Each line comprises the quantities at
                                        ///a gives time. This allows for
                                        ///time-dependent modeling. Format:
                                        ///time(yrs) - e_dot( erg/s, luminosity
                                        ///of source put into acc.particles -
                                        ///radius (cm) - speed (cm/s) -
                                        ///ambient_density (cm^-3) -
                                        ///ambient_density2 (cm^-3) - magn.Field
                                        ///(G) - max. energy of injected protons
                                        ///(erg) - max. energy of injected
                                        ///electrons (erg) - escape time

  vector<vector<double> > ICLossVector, LumVector, NVector, BVector, eMaxVector,
      escapeTimeVectorEdep, escapeTimeVectorTdep, RVector, VVector,tempspec,
      CustomInjectionSpectrumTimeEvolutionVector,
      EscapeTimeEnergyTimeEvolutionVector;
  gsl_spline *ICLossLookup, *LumLookup, *NLookup, *BFieldLookup, *eMaxLookup,
      *escapeTimeLookupEdep, *escapeTimeLookupTdep, *RLookup, *VLookup, *energyTrajectory,
      *energyTrajectoryInverse,*CustomInjectionSpectrum;
  interp2d_spline *CustomInjectionSpectrumTimeEvolution,
      *EscapeTimeEnergyTimeEvolution;
  gsl_interp_accel *accIC, *accLum, *accN, *accBField, *acceMax, *accescapeTimeEdep,
      *accescapeTimeTdep,*accR, *accV, *accTr, *accTrInv,*accCustInj,
      *taccsp,*eaccsp,*taccesc,*eaccesc;
  double SourceSpectrum(double e);  ///< particle injection spectrum
  double EscapeTime(double e, double t);  ///< particle escape time spectrum
  vector<double> energyAxis;        ///< energy axis for numerical integrator
  vector<vector<double> > grid;     ///< 2D grid on which solving takes place
  vector<vector<double> > diffusegrid;  ///< 2D grid of the streaming particles
  double ebins;                         ///< energy bins of grid (default: 100)
  void DetermineTMin(double emin,
                     double &tmin);  ///< determine the minimum time from where
                                     ///to start the calculation. Electrons
                                     ///before that time are injected as a
                                     ///single 'blob'. This time is derived from
                                     ///the requirement that the blob has slid
                                     ///down to energies E<1GeV at t=Age.
  void PrepareAndRunNumericalSolver(vector<vector<double> > &particlespectrum,
                                    bool onlyprepare = false,
                                    bool dontinitialise =
                                        false);  ///< grid solver to calculate
                                                 ///the final particle spectrum
  void GetAxis(double min, double max, int steps, vector<double> &Axis,
               bool logarithmic);  ///< create an axis that attributes each bin
                                   ///with a real value (e.g. time axis, energy
                                   ///axis).
  void CreateGrid();  ///< create the 2D propagation grid out of 2 1D-vectors
  void SetInitialCondition(vector<double> EnergyAxis,
                           double startTime);  ///< set initial condition
                                               ///(a.k.a. set the first time
                                               ///spectrum in the grid)
  void ComputeGrid(vector<double> EnergyAxis,
                   double startTime, double Age, double minTimeBin =
                       1. * yr_to_sec);  ///< Iterate through the time-energy
                                         ///grid. Implemented is a piece-wise
                                         ///linear advection algorithm,
                                         ///optionally with the 'Superbee'-Slope
                                         ///limiter
  void FillEnergyTrajectoryLookup();
  double CalcSpecSemiAnalyticConstELossIntegrand(double *x, double *par);
  void CalcSpecSemiAnalyticConstELoss();
  void CalcSpecSemiAnalyticNoELoss();
  double MaxMod(double a,
                double b);  ///< Mod function required for slope delimiters
  double MinMod(double a,
                double b);  ///< Mod function required for slope delimiters
  double GetSuperBeeSlope(int i, double deltaX,int i0);  ///< slope for
                                                           ///superbee slope
                                                           ///limiter method
  double GetMinModSlope(int i, double deltaX,int i0);  ///< slope for the
                                                         ///minmod slope limiter
                                                         ///method
  double LaxWendroffSlope(int i, double deltaX, int i0);  ///< slope for the
                                                          // Lax-Wendroff slope limiter
  vector<vector<double> > ParticleSpectrum;  ///< vector containing the final
                                             ///source particle spectrum (i.e.
                                             ///at time=Age)
  vector<vector<double> > EscapeVector;      ///< vector containing escape delta
                                         ///functions for later diffusion (done
                                         ///by the 'Diffusion' class)
  double DetermineEmax(double tmin);  ///< maximum allowed particle energy over
                                      ///the course of the source's history
  bool DEBUG;                         ///< Debug flag
  double adLossCoeff;                 ///< coefficient for adiabatic losses
  double escapeTime;                  ///< particle escape time
  double escapeTimeConstant;  ///< particle escape time, Constantly set to a
                              ///single value
  double EnergyAxisLowerBoundary;  ///< might be superfluous (compare to
                                   ///EminConstant)
  double EnergyAxisUpperBoundary;  ///< might be superfluous (compare to
                                   ///EmaxConstant)
  bool QUIETMODE;  ///< Quietmode boolean: no output text at all if set to true
                   ///. Default is false.
  double Polynomial(double x, vector<double>);  ///< polynomial function, used
                                                ///for deriving systematic
                                                ///errors due to numerical
                                                ///wiggling.
  void CalculateConstants();  ///< speed hack where often used constants are
                              ///calculated ahead of the grid solving
  double PowerLawInjectionSpectrum(double e, double ecut,
                                   double emax);  ///< power law injection
                                                  ///spectrum. If eBreak and
                                                  ///SpectralIndex2 are
                                                  ///specified, returns a broken
                                                  ///power law.
  double InverseLossRate(double *x, double *par);
  void SetLookup(vector<vector<double> > v, string LookupType);  ///<Sets lookup for evolution
  void ExtendLookup(vector<vector<double> > v, string LookupType);
  double SemiAnalyticConstELossIntegrand(double T, void *par);
  double SourceSpectrumWrapper(double E, void *par);
  double integratorTolerance;
  double Integrate(fPointer f, double *x, double emin, double emax,
                   double tolerance, int kronrodrule);
  int gslmemory; ///< memory of the GSL workspace when integrating.
                 /// Default = 5000
  int kronrodrule;
  vector<double> Constants;  ///< Auxiliary vector for speed-up
  vector<gsl_spline *> splines;
  vector<gsl_interp_accel *> accs;
  vector<double *> vals;
  vector<vector<double> > vETrajectory;
  vector<vector<vector<double> > > vs;
  void CalculateEnergyTrajectory(double TExt = 0.);
  void DetermineLookupTimeBoundaries();
  void ComputeGridInTimeInterval(double T1, double T2, string type, int bins,
                                 bool ICRESET);
                                              ///< wrapper function to calculate
                                              /// the grid only in a specified
                                              /// time interval dT = T2-T2 (yrs)
  vector< vector<double> > GetEnergyLossRateVector(vector<double> epoints,
                                                   string type, double age, 
                                                        bool TIMESCALE); ///< Returns either
                                                                         /// energy loss rate or
                                                                         /// cooling time
  vector< vector<double> > GetInjectionSpectrumVector(vector<double> epoints,
                                                      double age, bool SED);
  void CalculateParticleSpectrum(string type = "electrons", int bins = 100,
                                 bool onlyprepare = false,
                                 bool dontinitialise = false); ///< fill the
                                                        ///lookup that holds
                                                    ///the particle spectrum.
 public:
  Particles();
  ~Particles();
  void SetMembers(double t);  ///< set the values for the class variables
                              ///BField,eElectronMax and Ecr at a given time t

  void CalculateProtonSpectrum(int bins = 100) {CalculateParticleSpectrum("protons",bins);}
  void CalculateElectronSpectrum(int bins = 100) {CalculateParticleSpectrum("electrons",bins);}
  void CalculateProtonSpectrumInTimeInterval(double T1, double T2,
                                             int bins = 100) {
    ComputeGridInTimeInterval(T1,T2,"protons",bins,false);}
  void CalculateElectronSpectrumInTimeInterval(double T1, double T2,
                                               int bins = 100, bool ICRESET=false) {
    ComputeGridInTimeInterval(T1,T2,"electrons",bins,ICRESET);}
  void CalculateElectronSpectrumInTimeInterval(int bins = 100) {CalculateParticleSpectrum("electrons",bins);}
  void SetType(string type);  ///< Set the particle type ("protons" or "electrons")
  double EnergyLossRate(double E);  ///< total energy loss rate of particles
  void SetIonization() {
	  IONIZATION = true;
  }  ///< Set the calculation of the ionization. This is because normally we consider fully ionized plasma
     /// if you really want to compute ionization losses, call this function
  void UnsetIonization() {
	  IONIZATION = false;
  }  ///< Unset accounting for ionization losses.
  void SetAge(double age) {
    Age = age;
    SetMembers(Age);
  }                                      ///< set source age
  double GetAge() const { return Age; }  ///< get source age
  double GetLuminosity(double age=0.) {
    if(age) SetMembers(age);
    CalculateConstants();
    return Lum;
  }  ///< get particle luminosity of the source
  double GetBField(double age=0.) {
    if(age) SetMembers(age);
    CalculateConstants();
    return BField;
  }  ///< get the BField at the acceleration site of the astrophysical particle
     ///accelerator.
  double GetAmbientDensity(double age=0.) {
    if(age) SetMembers(age);
    CalculateConstants();
    return N;
  }  ///< get the ambient density at the astrophysical particle accelerator.
  double GetRadius(double age=0.) {
    if(age) SetMembers(age);
/*    CalculateConstants();*/
    return R/pc_to_cm;
  }  ///< get the extension of the astrophysical particle accelerator (cm).
  double GetExpansionVelocity(double age=0.) {
    if(age) SetMembers(age);
    CalculateConstants();
    return V;
  }  ///< get the extension speed of the astrophysical particle accelerator
     ///(cm/s).
  double GetEmax(double age=0.) {
    if(age) SetMembers(age);
    CalculateConstants();
    return eMax; }  ///< get the max. particle energy (erg)
  double GetEscapeTime(double e, double age=0.) {
    if(age) SetMembers(age);
    CalculateConstants();
    return EscapeTime(e,age) / yr_to_sec;
  }  ///< the particle escape time scale (s)
  vector < vector<double> > GetEscapeTime(vector<double> epoints, double age=0.) {
      return GetQuantityVector(epoints,age,"escape_time");
  } ///< the particle escape time scale (s), returns spectral shape
  void ToggleDebugging() {
    DEBUG = true;
  }  ///< switch on Debugging/Testing mode
  void SetICLossLookup(vector<vector<double> > ICLOSSLOOKUP) {
    SetLookup(ICLOSSLOOKUP, "ICLoss");
  }  ///< set the lookup holding energy-dependent IC cooling rate {E-ICLossRate}

  void SetAmbientDensity(vector<vector<double> > NLOOKUP) {
    NConstant = NAN;
    SetLookup(NLOOKUP, "AmbientDensity");
  }  ///< Set ambient density evolution
  void SetBField(vector<vector<double> > BFIELDLOOKUP) {
    BConstant = NAN;
    SetLookup(BFIELDLOOKUP, "BField");
  }  ///< Set BField evolution

  void SetRadius(vector<vector<double> > RADIUSLOOKUP) {
    RConstant = NAN;
    SetLookup(RADIUSLOOKUP, "Radius");
  }  ///< Set radius evolution

  void SetExpansionVelocity(vector<vector<double> > VELOCITYLOOKUP) {
    VConstant = NAN;
    SetLookup(VELOCITYLOOKUP, "Speed");
  }  ///< Set expansion velocity evolution

  void SetLuminosity(vector<vector<double> > LUMLOOKUP) {
    LumConstant = NAN;
    SetLookup(LUMLOOKUP, "Luminosity");
  }  ///< DEPRECATED

  void ExtendICLossLookup(vector<vector<double> > ICLOSSLOOKUP) {
    ExtendLookup(ICLOSSLOOKUP, "ICLoss");
  }  ///< set the lookup holding energy-dependent IC cooling rate {E-ICLossRate}
  void ExtendLuminosityLookup(vector<vector<double> > LUMLOOKUP) {
    ExtendLookup(LUMLOOKUP, "Luminosity");
  }  ///< Set BField evolution
  void ExtendAmbientDensityLookup(vector<vector<double> > NLOOKUP) {
    ExtendLookup(NLOOKUP, "AmbientDensity");
  }  ///< Set ambient density evolution
  void ExtendBFieldLookup(vector<vector<double> > BFIELDLOOKUP) {
    ExtendLookup(BFIELDLOOKUP, "BField");
  }  ///< Set BField evolution
  void ExtendEmaxLookup(vector<vector<double> > EMAXLOOKUP) {
    ExtendLookup(EMAXLOOKUP, "Emax");
  }  ///< Set max. particle evolution
  void ExtendEscapeTimeLookup(vector<vector<double> > ESCTIMELOOKUP) {
    ExtendLookup(ESCTIMELOOKUP, "EscapeTime");
  }  ///< Set escape time evolution
  void ExtendRadiusLookup(vector<vector<double> > RLOOKUP) {
    ExtendLookup(RLOOKUP, "Radius");
  }  ///<
  void ExtendVelocityLookup(vector<vector<double> > VLOOKUP) {
    ExtendLookup(VLOOKUP, "Speed");
  }  ///<
/*  vector<vector<double> > GetICLossLookup() { return ICLossVector; }*/
  vector<vector<double> > GetLuminosityLookup() {
    return fUtils->VectorAxisPow10(LumVector,1);
  }
  vector<vector<double> > GetEmaxLookup() {
    return fUtils->VectorAxisPow10(eMaxVector,1);
  }
  vector<vector<double> > GetBFieldLookup() {return BVector; }
  vector<vector<double> > GetRadiusLookup() {return RVector; }
  vector<vector<double> > GetVelocityLookup() {return VVector; }
  double GetEnergyLossRate(double E, string type,double age=0.);
  double GetCoolingTimeScale(double E, string type="sum", double age=0.) {//< this gives the cooling time scale at E and time age in [years]
      return  E / GetEnergyLossRate(E,type) / yr_to_sec; }

  vector< vector<double> > GetEnergyLossRate(vector<double> epoints, string type="sum",
                                                            double age=0.) {
      return GetEnergyLossRateVector(epoints,type,age,false);
  }
  vector< vector<double> > GetCoolingTimeScale(vector<double> epoints, string type="sum",
                                               double age=0.) {
      return GetEnergyLossRateVector(epoints,type,age,true);
  }
  vector< vector<double> > GetInjectionSpectrum(vector<double> epoints, double age=0.) {
      return GetInjectionSpectrumVector(epoints,age,false);
  }
  vector< vector<double> > GetInjectionSED(vector<double> epoints, double age=0.) {
      return GetInjectionSpectrumVector(epoints,age,true);
  }
  vector< vector<double> > GetQuantityVector(vector<double> epoints,
                                                      double age, string mode);
  double GetEnergyBins() {
    return ebins;
  }  ///< get energy binning of the numerical solution
  vector<vector<double> > GetParticleSpectrum() {
    return ParticleSpectrum;
  }  ///< get vector that holds source particle spectrum
  vector<vector<double> > GetEscapeVector() {
    return EscapeVector;
  }  ///< get vector that holds the escape particles.
  void SetTmin(double TMIN) {
    TminExternal = TMIN;
  }  ///< Constanty set minimal time of injected particles
  void SetBField(double BEXT) {
    BVector.clear();
    BConstant = BEXT;
    SetMembers(TminInternal);
  }  ///< Constanty set B-Field value
  void SetAmbientDensity(double NCONSTANT) {
    NVector.clear();
    NConstant = NCONSTANT;
    SetMembers(TminInternal);
  }  ///< Constanty set value of ambient density
  void SetRadius(double r) {
    RVector.clear();
    RConstant = pc_to_cm*r;
    SetMembers(TminInternal);
  }  ///< Constanty set value of source extension (pc)
  void SetExpansionVelocity(double v) {
    VVector.clear();
    VConstant = v;
    SetMembers(TminInternal);
  }  ///< Constanty set value of source expansion speed (cm/s)
  void SetCutOffFactor(double CUTOFFFACTOR) {
    CutOffFactor = CUTOFFFACTOR;
  }  ///< set shape of the Electron spectrum cut-off
  double GetCutOffFactor() {
    return CutOffFactor;
  }  ///< get shape of the Electron spectrum cut-off
  void SetEnergyMarginFactor(double ENERGYMARGINFACTOR) {
    energyMarginFactor = ENERGYMARGINFACTOR;
  }  ///< set safety margin on the upper boundary of the particle spectrum
  void ToggleSharpEnergyCut() {
    sharpEnergyCut = true;
  }  ///< toggle (default: off) a sharp energy cut in the particle spectrum
     ///(position given by Emax)
  void UnToggleSharpEnergyCut() {
    sharpEnergyCut = false;
  }  ///< untoggle a sharp energy cut in the particle spectrum (position given
     ///by Emax)
  void SetTminInternal(double TMININTERNAL) {
    TminInternal = TMININTERNAL;
  }  ///< set minimal time from where to start the iteration (default: 1yr).
  double GetTminInternal() {
    return TminInternal;
  }  ///< set minimal time from where to start the iteration (default: 1yr).
  vector<vector<double> > GetEnergyTrajectoryVector() { 
    vector< vector<double> > v;
    v = fUtils-> VectorAxisPow10(vETrajectory,0);
    v = fUtils-> VectorAxisPow10(v,1);
    return v; }
  vector<vector<double> > GetParticleSED();
  vector<vector<double> > GetGrid(){return grid;}
  double GetParticleEnergyContent(double E1=0., double E2=0.);

  void SetLuminosity(double LUMConstant) {
    LumVector.clear();
    LumConstant = LUMConstant;
  }
  void SetEmax(vector<vector<double> > EMAXLOOKUP) {
    EmaxConstant = NAN;
    SetLookup(EMAXLOOKUP, "Emax");
  } 
  void SetCriticalMinEnergyForGridSolver(double eminint) {EminInternal=eminint;}
  void SetIntegratorMemory(string mode);
  void SetInterpolationMethod(string intermeth)
    {fUtils->SetInterpolationMethod(intermeth);}
  void SetCustomTimeEnergyLookup(vector< vector<double> > vCustom, int mode);
  void SetCustomEnergylookup(vector< vector<double> > vCustom,int mode);

  /* methods to set custom source injection spectra */
  void SetCustomInjectionSpectrum(vector< vector<double> > vSpectrum) {
    CustomInjectionSpectrumTimeEvolution = NULL;
    LumConstant = NAN;
    fUtils->Clear2DVector(LumVector);
    SetCustomEnergylookup(vSpectrum,0);} ///< custom injection spectrum

  void SetCustomInjectionSpectrumTimeEvolution(vector< vector<double> > vCustomSpectrum) {
    CustomInjectionSpectrum = NULL;
    SetCustomTimeEnergyLookup(vCustomSpectrum,0);}///< custom injection spectrum, changing over time

  void SetCustomInjectionSpectrumTimeEvolution(vector<double> t, vector<double> e, 
                                                vector< vector<double> > mesh) {
    CustomInjectionSpectrum = NULL;
    vector< vector<double> > vCustomSpectrum = fUtils->MeshgridToTwoDVector(t,e,mesh);
    SetCustomTimeEnergyLookup(vCustomSpectrum,0);}///< custom injection spectrum, changing over time, input numpy-style meshgrid
 
  /* methods to set particle escape */
  void SetConstantEscapeTime(double EscapeTime) {
    if(EscapeTime>0.) ESCTIME = true;
    escapeTimeLookupTdep = NULL;
    escapeTimeLookupEdep = NULL;
    EscapeTimeEnergyTimeEvolution = NULL;
    escapeTimeConstant = EscapeTime;
  }  ///< set the time scale of particle escape

  void SetTimeDependentEscapeTime(vector< vector<double> > vEsc) {
    ESCTIME = true;
    escapeTimeConstant = NAN;
    escapeTimeLookupEdep = NULL;
    EscapeTimeEnergyTimeEvolution = NULL;
    SetLookup(vEsc, "EscapeTimeTdep");} ///< Set escape time evolution

  void SetEnergyDependentEscapeTime(vector<vector<double> > ESCTIMELOOKUP) {
    ESCTIME = true;
    escapeTimeConstant = NAN;
    escapeTimeLookupTdep = NULL;
    EscapeTimeEnergyTimeEvolution = NULL;
    SetCustomEnergylookup(ESCTIMELOOKUP, 1);
  } ///< set energy-dependent escape time, constant 

  void SetTimeAndEnergyDependentEscapeTime(vector< vector<double> > vEsc) {
    ESCTIME = true;
    escapeTimeConstant = NAN;
    escapeTimeLookupEdep = NULL;
    escapeTimeLookupTdep = NULL;
    SetCustomTimeEnergyLookup(vEsc,1);}///< set energy-dependent escape time, dynamic

  void SetTimeAndEnergyDependentEscapeTime(vector<double> t, vector<double> e, 
                                                vector< vector<double> > mesh) {
    ESCTIME = true;
    escapeTimeConstant = NAN;
    escapeTimeLookupEdep = NULL;
    escapeTimeLookupTdep = NULL;
    vector< vector<double> > vEsc = fUtils->MeshgridToTwoDVector(t,e,mesh);
    SetCustomTimeEnergyLookup(vEsc,1);}///< set energy-dependent escape time, dynamic. input numpy-style meshgrid.

/*  Radiation *GetSSCEquilibrium(Radiation *fr,double t, double tolerance=1e-2);*/
  void SetSolverMethod(int method);
  void ToggleQuietMode() { QUIETMODE = QUIETMODE == true ? false : true; fRadiation->ToggleQuietMode();}  ///< toggle quiet mode on or off ( if on, no progress printout on the console)
  bool GetQuietMode() {return QUIETMODE;}
  void AddThermalTargetPhotons(double T, double edens, int steps=200);// wrapped from Radiation class. See Docu there.
  void ResetWithThermalTargetPhotons(int i, double T, double edens, int steps=200);// wrapped from Radiation class. See Docu there.
  void AddArbitraryTargetPhotons(vector<vector<double> > PhotonArray);// wrapped from Radiation class. See Docu there.
  void ResetWithArbitraryTargetPhotons(int i,vector<vector<double> > PhotonArray);
  void ImportTargetPhotonsFromFile(const char *phFile); // wrapped from Radiation class. See Docu there.
  void ResetWithTargetPhotonsFromFile(int i,const char *phFile);
  void AddSSCTargetPhotons(int steps=100); // wrapped from Radiation class. See Docu there. Radius in pc
  void ResetWithSSCTargetPhotons(int i, int steps=100);
  vector<vector<double> > GetICLossLookup(int i=-1) {
    return fRadiation->GetICLossLookup(i);
  };///< return TotalTargetPhotonVector
  double GetTargetPhotonFieldEnergyDensity(unsigned int i) {
    return fRadiation->GetTargetPhotonFieldEnergyDensity(i);
  }
  void CheckSanityOfTargetPhotonLookup(); // wrapped from Radiation class. See Docu there.
  vector<vector<double> > GetTargetPhotons(int i=-1); // wrapped from Radiation class. See Docu there.
  void ClearTargetPhotons(); // wrapped from Radiation class. See Docu there.

  void SetSynchrotronEmissionModel(int SYNCHMODEL) {
    fRadiation->SetSynchrotronEmissionModel(SYNCHMODEL);
  }  ///< externally switch the parameterisation of the synchrotron emission
     ///model. See SynchModel docu for options.
  vector<double> CalculateSSCEquilibrium(double tolerance = 5e-2, int bins = 100);
  void SetFastIteration() {FASTMODE = true;}
  void SetPreciseIteration() {FASTMODE = false;}
  int GetTargetFieldCount(){ return fRadiation->GetTargetFieldCount();}
  /*********************************************************/
  /* DEPRECATED FUNCTIONS KEPT FOR BACKWARDS COMPATIBILITY */
  void SetEmin(double EMIN, bool ONLYFORNORMALISATION =
                                false) {  ///< Constanty set minimal energy of
                                          ///particle spectrum
    cout<< "SetEmin: This way of specifying the minimum energy of the source spectrum is DEPRECATED. "
           "This is done now by specifying a 2D-spectrum as input. See the documentation "
           "website on how to do it now!"<<endl;
    EminConstant = EMIN;
    if (ONLYFORNORMALISATION) EminConstantForNormalisationOnly = true;
  } ///< DEPRECATED
  void SetEmax(double EMAX) {
    cout<< "SetEmax: This way of specifying the maximum energy of the source spectrum is DEPRECATED. "
           "This is done now by specifying a 2D-spectrum as input. See the documentation "
           "website on how to do it now!"<<endl;
    eMaxVector.clear();
    EmaxConstant = EMAX;
  }  ///< DEPRECATED
  void SetSpectralIndex(double spectralindex) {
    cout<< "SetSpectralIndex: This way of specifying the index of the source spectrum is DEPRECATED. "
           "This is done now by specifying a 2D-spectrum as input. See the documentation "
           "website on how to do it now!"<<endl;
    SpectralIndex = spectralindex;
    if (SpectralIndex == 2.) SpectralIndex += 0.0000001;
  }  ///< DEPRECATED
  void SetLowSpectralIndex(double spectralindex2) {
    cout<< "SetLowSpectralIndex: This way of specifying the lower index of the source spectrum is DEPRECATED. "
           "This is done now by specifying a 2D-spectrum as input. See the documentation "
           "website on how to do it now!"<<endl;
    SpectralIndex2 = spectralindex2;
    if (SpectralIndex2 == 2.) SpectralIndex2 += 0.0000001;
  }  ///< DEPRECATED
  void SetBreakEnergy(double ebreak) {
    cout<< "SetBreakEnergy: This way of specifying the break energy of the source spectrum is DEPRECATED. "
           "This is done now by specifying a 2D-spectrum as input. See the documentation "
           "website on how to do it now!"<<endl;
    eBreak = ebreak;
  }  ///< DEPRECATED
  void SetSourceExpansionSpeed(double v) {
    cout << "SetSourceExpansionSpeed: DEPRECATED! USE Particles::SetExpansionVelocity(double v) instead!"<<endl;
    SetExpansionVelocity(v);
  }  ///< DEPRECATED
  double GetSpeed() {
    cout<<"GetSpeed: DEPRECATED! Use Particles::GetExpansionVelocity() instead! " <<endl; 
    return V;
  }  ///< DEPRECATED
  void SetEmaxLookup(vector<vector<double> > EMAXLOOKUP) {
    cout<< "SetEmaxLookup: This way of specifying the maximum energy is DEPRECATED. "
           "this is done now by a 3D-spectrum lookup. See the documentation "
           "website on how to do it now!"<<endl;
    EmaxConstant = NAN;
    SetLookup(EMAXLOOKUP, "Emax");
  }  ///< DEPRECATED
  void SetAmbientDensityLookup(vector<vector<double> > NLOOKUP) {
    cout<< "SetAmbientDensityLookup: DEPRECATED! USE Particles::SetAmbientDensity(<vector<vector<double>>v) instead!"<<endl;
    SetAmbientDensity(NLOOKUP);
  }  ///< DEPRECATED
  void SetLuminosityLookup(vector<vector<double> > LUMLOOKUP) {
    LumConstant = NAN;
    SetLookup(LUMLOOKUP, "Luminosity");
  }  ///< DEPRECATED
  void SetBFieldLookup(vector<vector<double> > BFIELDLOOKUP) {
    cout<< "SetBFieldLookup: DEPRECATED! USE Particles::SetBField(<vector<vector<double>>v) instead!"<<endl;
    SetBField(BFIELDLOOKUP);
  }  ///< DEPRECATED!
  void SetEscapeTimeLookup(vector<vector<double> > ESCTIMELOOKUP) {
    cout<< "SetEscapeTimeLookup: DEPRECATED! USE Particles::SetEscapeTime(<vector<vector<double>>v) instead!"<<endl;
    ESCTIME = true;
    SetTimeDependentEscapeTime(ESCTIMELOOKUP);
  }  ///< DEPRECATED
  void SetRadiusLookup(vector<vector<double> > RADIUSLOOKUP) {
    cout<< "SetRadiusLookup: DEPRECATED! USE Particles::SetRadius(<vector<vector<double>>v) instead!"<<endl;
    SetRadius(RADIUSLOOKUP);
  }  ///< DEPRECATED
  void SetVelocityLookup(vector<vector<double> > VELOCITYLOOKUP) {
    cout<< "SetVelocityLookup: DEPRECATED! USE Particles::SetExpansionVelocity(<vector<vector<double>>v) instead!"<<endl;
    SetExpansionVelocity(VELOCITYLOOKUP);
  }  ///< DEPRECATED
  void SetEnergyBins(double EBINS) {
    ebins = EBINS;
  }   ///< DEPRECATED
};
#endif
