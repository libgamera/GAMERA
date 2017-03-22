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

  static void initialise(fPointer funcPtr, Particles *partPtr) {
    _funcPtr = funcPtr;
    _partPtr = partPtr;
  }
  static double evaluate(double x, void* params) {
    return (_partPtr->*_funcPtr)(x, params);
  }
  static fPointer _funcPtr;
  static Particles *_partPtr;

struct timespec time0, time1, time2, time3;
 private:
  Utils *fUtils;
  int Type;  ///< integer indicating particle type. supported: 0-electrons and
             ///1-protons
  int METHOD; /// calculation method. This is determined automatically, but can
              /// be force set via SetSolverMethod(). However, if an inappropriate
              /// method is chosen manually, the result might be wrong or 
              /// the program might crash.
              /// 0: grid solver,
              /// 1: semi-analytical method in the presence of constant losses
              /// 2: no losses (simply adds up according to luminosity)
                
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
  double TminConstant;  ///< Constantly set lower time boundary at which the
                        ///particle iteration is started
  double BConstant;     ///< Constantly set constant B-Field
  double NConstant;     ///< Constantly set constant ambient density
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
  double SACELI_Told;  ///< Helper quantity that saves the last time step in the
                       ///integration in CalcSpecSemiAnalyticConstELoss

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
      escapeTimeVector, RVector, VVector,tempspec,
      CustomInjectionSpectrumTimeEvolutionVector,
      EscapeTimeEnergyTimeEvolutionVector;
  gsl_spline *ICLossLookup, *LumLookup, *NLookup, *BFieldLookup, *eMaxLookup,
      *escapeTimeLookup, *RLookup, *VLookup, *energyTrajectory,
      *energyTrajectoryInverse,*CustomInjectionSpectrum;
  interp2d_spline *CustomInjectionSpectrumTimeEvolution,
      *EscapeTimeEnergyTimeEvolution;
  gsl_interp_accel *accIC, *accLum, *accN, *accBField, *acceMax, *accescapeTime,
      *accR, *accV, *accTr, *accTrInv,*accCustInj,
      *taccsp,*eaccsp,*taccesc,*eaccesc;
  double SourceSpectrum(double e);  ///< particle injection spectrum
  double EscapeTime(double e);  ///< particle escape time spectrum
  vector<double> timeAxis;          ///< time axis for numerical integrator
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
  void SetInitialCondition(vector<vector<double> > &Grid,
                           vector<double> EnergyAxis,
                           double startTime);  ///< set initial condition
                                               ///(a.k.a. set the first time
                                               ///spectrum in the grid)
  void ComputeGrid(vector<vector<double> > &Grid, vector<double> EnergyAxis,
                   double startTime, double Age, vector<double> &TimeAxis,
                   double minTimeBin =
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
  double GetSuperBeeSlope(int i, double deltaX,
                          vector<vector<double> > *Grid);  ///< slope for
                                                           ///superbee slope
                                                           ///limiter method
  double GetMinModSlope(int i, double deltaX,
                        vector<vector<double> > *Grid);  ///< slope for the
                                                         ///minmod slope limiter
                                                         ///method
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
  void SetLookup(vector<vector<double> > v, string LookupType);  ///<
  void ExtendLookup(vector<vector<double> > v, string LookupType);
  double SemiAnalyticConstELossIntegrand(double T, void *par);
  double SourceSpectrumWrapper(double E, void *par);
  double integratorTolerance;
  double Integrate(fPointer f, double *x, double emin, double emax,
                   double tolerance, int kronrodrule);
  int gslmemory; ///< memory of the GSL workspace when integrating.
                 /// Default = 5000
  int kronrodrule;
  vector<double> Constants;
  vector<gsl_spline *> splines;
  vector<gsl_interp_accel *> accs;
  vector<double *> vals;
  vector<vector<double> > vETrajectory;
  vector<vector<vector<double> > > vs;
  void CalculateEnergyTrajectory(double TExt = 0.);
  void DetermineLookupTimeBoundaries();
  void SetSolverMethod(int method);
  void ComputeGridInTimeInterval(double T1, double T2, string type, int bins);  
                                              ///< wrapper function to calculate
                                              ///the grid only in a specified
                                              ///time interval dT = T2-T2 (yrs)

 public:
  Particles();
  ~Particles();
  void SetMembers(double t);  ///< set the values for the class variables
                              ///BField,eElectronMax and Ecr at a given time t
  void CalculateParticleSpectrum(string type = "electrons", int bins = 100,
                                 bool onlyprepare = false,
                                 bool dontinitialise = false); ///< fill the 
                                                        ///lookup that holds
                                                    ///the particle spectrum.
  void CalculateProtonSpectrum(int bins = 100) {CalculateParticleSpectrum("protons",bins);}
  void CalculateElectronSpectrum(int bins = 100) {CalculateParticleSpectrum("electrons",bins);}
  void CalculateProtonSpectrumInTimeInterval(double T1, double T2, 
                                             int bins = 100) {
    ComputeGridInTimeInterval(T1,T2,"protons",bins);}
  void CalculateElectronSpectrumInTimeInterval(double T1, double T2, 
                                               int bins = 100) {
    ComputeGridInTimeInterval(T1,T2,"electrons",bins);}
  void CalculateElectronSpectrumInTimeInterval(int bins = 100) {CalculateParticleSpectrum("electrons",bins);}
  void SetType(string type);
  double EnergyLossRate(double E);  ///< total energy loss rate of particles

  void SetAge(double age) {
    Age = age;
    SetMembers(Age);
  }                                      ///< set source age
  double GetAge() const { return Age; }  ///< get source age
  double GetLuminosity() {
    return Lum;
  }  ///< get particle luminosity of the source
  double GetBField() {
    return BField;
  }  ///< get the BField at the acceleration site of the astrophysical particle
     ///accelerator.
  double GetAmbientDensity() {
    return N;
  }  ///< get the ambient density at the astrophysical particle accelerator.
  double GetRadius() {
    return R;
  }  ///< get the extension of the astrophysical particle accelerator (cm).
  double GetSpeed() {
    return V;
  }  ///< get the extension speed of the astrophysical particle accelerator
     ///(cm/s).
  double GetEmax() { return eMax; }  ///< get the max. particle energy (erg)
  double GetEscapeTime() {
    return escapeTime;
  }  ///< the particle escape time scale (s)
  void SetConstantEscapeTime(double EscapeTime) {
    escapeTimeConstant = EscapeTime;
  }  ///< set the time scale of particle escape
  void SetSpectralIndex(double spectralindex) {
    SpectralIndex = spectralindex;
    if (SpectralIndex == 2.) SpectralIndex += 0.0000001;
  }  ///< set spectral index of the particle injection spectrum
  void SetLowSpectralIndex(double spectralindex2) {
    SpectralIndex2 = spectralindex2;
    if (SpectralIndex2 == 2.) SpectralIndex2 += 0.0000001;
  }  ///< set low-energy spectral index of the broken power-law particle
     ///injection spectrum. If this is not set, a single power-law injection
     ///spectrum will be assumed.
  void SetBreakEnergy(double ebreak) {
    eBreak = ebreak;
  }  ///< set break energy of broken powerlaw. If this is not set, a single
     ///power-law injection spectrum will be assumed.
  void ToggleDebugging() {
    DEBUG = true;
  }  ///< switch on Debugging/Testing mode
  void SetICLossLookup(vector<vector<double> > ICLOSSLOOKUP) {
    SetLookup(ICLOSSLOOKUP, "ICLoss");
  }  ///< set the lookup holding energy-dependent IC cooling rate {E-ICLossRate}
  void SetLuminosityLookup(vector<vector<double> > LUMLOOKUP) {
    SetLookup(LUMLOOKUP, "Luminosity");
  }  ///< Set BField evolution
  void SetAmbientDensityLookup(vector<vector<double> > NLOOKUP) {
    SetLookup(NLOOKUP, "AmbientDensity");
  }  ///< Set ambient density evolution
  void SetBFieldLookup(vector<vector<double> > BFIELDLOOKUP) {
    SetLookup(BFIELDLOOKUP, "BField");
  }  ///< Set BField evolution
  void SetEmaxLookup(vector<vector<double> > EMAXLOOKUP) {
    SetLookup(EMAXLOOKUP, "Emax");
  }  ///< Set max. particle evolution
  void SetEscapeTimeLookup(vector<vector<double> > ESCTIMELOOKUP) {
    SetLookup(ESCTIMELOOKUP, "EscapeTime");
  }  ///< Set escape time evolution
  void SetRadiusLookup(vector<vector<double> > RADIUSLOOKUP) {
    SetLookup(RADIUSLOOKUP, "Radius");
  }  ///< Set radius evolution
  void SetVelocityLookup(vector<vector<double> > VELOCITYLOOKUP) {
    SetLookup(VELOCITYLOOKUP, "Speed");
  }  ///< Set expansion velocity evolution
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
  vector<vector<double> > GetICLossLookup() { return ICLossVector; }
  vector<vector<double> > GetLuminosityLookup() { 
    return fUtils->VectorAxisPow10(LumVector,1); 
  }
  vector<vector<double> > GetEmaxLookup() { 
    return fUtils->VectorAxisPow10(eMaxVector,1); 
  }
  vector<vector<double> > GetBFieldLookup() {return BVector; }
  vector<vector<double> > GetRadiusLookup() {return RVector; }
  vector<vector<double> > GetVelocityLookup() {return VVector; }
  double GetEnergyLossRate(double E) { return EnergyLossRate(E); }
  void SetEnergyBins(double EBINS) {
    ebins = EBINS;
  }  ///< set energy binning of the numerical solution
  double GetEnergyBins() {
    return ebins;
  }  ///< get energy binning of the numerical solution
  vector<vector<double> > GetParticleSpectrum() {
    return ParticleSpectrum;
  }  ///< get vector that holds source particle spectrum
  vector<vector<double> > GetEscapeVector() {
    return EscapeVector;
  }  ///< get vector that holds the escape particles.
  void SetEmin(double EMIN, bool ONLYFORNORMALISATION =
                                false) {  ///< Constanty set minimal energy of
                                          ///particle spectrum
    EminConstant = EMIN;
    if (ONLYFORNORMALISATION) EminConstantForNormalisationOnly = true;
  }
  void SetTmin(double TMIN) {
    TminConstant = TMIN;
  }  ///< Constanty set minimal time of injected particles
  void SetEmax(double EMAX) {
    eMaxVector.clear();
    EmaxConstant = EMAX;
  }  ///< Constanty set maximal energy of particle spectrum
  void SetBField(double BEXT) {
    BVector.clear();
    BConstant = BEXT;
  }  ///< Constanty set B-Field value
  void SetAmbientDensity(double NCONSTANT) {
    NVector.clear();
    NConstant = NCONSTANT;
  }  ///< Constanty set value of ambient density
  void SetLuminosity(double LUMConstant) {
    LumVector.clear();
    LumConstant = LUMConstant;
  }  ///< Constantly set value of source luminosity
  void SetRadius(double r) {
    RVector.clear();
    RConstant = pc_to_cm*r;
  }  ///< Constanty set value of source extension (pc)
  void SetSourceExpansionSpeed(double v) {
    VVector.clear();
    VConstant = v;
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
  void ToggleQuietMode() {
    QUIETMODE = true;
  }  ///< untoggle quiet mode (no progress printout on the console)
  void SetEnergyAxisLowerBoundary(double BOUND) {
    EnergyAxisLowerBoundary = BOUND;
  }  ///< might be superflous
  void SetEnergyAxisUpperBoundary(double BOUND) {
    EnergyAxisUpperBoundary = BOUND;
  }  ///< might be superflous
  void SetTminInternal(double TMININTERNAL) {
    TminInternal = TMININTERNAL;
  }  ///< set minimal time from where to start the iteration (default: 1yr).
  double GetTminInternal() {
    return TminInternal;
  }  ///< set minimal time from where to start the iteration (default: 1yr).
  vector<vector<double> > GetEnergyTrajectoryVector() { return vETrajectory; }
  vector<vector<double> > GetParticleSED();
  double GetParticleEnergyContent(double E1, double E2);
  void SetCriticalMinEnergyForGridSolver(double eminint) {EminInternal=eminint;}
  void SetIntegratorMemory(string mode);
  void SetInterpolationMethod(string intermeth)
    {fUtils->SetInterpolationMethod(intermeth);}
  void SetCustomTimeEnergyLookup(vector< vector<double> > vCustom, int mode);
  void SetCustomEnergylookup(vector< vector<double> > vCustom,int mode);
  void SetCustomInjectionSpectrumTimeEvolution(vector< vector<double> > vCustomSpectrum) {
    SetCustomTimeEnergyLookup(vCustomSpectrum,0);}
  void SetCustomInjectionSpectrum(vector< vector<double> > vSpectrum) {
    SetCustomEnergylookup(vSpectrum,0);}
  void SetTimeAndEnergyDependentEscapeTime(vector< vector<double> > vEsc) {
    SetCustomTimeEnergyLookup(vEsc,1);}
  void SetEnergyDependentEscapeTime(vector< vector<double> > vEsc) {
    SetCustomEnergylookup(vEsc,1);}
  Radiation *GetSSCEquilibrium(Radiation *fr,double t, double tolerance=1e-2);
};
#endif
