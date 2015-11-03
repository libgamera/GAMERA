#ifndef _RADIATION_
#define _RADIATION_

#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

/* TeV->erg */
#define TeV_to_erg 1.602
/* GeV->erg */
#define GeV_to_erg 1.602e-3
/* Thomson cross section */
#define sigma_T 6.6524e-25
/* electron mass in erg */
#define m_e 8.187e-7
/* boltzmann constant (erg/K) */
#define kb 1.380658e-16
/* proton mass in erg */
#define m_p 1.50310854e-3
/* pi0 mass in erg */
#define m_pi 2.1622194e-4
/* parsec to cm */
#define pc_to_cm 3.0857e18
/* well... pi! */
#define pi 3.1416
/* year in seconds */
#define yr_to_sec 3.15576e7
/* solar mass */
#define mSol 1.9891e33
/* speed of light cm/s */
#define c_speed 29979245800.
/* elementary charge */
#define el_charge 4.80320427e-10
/* classical electron radius (cm) */
#define e_radius 2.8179e-13
/* planck's constant */
#define hp 6.62606896e-27
/* fine structure constant */
#define fineStructConst 7.2974e-3

using namespace std;

typedef double (*fGSLPointer)(double, void *);

template <typename F>
class GSLfuncRad : public gsl_function {
 public:
  GSLfuncRad(const F &func) : _func(func) {
    function = &GSLfuncRad::Call;
    params = this;
  }

 private:
  const F &_func;
  static double Call(double x, void *params) {
    return static_cast<GSLfuncRad *>(params)->_func(x);
  }
};

/**  @file Radiation.C
 *   @author Joachim Hahn
 *   @short Class that calculates broad band emission from particle spectra.
 */

class Radiation {
  typedef double (Radiation::*fPointer)(double, void *);

 private:
  void CalculateLuminosityAndFlux(string mechanism, double e, double &l, 
                                  double &f); 
  double DifferentialEmissionComponent(double e, void *par);
  double GreyBody(double ephoton, double temp, double edens);
  double ICEmissivityRadFieldIntegrated(double x, void *par);
  double ICEmissivity(double x, void *par);
  double K(double nu, double x);             ///< modified Bessel function
  double K_53(double x, void *par);  ///< modified Bessel function of order 5/3
  double SynchEmissivity(double x, void *par);  ///< Synchrotron emission from a
                                                ///number of electrons with
                                                ///energy eElectron at frequency
                                                ///nu
  double SynchEmissivityExplicit(double e,
                                 void *par);  ///< Synchrotron emission from a
                                              ///number of electrons with energy
                                              ///eElectron at frequency nu,
                                              ///following Blumenthal&Gould
  double BremsEmissivity(double x, void *par);  ///< Bremsstrahlung emissivity,
                                                ///e-p and e-e, follwing From
                                                ///Baring 1999, ApJ, 513,
                                                ///311-338, Eq. (27)
  double A(double g,
           double e);  ///< Equation (A4) from Baring 1999, ApJ, 513, 311-338
  double sigma1(double g, double e);   ///< Equation (A2) from Baring 1999, ApJ,
                                       ///513, 311-338
  double sigma2(double g, double e);   ///< Equation (A3) from Baring 1999, ApJ,
                                       ///513, 311-338
  double sigmaNR(double g, double e);  ///< Equation (A5) from Baring 1999, ApJ,
                                       ///513, 311-338
  double Fbr(double x, double g);  ///< Equations (A6,A7) from Baring 1999, ApJ,
                                   ///513, 311-338
  double ParticleSpectrum(double E);  ///< particle spectrum emitting the
                                      ///emission. This is basically reading in
                                      ///a 2D vector(i.e. 'ParticleLookup')
                                      ///holding E-N(particles).
  vector<vector<double> > ParticleVector;  ///< 2D vector holding the particle
                                           ///spectrum that emits the radiation.
                                           ///Format: E(erg)-N(number). This can
                                           ///be either electrons or protons,
                                           ///and is set to either species
                                           ///automatically depending on the
                                           ///radiation mechanism that needs to
                                           ///be calculated.
  vector<vector<double> > ElectronVector;  ///< format of 'ParticleLookup',
                                           ///holding the electron spectrum
  vector<vector<double> > ProtonVector;    ///< format of 'ParticleLookup',
                                           ///holding the proton spectrum
  vector<vector<double> > TargetPhotonVector;  ///< 2D format
                                               ///E(erg)-dN/dE(differential
                                               ///number) vector holding the
                                               ///total target field for IC
                                               ///scattering. This is the sum of
                                               ///individual components
                                               ///specified by public functions
  vector<vector<double> > TargetPhotonVectorOld;  ///< 2D format
                                                  ///E(erg)-dN/dE(differential
                                                  ///number) vector holding the
                                                  ///total target field for IC
                                                  ///scattering. This is the sum
                                                  ///of individual components
                                                  ///specified by public
                                                  ///functions
  vector<vector<double> > SSCTargetPhotons;  ///< 2D-vector object holding the
                                             ///SSC target photon spectrum  {
                                             ///E(erg) - Edens(erg cm^-3) }
  double n;  ///< ambient density for Bremsstrahlung and inelastic p-p emission
             ///mechanisms (cm^-3)
  double fintbrems;   ///< integrated flux above specified energy E (in
                      ///CalculateIntegratedGammaEmission) due to Bremsstrahlung
  double lintbrems;   ///< integrated luminosity above specified energy E (in
                      ///CalculateIntegratedGammaEmission) due to Bremsstrahlung
  double fintpp;      ///< integrated flux above specified energy E (in
                      ///CalculateIntegratedGammaEmission) due to inelastic p-p
                      ///collision
  double lintpp;      ///< integrated luminosity above specified energy E (in
                      ///CalculateIntegratedGammaEmission) due to inelastic p-p
                      ///collision
  double fintic;      ///< integrated flux above specified energy E (in
                      ///CalculateIntegratedGammaEmission) due to inelastic IC
                      ///radiation
  double lintic;      ///< integrated luminosity above specified energy E (in
                      ///CalculateIntegratedGammaEmission) due to inelastic IC
                      ///radiation
  double fdiffbrems;  ///< differential flux at specified energy E (in
                      ///CalculateDifferentialGammaEmission) due to
                      ///Bremsstrahlung
  double ldiffbrems;  ///< differential luminosity at specified energy E (in
                      ///CalculateDifferentialGammaEmission) due to
                      ///Bremsstrahlung
  double fdiffpp;     ///< differential flux at specified energy E (in
                   ///CalculateDifferentialGammaEmission) due to inelastic p-p
                   ///collision
  double ldiffpp;  ///< differential luminosity at specified energy E (in
                   ///CalculateDifferentialGammaEmission) due to inelastic p-p
                   ///collision
  double fdiffic;  ///< differential flux at specified energy E (in
                   ///CalculateDifferentialGammaEmission) due to IC emission
  double ldiffic;  ///< differential luminosity at specified energy E (in
                   ///CalculateDifferentialGammaEmission) due to IC emission
  double fdiffsynch;  ///< differential flux at specified energy E (in
                      ///CalculateDifferentialGammaEmission) Synchrotron
                      ///emission
  double ldiffsynch;  ///< differential luminosity at specified energy E (in
                      ///CalculateDifferentialGammaEmission) Synchrotron
                      ///emission
  double distance;    ///< source distance (cm)
  double targetphotonenergymin;  ///< lower energy boundary for the target
                                 ///photon spectrum (erg, dynamically
                                 ///determined)
  double targetphotonenergymax;  ///< upper energy boundary for the target
                                 ///photon spectrum (erg, dynamically
                                 ///determined)
  void FillTargetPhotonVectorAndGraph(int steps);  ///< function that adds
                                                   ///individual target photon
                                                   ///spectrum to
                                                   ///"TotalTargetPhotonGraph"
                                                   ///and
                                                   ///"TotalTargetPhotonVector"
  int TargetPhotonSteps;  ///< binning of the target photon spectrum
  bool DEBUG;             ///< debugging boolean
  double integratorTolerance;///
  string radiationMechanism;  ///< string query the radiation mechanism that is
                              ///to be calculated
  double hadronicAmpFactor;   ///< amplification factor of the p-p emission (due
                              ///to heavier species in the ISM)
  double BField;              ///< BField value (G)
  bool INTEGRATEOVERGAMMAS;   ///< boolean that switches between calculation IC
                              ///emission loss rate and emission
  bool QUIETMODE;  ///< boolean that toggles quiet output mode if set to true
  vector<vector<double> > diffSpec;  ///< vector holding all the individual
                                     ///differential spectra in erg -
                                     ///erg^-1s^-1cm^-2
  int PiModel;  ///< indicates which parameterisation to use in the Kafexhiu pi0
                ///model. 0 - Geant4, 1 - Pythia8, 2 - SIBYLL2.1, 3 - QGSJET-I.
                ///DEFAULT = 1
  int SynchModel;  ///< indicates which synchrotron emissivity model should be
                   ///used. 0 - random B-field (Gisellini 1988) 1 - regular
                   ///B-Field, with perpendicularly spiraling electron around
                   ///it. DEFAULT = 0
  void Clear2DVector(vector< vector<double> > &v);
  double PPEmissivity(double x, void *par);
  double InelasticPPXSectionKaf(double Tp);
  double InclusivePPXSection(double Tp);
  double DiffPPXSection(double Tp, double Eg);
  double MeanMultiplicity(double Tp);
  double Amax(double Tp);
  double F(double Tp, double Eg);
  double SigmaOnePi(double Tp);
  double SigmaTwoPi(double Tp);
  double GetMaximumGammaEnergy(double Tp);
  double GetMinimumGammaEnergy(double Tp);
  double Epilabmax(double Tp);
  double NuclearEnhancementFactor(double Tp);
  void GetABGParams(double Tp, double &alpha, double &beta, double &gamma,
                    double &lambda);
  void GetAParams(double Tp, double &a1, double &a2, double &a3, double &a4,
                  double &a5);
  void GetBParams(double Tp, double &b1, double &b2, double &b3);
  static double GSLFunctionWrapper(double x, void *params);
  /*    Utils *fUtils;*/
  double Integrate(fPointer f, double *x, double emin, double emax,
                   double tolerance);  ///< generic integration routine
  gsl_interp_accel *acc;  ///< gsl accelerator object for interpolation
  gsl_spline *ElectronLookup, *ProtonLookup, *TargetPhotonLookup,
      *TargetPhotonLookupEdens;
  void SetParticles(vector<vector<double> > PARTICLES, int type);
  void AddToTargetPhotonVector(gsl_spline *Spl, double logEminSpl,
                               double logEmaxSpl, int stepsSpl);
  void SetTargetPhotonVectorLookup();
  const gsl_interp_type *interp_type;
  vector<vector<double> > ICLossLookup;
  double TargetPhotonEdens;
  vector<vector<double> > GetParticleSED(string type);

 public:
  Radiation();                                       ///< standard constructor
  ~Radiation();                                      ///< standard destructor
  void SetProtons(vector<vector<double> > PROTONS);  ///< set the proton
                                                     ///spectrum (e.g.
                                                     ///calculated in the
                                                     ///"Particles" class, but
                                                     ///also arbitrary spectra).
                                                     ///Input format: 2D vector,
                                                     ///E(erg)-N(number)
  void SetElectrons(vector<vector<double> > ELECTRONS);  ///< set the electron
                                                         ///spectrum (e.g.
                                                         ///calculated in the
                                                         ///"Particles" class,
                                                         ///but also arbitrary
                                                         ///spectra). Input
                                                         ///format: 2D vector,
                                                         ///E(erg)-N(number)
  vector<vector<double> > GetProtonVector() {
    return ProtonVector;
  }  ///< return proton spectrum return format: 2D vector
  vector<vector<double> > GetElectronVector() {
    return ElectronVector;
  }  ///< return electron spectrum return format: 2D vector
  void SetAmbientDensity(double N) {
    n = N;
  }  ///< set ambient number density (cm^-3)
  double GetAmbientDensity() {
    return n;
  }  ///< get ambient number density (cm^-3)
  void CalculateIntegralGammaEmission(double e, int particletype);
  void CalculateDifferentialGammaEmission(double e, int particletype);
  void CalculateDifferentialPhotonSpectrum(int steps = 100, double emin = 0.,
                                           double emax = 0.);
  vector<vector<double> > ReturnDifferentialPhotonSpectrum(int i,
                                                           double emin = 0.,
                                                           double emax = 0.);
  vector<vector<double> > ReturnSED(int i, double emin = 0.,
                                    double emax = 0.);  ///< returns SED for
                                                        ///emission component i
                                                        ///as 2D vector
  void SetBField(double BFIELD) {
    BField = BFIELD;
  }  ///< set the source B-Field (G)
  void SetDistance(double d) {
    distance = d;
  }  ///< set the distance to the source (cm)
  double GetDifferentialICFlux() { return fdiffic; }        ///< get FdiffIC
  double GetDifferentialSynchFlux() { return fdiffsynch; }  ///< get fdiffsynch
  double GetDifferentialBremsFlux() { return fdiffbrems; }  ///< get fdiffbrems
  double GetDifferentialPPFlux() { return fdiffpp; }        ///< get fdiffpp
  double GetIntegralICFlux() { return fintic; }             ///< get fintic
  double GetIntegralBremsFlux() { return fintbrems; }       ///< get fintbrems
  double GetIntegralPPFlux() { return fintpp; }             ///< get fintpp
  double GetICLuminosity() { return lintic; }               ///< get lintic
  double GetBremsLuminosity() { return lintbrems; }         ///< get lintbrems
  double GetPPLuminosity() { return lintpp; }               ///< get lintpp
  void CreateICLossLookup(int bins = 200);  ///< creates a 2D vector holding the
                                            ///energy-dependent energy loss rate
                                            ///due to IC cooling. Useful to
                                            ///apply in spectral iterations in
                                            ///conjunction in the "Particles"
                                            ///class. Format: E(erg) -
                                            ///-1.*LossrateIC (erg/s)
  void AddThermalTargetPhotons(double T, double energydens,
                               int steps = 200);  ///< add a thermal target
                                                  ///radiation field component
                                                  ///to the total local
                                                  ///radiation field
  void AddArbitraryTargetPhotons(
      vector<vector<double> > PhotonArray);  ///< add an arbitrary target photon
                                             ///field. input is a 2D-vector of
                                             ///format E(erg)
                                             ///photon_density(erg^-1cm^-3)
  void ImportTargetPhotonsFromFile(
      const char *phFile);  ///< add a target photon density from an ASCII file
                            ///of format E(eV) photon_density(eV^-1cm^-3) !!!
                            ///here eV, as this is what is typically used in the
                            ///literature !!!
  void RemoveLastICTargetPhotonComponent();  ///< remove the latest component in
                                             ///TotalTargetPhotonVector and
                                             ///recompute the total target
                                             ///photon spectrum
  void AddSSCTargetPhotons(double R, int steps = 100);  ///< add target photons
                                                  ///resulting from Synchrotron
                                                  ///radiation due to current
                                                  ///electron spectrum and
                                                  ///B-Field.
  vector<vector<double> > GetTotalTargetPhotonVector() {
    return TargetPhotonVector;
  }              ///< return TotalTargetPhotonVector
  void Reset();  ///< reset ParticleLookup, Electrons, Protons, fintbrems,
                 ///lintbrems, fintpp, lintpp, fintic, lintic
  void ToggleQuietMode() {
    QUIETMODE = true;
  }  ///< enable quiet mode (very little cout output)
  vector<vector<double> > GetICLossLookup() {
    return ICLossLookup;
  }  ///< return TotalTargetPhotonVector
  void SetHadronicAmplificationFactor(double HADAMPFAC) {
    hadronicAmpFactor = HADAMPFAC;
  }  ///< set hadronicAmpFactor
  vector<vector<double> > GetTotalSpectrum(double emin = 0., double emax = 0.) {
    return ReturnDifferentialPhotonSpectrum(1, emin, emax);
  }  ///< return total spectrum
  vector<vector<double> > GetPPSpectrum(double emin = 0., double emax = 0.) {
    return ReturnDifferentialPhotonSpectrum(2, emin, emax);
  }  ///< return pi0 decay spectrum
  vector<vector<double> > GetICSpectrum(double emin = 0., double emax = 0.) {
    return ReturnDifferentialPhotonSpectrum(3, emin, emax);
  }  ///< return pi0 decay spectrum
  vector<vector<double> > GetBremsstrahlungSpectrum(double emin = 0.,
                                                    double emax = 0.) {
    return ReturnDifferentialPhotonSpectrum(4, emin, emax);
  }  ///< return Bremsstrahlung spectrum
  vector<vector<double> > GetSynchrotronSpectrum(double emin = 0.,
                                                 double emax = 0.) {
    return ReturnDifferentialPhotonSpectrum(5, emin, emax);
  }  ///< return Bremsstrahlung spectrum
  vector<vector<double> > GetTotalSED(double emin = 0., double emax = 0.) {
    return ReturnSED(1, emin, emax);
  }  ///< return total SED
  vector<vector<double> > GetPPSED(double emin = 0., double emax = 0.) {
    return ReturnSED(2, emin, emax);
  }  ///< return pi0 decay SED
  vector<vector<double> > GetICSED(double emin = 0., double emax = 0.) {
    return ReturnSED(3, emin, emax);
  }  ///< return pi0 decay SED
  vector<vector<double> > GetBremsstrahlungSED(double emin = 0.,
                                               double emax = 0.) {
    return ReturnSED(4, emin, emax);
  }  ///< return Bremsstrahlung SED
  vector<vector<double> > GetSynchrotronSED(double emin = 0.,
                                            double emax = 0.) {
    return ReturnSED(5, emin, emax);
  }  ///< return Bremsstrahlung sed
  Radiation *Clone() { return this; }
  void SetPPEmissionModel(int PIMODEL) {
    PiModel = PIMODEL;
  }  ///< externally switch the parameterisation of the pp emission model. See
     ///PiModel docu for options.
  int GetPPEmissionModel() {
    return PiModel;
  }  ///< return the parameterisation of the pp emission model. See PiModel docu
     ///for options.
  void SetSynchrotronEmissionModel(int SYNCHMODEL) {
    SynchModel = SYNCHMODEL;
  }  ///< externally switch the parameterisation of the synchrotron emission
     ///model. See SynchModel docu for options.
  int GetSynchrotronEmissionModel() {
    return SynchModel;
  }  ///< return the parameterisation of the synchrotron emission model. See
     ///SynchModel docu for options.
  vector< vector<double> > GetProtonSED() {return GetParticleSED("protons");}
  vector< vector<double> > GetElectronSED() {return GetParticleSED("electrons");}
};
#endif
