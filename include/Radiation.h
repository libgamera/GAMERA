#ifndef _RADIATION_
#define _RADIATION_

#include "Utils.h"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>



using namespace std;

typedef double (*fGSLPointer)(double, void *);

/**  @file Radiation.C
 *   @author Joachim Hahn
 *   @short Class that calculates broad band emission from particle spectra.
 */

class Radiation {
  typedef double (Radiation::*fPointer)(double, void *);

  static double evaluate(double x, void* params);
  static fPointer _funcPtr;
  static Radiation *_radPtr;

 private:
  void CalculateLuminosityAndFlux(string mechanism, double e, double &l,
                                  double &f);

  
  double lumtoflux;///< Conversion of differential photon rate to flux
                   /// [ph/(erg*s) -> ph/(erg*s*cm^2)]
  bool LUMFLAG; ///< this boolean is set to TRUE if no distance is given. 
                /// In this case, the luminosity is calculated.
  bool FASTMODE_IC; ///< speed up IC calculation by not calculating emission on
                   /// isotropic target fields individually but only sum
  bool IC_CALCULATED;
  bool FASTMODE_IC_LOSSLOOK; ///< speed up IC calculation by not calculating emission on
                   /// isotropic target fields individually but only sum
  bool IC_LOSSLOOK_CALCULATED;
  double DifferentialEmissionComponent(double e, void *par);
  double GreyBody(double ephoton, double temp, double edens);
  double ICEmissivityRadFieldIntegrated(double x, void *par);
  double ICEmissivity(double x, void *par);
  double ICEmissivityAnisotropic(double x, void *par);
  double ICAnisotropicAuxFunc(double phi_p, double theta_p, double ephoton,
                              double egamma, double lorentz, double beta, 
                              double cos_zeta);

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
  double SynchAngle;///< inclination angle between electrons spiraling along the B-Field and
                    /// observer as used in SynchEmissivityExplicit(). Default: 90degrees.
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
  double GetIntegratedFlux(int i, double emin, double emax, bool ENERGYFLUX=false); /// integrate spectrum i between emin and emax. if ENERGYFLUX == true, the integrated energy flux instead of the integrated flux will be calculated.
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
  vector<vector<double> > TargetPhotonVectorSumAll;  ///< 2D format
                                               ///E(erg)-dN/dE(differential
                                               ///number) vector holding the
                                               ///total target field for IC
                                               ///scattering. This is the sum of
                                               ///individual components
                                               ///specified by public functions
  vector<vector<double> > TargetPhotonVectorSumIso;///< same as TargetPhotonVectorSum
                                                  /// but only summed over isotropic fields
  vector<vector<double> > *TargetPhotonVectorCurrent;  ///< 2D format
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
                   ///CalculateDifferentialGammaEmission) due to IC emission.
                   ///SUM OF ALL TARGET FIELD CONTRIBUTIONS.
  vector<double> fdiffics;///< differential flux at specified energy E (in
                   ///CalculateDifferentialGammaEmission) due to IC emission.
                   ///HOLDS INDIVIDUAL TARGET FIELD CONTRIBUTIONS.
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
  void FillCosZetaLookup(int i);///< fill lookup for angular distance btw. photon and electrons 
                           ///< for anisotropic IC scattering following MoskalenkoStrong1999 
  double d_theta, d_phi; ///< binning of target photon anisotropy map (For MoskalenkoStrong 1999)
  int TargetPhotonSteps;  ///< binning of the target photon spectrum
  bool DEBUG;             ///< debugging boolean
  double integratorTolerance;///
  int integratorKronrodRule;///
  string radiationMechanism;  ///< string query the radiation mechanism that is
                              ///to be calculated
  double hadronicAmpFactor;   ///< amplification factor of the p-p emission (due
                              ///to heavier species in the ISM)
  double BField;              ///< BField value (G)
  bool INTEGRATEOVERGAMMAS;   ///< boolean that switches between calculation IC
                              ///emission loss rate and emission
  bool QUIETMODE;  ///< boolean that toggles quiet output mode if set to true
  bool VERBOSEMODE;
  vector<vector<double> > diffSpec;  ///< vector holding all the individual
                                     ///differential spectra in erg -
                                     ///erg^-1s^-1cm^-2
  vector<vector<double> > diffSpecICComponents;  ///< vector holding all 
                                     /// the individual components of the IC
                                     ///differential spectrum from different 
                                     /// target fields in erg - erg^-1s^-1cm^-2
  int PiModel;  ///< indicates which parameterisation to use in the Kafexhiu pi0
                ///model. 0 - Geant4.10, 1 - Pythia8.1, 2 - SIBYLL2.1, 3 - QGSJET-I.
                ///DEFAULT = 1
  int SynchModel;  ///< indicates which synchrotron emissivity model should be used.
                   /// 0 - isotropic pitch angle distribution, following Ghisellini et al. 1988 
                   /// 1 - fixed pitch angle, default 90 degrees, following Blumenthal&Gould 1970
                   /// if set to 1, angle to observer is controled by'SynchAngle' (default angle 90 degrees.). 
                   ///DEFAULT = 0
  double PPEmissivity(double x, void *par);
  void GetABGParams(double Tp, double &alpha, double &beta, double &gamma,
                    double &lambda);
  void GetAParams(double Tp, double &a1, double &a2, double &a3, double &a4,
                  double &a5);
  void GetBParams(double Tp, double &b1, double &b2, double &b3);
  static double GSLFunctionWrapper(double x, void *params);
  /*    Utils *fUtils;*/
  double Integrate(fPointer f, double *x, double emin, double emax,
                   double tolerance, int pointslevel);  ///< generic integration routine
  gsl_interp_accel *acc, *acciso, *accall,  
                   *loraccesc, *edaccesc, *ICLossLookupAccIso,*ICLossLookupAccAll; ///< gsl accelerator objects 
  gsl_interp_accel **TargetAccCurrent,**phiaccescCurrent, **thetaaccescCurrent,
                   **phiaccesc_zetaCurrent,**thetaaccesc_zetaCurrent,
                   **ICLossLookupAccCurrent;
                                                   ///for interpolation
  gsl_spline *ElectronLookup, *ProtonLookup, *TargetPhotonLookupSumIso,*TargetPhotonLookupSumAll,
       *ICLossLookupSumIso,*ICLossLookupSumAll;
  gsl_spline **TargetPhotonLookupCurrent,**ICLossLookupCurrent;
  interp2d_spline **TargetPhotonAngularDistrCurrent,**CosZetaLookupCurrent;
  unsigned int RADFIELDS_MAX,RADFIELD_CURRENT,RADFIELD_COUNTER;
  vector<double> epoints_temp, TargetPhotonEdensities;
  vector<bool> ANISOTROPY;
  bool ANISOTROPY_CURRENT;
  vector<double> *TargetPhotonAngularBoundsCurrent;
  vector< vector< vector<double> > > TargetPhotonVectors,ICLossVectors,TargetPhotonAngularDistrsVectors;
  vector< gsl_interp_accel * > TargetPhotonAccs, phiaccescs, thetaaccescs,phiaccesc_zetas,thetaaccesc_zetas,ICLossLookupAccs;
  vector<gsl_spline *> TargetPhotonLookups,ICLossLookups;
  vector<interp2d_spline*> TargetPhotonAngularDistrs, CosZetaLookups;
  vector< vector<double> > TargetPhotonAngularBounds,ICLossVectorSumIso,ICLossVectorSumAll,TargetPhotonAngularPhiVectors,TargetPhotonAngularThetaVectors;
  vector< vector<double> > *ICLossVectorCurrent;
  void SetParticles(vector<vector<double> > PARTICLES, int type);
  void SetTargetPhotonVectorLookup(vector< vector<double> > v, int i);
  void SetICLookups(int i);
  const gsl_interp_type *interp_type;
  double TargetPhotonEdensSumIso,TargetPhotonEdensSumAll;
  double *TargetPhotonEdensCurrent;
  vector<vector<double> > GetParticleSED(string type);
  bool SSCSET; ///< boolean that states if SSC target photons have already been added.
  double phi_min,phi_max,theta_min, theta_max;///< bounds of target photon anisotropy map
  double phi_e,theta_e;///< offset to observer-source direction (context of anisotropic IC)
  double ani_minval, ani_maxval; ///< extreme values of target field anisotropy dist.
  double sin_phi_e, cos_phi_e, sin_theta_e, cos_theta_e;
  void SetThermalTargetPhotons(double T, double edens, int steps, int i);
  void SetTargetPhotonsFromFile(
      const char *phFile, int i);  ///< add a target photon density from an ASCII file
                            ///of format E(eV) photon_density(eV^-1cm^-3) !!!
                            ///here eV, as this is what is typically used in the
                            ///literature !!!
  void CreateICLossLookupIndividual(int i=-1, int bins = 100); 
  void SetArbitraryTargetPhotons(
      vector<vector<double> > PhotonArray, int i);  ///< add an arbitrary target photon
                                             ///field. input is a 2D-vector of
                                             ///format E(erg)
                                             ///photon_density(erg^-1cm^-3)

  void SetSSCTargetPhotons(double R, int steps, int i);
  
  vector <double> sizephfield; //size of the photon field in which gamma-gamma absoption has to be considered (allows diffrent size for each field)

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
  void CalculateDifferentialGammaEmission(double e, int particletype);
  void CalculateDifferentialPhotonSpectrum(int steps = 100, double emin = 0.,
                                           double emax = 0.);
  void CalculateDifferentialPhotonSpectrum(vector<double> points);
  vector<vector<double> > ReturnDifferentialPhotonSpectrum(int i,
                                                           double emin,
                                                           double emax,
                                                  vector< vector<double> > vec);
  vector<vector<double> > ReturnSED(int i, double emin,double emax, 
                                      vector< vector<double> > vec);  ///< returns SED for
                                                        ///emission component i
                                                        ///as 2D vector
  void SetBField(double BFIELD) {
    BField = BFIELD;
  }  ///< set the source B-Field (G)
  double GetBField() {return BField;}
  void SetDistance(double d) {
    distance = d*pc_to_cm;
  }  ///< set the distance to the source (cm)
  double GetDifferentialICFlux(int i=-1) { 
    if(i<-1) return 0.;
    else if(i==-1) return fdiffic;
    else return fdiffics[i]; }        ///< get FdiffIC
  double GetDifferentialSynchFlux() { return fdiffsynch; }  ///< get fdiffsynch
  double GetDifferentialBremsFlux() { return fdiffbrems; }  ///< get fdiffbrems
  double GetDifferentialPPFlux() { return fdiffpp; }        ///< get fdiffpp
  double GetIntegralTotalFlux(double emin, double emax) {  ///< get integrated
    return GetIntegratedFlux(1,emin,emax); }             /// gamma-ray flux between
                                                         /// emin and emax (erg) 
                                                         /// summed over all 
                                                         /// radiation processes
                                                         
  double GetIntegralPPFlux(double emin, double emax) { ///< get integrated 
    return GetIntegratedFlux(2,emin,emax); }             /// gamma-ray flux between
                                                         /// emin and emax (erg) 
                                                         /// from proton-proton
                                                         /// interaction
                                                         
  double GetIntegralICFlux(double emin, double emax) { ///< get integrated 
    return GetIntegratedFlux(3,emin,emax); }             /// gamma-ray flux between
                                                         /// emin and emax (erg) 
                                                         /// from IC-mechanism
                                                         
  double GetIntegralBremsstrahlungFlux(double emin, double emax) { ///< get integrated 
    return GetIntegratedFlux(4,emin,emax); }             /// gamma-ray flux between
                                                         /// emin and emax (erg) 
                                                         /// from Bremsstrahlung
                                                         
  double GetIntegralSynchrotronFlux(double emin, double emax) { ///< get integrated 
    return GetIntegratedFlux(5,emin,emax); }             /// flux between
                                                         /// emin and emax (erg) 
                                                         /// from synchrotron
                                                         /// process
                                                         
  double GetIntegralTotalEnergyFlux(double emin, double emax) { ///< get integrated 
    return GetIntegratedFlux(1,emin,emax,true); }        /// energy flux between
                                                         /// emin and emax (erg) 
                                                         /// summed over all 
                                                         /// radiation processes
                                                         
  double GetIntegralPPEnergyFlux(double emin, double emax) { ///< get integrated 
    return GetIntegratedFlux(2,emin,emax,true); }        /// energy flux between
                                                         /// emin and emax (erg) 
                                                         /// from proton-proton
                                                         /// interaction
                                                         
  double GetIntegralICEnergyFlux(double emin, double emax) { ///< get integrated
    return GetIntegratedFlux(3,emin,emax,true); }        /// energy flux between
                                                         /// emin and emax (erg) 
                                                         /// from IC-mechanism
                                                         
  double GetIntegralBremsstrahlungEnergyFlux(double emin, double emax) { ///< get integrated 
    return GetIntegratedFlux(4,emin,emax,true); }        /// energy flux between
                                                         /// emin and emax (erg) 
                                                         /// from Bremsstrahlung
                                                         
  double GetIntegralSynchrotronEnergyFlux(double emin, double emax) { ///< get integrated 
    return GetIntegratedFlux(5,emin,emax,true); }        /// energy flux between
                                                         /// emin and emax (erg) 
                                                         /// from synchrotron
                                                         /// process
                                                         
  void CreateICLossLookup(int bins = 100); ///< creates a 2D vector holding the
                                            ///energy-dependent energy loss rate
                                            ///due to IC cooling. Useful to
                                            ///apply in spectral iterations in
                                            ///conjunction in the "Particles"
                                            ///class. Format: E(erg) -
                                            ///-1.*LossrateIC (erg/s)
  void AddThermalTargetPhotons(double T, double energydens,
                               int steps = 1000);  ///< add a thermal target
                                                  ///radiation field component
                                                  ///to the total local
                                                  ///radiation field
  void ResetWithThermalTargetPhotons(int i, double T, double energydens,
                               int steps = 1000);  ///< add a thermal target
                                                  ///radiation field component
                                                  ///to the total local
                                                  ///radiation field
  void AddArbitraryTargetPhotons(
      vector<vector<double> > PhotonArray);  ///< add an arbitrary target photon
                                             ///field. input is a 2D-vector of
                                             ///format E(erg)
                                             ///photon_density(erg^-1cm^-3)
  void ResetWithArbitraryTargetPhotons(int i,vector<vector<double> > PhotonArray);
  void ImportTargetPhotonsFromFile(const char *phFile);
  void ResetWithTargetPhotonsFromFile(int i,const char *phFile);
  void AddSSCTargetPhotons(double R, int steps = 200);  ///< add target photons
                                                  ///resulting from Synchrotron
                                                  ///radiation due to current
                                                  ///electron spectrum and
                                                  ///B-Field. R (source extension) in pc
  void ResetWithSSCTargetPhotons(int i, double R, int steps = 200); 
  vector<vector<double> > GetTargetPhotons(int i=-1);///< return TotalTargetPhotonVector
  void ClearTargetPhotons(); ///< remove all previously set IC target photons
  void Reset();  ///< reset ParticleLookup, Electrons, Protons, fintbrems,
                 ///lintbrems, fintpp, lintpp, fintic, lintic
  vector<vector<double> > GetICLossLookup(int i=-1);///< return TotalTargetPhotonVector
  double GetDistance() {return distance/pc_to_cm;}
  vector<vector<double> > GetTotalSpectrum(double emin = 0., double emax = 0.) {
    return ReturnDifferentialPhotonSpectrum(1, emin, emax, diffSpec);
  }  ///< return total spectrum
  vector<vector<double> > GetPPSpectrum(double emin = 0., double emax = 0.) {
    return ReturnDifferentialPhotonSpectrum(2, emin, emax, diffSpec);
  }  ///< return pi0 decay spectrum
  vector<vector<double> > GetICSpectrum(double emin = 0., double emax = 0.) {
    return ReturnDifferentialPhotonSpectrum(3, emin, emax, diffSpec);
  }  ///< return total Inverse-Compton spectrum
  vector<vector<double> > GetBremsstrahlungSpectrum(double emin = 0.,
                                                    double emax = 0.) {
    return ReturnDifferentialPhotonSpectrum(4, emin, emax, diffSpec);
  }  ///< return Bremsstrahlung spectrum
  vector<vector<double> > GetSynchrotronSpectrum(double emin = 0.,
                                                 double emax = 0.) {
    return ReturnDifferentialPhotonSpectrum(5, emin, emax, diffSpec);
  }  ///< return Bremsstrahlung spectrum
  vector<vector<double> > GetTotalSED(double emin = 0., double emax = 0.) {
    return ReturnSED(1, emin, emax, diffSpec);
  }  ///< return total SED
  vector<vector<double> > GetPPSED(double emin = 0., double emax = 0.) {
    return ReturnSED(2, emin, emax, diffSpec);
  }  ///< return pi0 decay SED
  vector<vector<double> > GetICSED(double emin = 0., double emax = 0.) {
    return ReturnSED(3, emin, emax, diffSpec);
  }  ///< return pi0 decay SED
  vector<vector<double> > GetBremsstrahlungSED(double emin = 0.,
                                               double emax = 0.) {
    return ReturnSED(4, emin, emax, diffSpec);
  }  ///< return Bremsstrahlung SED
  vector<vector<double> > GetSynchrotronSED(double emin = 0.,
                                            double emax = 0.) {
    return ReturnSED(5, emin, emax, diffSpec);
  }  ///< return Bremsstrahlung sed


  vector<vector<double> > GetICSpectrum(unsigned int i, double emin = 0., double emax = 0.);
  vector<vector<double> > GetICSED(unsigned int i, double emin = 0., double emax = 0.); ///< return pi0 decay spectrum


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
  void CheckSanityOfTargetPhotonLookup();
  vector< vector<double> > GetProtonSED() {return GetParticleSED("protons");}
  vector< vector<double> > GetElectronSED() {return GetParticleSED("electrons");}
  Utils *fUtils;
  void SetTargetPhotonAnisotropy(int i, vector<double> obs_angle, 
                                 vector<double> phi, vector<double> theta, 
                                 vector< vector<double> > mesh);
  vector< vector<double> > GetTargetPhotonAnisotropy(int i, 
                                     vector<double> phi, vector<double> theta);
  void SetSynchrotronPitchAngle(double synchangle) {SynchAngle = synchangle;}
  void SetInterpolationMethod(string intermeth)
    {fUtils->SetInterpolationMethod(intermeth);}
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
  double InelasticPPXSectionKaf(double Tp);
  double InclusivePPXSection(double Tp);
/*  void SetObserverOffsetAngle(double phi, double theta) {*/
/*    phi_e = phi; theta_e = theta; */
/*    sin_phi_e = sin(phi_e); */
/*    cos_phi_e = cos(phi_e);*/
/*    sin_theta_e = sin(theta_e); */
/*    cos_theta_e = cos(theta_e);*/
/*    return;}*/
  void ToggleQuietMode() { QUIETMODE = QUIETMODE == true ? false : true; }  ///< enable quiet mode (very little cout output)
  bool GetQuietMode() {return QUIETMODE;}
  void ToggleVerboseMode() { VERBOSEMODE = VERBOSEMODE == true ? false : true; }  ///< enable quiet mode (very little cout output)
  bool GetVerboseMode() {return VERBOSEMODE;}
  double ICEmissivityWrapper(double e_ph, double e_e, double e_g);
  double ICEmissivityAnisotropicWrapper(double e_ph, double e_e, double e_g);
/*  vector< vector<double> > GetTargetPhotonFieldVector(unsigned int i){*/
/*    vector< vector<double> >v = fUtils->VectorAxisPow10(TargetPhotonVectors[i],-1);*/
/*    return v;}*/
  double GetTargetPhotonFieldEnergyDensity(unsigned int i) {
    return TargetPhotonEdensities[i];
  }
  vector< vector<double> > GetTargetFieldAnisotropyMap(int i) {
    if (i>(int)RADFIELDS_MAX) cout<<"Radiation::GetTargetFieldAnisotropyMap: invalid index "<<i<<". Returning"
          "empty vector." <<endl;
    return TargetPhotonAngularDistrsVectors[i];}

  void ClearTargetPhotonField(int i); ///< reset values for target photon 
                                          ///< field component i
  void RemoveLastICTargetPhotonComponent() {
    ClearTargetPhotonField(--RADFIELD_COUNTER);
  }
  vector< vector<double> > SumTargetFields(int bins,bool ISO=true);
  void SumTargetFieldsAll(int bins=1000);
  void SumTargetFieldsIsotropic(int bins=1000);
  unsigned int GetTargetFieldCount(){return RADFIELD_COUNTER;}
  void SetICFastMode() {FASTMODE_IC = true;}
  
  /* 
   * Adding functions for the gamma-gamma absorption calculations
   */
  void SetSizePhotonField(vector <double> size){
    // if not inizialized, the default value is 1 pc
    // size to be given in [pc]
      for (unsigned int i=0;i<size.size();i++){
          sizephfield.push_back(size[i] * pc_to_cm);
      }
  } // Sets the size of the photon field for which gamma-gamma absorption has to be considered
  
  void ClearPhotonFieldSize();
  
  vector <double> GetSizePhotonField(); //get the photon filed
  
  double AverageSigmaGammaGamma(double Eph1, double Eph2);             // Average cross section for isotropic and homogeneous case
  double SigmaGammaGamma(double Eph1, double Eph2, double costheta);      // Full gamma-gamma cross section
  double ComputeAbsCoeff(double Egamma, int target);  //Auxiliary function to compute only the absorption coefficient, no spatial integration
  double ComputeOptDepth(double Egamma, int target, double phsize);
  double ComputeOptDepthIsotropic(double Egamma, int target, double phsize); // Computation of the optical depth parameter isotropic only
  vector< vector<double> > ReturnAbsorbedSEDonFields(double emin, double emax, 
                                                         vector <int> fields, vector <double> size); // Wrapper around function ReturnSED to return the gamma-gamma absorbed values
  vector< vector<double> > ReturnAbsorbedSpectrumOnFields(double emin, double emax, 
                                                         vector <int> fields, vector <double> size); // Wrapper around function ReturnDifferentialPhotonSpectrum to return the gamma-gamma absorbed values
  vector<vector<double> > GetTotalAbsorbedSpectrum(vector <int> fields, double emin = 0., double emax = 0.) {
    return ReturnAbsorbedSpectrumOnFields(emin, emax, fields, sizephfield);
  }
  vector<vector<double> > GetTotalAbsorbedSED(vector <int> fields, double emin = 0., double emax = 0.) {
    return ReturnAbsorbedSEDonFields(emin, emax, fields, sizephfield);
  }  ///< return total absorbed SED on selected fields, default behaviour should be to use all of them. TODO
  double ReturnAbsorbedIntergratedFlux(double emin, double emax, bool ENERGYFLUX,vector <int> fields, vector <double> size);
  
  double GetIntegralTotalAbsEnergyFlux(double emin, double emax, vector <int> fields) { ///< get integrated 
    return ReturnAbsorbedIntergratedFlux(emin,emax,true,fields,sizephfield); }                      /// absorbed energy flux between
                                                                                          /// emin and emax (erg) 
                                                                                          /// summed over all 
                                                                                          /// radiation processes
  double GetIntegralTotalAbsFlux(double emin, double emax, vector <int> fields) { ///< get integrated 
    return ReturnAbsorbedIntergratedFlux(emin,emax,false,fields,sizephfield); }                      /// absorbed energy flux between
                                                                                          /// emin and emax (erg) 
                                                                                          /// summed over all 
                                                                                          /// radiation processes
};
#endif
