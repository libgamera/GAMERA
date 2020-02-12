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
                   ///< [ph/(erg*s) -> ph/(erg*s*cm^2)]
  bool LUMFLAG; ///< this boolean is set to TRUE if no distance is given. 
                ///< In this case, the luminosity is calculated.
  bool FASTMODE_IC; ///< speed up IC calculation by not calculating emission on
                   ///< isotropic target fields individually but only sum
  bool IC_CALCULATED;
  bool ISOTROPIC_ELECTRONS; ///< Set to true for isotropic electrons in IC scattering

  bool FASTMODE_IC_LOSSLOOK; ///< speed up IC calculation by not calculating emission on
                             ///< isotropic target fields individually but only sum
  bool IC_LOSSLOOK_CALCULATED;  ///< Store if the Inverse Compton Lookup has been computed
  double DifferentialEmissionComponent(double e, void *par); ///< Calculates the differential photon rate
  double GreyBody(double ephoton, double temp, double edens); ///< Thermal grey body
  double ICEmissivityRadFieldIntegrated(double x, void *par); ///< IC emission from electrons integrated over the target photon population.
  double ICEmissivity(double x, void *par);           ///< Emissivity from IC scattering
  double ICEmissivityAnisotropic(double x, void *par);  ///< Anisitropic IC emissivity
  double ICAnisotropicAuxFunc(double phi_p, double theta_p, double ephoton,
                              double egamma, double lorentz, double beta, 
                              double cos_zeta); ///< Auxiliary function for IC scattering 
  double ICEmissivityAnisotropicIsotropicElectrons(double x, void *par);  ///< Anisotropic IC emissivity with isotropic electrons
  double ICEmissivityAnisotropicIsotropicElectronsSecondIntegral(double x, void *par);
  double ICEmissivityAnisotropicIsotropicElectronsFirstIntegral(double x, void *par);
  double ICEmissivityAnisotropicSecondIntegral(double x, void *par);
  double ICEmissivityAnisotropicFirstIntegral(double x, void *par);
  double ICEmissivityAnisotropicIsotropicElectronsAharonian(double x, void *par);
  double ICEmissivityAnisotropicIsotropicElectronsAharonianSecondIntegral(double x, void *par);
  double ICEmissivityAnisotropicIsotropicElectronsAharonianFirstIntegral(double x, void *par);
  double K(double nu, double x);             ///< modified Bessel function
  double K_53(double x, void *par);          ///< modified Bessel function of order 5/3
  double SynchEmissivity(double x, void *par);  ///< Synchrotron Emissivity
  double SynchEmissivityExplicit(double e,
                                 void *par);  ///< Synchrotron Emissivity with pitch angle
  double SynchAngle;///< inclination angle between electrons spiraling along the B-Field and
                    ///< observer as used in Radiation::SynchEmissivityExplicit. DEFAULT: 90degrees.
  double BremsEmissivity(double x, void *par);  ///< Bremsstrahlung emissivity for e-p and e-e  

  double epsilonc; ///< epsilon parameters from Kafhexiu paper
  double epsilon1;
  double epsilon2;

  
  double A(double g,
           double e);  ///< Equation (A4) from Baring 1999, ApJ, 513, 311-338
  double sigma1(double g, double e);   ///< Equation (A2) from Baring 1999, ApJ,
                                       ///< 513, 311-338
  double sigma2(double g, double e);   ///< Equation (A3) from Baring 1999, ApJ,
                                       ///< 513, 311-338
  double sigmaNR(double g, double e);  ///< Equation (A5) from Baring 1999, ApJ,
                                       ///< 513, 311-338
  double Fbr(double x, double g);  ///< Equations (A6,A7) from Baring 1999, ApJ,
                                   ///< 513, 311-338
  double ParticleSpectrum(double E);  ///< DEPRECATED: particle spectrum emitting the
                                      ///< emission. This is basically reading in
                                      ///< a 2D vector(i.e. 'ParticleLookup')
                                      ///< holding E-N(particles).
  double GetIntegratedFlux(int i, double emin, double emax, bool ENERGYFLUX=false); ///< Integrate spectrum "i" between emin and emax. if ENERGYFLUX == true, the integrated energy flux instead of the integrated flux will be calculated.
  
  /** \brief 2D vector holding the particle spectrum that emits the radiation. */
  vector<vector<double> > ParticleVector;  ///< Format: E (erg)-dN/dE (number). This can
                                           ///< be either electrons or protons,
                                           ///< and is set to either species
                                           ///< automatically depending on the
                                           ///< radiation mechanism that needs to
                                           ///< be calculated.
  vector<vector<double> > ElectronVector;  ///< format of Radiation::ParticleVector, holding the electron spectrum
  vector<vector<double> > ProtonVector;    ///< format of Radiation::ParticleVector, holding the proton spectrum
  /**\brief 2D Vector holding the total target field for IC scattering*/
  
  
  vector<double> ProjHadronMassVector;            ///holding the hadron mass vector of projectile spectrum
  vector<vector<double> > ProjHadronsVector; /// holding the relative hadron abundances in projectile spectrum

  vector<double> MedHadronMassVector;            ///holding the hadron mass vector in medium of propagation
  vector<double> MedRelativeHadronAbundancesVector; /// holding the relative hadron abundances in medium of propagation

  vector<double> ProjRelativeHadronsAbunancesVector;
  void CalculateProjRelativeAbundancesLookup(int bins=100);
  void CalculateEpsilonLookups(int bins=100);
  vector<gsl_spline *> ProjRelativeAbundanceLookups;
  vector<gsl_interp_accel *> ProjRelativeAbundanceLookupsAcc;
  vector< vector<double> > ProjRelativeAbundanceLookupsRanges;
  vector< gsl_spline * > EpsilonLookups;
  vector< gsl_interp_accel * > EpsilonLookupsAcc;
  vector< vector<double> > EpsilonLookupsRanges;
  vector< vector<double> > CreateVectorForProjRelativeAbundances(vector< double > x,
                                                      vector< double > y);
  
  vector<vector<double> > TargetPhotonVectorSumAll;  ///< Format: 
                                               ///< E(erg)-dN/dE(differential
                                               ///< number). The vector holds the
                                               ///< total target field for IC
                                               ///< scattering. This is the sum of
                                               ///< individual components
                                               ///< specified by public functions
  vector<vector<double> > TargetPhotonVectorSumIso;///< same as Radiation::TargetPhotonVectorSumAll
                                                   ///< but only summed over isotropic fields
  vector<vector<double> > *TargetPhotonVectorCurrent;  ///< Pointer to the current target field in use
  /**\brief Target photon field for Synchrotron Self Compton mechanism */
  vector<vector<double> > SSCTargetPhotons;  ///< 2D-vector object holding the
                                             ///< SSC target photon spectrum  {
                                             ///< E(erg) - Edens(erg cm^-3) }
  double n;  ///< ambient density for Bremsstrahlung and inelastic p-p emission mechanisms (cm^-3)
  double fdiffbrems;  ///< differential flux at specified energy E (in
                      ///< Radiation::CalculateDifferentialGammaEmission) due to
                      ///< Bremsstrahlung
  double ldiffbrems;  ///< differential luminosity at specified energy E (in
                      ///< Radiation::CalculateDifferentialGammaEmission) due to
                      ///< Bremsstrahlung
  double fdiffpp;     ///< differential flux at specified energy E (in
                   ///< Radiation::CalculateDifferentialGammaEmission) due to inelastic p-p
                   ///< collision
  double ldiffpp;  ///< differential luminosity at specified energy E (in
                   ///< Radiation::CalculateDifferentialGammaEmission) due to inelastic p-p
                   ///< collision
  double fdiffic;  ///< differential flux at specified energy E (in
                   ///< Radiation::CalculateDifferentialGammaEmission) due to IC emission.
                   ///< SUM OF ALL TARGET FIELD CONTRIBUTIONS.
  vector<double> fdiffics;///< differential flux at specified energy E (in
                   ///< Radiation::CalculateDifferentialGammaEmission) due to IC emission.
                   ///< HOLDS INDIVIDUAL TARGET FIELD CONTRIBUTIONS.
  double ldiffic;  ///< differential luminosity at specified energy E (in
                   ///< Radiation::CalculateDifferentialGammaEmission) due to IC emission
  double fdiffsynch;  ///< differential flux at specified energy E (in
                      ///< Radiation::CalculateDifferentialGammaEmission) Synchrotron
                      ///< emission
  double ldiffsynch;  ///< differential luminosity at specified energy E (in
                      ///< Radiation::CalculateDifferentialGammaEmission) Synchrotron
                      ///< emission
  double distance;    ///< source distance (cm)
  double targetphotonenergymin;  ///< lower energy boundary for the target
                                 ///< photon spectrum (erg, dynamically
                                 ///< determined)
  double targetphotonenergymax;  ///< upper energy boundary for the target
                                 ///< photon spectrum (erg, dynamically
                                 ///< determined)
  void FillTargetPhotonVectorAndGraph(int steps);  ///< DEPRECATED: function that adds
                                                   ///< individual target photon
                                                   ///< spectrum to
                                                   ///< "TotalTargetPhotonGraph"
                                                   ///< and
                                                   ///< "TotalTargetPhotonVector"
  void FillCosZetaLookup(int i);///< fill lookup for angular distance btw. photon and electrons for anisotropic IC scattering
  double d_theta, d_phi; ///< binning of target photon anisotropy map (For MoskalenkoStrong 1999)
  int TargetPhotonSteps;  ///< binning of the target photon spectrum
  bool DEBUG;             ///< debugging boolean
  double integratorTolerance; ///< Tolerance for gsl integration DEFAULT = 0.1
  int integratorKronrodRule;  ///< Integration rule for gsl integration DEFAULT = 2
  string radiationMechanism;  ///< string query the radiation mechanism that is to be calculated
  double hadronicAmpFactor;   ///< amplification factor of the p-p emission (due to heavier species in the ISM)
  double BField;              ///< BField value (in Gauss)
  bool INTEGRATEOVERGAMMAS;   ///< boolean that switches between calculation IC emission loss rate and emission
  bool QUIETMODE;  ///< boolean that toggles quiet output mode if set to true
  bool VERBOSEMODE;
  bool CALCULATEHADRONMIX; ///< boolean that indicates if hadron-mix related things have been calculated in  NuclearEnhancementFactor()
  /**\brief Vector holding all the individual differential spectra*/
  vector<vector<double> > diffSpec;  ///< Format: (erg) vs. (erg^-1 s^-1 cm^-2)
  vector<vector<double> > diffSpecICComponents;  ///< vector holding all 
                                     /// the individual components of the IC
                                     ///differential spectrum from different 
                                     /// target fields in erg - erg^-1s^-1cm^-2
 /**
  * \brief hadronic model for pi0 emission
  * 
  * Indicates which parameterisation to use in the Kafexhiu pi0 model. 
  * 
  * \li 0 - Geant4.10, 
  * \li 1 - Pythia8.1, 
  * \li 2 - SIBYLL2.1, 
  * \li 3 - QGSJET-I.
  * 
  * DEFAULT = 1
  */
  int PiModel;  
  
  /** 
   * \brief Define Synchrotron model being used
   * 
   * Indicates which synchrotron emissivity model should be used.
   * \li 0 = isotropic pitch angle distribution, following Ghisellini et al. 1988 
   * \li 1 = fixed pitch angle, default 90 degrees, following Blumenthal&Gould 1970
   * 
   * If set to 1, the angle to the observer is controled by Radiation::SynchAngle
   * (DEFAULT angle 90 degrees)
   */
  int SynchModel;  
  double PPEmissivity(double x, void *par);  ///< Compute pi0 emissivity
  void GetABGParams(double Tp, double &alpha, double &beta, double &gamma,
                    double &lambda); ///< auxiliary function for Radiation::PPEmissivity
  void GetAParams(double Tp, double &a1, double &a2, double &a3, double &a4,
                  double &a5); ///< auxiliary function for Radiation::PPEmissivity
  void GetBParams(double Tp, double &b1, double &b2, double &b3); ///< auxiliary function for Radiation::PPEmissivity
  static double GSLFunctionWrapper(double x, void *params);  
  /*    Utils *fUtils;*/
  double Integrate(fPointer f, double *x, double emin, double emax,
                   double tolerance, int pointslevel);  ///< Generic integration routine
  gsl_interp_accel *acc, *acciso, *accall,  
                   *loraccesc, *edaccesc, *ICLossLookupAccIso,*ICLossLookupAccAll; ///< gsl accelerator objects 
  gsl_interp_accel **TargetAccCurrent,**phiaccescCurrent, **thetaaccescCurrent,
                   **phiaccesc_zetaCurrent,**thetaaccesc_zetaCurrent,
                   **ICLossLookupAccCurrent;  ///< gsl accelerator objects for interpolation
                                                   
  gsl_spline *ElectronLookup, *ProtonLookup, *TargetPhotonLookupSumIso,*TargetPhotonLookupSumAll,
       *ICLossLookupSumIso,*ICLossLookupSumAll; ///< Lookup tables for particles and IC
  gsl_spline **TargetPhotonLookupCurrent,**ICLossLookupCurrent;  ///< Lookup tables for IC scattering on current field
  interp2d_spline **TargetPhotonAngularDistrCurrent,**CosZetaLookupCurrent; ///< 2D lookups for IC anisotropy calculations
  unsigned int RADFIELDS_MAX,RADFIELD_CURRENT,RADFIELD_COUNTER; ///< Counters for radiation fields
  vector<double> epoints_temp;  ///< temp variable for the energy values of the resulting emission spectrum
  vector<double> TargetPhotonEdensities;  ///< Vector of energy density of the target photon fields
  vector<bool> ANISOTROPY;  ///< Vector to store if atarget field is anisotripic or not
  bool ANISOTROPY_CURRENT;  ///< To state if current vector is anisotropic or not
  vector<double> *TargetPhotonAngularBoundsCurrent;  ///< Angular bounds for current target photon field
  
  /** \brief Vector storing all the target photon fields */
  vector< vector< vector<double> > > TargetPhotonVectors; ///< Format: vector of vectors of vectors of (Ephoton [erg] vs. differential photon density [N/cm^3/erg])
  vector< vector< vector<double> > > ICLossVectors; ///< Vector storing the IC loss rate for each field
  vector< vector< vector<double> > > TargetPhotonAngularDistrsVectors;  ///< Vector storing the angular distribution for each target photon field
  vector< gsl_interp_accel * > TargetPhotonAccs, phiaccescs, thetaaccescs,phiaccesc_zetas,thetaaccesc_zetas,ICLossLookupAccs; ///< vectors of gsl accelerator objects
  vector<gsl_spline *> TargetPhotonLookups; ///< Vectors of Lookups for the every target photon field
  vector<gsl_spline *> ICLossLookups; ///< Vectors of Lookups for the IC loos rates
  vector<interp2d_spline*> TargetPhotonAngularDistrs, CosZetaLookups;
  vector< vector<double> > TargetPhotonAngularBounds,ICLossVectorSumIso,ICLossVectorSumAll,TargetPhotonAngularPhiVectors,TargetPhotonAngularThetaVectors;
  vector< vector<double> > *ICLossVectorCurrent; ///< Vector for current ICLossVector filled in Radiation::CreateICLossLookupIndividual
  
  void SetParticles(vector<vector<double> > PARTICLES, int type);   ///< Fill the lookup for the particle spectrum. 0 for electrons, 1 for protons
  void SetTargetPhotonVectorLookup(vector< vector<double> > v, int i);  ///< Fill the lookup for the target photon field
  void SetICLookups(int i);  ///< Fill in all the IC related Lookups for field \a i
  const gsl_interp_type *interp_type;
  double TargetPhotonEdensSumIso,TargetPhotonEdensSumAll;
  double *TargetPhotonEdensCurrent;
  vector<vector<double> > GetParticleSED(string type);
  bool SSCSET; ///< boolean that states if SSC target photons have already been added.
  //@{
  double phi_min,phi_max,theta_min, theta_max; ///< bounds of target photon anisotropy map
  //@}
  //@{
  double phi_e,theta_e;///< offset to observer-source direction (context of anisotropic IC)
  //@}
  //@{
  double ani_minval, ani_maxval; ///< extreme values of target field anisotropy dist.
  //@}
  //@{
  double sin_phi_e, cos_phi_e, sin_theta_e, cos_theta_e; ///< sin and cos of offset to observer-source direction (context of anisotropic IC)
  //@}
  void SetThermalTargetPhotons(double T, double edens, int steps, int i); ///< Set a grey body photon field
  void SetTargetPhotonsFromFile(
      const char *phFile, int i);  ///< add a target photon density from an ASCII file
                            
  void CreateICLossLookupIndividual(int i=-1, int bins = 100); ///< Fills a lookup table for the IC loos rate
  void SetArbitraryTargetPhotons(
      vector<vector<double> > PhotonArray, int i);  ///< add an arbitrary target photon field.

  void SetSSCTargetPhotons(double R, int steps, int i); ///< Add SSC target photons.
  
  vector <double> sizephfield; ///< size of the photon field in which gamma-gamma absorption has to be considered
                               ///(allows different size for each field)

  vector < vector <vector<double> > > SpatialDep;  ///< Triple vector to store the spatial dependency for each photon field
  //gsl_spline *SpatialDepLookup; ///< Lookup for the spatial dependency for the gamma-gamma absorption
  //gsl_interp_accel *accsp; ///< gsl accelerator for the spatial lookup
  bool SPATIALDEP_CURRENT; ///< boolean state if the absorption target field has a homogeneous density or not
  vector <bool> SPATIALDEP; ///< vector of booleans to state which photon field has a spatial dependency


 public:
  Radiation();                                       ///< standard constructor
  ~Radiation();                                      ///< standard destructor
  void SetProtons(vector<vector<double> > PROTONS);  ///< Set Protons
  void SetElectrons(vector<vector<double> > ELECTRONS);  ///< Set Electrons
  void SetElectronsIsotropic(void);    ///<Setting isotropic electrons variable to true for anisotropic IC scattering
  
  void SetProjHadronMass(vector<double> PROJ_HADRON_MASS);
  void SetMedHadronMass(vector<double> MED_HADRON_MASS);
  
  vector<double> GetProjHadronMassVector(){
    return ProjHadronMassVector;  ///return hadron mass vector of the projectile spectrum
  }

  vector<double> GetMedHadronMassVector(){
    return MedHadronMassVector;  ///return hadron mass vector in medium of propagation
  }
   
  void SetProjHadrons(vector<vector<double> > PROJ_HADRONS, int bins = 100); ///Set the relative abundances vector for different species in comparison to proton spectrum in projectile spectrum
  
  vector<vector<double> > GetProjHadronsVector(){
    return ProjHadronsVector;
  }

  void SetMedRelativeHadronAbundances(vector<double> MED_REL_HADRONS_ABUN); ///Set the relative abundances vector for different species in comparison to proton spectrum in medium of propagation
  
  vector<double> GetMedRelativeHadronAbundancesVector(){
    return MedRelativeHadronAbundancesVector;
  }

  void CalculateProjRelativeAbundancesVectors(double Tp);

  vector<double> GetProjRelativeHadronAbundancesVector(){
    return ProjRelativeHadronsAbunancesVector;
  }
  
  double NuNuXSection(double ProjMass, double TargetMass);
  double CalculateEpsilonc(vector<double> ProjMass, vector<double> ProjRelAbun, vector<double> MedMass, vector<double> MedRelAbun);
  double CalculateEpsilon1(vector<double> ProjMass, vector<double> ProjRelAbun, vector<double> MedMass, vector<double> MedRelAbun);
  double CalculateEpsilon2(vector<double> ProjMass, vector<double> ProjRelAbun, vector<double> MedMass, vector<double> MedRelAbun);
  vector<double> CalculateEpsilons(double Tp);
  
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
  }  ///< set the distance to the source (\a d to be given in units of parsecs)
  double GetDifferentialICFlux(int i=-1) { 
    if(i<-1) return 0.;
    else if(i==-1) return fdiffic;
    else return fdiffics[i]; }        ///< get Radiation::fdiffIC. \a i to retrieve different components. With \a i = -1 (DEFAULT), returns the total IC components
  double GetDifferentialSynchFlux() { return fdiffsynch; }  ///< get Radiation::fdiffsynch
  double GetDifferentialBremsFlux() { return fdiffbrems; }  ///< get Radiation::fdiffbrems
  double GetDifferentialPPFlux() { return fdiffpp; }        ///< get Radiation::fdiffpp
  double GetIntegralTotalFlux(double emin, double emax) {  
    return GetIntegratedFlux(1,emin,emax); }               ///< get integrated gamma-ray flux between
                                                           ///< emin and emax (erg) 
                                                           ///< summed over all 
                                                           ///< radiation processes
                                                           ///< Using Radiation::GetIntegratedFlux
                                                         
  double GetIntegralPPFlux(double emin, double emax) {  
    return GetIntegratedFlux(2,emin,emax); }           ///< get integrated gamma-ray flux between
                                                       ///< emin and emax (erg) 
                                                       ///< from proton-proton
                                                       ///< interaction
                                                       ///< Using Radiation::GetIntegratedFlux
                                                         
  double GetIntegralICFlux(double emin, double emax) { 
    return GetIntegratedFlux(3,emin,emax); }           ///< get integrated gamma-ray flux between
                                                       ///< emin and emax (erg) 
                                                       ///< from IC-mechanism
                                                       ///< Using Radiation::GetIntegratedFlux
                                                         
  double GetIntegralBremsstrahlungFlux(double emin, double emax) { 
    return GetIntegratedFlux(4,emin,emax); }             ///< get integrated 
                                                         ///< gamma-ray flux between
                                                         ///< emin and emax (erg) 
                                                         ///< from Bremsstrahlung
                                                         ///< Using Radiation::GetIntegratedFlux
                                                         
  double GetIntegralSynchrotronFlux(double emin, double emax) { 
    return GetIntegratedFlux(5,emin,emax); }               ///< get integrated 
                                                           ///< flux between
                                                           ///< emin and emax (erg) 
                                                           ///< from synchrotron
                                                           ///< process.
                                                           ///< Using Radiation::GetIntegratedFlux
                                                         
  double GetIntegralTotalEnergyFlux(double emin, double emax) { 
    return GetIntegratedFlux(1,emin,emax,true); }        ///< get integrated energy flux between
                                                         ///< emin and emax (erg) 
                                                         ///< summed over all 
                                                         ///< radiation processes
                                                         ///< Using Radiation::GetIntegratedFlux
                                                         
  double GetIntegralPPEnergyFlux(double emin, double emax) { 
    return GetIntegratedFlux(2,emin,emax,true); }        ///< get integrated energy flux between
                                                         ///< emin and emax (erg) 
                                                         ///< from proton-proton
                                                         ///< interaction
                                                         ///< Using Radiation::GetIntegratedFlux
                                                         
  double GetIntegralICEnergyFlux(double emin, double emax) { 
    return GetIntegratedFlux(3,emin,emax,true); }        ///< get integrated energy flux between
                                                         ///< emin and emax (erg) 
                                                         ///< from IC-mechanism
                                                         ///< Using Radiation::GetIntegratedFlux
                                                         
  double GetIntegralBremsstrahlungEnergyFlux(double emin, double emax) {  
    return GetIntegratedFlux(4,emin,emax,true); }        ///< get integrated energy flux between
                                                         ///< emin and emax (erg) 
                                                         ///< from Bremsstrahlung
                                                         ///< Using Radiation::GetIntegratedFlux
                                                         
  double GetIntegralSynchrotronEnergyFlux(double emin, double emax) {  
    return GetIntegratedFlux(5,emin,emax,true); }        ///< get integrated energy flux between
                                                         ///< emin and emax (erg) 
                                                         ///< from synchrotron
                                                         ///< process
                                                         ///< Using Radiation::GetIntegratedFlux
                                                         
  void CreateICLossLookup(int bins = 100); ///< Creates a 2D vector holding the energy-dependent energy loss rate due to IC cooling.
 
  
  void AddThermalTargetPhotons(double T, double energydens,
                               int steps = 1000);  ///< add a thermal target radiation field component
                                                  
  void ResetWithThermalTargetPhotons(int i, double T, double energydens,
                               int steps = 1000);  ///< Resets photon field i with another thermal field
  void AddArbitraryTargetPhotons(
      vector<vector<double> > PhotonArray);  ///< add an arbitrary target photon field.
  void ResetWithArbitraryTargetPhotons(int i,vector<vector<double> > PhotonArray);  ///< Resets photon field i with another arbitrary field
  void ImportTargetPhotonsFromFile(const char *phFile);  ///< Add photon field from file
  void ResetWithTargetPhotonsFromFile(int i,const char *phFile);  ///< Reset field i with field from file
  void AddSSCTargetPhotons(double R, int steps = 200);  ///< Add target photons resulting from Synchrotron radiation
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
  }  ///< return Synchrotron spectrum
  vector<vector<double> > GetTotalSED(double emin = 0., double emax = 0.) {
    return ReturnSED(1, emin, emax, diffSpec);
  }  ///< return total SED
  vector<vector<double> > GetPPSED(double emin = 0., double emax = 0.) {
    return ReturnSED(2, emin, emax, diffSpec);
  }  ///< return pi0 decay SED
  vector<vector<double> > GetICSED(double emin = 0., double emax = 0.) {
    return ReturnSED(3, emin, emax, diffSpec);
  }  ///< return Inverse Compton SED
  vector<vector<double> > GetBremsstrahlungSED(double emin = 0.,
                                               double emax = 0.) {
    return ReturnSED(4, emin, emax, diffSpec);
  }  ///< return Bremsstrahlung SED
  vector<vector<double> > GetSynchrotronSED(double emin = 0.,
                                            double emax = 0.) {
    return ReturnSED(5, emin, emax, diffSpec);
  }  ///< return Synchrotron SED


  vector<vector<double> > GetICSpectrum(unsigned int i, double emin = 0., double emax = 0.); ///< return Inverse Compton spectrum
  vector<vector<double> > GetICSED(unsigned int i, double emin = 0., double emax = 0.); ///< return Inverse Compton SED


  Radiation *Clone() { return this; }
  void SetPPEmissionModel(int PIMODEL) {
    PiModel = PIMODEL;
  }  ///< externally switch the parameterisation of the pp emission model. See Radiation::PiModel documentation for options.
  int GetPPEmissionModel() {
    return PiModel;
  }  ///< return the parameterisation of the pp emission model. See Radiation::PiModel docu for options.
  void SetSynchrotronEmissionModel(int SYNCHMODEL) {
    SynchModel = SYNCHMODEL;
  }  ///< externally switch the parameterisation of the synchrotron emission model. See Radiation::SynchModel docu for options.
  int GetSynchrotronEmissionModel() {
    return SynchModel;
  }  ///< return the parameterisation of the synchrotron emission model. See Radiation::SynchModel docu for options.
  void CheckSanityOfTargetPhotonLookup();
  vector< vector<double> > GetProtonSED() {return GetParticleSED("protons");}
  vector< vector<double> > GetElectronSED() {return GetParticleSED("electrons");}
  Utils *fUtils;
  void SetTargetPhotonAnisotropy(int i, vector<double> obs_angle, 
                                 vector<double> phi, vector<double> theta, 
                                 vector< vector<double> > mesh);
  void SetTargetPhotonAnisotropy(int i, vector<double> phi, vector<double> theta, 
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
  void ToggleVerboseMode() { VERBOSEMODE = VERBOSEMODE == true ? false : true; }  ///< enable verbose mode (more cout output)
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

  void ClearTargetPhotonField(int i); ///< reset values for target photon field component i
  void RemoveLastICTargetPhotonComponent() {
    ClearTargetPhotonField(--RADFIELD_COUNTER);
  }
  vector< vector<double> > SumTargetFields(int bins,bool ISO=true);
  void SumTargetFieldsAll(int bins=1000);
  void SumTargetFieldsIsotropic(int bins=1000);

  unsigned int GetTargetFieldCount(){return RADFIELD_COUNTER;}
  void SetICFastMode() {FASTMODE_IC = true;}
  void UnsetICFastMode() {FASTMODE_IC = false;}
  
  void SetSizePhotonField(int i, double size);
  
  void ClearPhotonFieldSize();
  
  vector <double> GetSizePhotonField(); //get the photon field
  double GetSizePhotonField(int i);
  void SetTargetFieldSpatialDep(int i, vector< vector<double> > SpatialDep); ///< Set the spatial dependence for the Target field
  vector< vector<double> > GetTargetFieldSPatialDep(int i);  ///< Return the spatial dependency of the photon field (space in cm)
  double AverageSigmaGammaGamma(double Eph1, double Eph2);             // Average cross section for isotropic and homogeneous case
  double SigmaGammaGamma(double Eph1, double Eph2, double costheta);      // Full gamma-gamma cross section
  double ComputeAbsCoeff(double Egamma, int target);  ///< Auxiliary function to compute only the absorption coefficient, no spatial integration
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
