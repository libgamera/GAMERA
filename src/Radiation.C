#include "Radiation.h"
/**
 * Constructor. Nothing fancy, just initialises some stuff.
 */
Radiation::Radiation() {
  /* Default values */
  fUtils = new Utils();
  ParticleVector.clear();
  DEBUG = false;
  FASTMODE_IC = true;
  FASTMODE_IC_LOSSLOOK = true;
  IC_LOSSLOOK_CALCULATED = false;
  IC_CALCULATED = false;
  SSCSET = false;
  ANISOTROPY_CURRENT = false;
  ISOTROPIC_ELECTRONS = false;  //If true, calculate anisotropip IC scattering with isotropic electrons
  lumtoflux = 0.;
  ldiffbrems = fdiffbrems = ldiffsynch = fdiffsynch = 0.;
  ldiffic = fdiffic = ldiffpp = fdiffpp = 0.;
  distance = BField = 0.;
  phi_min = phi_max = theta_min = theta_max = phi_e = theta_e = 0.;
  sin_phi_e = cos_phi_e = sin_theta_e = cos_theta_e = 0.;
  ani_minval = ani_maxval = 0.;
  ElectronLookup = NULL;
  ProtonLookup = NULL;
  TargetPhotonLookupSumIso = NULL;
  TargetPhotonLookupSumAll = NULL;
  ICLossLookupSumIso = NULL;
//  TargetPhotonLookupCurrent = NULL;
//  TargetPhotonAngularDistrCurrent = NULL;
  RADFIELD_COUNTER = 0;
  RADFIELDS_MAX = 1000;
  RADFIELD_CURRENT = -1;
  TargetPhotonVectors.resize(RADFIELDS_MAX);
  TargetPhotonLookups.resize(RADFIELDS_MAX);
  ICLossVectors.resize(RADFIELDS_MAX);
  ICLossLookups.resize(RADFIELDS_MAX);
  fdiffics.resize(RADFIELDS_MAX);
  TargetPhotonAccs.resize(RADFIELDS_MAX);
  TargetPhotonEdensities.resize(RADFIELDS_MAX);
  TargetPhotonAngularDistrs.resize(RADFIELDS_MAX);
  TargetPhotonAngularBounds.resize(RADFIELDS_MAX);
  CosZetaLookups.resize(RADFIELDS_MAX);
  phiaccescs.resize(RADFIELDS_MAX);
  thetaaccescs.resize(RADFIELDS_MAX);
  phiaccesc_zetas.resize(RADFIELDS_MAX);
  thetaaccesc_zetas.resize(RADFIELDS_MAX);
  ICLossLookupAccs.resize(RADFIELDS_MAX);
  ANISOTROPY.resize(RADFIELDS_MAX);
  TargetPhotonAngularDistrsVectors.resize(RADFIELDS_MAX);
  TargetPhotonAngularPhiVectors.resize(RADFIELDS_MAX);
  TargetPhotonAngularThetaVectors.resize(RADFIELDS_MAX);
  for(unsigned int i=0;i<RADFIELDS_MAX;i++) {
    fdiffics[i] = NAN;
    TargetPhotonEdensities[i] = 0.;
    ANISOTROPY[i] = false;
    TargetPhotonLookups[i] = NULL;
    ICLossLookups[i] = NULL;
    ICLossLookupAccs[i] = NULL;
    TargetPhotonAccs[i] = NULL;
    TargetPhotonAngularDistrs[i] = NULL;
    CosZetaLookups[i] = NULL;
    phiaccescs[i] = NULL;
    thetaaccescs[i] = NULL;
    phiaccesc_zetas[i] = NULL;
    thetaaccesc_zetas[i] = NULL;
  }
  TargetPhotonEdensSumIso = 0.;
  fUtils->Clear2DVector(TargetPhotonVectorSumAll);
  fUtils->Clear2DVector(TargetPhotonVectorSumIso);
//  TargetPhotonVectorSumIso.clear();
  LUMFLAG = false;
  INTEGRATEOVERGAMMAS = false;
  QUIETMODE = false;
  VERBOSEMODE = false;
  PiModel = 1;
  SynchModel = 0;
  integratorTolerance = 1.e-1;
  integratorKronrodRule = 2;
  n = 0.;
  SynchAngle = 90.;
  acc = gsl_interp_accel_alloc();
  acciso = gsl_interp_accel_alloc();
  ICLossLookupAccIso = gsl_interp_accel_alloc();
  ICLossLookupAccAll = gsl_interp_accel_alloc();
  loraccesc = gsl_interp_accel_alloc();
  edaccesc = gsl_interp_accel_alloc();
}

/**
 * Standard destructor.
 */
Radiation::~Radiation() {}

//FIXME make me a nice function!
void Radiation::Reset() {
  ParticleVector.clear();
  ElectronVector.clear();
  ProtonVector.clear();
  gsl_spline_free(ElectronLookup);
  gsl_spline_free(ProtonLookup);
  gsl_spline_free(TargetPhotonLookupSumIso);
  BField = 0.;
  n = 0.;
  return;
}

/*
 * remove all previously set IC target photons
 */
void Radiation::ClearTargetPhotons() {
    for(int i=-2;i<(int)RADFIELDS_MAX;i++) {
        ClearTargetPhotonField(i);
    }
    RADFIELD_COUNTER = 0;
    return;
}

/**
 * Calculates the differential photon emission at energy 'e' [erg]. This results
 *  in
 * - differential fluxes, fdiff* [(no. of photons)/(erg*s*cm^2)]
 * - differential photon rate, fdiff* [(no. of photons)/(erg*s)]
 *
 * The * stands for Brems(-strahlung), IC and pp (i.e. pi^0 decay) components.
 * These quantities are private members.
 *
 * It calls the 'DifferentialEmissionComponent'
 * function for the relevant radiation mechanisms:
 * - protons ('particletype'=1)
 *   + inelastic p-p scattering
 * - electrons ('particletype'=0)
 *   + IC emission
 *   + Bremsstrahlung
 *   + Synchrotron radiation
 *
 * This function calculates the appropriate (depending on the particle species)
 * radiation mechanism and their gamma-ray flux by calling
 * the method 'DifferentialEmissionComponent'.
 */
void Radiation::CalculateDifferentialGammaEmission(double e, int particletype) {
  ldiffbrems = fdiffbrems = ldiffsynch = fdiffsynch = 0.;
  ldiffic = fdiffic = ldiffpp = fdiffpp = 0.;

  void *p = NULL;
  if(!lumtoflux) {
    if (!distance) {
      if(LUMFLAG == false) {
        cout << "### Radiation::CalculateDifferentialGammaEmission: Distance to "
                "particles not specified -> Flux equals now the luminosity! ###"
             << endl;
      }
      LUMFLAG = true;
      lumtoflux = 1.;
    } else
      lumtoflux = 1. / (4. * pi * distance * distance);
  }
  if (particletype != 0. && particletype != 1) {
    cout << "### Radiation::CalculateDifferentialGammaEmission: Please provide "
            "proper particle spectrum and particle type!  ###" << endl;
    return;
  }
  if ((!particletype && !ElectronVector.size()) ||
      (particletype == 1 && !ProtonVector.size())) {
    cout << "### Radiation::CalculateDifferentialGammaEmission: No accelerated "
            "particles! Exiting... ###" << endl;
    return;
  } else if (!particletype) {
    ParticleVector = ElectronVector;
    if(n) {
      radiationMechanism = "Bremsstrahlung";
      ldiffbrems = DifferentialEmissionComponent(e, p);
      fdiffbrems = lumtoflux * ldiffbrems;
    }
    if(RADFIELD_COUNTER) {
        fdiffic = 0.;
        double ldiffic_sum = 0.;
        radiationMechanism = "InverseCompton";
        for(unsigned int i = 0;i<RADFIELDS_MAX;i++) {
//            std::cout<<i<<" "<<TargetPhotonAngularDistrs[i]<<" "<<FASTMODE_IC<<std::endl;
            if(TargetPhotonLookups[i]!=NULL) {
                if(FASTMODE_IC == true && TargetPhotonAngularDistrs[i] == NULL) 
                    continue;
                if(FASTMODE_IC == false && TargetPhotonAngularDistrs[i] != NULL) 
                    continue;
                SetICLookups(i);
                ldiffic = DifferentialEmissionComponent(e, p);
                fdiffics[i] = lumtoflux * ldiffic;
                ldiffic_sum += ldiffic;
            }
        }
        if(FASTMODE_IC == true && TargetPhotonVectorSumIso.size()) {
                SetICLookups(-1);
                ldiffic = DifferentialEmissionComponent(e, p);
                ldiffic_sum += ldiffic;
                IC_CALCULATED = true;
        }
        fdiffic = lumtoflux * ldiffic_sum;
    }
    if(BField) {
      radiationMechanism = "Synchrotron";
      ldiffsynch = DifferentialEmissionComponent(e, p);
      fdiffsynch = lumtoflux * ldiffsynch;
    }
  } else if (particletype == 1) {
    ParticleVector = ProtonVector;
    radiationMechanism = "ppEmission";
    ldiffpp = DifferentialEmissionComponent(e, p);
    fdiffpp = lumtoflux * ldiffpp;
  } else {
    cout << "### Radiation::CalculateDifferentialGammaEmission: WTF?! That is "
            "not possible!!" << endl;
  }
  return;
}

void Radiation::SetICLookups(int i) {
    
    if(i==-1) {
        vector<double> vec_null;
        TargetPhotonVectorCurrent = &TargetPhotonVectorSumIso;
        TargetPhotonLookupCurrent = &TargetPhotonLookupSumIso;
        TargetPhotonEdensCurrent = &TargetPhotonEdensSumIso;
        TargetAccCurrent = &acciso;

        ICLossVectorCurrent = &ICLossVectorSumIso;
        ICLossLookupCurrent = &ICLossLookupSumIso;
        ICLossLookupAccCurrent = &ICLossLookupAccIso;
        ANISOTROPY_CURRENT = false;
    }
    else {
        i = (int)i;
        TargetPhotonVectorCurrent = &TargetPhotonVectors[i];
        TargetPhotonLookupCurrent = &TargetPhotonLookups[i];
        TargetPhotonEdensCurrent = &TargetPhotonEdensities[i];

        TargetPhotonAngularDistrCurrent = &TargetPhotonAngularDistrs[i];
        CosZetaLookupCurrent = &CosZetaLookups[i];
        TargetAccCurrent = &TargetPhotonAccs[i];

        ICLossVectorCurrent = &ICLossVectors[i];
        ICLossLookupCurrent = &ICLossLookups[i];
        ICLossLookupAccCurrent = &ICLossLookupAccs[i];


        phiaccescCurrent = &phiaccescs[i];
        thetaaccescCurrent = &thetaaccescs[i];
        phiaccesc_zetaCurrent = &phiaccesc_zetas[i];
        thetaaccesc_zetaCurrent = &thetaaccesc_zetas[i];
        TargetPhotonAngularBoundsCurrent = &TargetPhotonAngularBounds[i];
        ANISOTROPY_CURRENT = ANISOTROPY[i];
    }
    return;
}

/**
 * Calculates the differential photon rate [(no. of photons)/(erg*s)]
 * at energy 'e' [erg] resulting from radiation mechanism
 * 'radiationMechanism' that has been specified before in
 * 'CalculateDifferentialGammaEmission' or 'CalculateIntegralGammaEmission'.
 */
double Radiation::DifferentialEmissionComponent(double e, void *par) {
  if (radiationMechanism.compare("Synchrotron") &&
      radiationMechanism.compare("Bremsstrahlung") &&
      radiationMechanism.compare("InverseCompton") &&
      radiationMechanism.compare("ppEmission")) {
    cout << "### Radiation::DifferentialEmissionComponent: no valid emission mechanism "
            "specified! Returning 0 value ... ###" << endl;
    return 0.;
  }
  if (!ParticleVector.size()) {
    cout << "### Radiation::DifferentialEmissionComponent: No accelerated particles! "
            "Exiting... ###" << endl;
    return 0.;
  }
  e = log10(e);
  double emax = log10(ParticleVector[ParticleVector.size() - 1][0]);
  double emin = log10(ParticleVector[0][0]);
  if (e > emax) return 0.;
  double egamma = e;
  if(e<emin) e = emin;
  fPointer IntFunc = NULL;
  if (!radiationMechanism.compare("Synchrotron")) {
    if (!BField) {
      if(!QUIETMODE) cout << "Radiation::DifferentialEmissionComponent: No "
                             "BField value set for Synchrotron radiation. "
                             "Returning zero value." << endl;
      return 0.;
    }
    if (!SynchModel) IntFunc = &Radiation::SynchEmissivity;
    else if (SynchModel == 1) IntFunc = &Radiation::SynchEmissivityExplicit;
    else {
      cout << "Radiation::DifferentialEmissionComponent: Specify valid "
              "Synchrotron emission model. 0 - random B-Field, 1 - regular "
              "B-Field with 90 degree electron-BField pitch angle. "
              "Returning zero value."
           << endl;
      return 0.;
    }
  } else if (!radiationMechanism.compare("Bremsstrahlung")) {
    if (!n) {
      if(!QUIETMODE) cout << "Radiation::DifferentialEmissionComponent: "
                             "No ambient density value set "
                             "for Bremsstrahlung. Returning zero value." << endl;
      return 0.;
    }
    IntFunc = &Radiation::BremsEmissivity;
  } else if (!radiationMechanism.compare("InverseCompton")) {
    
    IntFunc = &Radiation::ICEmissivityRadFieldIntegrated;
  } else if (!radiationMechanism.compare("ppEmission")) {
    if (!n) {
      if(!QUIETMODE) cout << "Radiation::DifferentialEmissionComponent:"
                             "No ambient density value set for"
                             "p-p scattering. Returning zero value." << endl;
      return 0.;
    }
    if(PiModel<0 || PiModel>3)  {
      cout << "Radiation::DifferentialEmissionComponent:"
              "Please specify valid p-p interaction model." << endl;
      cout << "Options are: " << endl;
      cout << "  0 - Geant4"  << endl;
      cout << "  1 - Pythia8" << endl;
      cout << "  2 - SIBYLL2" << endl;
      cout << "  3 - QGSJET I"  << endl;
      cout << "Set it via the SetPPEmissionModel(<OPTION>) method. " << endl;
      return 0.;
    }
    IntFunc = &Radiation::PPEmissivity;
  } else
    return 0.;
  gsl_interp_accel_reset(acc);
  double gammas = Integrate(IntFunc, &egamma, e, emax, 0.5*integratorTolerance,
                            integratorKronrodRule+1);
  if (std::isnan(gammas)) return 0.;
  if (gammas < 0.) return 0.;
  return gammas;
}


/********* RADIATION MECHANISMS ***********************************/

/* Inverse Compton part */

/**
 * Describes grey body with temperature 'temp' [K] and energy density
 * 'edens' [erg/cm^3] and returns differential photon
 * density [(no of photons)/(cm^3*erg)] at energy 'ephoton' [erg].
 */
double Radiation::GreyBody(double ephoton, double temp, double edens) {
  double val = 15. * edens * pow(pi, -4.) * pow(kb * temp, -4.) *
                       pow(ephoton, 2.) / (exp(ephoton / (kb * temp)) - 1.);
  if(fabs(log10(val))>30.) val = 0.;
  return val;
}

/**
 * IC emission from electrons integrated over the target photon population.
 * This method has the switch INTEGRATEOVERGAMMAS, which determines its
 * output. If INTEGRATEOVERGAMMAS is set to
 *
 * - TRUE: method returns a photon production rate [(no. of photons)/s] at
 *         energy 'egamma'. It is used to calculate the total loss rate due to
 *         IC emission for single electrons of energy 'eelectrons' by
 *         integrating this production rate from a (hardcoded) minimum energy of
 *         1.e-10 eV  up to energy 'eelectrons'.
 *         This is done in the method 'CreateICLossLookup'.
 *
 *
 * - FALSE: method returns a differential production rate of photons
 *          [(no. of photons)/(erg*s)] with energy 'egamma' from the
 *          differential number of electrons at energy 'eelectron', which is
 *          then used to calculate the total IC gamma-ray emission at photon
 *          energy egamma from the total electron distribution. The
 *          corresponding integration is performed from the minimum electron
 *          energy to egamma and is implemented in the method
 *          'DifferentialEmissionComponent'
 */
double Radiation::ICEmissivityRadFieldIntegrated(double x, void *par) {
  /* energy of scattering electron */
  double eelectron = 0.;
  /* energy of gamma ray */
  double egamma = 0.;
  /* number of gamma rays @ egamma */
  double icgammas = 0.;

  /* change mode of this function by switching par and x. */
  if (INTEGRATEOVERGAMMAS == true) {
    eelectron = pow(10.,*(double *)par);
    egamma = pow(10.,x);
  } else {
    eelectron = pow(10.,x);
    egamma = pow(10.,*(double *)par);
  }

  gsl_interp_accel_reset(acc);
  fPointer IntFunc = &Radiation::ICEmissivity;
  double integratorTolerance_IC = integratorTolerance;
  int kronrod = 6; 
  if(ANISOTROPY_CURRENT == true)  {
    if(ISOTROPIC_ELECTRONS){
        IntFunc = &Radiation::ICEmissivityAnisotropicIsotropicElectrons;
    }
    else {
        IntFunc = &Radiation::ICEmissivityAnisotropic;
    }
    integratorTolerance_IC = integratorTolerance*5;
    kronrod = 2;
  }

  double xpars[2] = {eelectron, egamma};

  /* detemine integration boundaries for the target photon energy from Eq. 2.50
     in Blumenthal&Gould (Reviews of Modern Physics, vol. 42, no. 2, 1970) */
  double lorentz = eelectron / m_e;
  double edash = egamma / (1. - egamma / (lorentz * m_e));
  double k = 1. / (4. * lorentz * lorentz);
  double boundmin, boundmax;
  vector<double> minmax = fUtils->GetVectorMinMax(*TargetPhotonVectorCurrent,0);//FIXME:Slow!
  double tph_min = pow(10.,minmax[0]);
  double tph_max = pow(10.,minmax[1]);
  (edash *k < tph_min) ? (boundmin = tph_min)
                                     : boundmin = edash * k;
  (edash > tph_max) ? (boundmax = tph_max)
                                  : boundmax = edash;

   boundmin = log10(boundmin);
   boundmax = log10(boundmax);

  if (k > 0.1 || boundmin >= boundmax) return 0.;
  icgammas = ln10*pow(10.,x) * Integrate(IntFunc, xpars, boundmin, boundmax, integratorTolerance_IC,kronrod);
  if (std::isnan(icgammas) || std::isinf(icgammas)) return 0.;
  if (INTEGRATEOVERGAMMAS == true)
    return icgammas * egamma;
  else {
    gsl_interp_accel_reset(acc);
    double elnumber = fUtils->EvalSpline(log10(eelectron),ElectronLookup,acc
                                  ,__func__,__LINE__);
    icgammas *= pow(10., elnumber);
    return icgammas;
  }
}

/**
 * emissivity from IC scattering. Taken from Blumenthal and Gould 1970,
 * Eq(2.48).
 */
double Radiation::ICEmissivity(double x, void *par) {
  double ephoton = pow(10.,x);  ///< energy of the target photon
  double *p = (double *)par;
  double lorentz = p[0] / m_e;
  double egamma = p[1];  ///< energy of the resulting gamma photon
  double e1 = egamma / (lorentz * m_e);  ///< gamma-ray energy in units of the
                                         ///electron energy
  double gamma = 4. * ephoton * lorentz / m_e;  ///< parameter that describes
                                                ///the regime of the scattering
                                                ///process. Small value: Thomson
                                                ///regime, large: KN-regime
  double q = e1 / (gamma * (1. - e1));  ///< yet another parameter telling us
                                        ///the scattering domain
  double f = 1./(4.*lorentz*lorentz);
  if(f > 0.1) return 0.;
  if(q>1. ||  q<f) return 0.;
  /// Eq(2.48):
  double bracket = 2. * q * log(q) + (1. + 2. * q) * (1. - q)
                   + 0.5 * (1. - q) * gamma * q * gamma * q / (1. + gamma * q);
  double targetphotons = fUtils->EvalSpline(x,
                                            *TargetPhotonLookupCurrent,
                                            *TargetAccCurrent,__func__,__LINE__);
  double integrand = 2. * pi * e_radius * e_radius * m_e * c_speed / lorentz *
                     pow(10., targetphotons) / ephoton * bracket;

  integrand /= lorentz*m_e;
  integrand *= ln10*ephoton;
  return integrand;
}

double Radiation::ICEmissivityWrapper(double e_ph, double e_e, double e_g) {

    double p[2] = {e_e,e_g};

    return ICEmissivity(e_ph,p);
}

double Radiation::ICEmissivityAnisotropicWrapper(double e_ph, double e_e, double e_g) {

    double p[2] = {e_e,e_g};

    return ICEmissivityAnisotropic(e_ph,p);
}

/**
 * IC emissivity in anisotropic radiation field. From Moskalenko & Strong,
 * Astrophys.J. 528 (2000) 357-367.
 */
double Radiation::ICEmissivityAnisotropic(double x, void *par) {
    double ephoton = pow(10.,x);  ///< energy of the target photon
    double *p = (double *)par;
    double lorentz = p[0] / m_e;
    double egamma = p[1];  ///< energy of the resulting gamma photon

    double beta = sqrt(1. - 1. / (lorentz*lorentz));
    double Q = 0.; double F = 0.; double cos_zeta = 0.;
    double cos_zeta_min = egamma/(2.*ephoton*lorentz*(lorentz-egamma/c_speed))-1.; 

//    double cos_zeta_min2 = (egamma/m_e)/(2.*ephoton/m_e*lorentz*(lorentz-egamma/m_e))-1.; 
    double zeta_min = acos(cos_zeta_min);
    double integral = 0.;
    double cos_kappa = cos(theta_e); double sin_kappa = sin(theta_e);   
    double pi_min_ka = pi - theta_e; 

    double phi_min = (*TargetPhotonAngularBoundsCurrent)[0];
    double phi_max = (*TargetPhotonAngularBoundsCurrent)[1];
    double theta_min = (*TargetPhotonAngularBoundsCurrent)[2];
    double theta_max = (*TargetPhotonAngularBoundsCurrent)[3];
    for (double phi = phi_min; phi <= phi_max; phi += d_phi) {
//        std::cout<<phi<<"("<<phi_min<<" - "<<phi_max<<")"<<std::endl;
        for (double theta = theta_min; theta <= theta_max; theta += d_theta) {
//            cos_zeta = interp2d_spline_eval(*CosZetaLookupCurrent, phi, theta, 
//                                     *phiaccesc_zetaCurrent,*thetaaccesc_zetaCurrent);

            cos_zeta = -cos_kappa *cos(theta) + sin_kappa * sin(theta) * cos(phi - phi_e);
////            std::cout<<"cosz,cosz_min = "<<cos_zeta<<" "<<cos_zeta_min<<" "<<cos_zeta_min2<<std::endl;
//            cos_zeta  = cos(phi) * cos(phi_e) * sin(theta) * sin(theta_e);
//            cos_zeta += sin(phi) * sin(phi_e) * sin(theta) * sin(theta_e);
//            cos_zeta += cos(theta) * cos(theta_e);
//            std::cout<<"cos_zeta1,cos_zeta "<<cos_zeta1<<" "<<cos_zeta<<std::endl;
//            cos_zeta = cos(phi) * sin(theta);
//            std::cout<<phi<<" "<<theta<<" "<<Q<<" "<<F<<std::endl;
            if (cos_zeta < cos_zeta_min) continue;

            // Kinematic limits from sec. 2.2.1.
            double upsilon = cos_zeta_min+cos_kappa*cos(theta) / (sin_kappa*sin(theta));
            if (std::isinf(upsilon) || std::isnan(upsilon)) continue;
            if(fabs(theta - pi_min_ka) > zeta_min) { continue;}
            if (fabs(upsilon) < 1.) {
                
                if(fabs(phi - phi_e) > acos(upsilon)) continue;
            }            
            else {
                if(upsilon > 0.) continue;
                if(upsilon < 0. && fabs(phi-phi_e) > pi) continue;
            }
        
            Q = interp2d_spline_eval(*TargetPhotonAngularDistrCurrent, 
                                     phi, theta, *phiaccescCurrent,*thetaaccescCurrent);
            double eph_d = ephoton/m_e * lorentz * (1. + beta*cos_zeta);    // M&S equ. 9
            if (egamma/m_e > 2. * lorentz * eph_d / (1. + 2.*eph_d)) continue; // see M&S equ.9

            F = ICAnisotropicAuxFunc(phi,theta,ephoton,egamma,lorentz,beta,cos_zeta);
            integral += sin(theta) * d_phi * d_theta * F * Q;
//                std::cout<<"Q,F,int,cz,dph,dth: "<<Q<<","<<F<<","<<integral<<","<<cos_zeta<<","<<d_phi<<","<<d_theta<<std::endl;
            
        }
    }

    double targetphotons = fUtils->EvalSpline(x,
                                              *TargetPhotonLookupCurrent,
                                              *TargetAccCurrent,__func__,__LINE__);
    double integrand = pi * e_radius * e_radius * c_speed;
    integrand /= ephoton * (lorentz-egamma/m_e) * (lorentz-egamma/m_e);
    integrand *= integral * pow(10.,targetphotons);
    integrand *= ln10 * ephoton;
//    std::cout<<"-------------------------------------"<<std::endl;
    return integrand;
}







/*----------------------------------------------------------------------------------------------------
 *      This Block calculates the Anisotropic IC emission for isotropic electrons
 *      with equation 20 from Aharonian & Atoyan 1981
 *--------------------------------------------------------------------------------------------------*/

double Radiation::ICEmissivityAnisotropicIsotropicElectrons(double x, void *par){
    double ephoton = pow(10.,x);  ///< energy of the target photon
    double *p = (double *)par;
    double eelectron = p[0];
    double egamma = p[1];  ///< energy of the resulting gamma photon
    
    double phi_min = (*TargetPhotonAngularBoundsCurrent)[0];
    double phi_max = (*TargetPhotonAngularBoundsCurrent)[1];
    double theta_min = (*TargetPhotonAngularBoundsCurrent)[2];
    double theta_max = (*TargetPhotonAngularBoundsCurrent)[3];
    
    fPointer theta_int_function = &Radiation::ICEmissivityAnisotropicIsotropicElectronsSecondIntegral;
    
    
    
    double integralinput[6];
    
    integralinput[0] = ephoton;
    integralinput[1] = egamma;
    integralinput[2] = eelectron;
    integralinput[3] = phi_min;
    integralinput[4] = phi_max;
    
    
    double result;
    
    result = Integrate(theta_int_function, integralinput, theta_min, theta_max,integratorTolerance*5.,integratorKronrodRule);
    
    double targetphotons = fUtils->EvalSpline(x,
                                              *TargetPhotonLookupCurrent,
                                              *TargetAccCurrent,__func__,__LINE__);

    result = e_radius*e_radius/(2.*ephoton*eelectron*eelectron)*result*pow(10.,targetphotons)*ln10*ephoton;

    
    result *= m_e*m_e*c_speed;      // To Get the right units
    result *= 4.*pi;                // Have to multiply with 4*pi, because otherwise the flux at a specific
                                    // distance is wrong (see Radiation::CalculateDifferentialGammaEmission)
    return result;
}


double Radiation::ICEmissivityAnisotropicIsotropicElectronsSecondIntegral(double x, void *par){
    double theta = x;
    double *p = (double *)par;
    
    double phi_min = p[3];
    double phi_max = p[4];
    
    double integralinput[4];
    
    fPointer phi_int_function = &Radiation::ICEmissivityAnisotropicIsotropicElectronsFirstIntegral;
    
    integralinput[0] = p[0];
    integralinput[1] = p[1];
    integralinput[2] = p[2];
    integralinput[3] = theta;
    
    double result;
    result = Integrate(phi_int_function, integralinput, phi_min, phi_max,integratorTolerance*5., integratorKronrodRule);
    
    return (result*sin(theta));
}

double Radiation::ICEmissivityAnisotropicIsotropicElectronsFirstIntegral(double x, void *par){
    double phi = x;
    double *p = (double *)par;
    
    double ephoton = p[0];
    double egamma = p[1];
    double eelectron = p[2];
    double theta = p[3];
    
    
    double cos_angle = sin(theta)*cos(phi);

    
    double z = egamma/eelectron;
    
    double b_theta = 2.*(1.- cos_angle)*ephoton*eelectron/(m_e*m_e);
    
    
    if(egamma/m_e > (b_theta/(1.+b_theta) * eelectron/m_e)) return 0.0;
    if((10.*ephoton) > egamma) return 0.0;
    //if( eelectron < 10.*m_e ) return 0.0;
    
    double oneminusz = 1.-z;
    
    double Q = interp2d_spline_eval(*TargetPhotonAngularDistrCurrent, 
                                     phi, theta, *phiaccescCurrent,*thetaaccescCurrent);

    double result = 1. + z*z/(2.*oneminusz) -
                    2.*z/(b_theta*oneminusz) + 2.*z*z/(b_theta*b_theta*oneminusz*oneminusz);
    

                    
    return (result*Q);
    
    
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------








void Radiation::FillCosZetaLookup(int i) {
    vector <vector<double> > v;

    double phi_min = TargetPhotonAngularBounds[i][0];
    double phi_max = TargetPhotonAngularBounds[i][1];
    double theta_min = TargetPhotonAngularBounds[i][2];
    double theta_max = TargetPhotonAngularBounds[i][3];

    double ph_ma = phi_max+0.1*fabs(phi_max); double ph_mi = phi_min-0.1*fabs(phi_min);
    double th_ma = theta_max+0.1*fabs(theta_max); double th_mi = theta_min-0.1*fabs(theta_min);
    double d_ph = 0.25*d_phi;
    double d_th = 0.25*d_theta;

    for (double theta = th_mi; theta <= th_ma; theta += d_th) {
        for (double phi = ph_mi; phi <= ph_ma; phi += d_ph) {
            double sin_phi = sin(phi); double cos_phi = cos(phi);
            double sin_theta = sin(theta); double cos_theta = cos(theta);

            double cos_zeta = (cos_phi_e * cos_phi + sin_phi_e * sin_phi);
            cos_zeta = sin_theta_e * sin_theta * cos_zeta + cos_theta_e * cos_theta;

            fUtils->TwoDVectorPushBack(phi,theta,cos_zeta,v);
        }
    }
    double a,b,c,d;
    CosZetaLookups[i] =
      fUtils->TwoDsplineFromTwoDVector(v,a,b,c,d);  
    return;
}


/** 
 * Moskalenko and Strong 1999 Eq. 8 integrand expression in brackets
 **/
double Radiation::ICAnisotropicAuxFunc(double phi_p, double theta_p, 
                                       double ephoton, double egamma,
                                       double lorentz, double beta, 
                                       double cos_zeta) {

    ephoton /= m_e;
    egamma  /= m_e;
    if (egamma > lorentz) return 0.;

    double one_lor = 1./lorentz;

    double eph_d = ephoton * lorentz * (1. + beta*cos_zeta);
    double one_eph_d = 1. / eph_d;
    double eg_lo = egamma *one_lor;
    
    double val = 0.;
    if (egamma <= 2. * lorentz * eph_d / (1. + 2.*eph_d)) {
        val = 2.-2.*eg_lo*(one_eph_d+2.);
        val += eg_lo*eg_lo*(one_eph_d*one_eph_d+2.*one_eph_d+3.);
        val -= eg_lo*eg_lo*eg_lo;
    }
    return val;
}



void Radiation::CreateICLossLookup(int bins) {
    if (FASTMODE_IC_LOSSLOOK == true && RADFIELD_COUNTER) {
        SumTargetFieldsIsotropic();
        SumTargetFieldsAll();
    }
    for(unsigned int i = 0;i<RADFIELDS_MAX;i++) {
        if(ANISOTROPY[i] == false || TargetPhotonLookups[i]==NULL) 
            continue;
        CreateICLossLookupIndividual(i,bins);
    }
    if(TargetPhotonVectorSumIso.size()) CreateICLossLookupIndividual(-1,bins);

    vector<double> minmax;
    double min,max;
    if(ICLossVectorSumIso.size()) {
        minmax = fUtils->GetVectorMinMax(ICLossVectorSumIso,0);
        min = minmax[0];
        max = minmax[1];
    }
    else {
        min = 1e100;
        max = -1e100;
    }
    minmax.clear();
    for(unsigned int i=0;i<RADFIELDS_MAX;i++) {
        if(!ICLossVectors[i].size()) continue;
        minmax = fUtils->GetVectorMinMax(ICLossVectors[i],0);
        if(minmax[0] < min) min = minmax[0];
        if(minmax[1] > max) max = minmax[1];
        minmax.clear();
    }


    fUtils->Clear2DVector(ICLossVectorSumAll);
    ICLossLookupSumAll = NULL;

    double logestep = ( max-min ) / bins;
    for(double loge=min;loge<max;loge+=logestep) {
        double sum = 0.;
        if(ICLossVectorSumIso.size()) {
            minmax = fUtils->GetVectorMinMax(ICLossVectorSumIso,0);
            if(loge > minmax[0] && loge < minmax[1]) {
                sum += pow(10.,fUtils->EvalSpline(loge,ICLossLookupSumIso,
                                      ICLossLookupAccIso,__func__,__LINE__));
            }
            minmax.clear();
        }
        for(unsigned int i=0;i<RADFIELDS_MAX;i++) {
            if(ANISOTROPY[i] == true ) {
                minmax = fUtils->GetVectorMinMax(ICLossVectors[i],0);
                if(loge > minmax[0] && loge < minmax[1]) {
                    sum += pow(10.,fUtils->EvalSpline(loge,ICLossLookups[i],
                                      ICLossLookupAccs[i],__func__,__LINE__));
                }
                minmax.clear();
            }
        }
        if(sum) fUtils->TwoDVectorPushBack(loge,log10(sum),ICLossVectorSumAll);
    }

    int size = (int)ICLossVectorSumAll.size();
    double e[size]; double l[size];
    for (int i=0;i<size;i++) {
        e[i] = ICLossVectorSumAll[i][0];
        l[i] = ICLossVectorSumAll[i][1];  
    }
    ICLossLookupSumAll = gsl_spline_alloc(gsl_interp_linear, size);
    gsl_spline_init(ICLossLookupSumAll, e, l, size);

    return;
}

/** return a lookup table holding the differential electron energy loss rate due
 * to inverse-Compton
 *  scattering. The format of the lookup is: { Energy(erg) - Energy Loss Rate
 * from IC scattering(erg/s) }
 */
void Radiation::CreateICLossLookupIndividual(int i, int bins) {

  if(i<-1 || i >= (int)RADFIELDS_MAX) {
    cout<<"Radiation::CreateICLossLookupIndividual: field "<<i<<" not there. Exiting"<<endl;
    return;
  }

  if( (i>=0 && i< (int)RADFIELDS_MAX && ICLossVectors[i].size()) ||
      (i==-1 && ICLossVectorSumIso.size()) ) {
    
    cout<<"Radiation::CreateICLossLookupIndividual: loss lookup for field "<<i<<" already "
          "computed. Exiting."<<endl;
    return;
  }
  SetICLookups(i);
  if(!(*TargetPhotonVectorCurrent).size()) {
    cout << "Radiation::CreateICLossLookupIndividual: No target photons! Exiting." <<endl;
    return;
  }

  fUtils->Clear2DVector(*ICLossVectorCurrent);

  double av_cos_xi = 0.; double area = 0.;
  if (ANISOTROPY_CURRENT == true) {
      phi_min = (*TargetPhotonAngularBoundsCurrent)[0];
      phi_max = (*TargetPhotonAngularBoundsCurrent)[1];
      theta_min = (*TargetPhotonAngularBoundsCurrent)[2];
      theta_max = (*TargetPhotonAngularBoundsCurrent)[3];
      for (double phi = phi_min; phi <= phi_max; phi += d_phi) {
          for (double theta = theta_min; theta <= theta_max; theta += d_theta) {
              double zeta = acos(interp2d_spline_eval(*CosZetaLookupCurrent, phi, theta, 
                                     *phiaccesc_zetaCurrent,*thetaaccesc_zetaCurrent));
              double Q = interp2d_spline_eval(*TargetPhotonAngularDistrCurrent,
                                     phi, theta,*phiaccescCurrent,*thetaaccescCurrent);
              double xi = pi - zeta;
              double cos_xi = cos(xi);
              av_cos_xi += Q*(1.-cos_xi)*(1.-cos_xi)*d_theta*d_phi*sin(theta);
              area += d_theta * d_phi * sin(theta);
          }
      }
  }
  else av_cos_xi = 4. / 3.;

 
  INTEGRATEOVERGAMMAS = true;
  /* lower integration boundary over emitted (i.e. 'loss-') IC photons */
  double EGammaMin = 1.e-25 * TeV_to_erg;
  /* Upper integration boundary over emitted (i.e. 'loss-') IC photons */
  double EGammaMax = 1.e8 * TeV_to_erg;
  /* lower integration boundary for the IC loss lookup (i.e. here simply the
   * electron rest mass) */
  double logemin = log10(1.e-8 * m_e);
  /* upper integration boundary for the IC loss lookup */
  double logemax = log10(EGammaMax);
  double logestep = (double)(logemax - logemin) / bins;
  int tt = 1;
  int ii = 1;
  if (!QUIETMODE) {
    cout << ">> CALCULATING IC LOSS LOOKUP " << endl;
  }

  double phEmax = pow(10.,(*TargetPhotonVectorCurrent)[(*TargetPhotonVectorCurrent).size()-1][0]);
  /* gamma value that indicates Thomson regime (see Blumenthal&Gould) */
  double GammaLow = 1.e-1;
  /* transition energy to Thomson regime. Losses are then simply Edot~E*E*edens*/
  double Etrans = m_e * m_e * GammaLow / phEmax;

  for (double loge = logemin; loge < logemax; loge += logestep) {
    if ((double)ii / bins > 0.0001 * tt && QUIETMODE == false) {
      cout << "\r";
      cout << "    " << (int)(100. * ii / bins) - 1 << "\% done" << std::flush;
      tt++;
    }
    double LossRate = 0.;
    double Eelectron = pow(10., loge);
    if (Eelectron>Etrans) {
      fPointer IntFunc = &Radiation::ICEmissivityRadFieldIntegrated;
      
      LossRate =
          Integrate(IntFunc, &loge, log10(EGammaMin), log10(Eelectron),
                    integratorTolerance,integratorKronrodRule);
    } else {
      double gamma = (Eelectron + m_e) / m_e;
      LossRate =
          av_cos_xi * sigma_T * c_speed * (*TargetPhotonEdensCurrent) * gamma * gamma;
    }
    if (std::isnan(LossRate)) {
      cout << __func__ << ",l." << __LINE__ <<": LossRate is nan! Exiting."
           << endl;
      exit(1);
    }
    if(LossRate) fUtils->TwoDVectorPushBack(log10(Eelectron),log10(LossRate),*ICLossVectorCurrent);
    
    ii++;
  }
  
  unsigned int size = (*ICLossVectorCurrent).size();
  double e[size];
  double l[size];
  for(unsigned int g=0;g<size;g++){
     e[g] = (*ICLossVectorCurrent)[g][0];  
     l[g] = (*ICLossVectorCurrent)[g][1];       
  }


  *ICLossLookupCurrent = gsl_spline_alloc(gsl_interp_linear, size);
  gsl_spline_init(*ICLossLookupCurrent, e, l, size);
  *ICLossLookupAccCurrent = gsl_interp_accel_alloc();
  INTEGRATEOVERGAMMAS = false;
  if (QUIETMODE == false) {
    cout << endl;
    cout << "    -> DONE!   " << endl;
    cout << endl;
    cout << ">> COMPUTATION OF IC LOSS LOOKUP COMPLETE " << endl;
    cout << endl;
  }
  return;
}

vector<vector<double> > Radiation::GetICLossLookup(int i) {
    vector< vector<double> > v;
    if(!RADFIELD_COUNTER) return v;
    if(i>=0 && i<(int)RADFIELDS_MAX) {
        if(TargetPhotonLookups[i] == NULL) {
            cout<<"Radiation::GetICLossLookup: Field "<<i<<" not set. Returning "
                  "empty vector." <<endl;
            return v;
        }
        if(!ICLossVectors[(unsigned int)i].size()) {
            CreateICLossLookupIndividual(i);
        }
        for(unsigned int j=0;j<ICLossVectors[(unsigned int)i].size();j++) {
            double E = ICLossVectors[(unsigned int)i][j][0];
            double L = fUtils->EvalSpline(E, ICLossLookups[(unsigned int)i],
                                             ICLossLookupAccs[(unsigned int)i],
                                              __func__,__LINE__);
            fUtils->TwoDVectorPushBack(pow(10.,E),pow(10.,L),v);
        }
        return v;
    }
    else if (i==-1) {
        if(!ICLossVectorSumAll.size()) {
//            cout<< "Radiation::GetICLossLookup: Lookup not calculated yet. "
//                   "Please run Radiation::CreateICLossLookup first. Returning "
//                   "empty vector."<<endl;
            CreateICLossLookup();
//            return v;
        }
//        else {
            for(unsigned int j=0;j<ICLossVectorSumAll.size();j++) {
                
                double E = ICLossVectorSumAll[j][0];
                double L = fUtils->EvalSpline(E, ICLossLookupSumAll,ICLossLookupAccAll,
                                              __func__,__LINE__);
                
                fUtils->TwoDVectorPushBack(pow(10.,E),pow(10.,L),v);
            }
            return v;
//        }
    }
    else {
       cout<<"Radiation::GetICLossLookup: Index "<<i<<" not valid. Returning "
             "empty vector." << endl;
        return v;
    }
}



/* end of Inverse Compton part */

/* SYNCHROTRON PART */
/**
 * modified bessel functions
 */
double Radiation::K(double nu, double x) {
  if (x <= 0. || x > 700.)
    return 0.;
  else
    return gsl_sf_bessel_Knu(nu, x);
}

double Radiation::K_53(double x, void *par) {
//  double pp = *(double *)par;
  double K_4 = K(5. / 3., x);
  return K_4;
}

/**
 * emissivity of synchrotron radiation
 * adapted from galprop!ghisellini svensson 1988 'the synchrotron boiler'
 */
double Radiation::SynchEmissivity(double x, void *par) {
  /* frequency of emmited synchr. radiation */
  double nu = pow(10.,*(double *)par) / hp;
  /* electron energy */
  double eElectron = pow(10.,x);
  /* lorentz-factor of electrons */
  double gamma = eElectron / m_e;
  double nu_b = el_charge * BField * c_speed * pow(2. * pi * m_e, -1.);
  double j = nu / (3. * nu_b * pow(gamma, 2.));
  /* bessel fct. K_1/3 */
  double K_1 = K(1. / 3., j);
  /* bessel fct. K_4/3 */
  double K_4 = K(4. / 3., j);
  double value = 0.;
  if (nu < nu_b) {
    value = 0.;
  } else {
    double electrons = fUtils->EvalSpline(x,ElectronLookup,
                             acc,__func__,__LINE__);
    value = 4. * pi * sqrt(3.) * el_charge * el_charge * nu_b /
            (hp * hp * nu * c_speed);
    value *= pow(10., electrons) * pow(j, 2.);
    value *= (K_4 * K_1 - (3. / 5.) * j * (K_4 + K_1) * (K_4 - K_1));
  }
  return value * ln10 * eElectron;
}

/**
 * emissivity of synchrotron radiation
 * following Blumenthal&Gould, Eqs 4.44 and 4.48 !!!CHECK THAT!!!
 */
double Radiation::SynchEmissivityExplicit(double e, void *par) {

  double eElectron = pow(10.,e);
  double gamma = (eElectron + m_e) / m_e;
  double nu = pow(10.,*(double *)par) / hp;
  double geometry = sin(pi * SynchAngle / 180.);
  double norm = sqrt(3.) * geometry * el_charge * el_charge * el_charge * BField / m_e;
  double nu_c =
      3. * geometry * el_charge * BField * gamma * gamma * c_speed / (4. * pi * m_e);
  double x = nu / nu_c;
  fPointer IntFunc = &Radiation::K_53;
  double *v = NULL;
  double F = x * Integrate(IntFunc, v, x, 1.e3 * x, integratorTolerance,
                           integratorKronrodRule);
  if (!F) return 0.;
  double electrons = fUtils->EvalSpline(e,ElectronLookup,
                                        acc,__func__,__LINE__);
  double val = norm * F * pow(10., electrons) / (hp * hp * nu) / geometry / geometry;

  return val * ln10 * eElectron;
}

/* End of the Synchrotron part */

/* ---       BREMSSTRAHLUNG   --- */
/** emissivity of Bremsstrahlung,
 * proton-electron as well as electron-electron
 * From Baring 1999, ApJ, 513, 311-338
 */
double Radiation::BremsEmissivity(double x, void *par) {
  /* initial electron energy */
  double EI = pow(10.,x);
  /* bremsstrahlung photon energy */
  double EP = pow(10.,*(double *)par);
  /* threshold put by hand */
  if (EP < 1.e-8 * EI) return 0.;
  /* kinematic threshold */
  if (EI - EP <= m_e) return 0.;
  /* electron lorentz factor */
  double g = EI / m_e;
  /* electron velocity */
  double b = sqrt(1. - 1. / (g * g));
  /* photon energy in units of electron rest mass */
  double e = EP / m_e;
  /* equation (A1) */
  double sigma_e = (sigma1(g, e) + sigma2(g, e)) * A(g, e);
  /* 1.4 correction factor for 10% Helium */
  double S = n * 1.4;
  /* emissivity assuming a fully ionised plasma (n_e = S), Eq. (27) */
  double N;
  if (EI < 2.e-3 * GeV_to_erg)
    N = c_speed * b * S * (sigma1(g, e) + sigmaNR(g, e));
  else
    N = c_speed * b * S * (sigma1(g, e) + sigma_e);
  double electrons = fUtils->EvalSpline(log10(EI),ElectronLookup,
                                        acc,__func__,__LINE__);
  return N * pow(10., electrons) / m_e * ln10 * EI;
}

/** equation (A4)
  */
double Radiation::A(double g, double e) {
  return 1. - (8. / 3.) * pow(g - 1., 0.2) / (g + 1.) * pow(e / g, 0.3333);
}

/** equation (A2)
 */
double Radiation::sigma1(double g, double e) {
  return (4. * e_radius * e_radius * fineStructConst / e) *
         (1. + (0.3333 * -e / g) * (1. - e / g)) *
         (log(2. * g * (g - e) / e) - 0.5);
}

/** equation (A3)
 */
double Radiation::sigma2(double g, double e) {
  double k = e_radius * e_radius * fineStructConst / (3. * e);
  if (e <= .5) {
    return k * (16. * (1. - e + e * e) * log(g / e) - 1. / (e * e) + 3. / e -
                4. + 4. * e - 8. * e * e -
                2. * (1. - 2. * e) * log(1. - 2. * e) *
                    (1. / (4. * e * e * e) - 1. / (2. * e * e) + 3 / e - 2. +
                     4. * e));
  } else {
    return k * (2. / e) * ((4. - 1. / e + 1. / (4. * e * e)) * log(2. * g) -
                           2. + 2. / e - 5. / (8. * e * e));
  }
}

/** equation (A5), non-relativistic Bremsstrahlung
 */
double Radiation::sigmaNR(double g, double e) {
  if (e <= 0. || e >= 0.25 * (g * g - 1.)) return 0.;
  double sig = (4. * e_radius * e_radius * fineStructConst / (15. * e)) *
               Fbr(4. * e / (g * g - 1.), g);
  return sig;
}

/** equations (A6,A7)
 */
double Radiation::Fbr(double x, double g) {
  double B = 1. + 0.5 * (g * g - 1.);
  double b = sqrt(1. - 1. / (g * g));
  double C = 10. * x * g * b * (2. + g * b) / (1. + x * x * (g * g - 1.));
  double F = B * (17. - 3. * x * x / ((2. - x) * (2. - x)) - C) * sqrt(1. - x);
  F += (12. * (2. - x) - 7. * x * x / (2. - x) -
        3. * x * x * x * x / ((2. - x) * (2. - x) * (2. - x))) *
       log((1. + sqrt(1. - x)) / sqrt(x));
  return F;
}

/* end Bremsstrahlung part */

/* ---      pi0 decay     --- */

double Radiation::PPEmissivity(double x, void *par) {
  /* proton energy */
  double EP = pow(10.,x);
  /* pi0 decay photon energy */
  double Eg = pow(10.,*(double *)par);
  if (EP <= m_p) return 0.;
  double Tp = sqrt(EP * EP - m_p * m_p);
  double Tpth = 0.2797 * GeV_to_erg;
  if (Tp <= Tpth) return 0.;
  if (Eg <= GetMinimumGammaEnergy(Tp)) return 0.;
  if (Eg >= GetMaximumGammaEnergy(Tp)) return 0.;
  double N = DiffPPXSection(Tp, Eg);
  double logprotons = fUtils->EvalSpline(log10(EP),ProtonLookup,
                                      acc,__func__,__LINE__);
  return c_speed * n * N * pow(10., logprotons) * ln10 * EP;
}

/** differential cross section following Kafexhiu 2014
 *  (Eq. 8)
 */
double Radiation::DiffPPXSection(double Tp, double Eg) {
  double dsigmadE =
      NuclearEnhancementFactor(Tp) * Amax(Tp) * F(Tp, Eg) / GeV_to_erg;
  return dsigmadE;
}

/** inclusive pi0 production cross section over the full energy range.
 *  following Kafexhiu 2014, Paragragh II B 4
 */
double Radiation::InclusivePPXSection(double Tp) {
  double InclusiveXSection;
  if (Tp < 2. * GeV_to_erg)
    InclusiveXSection = SigmaOnePi(Tp) + SigmaTwoPi(Tp);
  else
    InclusiveXSection = InelasticPPXSectionKaf(Tp) * MeanMultiplicity(Tp);
  return InclusiveXSection;
}

/** inelastic proton-proton cross-section following
 *  Kafexhiu 2014 (Eq. 1)
 */
double Radiation::InelasticPPXSectionKaf(double Tp) {
  double Tpth = 0.2797 * GeV_to_erg;
  double L = Tp / Tpth;
  double logL = log(L);
  double sigma =
      (30.7 - 0.96 * logL + 0.18 * logL * logL) * pow(1. - pow(L, -1.9), 3.);
  return sigma * 1.e-27;
}

/** mean multiplicity, Eqs.6&7 in Kafexhiu 2014
 */
double Radiation::MeanMultiplicity(double Tp) {
  Tp /= GeV_to_erg;
  double meanMult, a1, a2, a3, a4, a5, TpTrans;
  meanMult = a1 = a2 = a3 = a4 = a5 = TpTrans = 0.;
  if (PiModel == 0)
    TpTrans = 0.;
  else if (PiModel == 1)
    TpTrans = 50.;
  else if (PiModel == 2 || PiModel == 3)
    TpTrans = 100.;
  else
    cout << "Radiation::MeanMultiplicity: WTF, specify supported "
            "parameterisation!" << endl;
  /* Eq. 6 */
  if (Tp < 5.) {
    double Tpth = 0.2797;
    double Q = (Tp - Tpth) / (m_p / GeV_to_erg);
    meanMult = -6.e-3 + 0.237 * Q - 0.023 * Q * Q;
  } else {
    double Xi = (Tp - 3.) / (m_p / GeV_to_erg);
    if (Tp < TpTrans) {
      double PiModelprev = PiModel;
      PiModel = 0;
      GetAParams(Tp, a1, a2, a3, a4, a5);
      PiModel = PiModelprev;
    } else
      GetAParams(Tp, a1, a2, a3, a4, a5);
    meanMult = a1 * pow(Xi, a4) * (1. + exp(-a2 * pow(Xi, a5))) *
               (1. - exp(-a3 * pow(Xi, 0.25)));
  }
  return meanMult;
}

/** Eq. 2 Kafexhiu 2014, needed for lowest energies
 */
double Radiation::SigmaOnePi(double Tp) {
  double sigma0 = 7.66e-3;
  double s = 2. * m_p * (Tp + 2. * m_p);
  double Mres = 1.1883 * GeV_to_erg;
  double Gres = 0.2264 * GeV_to_erg;
  double g = sqrt(Mres * Mres * (Mres * Mres + Gres * Gres));
  double kk = s - m_pi * m_pi - 4. * m_p * m_p;
  /* Eq. 3 */
  double eta =
      sqrt(kk * kk - 16. * m_pi * m_pi * m_p * m_p) / (2. * m_pi * sqrt(s));
  double K = sqrt(8.) * Mres * Gres * g / (pi * sqrt(Mres * Mres + g));
  double ll = (sqrt(s) - m_p) * (sqrt(s) - m_p) - Mres * Mres;
  /* Eq. 4 */
  double fBW = m_p * K / (ll * ll + Mres * Mres * Gres * Gres);
  double SigmaOnePi = sigma0 * pow(eta, 1.95) *
                      (1. + eta + eta * eta * eta * eta * eta) * pow(fBW, 1.86);
  return SigmaOnePi * 1.e-27;
}

/** Eq. 2 Kafexhiu 2014, needed for lowest energies
 */
double Radiation::SigmaTwoPi(double Tp) {
  double SigmaTwoPi = 5.7 / (1. + exp(-9.3 * (Tp / GeV_to_erg - 1.4)));
  return SigmaTwoPi * 1.e-27;
}

/** X-section normalisation following Eq. 12 in Kafexhiu 2014
 */
double Radiation::Amax(double Tp) {
  double b0, b1, b2, b3, thetap, logthetap, amax;
  if (Tp < 1. * GeV_to_erg) {
    b0 = 5.9;
    amax = b0 * InclusivePPXSection(Tp) / Epilabmax(Tp);
  } else {
    GetBParams(Tp, b1, b2, b3);
    thetap = Tp / m_p;
    logthetap = log(thetap);
    amax = b1 * pow(thetap, -b2) * exp(b3 * logthetap * logthetap) *
           InclusivePPXSection(Tp) / m_p;
  }
  amax *= GeV_to_erg;
  return amax;
}

/** shape of the gamma spectrum per pion decay from Kafexhiu 2014 Eq. 11
 */
double Radiation::F(double Tp, double Eg) {
  double alpha, beta, gamma, lambda;
  GetABGParams(Tp, alpha, beta, gamma, lambda);
  double Egmax = GetMaximumGammaEnergy(Tp);
  double Yg = Eg + m_pi * m_pi / (4. * Eg);
  double Ygmax = Egmax + m_pi * m_pi / (4. * Egmax);
  double Xg = (Yg - m_pi) / (Ygmax - m_pi);
  double C = lambda * m_pi / Ygmax;
  double f = pow(1. - pow(Xg, alpha), beta) / pow(1. + Xg / C, gamma);
  return f;
}

/** maximum gamma energy from decay of pion with kinetic energy Tp.
 *  from Kafexhiu 2014, Eq. 10
 */
double Radiation::GetMaximumGammaEnergy(double Tp) {
  double gammapilab = Epilabmax(Tp) / m_pi;
  double betapilab = sqrt(1. - 1. / (gammapilab * gammapilab));
  return 0.5 * m_pi * gammapilab * (1. + betapilab);
}
/** minimum gamma energy from decay of pion with kinetic energy Tp.
 *  from Kafexhiu 2014, Eq. 10
 */
double Radiation::GetMinimumGammaEnergy(double Tp) {
  double gammapilab = Epilabmax(Tp) / m_pi;
  double betapilab = sqrt(1. - 1. / (gammapilab * gammapilab));
  return (m_pi / 2.) * gammapilab * (1. - betapilab);
}

/** maximum pion energy in the lab frame (Kafexhiu 2014, Eq. 10)
 */
double Radiation::Epilabmax(double Tp) {
  double s = 2. * m_p * (Tp + 2. * m_p);
  double EpiCM = (s - 4. * m_p * m_p + m_pi * m_pi) / (2. * sqrt(s));
  double gammaCM = (Tp + 2. * m_p) / sqrt(s);
  double PpiCM = sqrt(EpiCM * EpiCM - m_pi * m_pi);
  double betaCM = sqrt(1. - 1. / (gammaCM * gammaCM));
  double Epilabmax = gammaCM * (EpiCM + PpiCM * betaCM);
  return Epilabmax;
}

/** correction factor that accounts for the abundance of heavier nuclei in
 *  the medium.(Eq. 24)
 */
double Radiation::NuclearEnhancementFactor(double Tp) {
  double Tp0 = 1.e3 * GeV_to_erg;
  double G;
  double sigmaTp = InelasticPPXSectionKaf(Tp);
  double sigmaTp0 = InelasticPPXSectionKaf(Tp0);
  /* Eq. 19 */
  (sigmaTp / sigmaTp0 > 1.) ? (G = 1. + log(sigmaTp / sigmaTp0)) : G = 1.;
  double sigmaRpp = 31.4e-27;
  /* these values are derived from local galactic ISM values, see Kafexhiu paper
   * p.13*/
  double epsilonc = 1.37;
  double epsilon1 = 0.29;
  double epsilon2 = 0.1;
  double epsilon = epsilonc + (epsilon1 + epsilon2) * (sigmaRpp * G) / sigmaTp;
  return epsilon;
}

/** get the alpha beta and gamma parameters used in Eq. 11, Kafexhiu 2014
 *  this is basically an implementation of Table V in the paper.
 *  is lambda right for Tp<1GeV? in the paper, a dash is listed...
 */
void Radiation::GetABGParams(double Tp, double &alpha, double &beta,
                             double &gamma, double &lambda) {
  Tp /= GeV_to_erg;
  double q = (Tp - 1.) / (m_p / GeV_to_erg);
  double mu = 1.25 * pow(q, 1.25) * exp(-1.25 * q);
  double kappa = 3.29 - 0.2 * pow(Tp / (m_p / GeV_to_erg), -1.5);
  if (Tp < 1.)
    lambda = 1., alpha = 1., beta = kappa, gamma = 0.;
  else if (Tp < 4.)
    lambda = 3., alpha = 1., beta = mu + 2.45, gamma = mu + 1.45;
  else if (Tp < 20.)
    lambda = 3., alpha = 1., beta = 1.5 * mu + 4.95, gamma = mu + 1.5;
  else if (Tp < 100.)
    lambda = 3., alpha = 0.5, beta = 4.2, gamma = 1.;
  else if (PiModel == 0)
    lambda = 3., alpha = 0.5, beta = 4.9, gamma = 1.;
  else if (PiModel == 2)
    lambda = 3.55, alpha = 0.5, beta = 3.6, gamma = 1.;
  else if (PiModel == 3)
    lambda = 3.55, alpha = 0.5, beta = 4.5, gamma = 1.;
  if (Tp > 50. && PiModel == 1)
    lambda = 3.5, alpha = 0.5, beta = 4., gamma = 1.;
  return;
}

/** get the a parameters used in Eq. 7, Kafexhiu 2014
 *  this is basically an implementation of Table IV in the paper.
 */
void Radiation::GetAParams(double Tp, double &a1, double &a2, double &a3,
                           double &a4, double &a5) {
  Tp /= GeV_to_erg;
  if (Tp < 100.)
    a1 = 0.728, a2 = 0.596, a3 = 0.491, a4 = 0.2503, a5 = 0.117;
  else if (PiModel == 0)
    a1 = 0.728, a2 = 0.596, a3 = 0.491, a4 = 0.2503, a5 = 0.117;
  else if (PiModel == 2)
    a1 = 5.436, a2 = 0.254, a3 = 0.072, a4 = 0.075, a5 = 0.166;
  else if (PiModel == 3)
    a1 = 0.908, a2 = 9.e-4, a3 = 6.089, a4 = 0.176, a5 = 0.448;
  if (Tp > 50. && PiModel == 1)
    a1 = 0.652, a2 = 1.6e-3, a3 = 0.488, a4 = 0.1928, a5 = 0.483;
  return;
}

/** get the b parameters used in Eq. 12, Kafexhiu 2014
 *  this is basically an implementation of Table VII in the paper.
 */
void Radiation::GetBParams(double Tp, double &b1, double &b2, double &b3) {
  Tp /= GeV_to_erg;
  if (Tp < 5.)
    b1 = 9.53, b2 = 0.52, b3 = 0.054;
  else if (Tp < 100.)
    b1 = 9.13, b2 = 0.35, b3 = 9.7e-3;
  else if (PiModel == 0)
    b1 = 9.13, b2 = 0.35, b3 = 9.7e-3;
  else if (PiModel == 2)
    b1 = 10.77, b2 = 0.412, b3 = 1.264e-2;
  else if (PiModel == 3)
    b1 = 13.16, b2 = 0.4419, b3 = 1.439e-2;
  if (Tp > 50. && PiModel == 1) b1 = 9.06, b2 = 0.3795, b3 = 1.105e-2;
  return;
}

/* ----         END OF RADIATION MODELS        ----  */

/**
 * Set a gsl interpolation object for fast reading of the proton
 * spectrum. x = energy, y = differential number. x has to be strictly
 * ordered ascending in energy!
 * Units: [x]=erg, [y]=1/erg
 */
void Radiation::SetParticles(vector<vector<double> > PARTICLES, int type) {
  if (type && type != 1) {
    cout << "Radiation::SetParticles: particle type unknown! Either "
            "electrons(type=0) or protons(type=1). Exiting!" << endl;
  }
  if (!PARTICLES.size()) {
    if (!type)
      cout << "Radiation::SetParticles: electron vector empty. Exiting."
           << endl;
    else
      cout << "Radiation::SetParticles: proton vector empty. Exiting." << endl;
    return;
  }
  if(!type && ElectronLookup && !QUIETMODE) {
      cout << "Radiation::SetParticles: Overwriting existing electron lookup."
           << endl;
  }
  if(type && ProtonLookup && !QUIETMODE) {
      cout << "Radiation::SetParticles: Overwriting existing proton lookup."
           << endl;
  }
  int size = (int)PARTICLES.size();
  double x[PARTICLES.size()];
  double y[PARTICLES.size()];
  for (unsigned int i = 0; i < PARTICLES.size(); i++) {
    x[i] = PARTICLES[i][0] > 0. ? log10(PARTICLES[i][0]) : -100.;
    y[i] = PARTICLES[i][1] > 0. ? log10(PARTICLES[i][1]) : -100.;
  }
  if (!type) {
    ElectronLookup = gsl_spline_alloc(gsl_interp_linear, size);
    gsl_spline_init(ElectronLookup, x, y, size);
  } else {
    ProtonLookup = gsl_spline_alloc(gsl_interp_linear, size);
    gsl_spline_init(ProtonLookup, x, y, size);
  }
  return;
}

void Radiation::SetElectrons(vector<vector<double> > ELECTRONS) {
  vector<vector<double> > *eladr = &ElectronVector;
  *eladr = ELECTRONS;
  SetParticles(ElectronVector, 0);
  return;
}

void Radiation::SetProtons(vector<vector<double> > PROTONS) {
  vector<vector<double> > *pradr = &ProtonVector;
  *pradr = PROTONS;
  SetParticles(ProtonVector, 1);
  return;
}

/* Here comes code that defines the spectral distributions of target photons in
 * the IC process.
 * The general idea is that you can add different components to the
 * "TargetPhotonGraphs" vector
 * storing TGraphs, which are then in the end added up to
 * "TotalTargetPhotonGraph", which is the
 * final, total radiation field in the IC process.
 */


void Radiation::ClearTargetPhotonField(int i) {



    if(i==-2) {

        if (VERBOSEMODE == true) {
            cout << "Clearing sum of all photon target field spectra..."
              <<endl;}
        fUtils->Clear2DVector(TargetPhotonVectorSumAll);
        TargetPhotonLookupSumAll = NULL;
        TargetPhotonEdensSumAll = 0.;
        accall = NULL;

        fUtils->Clear2DVector(ICLossVectorSumAll);
        ICLossLookupSumAll = NULL;
        ICLossLookupAccAll = NULL;
    }
    else if(i==-1) {

        if (VERBOSEMODE == true) {
            cout << "Clearing sum of isotropic target photon field spectra..."
              <<endl;}
        fUtils->Clear2DVector(TargetPhotonVectorSumIso);
        TargetPhotonLookupSumIso = NULL;
        TargetPhotonEdensSumIso = 0.;
        acciso = NULL;

        fUtils->Clear2DVector(ICLossVectorSumIso);
        ICLossLookupSumIso = NULL;
        ICLossLookupAccIso = NULL;
    }
    else {
            
        if (VERBOSEMODE == true) {
            cout << "Clearing spectrum of target photon field "<<i<<" and any anisotropy for it..."
              <<endl;}
        i = (int)i;
        fUtils->Clear2DVector(TargetPhotonVectors[i]);
        TargetPhotonLookups[i] = NULL;
        TargetPhotonEdensities[i] = 0.;

        TargetPhotonAngularDistrs[i] = NULL;
        CosZetaLookups[i] = NULL;
        TargetPhotonAccs[i] = NULL;

        fUtils->Clear2DVector(ICLossVectors[i]);
        ICLossLookups[i] = NULL;
        ICLossLookupAccs[i] = NULL;

        phiaccescs[i] = NULL;
        thetaaccescs[i] = NULL;
        phiaccesc_zetas[i] = NULL;
        thetaaccesc_zetas[i] = NULL;
        TargetPhotonAngularBounds[i].clear();
    }
    return;

}


void Radiation::AddThermalTargetPhotons(double T, double edens, int steps) {
    SetThermalTargetPhotons(T,edens,steps,RADFIELD_COUNTER);
    RADFIELD_COUNTER++;
    return;
}

void Radiation::ResetWithThermalTargetPhotons(int i, double T, double edens, int steps) {
    if (i<0 || i>=(int)RADFIELDS_MAX) {
      cout<<"Radiation::ResetWithThermalTargetPhotons: Invalid index "<<i<<
            ". Exiting."<<endl;
      return;
    }
    if (!TargetPhotonVectors[i].size()) {
      cout<<"Radiation::ResetWithThermalTargetPhotons: Vector "<<i<<
            " not set before. Set it up first before resetting. Exiting."<<endl;
      return;
    }
    SetThermalTargetPhotons(T,edens,steps,i);
    return;
}

/** Add a greybody distribution of target photons to TotalTargetPhotonGraph,
 * which is used in the
 * IC emission process in this class, but which can also be used 'Particles'
 * class to calculate
 * IC cooling losses
 */
void Radiation::SetThermalTargetPhotons(double T, double edens, int steps, int i) {
  if (edens > 1.e-8)
    cout << "Radiation::AddThermalTargetPhotons: energy density of radiation "
            "field insane. Are you sure of this?" << endl;
  if(edens<=0.) {
    cout << "Radiation::AddThermalTargetPhotons: energy density of target "
            "radiation field negative or zero? Exiting." << endl;
    return;
  }

  double logemin, logemax, low_boundary, high_boundary;
  low_boundary = 1.e-12;
  high_boundary = 1.e6;
  logemin = log10(low_boundary * kb * T);
  logemax = log10(high_boundary * kb * T);

  double estep = (logemax - logemin) / steps;
  double ePhoton = 0.;
  double nPhoton = 0.;
  vector< vector<double> > vint;
  double loge;
  for (loge = logemin; loge < logemax; loge += estep) {
    ePhoton = pow(10., loge);
    nPhoton = GreyBody(ePhoton, T, edens);
    if(!nPhoton) continue;
    fUtils->TwoDVectorPushBack(loge,log10(nPhoton),vint);
  }
  SetTargetPhotonVectorLookup(vint,i);
  return;
}

/** Add an arbitray distribution of target photons to TotalTargetPhotonGraph,
 * which is used in the
 * IC emission process in this class, but which can also be used 'Particles'
 * class to calculate
 * IC cooling losses. This requires as input a 2D vector of format:
 *              ~~~    energy[erg] number_density   ~~~
 */
void Radiation::AddArbitraryTargetPhotons(vector<vector<double> > PhotonArray) {
    SetArbitraryTargetPhotons(PhotonArray,RADFIELD_COUNTER);
    RADFIELD_COUNTER++;
    return;
}

void Radiation::ResetWithArbitraryTargetPhotons(int i,vector<vector<double> > PhotonArray) {
    if (i<0 || i>=(int)RADFIELDS_MAX) {
      cout<<"Radiation::ResetWithThermalTargetPhotons: Invalid index "<<i<<
            ". Exiting."<<endl;
      return;
    }
    if (!TargetPhotonVectors[i].size()) {
      cout<<"Radiation::ResetWithThermalTargetPhotons: Vector "<<i<<
            " not set before. Set it up first before resetting. Exiting."<<endl;
      return;
    }
    SetArbitraryTargetPhotons(PhotonArray,i);
    return;
}

void Radiation::SetArbitraryTargetPhotons(vector<vector<double> > PhotonArray, int i) {
  vector< vector<double> > vint;
  for (unsigned int j = 1; j < PhotonArray.size() - 1; j++) {
    double E = PhotonArray[j][0];
    double N = PhotonArray[j][1];
    if(E <=0. || N <=0.) continue;
    fUtils->TwoDVectorPushBack(log10(E),log10(N),vint);
  }
  SetTargetPhotonVectorLookup(vint,i);
  return;
}

/** Import target photons from file. File has to be in ASCII format, namely:
 *              ~~~    energy[eV] number_density   ~~~
 * The photons will be added to TotalTargetPhotonGraph
 */
void Radiation::ImportTargetPhotonsFromFile(const char *phFile) {
    SetTargetPhotonsFromFile(phFile,RADFIELD_COUNTER);
    RADFIELD_COUNTER++;
    return;
}

void Radiation::ResetWithTargetPhotonsFromFile(int i,const char *phFile) {
    if (i<0 || i>=(int)RADFIELDS_MAX) {
      cout<<"Radiation::ResetWithTargetPhotonsFromFile: Invalid index "<<i<<
            ". Exiting."<<endl;
      return;
    }
    if (!TargetPhotonVectors[i].size()) {
      cout<<"Radiation::ResetWithTargetPhotonsFromFile: Vector "<<i<<
            " not set before. Set it up first before resetting. Exiting."<<endl;
      return;
    }
    SetTargetPhotonsFromFile(phFile,i);
    return;
}

void Radiation::SetTargetPhotonsFromFile(const char *phFile, int i) {
  ifstream PHfile(phFile);
  vector<vector<double> > v;
  while (1) {
    if (PHfile.eof()) break;
    double e = 0.;
    double n = 0.;
    PHfile >> e >> n;
    if (e <= 0. || n <= 0.) continue;
    fUtils->TwoDVectorPushBack(log10(TeV_to_erg * e * 1.e-12),
                               log10(1.e12 * n / TeV_to_erg),v);
  }
  vector< vector<double> > vint;
  for (unsigned int j = 0; j < v.size(); j++) {
    fUtils->TwoDVectorPushBack(v[j][0],v[j][1],vint);
  }
  PHfile.close();
  SetTargetPhotonVectorLookup(vint,i);
  return;
}

/** Add SSC target photons.
 * This function calls the Synchrotron code in this class.
 * The photons will be added to TotalTargetPhotonGraph
 * If 'UPDATE' is 'true'
 * then recalculate the synchroton target field and
 * replace the previous one by this updated field.
 * DANGER: for the 'UPDATE' option to work, the SSC field
 * must be the last entry in the 'TargetPhotonGraphs' vector!
 * It uses Atoyan&Aharonian1996: MNRAS, Volume 278, Issue 2, pp. 525-541
 */
void Radiation::AddSSCTargetPhotons(double R, int steps) {
    SetSSCTargetPhotons(R,steps,RADFIELD_COUNTER);
    RADFIELD_COUNTER++;
    return;
}

void Radiation::ResetWithSSCTargetPhotons(int i,double R, int steps) {
    if (i<0 || i>=(int)RADFIELDS_MAX) {
      cout<<"Radiation::ResetWithSSCTargetPhotons: Invalid index "<<i<<
            ". Exiting."<<endl;
      return;
    }
    if (!TargetPhotonVectors[i].size()) {
      cout<<"Radiation::ResetWithSSCTargetPhotons: Vector "<<i<<
            " not set before. Set it up first before resetting. Exiting."<<endl;
      return;
    }
    SetSSCTargetPhotons(R,steps,i);
    return;
}

void Radiation::SetSSCTargetPhotons(double R, int steps, int i) {
  if (R <= 0.) {
    cout
        << "Radiation::SetSSCTargetPhotons: Souce extension is <= 0... exiting!"
        << endl;
    return;
  }
  R *= pc_to_cm;
  ParticleVector = ElectronVector;
  if (!ParticleVector.size()) {
    cout << "Radiation::SetSSCTargetPhotons: No particles in spectrum... "
            "exiting!" << endl;
    return;
  }
  void *p = NULL;
  double logemin = log10(1.e-8 * eV_to_erg);
  double logemax = log10(1e-4*ParticleVector[ParticleVector.size()-1][0]);
//  if(logemax>-2.) logemax = -2.;
  double estep = (logemax - logemin) / steps;
  double E = 0.;
  double N = 0.;
  radiationMechanism = "Synchrotron";
  double U = 2.24;
  vector< vector<double> > vint;
  for (double loge = logemin; loge < logemax; loge += estep) {
    E = pow(10., loge);
    N = DifferentialEmissionComponent(E, p) * U / (4. * pi * R * R * c_speed);
    if(N <= 0.) continue;
    fUtils->TwoDVectorPushBack(loge,log10(N),vint);
  }
  SetTargetPhotonVectorLookup(vint,i);
  return;
}

vector< vector<double> > Radiation::SumTargetFields(int bins, bool ISO) {
    vector< vector<double> > vec;
    double eph_min = 100;
    double eph_max = -100;
    vector<double> minmax;
    for(unsigned int i=0;i<RADFIELDS_MAX;i++) {
        SetICLookups(i);
        if((*TargetPhotonVectorCurrent).size()) {
            if(ISO == true && ANISOTROPY_CURRENT == true) continue;
            minmax = fUtils->GetVectorMinMax(*TargetPhotonVectorCurrent,0);
            if(minmax[0] < eph_min) eph_min = minmax[0];
            if(minmax[1] > eph_max) eph_max = minmax[1];
            minmax.clear();
        }
    }
    double logestep = (eph_max - eph_min) / bins;

    for(double loge=eph_min;loge<eph_max;loge+=logestep){
        double sum = 0;
        for(unsigned int i=0;i<RADFIELDS_MAX;i++) {
            SetICLookups(i);
            if((*TargetPhotonVectorCurrent).size()) {
                if(ISO == true && ANISOTROPY_CURRENT == true) continue;
                minmax = fUtils->GetVectorMinMax(*TargetPhotonVectorCurrent,0);
                if(loge>=minmax[0] && loge<=minmax[1]) {
                    sum += pow(10.,fUtils->EvalSpline(loge,*TargetPhotonLookupCurrent,
                                                     *TargetAccCurrent,
                                                     __func__,__LINE__));
                minmax.clear();
                }
            }
        }
        fUtils->TwoDVectorPushBack(loge,log10(sum),vec);
    }
    return vec;
}



void Radiation::SumTargetFieldsAll(int bins) {
    vector< vector<double> > v = SumTargetFields(bins,false);
    if(v.size()) SetTargetPhotonVectorLookup(v,-2);
}

void Radiation::SumTargetFieldsIsotropic(int bins) {
    vector< vector<double> > v = SumTargetFields(bins,true);
    if(v.size()) SetTargetPhotonVectorLookup(v,-1);
}



/** Function that adds up all individual target photon contributions into
 * TargetPhotonVector,
 * which is what is then used by the code in the end.
 */
void Radiation::SetTargetPhotonVectorLookup(vector< vector<double> > v, int i) {


  if(!v.size()) {
    cout<< "Radiation::SetTargetPhotonVectorLookup: target photon vector empty. "
           "Returning." <<endl;
    return;
  }
  if(RADFIELD_COUNTER >= RADFIELDS_MAX) {
    cout << "Radiation::SetTargetPhotonVectorLookup: Maximum number of Radiation "
            "fields reached. If you want more, you can set the maximum number "
            "via the Radiation::SetRadfieldMaxNumber() function. Exiting." <<endl;
             exit(1);
  }

  ClearTargetPhotonField(-2);
  int size = (int)v.size();
  gsl_spline **spl;
  gsl_interp_accel **accspl;
  
  if(i>=0 && i<(int)RADFIELDS_MAX) {
    spl = &TargetPhotonLookups[i];
    accspl = &TargetPhotonAccs[i];
    if(TargetPhotonVectors[i].size()) ClearTargetPhotonField(i);
    TargetPhotonVectors[i] = v;
  }
  else if(i==-1) {
    spl = &TargetPhotonLookupSumIso;
    accspl = &acciso;
    if(TargetPhotonVectorSumIso.size()) ClearTargetPhotonField(i);
    TargetPhotonVectorSumIso = v;
  }
  else if(i==-2) {
    spl = &TargetPhotonLookupSumAll;
    accspl = &accall;
    if(TargetPhotonVectorSumAll.size()) ClearTargetPhotonField(i);
    TargetPhotonVectorSumAll = v;
  }
  else{
      cout << "Radiation::SetTargetPhotonVectorLookup: index " <<i<<
              " not valid. Exiting!" << endl;
      return;
  }

  double e[size];
  double n[size];
  double elin[size];
  double en[size];
  double logEOld = 0;
  for (unsigned int j = 0; j < (unsigned int)size; j++) {
    double logE = v[j][0];
    if (logE < logEOld && logEOld) {
      cout << "Radiation::SetTargetPhotonVectorLookup: Target field not "
              "ordered ascending in energy! Exiting!" << endl;
      return;
    }
    double logN = v[j][1];
    e[j] = logE;
    n[j] = logN;
    elin[j] = pow(10., e[j]);
    en[j] = pow(10., e[j]) * pow(10., n[j]);
    logE = logEOld;
  }
  double emin = pow(10., e[0]);
  double emax = pow(10., e[size - 1]);
  *spl = gsl_spline_alloc(gsl_interp_linear, size);
  gsl_spline_init(*spl, e, n, size);
  *accspl = gsl_interp_accel_alloc();

  gsl_spline *edens_int = gsl_spline_alloc(gsl_interp_linear, size);
  gsl_interp_accel *edens_int_acc = gsl_interp_accel_alloc();
  gsl_spline_init(edens_int, elin, en, size);

  double edens = 0.;
  if (gsl_spline_eval_integ_e(edens_int, emin,
                              emax, edens_int_acc, &edens))
    edens = 0.;

  
  if(i>=0 && i<(int)RADFIELDS_MAX) TargetPhotonEdensities[i] = edens; 
  else if(i==-1) TargetPhotonEdensSumIso = edens;
  else if(i==-2) TargetPhotonEdensSumAll = edens;
  gsl_spline_free(edens_int);
  gsl_interp_accel_free(edens_int_acc);
  return;
}


void Radiation::CheckSanityOfTargetPhotonLookup() {

//  gsl_interp_accel *acc = gsl_interp_accel_alloc();
//  for(unsigned int i=0;i<TargetPhotonVectorSum.size();i++) {
//    double e = TargetPhotonVectorSum[i][0];
//    double n = TargetPhotonVectorSum[i][1];
//    double nl = fUtils->EvalSpline(e,TargetPhotonLookupSum,acc,__func__,__LINE__);
//    cout << "rel. diff: " << nl/n - 1. << " " << nl << " " << n <<endl;
//  }
  return;
}

/** remove the latest component in
 * TotalTargetPhotonVector and recompute
 * the total target photon spectrum
 */
//void Radiation::RemoveLastICTargetPhotonComponent() {
//  TargetPhotonVector = TargetPhotonVectorOld;
//  SetTargetPhotonVectorLookup();
//  return;
//}


void Radiation::SetTargetPhotonAnisotropy(int i, vector<double> obs_angle, 
                                          vector<double> phi, vector<double> theta, 
                                          vector< vector<double> > mesh) {

    if(TargetPhotonAngularDistrs[i] != NULL || CosZetaLookups[i] != NULL) {
        TargetPhotonAngularDistrs[i] = NULL;
        CosZetaLookups[i] = NULL;
        TargetPhotonAngularBounds[i].clear();
        gsl_interp_accel_free(phiaccescs[i]);
        gsl_interp_accel_free(thetaaccescs[i]);
        fUtils->Clear2DVector(TargetPhotonAngularDistrsVectors[i]);
        TargetPhotonAngularPhiVectors[i].clear();
        TargetPhotonAngularPhiVectors[i].clear();
    }
    // first set the electron beam direction
    phi_e = obs_angle[0]; theta_e = obs_angle[1]; 
    sin_phi_e = sin(phi_e); 
    cos_phi_e = cos(phi_e);
    sin_theta_e = sin(theta_e); 
    cos_theta_e = cos(theta_e);
    d_phi = phi[1]-phi[0];
    d_theta = theta[1]-theta[0];

    
    // now set the target field anisotropy map
    vector< vector<double> > ani = fUtils->MeshgridToTwoDVector(phi,theta,mesh);

    vector<double> extrema = fUtils->GetVectorMinMax(ani,2);

    ani_minval = extrema[0];
    ani_maxval = extrema[1];
    TargetPhotonAngularDistrs[i] = 
      fUtils->TwoDsplineFromTwoDVector(ani,phi_min,phi_max,theta_min,theta_max);
    TargetPhotonAngularBounds[i].push_back(phi_min);
    TargetPhotonAngularBounds[i].push_back(phi_max);
    TargetPhotonAngularBounds[i].push_back(theta_min);
    TargetPhotonAngularBounds[i].push_back(theta_max);
    TargetPhotonAngularBounds[i].push_back(ani_minval);
    TargetPhotonAngularBounds[i].push_back(ani_maxval);

    TargetPhotonAngularDistrsVectors[i] = mesh;
    TargetPhotonAngularPhiVectors[i] = phi;
    TargetPhotonAngularThetaVectors[i] = theta;

    phiaccescs[i] = gsl_interp_accel_alloc();
    thetaaccescs[i] = gsl_interp_accel_alloc();
    phiaccesc_zetas[i] = gsl_interp_accel_alloc();
    thetaaccesc_zetas[i] = gsl_interp_accel_alloc();
    //FillCosZetaLookup(i);
    ANISOTROPY[i] = true;
    
    if ( distance == 0.0 ){
        SetDistance(1./pc_to_cm);
        if ( !QUIETMODE ) 
            cout << "\nRadiation::SetTargetPhotonAnisotropy: So far no distance specified. Since an anisotropic photon field was defined, not the total luminosity but the radiation in the observation direction is calculated.\n";
    }

    return;
}





void Radiation::SetTargetPhotonAnisotropy(int i, 
                                          vector<double> phi, vector<double> theta, 
                                          vector< vector<double> > mesh) {

    if(TargetPhotonAngularDistrs[i] != NULL || CosZetaLookups[i] != NULL) {
        TargetPhotonAngularDistrs[i] = NULL;
        CosZetaLookups[i] = NULL;
        TargetPhotonAngularBounds[i].clear();
        gsl_interp_accel_free(phiaccescs[i]);
        gsl_interp_accel_free(thetaaccescs[i]);
        fUtils->Clear2DVector(TargetPhotonAngularDistrsVectors[i]);
        TargetPhotonAngularPhiVectors[i].clear();
        TargetPhotonAngularPhiVectors[i].clear();
    }

    d_phi = phi[1]-phi[0];
    d_theta = theta[1]-theta[0];

    
    // now set the target field anisotropy map
    vector< vector<double> > ani = fUtils->MeshgridToTwoDVector(phi,theta,mesh);

    vector<double> extrema = fUtils->GetVectorMinMax(ani,2);

    ani_minval = extrema[0];
    ani_maxval = extrema[1];
    TargetPhotonAngularDistrs[i] = 
      fUtils->TwoDsplineFromTwoDVector(ani,phi_min,phi_max,theta_min,theta_max);
    TargetPhotonAngularBounds[i].push_back(phi_min);
    TargetPhotonAngularBounds[i].push_back(phi_max);
    TargetPhotonAngularBounds[i].push_back(theta_min);
    TargetPhotonAngularBounds[i].push_back(theta_max);
    TargetPhotonAngularBounds[i].push_back(ani_minval);
    TargetPhotonAngularBounds[i].push_back(ani_maxval);

    TargetPhotonAngularDistrsVectors[i] = mesh;
    TargetPhotonAngularPhiVectors[i] = phi;
    TargetPhotonAngularThetaVectors[i] = theta;

    phiaccescs[i] = gsl_interp_accel_alloc();
    thetaaccescs[i] = gsl_interp_accel_alloc();
    phiaccesc_zetas[i] = gsl_interp_accel_alloc();
    thetaaccesc_zetas[i] = gsl_interp_accel_alloc();

    ANISOTROPY[i] = true;
    ISOTROPIC_ELECTRONS = true;
    if ( distance == 0.0 ){
        SetDistance(1./pc_to_cm);
        if ( !QUIETMODE ) 
            cout << "\nRadiation::SetTargetPhotonAnisotropy: So far no distance specified. Since an anisotropic photon field was defined, not the total luminosity but the radiation in the observation direction is calculated.\n";
    }

    return;
}















/* Function to set an isotropic electron distribution for the case of   *
 * anisotropic inverse Compton scattering                               */
void Radiation::SetElectronsIsotropic(void){
    ISOTROPIC_ELECTRONS = true;
}


vector< vector<double> > Radiation::GetTargetPhotonAnisotropy(int i, 
                                     vector<double> phi, vector<double> theta) { 

    vector< vector< double > > mesh;
    if (TargetPhotonAngularDistrs[i] == NULL) {
        cout<<"Radiation::GetTargetPhotonAnisotropy: No anisotropy set for field "
            <<i<<". Returning empty vector." <<endl;
        return mesh;
    }
    phi_min = TargetPhotonAngularBounds[i][0];
    phi_max = TargetPhotonAngularBounds[i][1];
    theta_min = TargetPhotonAngularBounds[i][2];
    theta_max = TargetPhotonAngularBounds[i][3];

    double Q = 0.;
    for(unsigned int k=0;k<theta.size();k++) {
        mesh.push_back(vector<double>());
        for(unsigned int j=0;j<phi.size();j++) {
            if (phi[j]>=phi_min && phi[j] <= phi_max && 
                theta[k] >= theta_min && theta[k] <= theta_max) {
                Q = interp2d_spline_eval(TargetPhotonAngularDistrs[i],phi[j],
                                     theta[k],phiaccescs[i],thetaaccescs[i]);
            }
            else Q = 0.;
            mesh[mesh.size()-1].push_back(Q);
        }
    }
    return mesh;

}

/**
 * Wrapper around CalculateDifferentialPhotonSpectrum(vector<double> points) that
 * will create an log-evenly distributed set of energy points between emin and emax.
 */
void Radiation::CalculateDifferentialPhotonSpectrum(int steps, double emin,
                                                    double emax) {
  if (!emin || !emax) {
    cout << "Radiation::CalculateDifferentialPhotonSpectrum: Please specify non "
            "zero values for the boundaries of your spectrum. Exiting..." << endl;
    return;

  }
  if (emin > emax) {
    cout << "Radiation::CalculateDifferentialPhotonSpectrum: emin>emax! Check your "
            "boundaries. Exiting..." << endl;
    return;
  }
  if (!steps) {
    cout << "Radiation::CalculateDifferentialPhotonSpectrum: Requested 0 steps! "
            "Exiting..." << endl;
    return;
  }

  double estep = (log10(emax) - log10(emin)) / steps;
  vector<double> epoints;

  if (QUIETMODE == false)
    cout << "** Calculating differential gamma-ray emission:" << endl;
  double loge;
  int tt,jj;
  for (loge = log10(emin), jj = 1, tt = 1; loge < log10(emax);
       loge += estep, jj++) {
    if ((double)jj / steps > 0.01 * tt && QUIETMODE == false) {
      cout << "\r";
      cout << "    "
           << (int)(100. * (loge - log10(emin)) / (log10(emax) - log10(emin)))
           << "\% done" << std::flush;
      tt++;
    }
    epoints.push_back(pow(10., loge));
  }
  
  CalculateDifferentialPhotonSpectrum(epoints);
  return;

}

/** Calculate differential photon spectra for the different radiation processes.
 *  Spectral points will be calculated for energy points given in 'points'.
 *  These points have to be in units of 'erg'!
 *  They are stored in the 2D 'diffspec' vector and can be accessed via the
 *  Radiation::ReturnDifferentialSpectrum() member function.
 */
void Radiation::CalculateDifferentialPhotonSpectrum(vector<double> points) {
  if (!points.size()) {
    cout << "Radiation::ReturnDifferentialSpectrum: Emin>Emax! Check your "
            "boundaries. Exiting..." << endl;
    return;
  }
  if (diffSpec.size()) fUtils->Clear2DVector(diffSpec);
  if (diffSpecICComponents.size()) fUtils->Clear2DVector(diffSpecICComponents);
  if (!ElectronVector.size() && !ProtonVector.size()) {
    cout << "Radiation::ReturnDifferentialSpectrum: No particle spectra filled "
            "-> No gamma spectra to calculate. Exiting..." << endl;
    return;
  }
  if (!QUIETMODE) {
    cout << "_________________________________________" << endl;
    cout << ">> CALCULATING SED FROM PARENT PARTICLES " << endl;
  }

  double ICVal, SynchVal, BremsVal, ppVal, E, Emin, Emax;
  if (!ProtonVector.size()) {
    Emax = ElectronVector[ElectronVector.size() - 1][0];
    Emin = ElectronVector[0][0] * 1.e-6;
    if (BField)
      Emin *= 1.e-10;  // because in this case we have to go to radio energies
                       // (synchrotron)
  } else if (!ElectronVector.size()) {
    Emax = ProtonVector[ProtonVector.size() - 1][0];
    Emin = ProtonVector[0][0] * 1.e-6;
  } else {
    (ElectronVector[ElectronVector.size() - 1][0] >
     ProtonVector[ProtonVector.size() - 1][0])
        ? (Emax = ElectronVector[ElectronVector.size() - 1][0])
        : (Emax = ProtonVector[ProtonVector.size() - 1][0]);
    Emin = ElectronVector[0][0] * 1.e-6;
    if (BField)
      Emin *= 1.e-10;  // because in this case we have to go to radio energies
                       // (synchrotron)
  }
  if (QUIETMODE == false)
    cout << "** Calculating differential gamma-ray emission:" << endl;
  int size = (int)points.size();
  epoints_temp = points;
  if (ElectronVector.size() && FASTMODE_IC == true && RADFIELD_COUNTER) {
    SumTargetFieldsIsotropic();
    SumTargetFieldsAll();
    }
  for (int i = 0; i < size; i++) {
    if (QUIETMODE == false) {
      cout << "\r";
      cout << "    "
           << i+1 << " / " << size
           << " points calculated" << std::flush;
    }
    E = points[i];

    ICVal = SynchVal = BremsVal = ppVal = 0.;

    if (ElectronVector.size()) {
      CalculateDifferentialGammaEmission(E, 0);
      ICVal = GetDifferentialICFlux();
      SynchVal = GetDifferentialSynchFlux();
      BremsVal = GetDifferentialBremsFlux();
    } else
      ICVal = SynchVal = BremsVal = 0.;
    if (ProtonVector.size()) {
      CalculateDifferentialGammaEmission(E, 1);
      ppVal = GetDifferentialPPFlux();
    } else
      ppVal = 0.;

    diffSpec.push_back(vector<double>());
    diffSpec[diffSpec.size() - 1].push_back(E);
    diffSpec[diffSpec.size() - 1]
        .push_back(ppVal + ICVal + BremsVal + SynchVal);
    diffSpec[diffSpec.size() - 1].push_back(ppVal);
    diffSpec[diffSpec.size() - 1].push_back(ICVal);
    diffSpec[diffSpec.size() - 1].push_back(BremsVal);
    diffSpec[diffSpec.size() - 1].push_back(SynchVal);

    
    diffSpecICComponents.push_back(vector<double>());
    diffSpecICComponents[diffSpecICComponents.size() - 1].push_back(E);
    for(unsigned int i=0; i<RADFIELDS_MAX; i++) {
        diffSpecICComponents[diffSpecICComponents.size() - 1].push_back(fdiffics[i]);
    }
  }
  if (QUIETMODE == false) {
    cout << endl;
    cout << "    -> DONE!   " << endl;
    cout << endl;
    cout << ">> SED CALCULATION DONE. EXITING." << endl;
    cout << endl;
  }
  return;
}
/**
 * Return the differential photon spectra dN/dE [number of photons/(erg*s*cm^2)]
 * vs E [erg] in a 2D-Vector
 * between the energies #emin and #emax. The argument #i determines the spectral
 * component to return:
 * - i = 1 : total spectrum
 * - i = 2 : pi0 decay
 * - i = 3 : IC scattering
 * - i = 4 : Bremsstrahlung
 * - i = 5 : Synchrotron radiation
 */
vector<vector<double> > Radiation::ReturnDifferentialPhotonSpectrum(
    int i, double emin, double emax ,vector< vector<double> > vec) {
  vector<vector<double> > tempVec;
  if (!vec.size()) {
    cout << "Radiation::ReturnDifferentialSpectrum: Differential spectrum "
            "vector empty. Fill it via "
            "Radiation::CalculateDifferentialPhotonSpectrum() first! Returning empty "
            "vector." << endl;
    return tempVec;
  }
  double e, dNdE;
  for (unsigned int j = 0; j < vec.size(); j++) {
    dNdE = vec[j][i];
    if (dNdE <= 0. || std::isnan(dNdE)) continue;
    e = vec[j][0];
    if (e < emin && emin) continue;
    if (e > emax && emax) continue;
    fUtils->TwoDVectorPushBack(e,dNdE,tempVec);
  }
  return tempVec;
}

/** Return the photon SED (EdN/dE (erg/s/cm^2) vs E (TeV) in a 2D-Vector of
 *  i = 1 : total spectrum
 *  i = 2 : pi0 decay
 *  i = 3 : IC scattering
 *  i = 4 : Bremsstrahlung
 *  i = 5 : Synchrotron radiation
 */
vector<vector<double> > Radiation::ReturnSED(int i, double emin, double emax, 
                                                vector< vector<double> > vec) {
  vector<vector<double> > tempVec;
  if (!vec.size()) {
    cout << "Radiation::ReturnSED: Differential spectrum vector empty. Fill it "
            "via Radiation::CalculateDifferentialPhotonSpectrum() first! Returning "
            "empty vector." << endl;
    return tempVec;
  }
  double e, eTeV, dNdE;
  for (unsigned int j = 0; j < vec.size(); j++) {
    dNdE = vec[j][i];
    if (dNdE <= 0. || std::isnan(dNdE)) continue;
    e = vec[j][0];
    eTeV = e / TeV_to_erg;
    if (eTeV < emin && emin) continue;
    if (eTeV > emax && emax) continue;
    fUtils->TwoDVectorPushBack(eTeV,e * e * dNdE,tempVec);
  }
  return tempVec;
}

vector<vector<double> > Radiation::GetICSpectrum(unsigned int i, double emin, double emax) {

    if(IC_CALCULATED && FASTMODE_IC==true && epoints_temp.size()) {
        FASTMODE_IC=false;
        CalculateDifferentialPhotonSpectrum(epoints_temp);
        FASTMODE_IC=true;
    }
    return ReturnDifferentialPhotonSpectrum(i+1, emin, emax, diffSpecICComponents);
}  ///< return pi0 decay spectrum


vector<vector<double> > Radiation::GetICSED(unsigned int i, double emin, double emax) {
    if(IC_CALCULATED && FASTMODE_IC==true && epoints_temp.size()) {
        FASTMODE_IC=false;
        CalculateDifferentialPhotonSpectrum(epoints_temp);
        FASTMODE_IC=true;
    }
    return ReturnSED(i+1, emin, emax, diffSpecICComponents);
  }  ///< return pi0 decay spectrum

/**
 * Return a particle SED dN/dE vs E (erg vs TeV)
 */
vector<vector<double> > Radiation::GetParticleSED(string type) {
  vector<vector<double> > v;
  if(!type.compare("electrons")) {
    if(!ElectronVector.size()) {
      cout<<"Radiation::GetParticleSED: electron vector empty! Returning "
                  "empty vector!" << endl;
      return v;
    }
    else v = ElectronVector;
  }
  else if(!type.compare("protons")) {
    if(!ProtonVector.size()) {
      cout<<"Radiation::GetParticleSED: proton vector empty! Returning "
                  "empty vector!" << endl;
      return v;
    }
    else v = ProtonVector;
  }
  else {
    cout << "Radiation::GetParticleSED: unknown particle type >"<< type <<"<. "
            "Supported are 'electrons' and 'protons'. Returning empty vector"
         << endl;
    return v;
  }

  for(unsigned int i=0;i<v.size();i++) {
    double E = v[i][0];
    double E_in_TeV = E/TeV_to_erg;
    double N = v[i][1];

    v[i][0] = E_in_TeV;
    v[i][1] = E*E*N;
  }

  return v;
}

double Radiation::GetIntegratedFlux(int i, double emin, double emax, bool ENERGYFLUX) {

  if (!diffSpec.size()) {
    cout << "Radiation::GetIntegratedFlux: Differential spectrum "
            "vector empty. Fill it via "
            "Radiation::CalculateDifferentialSpectrum() first! Returning zero."
            << endl;
    return 0.;
  }
  if (!emin) emin = diffSpec[0][0];
  if (!emax) emax = diffSpec[diffSpec.size()-1][0];
  vector <vector <double> > tempVec;
  double e, dNdE, val;
  for (unsigned int j = 0; j < diffSpec.size(); j++) {
    dNdE = diffSpec[j][i];
    if (dNdE <= 0. ) continue;
    e = diffSpec[j][0];
    if (ENERGYFLUX == false) fUtils->TwoDVectorPushBack(e,dNdE,tempVec);
    else fUtils->TwoDVectorPushBack(e,e*dNdE,tempVec);
  }
  if (tempVec.size() < 3) return 0.;
  fUtils->ToggleQuietMode();
  val = fUtils->Integrate(tempVec,emin,emax);
  fUtils->ToggleQuietMode();
  return val;
}


/**
 * Integration function using the GSL QAG functionality
 *
 */
Radiation::fPointer Radiation::_funcPtr;
Radiation *Radiation::_radPtr;

double Radiation::evaluate(double x, void* params) {
  return (_radPtr->*_funcPtr)(x, params);
}


double Radiation::Integrate(fPointer f, double *x, double emin, double emax,
                            double tolerance, int kronrodrule) {
  double integral, error;

  fPointer ftemp = _funcPtr;
  gsl_function F;
  _funcPtr = f;
  _radPtr = this;
  F.function = &Radiation::evaluate;
  F.params = x;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);
  int val = gsl_integration_qag(&F, emin, emax, 0, tolerance, 10000, kronrodrule, w,
                          &integral, &error);
  gsl_integration_workspace_free(w);
  _funcPtr = ftemp;
  if (val)  return 0.;
  else return integral;
}

vector<vector<double> > Radiation::GetTargetPhotons(int i) {

    if(i<-1 || i > (int)RADFIELDS_MAX) {
        cout<<"Radiation::GetTargetPhotons: Index "
            << i << " not valid. Exiting."<<endl;
        exit(1);}
    if(i==-1) {
        SumTargetFieldsAll();
        TargetPhotonVectorCurrent = &TargetPhotonVectorSumAll;
        TargetPhotonLookupCurrent = &TargetPhotonLookupSumAll;
        TargetAccCurrent = &accall;
    }
    else SetICLookups(i);

    vector< vector<double> >  vint;

    if (!RADFIELD_COUNTER) {
        cout<<"Radiation::GetTargetPhotons: No target photons set. "
            <<"Returning empty vector."<<endl;
        return vint;
    }

    double logemin = fUtils->GetVectorMinMax(*TargetPhotonVectorCurrent,0)[0];
    double logemax = fUtils->GetVectorMinMax(*TargetPhotonVectorCurrent,0)[1];
    int bins = (int)(*TargetPhotonVectorCurrent).size();
    
    double logestep = (logemax - logemin) / bins;

    for (double loge=logemin;loge<logemax;loge+=logestep) {
        
      double ph = fUtils->EvalSpline(loge,*TargetPhotonLookupCurrent,
                             *TargetAccCurrent,__func__,__LINE__);
      fUtils->TwoDVectorPushBack(loge,ph,vint);
    }
    return fUtils->VectorAxisPow10(vint,-1);
}
