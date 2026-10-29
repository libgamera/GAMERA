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
  ISOTROPIC_ELECTRONS = false;  //If true, calculate anisotropic IC scattering with isotropic electrons
  SPATIALDEP_CURRENT = false;
  DEFAULT_HADRON_COMPOSITION = false; //TODO: Check if it is still needed, otherwise delete it
  USE_KELNER = false;
  USE_KAMAE = false;
  lumtoflux = 0.;
  ldiffbrems = fdiffbrems = ldiffsynch = fdiffsynch = 0.;
  ldiffic = fdiffic = ldiffpp = fdiffpp = 0.;
  fdiffhadr.clear();
  distance = BField = 0.;
  phi_min = phi_max = theta_min = theta_max = phi_e = theta_e = 0.;
  sin_phi_e = cos_phi_e = sin_theta_e = cos_theta_e = 0.;
  ani_minval = ani_maxval = 0.;
  ElectronLookup = NULL;
  ProtonLookup = NULL;
  TargetPhotonLookupSumIso = NULL;
  TargetPhotonLookupSumAll = NULL;
  ICLossLookupSumIso = NULL;
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
  SpatialDep.resize(RADFIELDS_MAX);
  SPATIALDEP.resize(RADFIELDS_MAX);
  sizephfield.resize(RADFIELDS_MAX);

  // This is used for calculation of hadronic interactions of nuclei
  HadronMasses.clear();
  HadronSpectra.clear();
  HadronSpectraLookups.clear();
  AmbientMediumComposition.clear();
  diffSpecHadronComponents.clear();
  
  // This are the vectors used for the neutrino calculations
  HadronicMuonNeutrinoVectors.clear();
  HadronicElectronNeutrinoVectors.clear();
  HadronicTotalNeutrinoVectors.clear();
  
  ProtonMuonNeutrinoVector.clear();
  ProtonElectronNeutrinoVector.clear();
  ProtonTotalNeutrinoVector.clear();
  
  MuonNeutrinoVector.clear();
  ElectronNeutrinoVector.clear();
  TotalNeutrinoVector.clear();
  
  
  for(unsigned int i=0;i<RADFIELDS_MAX;i++) {
    fdiffics[i] = NAN;
    TargetPhotonEdensities[i] = 0.;
    ANISOTROPY[i] = false;
    SPATIALDEP[i] = false;
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
    sizephfield[i]=0.;
  }
  TargetPhotonEdensSumIso = 0.;
  fUtils->Clear2DVector(TargetPhotonVectorSumAll);
  fUtils->Clear2DVector(TargetPhotonVectorSumIso);
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
  /// for spatial dependency of the gamma-gamma absorption target field
  // accsp = gsl_interp_accel_alloc();
}

/**
 * Destructor. All allocated space has to be freed.
 */
Radiation::~Radiation() {
    delete fUtils;
    gsl_interp_accel_free(acc);
    gsl_interp_accel_free(acciso);
    //if (accall != NULL){ gsl_interp_accel_free(accall);}
    gsl_interp_accel_free(ICLossLookupAccIso);
    gsl_interp_accel_free(ICLossLookupAccAll);
    gsl_interp_accel_free(loraccesc);
    gsl_interp_accel_free(edaccesc);
    
    //gsl_spline_free(ICLossLookupSumAll);
    //gsl_spline_free(ICLossLookupSumIso);
    //gsl_spline_free(ElectronLookup);
    //gsl_spline_free(ProtonLookup);
    
    //gsl_spline_free(TargetPhotonLookupSumAll);
    //gsl_spline_free(TargetPhotonLookupSumIso);
    
    for (int i = 0; i < (int)TargetPhotonLookups.size(); i++){
        gsl_spline_free(TargetPhotonLookups[i]);
    }
    for (int i = 0; i < (int)TargetPhotonAccs.size(); i++){
        gsl_interp_accel_free(TargetPhotonAccs[i]);
    }
    
    for (int i = 0; i < (int)ICLossLookups.size(); i++){
        gsl_spline_free(ICLossLookups[i]);
    }
    for (int i = 0; i < (int)ICLossLookupAccs.size(); i++){
        gsl_interp_accel_free(ICLossLookupAccs[i]);
    }
    
    for (int i = 0; i < (int)phiaccescs.size();i++){
        gsl_interp_accel_free(phiaccescs[i]);
    }
    for (int i = 0; i < (int)thetaaccescs.size();i++){
        gsl_interp_accel_free(thetaaccescs[i]);
    }
    for (int i = 0; i < (int)phiaccesc_zetas.size();i++){
        gsl_interp_accel_free(phiaccesc_zetas[i]);
    }
    for (int i = 0; i < (int)thetaaccesc_zetas.size();i++){
        gsl_interp_accel_free(thetaaccesc_zetas[i]);
    }
    for (int i=0; i< (int)CosZetaLookups.size(); i++){
        interp2d_spline_free(CosZetaLookups[i]);
    }
    
    for (int i = 0; i< (int)TargetPhotonAngularDistrs.size(); i++){
        interp2d_spline_free(TargetPhotonAngularDistrs[i]);
    }
    
    for (int i = 0; i< (int)HadronSpectraLookups.size(); i++){
        gsl_spline_free(HadronSpectraLookups[i]);
    }
    
}

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

/************************************************
 * Remove all previously defined hadrons
************************************************/
void Radiation::ClearHadrons() {
    HadronMasses.clear();
    HadronSpectra.clear();
    for (int i = 0; i< (int)HadronSpectraLookups.size(); i++){
        gsl_spline_free(HadronSpectraLookups[i]);
    }
    HadronSpectraLookups.clear();
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
 * - hadrons ('particletype'=2)
 *   + inelastic scattering from different nuclear species
 *     that can be intriduced by the user
 *
 * This function calculates the appropriate (depending on the particle species)
 * radiation mechanism and their gamma-ray flux by calling
 * the method 'DifferentialEmissionComponent'.
 */
void Radiation::CalculateDifferentialGammaEmission(double e, int particletype) {
  ldiffbrems = fdiffbrems = ldiffsynch = fdiffsynch = 0.;
  ldiffic = fdiffic = ldiffpp = fdiffpp = 0.;
  fdiffhadr.clear();

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
  if (particletype != 0. && particletype != 1 && particletype != 2) {
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
            // std::cout<<i<<" "<<TargetPhotonAngularDistrs[i]<<" "<<FASTMODE_IC<<std::endl;
            if(TargetPhotonLookups[i]!=NULL) {
                if(FASTMODE_IC == true && TargetPhotonAngularDistrs[i] == NULL) 
                    continue;
                // if(FASTMODE_IC == false && TargetPhotonAngularDistrs[i] != NULL)
                //     continue;
                SetICLookups(i);
                // cout << "Filling field number: "<<i<<endl;
                ldiffic = DifferentialEmissionComponent(e, p);
                fdiffics[i] = lumtoflux * ldiffic;
                // for (int k=0;k<3;k++){
                //     cout<<"fdiffics["<<k<<"]: "<<fdiffics[k]<<endl;
                // 	}
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
    current_mass_number = 1.0;
    radiationMechanism = "ppEmission";
    current_Hadron_lookup = ProtonLookup;
    ldiffpp = DifferentialEmissionComponent(e, p);
    fdiffpp = lumtoflux * ldiffpp;
  } else if (particletype == 2) {
    radiationMechanism = "hadronicEmission";
    double ldiffhadr = 0.;
    for(int i = 0; i < (int)HadronSpectra.size(); i ++){
        ParticleVector = HadronSpectra[i];
        current_mass_number = HadronMasses[i];
        current_Hadron_lookup = HadronSpectraLookups[i];
        ldiffhadr = DifferentialEmissionComponent(e, p);
        fdiffhadr.push_back(lumtoflux * ldiffhadr);
    }
  } else {
    cout << "### Radiation::CalculateDifferentialGammaEmission: "
    		"particletype value not recognized (must be either "
    		"0, 1 or 2)" << endl;
  }
  return;
}


/**
 * The function selects the field \a i as the 'current' field
 * updating all the needed quantities and lookups.
 */
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
 */  // TODO: the input parameter void *par became useless? if it is not used -> delete it to avoid confusion!!
double Radiation::DifferentialEmissionComponent(double e, void *par) {
  if (radiationMechanism.compare("Synchrotron") &&
      radiationMechanism.compare("Bremsstrahlung") &&
      radiationMechanism.compare("InverseCompton") &&
      radiationMechanism.compare("ppEmission") &&
      radiationMechanism.compare("hadronicEmission")) {
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
    if (!n && !AmbientMediumComposition.size()) {
      if(!QUIETMODE) cout << "Radiation::DifferentialEmissionComponent:"
                             "No ambient density value set for "
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
    if(USE_KELNER) {IntFunc = &Radiation::PPEmissivityKelner;}
    else if (USE_KAMAE) {IntFunc = &Radiation::PPEmissivityKamae;}
    else {IntFunc = &Radiation::PPEmissivity;}
  } else if (!radiationMechanism.compare("hadronicEmission")) {
    if(USE_KELNER) { IntFunc = &Radiation::PPEmissivityKelner;}
    else if (USE_KAMAE) {IntFunc = &Radiation::PPEmissivityKamae;}
    else {IntFunc = &Radiation::PPEmissivity;}
	  if (!n && !AmbientMediumComposition.size()) {
	    if(!QUIETMODE) cout << "Radiation::DifferentialEmissionComponent:"
	                           "No ambient density value set for "
	                           "p-p scattering. Returning zero value." << endl;
	    return 0.;
	  }
      
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

  /* determine integration boundaries for the target photon energy from Eq. 2.50
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
 * Emissivity from IC scattering.
 * Taken from Blumenthal and Gould 1970, Eq(2.48).
 * 
 * Arguments:
 * \arg double x = log10 of target photon energy
 * \arg void *par = parameter vector:
 *    - electron energy
 *    - resulting gamma-ray photon
 * 
 * \return Emissivity: dN(electronE,targetE_ph)/(dt/dscatteredE_ph)
 * with scatteredE_ph in units of the initial electron energy
 */
double Radiation::ICEmissivity(double x, void *par) {
  double ephoton = pow(10.,x);  // energy of the target photon
  double *p = (double *)par;
  double lorentz = p[0] / m_e;
  double egamma = p[1];  // energy of the resulting gamma photon
  double e1 = egamma / (lorentz * m_e);  // gamma-ray energy in units of the
                                         // electron energy
  double gamma = 4. * ephoton * lorentz / m_e;  // parameter that describes
                                                // the regime of the scattering
                                                // process. Small value: Thomson
                                                // regime, large: KN-regime
  double q = e1 / (gamma * (1. - e1));  // yet another parameter telling us
                                        // the scattering domain
  double f = 1./(4.*lorentz*lorentz);
  if(f > 0.1) return 0.;
  if(q>1. ||  q<f) return 0.;
  // Eq(2.48):
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
 * IC emissivity in anisotropic radiation field.
 * 
 * From Moskalenko & Strong,
 * Astrophys.J. 528 (2000) 357-367.
 * 
 * Arguments:
 * \arg double x = log10 of target photon energy
 * \arg void *par = parameter vector:
 *    - electron energy
 *    - resulting gamma-ray photon
 * 
 * \return Emissivity: dN(electronE,targetE_ph)/(dt/dscatteredE_ph)
 * with scatteredE_ph in units of the initial electron energy
 */
double Radiation::ICEmissivityAnisotropic(double x, void *par) {
    double ephoton = pow(10.,x);  // energy of the target photon
    double *p = (double *)par;
    double lorentz = p[0] / m_e;
    double egamma = p[1];  // energy of the resulting gamma photon

    double beta = sqrt(1. - 1. / (lorentz*lorentz));
    double Q = 0.; double F = 0.; double cos_zeta = 0.;
    double cos_zeta_min = egamma/(2.*ephoton*lorentz*(lorentz-egamma/c_speed))-1.; 

    // double cos_zeta_min2 = (egamma/m_e)/(2.*ephoton/m_e*lorentz*(lorentz-egamma/m_e))-1.;
    double zeta_min = acos(cos_zeta_min);
    double integral = 0.;
    double cos_kappa = cos(theta_e); double sin_kappa = sin(theta_e);   
    double pi_min_ka = pi - theta_e; 

    double phi_min = (*TargetPhotonAngularBoundsCurrent)[0];
    double phi_max = (*TargetPhotonAngularBoundsCurrent)[1];
    double theta_min = (*TargetPhotonAngularBoundsCurrent)[2];
    double theta_max = (*TargetPhotonAngularBoundsCurrent)[3];
    for (double phi = phi_min; phi <= phi_max; phi += d_phi) {
        // std::cout<<phi<<"("<<phi_min<<" - "<<phi_max<<")"<<std::endl;
        for (double theta = theta_min; theta <= theta_max; theta += d_theta) {
            // cos_zeta = interp2d_spline_eval(*CosZetaLookupCurrent, phi, theta,
            //                   *phiaccesc_zetaCurrent,*thetaaccesc_zetaCurrent);
            cos_zeta = -cos_kappa *cos(theta) + sin_kappa * sin(theta) * cos(phi - phi_e);
            // std::cout<<"cosz, cosz_min = "<<cos_zeta<<" "<<cos_zeta_min<<" "<<cos_zeta_min2<<std::endl;
            // cos_zeta  = cos(phi) * cos(phi_e) * sin(theta) * sin(theta_e);
            // cos_zeta += sin(phi) * sin(phi_e) * sin(theta) * sin(theta_e);
            // cos_zeta += cos(theta) * cos(theta_e);
            // std::cout<<"cos_zeta1,cos_zeta "<<cos_zeta1<<" "<<cos_zeta<<std::endl;
            // cos_zeta = cos(phi) * sin(theta);
            // std::cout<<phi<<" "<<theta<<" "<<Q<<" "<<F<<std::endl;
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
            // std::cout<<"Q,F,int,cz,dph,dth: "<<Q<<","<<F<<","<<integral<<","<<cos_zeta<<","<<d_phi<<","<<d_theta<<std::endl;
        }
    }

    double targetphotons = fUtils->EvalSpline(x,
                                              *TargetPhotonLookupCurrent,
                                              *TargetAccCurrent,__func__,__LINE__);
    double integrand = pi * e_radius * e_radius * c_speed;
    integrand /= ephoton * (lorentz-egamma/m_e) * (lorentz-egamma/m_e);
    integrand *= integral * pow(10.,targetphotons);
    integrand *= ln10 * ephoton;
    // std::cout<<"-------------------------------------"<<std::endl;
    return integrand;
}




/*----------------------------------------------------------------------------------------------------
 *      This Block calculates the Anisotropic IC emission for isotropic electrons
 *      with equation 20 from Aharonian & Atoyan 1981
 *--------------------------------------------------------------------------------------------------*/

/**
 * Anisotropic IC emissivity with isotropic electrons.
 * 
 * From Aharonian&Atoyan, Ap&SS 1981
 */
double Radiation::ICEmissivityAnisotropicIsotropicElectrons(double x, void *par){
    double ephoton = pow(10.,x);  // energy of the target photon
    double *p = (double *)par;
    double eelectron = p[0];
    double egamma = p[1];  // energy of the resulting gamma photon
    
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

//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------


/**
 * DEPRECATED: Filling lookup table for the cosine of the scattering angle between
 * photon and electron.
 * This function is not used at the moment!
 */
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
    for (unsigned int k=0;k<v.size();k++){for (unsigned int l=0;l<v[k].size();l++){ cout<<v[k][l]<<endl;}}
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


/** 
 * Useful to
 * apply in spectral iterations in
 * conjunction in the "Particles" class.
 * 
 * \param bins : number of bins in logarithmic energy
 * 
 * Produces a lookup table with format E[erg] ; -1.*LossrateIC [erg/s]
 * 
 * It fills the Radiation::ICLossLookupSumIso lookup for all the isotropic fields
 * and and single lookups for each of the anisotropic ones with Radiation::ICLossLookups
 */
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

/** 
 * Fills a lookup table holding the differential electron energy loss rate due
 * to inverse-Compton scattering.
 * The format of the lookup is: { Energy(erg) - Energy Loss Rate
 * from IC scattering(erg/s) }
 * 
 * The values are returned as the log10.
 *
 * The function integrates Radiation::ICEmissivityRadFieldIntegrated for
 * a given energy of the electron and stores the results in Radiation::ICLossVectorCurrent
 * 
 * \param i = index of the photon field
 * \param bins = number of bins for the electron energy range
 * 
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

/**
 * Returns the Inverse Compton lookup for losses for the field \a i
 *
 * @param i : index of the
 * @return vector of 2Dtuple with the energy loss for each energy
 *         log10(E/[erg]) and log10(-dE/dt / [erg/s])
 */
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
            // cout<< "Radiation::GetICLossLookup: Lookup not calculated yet. "
            //        "Please run Radiation::CreateICLossLookup first. Returning "
            //        "empty vector."<<endl;
            CreateICLossLookup();
            // return v;
        }
        for(unsigned int j=0;j<ICLossVectorSumAll.size();j++) {
            double E = ICLossVectorSumAll[j][0];
            double L = fUtils->EvalSpline(E, ICLossLookupSumAll,ICLossLookupAccAll,
                                          __func__,__LINE__);
                
            fUtils->TwoDVectorPushBack(pow(10.,E),pow(10.,L),v);
        }
        return v;
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
 * modified Bessel functions
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
 * Emissivity of synchrotron radiation.
 * Emission at frequency nu from a number of electrons with energy eElectron 
 * adapted from galprop!
 * Ghisellini&Svensson, 1988: 'the synchrotron boiler'
 */
double Radiation::SynchEmissivity(double x, void *par) {
  /* frequency of emitted synchr. radiation */
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
 * Emissivity of synchrotron radiation with explicit
 * pitch angle. Default value is 90 degrees.
 * 
 * Following Blumenthal&Gould, 1970. Eqs 4.44 and 4.48
 */
// !!!CHECK THAT!!!
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
/** Emissivity of Bremsstrahlung,
 * proton-electron as well as electron-electron
 * From Baring 1999, ApJ, 513, 311-338  - Eq. 27
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
  // double S = n * 1.4; // OLD STUFF
  double S = 0;
  double Se = 0;
  for (unsigned int i=0;i<AmbientMediumComposition.size();i++){
     S += AmbientMediumComposition[i][1] * AtomicNumber(AmbientMediumComposition[i][0]) * AtomicNumber(AmbientMediumComposition[i][0]);
     Se += AmbientMediumComposition[i][1] * AtomicNumber(AmbientMediumComposition[i][0]);
  }
  /* emissivity assuming a fully ionised plasma (n_e = S), Eq. (27) */
  double N;
  if (EI < 2.e-3 * GeV_to_erg)
    //N = c_speed * b * S * (sigma1(g, e) + sigmaNR(g, e));
	N = c_speed * b * (S * sigma1(g, e) + Se * sigmaNR(g, e));
  else
    //N = c_speed * b * S * (sigma1(g, e) + sigma_e);
	N = c_speed * b * (S * sigma1(g, e) + Se * sigma_e);
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

  double logprotons = fUtils->EvalSpline(log10(EP),current_Hadron_lookup,
                                      acc,__func__,__LINE__);
  return c_speed * N * pow(10., logprotons) * ln10 * EP;
}



/**
 * differential cross section following Kafexhiu 2014
 * (Eq. 8)
 */
double Radiation::DiffPPXSection(double Tp, double Eg) {
  double NuclearEnhancement = CalculateEpsilon(Tp, current_mass_number);
  double dsigmadE =
      NuclearEnhancement *Amax(Tp) * F(Tp, Eg) / GeV_to_erg;
  
  return dsigmadE;
}

/** inclusive pi0 production cross section over the full energy range.
 *  following Kafexhiu 2014, Paragragh II B 4
 */
double Radiation::InclusivePPXSection(double Tp) {
  double InclusiveXSection;
  if (Tp < 2. * GeV_to_erg){
    InclusiveXSection = SigmaOnePi(Tp) + SigmaTwoPi(Tp);}
  else{
    InclusiveXSection = InelasticPPXSectionKaf(Tp) * MeanMultiplicity(Tp);}
  return InclusiveXSection;
}

/** inelastic proton-proton cross-section following
 *  Kafexhiu 2014 (Eq. 1)
 */
double Radiation::InelasticPPXSectionKaf(double Tp) {
  double Tpth = 0.2797 * GeV_to_erg;
  double L = Tp / Tpth;
  if (L < 1.0) return 0.0;   // Return 0 if we are below the threshold
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
    cout << "Radiation::MeanMultiplicity: please, specify supported "
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
  if ( Tp/GeV_to_erg < 0.2797) return 0.0;
  if (Tp/GeV_to_erg > 2.0) return 0.0;
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

/** Eq. 5 Kafexhiu 2014, needed for lowest energies
 */
double Radiation::SigmaTwoPi(double Tp) {
  if ((Tp / GeV_to_erg) < 0.56) { return 0.0;}
  else {
  double SigmaTwoPi = 5.7 / (1. + exp(-9.3 * (Tp / GeV_to_erg - 1.4)));
  return SigmaTwoPi * 1.e-27;
  }
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
  if (Xg > 1.0) return 0.0;
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





/******************************************************************
 * With this function one can set the parametrization of
 * pp and hadronic collisions to the one by Kamae et al. 2006
 * \param kamae = Boolean, if true Kamae 2006 par. will be used
 ******************************************************************/
void Radiation::UseKamae(bool kamae){
    if (kamae) {
     USE_KELNER = false;
     USE_KAMAE = kamae;
    }
    else {
     USE_KAMAE = kamae;   
    }
    return;
}




double Radiation::PPEmissivityKamae(double x, void *par){
  /* proton energy */
  double EP = pow(10.,x);
  /* pi0 decay photon energy */
  double Eg = pow(10.,*(double *)par);
  if (EP <= m_p) return 0.;
  double Tp = sqrt(EP * EP - m_p * m_p); 
  double Tpth = 0.2797 * GeV_to_erg; 
  if (Tp <= Tpth) return 0.;
  
  
  // Calculation of the parameters in table 3 in Kamae et al. 2006
  double y = log(Tp * erg_to_TeV);
  
  double a0 = -0.51187*(y + 3.3) + 7.6179*(y + 3.3)*(y + 3.3) - 2.1332*pow(y + 3.3, 3) + 0.22184*pow(y + 3.3, 4);
  double a1 = -1.2592e-5 + 1.4439e-5 * exp(-0.29360 * ( y + 3.4)) + 5.9363e-5/( y + 4.1485) + 2.2640e-6 * y - 3.3723e-7 * y*y;
  double a2 = -1.0*(174.83 + 152.78) * log(1.5682 * ( y + 3.4)) -  808.74/(y + 4.6157);
  double a3 = 0.81177 + 0.56385*y + 0.0040031*y*y - 0.0057658*y*y*y + 0.00012057*pow(y,4);
  double a4 = 0.68631*( y + 3.32) + 10.145* pow( y + 3.32 , 2) - 4.6176*pow( y + 3.32 , 3) + 0.86824*pow( y + 3.32 , 4) - 0.053741*pow( y + 3.32 , 5);
  double a5 = 9.0466e-7 + 1.4539e-6 * log(0.015204*( y + 3.4)) + 1.3253e-4/((y + 4.7171)*(y + 4.7171)) - 4.1228e-7 * y + 2.2036e-7 *y*y;
  double a6 = -1.0*(339.45 + 618.73) *  log(0.31595*( y + 3.9)) + 250.20/(( y + 4.4395)*( y + 4.4395));
  double a7 = -35.105 + 36.167*y - 9.3575*y*y + 0.33717*y*y*y;
  double a8 = 0.17554 + 0.37300*y - 0.014938*y*y + 0.0032314*y*y*y + 0.0025579*pow(y, 4);
  double r = 1.01;
  if (Tp/GeV_to_erg <= 1.95){
   r = 3.05*exp( -107.0*pow((y + 3.25)/(1.0 + 8.08*(y + 3.25)), 2));   
  }
    
  double b0 = 60.142 * tanh(-0.37555*( y + 2.2)) - 5.9564*(y + 0.59913)*(y + 0.59913) + 6.0162e-3 *pow(y + 9.4773,4);
  double b1 = 35.322 + 3.8026 * tanh(-2.5979*(y + 1.9)) - 2.1870e-4 * (y + 369.13)*(y + 369.13);
  double b2 = -15.732 - 0.082064 * tanh(-1.9621*(y + 2.1)) + 2.3355e-4 * (y + 252.43) *(y + 252.43);
  double b3 = -0.086827 + 0.37646 * exp(-0.53053*pow(((y + 1.0444)/(1.0 + 0.27437*(y + 1.0444))),2));
  double b4 = 2.5982 + 0.39131*(y + 2.95)*(y + 2.95) - 0.0049693*pow( y + 2.95, 4) + 0.94131 * exp(-24.347*pow(y + 2.45 - 0.19717*(y + 2.45)*(y + 2.45),2));
  double b5 = 0.11198 - 0.64582*y + 0.16114*y*y + 2.2853 * exp(-0.0032432*pow(( y - 0.83562)/(1.0 + 0.33933*( y - 0.83562)),2) );
  double b6 = 1.7843 + 0.91914*y + 0.050118*y*y + 0.038096*y*y*y - 0.027334*pow(y, 4) - 0.0035556*pow(y, 5) + 0.0025742*pow(y, 6);
  double b7 = -0.19870 - 0.071003*y + 0.019328*y*y - 0.28321 * exp(-6.0516*(y + 1.8441)*(y + 1.8441));
  
  double c0 = 2.4316 * exp(-69.484*pow((y + 3.1301)/(1.0 + 1.24921*(y + 3.1301)),2)) - 6.3003 - 9.5349/y + 0.38121*y*y;
  double c1 = 56.872 + 40.627*y + 7.7528*y*y;
  double c2 = -5.4918 - 6.7872 * tanh(4.7128*(y + 2.1)) + 0.68048*y;
  double c3 = -0.36414 + 0.039777*y;
  double c4 = -0.72807 - 0.48828*y - 0.092876*y*y;
  
  double d0 = 3.2433 * exp(-57.133*pow((y + 2.9507)/(1.0 + 1.2912*(y + 2.9507)),2)) - 1.0640 - 0.43925*y;
  double d1 = 16.901 + 5.9539*y - 2.1257*y*y - 0.92057*y*y*y;
  double d2 = -6.6638 - 7.5010 * tanh(30.322*(y + 2.1)) + 0.54662*y;
  double d3 = -1.50648 - 0.87211*y - 0.17097*y*y;
  double d4 = 0.42795 + 0.55136*y + 0.20707*y*y + 0.027552*y*y*y;


  double Wndl = 15.0;
  double Wndh = 44.0;
  double Lmax = 0.96*log(Tp/GeV_to_erg);
  
  double X = log(Eg/GeV_to_erg);
  
  
  double Fnd = FndKamae(X, a0, a1, a2, a3, a4, a5, a6, a7, a8) * FndklKamae(X,Wndl, Wndh, Lmax)*r;
  
  double Fdiff = FdiffKamae(X, b0,b1,b2,b3,b4,b5,b6,b7) * FklKamae(X, Tp/GeV_to_erg);
  
  double FresDelta = FresKamae(X,c0,c1,c2,c3,c4) * FklKamae(X, Tp/GeV_to_erg);
  double Fres1600 = FresKamae(X,d0,d1,d2,d3,d4) * FklKamae(X, Tp/GeV_to_erg);
  
  double N = Fnd * SigmaNdKamae(Tp) + Fdiff * SigmaDiffKamae(Tp) + FresDelta * SigmaDeltaKamae(Tp) + Fres1600 * SigmaResKamae(Tp);
  
  double NuclearEnhancement = CalculateEpsilon(Tp, current_mass_number);
  
  double logprotons = fUtils->EvalSpline(log10(EP),current_Hadron_lookup,
                                      acc,__func__,__LINE__);
  return c_speed * NuclearEnhancement * N * pow(10., logprotons) * ln10 * EP * Eg;

    
}




/**************************************************************
 * Equation 6 from Kamae et al. 2006 for the non-diffractive
 * channel
 * ***********************************************************/
double Radiation::FndKamae(double x, double a0, double a1, double a2, double a3, double a4, double a5, double a6, double a7, double a8){
    //return 0.0;
    //cout << " a0: " << a0 << " a1: " << a1 << " a5: " << a5 << " arg_exp1: " << -a1* pow((x-a3 + a2*(x-a3)*(x-a3)),2.0) << " arg_exp2: " << -a5* pow((x-a8 + a6*(x-a8)*(x-a8)+ a7*pow((x-a8),3.0)) , 2.0);
    double result = a0 * exp(-a1* pow((x-a3 + a2*(x-a3)*(x-a3)),2.0))
                    + a4 * exp(-a5* pow((x-a8 + a6*(x-a8)*(x-a8)+ a7*pow((x-a8),3.0)) , 2.0));
    //cout << result << "\n";
    if (not isinf(result)) return result;
    else return 0.0;
}

/**************************************************************
 * Equation 7 from Kamae et al. 2006 for the non-diffractive
 * channel
 * ***********************************************************/
double Radiation::FndklKamae(double x, double Wndl, double Wndh, double Lmax){
    //cout << "In FndklKamae\n";
    //return 0.0;
    double Lmin = -2.6;
    double result = 1.0/(exp(Wndl*(Lmin - x)) + 1.0) * 1.0/(exp(Wndh*(x - Lmax)) + 1.0);
    
    return result;
}


/**************************************************************
 * Equation 9 from Kamae et al. 2006 for the non-diffractive
 * channel
 * ***********************************************************/
double Radiation::FdiffKamae(double x, double b0, double b1, double b2, double b3, double b4, double b5, double b6, double b7){
    double result = b0 * exp(-1.0*b1 * pow((x-b2)/(1.0 + b3*(x-b2)),2))
                    + b4*exp(-1.0*b5*pow((x-b6)/(1.0+b7*(x-b6)),2)); //DANGER: There is a typo in Kamae in equ. 9
    if (not isinf(result)) return result;
    else return 0.0;
}
// Equation 10 in Kamae
double Radiation::FklKamae(double x, double Tp){
    //cout << "In FklKamae\n";
    //return 0.0;
    double Lmax = log(Tp);
    double Wdiff = 75.0;
    double result = 1.0/(exp(Wdiff*(x-Lmax)) + 1.0);
    //if (not isinf(result)) return result;
    return result;
}

// Equation 12 in Kamae
double Radiation::FresKamae(double x, double c0, double c1, double c2, double c3, double c4){
    //cout << "In FresKamae\n";
    //return 0.0;
    double result = c0*exp(-1.0*c1* pow(((x-c2)/(1.0 + c3*(x-c2) + c4*(x-c2)*(x-c2))),2));
    return result;
}


/*Cross sections for Kamae. NOTE: There is an erratum to the original paper, the equations
 were taken from the erratum. */
// Equation 1
double Radiation::SigmaNdKamae(double Tp){
    //cout << "In SigmaNdKamae\n";
    //return 0.0;
    double x = log( Tp/GeV_to_erg/c_speed );
    double a0 = 0.1176, a1 = 0.3829, a2 = 23.10, a3 = 6.454, a4 = -5.764, a5 = -23.63, a6 = 94.75, a7 = 0.02667;
    double b0 = 11.34, b1 = 23.72;
    double c0 = 28.5, c1 = -6.133, c2 = 1.464;
    
    double sigma;
    if (x < 1.0) sigma = 0.0;
    else if (x < 1.3) sigma = 0.57*pow(x/a0,1.2)*(a2 + a3*x*x + a4*x*x*x + a5*exp(-a6*(x+a7)*(x+a7)));
    else if (x < 2.4) sigma = (b0*abs(a1-x) + b1*abs(a0-x) - x)/(a1-a0);
    else if (x < 10.0) sigma = a2 + a3*x*x + a4*x*x*x + a5*exp(-a6*(x+a7)*(x+a7));
    else sigma = c0 + c1*x + c2*x*x;
        
    return sigma*1.0e-27;
}

//Equation 2
double Radiation::SigmaDiffKamae(double Tp){
    //cout << "In SigmaDiffKamae\n";
    //return 0.0;
    double x = log( Tp/GeV_to_erg/c_speed );
    double d0 = 0.3522, d1 = 0.1530, d2 = 1.498, d3 = 2.0, d4 = 30.0, d5 = 3.155, d6 = 1.042;
    double e0 = 5.922, e1 = 1.632;
    
    double sigma;
    if (x<2.25) sigma = 0.0;
    else if (x< 3.2) sigma = sqrt((x-d0)/d1) * (d2 + d3*log(d4*(x-0.25)) + d5*x*x - d6*x*x*x);
    else if(x < 100.0) sigma = d2 + d3 * log(d4*(x-0.25)) + d5*x*x - d6*x*x*x;
    else sigma = e0 + e1*x;
    
    
    return sigma*1.0e-27;
}

// Equation 3
double Radiation::SigmaDeltaKamae(double Tp){
    //cout << "In SigmaDeltaKamae\n";
    //return 0.0;
    double Ep = sqrt(Tp*Tp + m_p*m_p)/GeV_to_erg;
    double f0 = 0.0834, f1 = 9.5, f2 = -5.5, f3 = 1.68, f4 = 3134.0;
    
    double sigma;
    if (Ep < 1.4) sigma = 0.0;
    else if (Ep < 1.6) sigma = f0 * pow(Ep,10);
    else if (Ep < 1.8) sigma = f1*exp(-f2*(Ep - f3)*(Ep - f3));
    else if (Ep < 10.0) sigma = f4* pow(Ep, -10);
    else sigma = 0.0;
    
    return sigma*1.0e-27;
}

// Equation 4
double Radiation::SigmaResKamae(double Tp){
    //cout << "In SigmaResKamae\n";
    //return 0.0;
    double Ep = sqrt(Tp*Tp + m_p*m_p)/GeV_to_erg;
    double g0 = 0.0004257, g1 = 4.5, g2 = -7.0, g3 = 2.1, g4 = 503.5;
    
    double sigma;
    if (Ep < 1.6) sigma = 0.0;
    else if (Ep < 1.9) sigma = g0* pow(Ep,14);
    else if (Ep < 2.3) sigma = g1 * exp(-g2*(Ep - g3)*(Ep - g3));
    else if (Ep < 20.0) sigma = g4 * pow(Ep, -6);
    else sigma = 0.0;
    
    return sigma*1.0e-27;
}










/******************************************************************
 * With this function one can set the parametrization of
 * pp and hadronic collisions to the one by Kelner et al. 2006
 * \par  - Boolean, if true Kelner 2006 par. will be used
 ******************************************************************/
void Radiation::UseKelner(bool kelner){
    if (kelner) {
     USE_KAMAE = false;
     USE_KELNER = kelner;
    }
    else {
     USE_KELNER = kelner;   
    }
    return;
}


/**************************************************************
 * PP emissivity according to Kelner et al. 2006
 * ***********************************************************/
double Radiation::PPEmissivityKelner(double x, void *par) {
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

  double N = Fgamma(Eg, EP); 

  double NuclearEnhancement = CalculateEpsilon(Tp, current_mass_number);
  
  double sigma = InelasticPPXSectionKaf(Tp);   // calculate the cross section
  //double sigma = InelasticPPXSectionKelner(EP);  // Here, the cross section from Kelner is used

  
  double logprotons = fUtils->EvalSpline(log10(EP),current_Hadron_lookup,
                                      acc,__func__,__LINE__);
  
  return (c_speed * NuclearEnhancement * sigma * pow(10., logprotons) * N * ln10);
}



/**
 * Inelastic PP interaction corss section according to Kelner et al. 2006.
 * \param EP = Total proton energy in [erg]
 * \return Cross section in [barn]
 * */
double Radiation::InelasticPPXSectionKelner(double EP){
  double L = log(EP*erg_to_TeV);
  double Eth = 1.22*GeV_to_erg;
  double sigma = (34.3 + 1.88*L + 0.25*L*L) * (1.0 - Eth*Eth*Eth*Eth/(EP*EP*EP*EP))*(1.0 - Eth*Eth*Eth*Eth/(EP*EP*EP*EP)) * 1.0e-27;
  
  return sigma;
}



/*******************************************************
 * Equations 58 to 61 in Kelner et al. 2006
 * Input:   - Energy of the gamma ray in [erg]
 *          - Energy of the proton in [erg]
 ******************************************************/
double Radiation::Fgamma(double Egamma, double Eproton){
    double x = Egamma/Eproton;
    double L = log(Eproton*erg_to_TeV);
    
    double Bgamma = 1.30 + 0.14*L + 0.011*L*L;
    double betagamma = 1.0/(1.79 + 0.11*L + 0.008*L*L);
    double kgamma = 1.0/(0.801 + 0.049*L + 0.014*L*L);
    
    double xbeta = pow(x, betagamma);
    double Fg = Bgamma * log(x)/x * pow(((1.0-xbeta)/(1.0 + kgamma*xbeta*(1.0-xbeta))),4.0)
                *( 1.0/log(x) - 4.0*betagamma*xbeta/(1.0-xbeta) - 
                4.0*kgamma*betagamma*xbeta*(1.0-2.0*xbeta) / (1.0 + kgamma*xbeta*(1.0-xbeta)));
    
    
    
    return Fg;
}



/*---- Calculation of neutrinos following Kelner et al 2006 --- */

/*************************************************************************
 * Calculation of the Neutrino spectrum. It fills the vectors
 * MuonNeutrinoVector, ElectronNeutrinoVector and TotalNeutrinoVector.
 * Input:   - Energy points where to calculate the spectrum in [erg]
 * **********************************************************************/
void Radiation::CalculateNeutrinoSpectrum(vector<double> points){
  if (MuonNeutrinoVector.size()) fUtils->Clear2DVector(MuonNeutrinoVector);
  if (ElectronNeutrinoVector.size()) fUtils->Clear2DVector(ElectronNeutrinoVector);
  if (TotalNeutrinoVector.size()) fUtils->Clear2DVector(TotalNeutrinoVector);
  
  if (ProtonMuonNeutrinoVector.size()) fUtils->Clear2DVector(ProtonMuonNeutrinoVector);
  if (ProtonElectronNeutrinoVector.size()) fUtils->Clear2DVector(ProtonElectronNeutrinoVector);
  if (ProtonTotalNeutrinoVector.size()) fUtils->Clear2DVector(ProtonTotalNeutrinoVector);
  
  if (HadronicMuonNeutrinoVectors.size()) fUtils->Clear3DVector(HadronicMuonNeutrinoVectors);
  if (HadronicElectronNeutrinoVectors.size()) fUtils->Clear3DVector(HadronicElectronNeutrinoVectors);
  if (HadronicTotalNeutrinoVectors.size()) fUtils->Clear3DVector(HadronicTotalNeutrinoVectors);
  
  if (!ProtonVector.size() && !HadronSpectra.size()) {
    cout << "Radiation::CalculateNeutrinoSpectrum: No proton or hadron spectra filled "
            "-> No neutrino spectra to calculate. Exiting..." << endl;
    return;
  }
  
  
  if( !AmbientMediumComposition.size()){
      vector<double> temp0; temp0.resize(2);
      vector<vector<double> > temp;
      if (!n) { temp0[0] = 1.0; temp0[1] = 0.0;
          temp.push_back(temp0);
          if(!QUIETMODE) cout << "Warning in Radiation::CalculateNeutrinoSpectrum: "
                             "No ambient medium defined "
                             "for p-p scattering. Values will be zero." << endl;
          
      }
      else { temp0[0] = 1.0; temp0[1] = n;
          temp.push_back(temp0);
      }
      SetAmbientMediumComposition(temp);
  }
  
  if (HadronSpectra.size()){
    HadronicMuonNeutrinoVectors.resize(HadronSpectra.size());
    HadronicElectronNeutrinoVectors.resize(HadronSpectra.size());
    HadronicTotalNeutrinoVectors.resize(HadronSpectra.size());
  }

  double distance_normalisation;
  if (!distance) {
        cout << "### Radiation::CalculateNeutrinoSpectrum: Distance to "
                "particles not specified -> Flux equals now the luminosity! ###"
             << endl;
      distance_normalisation = 1.;
    } else
      distance_normalisation = 1. / (4. * pi * distance * distance);
  
  
  
  if (!QUIETMODE) {
    cout << "___________________________________________________________" << endl;
    cout << ">> CALCULATING NEUTRINO SED FROM PARENT PROTONS AND HADRONS " << endl;
  }
  

  
  double E, nu1, nu2, total_flux;   // will contain the summed contribution from protons and all hadrons
  double nu1_protons = 0.0, nu2_protons = 0.0, total_flux_protons = 0.0;
  double nu1_hadrons = 0.0, nu2_hadrons = 0.0, total_flux_hadrons = 0.0;
  
  int size = (int)points.size();
 for(int i = 0; i < size; i++ ){
     if (QUIETMODE == false) {
      cout << "\r";
      cout << "    "
           << i+1 << " / " << size
           << " points calculated" << std::flush;
     }
     E = points[i];

     // If there are protons defined, calculate the neutrinos from protons
     if (ProtonVector.size()){
         ParticleVector = ProtonVector;
        current_Hadron_lookup = ProtonLookup;
        current_mass_number = 1.0;
        nu1_protons = CalculateNeutrinoFlux(E, 1)*distance_normalisation;
        nu2_protons = CalculateNeutrinoFlux(E, 0)*distance_normalisation;
        total_flux_protons = nu1_protons + 2.*nu2_protons;
        fUtils->TwoDVectorPushBack(E, nu1_protons, ProtonMuonNeutrinoVector);
        fUtils->TwoDVectorPushBack(E, nu2_protons, ProtonElectronNeutrinoVector);
        fUtils->TwoDVectorPushBack(E, total_flux_protons, ProtonTotalNeutrinoVector);
     }

     // If there are hadrons defined, calculate the neutrinos from all hadron species
     if (HadronSpectra.size()){
         nu1_hadrons = 0.0; nu2_hadrons = 0.0, total_flux_hadrons = 0.0;
         for(int j = 0; j < (int)HadronSpectra.size(); j++){
             ParticleVector = HadronSpectra[j];
             current_Hadron_lookup = HadronSpectraLookups[j];
             current_mass_number = HadronMasses[j];
             nu1 = CalculateNeutrinoFlux(E, 1)*distance_normalisation;
             nu2 = CalculateNeutrinoFlux(E, 0)*distance_normalisation;
             total_flux = nu1 + 2.*nu2;
             fUtils->TwoDVectorPushBack(E, nu1, HadronicMuonNeutrinoVectors[j]);
             fUtils->TwoDVectorPushBack(E, nu2, HadronicElectronNeutrinoVectors[j]);
             fUtils->TwoDVectorPushBack(E, total_flux, HadronicTotalNeutrinoVectors[j]);    
             
             nu1_hadrons += nu1;
             nu2_hadrons += nu2;
             total_flux_hadrons += total_flux;
         }
         
     }
     
     // Save the total flux from protons and all hadron species
     nu1 = nu1_protons + nu1_hadrons;
     nu2 = nu2_protons + nu2_hadrons;
     total_flux = total_flux_protons + total_flux_hadrons;
     fUtils->TwoDVectorPushBack(E, nu1, MuonNeutrinoVector);
     fUtils->TwoDVectorPushBack(E, nu2, ElectronNeutrinoVector);
     fUtils->TwoDVectorPushBack(E, total_flux, TotalNeutrinoVector);
 }
 if (QUIETMODE == false) {
    cout << endl;
    cout << "    -> DONE!   " << endl;
    cout << endl;
    cout << ">> NEUTRINO CALCULATION DONE. EXITING." << endl;
    cout << endl;
 }
 return;
}



/*****************************************************************************************
 *  Integrates and returns the result, used by the function CalculateNeutrinoSpectrum
 *   Input:  - Energy of the neutrinos (or electrons) in [erg]
 *           - leptontype: '0' for electron/muon neutrinos from decay of muons and '1' 
 *              for muons from direct charged pion decay
 *  Output: - Flux
 ****************************************************************************************/
double Radiation::CalculateNeutrinoFlux(double energy, int leptontype) {
    if (!n) {
        if ((AmbientMediumComposition.size() == 1) && (AmbientMediumComposition[0][1] == 0.0)) {
        return 0.;
        }
    }
    double integral = 0;
    const double ethresh = 1.22e-3 * TeV_to_erg;  // Threshold for pion production
    double emax = ParticleVector[ParticleVector.size() - 1][0];
    double emin = energy;
    if (emin < ParticleVector[0][0]) emin = ParticleVector[0][0];
    if (emin < ethresh) emin = ethresh; // Don't start integration below the threshold
    if (emin >= emax) return 0.; // We can not produce neutrinos more energetic than the protons
    

    // Contains the energy of the muon
    double integralinput[1];
    integralinput[0] = energy;    
    
    // Convert everything to logarithmic values
    double emin_log = log10(emin);
    double emax_log = log10(emax);
    
    // Set the function for the spectrum to electron neutrinos (or electrons) or muon neutrinos
    if (leptontype == 0) {
        fPointer Fnu = &Radiation::NeutrinoFlux1;
        integral = Integrate(Fnu, integralinput, emin_log, emax_log,integratorTolerance*5., integratorKronrodRule);
    }
    
    else if (leptontype == 1) {
        double emin_new = energy/0.427;
        if(emin_new > emin) emin_log = log10(emin_new);
        if(emax_log < emin_log) return 0.0;
        fPointer Fnu = &Radiation::NeutrinoFlux2;
        // integrate the function
        integral = Integrate(Fnu, integralinput, emin_log, emax_log,integratorTolerance*5., integratorKronrodRule);
    }
    else {
        cout << "Radiation::CalculateNeutrinoFlux: No valid lepton type specified.";
        return 0.;
    }

    double result = c_speed * integral * ln10;
    return result;
}



/**********************************************************************
 * Function to be integrated to get the neutrino emission for
 * electron- and muon-Neutrino produced in the muon decay,
 * as well as the electron spectrum.
 * Input:   - Energy of the proton (or hadron) in [erg]
 *          - energy of the muon/electron as a parameter in [erg]
 * Output:  - The integrand (see e.g. equ. 71 or 72 in Kelner et al.
 *              2006
 *********************************************************************/
double Radiation::NeutrinoFlux1(double energy_proton, void *par) {
    double *p = (double *)par;
    double energy_nu = p[0];
    
    double logprotons = fUtils->EvalSpline(energy_proton,current_Hadron_lookup,
                                      acc,__func__,__LINE__);
    double Jp = pow(10, logprotons);   // Proton flux, Jp in Kelner 2006
    
    // Convert the energy to the real value again
    energy_proton = pow(10, energy_proton);        
    double Tp = sqrt(energy_proton * energy_proton - m_p * m_p);
    double NuclearEnhancement = CalculateEpsilon(Tp, current_mass_number);
    
    double sigma = InelasticPPXSectionKaf(Tp);   // calculate the cross section
    //double sigma = InelasticPPXSectionKelner(energy_proton); // Cross section from Kelner
    double Fnu = Felectron(energy_nu, energy_proton);
    
    double result = sigma*Jp*Fnu*NuclearEnhancement;
    return result;
}



/**********************************************************************
 * Function to be integrated to get the muon-neutrino emission produced
 * in the pion decay.
 * Input:   - Energy of the proton in [erg]
 *          - energy of the muon as a parameter in [erg]
 * Output:  - The integrand (see e.g. equ. 71 or 72 in Kelner et al.
 *              2006
 *********************************************************************/
double Radiation::NeutrinoFlux2(double energy_proton, void *par) {
    double *p = (double *)par;
    double energy_nu = p[0];
    
    
    double logprotons = fUtils->EvalSpline(energy_proton,current_Hadron_lookup,
                                      acc,__func__,__LINE__);
    double Jp = pow(10, logprotons);   // Proton flux, Jp in Kelner 2006
    
    // Convert the energy to the real value again
    energy_proton = pow(10, energy_proton);        
    double Tp = sqrt(energy_proton * energy_proton - m_p * m_p);
    
    double sigma = InelasticPPXSectionKaf(Tp);   // calculate the cross section
    //double sigma = InelasticPPXSectionKelner(energy_proton);  // Cross section from Kelner
    double Fnu = Fnumu(energy_nu, energy_proton);
    
    double NuclearEnhancement = CalculateEpsilon(Tp, current_mass_number);
    
    double result = sigma*Jp*Fnu*NuclearEnhancement;
    return result;
}



/********************************************************************************* 
 * Function for electrons and muon neutrinos from pion and muon decay
 * Equations 62 to 65 from Kelner et al 2006
 * Input:
 *       - Energy of electrons or neutrinos, where the flux should
 *           be calculated in [erg]
 *       - Energy of protons in [erg]
 *   Output:
 *       - Flux
 ********************************************************************************/
double Radiation::Felectron(double Ee_erg, double Ep_erg) 
{

  double Ee = Ee_erg*erg_to_TeV;    // Electron/neutrino energy in TeV
  double Ep = Ep_erg*erg_to_TeV;


  double x = Ee/Ep; 
  //if (x<1.e-4) return 0.0; /// hack

  if (Ep < 0.04) Ep = 0.04; // HACK!!!!!!!!

  double L = log(Ep);
  double L2 = L*L; 
  double B = 1.0/(69.5 + 2.65*L + 0.3*L2);
  double beta = 1.0/pow(0.201 + 0.062*L+0.00042*L2,0.25); 
  double kappa = (0.279 + 0.141*L + 0.0172*L2)/(0.3 + pow(2.3 + L,2.0));
  
  double xb = pow(x,beta);

  double p1 = B;
  double p2 = pow(1+kappa*log(x)*log(x),3)/(x*(1 + 0.3/xb));
  double p3 = pow(-log(x),5);
  
  return p1*p2*p3;
}

/********************************************************************************* 
 * Function for muon neutrinos from pion decay
 * Equations 66 to 69 from Kelner et al 2006
 * Input:
 *       - Energy of neutrinos, where the flux should
 *           be calculated in [erg]
 *       - Energy of protons in [erg]
 *   Output:
 *       - Flux
 ********************************************************************************/
double Radiation::Fnumu(double Enu_erg, double Ep_erg)
{

  double Ep = Ep_erg*erg_to_TeV;
  double Enu = Enu_erg*erg_to_TeV;

  double x = Enu/Ep;
  double y = x/0.427;
  if (x > 0.427) return 0.0;

  double L = log(Ep);
  double L2 = L * L; 
  double B = 1.75 + 0.204 * L + 0.010 * L2;
  double beta = 1.0/(1.67 + 0.111 * L + 0.0038*L2); 
  double kappa = 1.07 - 0.086*L + 0.002*L2;

  double yb = pow(y,beta);

  double p1 = B*log(y)/y;
  double p2 = pow((1-yb)/(1+kappa*yb*(1-yb)),4);
  double p3 = (1/log(y)) 
    - (4*beta*yb/(1-yb)) 
    - 4*kappa*beta*yb*(1-2*yb)/(1+kappa*yb*(1-yb));
  
  return p1*p2*p3;
} 

/******** END of Neutrino calculations *******************************/






/* ----         END OF RADIATION MODELS        ----  */


/**
 * Set a gsl interpolation object for fast reading of the proton
 * spectrum. x = energy, y = differential number.
 * 
 * x has to be strictly ordered ascending in energy!
 * 
 * Units: [x]=erg, [y]=1/erg
 * 
 * \param PARTICLES = a vector of tuples (E; dN/dE)
 * \param type = 0 for electrons; 1 for protons
 *
 * NOTE: Does not work for the hadrons, use the function AddHadrons() instead
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


/** 
 * set the Electron spectrum (e.g. calculated in the "Particles" class,
 * but also arbitrary spectra).
 * 
 * \param ELECTRONS = vector of tuples (E[erg],N[erg^-1])
*/
void Radiation::SetElectrons(vector<vector<double> > ELECTRONS) {
  vector<vector<double> > *eladr = &ElectronVector;
  *eladr = ELECTRONS;
  SetParticles(ElectronVector, 0);
  return;
}

/** 
 * set the Proton spectrum (e.g. calculated in the "Particles" class,
 * but also arbitrary spectra).
 * 
 * \param PROTONS = vector of tuples (E[erg],N[erg^-1])
 *
 * TODO: This function is not really consistent with the AddHadrons
*/
void Radiation::SetProtons(vector<vector<double> > PROTONS) {
  vector<vector<double> > *pradr = &ProtonVector;
  *pradr = PROTONS;
  SetParticles(ProtonVector, 1);
  return;
}




/**
 * Function to add a spectrum of hadrons with a specific mass number. It also
 * calculates and saves a lookup for the spectrum.
 *
 * Inner functioning:
 *  - The energy used in the equations is the energy per nucleon.
 * Therefore in the function, the passed energy is divided by the
 * mass Mass_number before creating the lookup and saving the spectrum.
 * Furthermore, because we are dealing with a dN/dE quantity, we
 * multiply the dN/dE value that was given again for the Mass_number.
 *
 * @param Spectrum : spectrum of the hadrons (as vector of 2D-tuples (E,dN/dE) ([erg],[erg^-1]))
 * @param Mass_number : mass number of the hadrons (number of the nucleons)
 *
 */
void Radiation::AddHadrons(vector<vector<double> > Spectrum, double Mass_number){
  HadronMasses.push_back(Mass_number);
  // The energy used in the equations is the energy per nucleon. Therefore we divide with the Mass_number
  // before creating the lookup and saving the spectrum
  // Since we are dealing with dN/dE and not with dN, we have to multiply with the mass number
  for(unsigned int i = 0; i < Spectrum.size(); i++ ){
     Spectrum[i][0] = Spectrum[i][0]/Mass_number;
     Spectrum[i][1] = Spectrum[i][1]*Mass_number;
  }
  
  int size = (int)Spectrum.size();
  double x[Spectrum.size()];
  double y[Spectrum.size()];
  // Use the logarithmic values for the lookups:
  for (unsigned int i = 0; i < Spectrum.size(); i++) {
    x[i] = Spectrum[i][0] > 0. ? log10(Spectrum[i][0]) : -100.;
    y[i] = Spectrum[i][1] > 0. ? log10(Spectrum[i][1]) : -100.;
  }

  gsl_spline *HadronLookup;
  HadronLookup = gsl_spline_alloc(gsl_interp_linear, size);
  gsl_spline_init(HadronLookup, x, y, size);
  
  // If the lowest energy is lower than the rest mass energy, we have
  // to adapt the spectrum and the lookup
  if (Spectrum[0][0]*Spectrum[0][0] - m_p*m_p < 0.0){
    vector<vector<double> > tempVec;
    double value = fUtils->EvalSpline(m_p*1.001,HadronLookup,
                                      acc,__func__,__LINE__);
    
    for(unsigned int i = 0; i < Spectrum.size(); i++){
        if (Spectrum[i][0]*Spectrum[i][0] - m_p*m_p > 0.0) {
            if ((tempVec.size() == 0) && (Spectrum[i][0] > m_p*1.001)){
                tempVec.push_back({m_p*1.001, value});
            }
            tempVec.push_back({Spectrum[i][0], Spectrum[i][1]});
        }
    }
    double x2[tempVec.size()];
    double y2[tempVec.size()];
    for (unsigned int i = 0; i < tempVec.size(); i++) {
      x2[i] = tempVec[i][0] > 0. ? log10(tempVec[i][0]) : -100.;
      y2[i] = tempVec[i][1] > 0. ? log10(tempVec[i][1]) : -100.;
    }    
    gsl_spline *HadronLookup2;
    HadronLookup2 = gsl_spline_alloc(gsl_interp_linear, tempVec.size());
    gsl_spline_init(HadronLookup2, x2, y2, tempVec.size());
    
    HadronSpectra.push_back(tempVec);        // Save the spectrum
    HadronSpectraLookups.push_back(HadronLookup2);
    gsl_spline_free(HadronLookup);  // free the not needed allocated space of HadronLookup
  }
  else {
  HadronSpectra.push_back(Spectrum);        // Save the spectrum
  HadronSpectraLookups.push_back(HadronLookup);
  }
  
  return;
}


/**
 * Test the hadron lookup to see if the spectral returned is correct
 *
 * @param i : index of the hadron spectrum to consider
 * @param e : energy per nucleon at which test the lookup (erg)
 * @return : value of the spectrum as 1/(erg*cm^3)
 *
 * TODO: make private?
 */
double Radiation::TestHadronLookup(int i, double e){
 double value = fUtils->EvalSpline(e,HadronSpectraLookups[i],
                                      acc,__func__,__LINE__);
 return value;
}


/**
 * Set the default ambient composition with simply protons and 10% Helium.
 * This is for consistency with Bremsstrahlung and ionization processes
 * which use this default values for the composition of the ISM.
 *
 * @param N : proton number density in [cm^-3].
 */
void Radiation::SetAmbientDensity(double N){
	n = N; // number density for the protons
	vector<vector<double>> default_comps = {{1., n}, {4., 0.1*n}};
	AmbientMediumComposition = default_comps;
}
 
 
/**
 * Set the composition and the abundances of the ambient medium through a list of
 * mass numbers and number densities.
 *
 * @param composition : Vector of tuples containing the mass number (first entry) and
 *                      the density in [1/cm^3]
 */
void Radiation::SetAmbientMediumComposition(vector<vector< double > > composition){
    AmbientMediumComposition = composition;
    return;
}
 
 
/**
 * Function to return the Hadron spectrum number \a i.
 *
 * @param i : Number of the hadron component
 *
 * @return Spectrum, energy in [erg] and dN/dE in [1/(erg*cm^3)]
 */
vector<vector< double > > Radiation::GetHadrons(int i){ 
  vector<vector<double> > tempvec;
  if (i  >= (int)HadronSpectra.size()){
      cout << "In Radiation::GetHadrons: The hadron spectrum for i="
    		  << i << " does not exist, " << HadronSpectra.size()
			  << " different hadrons are defined so far. Returning empty vector. ";
      return tempvec;
  }
  tempvec = HadronSpectra[i];
  for(unsigned int j = 0; j < tempvec.size(); j++ ){
     tempvec[j][0] = tempvec[j][0]*HadronMasses[i];
     tempvec[j][1] = tempvec[j][1]/HadronMasses[i];
  }
  return tempvec;
}

/**
 * Function to return the mass numbers of the nuclei.
 *
 * @return Vector with the hadron masses in atomic units
 */
vector< double > Radiation::GetHadronMasses(void){
    return HadronMasses;
}

/**
 * Function to get back the composition and abundances of the ambient
 * medium.
 *
 * @return Vector of tuples containing the mass number (first entry)
 *         and the corresponding density in cm^-3
 */
vector<vector<double> > Radiation::GetAmbientMediumComposition(void){
    return AmbientMediumComposition;
}

/**
 * Implementation of equation 20 in Kafexhiu et al. 2014
 * for one projectile species with mass number 'Mass'.
 * @param Tp : Energy [erg]
 * @param Mass : mass number of the projectile
 *
 * @return Epsilon factor
 */
double Radiation::CalculateEpsilon(double Tp, double Mass){
    double epsilon = 0.;
    double density = 0.;
    double ambient_mass_number = 0.;
    double epsilon_component = 0.;
    double Tp0 = 1.e3 * GeV_to_erg;
    double G;
    double sigmaTp = InelasticPPXSectionKaf(Tp);
    double sigmaTp0 = InelasticPPXSectionKaf(Tp0);
    double sigma1 = 0.0;
    double sigma2 = 0.0;
    /* Eq. 19 */
    (sigmaTp / sigmaTp0 > 1.) ? (G = 1. + log(sigmaTp / sigmaTp0)) : G = 1.;
    //double sigmaRpp = 31.4e-27;
    for (int i = 0; i < (int)AmbientMediumComposition.size(); i++){
        ambient_mass_number = AmbientMediumComposition[i][0];
        density = AmbientMediumComposition[i][1];
        // If necessary use the pp inelastic cross section and not the nucleus nucleus one
        if(abs(ambient_mass_number - 1.0) <  0.5) sigma1 = sigmaTp;
        else sigma1 = NuNuXSection(1., ambient_mass_number)*G;
        if( abs(Mass - 1.0) < 0.5)  sigma2 = sigmaTp;
        else sigma2 = NuNuXSection(1.,Mass)*G;
        
        epsilon_component = density * (Mass*sigma1 + ambient_mass_number * sigma2);
        epsilon += epsilon_component;
    }
    epsilon = epsilon/(2.*sigmaTp);
    return epsilon;
}

/**
 * Calculates the nucleus-nucleus reaction cross section (equation 17 in Kafexhiu et al. 2014)
 * It is valid for projectile energies > 0.2 GeV for a proton projectile and > 0.1 GeV/nucleon for a nuclei
 * projectile different than a proton.
 *
 * @param ProjMass : Projectile mass number
 * @param TargetMass : Target mass number
 * 
 * @return Cross section in [mb]
 */
double Radiation::NuNuXSection(double ProjMass, double TargetMass){
  double sigma_r0 = 58.1*1.0e-27; //mb to cm
  double beta0 = 2.247 - 0.915*(1.0+pow(TargetMass,-1.0/3.)); //projectile is a proton
  //double beta0 = 1.581 - 0.876*(pow(ProjMass,-1/3.)+pow(TargetMass,-1/3.)); //projectile different than a proton. Not needed here
  double sigma_r = sigma_r0*pow((pow(ProjMass,1.0/3.)+pow(TargetMass,1.0/3.)-beta0*(pow(ProjMass,-1.0/3.)+pow(TargetMass,-1.0/3.))),2.0); // equation 17
  return sigma_r;
}

/**
 * Function to return the average atomic number corresponding to
 * a certain mass number.
 * This is the semi-empirical formula (least square coefficients)
 * The relation is an approximation!
 *
 * @param mass_number : mass number of the atomic species
 */
double Radiation::AtomicNumber(double mass_number){
	if (mass_number < 3) {return 1;}
	if (mass_number == 4) {return 2;}
	double aC = 0.714;
	double aA = 23.2;
	double Z = mass_number / (2 + (aC/(2*aA)) * pow(mass_number,2/3));
	return Z;
}

/****** End of Functions for hadronic interactions of nuclei *****************/


/* Here comes code that defines the spectral distributions of target photons in
 * the IC process.
 * The general idea is that you can add different components to the
 * "TargetPhotonGraphs" vector
 * storing TGraphs, which are then in the end added up to
 * "TotalTargetPhotonGraph", which is the
 * final, total radiation field in the IC process.
 */

/**
 * Remove the target field \a i cleaning all the related lookups.
 * Including those related to the sum of the fields.
 *
 * To be used with care because this function does not touch the
 * total counter of fields Radiation::RADFIELD_COUNTER
 *
 * @param i : index of target photon to be cleared
 */
void Radiation::ClearTargetPhotonField(int i) {
    if(i==-2) {
        if (VERBOSEMODE == true) {
            cout << "Clearing sum of all photon target field spectra..."
                 << endl;
        }
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
                 << endl;
        }
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
                 << endl;
        }
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

/**
 * Uses the method Radiation::SetThermalTargetPhotons.
 * The function updates the counter of radiation fields 
 * updating the variable Radiation::RADFIELD_COUNTER 
 */
void Radiation::AddThermalTargetPhotons(double T, double edens, int steps) {
    SetThermalTargetPhotons(T,edens,steps,RADFIELD_COUNTER);
    RADFIELD_COUNTER++;
    return;
}


/**
 * Uses the method Radiation::SetThermalTargetPhotons
 * and updates the field \a i with the new parameters.
 */
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

/** 
 * Add a greybody distribution of target photons
 * calling Radiation::GreyBody and Radiation::SetTargetPhotonVectorLookup,
 * which is used in the
 * IC emission process in this class, but which can also
 * be used in the Particles
 * class to calculate IC cooling losses
 * 
 * \param T = temperature [K]
 * \param edens = energy density in erg/cm3
 * \param steps = number of steps in the energy domain
 * \param i = counter for the photon field
 * 
 * The boundaries of the energy of the target photon field
 * are hard-coded to go from 1e-12 to 1e6 in units of (kb*T)
 * where kb is the Boltzmann constant.
 * 
 * The energy and the photon number density are passed as log10
 * of the quantity
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

/** 
 * Calls the method Radiation::SetArbitraryTargetPhotons
 * to add an arbitrary target photon to the vector of target photons.
 * The function updates the number of total fields via the variable
 * Radiation::RADFIELD_COUNTER
 * 
 * The field is used to calculate the IC emission process in this class,
 * but which can also be used 'Particles' class to calculate
 * IC cooling losses.
 * 
 * \param PhotonArray = vector of tuples (energy,number density) with units
 *                      \li energy -> erg
 *                      \li number density -> erg^-1 cm^-3
 * 
 */
void Radiation::AddArbitraryTargetPhotons(vector<vector<double> > PhotonArray) {
    SetArbitraryTargetPhotons(PhotonArray,RADFIELD_COUNTER);
    RADFIELD_COUNTER++;
    return;
}

/**
 * Uses the method Radiation::SetArbitraryTargetPhotons
 * and updates the field \a i with the new parameters.
 */
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

/**
 * Input is a vector of tuple of format (E,photon density)
 * in unitsE(erg) and photon_density(erg^-1cm^-3)
 * 
 * The values are passed to Radiation::SetTargetPhotonVectorLookup using
 * log10 of the quantities
 */
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

/** 
 * Calls the method Radiation::SetTargetPhotonsFromFile
 * to add an arbitrary target photon to the vector of target photons.
 * The function updates the number of total fields via the variable
 * Radiation::RADFIELD_COUNTER
 * 
 * Format of quantities in the file:
 * 
 * E(eV) vs. photon_density(eV^-1cm^-3)
 */
void Radiation::ImportTargetPhotonsFromFile(const char *phFile) {
    SetTargetPhotonsFromFile(phFile,RADFIELD_COUNTER);
    RADFIELD_COUNTER++;
    return;
}

/**
 * @short Reset photon field \a i with one imported from file
 *
 * Format of the file is:
 *  - E(eV) vs. photon_density(eV^-1cm^-3)
 *
 * See also Radiation::SetTargetPhotonsFromFile
 *
 * @param i : index of the field to replace
 * @param phFile : input file
 */
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


/**
 * Add an arbitrary target photon field from an ASCII file
 * 
 * Format: E(eV) vs. photon_density(eV^-1cm^-3)
 * 
 * Use of eV, as this is what is typically used in the literature.
 * 
 * The function calls back Radiation::SetTargetPhotonVectorLookup
 * 
 * \param phFile = filname
 * \param i = index for the target photon
 */
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

/** 
 * This function calls the Synchrotron code in this class.
 * The photons will be added to Radiation::TotalTargetPhotonVector
 * It uses Atoyan&Aharonian1996: MNRAS, Volume 278, Issue 2, pp. 525-541
 * 
 * The function calls the private method Radiation::SetSSCTargetPhotons
 * passing the parameters:
 * \param R = Size of the emission region (in units of parsec)
 * \param steps = steps in the energy domain for the synchrotron spectrum 
 * 
 * The field is added at the end of the list of photon fields
 */
void Radiation::AddSSCTargetPhotons(double R, int steps) {
    SetSSCTargetPhotons(R,steps,RADFIELD_COUNTER);
    RADFIELD_COUNTER++;
    return;
}

/**
 * Replace photon field \a i with the SSC field. See Radiation::AddSSCTargetPhotons
 *
 * @param i : index of the photon field to replace
 * @param R : size of the emission region (in units of parsecs)
 * @param steps : steps in the energy domain for the synchrotron spectrum
 */
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

/**
 * Set the target photon to be the Synchrotron photons produced
 * in the interaction with the magnetic field.
 *
 * Simple assumption that
 * the producing the synchrotron radiation does not affect too much
 * the particle spectrum so that the SSC uses the same electron distribution
 *
 * @param R : size of the emission regio (units of parsec)
 * @param steps : number of steps to sample the synchrotron spectrum
 * @param i : field number
 */
void Radiation::SetSSCTargetPhotons(double R, int steps, int i) {
  if (R <= 0.) {
    cout
        << "Radiation::SetSSCTargetPhotons: Source extension is <= 0... exiting!"
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
  // Factor U to account that the radiation is produced from the whole sphere
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

/**
 * Returns the sum of the target photon fields as vector of 2D tuples
 * (log10(E/(1 erg)), log10((dN/dE)/(erg-1/cm-3)))
 *
 * @param bins : bins for sampling the photon target field
 * @param ISO : true to account only for the isotropic fields. false
 *              to account for all the fields.
 *
 * @return vector of 2D tuple with the sum of the fields
 *
 * NOTE: this function could be private
 */
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



/** 
 * Function that adds up all individual target photon contributions into
 * Radiation::TargetPhotonVectors (or Radiation::TargetPhotonVectorSumIso,
 * or Radiation::TargetPhotonVectorSumAll),
 * which is what is then used by the code in the end.
 * 
 * Arguments:
 * \arg \a v = vector of the target photon field
 * \arg \a i = Target field being considered
 * 
 * If \a i = -1, the function works on the sum of all the 
 * ispotropic target photon fields
 * 
 * If \a i = -2, the function works on the sum of all fields.
 * 
 * The target field contains the photon field as log10 of the energy
 * and log10 of the photon number density.
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

// TODO: Remove this function? Keep it and make it do something?
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

/* remove the latest component in
 * TotalTargetPhotonVector and recompute
 * the total target photon spectrum
 */
//void Radiation::RemoveLastICTargetPhotonComponent() {
//  TargetPhotonVector = TargetPhotonVectorOld;
//  SetTargetPhotonVectorLookup();
//  return;
//}

/**
 * Set anisotropy for photon field \a i. The anisotropy is passed can be passed
 * as a python mesh grid. This function assumes that the electrons have a certain angle
 * with respect to the observer.
 *
 * See http://libgamera.github.io/GAMERA/docs/inverse_compton.html for a schematic
 * view of the geometry
 *
 * @param i : index of the target field
 * @param obs_angle : angle between the electron beam and x axis (observer)
 *                    as a vector of coordinate (phi, theta)
 * @param phi : vector of the phi coordinate of the anisotropy mesh grid of
 *              the target photon field
 * @param theta : vector of the theta coordinates of the anisotropy mesh grid
 *                of the target photon field
 * @param mesh : mesh grid of the anisotropy distribution of the target field
 */
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
    // FillCosZetaLookup(i);
    ANISOTROPY[i] = true;
    
    if ( distance == 0.0 ){
        SetDistance(1./pc_to_cm);
        if ( !QUIETMODE ) 
            cout << "\nRadiation::SetTargetPhotonAnisotropy: "
            		"So far no distance specified. Since an "
            		"anisotropic photon field was defined, "
            		"not the total luminosity but the radiation "
            		"in the observation direction is calculated."<<endl;
    }
    return;
}

/**
 * Set anisotropy for photon field \a i. The anisotropy is passed can be passed
 * as a python mesh grid. This function assumes that the electrons have a
 * isotropic distribution. This function sets Radiation::ISOTROPIC_ELECTRONS
 * to true. This setting is not reversible. Do not mix isotropic and anisotropic
 * distribuition of electrons.
 *
 * See http://libgamera.github.io/GAMERA/docs/inverse_compton.html for details.
 *
 * @param i : index of the target field
 * @param phi : vector of the phi coordinate of the anisotropy mesh grid of
 *              the target photon field
 * @param theta : vector of the theta coordinates of the anisotropy mesh grid
 *                of the target photon field
 * @param mesh : mesh grid of the anisotropy distribution of the target field
 */
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
            cout << "\nRadiation::SetTargetPhotonAnisotropy: "
            		"So far no distance specified. Since an "
            		"anisotropic photon field was defined, "
            		"not the total luminosity but the radiation "
            		"in the observation direction is calculated."<<endl;
    }

    return;
}



/**
 *  Function to set an isotropic electron distribution for the case of
 *  anisotropic inverse Compton scattering
 */
void Radiation::SetElectronsIsotropic(void){
    ISOTROPIC_ELECTRONS = true;
}

/**
 * Returns the anisotropy map for the target photon \a i in the phormat of a
 * mash grid, given the vectors of angular coordinates phi and theta.
 *
 * See http://libgamera.github.io/GAMERA/docs/inverse_compton.html for details.
 *
 * @param i : index of target field
 * @param phi : vector phi coordinates
 * @param theta : vector phi coordinates
 *
 * @return mesh grid style 2D array given the arrays of phi and theta
 */
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

/**
 * Calculate differential photon spectra for the different radiation processes.
 *  Spectral points will be calculated for energy points given in 'points'.
 *  These points have to be in units of 'erg'!
 *  They are stored in the 2D 'diffspec' vector and can be accessed via the
 *  Radiation::ReturnDifferentialSpectrum() member function.
 * TODO: Are the values Emin and Emax used at all, or just a relic from Joachim? If not, delete them.
 */
void Radiation::CalculateDifferentialPhotonSpectrum(vector<double> points) {
  if (!points.size()) {
    cout << "Radiation::ReturnDifferentialSpectrum: Emin>Emax! Check your "
            "boundaries. Exiting..." << endl;
    return;
  }
  if (diffSpec.size()) fUtils->Clear2DVector(diffSpec);
  if (diffSpecICComponents.size()) fUtils->Clear2DVector(diffSpecICComponents);
  if (!ElectronVector.size() && !ProtonVector.size() && !HadronSpectra.size()) {
    cout << "Radiation::ReturnDifferentialSpectrum: No particle spectra filled "
            "-> No gamma spectra to calculate. Exiting..." << endl;
    return;
  }
  
  if (diffSpecHadronComponents.size()) fUtils->Clear2DVector(diffSpecHadronComponents);
  if( !AmbientMediumComposition.size()){
      vector<double> temp0; temp0.resize(2);
      vector<vector<double> > temp;
      if (!n) { temp0[0] = 1.0; temp0[1] = 0.0;
          temp.push_back(temp0);
      }
      else { temp0[0] = 1.0; temp0[1] = n;
          temp.push_back(temp0);
      }
      SetAmbientMediumComposition(temp);
  }
  
  if (!QUIETMODE) {
    cout << "_________________________________________" << endl;
    cout << ">> CALCULATING SED FROM PARENT PARTICLES " << endl;
  }

  double ICVal, SynchVal, BremsVal, ppVal, E, Emin, Emax;
  vector<double> hadronVal;
  if (!ProtonVector.size()) {
    Emax = ElectronVector[ElectronVector.size() - 1][0];
    Emin = ElectronVector[0][0] * 1.e-6;
    if (BField)
      Emin *= 1.e-10;  // because in this case we have to go to radio energies
                       // (synchrotron)
  } else if (!ElectronVector.size()) {
    Emax = ProtonVector[ProtonVector.size() - 1][0];
    Emin = ProtonVector[0][0] * 1.e-6;      //DANGER: Does not include Hadrons!!!
  } else {
	//TODO: Adjust this to take into account all Hadron crontributions!!! ->
	//      but this might be not necessary, if Emin and Emax are not used somewhere else
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
  
  // Loop over all the energy points
  for (int i = 0; i < size; i++) {
    if (QUIETMODE == false) {
      cout << "\r";
      cout << "    "
           << i+1 << " / " << size
           << " points calculated" << std::flush;
    }
    E = points[i];

    ICVal = SynchVal = BremsVal = ppVal = 0.;
    hadronVal.clear();

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
    if (HadronSpectra.size()) {
      CalculateDifferentialGammaEmission(E, 2);
      hadronVal = GetDifferentialHadronFlux();
    } else
        hadronVal.push_back(0.);
    double total_hadron_emission = 0.;
    
    // Calculate the total contribution from hadronic interactions
    for(int j = 0; j < (int)hadronVal.size(); j++){
        total_hadron_emission += hadronVal[j];
    }

    
    diffSpec.push_back(vector<double>());
    diffSpec[diffSpec.size() - 1].push_back(E);
    diffSpec[diffSpec.size() - 1]
        .push_back(ppVal + ICVal + BremsVal + SynchVal + total_hadron_emission);
    diffSpec[diffSpec.size() - 1].push_back(ppVal);
    diffSpec[diffSpec.size() - 1].push_back(ICVal);
    diffSpec[diffSpec.size() - 1].push_back(BremsVal);
    diffSpec[diffSpec.size() - 1].push_back(SynchVal);
    diffSpec[diffSpec.size() - 1].push_back(total_hadron_emission);
    
    // Store the different values of the Hadronic interactions
    diffSpecHadronComponents.push_back(vector<double>());
    diffSpecHadronComponents[diffSpecHadronComponents.size() - 1].push_back(E);
    for(int j=0; j < (int)hadronVal.size(); j++) {
        diffSpecHadronComponents[diffSpecHadronComponents.size() - 1].push_back(hadronVal[j]);
    }    

    
    
    diffSpecICComponents.push_back(vector<double>());
    diffSpecICComponents[diffSpecICComponents.size() - 1].push_back(E);
    for(unsigned int j=0; j<RADFIELDS_MAX; j++) {
        diffSpecICComponents[diffSpecICComponents.size() - 1].push_back(fdiffics[j]);
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
 * between the energies \a emin and \a emax. The argument \a i determines the spectral
 * component to return:
 * - i = 1 : total spectrum
 * - i = 2 : pi0 decay
 * - i = 3 : IC scattering
 * - i = 4 : Bremsstrahlung
 * - i = 5 : Synchrotron radiation
 * - i - 6 : hadronic spectrum
 *
 * TODO: might need to be private. Needs the spectral lookup as argument
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

/**
 * Return the photon SED (EdN/dE (erg/s/cm^2) vs E (TeV) in a 2D-Vector of
 *  i = 1 : total spectrum
 *  i = 2 : pi0 decay
 *  i = 3 : IC scattering
 *  i = 4 : Bremsstrahlung
 *  i = 5 : Synchrotron radiation
 *  i = 6 : hadronic spectrum
 *
 *  TODO: might need to be private. Needs the spectral lookup as argument
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

/**
 * Returns the Inverse Compton spectrum for target field \a i
 * as vector of 2D tuples (E, dN/dE) in units of erg vs (1/cm2/s/erg)
 *
 * @param i : index of the target field for which to return the IC component
 * @param emin : minimum energy in erg
 * @param emax : maximum energy in erg
 * @return Inverse compton differential spectrum
 *
 * Wrapped around Radiation::ReturnDifferentialPhotonSpectrum
 */
vector<vector<double> > Radiation::GetICSpectrum(unsigned int i, double emin, double emax) {

    if(IC_CALCULATED && FASTMODE_IC==true && epoints_temp.size()) {
        FASTMODE_IC=false;
        CalculateDifferentialPhotonSpectrum(epoints_temp);
        FASTMODE_IC=true;
    }
    return ReturnDifferentialPhotonSpectrum(i+1, emin, emax, diffSpecICComponents);
}

/**
 * Returns the Inverse Compton SED for target field \a i
 * as vector of 2D tuples (E, E^2dN/dE) in units of TeV vs (erg/cm2/s)
 *
 * @param i : index of the target field for which to return the IC component
 * @param emin : minimum energy in TeV
 * @param emax : maximum energy in erg
 * @return Inverse compton SED
 *
 * Wrapped around Radiation::ReturnSED
 */
vector<vector<double> > Radiation::GetICSED(unsigned int i, double emin, double emax) {
    if(IC_CALCULATED && FASTMODE_IC==true && epoints_temp.size()) {
        FASTMODE_IC=false;
        CalculateDifferentialPhotonSpectrum(epoints_temp);
        FASTMODE_IC=true;
    }
    return ReturnSED(i+1, emin, emax, diffSpecICComponents);
  }

  
  
  
/**
 * Return pi0 decay spectrum from Hadron component \a i
 *
 * @param i : index of hadronic component
 * @param emin : minimum energy in erg
 * @param emax : maximum energy in erg
 *
 * @return Differential spectrum as vector of tuple E,dN/dE (erg, 1/erg/cm2/s)
 *
 * Wrapped around Radiation::ReturnDifferentialPhotonSpectrum
 */
vector<vector<double> > Radiation::GetHadronSpectrum(unsigned int i, double emin, double emax) {
    return ReturnDifferentialPhotonSpectrum(i+1, emin, emax, diffSpecHadronComponents);
}

/**
 * Return pi0 decay SED from Hadron component \a i
 *
 * @param i : index of hadronic component
 * @param emin : minimum energy in TeV
 * @param emax : maximum energy in TeV
 *
 * @return Differential spectrum as vector of tuple E,E^2dN/dE (TeV, 1/erg/cm2/s)
 *
 * Wrapped around Radiation::ReturnSED
 */
vector<vector<double> > Radiation::GetHadronSED(unsigned int i, double emin, double emax) {
    return ReturnSED(i+1, emin, emax, diffSpecHadronComponents);
  }
  
  
/**
 * \brief Return a particle SED E^2dN/dE vs E (erg vs TeV)
 * 
 * The particle vector is the one given by
 * Radiation::ElectronVector or Radiation::ProtonVector
 * 
 * @param type = "electron" or "proton"
 *
 * @return vector of 2D tuple (E, E^2dN/dE) ([TeV], [erg])
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


/**
 * Core method for the integrated spectrum.
 * Performs the integral between \a emin and \a emax. If ENERGYFLUX == true
 * then the integrated energy flux is computed
 * 
 * The parameter \a i is to specify which radiation mechanism to consider
 * \li \a i = 1 total flux
 * \li \a i = 2 pp interaction
 * \li \a i = 3 IC mechanism
 * \li \a i = 4 Synchrotron process
 * \li \a i = 5 Bremmstrahlung
 */
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



vector<vector<double> > Radiation::ReturnNeutrinoSpectrum(vector< vector<double> > vec){
 int size = (int)vec.size();
 if (size == 0){
     cout << "Warning in Radiation::ReturnNeutrinoSpectrum: No neutrinos calculated. "
            "Returning empty vector. Add missing protons or hadrons and/or fill it via "
            "Radiation::CalculateNeutrinoSpectrum() first!" << endl;
 }
 return vec;    
}

vector<vector<double> > Radiation::ReturnNeutrinoSED(vector< vector<double> > vec){
 int size = (int)vec.size();
 vector<vector<double> > tempVec;
 if (size == 0){
     cout << "Warning in Radiation::ReturnNeutrinoSED: No neutrinos calculated. "
            "Returning empty vector. Add missing protons or hadrons and/or fill it via "
            "Radiation::CalculateNeutrinoSpectrum() first!" << endl;
    return tempVec;
 }
 double E, E_TeV, flux, sed_value;
 for ( int i = 0; i < size; i++ ){
     E = vec[i][0];
     flux = vec[i][1];
     E_TeV = E*erg_to_TeV;
     sed_value = flux*E*E;
     fUtils->TwoDVectorPushBack(E_TeV, sed_value, tempVec);
 }
 return tempVec;
}


vector<vector<double> > Radiation::ReturnHadronNeutrinoSpectrum(int i, vector< vector< vector<double> > > vec){
    vector<vector<double> > tempvec;
 if (i >= (int)vec.size()){
     cout << "In Radiation::ReturnHadronNeutrinoSpectrum: Index value is to high, returning empty vector. Please define more hadron species with the function Radiation::AddHadrons and/or fill it via Radiation::CalculateNeutrinoSpectrum() first!" << endl;
     return tempvec;
 }
 tempvec = ReturnNeutrinoSpectrum(vec[i]);
 return tempvec;
}

vector<vector<double> > Radiation::ReturnHadronNeutrinoSED(int i, vector< vector< vector<double> > > vec){
    vector<vector<double> > tempvec;
 if (i >= (int)vec.size()){
     cout << "In Radiation::ReturnHadronNeutrinoSED: Index value is to high, returning empty vector. please define more hadron species with the function Radiation::AddHadrons." << endl;
     return tempvec;
 }
 tempvec = ReturnNeutrinoSED(vec[i]);
 return tempvec;
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
	// Now throws error when the field has not been specified
	// Otherwise crash without explanation given
    if(i<-1 || i >= (int)RADFIELD_COUNTER) {
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


/**
 * Set the photon field size for the i-th field
 * size to be passed in units of pc
 */
void Radiation::SetSizePhotonField(int i,double size){
	sizephfield[i] = size * pc_to_cm;
}


/**
 * Returns the size for all the photon fields introduced.
 * Photon fields sizes returned in units of cm
 * The function returns only non 0 values of the field sizes
 */
vector <double> Radiation::GetSizePhotonField(){
	vector <double> tempVec;
    for (unsigned int i=0; i<sizephfield.size(); i++){
        if (sizephfield[i]==0.){continue;}
            std::cout << "Photon field "<<i<<", size: "<< sizephfield[i]<<" cm"<<std::endl;
            tempVec.push_back(sizephfield[i]);
        }
    return tempVec;
}

/**
 * Return size photon field only for the i-th field.
 *
 * @param i : target photon fiels
 *
 * @return size of the photon field \a i [cm]
 */
double Radiation::GetSizePhotonField(int i){
	if (sizephfield[i]==0.){
		cout<<"Photon field size not set or set to 0!"<<endl;
	}
	return sizephfield[i];
}

/**
 * Clears the vector of photon field sizes
 */
void Radiation::ClearPhotonFieldSize(){
    if (!sizephfield.size()){
          std::cout << "Already empty" <<std::endl;
      }
      sizephfield.clear();
      // reinitialization to 0 of the vector
      // This is to avoid problems with the conditional statements
      for (unsigned int i=0;i<RADFIELDS_MAX;i++){
    	  sizephfield[i]=0.;
      }
      return;
}

/**
 * Sets the spatial dependency of the photon field along the line of sight.
 * The user provides a vector of tuples composed of 2 elements:
 *  * distance from source along the line of sight (in parsecs)
 *  * fractional value of the intensity value given when setting the photon field
 *
 *  @param i = index of the target photon index
 *  @param SpDp = Vector of the spatial dependency
 *
 *  The function sets the boolean variable of the spatial dependency to true
 *  So to properly compute the integral of the absorption coefficient along
 *  the line of sight.
 *
 *  The function converts the spatial quantity given in parsec to a quantity in cm
 */
void Radiation::SetTargetFieldSpatialDep(int i,vector< vector<double> > SpDp){
	SPATIALDEP_CURRENT = true;
	SPATIALDEP[i] = true;
    vector <vector<double> > tempVec;
    tempVec = SpDp;
	for (unsigned int j = 0;j<SpDp.size();j++){
	    tempVec[j][0] = SpDp[j][0]*pc_to_cm;
	}
    SpatialDep[i] = tempVec;
	return;
}



/**
 *  Return the spatial dependence of the target field \a i
 *  The return gives the spatial size in cm
 *  If the target does not exist or has no spatial dependency
 *  returns an empty vector
 *
 *  @param i = target photon field
 *
 *  @return Spatial dependency array ( r,normalization(r) ) distance in units of cm
*/
vector< vector<double> > Radiation::GetTargetFieldSPatialDep(int i){
	if (SpatialDep[i].size()){
		cout<<"Spatial dependence for photon field "<<i<<" not set\n"
		    <<"Return empty vector!"<<endl;
	}
	return SpatialDep[i];
}

/**
 * Function for the gamma gamma absorption cross section
 * Analytical average over solid angle
 * Eq. 5 from (Eungwanichayapant & Aharonian, 2009) https://arxiv.org/pdf/0907.2971.pdf
 * This is approximated and good within 3%
 * 
 * @param Eph1 = energy of the first photon (in ergs)
 * @param Eph2 = energy of the second photon (in ergs)
 *
 * @return the averaged gamma-gamma cross section [cm^2]
 * !!!UNDER CONSTRUCTION!!!
 * 
 * To be consistent, the energy of the photons must be given in ergs
 * 
 */
double Radiation::AverageSigmaGammaGamma(double Eph1, double Eph2) {
  double CMene = Eph1*Eph2/(m_e*m_e);
  if (CMene < 1.){
      //std::cout << "ERROR, you'll get a negative number in the square root!\n"
      //  "You are below threshold for pair production\n";
      return 0;}
  return 3./(2.*CMene*CMene)*sigma_T*((CMene+0.5*log(CMene)-1./6.+1./(2.*CMene))
		  *log(sqrt(CMene)+sqrt(CMene-1))-(CMene+4./9.-1./(9.*CMene))*sqrt(1.-(1./CMene)));
}


/**
 * Function for the full gamma gamma absorption cross section
 * Eq. 1 from Vernetto&Lipari 2016 (https://arxiv.org/pdf/1608.01587v2.pdf)
 * 
 * @param Eph1 : energy of the first photon (in ergs)
 * @param Eph2 : energy of the second photon (in ergs)
 * @param costheta : cosine of scattering angle  (in radiands)
 *
 * @return gamma-gamma cross section [cm^2]
 */
double Radiation::SigmaGammaGamma(double Eph1,double Eph2, double costheta) {
	double CMene = 2.*Eph1*Eph2*(1-costheta)/(4.*m_e*m_e);  // Centre Mass energy
	if (CMene < 1) {
		//below threshold for pair production return a cross section = 0
		return 0;
	}
	double beta = sqrt(1.-1./CMene);  //auxiliary variable
	double crosssec = sigma_T * 3. / 16. * (1.-beta*beta) *
			            (2.*beta*(beta*beta-2.) + (3.-beta*beta*beta*beta) *
			            		log((1+beta)/(1-beta)));
	return crosssec;
}


/**
 * Auxiliary function to compute only the absorption coefficient,
 * no spatial integration.
 *
 * The function takes into account the angular anisotropy of the photon field
 * through interpolation of the mesh grid of the angular dependencies.
 * Angular integration is a simple rectangular integration.
 *
 * @param Egamma : energy of the gamma-ray photon (in erg)
 * @param target : target photon field to compute the absorption
 *
 * @return Absorption coefficient. Units of 1/cm
 */
double Radiation::ComputeAbsCoeff(double Egamma, int target) {
	vector< vector<double> > targets = Radiation::GetTargetPhotons(target);
    vector< vector<double> > TempVect;
    double product=0;
    double integrand = 0;
    double integrand2 = 0;
    double integral=0;
    double Q=0; //the scaling factor due to the angular distribution of target field
    if ( ANISOTROPY_CURRENT ){
        double cos_zeta = 0.;
        // because the photon direction is aligned with the x-axis by definition
        // hard coded so that it is not influenced by the direction of the electron beam
        double cos_kappa = 0.; double sin_kappa = 1.;
        double phi_min = (*TargetPhotonAngularBoundsCurrent)[0];
        double phi_max = (*TargetPhotonAngularBoundsCurrent)[1];
        double theta_min = (*TargetPhotonAngularBoundsCurrent)[2];
        double theta_max = (*TargetPhotonAngularBoundsCurrent)[3];
        for (unsigned int j=0;j<targets.size();j++){
            integrand2 = 0;
        	for (double phi = phi_min; phi <= phi_max; phi += d_phi) {
        		integrand = 0;
        		for (double theta = theta_min; theta <= theta_max; theta += d_theta) {
        			// simplified cos_zeta: the photon direction is always along x-axis
        			cos_zeta = -cos_kappa *cos(theta) + sin_kappa * sin(theta) * cos(phi);

        			Q = interp2d_spline_eval(*TargetPhotonAngularDistrCurrent,
        			                                     phi, theta, *phiaccescCurrent,*thetaaccescCurrent);
        			product = Radiation::SigmaGammaGamma(Egamma,targets[j][0],cos_zeta)*targets[j][1]*Q*(1-cos_zeta);
        			integrand += product*sin(theta)*d_theta;
        		}
                integrand2 += integrand*d_phi;
        	}
        	fUtils->TwoDVectorPushBack(targets[j][0],integrand2,TempVect);
        }
    }
    else {
    	for (unsigned int j=0;j<targets.size();j++){
    	        product = Radiation::AverageSigmaGammaGamma(Egamma,targets[j][0])*targets[j][1];
    	        fUtils->TwoDVectorPushBack(targets[j][0],product,TempVect);
    	    }
    }
    integral = fUtils->Integrate(TempVect,targets[0][0],targets[targets.size()-1][0]);  //Make this process smarter
    return integral;  // returned value is in units of 1/cm
}


/**
 * Functions for the calculation of the optical depth
 * Takes as argument just the energy f the gamma ray photon and computes the optical
 * depth for the chosen target photon
 *
 * @param Egamma : energy of the gamma ray [erg]
 * @param target : index of the photon field target
 * @param phsize : patial size of the photon field [cm]
 *
 * For the non-homogeneous field case, it assumes that Radiation::SpatialDep
 * has been initialized via Radiation::SetTargetFieldSpatialDep
 *
 * !!! CAREFUL UNDER CONSTRUCTION !!!
 * !!! Implemented simple spatial variability of the target photons
 * !!! The SpatialDep array is a scaling relation.
 */
double Radiation::ComputeOptDepth(double Egamma, int target, double phsize){
	double tauval = 0;
    double integral =0;
    integral = ComputeAbsCoeff(Egamma,target);
    if (SPATIALDEP[target]) {
    	// integrate over the real spatial dependence along the line of sight
    	// TEST: Use the simple multiplication (then will use the lookup)
    	// Assume a fixed step in the definition of the spatial range
    	double deltar = SpatialDep[target][1][0]-SpatialDep[target][0][0];
    	for (unsigned int i=0;i<SpatialDep[target].size();i++){
    	    tauval += (integral*SpatialDep[target][i][1]*deltar);
    	}
    }
    else {
            tauval = (integral*phsize); // assumes a homogeneous field.
        }
    return tauval;
}


/**
 * Functions for the calculation of the optical depth
 * Takes as argument just the energy f the gamma ray photon and computes the optical
 * depth for the sum of the target photons that are present
 * At the moment it works only for isotropic and homogeneous case
 *
 * @param Egamma : energy of the gamma ray [erg]
 * @param target : index of the photon field target
 * @param phsize : spatial size of the photon field [cm]
 */
double Radiation::ComputeOptDepthIsotropic(double Egamma, int target, double phsize){
    double tauval=0;
    vector< vector<double> > targets = Radiation::GetTargetPhotons(target);
    vector< vector<double> > TempVect;
    double product=0;
    double integral=0;
    for (unsigned int j=0;j<targets.size();j++){
        product = Radiation::AverageSigmaGammaGamma(Egamma,targets[j][0])*targets[j][1];
        fUtils->TwoDVectorPushBack(targets[j][0],product,TempVect);
    }
    integral = fUtils->Integrate(TempVect,targets[0][0],targets[targets.size()-1][0]);  //Make this process smarter
    tauval= (integral*phsize); //assumes a homogeneous field. Would need another integral if the density varies
    return tauval;
}


/**
 * Return the total absorbed differential spectrum (for the moment using the homogeneous and isotropic case)
 * This is a wrapper around the function ReturnDifferentialPhotonSpectrum to which it was added the distance
 * parameter needed to compute the right absorption and the target photon field for the absorption
 * 
 * @param emin : minimum energy [erg]
 * @param emax : maximum energy [erg]
 * @param fields : vector photon field counter to be used in the order the fields were added
 * @param size : size of the photon field integration [cm]
 *
 * @return Absorbed Spectrum : vector of 2D tuple (E, dN/dE) ([erg],[1/erg-1/cm2/s])
 */
vector< vector<double> > Radiation::ReturnAbsorbedSpectrumOnFields(double emin, double emax, vector <int> fields, vector <double> size) {
    double tauval=0;
    double singletau;
    vector< vector<double> > tempVec;
    //returns the total Differential Photon Spectrum
    tempVec = Radiation::ReturnDifferentialPhotonSpectrum(1, emin, emax, diffSpec);
    if (!fields.size()){
        cout<<"You have not added any photon field\n"
                 <<"Add photon fields before proceeding with this step\n"
                 <<"Returning unmodified spectrum!"<<endl;
        return tempVec;
    }
    for (unsigned int k=0; k < fields.size(); k++){
    	if (size[fields[k]]==0.) {
    	cout<<"You have not set the size of field "<<fields[k]<<endl;
    	cout<<"Add photon fields before proceeding with this step\n"
            <<"Returning unmodified spectrum!"<<endl;
    	return tempVec;
    	}
    }
    for (unsigned int j=0; j < tempVec.size(); j++){
        tauval = 0;
        for (unsigned int k=0; k < fields.size(); k++){
        	// energy units of tempVec is already in ergs
            singletau = Radiation::ComputeOptDepth(tempVec[j][0],fields[k],size[fields[k]]);
            tauval += singletau;
        }
        tempVec[j][1] = tempVec[j][1]*exp(-tauval);
    }
    return tempVec;
}


/**
 * Return the total absorbed SED (for the moment using the homogenous and isotropic case)
 * This is a wrapper around the function ReturnSED to which it was added the distance
 * parameter needed to compute the right absoption and the target photon field for the absorption
 * 
 * @param emin : minimum energy [erg]
 * @param emax : maximum energy [erg]
 * @param fields : vector photon field counter to be used in the order the fields were added
 * @param size : size of the photon field integration [cm]
 *
 * @return Absorbed SED : vector of 2D tuple (E, dN/dE) ([erg],[erg/cm2/s])
 */
vector< vector<double> > Radiation::ReturnAbsorbedSEDonFields(double emin, double emax,
                                                        vector <int> fields, vector <double> size){
    double tauval=0;
    double singletau;
    vector< vector<double> > tempVec;
    tempVec = Radiation::ReturnSED(1, emin, emax, diffSpec);  //returns the total SED
    if (!fields.size()){
        std::cout<<"You have not added any photon field\n"
                 <<"Add photon fields before proceeding with this step\n"
                 <<"Returning unmodified SED!"<<std::endl;
        return tempVec;
    }
    for (unsigned int k=0; k < fields.size(); k++){
        	if (size[fields[k]]==0.) {
        	cout<<"You have not set the size of field "<<fields[k]<<endl;
        	cout<<"Add photon fields before proceeding with this step\n"
                <<"Returning unmodified spectrum!"<<endl;
        	return tempVec;
        	}
        }
    for (unsigned int j=0; j < tempVec.size(); j++){
        tauval = 0;
        for (unsigned int k=0; k < fields.size(); k++){
        	// need to convert the energy of the photon field because it is returned in TeV
            singletau = Radiation::ComputeOptDepth(tempVec[j][0]*TeV_to_erg,fields[k],size[fields[k]]);
            tauval += singletau;
        }
        tempVec[j][1] = tempVec[j][1]*exp(-tauval);
    }
    return tempVec;
}


/**
 * Returns the absorbed integrated flux. The absorption is computed using all the photon fields that enter in the
 * absorption calculation
 *
 * @param emin : minimum energy for the integration [erg]
 * @param emax : maximum energy for the integration [erg]
 * @param ENERGYFLUX : if true intgrates Edn/dE if false dn/dE
 * @param fields : vector of indices of fields used to compute the absorption
 * @param size : vector of sizes of the fields used to compute absorption (sizes in cm)
 *
 * @return val : integral flux (1/cm2/s) or integral energy flux (erg/cm2/s)
 */
double Radiation::ReturnAbsorbedIntergratedFlux(double emin, double emax, bool ENERGYFLUX, vector <int> fields, vector <double > size){
    if (!diffSpec.size()){
        cout << "Radiation::GetIntegratedFlux: Differential spectrum "
            "vector empty. Fill it via "
            "Radiation::CalculateDifferentialSpectrum() first! Returning zero."
            << endl;
        return 0.;
    }
    double tauval=0;
    double singletau=0;
    if (!emin) emin = diffSpec[0][0];
    if (!emax) emax = diffSpec[diffSpec.size()-1][0];
    vector <vector <double> > tempVec;
    double e, dNdE, val;
    for (unsigned int j = 0; j < diffSpec.size(); j++) {
        tauval =0;
        e = diffSpec[j][0];
        for (unsigned int k=0; k < fields.size(); k++){
            singletau = Radiation::ComputeOptDepth(e,fields[k],size[fields[k]]);
            tauval += singletau;
        }
        dNdE = diffSpec[j][1]*exp(-tauval);
        if (dNdE <= 0. ) continue;
        if (ENERGYFLUX == false) fUtils->TwoDVectorPushBack(e,dNdE,tempVec);
        else fUtils->TwoDVectorPushBack(e,e*dNdE,tempVec);
    }
    if (tempVec.size() < 3) return 0.;
    fUtils->ToggleQuietMode();
    val = fUtils->Integrate(tempVec,emin,emax);
    fUtils->ToggleQuietMode();
    return val;
}

/*
 * Functions to directly retrieve the absorbed gamma ray spectrum and
 * the absorbed gammaray SED
 * using the above calculated optical depth
 * It makes things easier when scripting
 
vector <vector <double> > Radiation::RetAbsDiffSpectrum(double emin, double emax){
	vector <vector<double> > gammaspec = Radiation::ReturnDifferentialPhotonSpectrum(1,emin,emax,diffSpecICComponents);
    vector <vector<double> > TempVect;
    for (unsigned int i=0;i<gammaspec.size();i++){
        double tau = Radiation::ComputeOptDepth(gammaspec[i][0]);
    	fUtils->TwoDVectorPushBack(gammaspec[i][0],gammaspec[i][1]*exp(-tau),TempVect);
    }
    return TempVect;
}

vector <vector <double> > Radiation::RetAbsSED(double emin, double emax){
	vector <vector<double> > gammaspec = Radiation::ReturnSED(1,emin,emax,diffSpecICComponents);
    vector <vector<double> > TempVect;
    for (unsigned int i=0;i<gammaspec.size();i++){
        double tau = Radiation::ComputeOptDepth(gammaspec[i][0]);
    	fUtils->TwoDVectorPushBack(gammaspec[i][0],gammaspec[i][1]*exp(-tau),TempVect);
    }
    return TempVect;
}
*/





/**
 * Function to set the ambient medium composition to the one in the local galactic ISM
 * composition. It uses the same relative abundances as in Kafexhiu et al. 2014
 * 
 * \param density: Density of hydrogen in 1/cm^3
 * **/
void Radiation::SetLocalAmbientMediumComposition(double density, string norm){
    vector < double > mass_numbers = {1.0,4.0,12.0,14.0,16.0,20.0,24.0,28.0,32.0,56.0};
    double normalisation;
    vector < double > abundances = {1.0, 9.59e-2, 4.65e-4,8.3e-5, 8.3e-4,1.2e-4, 3.87e-5, 3.69e-5, 1.59e-5,3.25e-5};
    if (norm == "Hydrogen" || norm == "hydrogen"){
        //vector < double > abundances = {1.0*density, 9.59e-2*density, 4.65e-4*density,8.3e-5*density, 8.3e-4*density,1.2e-4*density, 3.87e-5*density, 3.69e-5*density, 1.59e-5*density,3.25e-5*density};
        normalisation = 1.0;
        //for(int i=0;i<abundances.size(); i++){
        //    abundances[i] *= density;
        //}
    }
    
    else if (norm == "number density"){
        normalisation = 0.0;
        for(unsigned int i=0;i<abundances.size(); i++){
            normalisation += abundances[i];
        }
    }
    
    else if (norm == "mass density"){
        normalisation = 0.0;
        for(unsigned int i=0;i<abundances.size(); i++){
            normalisation += (abundances[i]*mass_numbers[i]);
        }        
    }
    
    else{
     cout << "In function Radiation::SetLocalAmbientMediumComposition: ";
     cout << "No valid input string for parameter 'norm'.\n";
     cout << "Options are:  'hydrogen' (default)\n";
     cout << "              'number density'\n";
     cout << "              'mass density'\n";
     cout << "Exiting without setting the ambient medium composition.\n";
     return;
    }
    
    double density_normed = density/normalisation;
    for(unsigned int i=0;i<abundances.size(); i++){
        abundances[i] *= density_normed;
    }
    
    vector< vector <double> > zipped_vector = fUtils->ZipTwoOneDVectors(mass_numbers, abundances);
    SetAmbientMediumComposition(zipped_vector);
    
    return;
}





/**
 * Function to use the local cosmic ray model above 10 TeV by ... TODO: Add paper reference when it is published
 * 
 * \param extragalactic (optional): Boolean, if true use the extragalctic component,
 *                                  if false not.
 * **/
void Radiation::AddLocalCosmicRayModel(bool extragalactic, double slope){
    double N_H=1.0e-4, N_He=1.22e-4, N_C=2.45e-5, N_O=3.5e-5, N_Mg=0.87e-5, N_Si=1.2e-5, N_Fe=2.17e-5, alpha=2.65, ecut=3.0e6, alpha2=2.4, N_downscale=1.8e-2, ecut2 = 70.0e6, N_extragal = 2.5e-8, ecut_extragal = 0.45e11, alpha_extragal = 2.3;
    
    alpha = slope;
    
    double Ns[] = {N_H, N_He, N_C, N_O, N_Mg, N_Si, N_Fe};
    double mass_numbers_CR[] = {1.0,4.0,12.0, 16.0, 24.0, 28.0, 56.0, 1.0};
    double charges[] = {1.0,2.0,6.0,8.0,12.0,14.0,26.0};
    double ecuts[] = {1.0*ecut,2.0*ecut,6.0*ecut,8.0*ecut,12.0*ecut,14.0*ecut,26.0*ecut};
    
    vector < double > e_values;
    e_values.push_back(2.0);                    // These values are log10(E/GeV)
    double difference = (11.5 - 2.0)/300.0;
    for (int i = 0; i < 300; i++) {
        e_values.push_back(e_values[i] + difference);
    }
    
    
    double conversion = 4.0*pi/1.0e4/c_speed;
    
    // Here, the Hadron species of the first galactic component are initialized
    double N, mass, charge, ecut_value, e;
    vector < vector<double> > tempvec;
    for (int i = 0; i<7; i++) {
        tempvec.clear();
        N = Ns[i]; mass = mass_numbers_CR[i]; charge = charges[i]; ecut_value = ecuts[i];
        for (int j = 0; j< (int) e_values.size(); j++){
            e = pow(10.0,e_values[j]);
            fUtils->TwoDVectorPushBack(e*GeV_to_erg, N*pow(e/1.0e3, -1.0*alpha) * exp(-e/ecut_value)/GeV_to_erg * conversion, tempvec);
        }
        
        AddHadrons(tempvec, mass);
    }
    
    
    // Here, the Hadron species of the second galactic component are initialized
    for (int i = 0; i<7; i++) {
        tempvec.clear();
        N = Ns[i]; mass = mass_numbers_CR[i]; charge = charges[i]; ecut_value = charge*ecut2;
        for (int j = 0; j< (int) e_values.size(); j++){
            e = pow(10.0,e_values[j]);
            fUtils->TwoDVectorPushBack(e*GeV_to_erg, N*N_downscale*pow(e/1.0e3, -1.0*alpha2) * exp(-e/ecut_value)/GeV_to_erg*conversion, tempvec);
        }
        AddHadrons(tempvec, mass);
    }
    
    
    // If needed, add the extragalactic component. It is assumed, they are solely protons
    if(extragalactic) {
        tempvec.clear();
        for (int j = 0; j< (int) e_values.size(); j++){
            e = pow(10.0,e_values[j]);
            fUtils->TwoDVectorPushBack(e*GeV_to_erg, N_extragal*pow(e/1.0e3, -1.0*alpha_extragal) * exp(-e/ecut_extragal)/GeV_to_erg * conversion, tempvec);
        }
        AddHadrons(tempvec, 1.0);
    }
    
    return;
}


void Radiation::AddCRModel(bool extragalactic, bool CREAM){
    
    
    double mass_numbers_CR[] = {1.0,4.0,12.0, 16.0, 24.0, 28.0, 56.0, 1.0};
    double charges[] = {1.0,2.0,6.0,8.0,12.0,14.0,26.0};
    
    vector < double > e_values;
    e_values.push_back(2.0);                    // These values are log10(E/GeV)
    double difference = (11.5 - 2.0)/300.0;
    for (int i = 0; i < 300; i++) {
        e_values.push_back(e_values[i] + difference);
    }
    
    double conversion = 4.0*pi/1.0e4/c_speed;
    
    double N, mass, charge, ecut_value, e;
    vector < vector<double> > tempvec;
    
    double N0_H = 28.2, E1_H = 2.0e3, alpha1_H = 2.815, alpha2_H = 2.45, alpha3_H = 2.72, E2_H = 2.5e4, ecut_H = 9.0e6;
    if( CREAM) { N0_H = 28.2; E1_H = 2.0e3; alpha1_H = 2.815; alpha2_H = 2.4; alpha3_H = 3.0; E2_H = 1.5e4; ecut_H = 9.0e6;}
    
    double N0_He = 16.0, E1_He = 4.0e3, alpha1_He = 2.71, alpha2_He = 2.32, alpha3_He = 2.66, E2_He = 4.0e4, ecut_He = 5.0e6;
    if( CREAM ){ N0_He = 16.0; E1_He = 3.0e3; alpha1_He = 2.71; alpha2_He = 2.2; alpha3_He = 2.87; E2_He = 9.0e3; ecut_He = 5.0e6;}
    
    // Parameters only needed for H and He because of the double break spectrum:
    double alpha3s[] = {alpha3_H, alpha3_He};
    double E2s[] = {E2_H, E2_He};
    
    // Parameters used for all species:
    double N0_C = 3.5, E1_C = 7.0e3, alpha1_C = 2.7, alpha2_C = 2.45, ecut_C = 3.0e6*3.0;
    double N0_O = 6.0, E1_O = 7.0e3, alpha1_O = 2.7, alpha2_O = 2.5, ecut_O = 3.0e6*2.0;
    double N0_Mg = 3.0, E1_Mg = 1.0e4, alpha1_Mg = 2.75, alpha2_Mg = 2.55, ecut_Mg = 3.0e6*14.0;
    double N0_Si = 2.5, E1_Si = 1.0e4, alpha1_Si = 2.7, alpha2_Si = 2.5, ecut_Si = 3.0e6*13.0;
    double N0_Fe = 10.0, E1_Fe = 1.0e4, alpha1_Fe = 2.8, alpha2_Fe = 2.55, ecut_Fe = 3.0e6*12.0;
    
    double N0s[] = {N0_H, N0_He, N0_C, N0_O, N0_Mg, N0_Si, N0_Fe};
    double E1s[] = {E1_H, E1_He, E1_C, E1_O, E1_Mg, E1_Si, E1_Fe};
    double alpha1s[] = {alpha1_H, alpha1_He, alpha1_C, alpha1_O, alpha1_Mg, alpha1_Si, alpha1_Fe};
    double alpha2s[] = {alpha2_H, alpha2_He, alpha2_C, alpha2_O, alpha2_Mg, alpha2_Si, alpha2_Fe};
    double ecuts[] = {ecut_H, ecut_He, ecut_C, ecut_O, ecut_Mg, ecut_Si, ecut_Fe};
    
    
    //******************************************************************
    // First component for H and He
    double E0 = 10.0, w = 0.5;
    double alpha1, alpha2, alpha3, E1, E2;
    
    
        for (int i = 0; i<2; i++) {
        tempvec.clear();
        N = N0s[i]; mass = mass_numbers_CR[i]; charge = charges[i]; ecut_value = ecuts[i]; alpha1 = alpha1s[i]; alpha2 = alpha2s[i]; E1 = E1s[i];
        alpha3 = alpha3s[i]; E2 = E2s[i];
        for (int j = 0; j< (int) e_values.size(); j++){
            e = pow(10.0,e_values[j]);
            fUtils->TwoDVectorPushBack(e*GeV_to_erg, N * pow((e/E0), -alpha1) * pow((1.0 + (e*e/(E1*E1))), (-(alpha2-alpha1)*w)) * pow((1.0 + (e*e/(E2*E2))),(-(alpha3-alpha2)*w)) * exp(-e/ecut_value)/GeV_to_erg*conversion, tempvec);
        }
        AddHadrons(tempvec, mass);
    }
    
    
    
    //*********************************************************
    // First component for Carbon until Iron
    
    
    
    for (int i = 2; i<7; i++) {
        tempvec.clear();
        N = N0s[i]; mass = mass_numbers_CR[i]; charge = charges[i]; ecut_value = ecuts[i]; alpha1 = alpha1s[i]; alpha2 = alpha2s[i]; E1 = E1s[i];
        for (int j = 0; j< (int) e_values.size(); j++){
            e = pow(10.0,e_values[j]);
            fUtils->TwoDVectorPushBack(e*GeV_to_erg, N*pow(e/E0, -1.0*alpha1) *  pow( 1. + (e*e/(E1*E1)), -1.*(alpha2-alpha1)*w) * exp(-e/ecut_value)/GeV_to_erg*conversion, tempvec);
        }
        AddHadrons(tempvec, mass);
    }
    
    
    
    
    //**********************************************************
    // Adding the second galactic component for each species:
    
    double alpha_H2 = 2.5, alpha_He2 = 2.5, alpha_C2 = 2.45, alpha_O2 = 2.45, alpha_Mg2 = 2.45, alpha_Si2 = 2.45, alpha_Fe2 = 2.4;
    double N_H2 = 6.0e-6, N_He2 = 4.0e-6, N_C2 = 6.0e-7, N_O2 = 6.0e-7, N_Mg2 = 5.0e-7, N_Si2 = 3.0e-7, N_Fe2 = 8.0e-7;
    double ecut_H2 = 70.0e6, ecut_He2 = 70.0e6*2.0, ecut_C2 = 70.0e6*6.0, ecut_O2 = 70.0e6*8.0, ecut_Mg2 = 70.0e6*12.0, ecut_Si2 = 70.0e6*14.0, ecut_Fe2 = 70.0e6*12.0;
    
    double alphas2[] = {alpha_H2, alpha_He2, alpha_C2, alpha_O2, alpha_Mg2, alpha_Si2, alpha_Fe2};
    double Ns2[] = {N_H2, N_He2, N_C2, N_O2, N_Mg2, N_Si2, N_Fe2};
    double ecuts2[] = {ecut_H2, ecut_He2, ecut_C2, ecut_O2, ecut_Mg2, ecut_Si2, ecut_Fe2};
    
    // Here, the Hadron species of the second galactic component are initialized
    for (int i = 0; i<7; i++) {
        tempvec.clear();
        N = Ns2[i]; mass = mass_numbers_CR[i]; charge = charges[i]; ecut_value = ecuts2[i];
        for (int j = 0; j< (int) e_values.size(); j++){
            e = pow(10.0,e_values[j]);
            fUtils->TwoDVectorPushBack(e*GeV_to_erg, N*pow(e/1.0e3, -1.0*alphas2[i]) * exp(-e/ecut_value)/GeV_to_erg*conversion, tempvec);
        }
        AddHadrons(tempvec, mass);
    }
    
    
    
    
        
    // If needed, add the extragalactic component. It is assumed, they are solely protons
    double N_extragal = 1.4e-7, alpha_extragal = 2.4, ecut_extragal = 0.45e11;
    if(extragalactic) {
        tempvec.clear();
        for (int j = 0; j< (int) e_values.size(); j++){
            e = pow(10.0,e_values[j]);
            fUtils->TwoDVectorPushBack(e*GeV_to_erg, N_extragal*pow(e/1.0e3, -1.0*alpha_extragal) * exp(-e/ecut_extragal)/GeV_to_erg * conversion, tempvec);
        }
        AddHadrons(tempvec, 1.0);
    }
    
    
    
    return;
}



