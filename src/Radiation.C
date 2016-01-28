#include "Radiation.h"
/**
 * Constructor. Nothing fancy, just initialises some stuff.
 */
Radiation::Radiation() {
  /* Default values */
  fUtils = new Utils();
  ParticleVector.clear();
  DEBUG = false;
  fintbrems = lintbrems = fintpp = lintpp = fintic = lintic = 0.;
  ldiffbrems = fdiffbrems = ldiffsynch = fdiffsynch = 0.;
  ldiffic = fdiffic = ldiffpp = fdiffpp = 0.;
  TargetPhotonEdens = 0.;
  TargetPhotonVector.clear();
  TargetPhotonVectorOld.clear();
  INTEGRATEOVERGAMMAS = false;
  QUIETMODE = false;
  PiModel = 1;
  SynchModel = 0;
  integratorTolerance = 1.e-1;
  integratorKronrodRule = 2;
  acc = gsl_interp_accel_alloc();
}

/**
 * Standard destructor.
 */
Radiation::~Radiation() {}

void Radiation::Reset() {
  ParticleVector.clear();
  ElectronVector.clear();
  ProtonVector.clear();
  gsl_spline_free(ElectronLookup);
  gsl_spline_free(ProtonLookup);
  gsl_spline_free(TargetPhotonLookup);
  BField = 0.;
  n = 0.;
  fintbrems = 0.;
  lintbrems = 0.;
  fintpp = 0.;
  lintpp = 0.;
  fintic = 0.;
  lintic = 0.;
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

  /* Conversion of differential photon rate to flux
   [ph/(erg*s) -> ph/(erg*s*cm^2)] */
  double lumtoflux = 0.;
  void *p = NULL;
  if (!distance) {
    cout << "### Radiation::CalculateDifferentialGammaEmission: Distance to "
            "particles not specified -> Flux equals now the luminosity! ###"
         << endl;
    lumtoflux = 1.;
  } else
    lumtoflux = 1. / (4. * pi * distance * distance);

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
    radiationMechanism = "Bremsstrahlung";
    ldiffbrems = DifferentialEmissionComponent(e, p);
    fdiffbrems = lumtoflux * ldiffbrems;

    radiationMechanism = "InverseCompton";
    ldiffic = DifferentialEmissionComponent(e, p);
    fdiffic = lumtoflux * ldiffic;

    radiationMechanism = "Synchrotron";
    ldiffsynch = DifferentialEmissionComponent(e, p);
    fdiffsynch = lumtoflux * ldiffsynch;
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

/**
 * Calculates the integral photon Emission above energy 'e' [erg]. This gives
 * - integral fluxes, fint* [(no. of photons)/(s*cm^2)]
 * - photon rate, lint* [(no. of photons)/s]
 *
 * The * stands for Brems(-strahlung), IC and pp (i.e. pi^0 decay) components.
 * These quantities are private members.
 *
 * It calls the 'DifferentialEmissionComponent'
 * function for the relevant radiation mechanisms:
 * - protons
 *   + inelastic p-p scattering
 * - electrons
 *   + IC emission
 *   + Bremsstrahlung
 *   + Synchrotron radiation
 *
 * This function calculates the appropriate (depending on the particle species)
 * radiation mechanism and their gamma-ray flux by calling
 * the method 'CalculateLuminosityAndFlux' which integrates
 * 'DifferentialEmissionComponent'.
 */
void Radiation::CalculateIntegralGammaEmission(double e, int particletype) {

  lintbrems = fintbrems = lintic = fintic = lintpp = fintpp = 0.;

  if (particletype != 0. && particletype != 1) {
    cout << "### Radiation::CalculateIntegratedGammaEmission: Please provide "
            "proper particle spectrum and particle type! ###" << endl;
    return;
  }

  if ((!particletype && !ElectronVector.size()) ||
      (particletype == 1 && !ProtonVector.size())) {
    cout << "### Radiation::CalculateIntegratedGammaEmission: No accelerated "
            "particles! Exiting... ###" << endl;
    return;
  } else if (!particletype) {
    ParticleVector = ElectronVector;
    if (!n) {
      cout << "Radiation::CalculateIntegratedGammaEmission: No ambient density "
              "value set for Bremsstrahlung. Returning zero value." << endl;
      lintbrems = 0.;
      fintbrems = 0.;
    } else {
      CalculateLuminosityAndFlux("Bremsstrahlung", e, lintbrems, fintbrems);
    }
    CalculateLuminosityAndFlux("InverseCompton", e, lintic, fintic);
  } else if (particletype == 1) {
    ParticleVector = ProtonVector;
    CalculateLuminosityAndFlux("ppEmission", e, lintpp, fintpp);
  } else {
    cout << "### Radiation::CalculateIntegratedGammaEmission: WTF?! That is "
            "not possible!!" << endl;
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
  if (e > ParticleVector[ParticleVector.size() - 1][0]) return 0.;
  double egamma = e;
  double emax = ParticleVector[ParticleVector.size() - 1][0];
  double emin = ParticleVector[0][0];
  if(e<emin) e = emin;
  fPointer IntFunc = NULL;
  if (!radiationMechanism.compare("Synchrotron")) {
    if (!BField) {
      if(!QUIETMODE) cout << "Radiation::DifferentialEmissionComponent: No "
                             "BField value set for Synchrotron radiation. "
                             "Returning zero value." << endl;
      return 0.;
    }
    if (!SynchModel)
      IntFunc = &Radiation::SynchEmissivity;
    else if (SynchModel == 1)
      IntFunc = &Radiation::SynchEmissivityExplicit;
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
    if (!TargetPhotonVector.size()) {
      if(!QUIETMODE) cout << "Radiation::DifferentialEmissionComponent: "
                             "No radiation fields set for IC "
                             "emission. Returning zero value." << endl;
      return 0.;
    }
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
  double gammas = Integrate(IntFunc, &egamma, e, emax, integratorTolerance,
                            integratorKronrodRule);
  if (std::isnan(gammas)) return 0.;

  return gammas;
}

/**
 * Calculates
 * - f: integrated flux [(no. of photons)/(s*cm^2)] and
 * - l: photon rate [(no. of photons)/s]
 *
 * above energy 'e' [erg].
 *
 * This method integrates the 'DifferentialEmissionComponent' method
 * for the relevant radiation mechanisms:
 * - protons
 *   + inelastic p-p scattering
 * - electrons
 *   + IC emission
 *   + Bremsstrahlung
 *   + Synchrotron radiation
 *
 * The integration boundaries are {e,maximum particle energy}, where the
 * maximum particle energy is automatically determined from the
 * 'ParticleVector' vector.
 */
void Radiation::CalculateLuminosityAndFlux(string mechanism, double e,
                                           double &l, double &f) {

  double lumtoflux = 0.;
  double emingamma = e;
  /* if no distance is provided, Luminosity and Flux are treated as
   * being the same */
  if (!distance) {
    cout << "### Radiation::CalculateLuminosityAndFlux: Distance to particles "
            "not specified -> Flux equals now the luminosity! ###" << endl;
    lumtoflux = 1.;
  } else
    lumtoflux = 1. / (4. * pi * distance * distance);

  if (!ParticleVector.size()) {
    cout << "### Radiation::CalculateLuminosityAndFlux: No accelerated "
            "particles! Exiting... ###" << endl;
    l = 0.;
    f = 0.;
    return;
  }

  double emax = ParticleVector[ParticleVector.size() - 1][0];
  if (e >= 0.99999 * emax) {
    l = 0.;
    f = 0.;
    return;
  }
  radiationMechanism = mechanism;
  fPointer IntFunc = &Radiation::DifferentialEmissionComponent;

  l = Integrate(IntFunc, &emingamma, e, emax, integratorTolerance,
                integratorKronrodRule);
  f = lumtoflux * l;

  return;
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
    eelectron = *(double *)par;
    egamma = x;
  } else {
    eelectron = x;
    egamma = *(double *)par;
  }

  gsl_interp_accel_reset(acc);
  fPointer IntFunc = &Radiation::ICEmissivity;
  double xpars[2] = {eelectron, egamma};

  /* detemine integration boundaries for the target photon energy from Eq. 2.50
     in Blumenthal&Gould (Reviews of Modern Physics, vol. 42, no. 2, 1970) */
  double lorentz = eelectron / m_e;
  double edash = egamma / (1. - egamma / (lorentz * m_e));
  double k = 1. / (4. * lorentz * lorentz);
  double boundmin, boundmax;
  (edash *k < targetphotonenergymin) ? (boundmin = targetphotonenergymin)
                                     : boundmin = edash * k;
  (edash > targetphotonenergymax) ? (boundmax = targetphotonenergymax)
                                  : boundmax = edash;
  if (k > 0.1 || boundmin >= boundmax) return 0.;
  icgammas = Integrate(IntFunc, xpars, boundmin, boundmax, integratorTolerance,6);
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

  double ephoton = x;  ///< energy of the target photon
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
  double targetphotons = fUtils->EvalSpline(log10(ephoton),TargetPhotonLookup,
                                            acc,__func__,__LINE__);
  double integrand = 2. * pi * pow(e_radius, 2.) * m_e * c_speed / lorentz *
                     pow(10., targetphotons) / ephoton * bracket;

  integrand /= lorentz*m_e;
  return integrand;
}

/** return a lookup table holding the differential electron energy loss rate due
 * to inverse-Compton
 *  scattering. The format of the lookup is: { Energy(erg) - Energy Loss Rate
 * from IC scattering(erg/s) }
 */
void Radiation::CreateICLossLookup(int bins) {

  if(!TargetPhotonVector.size()) {
    cout << "Radiation::CreateICLossLookup: No target photons! Exiting." <<endl;
    return;
  }
  INTEGRATEOVERGAMMAS = true;
  fUtils->Clear2DVector(ICLossLookup);
  /* lower integration boundary over emitted (i.e. 'loss-') IC photons */
  double EGammaMin = 1.e-25 * TeV_to_erg;
  /* Upper integration boundary over emitted (i.e. 'loss-') IC photons */
  double EGammaMax = 1.e8 * TeV_to_erg;
  /* lower integration boundary for the IC loss lookup (i.e. here simply the
   * electron rest mass) */
  double logemin = log10(0.1 * m_e);
  /* upper integration boundary for the IC loss lookup */
  double logemax = log10(EGammaMax);
  double logestep = (double)(logemax - logemin) / bins;
  int tt = 1;
  int ii = 1;
  if (!QUIETMODE) {
    cout << ">> CALCULATING IC LOSS LOOKUP " << endl;
  }

  double phEmax = pow(10.,TargetPhotonVector[TargetPhotonVector.size()-1][0]);
  /* gamma value that indicates Thomson regime (see Blumenthal&Gould) */
  double GammaLow = 1.e-1;
  /* transition energy to Thomson regime. Losses are then simply Edot~E*E*edens*/
  double Etrans = m_e * m_e * GammaLow / phEmax;
  for (double loge = logemin; loge < logemax; loge += logestep) {

    if ((double)ii / bins > 0.0001 * tt && QUIETMODE == false) {
      cout << "\r";
      cout << "    " << (int)(100. * ii / bins) << "\% done" << std::flush;
      tt++;
    }
    ii++;
    double LossRate = 0.;
    double Eelectron = pow(10., loge);
    //double Eelectron = Emax;
    if (Eelectron>Etrans) {
      fPointer IntFunc = &Radiation::ICEmissivityRadFieldIntegrated;
      LossRate =
          Integrate(IntFunc, &Eelectron, EGammaMin, Eelectron,
                    integratorTolerance,integratorKronrodRule);
    } else {
      double gamma = (Eelectron + m_e) / m_e;
      LossRate =
          (4. / 3.) * sigma_T * c_speed * TargetPhotonEdens * gamma * gamma;
    }
    if (std::isnan(LossRate)) {
      cout << __func__ << ",l." << __LINE__ <<": LossRate is nan! Exiting."
           << endl;
      exit(1);
    }
    fUtils->TwoDVectorPushBack(Eelectron,LossRate,ICLossLookup);
  }

  INTEGRATEOVERGAMMAS = false;
  if (QUIETMODE == false) {
    cout << endl;
    cout << "    -> DONE!   " << endl;
    cout << endl;
    cout << ">> CALCULATING OF IC LOSS LOOKUP COMPLETE " << endl;
    cout << endl;
  }
  return;
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
  double nu = *(double *)par / hp;
  /* electron energy */
  double eElectron = x;
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
    double electrons = fUtils->EvalSpline(log10(eElectron),ElectronLookup,
                             acc,__func__,__LINE__);
    value = 4. * pi * sqrt(3.) * el_charge * el_charge * nu_b /
            (hp * hp * nu * c_speed);
    value *= pow(10., electrons) * pow(j, 2.);
    value *= (K_4 * K_1 - (3. / 5.) * j * (K_4 + K_1) * (K_4 - K_1));
  }
  return value;
}

/**
 * emissivity of synchrotron radiation
 * adapted from galprop!ghisellini svensson 1988 'the synchrotron boiler'
 */
double Radiation::SynchEmissivityExplicit(double e, void *par) {

  double eElectron = e;
  double gamma = (eElectron + m_e) / m_e;
  double nu = *(double *)par / hp;

  double norm = sqrt(3.) * el_charge * el_charge * el_charge * BField / m_e;
  double nu_c =
      3. * el_charge * BField * gamma * gamma * c_speed / (4. * pi * m_e);
  double x = nu / nu_c;

  fPointer IntFunc = &Radiation::K_53;
  double *v = NULL;
  double F = x * Integrate(IntFunc, v, x, 1.e2 * x, integratorTolerance,
                           integratorKronrodRule);
  double electrons = fUtils->EvalSpline(log10(eElectron),ElectronLookup,
                                        acc,__func__,__LINE__);
  double val = norm * F * pow(10., electrons) / (hp * hp * nu);

  return val;
}

/* End of the Synchrotron part */

/* ---       BREMSSTRAHLUNG   --- */
/** emissivity of Bremsstrahlung,
 * proton-electron as well as electron-electron
 * From Baring 1999, ApJ, 513, 311-338
 */
double Radiation::BremsEmissivity(double x, void *par) {
  /* initial electron energy */
  double EI = x;
  /* bremsstrahlung photon energy */
  double EP = *(double *)par;
  /* threshold put by hand */
  if (EP < 1.e-12 * EI) return 0.;
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
  return N * pow(10., electrons) / m_e;
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
  double EP = x;
  /* pi0 decay photon energy */
  double Eg = *(double *)par;
  if (EP <= m_p) return 0.;
  double Tp = sqrt(EP * EP - m_p * m_p);
  double Tpth = 0.2797 * GeV_to_erg;
  if (Tp <= Tpth) return 0.;
  if (Eg <= GetMinimumGammaEnergy(Tp)) return 0.;
  if (Eg >= GetMaximumGammaEnergy(Tp)) return 0.;
  double N = DiffPPXSection(Tp, Eg);
  double logprotons = fUtils->EvalSpline(log10(EP),ProtonLookup,
                                      acc,__func__,__LINE__);
  return c_speed * n * N * pow(10., logprotons);
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
 *  the medium.
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
  int size = (int)PARTICLES.size();
  double x[PARTICLES.size()];
  double y[PARTICLES.size()];
  for (unsigned int i = 0; i < PARTICLES.size(); i++) {
    x[i] = log10(PARTICLES[i][0]);
    y[i] = log10(PARTICLES[i][1]);
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

/** Add a greybody distribution of target photons to TotalTargetPhotonGraph,
 * which is used in the
 * IC emission process in this class, but which can also be used 'Particles'
 * class to calculate
 * IC cooling losses
 */
void Radiation::AddThermalTargetPhotons(double T, double edens, int steps) {
  if (edens > 1.e-8)
    cout << "Radiation::AddThermalTargetPhotons: energy density of radiation "
            "field insane. Are you sure of this?" << endl;
  if(edens<=0.) {
    cout << "Radiation::AddThermalTargetPhotons: energy density of target "
            "radiation field negative or zero? Exiting." << endl;
    return;
  }
  double logemin, logemax, low_boundary, high_boundary, low, high, lowtp,
      hightp;
  low_boundary = 1.e-12;
  high_boundary = 1.e6;
  low = log10(low_boundary * kb * T);
  high = log10(high_boundary * kb * T);
  lowtp = log10(targetphotonenergymin);
  hightp = log10(targetphotonenergymax);
  if (!TargetPhotonVector.size()) {
    logemin = low;
    logemax = high;
  } else {
    (lowtp < low) ? logemin = lowtp : logemin = low;
    (hightp > high) ? logemax = hightp : logemax = high;
  }
  double estep = (logemax - logemin) / steps;
  double ePhoton = 0.;
  double nPhoton = 0.;
  vector< vector<double> > vint;
  int i;
  double loge;
  for (loge = logemin, i = 0; loge < logemax; loge += estep, i++) {
    ePhoton = pow(10., loge);
    nPhoton = GreyBody(ePhoton, T, edens);
    if(!nPhoton) continue;
    fUtils->TwoDVectorPushBack(loge,log10(nPhoton),vint);
  }
  AddToTargetPhotonVector(vint);
  return;
}

/** Add an arbitray distribution of target photons to TotalTargetPhotonGraph,
 * which is used in the
 * IC emission process in this class, but which can also be used 'Particles'
 * class to calculate
 * IC cooling losses. This requires as input a 2D vector of format:
 *              ~~~    energy[erg] number_density   ~~~
 * The photons will be added to TotalTargetPhotonGraph
 */
void Radiation::AddArbitraryTargetPhotons(vector<vector<double> > PhotonArray) {
  vector< vector<double> > vint;
  for (unsigned int i = 1; i < PhotonArray.size() - 1; i++) {
    double E = PhotonArray[i][0];
    double N = PhotonArray[i][1];
    if(E <=0. || N <=0.) continue;
    fUtils->TwoDVectorPushBack(log10(E),log10(N),vint);
  }
  AddToTargetPhotonVector(vint);
  return;
}

/** Import target photons from file. File has to be in ASCII format, namely:
 *              ~~~    energy[eV] number_density   ~~~
 * The photons will be added to TotalTargetPhotonGraph
 */
void Radiation::ImportTargetPhotonsFromFile(const char *phFile) {
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
  for (unsigned int i = 0; i < v.size(); i++) {
    fUtils->TwoDVectorPushBack(v[i][0],v[i][1],vint);
  }
  AddToTargetPhotonVector(vint);
  PHfile.close();
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
  double eminxray = 1.e-19;
  double logeminxray = log10(eminxray);
  double logemaxxray =
      log10(1.e-3 * ParticleVector[ParticleVector.size() - 1][0]);
  double estep = (logemaxxray - logeminxray) / steps;
  double E = 0.;
  double N = 0.;
  radiationMechanism = "Synchrotron";
  double U = 2.24;
  vector< vector<double> > vint;
  for (double loge = logeminxray; loge < logemaxxray; loge += estep) {
    E = pow(10., loge);
    N = DifferentialEmissionComponent(E, p) * U / (4. * pi * R * R * c_speed);
    if(N <= 0.) continue;
    fUtils->TwoDVectorPushBack(loge,log10(N),vint);
  }
  AddToTargetPhotonVector(vint);
  return;
}

void Radiation::AddToTargetPhotonVector(vector< vector<double> > vint) {

  gsl_interp_accel *accT1 = gsl_interp_accel_alloc();
  gsl_interp_accel *accT2 = gsl_interp_accel_alloc();
  gsl_spline *Spl = fUtils->GSLsplineFromTwoDVector(vint);
  double logEminSpl = vint[0][0];
  double logEmaxSpl = vint[vint.size()-1][0];
  double stepsSpl = (double)vint.size();
  if (!TargetPhotonVector.size()) {
    gsl_interp_accel_reset(accT1);
    double logdE = (logEmaxSpl - logEminSpl) / stepsSpl;
    for (double logE = logEminSpl; logE < logEmaxSpl; logE += logdE) {
      double val = fUtils->EvalSpline(logE,Spl,accT1,__func__,__LINE__);
      fUtils->TwoDVectorPushBack(logE,val,TargetPhotonVector);
    }
  } else {
    gsl_interp_accel_reset(accT1);
    gsl_interp_accel_reset(accT2);
    /* safe the old vector */
    fUtils->Clear2DVector(TargetPhotonVectorOld);
    for (unsigned int i = 0; i < TargetPhotonVector.size(); i++)
      TargetPhotonVectorOld.push_back(TargetPhotonVector[i]);

    double logEmin, logEmax, logEminOld, logEmaxOld;
    logEminOld = TargetPhotonVector[0][0];
    logEmaxOld = TargetPhotonVector[TargetPhotonVector.size() - 1][0];
    int stepsOld = TargetPhotonVector.size();
    (logEminOld <= logEminSpl) ? (logEmin = logEminOld) : logEmin = logEminSpl;
    (logEmaxOld >= logEmaxSpl) ? (logEmax = logEmaxOld) : logEmax = logEmaxSpl;

    fUtils->Clear2DVector(TargetPhotonVector);
    int steps =
        (int)(stepsOld * (logEmax - logEmin) / (logEmaxOld - logEminOld));
    double logdE = (logEmax - logEmin) / steps;
    for (double logE = logEmin; logE < logEmax; logE += logdE) {
      double val = 0.;
      double valOld = 0.;
      double valSpl = 0.;
      if (logE > logEminOld && logE < logEmaxOld) {
        valOld = fUtils->EvalSpline(logE,TargetPhotonLookup,
                                         accT1,__func__,__LINE__);
        valOld = pow(10., valOld);
      }
      if (logE > logEminSpl && logE < logEmaxSpl) {
        valSpl = fUtils->EvalSpline(logE,Spl,accT2,__func__,__LINE__);
        valSpl = pow(10., valSpl);
      }
      val = valOld + valSpl;
      if (!val) continue;
      val = log10(val);
      fUtils->TwoDVectorPushBack(logE,val,TargetPhotonVector);
    }
  }
  SetTargetPhotonVectorLookup();
  return;
}

/** Function that adds up all individual target photon contributions into
 * TotalTargetPhotonGraph,
 * which is what is then used by the code in the end.
 */
void Radiation::SetTargetPhotonVectorLookup() {
  int size = (int)TargetPhotonVector.size();
  double e[size];
  double n[size];
  double elin[size];
  double en[size];
  double logEOld = 0;
  for (unsigned int i = 0; i < TargetPhotonVector.size(); i++) {
    double logE = TargetPhotonVector[i][0];
    if (logE < logEOld && logEOld) {
      cout << "Radiation::SetTargetPhotonVectorLookup: Target field not "
              "ordered ascending in energy! Exiting!" << endl;
      return;
    }
    double logN = TargetPhotonVector[i][1];
    e[i] = logE;
    n[i] = logN;
    elin[i] = pow(10., e[i]);
    en[i] = pow(10., e[i]) * pow(10., n[i]);
    logE = logEOld;
  }
  targetphotonenergymin = pow(10., e[0]);
  targetphotonenergymax = pow(10., e[size - 1]);
  if(TargetPhotonVectorOld.size()) {
    gsl_spline_free(TargetPhotonLookup);
    gsl_spline_free(TargetPhotonLookupEdens);
  }
  TargetPhotonLookup = gsl_spline_alloc(gsl_interp_linear, size);
  TargetPhotonLookupEdens = gsl_spline_alloc(gsl_interp_linear, size);
  gsl_spline_init(TargetPhotonLookup, e, n, size);
  gsl_spline_init(TargetPhotonLookupEdens, elin, en, size);
  gsl_interp_accel_reset(acc);

  if (gsl_spline_eval_integ_e(TargetPhotonLookupEdens, targetphotonenergymin,
                              targetphotonenergymax, acc, &TargetPhotonEdens))
    TargetPhotonEdens = 0.;
  return;
}

void Radiation::CheckSanityOfTargetPhotonLookup() {

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  for(unsigned int i=0;i<TargetPhotonVector.size();i++) {
    double e = TargetPhotonVector[i][0];
    double n = TargetPhotonVector[i][1];
    double nl = fUtils->EvalSpline(e,TargetPhotonLookup,acc,__func__,__LINE__);
    cout << "rel. diff: " << nl/n - 1. << " " << nl << " " << n <<endl;
  }
  return;
}

/** remove the latest component in
 * TotalTargetPhotonVector and recompute
 * the total target photon spectrum
 */
void Radiation::RemoveLastICTargetPhotonComponent() {
  TargetPhotonVector = TargetPhotonVectorOld;
  SetTargetPhotonVectorLookup();
  return;
}

/** Calculate differential photon spectra for the different radiation processes.
 *  They are stored in the 2D 'diffspec' vector and can be accessed via the
 *  Radiation::ReturnDifferentialSpectrum() member function.
 */
void Radiation::CalculateDifferentialPhotonSpectrum(int steps, double emin,
                                                    double emax) {
  if (emin > emax) {
    cout << "Radiation::ReturnDifferentialSpectrum: Emin>Emax! Check your "
            "boundaries. Exiting..." << endl;
    return;
  }
  if (!steps) {
    cout << "Radiation::ReturnDifferentialSpectrum: Requested 0 steps! "
            "Exiting..." << endl;
    return;
  }
  fUtils->Clear2DVector(diffSpec);
  if (!ElectronVector.size() && !ProtonVector.size()) {
    cout << "Radiation::ReturnDifferentialSpectrum: No particle spectra filled "
            "-> No gamma spectra to calculate. Exiting..." << endl;
    return;
  }
  if (!QUIETMODE) {
    cout << "_________________________________________" << endl;
    cout << ">> CALCULATING SED FROM PARENT PARTICLES " << endl;
  }
  int tt, jj;
  double ICVal, SynchVal, BremsVal, ppVal, E, loge, Emin, Emax, estep;
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
  if (emin) Emin = emin;
  if (emax) Emax = emax;

  estep = (log10(Emax) - log10(Emin)) / steps;

  if (QUIETMODE == false)
    cout << "** Calculating differential gamma-ray emission:" << endl;

  for (loge = log10(Emin), jj = 1, tt = 1; loge < log10(Emax);
       loge += estep, jj++) {
    if ((double)jj / steps > 0.01 * tt && QUIETMODE == false) {
      cout << "\r";
      cout << "    "
           << (int)(100. * (loge - log10(Emin)) / (log10(Emax) - log10(Emin)))
           << "\% done" << std::flush;
      tt++;
    }
    E = pow(10., loge);

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
  fUtils->Clear2DVector(diffSpec);
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
    int i, double emin, double emax) {
  vector<vector<double> > tempVec;
  if (!diffSpec.size()) {
    cout << "Radiation::ReturnDifferentialSpectrum: Differential spectrum "
            "vector empty. Fill it via "
            "Radiation::CalculateDifferentialSpectrum() first! Returning empty "
            "vector." << endl;
    return tempVec;
  }
  double e, dNdE;
  for (unsigned int j = 0; j < diffSpec.size(); j++) {
    e = diffSpec[j][0];
    dNdE = diffSpec[j][i];
    if (e < emin && emin) continue;
    if (e > emax && emax) continue;
    if (dNdE < 0.) dNdE = 0.;
    if (!dNdE) continue;
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
vector<vector<double> > Radiation::ReturnSED(int i, double emin, double emax) {
  vector<vector<double> > tempVec;
  if (!diffSpec.size()) {
    cout << "Radiation::ReturnSED: Differential spectrum vector empty. Fill it "
            "via Radiation::CalculateDifferentialSpectrum() first! Returning "
            "empty vector." << endl;
    return tempVec;
  }
  double e, eTeV, dNdE;
  for (unsigned int j = 0; j < diffSpec.size(); j++) {
    e = diffSpec[j][0];
    eTeV = e / TeV_to_erg;
    dNdE = diffSpec[j][i];
    if (eTeV < emin && emin) continue;
    if (eTeV > emax && emax) continue;
    if (dNdE < 0.) dNdE = 0.;
    if (!dNdE) continue;
    fUtils->TwoDVectorPushBack(eTeV,e * e * dNdE,tempVec);
  }
  return tempVec;
}


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
/**
 * Integration function using the GSL QAG functionality
 *
 */
double Radiation::Integrate(fPointer f, double *x, double emin, double emax,
                            double tolerance, int kronrodrule) {
  double integral, error;
  auto ptr = [=](double xx)->double {
    return (this->*f)(xx, (void *)x);
  };
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);
  GSLfuncRad<decltype(ptr)> Fp(ptr);
  gsl_function F = *static_cast<gsl_function *>(&Fp);
  if (gsl_integration_qag(&F, emin, emax, 0, tolerance, 10000, kronrodrule, w,
                          &integral, &error))
    return 0.;
  gsl_integration_workspace_free(w);
  return integral;
}

vector<vector<double> > Radiation::GetTargetPhotons() {
  vector< vector<double> >  vint;
  for(double i = targetphotonenergymin ; i < targetphotonenergymax ; i*=1.01) {
    double targetphotons = pow(10.,fUtils->EvalSpline(log10(i),TargetPhotonLookup,
                           acc,__func__,__LINE__));
     fUtils->TwoDVectorPushBack(i,targetphotons,vint);
  }
  return vint;
}
