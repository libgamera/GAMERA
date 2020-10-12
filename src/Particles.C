#include "Particles.h"

Particles::Particles() {
  /* Default values */
  SpectralIndex = 2.000;
  METHOD = -1;
  Tmin = 1.;
  TIterationStart = TminInternal;
  Type = 0;
  CutOffFactor = 1000.;
  TminInternal = 1./yr_to_sec; // 1 sec
  TmaxInternal = 1.e100;
  TminLookup = 0.;
  EminInternal = 1.e-3;
  energyMarginFactor = 1.e-3;
  energyMargin = 1.;
  ebins = 100;
  BField = 0.;
  N = 0.;
  R = -1.;
  V = -1.;
  TActual  = 0.;
  gslmemory=5000;
  kronrodrule = 4;
  QUIETMODE = false;
  EminConstantForNormalisationOnly = false; 
  DEBUG = false;
  FASTMODE = false;
  ESCTIME = false;
  IONIZATION = false;
  ICLossLookup = NULL;
  iclosslookupbins = 200;
  LumLookup = NULL;
  NLookup = NULL;
  eMaxLookup = NULL;
  BFieldLookup = NULL;
  RLookup = NULL;
  VLookup = NULL;
  escapeTimeLookupEdep = NULL;
  escapeTimeLookupTdep = NULL;
  EscapeTimeEnergyTimeEvolution = NULL;
  CustomInjectionSpectrum = NULL;
  CustomInjectionSpectrumTimeEvolution = NULL;
  energyTrajectory = NULL;
  energyTrajectoryInverse = NULL;
  integratorTolerance = 5.e-2;
  ParticleSpectrum.clear();
  ParticleSpectrum.resize(0);
  ICLossVector.resize(0);
  LumVector.resize(0);
  NVector.resize(0);
  BVector.resize(0);
  eMaxVector.resize(0);
  escapeTimeVectorEdep.resize(0);
  escapeTimeVectorTdep.resize(0);
  RVector.resize(0);
  VVector.resize(0);
  Constants.resize(4,NAN);
  splines.resize(4);
  vals.resize(4);
  accs.resize(4);
  vs.resize(4);
  CustomInjectionSpectrumTimeEvolutionVector.resize(0);
  EscapeTimeEnergyTimeEvolutionVector.resize(0);
  SpectralIndex2 = eBreak = TminExternal = adLossCoeff = EminConstant = 0.;
  
  escapeTime = 0.;
  eBreakS2 = eBreak2mS2 = eBreakS = eBreak2mS = emin2mS2 = emin2mS =
      emineBreak2mS2 = 0.;
  eBreak2mSInd2 = emin2mSInd2 = emineBreak2mSInd2 = fS = fS2 = bremsl_epf =
      bremsl_eef = 0.;
  LumConstant = BConstant = NConstant = VConstant = RConstant = EmaxConstant = 
      escapeTimeConstant = CustomInjectionNorm = Age = NAN;
  MinELookup = MaxELookup = 0.;
  accIC = gsl_interp_accel_alloc();
  accLum = gsl_interp_accel_alloc();
  accN = gsl_interp_accel_alloc();
  accBField = gsl_interp_accel_alloc();
  acceMax = gsl_interp_accel_alloc();
  accescapeTimeEdep = gsl_interp_accel_alloc();
  accescapeTimeTdep = gsl_interp_accel_alloc();
  accR = gsl_interp_accel_alloc();
  accV = gsl_interp_accel_alloc();
  accTr = gsl_interp_accel_alloc();
  accTrInv = gsl_interp_accel_alloc();
  accCustInj = gsl_interp_accel_alloc();
  taccsp = gsl_interp_accel_alloc();
  eaccsp = gsl_interp_accel_alloc();
  taccesc = gsl_interp_accel_alloc();
  eaccesc = gsl_interp_accel_alloc();


  // grid solver seems to be the best solution in all but the highest energy-loss
  // cases. Setting it as default. User can manually specify semi-analytical method.
  METHOD = 0;

  fUtils = new Utils();
  fRadiation = new Radiation();
}
Particles::~Particles() {
  for (unsigned int i = 0; i < grid.size(); i++) grid[i].clear();
  grid.clear();
}



/** Speed Hack: calculate some constant that are repeatedly needed in nested
 * loops
 */
void Particles::CalculateConstants() {
  eBreakS2 = pow(eBreak, SpectralIndex2);
  eBreak2mS2 = pow(eBreak, 2. - SpectralIndex2);
  eBreak2mSInd2 = log(eBreak);
  eBreakS = pow(eBreak, SpectralIndex);
  eBreak2mS = pow(eBreak, 2. - SpectralIndex);
  emin2mS2 = pow(Emin, 2. - SpectralIndex2);
  emin2mS = pow(Emin, 2. - SpectralIndex);
  emin2mSInd2 = log(Emin);
  emineBreak2mS2 = pow(Emin / eBreak, 2. - SpectralIndex2);
  emineBreak2mSInd2 = log(Emin / eBreak);
  fS = 1. / (2. - SpectralIndex);
  fS2 = 1. / (2. - SpectralIndex2);

  bremsl_epf = 3. * fineStructConst * sigma_T * c_speed * m_e / pi;
  bremsl_eef = (3. * fineStructConst * sigma_T * c_speed * m_e / (2. * pi));
  if(!Type && (fRadiation->GetICLossLookup()).size()) {
	// Take the lookup of the IC losses from the Radiation class
    SetLookup(fRadiation->GetICLossLookup(), "ICLoss");
  }
  return;
}

void Particles::SetSolverMethod(int method) {
  if(method != 0 && method != 1 && method != 2) {
      cout << "Particles::SetSolverMethod: Unsupported solver method ("
           << method << "). Options: 0-grid solver, 1-semianalytic with "
           << "const. losses, 2-semianalytic with no losses. Ignoring this "
           << "foolishness." << endl;
      exit(1);
  }
  METHOD = method;
  return;
}

/** fill the lookup holding the particle spectrum {E(erg) - N(erg^-1) */
void Particles::CalculateParticleSpectrum(string type, int bins, bool onlyprepare,
                                          bool dontinitialise) {

  if (!type.compare("electrons")) {
    Type = 0;
  } else if (!type.compare("protons")) {
    Type = 1;
  } else {
    cout << "Particles::CalculateParticleSpectrum: Particle type not "
    		"recognized. Please specify "
            "supported particle species! " << endl;
  }
  if (!QUIETMODE) {
    cout << "___________________________________" << endl;
    cout << ">> STARTING NEW PARTICLE EVOLUTION " << endl;
    if (Type == 0) cout << "   (-> electrons)     " << endl;
    if (Type == 1) cout << "   (-> protons)     " << endl;
  }

  if(std::isnan(Age)) {
    cout << "Particles::CalculateParticleSpectrum: No age set! Exiting"
         << endl; 
    return;
  }
    
  if (VVector.size() && !RVector.size()) {
    cout << "Particles::CalculateParticleSpectrum: You set a velocity time "
         << "lookup but no radius time lookup. Both need to be set via "
         << "Particles::SetRadiusLookup() and Particles::SetVelocityLookup(), "  
         << "otherwise I exit! Exiting. "<<endl;
    exit(1);
  }

  if(!LumVector.size() && std::isnan(LumConstant) 
     && CustomInjectionSpectrum == NULL 
     && CustomInjectionSpectrumTimeEvolution == NULL) {
    cout << "Particles::CalculateParticleSpectrum: Particle vector empty! Exiting"
         << endl;
    return;
  }

  /* reset Particle Lookup */
  fUtils->Clear2DVector(ParticleSpectrum);

  ebins = bins;

  if (Type == 1)
    Emin = m_p;  // norm for protons or for electrons if defined per e/p
                 // fraction
  else
    Emin = m_e;
  /* if lower energy bound is externally set
   * use this value for the normalisation of the spectrum
   */
  if (EminConstant) {
    Emin = EminConstant;
  }
  if (EminConstant>EminInternal) {
    EminInternal = 1.1*EminConstant;
  }
  CalculateConstants();
  DetermineLookupTimeBoundaries();

  /* determine time from where to start the iteration. Particles that would have
   * been injected before that time are injected as a blob at Tmin. This can
   * lead to bumps at low energies, depending on cooling.
   */
  Tmin = 1; // standard start of iteration time , 1 year.
  if (TminExternal) Tmin = TminExternal; // replace with externally set TminExternal
  if (Tmin<TminInternal) Tmin=TminInternal; // if external time is smaller than internal minimum time (standard 1s), set to latter

  /* get the upper energy boundary of the spectrum. This can either be
   * externally
   * set (if) or dynamically determined by the code (else).
   */
  if (!std::isnan(EmaxConstant))
    eMax = EmaxConstant;
  else
    eMax = DetermineEmax(Tmin);

  /* apply a safe margin to the upper energy boundary
   * in order to prevent numerical effects at the upper end of the spectrum.
   */
  energyMargin = pow(-log(energyMarginFactor), 1. / CutOffFactor);
  eMax *= energyMargin;

  /* if emax falls below emin, return dummy vector with zeroes */
  if (eMax <= Emin) {
    cout << "Particles::CalculateParticleSpectrum Whaat? eMax lower than "
            "Emin: eMax = " << eMax << " Emin = " << Emin << endl;
    return;
  }

  if (!METHOD) {
    PrepareAndRunNumericalSolver(ParticleSpectrum, onlyprepare, dontinitialise);
  }
  else if (METHOD==1) {
	// semi-analytical
    CalculateEnergyTrajectory();
    CalcSpecSemiAnalyticConstELoss();
  }
  else {
    CalcSpecSemiAnalyticNoELoss();
  }   
  if (!QUIETMODE) {
    cout << ">> PARTICLE EVOLUTION DONE. EXITING." << endl;
    cout << endl;
  }
  
  return;
}

/**
 * Particle injection spectrum, provided by a 2D vector (1D-lookup) with tuples
 * energy - differential luminosity.
 * Can also handle time-changing spectral shapes of injection spectra by providing
 * a 2D lookup in time and energy.
 * A historic option is to provide spectral index (indices) and high energy cutoff
 * value, this will result in a power-law (broken power-law) injection.
 */
double Particles::SourceSpectrum(double e) {
  double val = 0.;
  if(CustomInjectionSpectrumTimeEvolution != NULL) {
    if(e < MinELookup || e > MaxELookup) return 0.;
    val = pow(10.,interp2d_spline_eval(
              CustomInjectionSpectrumTimeEvolution,TActual,log10(e),taccsp,eaccsp));
  }
  else if(CustomInjectionSpectrum != NULL) {
    if(e < MinELookup || e > MaxELookup || !Lum ) return 0.;
    val = Lum * 
      pow(10.,gsl_spline_eval(CustomInjectionSpectrum, log10(e), accCustInj)) 
      / CustomInjectionNorm;
  }
  // this is perhaps outdated, as custom injection functions also cover this 
  // special case. FIXME: Remove this?
  else val = PowerLawInjectionSpectrum(e, eMax, 10. * eMax);

  if (val <1.e-99)  return 0.;
  else return val;
}

/**
 * Escape time at energy \a e (erg) and time \a t (yrs).
 *
 * @param e : energy of the particle [erg]
 * @param t : time [yr]
 * @return escape time of the particle in yr
 */
double Particles::EscapeTime(double e, double t) {
  // energy and time dependent escape -> 2D intepolation in t,e lookup
  if(EscapeTimeEnergyTimeEvolution != NULL) {
    if(e < MinELookup || e > MaxELookup) return 0.;
    return pow(10.,interp2d_spline_eval(
              EscapeTimeEnergyTimeEvolution, t, log10(e),taccesc,eaccesc));
  }
  // energy dependent, but time constant escape -> 1D interpolation in e lookup
  else if(escapeTimeLookupEdep != NULL) {
    if(e < MinELookup || e > MaxELookup) return 0.;
    return pow(10.,
               gsl_spline_eval(escapeTimeLookupEdep, log10(e), accescapeTimeEdep));
  }
  // time dependent, but energy constant escape -> 1D interpolation in t lookup
  else if(escapeTimeLookupTdep != NULL) {
    return pow(10.,
               gsl_spline_eval(escapeTimeLookupTdep, t, accescapeTimeTdep));
  }
  // overall constant for escape time
  else if(!std::isnan(escapeTimeConstant)) return escapeTimeConstant;
  else return 0.;
}

/**
 * Power law injection spectrum. If eBreak and SpectralIndex2 are specified,
 * returns a broken power law.
 *
 * This is perhaps outdated, as custom injection functions also cover this
 * special case. FIXME: Remove this?
 * Advised to use the Particles::CustomInjectionSpectrum function*
 *
 */
double Particles::PowerLawInjectionSpectrum(double e, double ecut,
                                            double emax) {

  /* exponential cut factor */
  double cutOffFactorInternal = exp(-pow(e / ecut, CutOffFactor));
  /* calculate norm */
  double integral = 0.;
  /* Broken power-law. Used, if both SpectralIndex2 and a break energy are
   * specified */
  // FIXME: special case SpectralIndex,SpectralIndex2 = 2
  if (SpectralIndex2 && eBreak) {
    if (SpectralIndex != 2. && SpectralIndex2 != 2.) {
      integral = eBreakS2 * (eBreak2mS2 - emin2mS2) * fS2;
      integral += eBreakS * (pow(emax, 2. - SpectralIndex) - eBreak2mS) * fS;
      if (emax <= eBreak) {
        integral =
            (pow(emax / eBreak, 2. - SpectralIndex2) - emineBreak2mS2) * fS2;
      }
    } else {
      cout << "Particles::PowerLawInjectionSpectrum: "
    		  "SpectralIndex,SpectralIndex2 = 2 case not implemented!"
    		  << endl;
    }
  }
  /* else, a single power-law is used in the integrated flux (introduces a small
   * error if
   * depending on "CutOffFactor":
   */
  else {
    if (SpectralIndex != 2.)
      integral = (pow(emax, 2. - SpectralIndex) - emin2mS) *
                 fS;  // FIXME replace with errfunc for exp cutof. (minor error)
    else
      integral = (log(emax) - emin2mSInd2);
  }
  /* calculate the normalisation */
  double SourceSpectrumNorm = Lum / integral;
  double J = 0.;
  /* Broken power-law. If both SpectralIndex2 and a break energy was specified */
  if (SpectralIndex2 && eBreak) {
    if (e > eBreak)
      J = SourceSpectrumNorm * pow(e / eBreak, -SpectralIndex) *
          cutOffFactorInternal;
    else if (e <= eBreak && emax <= eBreak)
      J = SourceSpectrumNorm * pow(e / eBreak, -SpectralIndex2) *
          cutOffFactorInternal;
    else
      J = SourceSpectrumNorm * pow(e / eBreak, -SpectralIndex2);
  }
  /* super-exponential cutoff as ~exp[-(e/emax)^cutOffFactorInternal] */
  else {
    J = SourceSpectrumNorm * pow(e, -SpectralIndex) * cutOffFactorInternal;
  }
  if (DEBUG == true) {
    J = 0.;
  }
  return J;
}

/** 
 * Set the values for the class variables and lookups like:
 * Bfield, Luminosity, Maximum energy, density.
 *
 * @param t : time to be given in units of years
 */
void Particles::SetMembers(double t) {
  if (t < TminInternal || t > TmaxInternal) {
    cout << "Particles::SetMembers: Time (" << t << "yrs vs {" << TminInternal
         << "," << TmaxInternal << "}) outside boundaries." << endl;
    cout << "  -> Using values from previous time step. If you don't want "
            "this, you have to extend " << endl;
    cout << "     your Lookups (Luminosity,R,V,B,N etc) or set constant values "
            "(via e.g. SetBField()). " << endl;
    return;
  }
  if (t == TminInternal) t *= 1. + 1.e-10;
  if (t == TmaxInternal) t /= 1. + 1.e-10;
  TActual = t;

  Constants[0] = LumConstant;  Constants[1] = NConstant;
  Constants[2] = BConstant;  Constants[3] = EmaxConstant;

  splines[0] = LumLookup;  splines[1] = NLookup;
  splines[2] = BFieldLookup;  splines[3] = eMaxLookup;

  accs[0] = accLum;  accs[1] = accN;
  accs[2] = accBField;  accs[3] = acceMax;

  vals[0] = &Lum;  vals[1] = &N;
  vals[2] = &BField;  vals[3] = &eMax;

  vs[0] =  LumVector;  vs[1] =  NVector; //FIXME: are they necessary here?
  vs[2] =  BVector;  vs[3] =  eMaxVector; //FIXME: are they necessary here?
  
  // check whether constants or lookups are present and set value accordingly.
  for (unsigned int i = 0; i < Constants.size(); i++) {
    *vals[i] = 0.;
    if (!std::isnan(Constants[i])) *vals[i] = Constants[i];
    else if (splines[i] == NULL) continue;
    else {
       double value = gsl_spline_eval(splines[i], t, accs[i]);
      if(i == 0 || i == 3) value = pow(10.,value);
      *vals[i] = value;
    }
  }

  // the following is important for adiabatic losses
  R = V = adLossCoeff = 0.;
  if (!std::isnan(VConstant)) V = VConstant;
  if (!std::isnan(RConstant)) R = RConstant;
  R += yr_to_sec*(t-TIterationStart)*V;
  if (R < 0.) {
    cout << "Particles::SetMembers: Radius (" << R << " cm) "
         << " is negative! Setting it to 0."  << endl;
    R = 0.;
  }

  if(RVector.size() && t > RVector[0][0] &&
          t < RVector[RVector.size() - 1][0])
    R = gsl_spline_eval(RLookup, t, accR) * pc_to_cm;

  if(VVector.size() && t > VVector[0][0] &&
          t < VVector[VVector.size() - 1][0])
    V = gsl_spline_eval(VLookup, t, accV);
  if (R && V) adLossCoeff = V / R;
  return;
}

/**
 * Sets up lookup for the changed of various quantities:
 * "ICLoss", "Luminosity", "AmbientDensity", "BField",
   "Emax", "Radius","Speed","EscapeTimeTdep"
 */
void Particles::SetLookup(vector<vector<double> > v, string LookupType) {
  int size = (int)v.size();
  if (!size) {
    cout << "Particles::SetLookup: lookup vector empty. Exiting." << endl;
    return;
  }
  vector< vector< double > > lookup;
  if (!LookupType.compare("Luminosity") || !LookupType.compare("Emax")
      || !LookupType.compare("EscapeTimeTdep")) {
    lookup = fUtils->VectorAxisLogarithm(v,1);  //transforms in log10
  }
  else {
    lookup = v;
  }
  gsl_spline *ImportLookup = fUtils->GSLsplineFromTwoDVector(lookup);
  std::string stArr[] = {"ICLoss", "Luminosity", "AmbientDensity", "BField",
                    "Emax", "Radius","Speed","EscapeTimeTdep"};
  vector<string> st( stArr, stArr + ( sizeof ( stArr ) /  sizeof ( stArr[0] ) ) );
  gsl_spline **splArr[] = {&ICLossLookup, &LumLookup,  &NLookup,
                            &BFieldLookup, &eMaxLookup,
                            &RLookup, &VLookup, &escapeTimeLookupTdep};
  vector<gsl_spline **> spl( splArr, splArr + ( sizeof ( splArr ) /  sizeof ( splArr[0] ) ) );

  vector<vector<double> > * vsArr[] = {
      &ICLossVector, &LumVector, &NVector, &BVector,
      &eMaxVector, &RVector, &VVector, &escapeTimeVectorTdep};
  vector<vector<vector<double> > *> vs( vsArr, vsArr + ( sizeof ( vsArr ) /  sizeof ( vsArr[0] ) ) );
  for (unsigned int i = 0; i < st.size(); i++) {
    if (!LookupType.compare(st[i])) {
      if (*spl[i] != NULL) {
        gsl_spline_free(*spl[i]);
      }
      *spl[i] = ImportLookup;
      *vs[i] = lookup;   
      DetermineLookupTimeBoundaries();
      return;
    }
  }
  return;
}

/** Append lookup to existing lookup. Time order has to be right, so the
 * youngest
 * appended entry has to be older than the oldest already existing one.
 * This function is useful when iteratively determining B-Field, CR Lum, source
 * speed and radius etc... That is, when the particle distribution influences
 * shock dynamics and B-Field.
 */
void Particles::ExtendLookup(vector<vector<double> > v, string LookupType) {
  if (!v.size()) {
    cout << "Particles::ExtendLookup: Input vector empty. Exiting." << endl;
    return;
  }
  string stArr[] = {"Luminosity", "AmbientDensity", "BField", "Emax",
                       "Radius",         "Speed"};
  vector<string> st( stArr, stArr + ( sizeof ( stArr ) /  sizeof ( stArr[0] ) ) );
  vector<vector<double> > vsArr[] = {LumVector,  NVector,          BVector,
                                         eMaxVector, RVector,  VVector};
  vector<vector<vector<double> > > vs( vsArr, vsArr + ( sizeof ( vsArr ) /  sizeof ( vsArr[0] ) ) );
  vector< vector< double > > lookup;
  if (!LookupType.compare("Luminosity") || !LookupType.compare("Emax")) {
    lookup = fUtils->VectorAxisLogarithm(v,1);
  }
  else {
    lookup = v;
  }
  for (unsigned int i = 0; i < st.size(); i++) {
    if (!LookupType.compare(st[i])) {
      if (vs[i][vs[i].size() - 1][0] > lookup[0][0]) {
        cout << "Particles::ExtendCRLumLookup - WTF, the vector which to add ("
             << st[i] << ") to the existing one starts at earlier times than "
                         "the existing ones ends. Please keep time order in "
                         "the vector! exiting." << endl;
        return;
      }
      if(vs[i][vs[i].size() - 1][0] == lookup[0][0])
        vs[i].insert(vs[i].end(), lookup.begin()+1, lookup.end());
      else vs[i].insert(vs[i].end(), lookup.begin(), lookup.end());
      SetLookup(vs[i], LookupType);      
      DetermineLookupTimeBoundaries();
      return;
    }
  }
  return;
}

/**
 * Calculates the energy loss rate at a given time t from the shock dynamics,
 * ambient photon and B-fields.
 *
 * Arguments:
 *   double E : energy of the particle [erg]
 * Returns:
 *   double loss_rate (in seconds)
 */
double Particles::EnergyLossRate(double E) {
  synchl = icl = adl = bremsl = ppcol = ionization = 0.;
  double bremsl_ep = 0.;
  double bremsl_ee = 0.;
  double gamma;
  if (Type) {
    gamma = (E + m_p) / m_p;
  }
  else{
    gamma = (E + m_e) / m_e;
  }
  if (gamma < 1.) gamma = 1.;
  double gamma2 = gamma * gamma;
  double beta=sqrt(1-1/gamma2);
  double p = sqrt(gamma2 - 1.);
  /* S-parameter. See Haug2004 (eq.2)
   * Here assumed as default hydrogen and 10% Helium
   * TODO: use the new formulation in the Radiation class
   * to allow different compositions.
   */
  double S = N * 1.4;
  /* synchrotron losses */
  synchl = (4. / 3.) * sigma_T * c_speed * BField * BField * gamma * gamma /
           (8. * pi);
  
  if(!ICLossVector.size()) icl=0.;
  else {
    /* IC losses from the lookup table (ICLossLookup) */
    icl = gsl_spline_eval(ICLossLookup, E, accIC);
  }

  /* adiabatic losses (adlossCoeff = V/R) */
  adl = adLossCoeff * E;

  /** Bremsstrahlung losses Haug+2004 */
  /* electron-proton bremsstrahlung */
  /* TODO: to be super self-consistent, calculate losses in a lookup, analog
   * to the IC lookup.
   */
  bremsl_ep = ((2. * gamma2 / 9. - 19. * gamma * p * p / 675. -
                0.06 * p * p * p * p / gamma) *
                   p * p * p / (gamma * gamma * gamma * gamma * gamma * gamma) +
               gamma * log(gamma + p) - p / 3.) *
              bremsl_epf * S * gamma2 / (gamma2 + p * p);

  /*
   * Electron-electron bremsstrahlung
   * Assuming also here the 10% Helium
   */
  bremsl_ee =
      N * (1.2) * bremsl_eef * (p * (gamma - 1.) / gamma) * (log(2. * gamma) - 1. / 3.);
  bremsl = bremsl_ep + bremsl_ee;

  /* in case of protons, take adiabatic and ionization losses into account */
  // TODO: for protons, define analog to loss rate from (catastrophic)
  // proton-proton interaction time scale?
  if (Type) {
    synchl = 0.;
    icl = 0.;
    bremsl = 0.;
    //Collision pp
    double Tpth=2*m_pi+(m_pi*m_pi/(2*m_p));
    if(E>Tpth) {
      double e=E/Tpth;
      double k=0.3;
      double Cbarn_cm=1e-24;
      double sigma;
      double sigma1=30.7-0.96*log(e)+0.18*log(e)*log(e);
      double sigma2=1-pow(1/e,1.9);
      sigma=sigma1*pow(sigma2,3)*1e-3*Cbarn_cm;
      ppcol=N*E*k*c_speed*sigma;
    }
    else {
      ppcol=0;
      }
    /* ionization losses */
    // Mean value in the interstellar medium if you consider 10% He
    if (IONIZATION) {
      double z=1.2;  // Includes the effect of 10% Helium in the target.
      double Z=1;    // The ions are the protons
      double I=15*eV_to_erg; //erg
      // Ionization should not play a role if the plasma is fully ionized!
      // If you want to apply these losses, then need to activate the IONIZATION flag
      // Ionization computed here is valid only for protons
      ionization=N*z*4.*pi*pow(el_charge,4)*(Z*Z)/(m_e*beta/c_speed)*(log((2*m_e*pow(beta,2)/I)*gamma2)-pow(beta,2));
    }
  }
  if (DEBUG == true) {
    synchl = 0.;
    bremsl = 0.;
    adl = 0.;
  }
  return synchl + icl + adl + bremsl + ionization + ppcol;
}

/**
 * Prepare axes for the numerical solver (and if onlyprepare==false) call the
 * solver
 *
 * @param particlespectrum : spectrum of the particles
 * @param onlyprepare : if false, calls the solver and initialise
 *                      the grid (default false)
 * @param dontinitialise : if false, sets the initial conditions (default false)
 */
void Particles::PrepareAndRunNumericalSolver(
    vector<vector<double> > &particlespectrum, bool onlyprepare,
    bool dontinitialise) {
  GetAxis(Emin, eMax, ebins, energyAxis, true);
  CreateGrid();
  if (dontinitialise == false) SetInitialCondition(energyAxis, Tmin);

  /* if onlyprepare=true, only the axes are initialised, but the grid is not
   * computed
   * until t=age. This is useful in conjunction with
   * ComputeGridInTimeInterval().
   */
  if (onlyprepare == false) ComputeGrid(energyAxis, Tmin, Age, Age*yr_to_sec/1e3);
  return;
}

/**
 * create an axis that attributes each bin with a real value. (e.g. time axis, energy axis).
 */
void Particles::GetAxis(double min, double max, int steps, vector<double> &Axis,
                        bool logarithmic) {

  double min_mass = Type ? m_p : m_e;
  double quot = 1.;
  if (logarithmic == true) {
    min = log10(min); max = log10(max);
    quot = (max - log10(min_mass)) / (max - min);   
    min = log10(min_mass);
  }
  else {
    quot = (max - min_mass) / (max - min);    
    min = min_mass;
  } 

  steps = int(quot * steps);
  double binsize = (max - min) / steps; //FIXME: replace by ebins!
  ebins = steps +1;
  Axis.resize(steps+1);
  for (unsigned int i = 0; i < Axis.size(); i++) {
    Axis[i] = min + i * binsize;
  }

  return;
}

/**
 * create the propagation grid out of 2 1D-vectors
 * energy in x-direction, time in y-direction
 */
void Particles::CreateGrid() {
  fUtils->Clear2DVector(grid);
  grid.push_back(vector<double>());
  for (int i = 0; i < ebins; i++) {
    grid[grid.size() - 1].push_back(0.);
  }
  return;
}

/**
 * Determine the minimum time from where to start
 * the calculation. Electrons before that time are injected as
 * a single 'blob'. This time is derived from the requirement
 * that the blob has slid down to energies e.g. E<1GeV(EMIN) at t=Age.
 *
 * @param emin : lost energy at time t=age
 * @param tmin : pointer to the time
 *
 */
void Particles::DetermineTMin(double emin, double &tmin) {
  double logt, logtmin, logtmax, logdt, logsteps, TMIN;
  if (TminInternal > 0.)
    logtmin = log10(TminInternal);
  else
    logtmin = log10(Age) - 5.;
  TMIN = 0.;
  logtmax = log10(Age);
  if (logtmin > logtmax) logtmin = logtmax - 3.;
  logsteps = 30.;
  logdt = (logtmax - logtmin) / logsteps;
  for (logt = logtmin; logt < logtmax; logt += logdt) {
    CalculateEnergyTrajectory(pow(10., logt));
    gsl_interp_accel_reset(accTrInv);
    if (energyTrajectoryInverse == NULL) continue;
    if (gsl_spline_eval_e(energyTrajectoryInverse, log10(emin), accTrInv,
                          &TMIN))
      continue;
    TMIN = pow(10., TMIN);
    if (TMIN > Age) break;
    tmin = pow(10., logt);
  }
  return;
}

/**
 * Determine the maximum energy of particles between tmin and Age.
 * This energy is then used as upper boundary for the energy dimension of the
 * grid.
 */
double Particles::DetermineEmax(double tmin) {
  if(eMaxLookup==NULL) {
    cout << "Particles::DetermineEmax: Neither constant value for maximum "
            "particle energy set nor time evolution lookup... Exiting." << endl;
    exit(1);
  }
  
  double t = 0.;
  double dt = 0.;
  double tt = 0.;
  double eMaxHistory = -1.;
  t = tmin;
  dt = (Age - tmin) / 10000.;
  while (t <= Age) {
    if (gsl_spline_eval_e(eMaxLookup, t, acceMax, &tt)) continue;
    if (tt > eMaxHistory) eMaxHistory = tt;
    t += dt;
  }
  return pow(10.,eMaxHistory);
}

/** set initial condition (a.k.a. set the first energy vector at t=tmin of the
 * grid) */
void Particles::SetInitialCondition(vector<double> EnergyAxis,
                                    double startTime) {

  double t0 = startTime;
  SetMembers(t0);
  for (unsigned int i = 0; i < EnergyAxis.size()-1; i++) {//last bin not filled, it is just as an edge. it wont hold any information later and will be removed. its function is to prevent the spectrum edge from oscillating. 
    double e =pow(10., EnergyAxis[i]);
    if (e>Emin) grid[0][i] = t0 * yr_to_sec * SourceSpectrum(e);
    else grid[0][i] = 0.;
  }

  return;
}

/**
* Center piece of the class: Numerical solver of the particle spectrum.
* It treats cooling as an advective flow in energy space, and uses a piece-wise
* linear numerical scheme to transport particles from one energy bin in the
* next,
* where the (energy-dependent) cooling rate is treated as the flow velocity of
* the 'fluid'. To realise sharp edges in the resulting spectra, the 'SuperBee'
* and 'MinMod' slope
* limiters are available (default: MinMod). This slows down the procedure, and
* can manually be
* disabled.
* TODO: write option to choose between slope limiter methods / disable them
*/
void Particles::ComputeGrid(vector<double> EnergyAxis, double startTime,
                            double Age, double minTimeBin) {

  double value = 0.;
  double quot = 0.;
  double t = 0.;
  double e0 = 0.;
  double e1 = 0.;
  double e2 = 0.;
  double ebin = 0.;
  double tbin = 0.;
  double deltaE1 = 0.;
  double deltaE2 = 0.;
  double ElossRate_e1 = 0.;
  double ElossRate_e2 = 0.;
//  unsigned int Esize = EnergyAxis.size() - 1;
  int tt = 1;
  int largestFilledBin = EnergyAxis.size() - 1;;
//  long int count = 0;
  if(minTimeBin > Age * yr_to_sec) minTimeBin = TminInternal * 10. *yr_to_sec;//Age*yr_to_sec/1000.;
  /* append a new energy vector that will always hold the energy spectrum at the
   *  next time step and initialise it with zeroes.
   */


  vector<double> E;
  vector<double> Ecuts;
  vector<double> EscTime;
  TIterationStart = startTime;
  vector< vector<double> > lossrates;
  double eMax0 = 0.;
  grid.push_back(vector<double>());
  int i0 = 0;
  int i1 = 1;

  for (unsigned int i = 0; i < EnergyAxis.size(); i++) {
    grid[grid.size() - 1].push_back(0.);
    double e =  pow(10., EnergyAxis[i]);
    E.push_back(e);
  }

  /* info writeout. Disable it by using 'ToggleQuietMode()' */
  if (!Type && QUIETMODE == false)
    cout << "** Evolving Electron Spectrum:" << endl;
  else if (Type == 1 && QUIETMODE == false)
    cout << "** Evolving Proton Spectrum:" << endl;

  /* main loop over time  */
  for (double T = startTime; T <= Age; T += tbin / yr_to_sec) {

    /* Set Members (CR luminosity, B-field etc.) at time T */
    SetMembers(T);
    Ecuts.clear();
    EscTime.clear();
    double MinEscTime = 1.e100;
    if(eMax != eMax0) {
        for (unsigned int i = 0; i < E.size(); i++) {
          Ecuts.push_back(exp(-pow( E[i] / eMax, CutOffFactor)));
        }      
    }
    if (ESCTIME == true) {
        for (unsigned int i = 0; i < E.size(); i++) {
            double esctime = EscapeTime(E[i],T);
            if(esctime && esctime < MinEscTime) MinEscTime = esctime;
            EscTime.push_back(esctime);
        }
    }
    /* dynamically determine tbin size. This is a critical step
     * for the speed of the algorithm. Since the time step size is
     * proportional to Ebinsize(eMax)/Edot(eMax) and Edot ~ E^2,
     * eMax - or the largest relevant energy bin - should be chosen
     * as small as possible. Thus, don't use the energy at eMax, but
     * but rather of the 
     * If an external emax is specified, always choose the highest energy bin
     * value.
     */

    tbin = 1e50;
    for(unsigned int i=1; i < E.size(); i++) {
        double t_temp = fabs((E[i]-E[i-1]) / EnergyLossRate(E[i]));
        if (t_temp < tbin) tbin = t_temp;
    }
    /* compare minimum timescale to that of particle escape 
     * and choose the smaller one
     */
    if (MinEscTime && MinEscTime < tbin) tbin = MinEscTime;

    /* if losses become small (e.g. degrading B-field, or low eMax), time bins
     * may become very large. This can become problematic if tbin << Age no
     * longer
     * holds. Then, replace tbin by minTimeBin (default: 1yr) if
     * tbin>minTimeBin.
     */
    if (tbin > minTimeBin) tbin = minTimeBin;


    /* info writeout */
    if (T > 0.01 * tt * (Age - startTime) && QUIETMODE == false) {
      cout << "\r";
      if (tbin / (yr_to_sec * Age) > 1.e-7) {
        cout << "                                                              "
                "         \r" << std::flush;
        cout << "    " << (int)(100. * (T - startTime) / (Age - startTime))
             << "\% done \r" << std::flush;
      } else {
        cout << "    " << (int)(100. * (T - startTime) / (Age - startTime))
             << "\% done (Energy losses are very high, iteration might take a "
                "while)" << std::flush;
      }
      tt++;
    }
    /* if eMax drops below the lower energy bound of the grid, exit. */
    if (largestFilledBin <= 0) break;

    /* This is the new time! (halve a time step forward) */
    t = T + 0.5 * tbin / yr_to_sec;
    /* if t is larger than age, exit! */
    if(t>Age) break;
    
    /* update the Members at the new time. */
    Lum = 0.;
    SetMembers(t);

    /* iterate over the previous spectrum, stored in grid[0] and calculate
     * grid[1] from it. This is done using a standard, piece-wise linear
     * advection
     * scheme. Per default, also a slope limiter is implemented (MinMod method),
     * that preserves sharp edges in the particle spectrum rather than smearing
     * it out as in a pure donor-cell algorithm.
     */

    /* just for the first step (for speed) */
    e0 = 0.;
    ElossRate_e2 = EnergyLossRate(pow(10., EnergyAxis[0]));
    double ElossRate_e0 = 0.;
    double particleCount = 0.;
    for (unsigned int i = 0; i < E.size()-1; i++) {

      value = 0.;
      quot = 0.;
      ebin = 0.;
      if(i) e0 = E[i-1];
      e1 = E[i]; 
      e2 = E[i+1]; 

      /* The following block calculates the streaming of particles in an out
       * of energy bin 'i' by cooling. This component of increase / decrease of
       * particles in bin 'i' is caused by particles already present in the
       * last time step (grid[0])
       */
      ebin = e2 - e1;
      quot = tbin / ebin;

      ElossRate_e1 = ElossRate_e2;
      ElossRate_e2 = EnergyLossRate(e2);
      ElossRate_e0 = i ? EnergyLossRate(e0): ElossRate_e1;
        

      deltaE1 = tbin * ElossRate_e1;
      deltaE2 = tbin * ElossRate_e2;


      /* Increase in particles in bin 'i' due to particle injection from the
       * source */
      value = tbin * SourceSpectrum(e1);
      value *= Ecuts[i];
      value += grid[i0][i];

      /* Donor-cell advection */
      if(ElossRate_e1>=0. && i) value -= quot * grid[i0][i] * ElossRate_e1;
      if(ElossRate_e2>=0.) value += quot * grid[i0][i+1] * ElossRate_e2;

      /* slope limiters (per default enabled, toggle on/off with FASTMODE flag) */    
      if(FASTMODE == false) {
        double mm,sb;
        /* Uncomment the following lines if you want to use the
         * Lax-Wendroff slope limiter. You also have to comment
         * the code after until the curly bracket. Normally you
         * don't want to do this, the MinMod slope limiter is better */
        /*mm =  quot * (LaxWendroffSlope(i, ebin, i0) * ElossRate_e1 *
                (ebin - deltaE1) -
                             LaxWendroffSlope(i + 1, ebin, i0) * ElossRate_e2 *
                                 (ebin - deltaE2));
        double lw = 0.5*quot*(LaxWendroffSlope(i,ebin,i0)*ElossRate_e1*(ebin-deltaE1)-LaxWendroffSlope(i+1,ebin,i0)*ElossRate_e2*(ebin-deltaE2));
        value -= lw;*/
        
        /*
         * Comment the next lines until the curly bracket if the
         * Lax-Wendroff slope limiter is used.
         */
        mm =  0.5*quot * (GetMinModSlope(i, ebin, i0) * ElossRate_e1 *
                                 (ebin - deltaE1) -
                             GetMinModSlope(i + 1, ebin, i0) * ElossRate_e2 *
                                 (ebin - deltaE2));
        sb =
            0.5*quot*(GetSuperBeeSlope(i,ebin,i0)*ElossRate_e1*(ebin-deltaE1)-GetSuperBeeSlope(i+1,ebin,i0)*ElossRate_e2*(ebin-deltaE2));
        value -= 0.5 * sqrt(mm*mm + sb*sb);
      }

      //TODO: get slope limiter right in case of adiabatic heating!
      if(ElossRate_e1<0.) {
                             value += quot * grid[i0][i] * ElossRate_e1;
      }
      if(ElossRate_e0 < 0. && i) {
                             value -= quot * grid[i0][i-1] * ElossRate_e0;
      }


      /* Decrease in particles in bin 'i' due to particle escape.*/
      if (ESCTIME == true) escapeTime = EscTime[i];
      if (escapeTime > 0.) value -= tbin * grid[i0][i] / escapeTime;

      /* if value is smaller 0 (can happen if escape time scale is low), set
       * to 0. */
      if (value < 0.) value = 0.;


      /* Determine the largest filled energy bin (needed for the efficient
       * calculation of the next iterative time bin.
       */
      if (value > 0.) largestFilledBin = i;

      /* set the particle number (value) in at bin 'i' in grid[i] */
      grid[i1][i] = value;
      particleCount += value * ebin;
    }
    /* set the just now calculated spectrum as base spectrum for the next step
     * This way, only 2 vectors are needed for the calculation of the spectrum
     */
    if(particleCount) {i0 = !i0; i1 = !i1;}

  }

  /* Fill the final lookup, holding the time evolved spectrum at time = Age.
   * Also, forego edge bins in order to avoid artifacts.
   */
  TIterationStart = 0.;
  fUtils->Clear2DVector(ParticleSpectrum);
  for (unsigned int j = 0; j < EnergyAxis.size() -1; j++) {
    e1 = pow(10., EnergyAxis[j]);
    double val = grid[i0][j];
    if(std::isnan(val) || std::isinf(val) || !val || val < 1.e-100)
      continue;
    if(val>0.) fUtils->TwoDVectorPushBack(e1,val,ParticleSpectrum);
  }
  /* Important for wrapper function 'ComputeGridInTimeInterval': remove the last
   * vector (grid[1]) so that grid has the same shape as in the beginning of
   * this function.
   */
  grid.pop_back();
  /* for the format of the info writeout */
  if (QUIETMODE == false) {
    cout << endl;
    cout << "    -> DONE!   " << endl;
    cout << endl;
  }
  return;
}

/** slope for MinMod slope limiter method */
double Particles::GetMinModSlope(int i, double deltaX, int i0) {
  double a = (grid[i0][i] - grid[i0][i + 1]) / deltaX;
  double b = (grid[i0][i - 1] - grid[i0][i]) / deltaX;
  double sigma = MinMod(a, b);
  return sigma;
}

/** slope for superbee slope limiter method */
double Particles::GetSuperBeeSlope(int i, double deltaX, int i0) {
  double a = (grid[i0][i - 1] - grid[i0][i]) / deltaX;
  double b = (grid[i0][i] - grid[i0][i + 1]) / deltaX;

  double sigma1 = MinMod(a, 2. * b);
  double sigma2 = MinMod(2. * a, b);

  double sigma = MaxMod(sigma1, sigma2);
  return sigma;
}


double Particles::LaxWendroffSlope(int i, double deltaX, int i0) {
  double sigma = (grid[i0][i] - grid[i0][i+1])/ deltaX;
  return sigma;
}


/** MaxMod function for slope limiters */
double Particles::MaxMod(double a, double b) {
  if (a * b > 0. && fabs(a) < fabs(b))
    return b;
  else if (a * b > 0. && fabs(a) >= fabs(b))
    return a;
  else
    return 0.;
}

/** MinMod function for slope limiters */
double Particles::MinMod(double a, double b) {
  if (a * b > 0. && fabs(a) < fabs(b))
    return a;
  else if (a * b > 0. && fabs(a) >= fabs(b))
    return b;
  else
    return 0.;
}

/**
 * Wrapper function to calculate the grid only in a specified time interval dT
 * = T2-T2
 * This is especially useful for the creation of time-series of spectra and
 * necessary
 * as the 'ComputeGrid' function only stores two spectra (the actual one and
 * that in the
 * time step before) due to memory storage reasons.
 */
void Particles::ComputeGridInTimeInterval(double T1, double T2, string type, 
                                          int bins, bool ICRESET) {
  
  if (TminExternal) Tmin = TminExternal;
  if (T1 <= Tmin) {
    cout << "Particles::ComputeGridInTimeInterval: T1<=internal min time (" << T1
         << "<=" << Tmin << "). set it artificially to " << Tmin << "yrs."
         << endl;
    T1 = Tmin * 1.001;
  }
  if (T1 >= T2) {
    cout << "Particles::ComputeGridInTimeInterval: T1>=T2 (" << T1
         << ">=" << T2 << "). Exiting."
         << endl;
    exit(1);
  }
  SetSolverMethod(0);
  if(!grid.size()) {
    Age = T1;
    CalculateParticleSpectrum(type, bins);
    fUtils->Clear2DVector(ParticleSpectrum);
  }
  // if ICRESET == true, this will re-create the IC loss lookup. This is important
  // when the radiation fields change during iteration. Danger: this is very slow!
  // TODO: implement 2D lookups for radiation fields, photon density vs. energy and 
  // time, to accommodate dynamical radiation fields.
  if(ICRESET==true) { 
    CalculateConstants();
  }
  ComputeGrid(energyAxis, T1, T2, yr_to_sec*(T2 - T1)/1e3);
  return;
}

/**
 * This calculates the particle spectrum without losses. It is essentially
 * N(E,t) = Q(E,t) * t
 */
void Particles::CalcSpecSemiAnalyticNoELoss() {
  if (Tmin >= Age) {
    cout << "CalcSpecSemiAnalyticNoELoss: Tmin is larger/equal "
            "than source age... Exiting" << endl;
    return;
  }
  double totallum = 0.;
  if(!std::isnan(LumConstant)) totallum = LumConstant*Age;
  else if(LumVector.size()) {
    gsl_interp_accel_reset(accLum);
    if(gsl_spline_eval_integ_e(LumLookup, Tmin, Age, accLum, &totallum))
      totallum = 0.;
  }
  else return;
  double LumConstantTemp = LumConstant;
  LumConstant = totallum*yr_to_sec;
  SetMembers(Age);

  double logstep = (log10(eMax) - log10(Emin)) / ebins;
  for (double e = Emin; e < eMax; e = pow(10., log10(e) + logstep)) {
    double val = SourceSpectrum(e);
    if(std::isnan(val) || std::isinf(val) || !val)
      continue;
    fUtils->TwoDVectorPushBack(e,val,ParticleSpectrum);
  }
  LumConstant = LumConstantTemp;
  return;
}

/**
 * Semi-analytical solution, only applicable for time-independent (i.e. constant)
 * energy losses. A nice derivation can be found in 
 * Atoyan & Aharonian MNRAS 302, 253-276 (1999). Particle escape is not implemented
 * here. TODO: Implement particle escape?
 */
void Particles::CalcSpecSemiAnalyticConstELoss() {
  if (Tmin >= Age) {
    cout << "Particles::CalcSpecSemiAnalyticConstELoss: Tmin is larger/equal "
            "than source age... Exiting" << endl;
    return;
  }

  /* info writeout. Disable it by using 'ToggleQuietMode()' */
  if (!Type && QUIETMODE == false)
    cout << "** Evolving Electron Spectrum:" << endl;
  else if (Type == 1 && QUIETMODE == false)
    cout << "** Evolving Proton Spectrum:" << endl;

  fPointer IntFunc = NULL;

  fUtils->Clear2DVector(ParticleSpectrum);

  double logstep = (log10(eMax) - log10(Emin)) / ebins;

  double e, lossrate, val = 0.;
  e = lossrate = val = 0.;
  int tt;

  Tmin = pow(10.,vETrajectory[0][0]);
  SetMembers(Age);
  // this is the time integral in Eq. (18) of the above mentioned paper
  for (e = Emin, tt = 0; e < eMax; e = pow(10., log10(e) + logstep), tt++) {
    if (QUIETMODE == false)
      cout << "    " << (int)(100. * tt / ebins) -1 << "\% done\r" << std::flush;
    lossrate = EnergyLossRate(e);
    if(lossrate <= 0.) continue;
    IntFunc = &Particles::SemiAnalyticConstELossIntegrand;
    val = Integrate(IntFunc, &e, log10(Tmin), log10(Age),
                    0.1*integratorTolerance,kronrodrule) / lossrate;
    if(val) fUtils->TwoDVectorPushBack(e,val,ParticleSpectrum);
  }
  return;
}

/**
 * Integrand of the time-integration of the semi-analytical solution. 
 * This is an integration over the time energy trajectory of the electron.
 *
 * @param T : time
 * @param par : additional parameters (energy)
 */ 
double Particles::SemiAnalyticConstELossIntegrand(double T, void *par) {
  if (energyTrajectoryInverse == NULL || energyTrajectory == NULL) {
    cout << "Particles::SemiAnalyticConstELossIntegrand: Couldn't calculate "
            "particle energy trajectory lookup. Exiting program." << endl;
    exit(1);
  }
  double tdash, E, Enow;
  Enow = *(double *)par;
  if(gsl_spline_eval_e(energyTrajectoryInverse, log10(Enow), accTrInv, &tdash)) 
    return 0.;
  T = pow(10.,T);
  tdash = pow(10.,tdash);

  if (tdash - T < Tmin || Age - T < Tmin) return 0.;
  if (gsl_spline_eval_e(energyTrajectory, log10(tdash - T), accTr, &E))
    return 0.;
  E = pow(10., E);
  if (E > eMax || E <= 0.) return 0.;
  SetMembers(Age - T);
  double dE = EnergyLossRate(E);
  double Sp = SourceSpectrum(E);
  // the number at the end is ln(10), needed because integration is in energy log space.
  return dE * Sp * T * yr_to_sec * 2.302585093; 
}

/**
 * Calculates the energy trajectory of an electron assuming time-constant energy losses.
 */
void Particles::CalculateEnergyTrajectory(double TExt) {
  if (TminInternal < 0.) {
    cout << "Particles::CalculateEnergyTrajectory: Calculate internal Tmin "
            "first by running DetermineLookupStartingTime(). Exiting." << endl;
    return;
  }
  bool UPDATE = false;
  if(vETrajectory.size()) UPDATE = true;
  fUtils->Clear2DVector(vETrajectory);

  double E, Edot, dt, T;
  T = (TExt > TminInternal) ? TExt : TminInternal;
  T *= 1.1;
  SetMembers(T);
  E = eMax;
  while (E >= Emin) {
    if(E<=0.) {
      cout << "Particles::CalculateEnergyTrajectory: Energy smaller 0. Exiting."
           << endl;
      return;
    }
    fUtils->TwoDVectorPushBack(log10(T),log10(E),vETrajectory);
    Edot = EnergyLossRate(E);
    dt = 1.e-3 * E / Edot;
    E -= dt * Edot;
    T += dt / yr_to_sec;
    if(T>TmaxInternal || Edot < 0.) break;
    SetMembers(T);
  }
  unsigned int size = vETrajectory.size();
  if(size<3) return;
  if(fabs(vETrajectory[0][0]/vETrajectory[size-1][0]-1.) < 1.e-3) return;
  // make an inverse of the energy loss trajectory, holding x=E(t),y=t
  vector< vector<double> > vETrajectoryInverse;
  for (unsigned int i = size-1; i > 0; i--)
    fUtils->TwoDVectorPushBack(vETrajectory[i][1],vETrajectory[i][0],
                               vETrajectoryInverse);
  if(UPDATE == true) {
    gsl_spline_free(energyTrajectory);
    gsl_spline_free(energyTrajectoryInverse);
  }
  energyTrajectory = fUtils->GSLsplineFromTwoDVector(vETrajectory);
  energyTrajectoryInverse = fUtils->GSLsplineFromTwoDVector(vETrajectoryInverse);
  return;
}



void Particles::SetType(string type) {
  if (!type.compare("electrons"))
    Type = 0;
  else if (!type.compare("protons"))
    Type = 1;
  else
    cout << "Particles::SetType: Particle type not recognized! "
    		"Please specify supported particle "
            "species! (electrons or protons)" << endl;
  return;
}

double Particles::SourceSpectrumWrapper(double E, void *par) {
  return SourceSpectrum(E);
}

/**
 * Determine the common time boundaries of the provided lookups for the
 * evolution of parameters (Source Luminosity,B-Field, Ambient density etc.).
 * If everything is set to constant values, i.e. in a stationary scenario,
 * these boundaries will be set to {TminInternal,TmaxInternal} = {1.e-4,1.e100}yrs.
 * TODO: Shorten this!!
 */
void Particles::DetermineLookupTimeBoundaries() {

  // determine starting time and energy
  double lumtmin, emaxtmin, ntmin, btmin, rtmin, vtmin, esctmin, lumtmax,
      emaxtmax, ntmax, btmax, rtmax, vtmax, esctmax,custspinjtmin,custspinjtmax;
  lumtmin = emaxtmin = ntmin = btmin = rtmin = vtmin = esctmin = custspinjtmin 
    = TminInternal; //std::cout ACHTUNG; hier stand 1e-10
  lumtmax = emaxtmax = ntmax = btmax = rtmax = vtmax = esctmax = custspinjtmax
    = TmaxInternal;//std::cout ACHTUNG; hier stand 1e100

  if (LumVector.size()) lumtmin = LumVector[0][0];
  if (eMaxVector.size()) emaxtmin = eMaxVector[0][0];
  if (NVector.size()) ntmin = NVector[0][0];
  if (BVector.size()) btmin = BVector[0][0];
  if (RVector.size()) rtmin = RVector[0][0];
  if (VVector.size()) vtmin = VVector[0][0];
  if (EscapeTimeEnergyTimeEvolutionVector.size()) 
        esctmin = EscapeTimeEnergyTimeEvolutionVector[0][0];
  if (escapeTimeVectorTdep.size()) 
        esctmin = escapeTimeVectorTdep[0][0];
  if (CustomInjectionSpectrumTimeEvolutionVector.size()) 
        custspinjtmin = CustomInjectionSpectrumTimeEvolutionVector[0][0];

  if (LumVector.size()) lumtmax = LumVector[LumVector.size() - 1][0];
  if (eMaxVector.size()) emaxtmax = eMaxVector[eMaxVector.size() - 1][0];
  if (NVector.size()) ntmax = NVector[NVector.size() - 1][0];
  if (BVector.size()) btmax = BVector[BVector.size() - 1][0];
  if (RVector.size()) rtmax = RVector[RVector.size() - 1][0];
  if (VVector.size()) vtmax = VVector[VVector.size() - 1][0];
  if (EscapeTimeEnergyTimeEvolutionVector.size()) esctmax = 
    EscapeTimeEnergyTimeEvolutionVector[EscapeTimeEnergyTimeEvolutionVector.size() - 1][0];
  if (escapeTimeVectorTdep.size()) 
        esctmax = escapeTimeVectorTdep[escapeTimeVectorTdep.size() - 1][0];
  if (CustomInjectionSpectrumTimeEvolutionVector.size()) custspinjtmax = 
    CustomInjectionSpectrumTimeEvolutionVector[CustomInjectionSpectrumTimeEvolutionVector.size() - 1][0];

  double T0, T1, T2, T3, T4, T5, T;

  // first determine highest common lower boundary time of the lookups
  if (!emaxtmin && !ntmin && !btmin && !rtmin && !vtmin && !esctmin && 
      !custspinjtmin)
    T = TminInternal;
  else {

    T0 = rtmin > vtmin ? rtmin : vtmin;
    T1 = ntmin > T0 ?  ntmin : T0;
    T2 = btmin > T1 ?  btmin : T1;
    T3 = esctmin > T2 ?  esctmin : T2;
    T4 = lumtmin > T3 ?  lumtmin : T3;
    T5 = emaxtmin > T4 ?  emaxtmin : T4;
    T = custspinjtmin > T5 ?  custspinjtmin : T5;

  }
  TminInternal = T;

  // now do the similar thing to get the lowest common higher lookup boundary
  if (!emaxtmin && !ntmin && !btmin && !rtmin && !vtmin)
    T = TmaxInternal;
  else {
    T0 = rtmax < vtmax ? rtmax : vtmax;
    T1 = ntmax < T0 ?  ntmax : T0;
    T2 = btmax < T1 ?  btmax : T1;
    T3 = esctmax < T2 ?  esctmax : T2;
    T4 = lumtmax < T3 ?  lumtmax : T3;
    T5 = emaxtmax < T4 ?  emaxtmax : T4;
    T = custspinjtmax < T5 ?  custspinjtmax : T5;

  }
  TmaxInternal = T;

  return;
}

/*
 * Return a particle SED E vs E^2dN/dE (TeV vs erg)
 */
vector<vector<double> > Particles::GetParticleSED() {
  vector<vector<double> > v;
  for (unsigned int i = 0; i < ParticleSpectrum.size(); i++) {
    double E = ParticleSpectrum[i][0];
    double ETeV = E / TeV_to_erg;
    double N = ParticleSpectrum[i][1];
//    if (!N) continue;
    fUtils->TwoDVectorPushBack(ETeV,E * E * N,v);
  }
  return v;
}

/*
 * Return total energy in particles. Input energy bounds in erg!
 * FIXME:This function sucks, have to think of something better!
 */
double Particles::GetParticleEnergyContent(double E1, double E2) {

  if(E2 <= E1 && E1 && E2) {
    cout << "Particles::GetParticleEnergyContent: upper energy bound "
            "equal/lower than lower energy bound: " << E2 << " <= " << E1 <<
            " Returning 0 value." << endl;
    return 0.;
  }
  
  vector<vector<double> > v;
  for (unsigned int i = 0; i < ParticleSpectrum.size(); i++) {
    double E = ParticleSpectrum[i][0];
    double N = ParticleSpectrum[i][1];
    if (!N) continue;
    fUtils->TwoDVectorPushBack(E,E * N,v);
  }
  
  if(E1 < ParticleSpectrum[0][0] || !E1)
    E1 = ParticleSpectrum[0][0];
  if(E2 > ParticleSpectrum[ParticleSpectrum.size()-1][0] || !E2)
    E2 = ParticleSpectrum[ParticleSpectrum.size()-1][0];
  return fUtils->Integrate(v,E1,E2);
}

void Particles::SetIntegratorMemory(string mode) {
  if(!mode.compare("light")) gslmemory=1000;
  else if(!mode.compare("normal")) gslmemory=5000;
  else if(!mode.compare("heavy")) gslmemory=10000;
  else {
    cout << "Particle::SetIntegratorMemory: Set valid mode. Possibilities: "
               "  light - normal - heavy. Default is 'normal'. " << endl;
  }
  return;
}

/*
 * Integration function using the GSL QAG functionality
 *
 */
Particles::fPointer Particles::_funcPtr;
Particles *Particles::_partPtr;

double Particles::evaluate(double x, void* params) {
  return (_partPtr->*_funcPtr)(x, params);
}
double Particles::Integrate(fPointer f, double *x, double emin, double emax,
                            double tolerance, int kronrodrule) {
  double integral, error;
  /* no comment */
  fPointer ftemp = _funcPtr;
  gsl_function F;
  _funcPtr = f;
  _partPtr = this;
  F.function = &Particles::evaluate;
  F.params = x;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(gslmemory);
  int val = gsl_integration_qag(&F, emin, emax, 0, tolerance, gslmemory, kronrodrule,
                          w, &integral,&error);
  gsl_integration_workspace_free(w);
  _funcPtr = ftemp;
  if (val)  return 0.;
  else return integral;
}

/*
 * Set a custom injection spectrum. Input is a 2D vector with 2 colums, first
 * column holds energy in erg. 
 * mode = 0 -> custom injection spectrum lookup:
 * second column holds differential rate in 
 * erg^-1 s^-1. If this is specified, otherwise defined spectral parameters 
 * like luminosity lookups and specified spectral index will be ignored.
 * NOTE: you still have to specify emax (either via lookup or constant value).
 * mode = 1 -> energy dependent particle escape time lookup:
 * second column holds particle escape time in seconds
 */
void Particles::SetCustomEnergylookup(vector< vector<double> > vCustom,
                                      int mode){
  if(!vCustom.size()) {
    cout << "Particles::SetCustomEnergylookup: Input vector empty."
            "Exiting." << endl;
    return;
  }
  if(vCustom[0].size() != 2) {
    cout << "Particles::SetCustomEnergylookup:This function supports "
            "only input vectors with 2 columns, energy vs. differential "
            "injection rate / escape time. The vector you specified has " << 
            vCustom[0].size() << " columns. Exiting." << endl;
    return;
  }

  if(std::isnan(EmaxConstant)) EmaxConstant = vCustom[vCustom.size()-1][0];
  if(!EminConstant) EminConstant = vCustom[0][0];
  
  vCustom = fUtils->RemoveZeroEntries(vCustom);
  if(!MinELookup) MinELookup = vCustom[0][0];
  if(!MaxELookup) MaxELookup = vCustom[vCustom.size()-1][0];
  if(!mode) {
    // if no explicit luminosity is specified, use the norm of the custom injection spectrum
    CustomInjectionNorm = fUtils->EnergyContent(vCustom); 
    if(std::isnan(LumConstant) && !LumVector.size()) {
       LumConstant = CustomInjectionNorm;
    }
  }
  vCustom = fUtils->VectorAxisLogarithm(vCustom,0);
  vCustom = fUtils->VectorAxisLogarithm(vCustom,1);
  if(!mode) {  
    CustomInjectionSpectrum = fUtils->GSLsplineFromTwoDVector(vCustom);
  }
  else if(mode == 1) {
    escapeTimeLookupEdep = fUtils->GSLsplineFromTwoDVector(vCustom);
  }
  else {
    cout << "Particles::SetCustomEnergylookup:"
            "You specified unsupported mode " << 
            mode << ". Supported options are mode = 0 -> custom time-dependent"
            " injection spectrum, mode = 1 -> time and energy dependent escape "
            " time. Doing nothing." << endl;
  }
  return;    
}

/*
 * Set a custom injection spectrum or escape time. Input is a 2D vector with 3 colums.
 * mode 0 = custom injection spectrum: first column holds time in yrs, second 
 * energy in erg, third holds differential rate in 
 * erg^-1 s^-1. If this is specified, otherwise defined spectral parameters 
 * like luminosity lookups and specified spectral index will be ignored.
 * NOTE: you still have to specify emax (either via lookup or constant value).
 * mode 1 = time and energy dependent escape time:
 * first column holds time in yrs, second 
 * energy in erg, third holds escape time in seconds 
 */
void Particles::SetCustomTimeEnergyLookup(vector< vector<double> > vCustom, int mode){

  if(!vCustom.size()) {
    cout << "Particles::SetCustomTimeEnergyLookup: Input vector empty."
            "Exiting." << endl;
    return;
  }
  if(vCustom[0].size() != 3) {
    cout << "Particles::SetCustomTimeEnergyLookup:"
            "This function supports "
            "only input vectors with 3 columns, time vs energy vs. differential "
            "injection rate / escape time. The vector you specified has " << 
            vCustom[0].size() << " columns. Exiting." << endl;
    return;
  }
  
  vCustom = fUtils->VectorAxisLogarithm(vCustom,1);
  vCustom = fUtils->VectorAxisLogarithm(vCustom,2);
  if(!mode) {
    CustomInjectionSpectrumTimeEvolution = 
                      fUtils->TwoDsplineFromTwoDVector(vCustom,
                                                       TminLookup,TmaxLookup,
                                                       MinELookup,MaxELookup);
    CustomInjectionSpectrumTimeEvolutionVector = vCustom;
  }
  else if(mode == 1) {
    EscapeTimeEnergyTimeEvolution = 
                      fUtils->TwoDsplineFromTwoDVector(vCustom,
                                                       TminLookup,TmaxLookup,
                                                       MinELookup,MaxELookup);
    EscapeTimeEnergyTimeEvolutionVector = vCustom;

  }
  else {
    cout << "Particles::SetCustomTimeEnergyLookup:"
            "You specified unsupported mode " << 
            mode << ". Supported options are mode = 0 -> custom time-dependent"
            " injection spectrum, mode = 1 -> time and energy dependent escape "
            " time. Exiting." << endl;
    return;
  }
  TminInternal = TminLookup;
  DetermineLookupTimeBoundaries();
  MinELookup = pow(10.,MinELookup);
  MaxELookup = pow(10.,MaxELookup);
  if(std::isnan(EmaxConstant)) {
    EmaxConstant = MaxELookup;
  }
  if(!EminConstant) {
    EminConstant = MinELookup;
  }

  return;    
}

/** Add a greybody distribution of target photons,
 * which is used in the
 * IC cooling process in this class.
 */
void Particles::AddThermalTargetPhotons(double T, double edens, int steps) {
  fRadiation->AddThermalTargetPhotons(T, edens, steps);
}

/** Reset photon target field i with a greybody distribution,
 * which is used in the
 * IC cooling process in this class.
 */
void Particles::ResetWithThermalTargetPhotons(int i, double T, double edens, int steps) {
  fRadiation->ResetWithThermalTargetPhotons(i, T, edens, steps);
}

/** Add an arbitray distribution of target photons,
 * which is used in the
 * IC cooling process in this class
 * This requires as input a 2D vector of format:
 *              ~~~    energy[erg] number_density   ~~~
 */
void Particles::AddArbitraryTargetPhotons(vector<vector<double> > PhotonArray) {
  fRadiation->AddArbitraryTargetPhotons(PhotonArray);
}

/** Reset photon target field i by a an arbitrary distribution of target photons,
 * which is used in the
 * IC cooling process in this class.
 */
void Particles::ResetWithArbitraryTargetPhotons(int i,vector<vector<double> > PhotonArray) {
  fRadiation->ResetWithArbitraryTargetPhotons(i, PhotonArray);
}

/** Import target photons from file. File has to be in ASCII format, namely:
 *              ~~~    energy[eV] number_density   ~~~
 */
void Particles::ImportTargetPhotonsFromFile(const char *phFile) {
  fRadiation->ImportTargetPhotonsFromFile(phFile);
}

/** Reset photon target field i by target photons from file,
 * which is used in the
 * IC cooling process in this class.
 */
void Particles::ResetWithTargetPhotonsFromFile(int i,const char *phFile) {
  fRadiation->ResetWithTargetPhotonsFromFile(i,phFile);
}


/** Add SSC target photons.
 * This function calls the Synchrotron code in this class.
 * If 'UPDATE' is 'true'
 * then recalculate the synchroton target field and
 * replace the previous one by this updated field.
 * DANGER: for the 'UPDATE' option to work, the SSC field
 * must be the last entry in the 'TargetPhotonGraphs' vector!
 * It uses Atoyan&Aharonian1996: MNRAS, Volume 278, Issue 2, pp. 525-541
 */
void Particles::AddSSCTargetPhotons(int steps) {
  fRadiation->SetElectrons(ParticleSpectrum);
  fRadiation->SetBField(BField);
  fRadiation->AddSSCTargetPhotons(R/pc_to_cm,steps);
}

/** Reset photon target field i with SSC target photons.
 *
 */
void Particles::ResetWithSSCTargetPhotons(int i,int steps) {
  fRadiation->SetElectrons(ParticleSpectrum);
  fRadiation->SetBField(BField);
  fRadiation->ResetWithSSCTargetPhotons(i, R/pc_to_cm, steps);
}

void Particles::ClearTargetPhotons() {
  fUtils->Clear2DVector(ICLossVector);
  fRadiation->ClearTargetPhotons();
}


/** remove the latest component in
 * TotalTargetPhotonVector and recompute
 * the total target photon spectrum
 */
//void Particles::RemoveLastICTargetPhotonComponent() {
//  fRadiation->RemoveLastICTargetPhotonComponent();
//  fRadiation->CreateICLossLookup(iclosslookupbins);
//  SetLookup(fRadiation->GetICLossLookup(), "ICLoss");
//}


void Particles::CheckSanityOfTargetPhotonLookup() {
  fRadiation->CheckSanityOfTargetPhotonLookup();
}


vector<vector<double> > Particles::GetTargetPhotons(int i) {
  return fRadiation->GetTargetPhotons(i);
}

/*
 * Returns the energy loss rate of a particle with energy E depending on
 * the loss process. Energy is given in erg.
 */
double Particles::GetEnergyLossRate(double E, string type, double age) {

    if(age) {
        SetMembers(age);
    }
    CalculateConstants();
    double sum = EnergyLossRate(E);
    if (!type.compare("sum")) return sum;     
    if (!type.compare("synchrotron") ) return synchl;
    else if (!type.compare("inverse_compton") ) return icl;
    else if (!type.compare("adiabatic_losses") ) return adl;
    else if (!type.compare("bremsstrahlung") ) return bremsl;
    else if (!type.compare("ionization") ) return ionization;
    else if (!type.compare("ppcol") ) return ppcol;
    
    else {
        cout<<"Particles::GetEnergyLossRate: Loss mechanism "<<type<<" not known."
              " returning zero."<<endl;
        return 0.;
    }
}

/**
 * Computes and returns either the energy loss rate of the particle or the
 * cooling time via setting the \a TIMESCALE boolean to false or true
 * respectively
 *
 * @param epoints : vector of energies of the particle
 * @param type : process for which you want the information
 * @param age : age at which you want to compute
 * @param TIMESCALE : 1 for cooling time; 0 for energy loss rate
 *
 * @return v as a vector of tuples (energy, value)
 */
vector< vector<double> > Particles::GetEnergyLossRateVector(vector<double> epoints,
                                                        string type, double age, 
                                                               bool TIMESCALE) {
    if(!epoints.size()) {
        cout << "Particles::GetEnergyLossRateVector: you input a vector of size 0. "
                "Exiting." << endl;
        exit(0);
    }
    if(age || !ICLossVector.size()) {
        SetMembers(age);
    }
    CalculateConstants();
    vector< vector<double> > v;
    for(unsigned int i=0; i < epoints.size(); i++) {
        double energy = epoints[i];
        double val = GetEnergyLossRate(energy,type);
        if(TIMESCALE == true) val = energy / val / yr_to_sec;
        fUtils->TwoDVectorPushBack(energy,val,v);
    }
    return v;
}

/**
 * Returns the injection spectrum:
 *
 * @param epoints : vector of energies to return the spectrum (erg units)
 * @param age : age at which to return the spectrum (if age not
 *                            given, then computes as the injection)
 * @param SED : true if the particle spectrum should be returned
 *                            as E^2dN/dE
 * @return v : vector of tuples (Energy, spectrum)
 *                               energy returned in TeV.
 */
vector< vector<double> > Particles::GetInjectionSpectrumVector(
                                                      vector<double> epoints,
                                                      double age, bool SED) {
    if(!epoints.size()) {
        cout << "Particles::GetInjectionSpectrumVector: you input a vector of "
                "size 0. Exiting." << endl;
        exit(0);
    }
    if(age) SetMembers(age);
    CalculateConstants();
    vector< vector<double> > v;
    for(unsigned int i=0; i < epoints.size(); i++) {
        double energy = epoints[i];
        double val = SourceSpectrum(energy);
        if(SED == true) {
            val *= energy * energy;
            energy /= TeV_to_erg;
        }
        fUtils->TwoDVectorPushBack(energy,val,v);
    }
    return v;
}

/**
 * Wrapper around other get functions
 */
vector< vector<double> > Particles::GetQuantityVector(vector<double> epoints,
                                                      double age, string mode) {
    if(!epoints.size()) {
        cout << "Particles::GetQuantityVector: you input a vector of size 0 when "
                "trying to extract "<<mode<<". Exiting." << endl;
        exit(0);
    }
    if(age) SetMembers(age);
    CalculateConstants();
    vector< vector<double> > v;

    for(unsigned int i=0; i < epoints.size(); i++) {
        double energy = epoints[i];
        double val = 0.;
        if(!mode.compare("energy_loss_rate")) {
            val = EnergyLossRate(energy);}
        else if(!mode.compare("cooling_time_scale")) {
            val = energy / EnergyLossRate(energy) / yr_to_sec;}
        else if(!mode.compare("injection_spectrum")) {
            val = SourceSpectrum(energy);}
        else if(!mode.compare("injection_sed")) {
            val = SourceSpectrum(energy)*energy*energy;
            energy /= TeV_to_erg;}
        else if(!mode.compare("escape_time")) {
            val = EscapeTime(energy, age) / yr_to_sec;}
        else {
            cout << "Particles::GetQuantityVector: '"<< mode <<"' is not a valid "
                "option. Allowed are 'energy_loss_rate','cooling_time_scale',"
                "'injection_spectrum','injection_sed','escape_time'. Exiting." <<endl;
            exit(0);}
    
        fUtils->TwoDVectorPushBack(energy,val,v);
    }
    return v;
}



/*
 * Function under construction! Use at own peril!
 * TODO: The handling of the SSC field should be done better. Not sure all the cases
 * are correctly taken into account.
 * To work the SSC must be the last target field added to list.
 */
vector<double> Particles::CalculateSSCEquilibrium(double tolerance, int bins) {
  if (TminExternal) Tmin = TminExternal;
  double r_orig = R;
  double r_start = 10.*R < 1. ? 1. : 10.*R;
  int steps = 10;
  double delta_log_r = (log10(r_start) - log10(r_orig)) / steps;
  double delta_t = (Age-Tmin) / steps;
  // returns Radiation::RADFIELD_COUNTER
  double target_field_count = fRadiation->GetTargetFieldCount();
  double e_p = 0.;
  double e_p_0 = 0.;
  double precision = 1e100;
  if (QUIETMODE == false) {
    cout << "\r";
    cout << "____________________________________" << endl;
    cout << ">> CALCULATING SSC EQUILIBRIUM STATE " << endl;
    cout << "   -> pre-iteration for speedup     " << endl;
  }
  // this first loop calculates spectrum on ever-decreasing radii until real
  // source radius is reached. The idea is to 'cool away' the highest energy 
  // electrons first in a lower intensity synchrotron target field (due to the
  // larger radii) before the actual iteration on the real radius starts. This
  // can considerably speed up the process.
  vector<double> iter;
  CalculateElectronSpectrum(bins);
  double log_r,t;
  int i;

  for(log_r = log10(r_start), t=Tmin, i=0; log_r >= log10(r_orig) && t<=Age && i<steps; log_r -= delta_log_r, t+= delta_t,i++) {
    
    if (QUIETMODE == false) {
     cout << "    " << (int)(100. * (i+1) / steps) << "\% done\r" << std::flush;
    }
    R = pow(10.,log_r);
    if (log_r == log10(r_start)) AddSSCTargetPhotons();
    // -1 because the counter starts from 1, while the indexing of the target
    // fields starts from 0
    else ResetWithSSCTargetPhotons(target_field_count-1);
    CalculateElectronSpectrum(bins);
    e_p = GetParticleEnergyContent();
    iter.push_back(e_p);
    e_p_0 = e_p;
  }   
  cout << "   -> precision reached / goal:     " << endl;
  R = r_orig;
  while(1){
    ResetWithSSCTargetPhotons(target_field_count-1);
    SetMembers(Age);
    CalculateElectronSpectrum(bins);
    e_p = GetParticleEnergyContent();
    iter.push_back(e_p);
    precision = fabs(e_p/e_p_0 - 1.);
        cout << " >> " << precision << " / " << tolerance <<"\r"<<std::flush;
    if (precision < tolerance) break;
    e_p_0 = e_p;
  }

    cout << "   -> DONE!     " << endl;
  return iter;

}

