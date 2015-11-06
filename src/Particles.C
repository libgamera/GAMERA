#include "Particles.h"

Particles::Particles() {
  /* Default values */
  eElectronMax = 100000.;
  SpectralIndex = 2.000;
  Tmin = 1.;
  CutOffFactor = 1000.;
  TminInternal = 1.e-10;
  TmaxInternal = 1.e100;
  EminInternal = 1.e-3;
  energyMarginFactor = 1.e-3;
  energyMargin = 1.;
  ebins = 100;
  SACELI_Told = 0.;
  R = -1.;
  V = -1.;
  gslmemory=5000;
  QUIETMODE = false;
  EminConstantForNormalisationOnly = false;
  DEBUG = false;
  ICLossLookup = NULL;
  LumLookup = NULL;
  NLookup = NULL;
  eMaxLookup = NULL;
  BFieldLookup = NULL;
  RLookup = NULL;
  VLookup = NULL;
  escapeTimeLookup = NULL;
  integratorTolerance = 5.e-2;
  ParticleSpectrum.clear();
  ParticleSpectrum.resize(0);
  ICLossVector.resize(0);
  LumVector.resize(0);
  NVector.resize(0);
  BVector.resize(0);
  eMaxVector.resize(0);
  escapeTimeVector.resize(0);
  RVector.resize(0);
  VVector.resize(0);

  SpectralIndex2 = eBreak = TminConstant = adLossCoeff = EminConstant =
      eMaxConstant = 0.;
  escapeTime = escapeTimeConstant = EnergyAxisLowerBoundary =
      EnergyAxisUpperBoundary = 0.;
  eBreakS2 = eBreak2mS2 = eBreakS = eBreak2mS = emin2mS2 = emin2mS =
      emineBreak2mS2 = 0.;
  eBreak2mSInd2 = emin2mSInd2 = emineBreak2mSInd2 = fS = fS2 = bremsl_epf =
      bremsl_eef = 0.;
  LumConstant = BConstant = NConstant = VConstant = RConstant = 0.;
  accIC = gsl_interp_accel_alloc();
  accLum = gsl_interp_accel_alloc();
  accN = gsl_interp_accel_alloc();
  accBField = gsl_interp_accel_alloc();
  acceMax = gsl_interp_accel_alloc();
  accescapeTime = gsl_interp_accel_alloc();
  accR = gsl_interp_accel_alloc();
  accV = gsl_interp_accel_alloc();
  accTr = gsl_interp_accel_alloc();
  accTrInv = gsl_interp_accel_alloc();
}
Particles::~Particles() {
  for (unsigned int i = 0; i < grid.size(); i++) grid[i].clear();
  grid.clear();
}

/*   ----    END OF BINARY SEARCHES     ----   */

/** Speed Hack: calculate some constant that are repeatedly needed in nested
 * loops
 */
void Particles::CalculateConstants() {
  double emin = Emin;
  if ((Kep && !Type) || Type == 1)
    emin = m_p;  // norm for protons or for electrons if defined per e/p
                 // fraction
  eBreakS2 = pow(eBreak, SpectralIndex2);
  eBreak2mS2 = pow(eBreak, 2. - SpectralIndex2);
  eBreak2mSInd2 = log(eBreak);
  eBreakS = pow(eBreak, SpectralIndex);
  eBreak2mS = pow(eBreak, 2. - SpectralIndex);
  emin2mS2 = pow(emin, 2. - SpectralIndex2);
  emin2mS = pow(emin, 2. - SpectralIndex);
  emin2mSInd2 = log(emin);
  emineBreak2mS2 = pow(emin / eBreak, 2. - SpectralIndex2);
  emineBreak2mSInd2 = log(emin / eBreak);
  fS = 1. / (2. - SpectralIndex);
  fS2 = 1. / (2. - SpectralIndex2);
  bremsl_epf = 3. * fineStructConst * sigma_T * c_speed * m_e / pi;
  bremsl_eef = (3. * fineStructConst * sigma_T * c_speed * m_e / (2. * pi));
  return;
}

/** fill the lookup holding the particle spectrum {E(erg) - N(erg^-1) */
void Particles::CalculateParticleSpectrum(string type, bool onlyprepare,
                                          bool dontinitialise) {

  if (!type.compare("electrons")) {
    Type = 0;
  } else if (!type.compare("protons")) {
    Type = 1;
  } else {
    cout << "Particles::FillParticleSpectrumLookup: What the fuck! Specify "
            "supported particle species! " << endl;
  }
  if (!QUIETMODE) {
    cout << "___________________________________" << endl;
    cout << ">> STARTING NEW PARTICLE EVOLUTION " << endl;
    if (Type == 0) cout << "   (-> electrons)     " << endl;
    if (Type == 1) cout << "   (-> protons)     " << endl;
  }
  /* reset Particle Lookup */
  Clear2DVector(ParticleSpectrum);

  if (Type == 1)
    Emin = m_p;  // norm for protons or for electrons if defined per e/p
                 // fraction
  else
    Emin = m_e;
  /* if lower energy bound is externally set
   * use this value for the normalisation of the spectrum
   */
  if (EminConstant>m_e) {
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

  if (TminConstant)
    Tmin = TminConstant;
  else if (Type == 1)
    Tmin = 1.e-3;
  else
    DetermineTMin(EminInternal, Tmin);
  if(Tmin<TminInternal) Tmin=TminInternal;
  /* get the upper energy boundary of the spectrum. This can either be
   * externally
   * set (if) or dynamically determined by the code (else).
   */
  if (eMaxConstant)
    eMax = eMaxConstant;
  else
    eMax = DetermineEmax(Tmin);

  /* apply a safe margin to the upper energy boundary
   * in order to prevent numerical effects at the upper end of the spectrum.
   */
  energyMargin = pow(-log(energyMarginFactor), 1. / CutOffFactor);
  eMax *= energyMargin;

  /* if emax falls below emin, return dummy vector with zeroes */
  if (eMax <= Emin) {
    cout << "Particles::FillParticleSpectrumLookup Whaat? eMax lower than "
            "Emin: eMax = " << eMax << " Emin = " << Emin << endl;
    return;
  }

  if (BVector.size() || NVector.size() || (RVector.size() && VVector.size())) {
    /* call numerical solver. */
    PrepareAndRunNumericalSolver(ParticleSpectrum, onlyprepare, dontinitialise);
  }
  else if (Type==0 && (BConstant || NConstant || VConstant)) {
    /* if B-Field, ambient density and speed are set externally thus constant.
     * This means constant energy losses due to Brems- and Synchrotronstrahlung
     * as well as adiabatic expansion (currently IC losses are always constant
     * in GAMERA.).
     */
    CalculateEnergyTrajectory();
    CalcSpecSemiAnalyticConstELoss();
  }
  else if (Type==1 && VConstant) {
    /* if B-Field, ambient density and speed are set externally thus constant.
     * This means constant energy losses due to Brems- and Synchrotronstrahlung
     * as well as adiabatic expansion (currently IC losses are always constant
     * in GAMERA.).
     */
    CalculateEnergyTrajectory();
    CalcSpecSemiAnalyticConstELoss();
  } else {
    CalcSpecSemiAnalyticNoELoss();
  }
  if (!QUIETMODE) {
    cout << ">> PARTICLE EVOLUTION DONE. EXITING." << endl;
    cout << endl;
  }
  return;
}

/** particle injection spectrum, assuming a power-law+super exponential Cut-off
 * supported injection Spectra:
 * - Power Law with exponential cut-off (pars: SpectralIndex,eMax - cut energy,
 * CutOffFactor - ~exp[-(e/emax)^cutOffFactorInternal])
 * - Broken power law with exponential cut-off (pars: SpectralIndex,eMax - cut
 * energy, eBreak - break energy, CutOffFactor -
 * ~exp[-(e/emax)^cutOffFactorInternal])
 * - exponential cut-offs can be mitigated by very high values of "CutOffFactor"
 */
double Particles::SourceSpectrum(double e) {
  return PowerLawInjectionSpectrum(e, eMax, 10. * eMax);
}

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
      cout << "Particles::PowerLawInjectionSpectrum: todo!" << endl;
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
  /* Broken power-law. If both SpectralIndex2 and a break energy was specified
   */
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

/** Set important class members to values at time t (as specified in
 * CRLumLookup). */
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
  Constants = {LumConstant,  NConstant,         BConstant,
               eMaxConstant, escapeTimeConstant};
  splines = {LumLookup, NLookup, BFieldLookup, eMaxLookup, escapeTimeLookup};
  accs = {accLum, accN, accBField, acceMax, accescapeTime};
  vals = {&Lum, &N, &BField, &eMax, &escapeTime};
  vs = {LumVector, NVector, BVector, eMaxVector, escapeTimeVector};
  for (unsigned int i = 0; i < Constants.size(); i++) {
    *vals[i] = -1.;
    if (Constants[i])
      *vals[i] = Constants[i];
    else if (splines[i] == NULL)
      continue;
    else
      *vals[i] = gsl_spline_eval(splines[i], t, accs[i]);
  }

  R = V = adLossCoeff = 0.;
  if(RConstant && VConstant) {
    R = VConstant*yr_to_sec + RConstant;
    V = VConstant;
  }

  else if(RConstant && !VConstant) {
    R = RConstant;
    V = 0.;
  }
  else {
    R = yr_to_sec*V;
    V = VConstant;
  }

  if(RVector.size() && t > RVector[0][0] &&
          t < RVector[RVector.size() - 1][0])
    R = gsl_spline_eval(RLookup, t, accR);

  if(VVector.size() && t > VVector[0][0] &&
          t < VVector[VVector.size() - 1][0])
    V = gsl_spline_eval(VLookup, t, accV);

  if (R && V) adLossCoeff = V / R;

  return;
}

void Particles::SetLookup(vector<vector<double> > v, string LookupType,
                          bool UPDATE) {
  int size = (int)v.size();
  if (!size) {
    cout << "Particles::SetLookup: lookup vector empty. Exiting." << endl;
    return;
  }
  double x[size];
  double y[size];
  for (int i = 0; i < size; i++) {
    x[i] = v[i][0];
    y[i] = v[i][1];
  }
  gsl_spline *ImportLookup = gsl_spline_alloc(gsl_interp_akima, size);
  gsl_spline_init(ImportLookup, x, y, size);

  vector<string> st{"ICLoss", "Luminosity", "AmbientDensity", "BField",
                    "Emax",   "EscapeTime", "Radius",         "Speed"};
  vector<gsl_spline **> spl{&ICLossLookup, &LumLookup,  &NLookup,
                            &BFieldLookup, &eMaxLookup, &escapeTimeLookup,
                            &RLookup,      &VLookup};
  vector<vector<vector<double> > *> vs{
      &ICLossVector, &LumVector,        &NVector, &BVector,
      &eMaxVector,   &escapeTimeVector, &RVector, &VVector};
  for (unsigned int i = 0; i < st.size(); i++) {
    if (!LookupType.compare(st[i])) {
      if (*spl[i] != NULL) {
        if (UPDATE == true) {
          gsl_spline_free(*spl[i]);
          *spl[i] = NULL;
        } else {
          cout << "Particles::SetLookup: " << st[i]
               << " lookup already set earlier. Exiting." << endl;
          return;
        }
      }
      *spl[i] = ImportLookup;
      *vs[i] = v;
      return;
    }
  }
  return;
}

/** Append CRLUMLOOKUP to CRLumLookup. Time order has to be right, so the
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
  vector<string> st = {"Luminosity", "AmbientDensity", "BField", "Emax",
                       "EscapeTime", "Radius",         "Speed"};
  vector<vector<vector<double> > > vs = {LumVector,  NVector,          BVector,
                                         eMaxVector, escapeTimeVector, RVector,
                                         VVector};
  for (unsigned int i = 0; i < st.size(); i++) {
    if (!LookupType.compare(st[i])) {
      if (vs[i][vs[i].size() - 1][0] >= v[0][0]) {
        cout << "Particles::ExtendCRLumLookup - WTF, the vector which to add ("
             << st[i] << ") to the existing one starts at earlier times than "
                         "the existing ones ends. Please keep time order in "
                         "the vector! exiting." << endl;
        return;
      }
      vs[i].insert(vs[i].end(), v.begin(), v.end());
      SetLookup(vs[i], LookupType, true);
      return;
    }
  }
  return;
}

/** calculates the energy loss rate at a given time t from the shock dynamics,
 * ambient photon and B-fields.
 */
double Particles::EnergyLossRate(double E) {
  double synchl = 0.;
  double icl = 0.;
  double adl = 0.;
  double bremsl = 0.;
  double bremsl_ep = 0.;
  double bremsl_ee = 0.;
  double gamma = (E + m_e) / m_e;
  if (gamma < 1.) gamma = 1.;
  double gamma2 = gamma * gamma;
  double p = sqrt(gamma2 - 1.);
  /* S-parameter. This is only the case in a pure hydrogen gas environment. For
   * more comlex mixture,
   * nuclear charge of the different gas species become important. See Haug2004
   */
  double S = N;
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

  /* Bremsstrahlung losses Haug+2004 */
  /* electron-proton bremsstrahlung */
  /* TODO: to be super self-consistent, calculate losses in a lookup, analogue
   * to the IC lookup */
  bremsl_ep = ((2. * gamma2 / 9. - 19. * gamma * p * p / 675. -
                0.06 * p * p * p * p / gamma) *
                   p * p * p / (gamma * gamma * gamma * gamma * gamma * gamma) +
               gamma * log(gamma + p) - p / 3.) *
              bremsl_epf * S * gamma2 / (gamma2 + p * p);

  /* electron-electron bremsstrahlung */
  bremsl_ee =
      N * bremsl_eef * (p * (gamma - 1.) / gamma) * (log(2. * gamma) - 1. / 3.);
  bremsl = bremsl_ep + bremsl_ee;

  /* in case of protons, only take adiabatic losses into account */
  if (Type) {
    synchl = 0.;
    icl = 0.;
    bremsl = 0.;
  }
  if (DEBUG == true) {
    synchl = 0.;
    bremsl = 0.;
    adl = 0.;
  }
  return synchl + icl + adl + bremsl;
}

/** Prepare axes for the numerical solver (and if onlyprepare==false) call the
 * solver
 */
void Particles::PrepareAndRunNumericalSolver(
    vector<vector<double> > &particlespectrum, bool onlyprepare,
    bool dontinitialise) {
  if (EnergyAxisUpperBoundary) {
    GetAxis(Emin, EnergyAxisUpperBoundary, ebins, energyAxis, true);
  } else if (EnergyAxisLowerBoundary) {
    GetAxis(EnergyAxisLowerBoundary, eMax, ebins, energyAxis, true);
  } else if (EnergyAxisLowerBoundary && EnergyAxisUpperBoundary) {
    GetAxis(EnergyAxisLowerBoundary, EnergyAxisUpperBoundary, ebins, energyAxis,
            true);
  } else {
    GetAxis(Emin, eMax, ebins, energyAxis, true);
  }
  CreateGrid();
  if (dontinitialise == false) SetInitialCondition(grid, energyAxis, Tmin);

  /* if onlyprepare=true, only the axes are initialised, but the grid is not
   * computed
   * until t=age. This is useful in conjunction with
   * ComputeGridInTimeInterval().
   */
  if (onlyprepare == false) ComputeGrid(grid, energyAxis, Tmin, Age, timeAxis);
  return;
}

/** create an axis that attributes each bin with a real value. */
void Particles::GetAxis(double min, double max, int steps, vector<double> &Axis,
                        bool logarithmic) {
  Axis.resize(steps);
  if (logarithmic == true) {
    min = log10(min);
    max = log10(max);
  }
  double binsize = (max - min) / steps;
  for (unsigned int i = 0; i < Axis.size(); i++) {
    Axis[i] = min + i * binsize;
  }

  return;
}

/** create the propagation grid out of 2 1D-vectors
 * energy in x-direction, time in y-direction
 */
void Particles::CreateGrid() {
  Clear2DVector(grid);
  grid.push_back(vector<double>());
  for (int i = 0; i < ebins; i++) {
    grid[grid.size() - 1].push_back(0.);
  }
  return;
}

/** determine the minimum time from where to start
 * the calculation. Electrons before that time are injected as
 * a single 'blob'. This time is derived from the requirement
 * that the blob has slid down to energies e.g. E<1GeV(EMIN) at t=Age.
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
    if (gsl_spline_eval_e(energyTrajectoryInverse, log10(emin), accTrInv,
                          &TMIN))
      continue;
    TMIN = pow(10., TMIN);
    if (TMIN > Age) break;
    tmin = pow(10., logt);
  }
  return;
}

/** Determine the maximum energy of particles between tmin and Age.
 * This energy is then used as upper boundary for the energy dimension of the
 * grid.
 */
double Particles::DetermineEmax(double tmin) {
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
  return eMaxHistory;
}

/** set initial condition (a.k.a. set the first energy vector at t=tmin of the
 * grid) */
void Particles::SetInitialCondition(vector<vector<double> > &Grid,
                                    vector<double> EnergyAxis,
                                    double startTime) {

  double t0 = startTime;
  SetMembers(t0);
  for (unsigned int i = 0; i < EnergyAxis.size(); i++) {
    double e = 0.5 * (pow(10, EnergyAxis[i]) + pow(10., EnergyAxis[i + 1]));
    Grid[0][i] = t0 * yr_to_sec * SourceSpectrum(e);
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
void Particles::ComputeGrid(vector<vector<double> > &Grid,
                            vector<double> EnergyAxis, double startTime,
                            double Age, vector<double> &TimeAxis,
                            double minTimeBin) {

  TimeAxis.push_back(startTime);
  TimeAxis.push_back(0.);
  double value = 0.;
  double quot = 0.;
  double t = 0.;
  double e1 = 0.;
  double e2 = 0.;
  double ebin = 0.;
  double tbin = 0.;
  double deltaE1 = 0.;
  double deltaE2 = 0.;
  double ElossRate_e1 = 0.;
  double ElossRate_e2 = 0.;
  unsigned int Esize = EnergyAxis.size() - 1;
  int tt = 1;
  int largestFilledBin = Esize;
  long int count = 0;
  /* append a new energy vector that will always hold the energy spectrum at the
   *  next time step and initialise it with zeroes.
   */
  Grid.push_back(vector<double>());
  for (unsigned int i = 0; i < EnergyAxis.size(); i++) {
    Grid[Grid.size() - 1].push_back(0.);
  }

  /* info writeout. Disable it by using 'ToggleQuietMode()' */
  if (!Type && QUIETMODE == false)
    cout << "** Evolving Electron Spectrum:" << endl;
  else if (Type == 1 && QUIETMODE == false)
    cout << "** Evolving Proton Spectrum:" << endl;

  /* main loop over time  */
  for (double T = startTime; T < Age; T += tbin / yr_to_sec) {
    /* Set Members (CR luminosity, B-field etc.) at time T */
    SetMembers(T);
    /* dynamically determine tbin size. This is a critical step
     * for the speed of the algorithm. Since the time step size is
     * proportional to Ebinsize(eMax)/Edot(eMax) and Edot ~ E^2,
     * eMax - or the largest relevant energy bin - should be chosen
     * as small as possible. Thus, don't use the energy at eMax, but
     * but rather of the largest filled energy bin.
     * This is 'largestFilledBin' which is time-dependent and is defined
     * in the next for-loop.
     * If an external emax is specified, always choose the highest energy bin
     * value
     */
    if (eMaxConstant) largestFilledBin = Esize - 1;
    e1 = pow(10., EnergyAxis[largestFilledBin - 1]);
    e2 = pow(10., EnergyAxis[largestFilledBin]);
    ebin = e2 - e1;

    /* the tbin size is then simply defined as deltaE/Edot_max */
    tbin = ebin / fabs(EnergyLossRate(e2));

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
    /* in case of protons, energy losses are negligible and the computation
     * is very fast. Thus, there is no reason not to use very fine time bins in
     * this case.
     */
    if (Type) tbin = 0.05 * ebin / EnergyLossRate(e2);
    /* if losses become small (e.g. degrading B-field, or low eMax), time bins
     * may become very large. This can become problematic if tbin << Age no
     * longer
     * holds. Then, replace tbin by minTimeBin (default: 1yr) if
     * tbin>minTimeBin.
     */
    if (tbin > minTimeBin) tbin = minTimeBin;

    /* negative time steps are not what we want! In this case shout out some
     * debug. */
    if (tbin < 0.) {
      cout << "Particles::ComputeGrid: ebin = " << ebin << " (e2,e1) = " << e2
           << "," << e1 << ") elossrate(" << e2 << ") = " << EnergyLossRate(e2)
           << " lastbin = " << largestFilledBin << endl;
    }

    /* if eMax drops below the lower energy bound of the grid, exit. */
    if (largestFilledBin <= 0) break;

    /* This is the new time! */
    t = T + 0.5 * tbin / yr_to_sec;
    /* update the Members at the new time. */
    Lum = 0.;
    SetMembers(t);

    /* iterate over the previous spectrum, stored in Grid[0] and calculate
     * Grid[1] from it. This is done using a standard, piece-wise linear
     * advection
     * scheme. Per default, also a slope limiter is implemented (MinMod method),
     * that preserves sharp edges in the particle spectrum rather than smearing
     * it out as in a pure donor-cell algorithm.
     */

    /* just for the first step (for speed) */
    ElossRate_e2 = EnergyLossRate(pow(10., EnergyAxis[0]));
    double particleCount = 0.;
    for (unsigned int i = 0; i < Esize; i++) {
      count++;
      value = 0.;
      quot = 0.;
      ebin = 0.;

      e1 = pow(10., EnergyAxis[i]);
      e2 = pow(10., EnergyAxis[i + 1]);

      /* The following block calculates the streaming of particles in an out
       * of energy bin 'i' by cooling. This component of increase / decrease of
       * particles in bin 'i' is caused by particles already present in the
       * last time step (Grid[0])
       */
      ebin = e2 - e1;
      quot = tbin / ebin;

      ElossRate_e1 = ElossRate_e2;
      ElossRate_e2 = EnergyLossRate(e2);

      deltaE1 = tbin * ElossRate_e1;
      deltaE2 = tbin * ElossRate_e2;
      /* Donor-cell advection */
      value = Grid[0][i] - quot * Grid[0][i] * ElossRate_e1 +
              quot * Grid[0][i + 1] * ElossRate_e2;
      if (!i) value += quot * Grid[0][i] * ElossRate_e1;

      /* these additional operations result in the superbee algorithm */
      //      value -=
      // 0.5*quot*(GetSuperBeeSlope(i,ebin,&Grid)*ElossRate_e1*(ebin-deltaE1)-GetSuperBeeSlope(i+1,ebin,&Grid)*ElossRate_e2*(ebin-deltaE2));
      /* these additional operations result in the minmod slope limiter
       * algorithm */
      value -= 0.5 * quot * (GetMinModSlope(i, ebin, &Grid) * ElossRate_e1 *
                                 (ebin - deltaE1) -
                             GetMinModSlope(i + 1, ebin, &Grid) * ElossRate_e2 *
                                 (ebin - deltaE2));

      value *= exp(-pow(e1 / eMax, CutOffFactor));

      /* Increase in particles in bin 'i' due to particle injection from the
       * source */
      value += tbin * SourceSpectrum(e1);

      /* Decrease in particles in bin 'i' due to particle escape.
       */
      if (escapeTime > 0.) value -= tbin * Grid[0][i] / escapeTime;

      /* Determine the largest filled energy bin (needed for the efficient
       * calculation of the next iterative time bin.
       */
      if (value > 0.) largestFilledBin = i;
      /* set the particle number (value) in at bin 'i' in Grid[i] */
      Grid[1][i] = value;
      particleCount += value * ebin;
    }
    /* put a "0" as the last element of this row in order to avoid edge effects.
     */
    Grid[1][Esize] = 0.;
    /* set the just now calculated spectrum as base spectrum for the next step
     * This way, only 2 vectors are needed for the calculation of the spectrum
     */
    double Econt = 0.;
    for (unsigned int ii = 0; ii < Grid[1].size(); ii++) Econt += Grid[0][ii];
    if (fabs(Econt) > 1.) Grid[0] = Grid[1];
  }

  /* Fill the final lookup, holding the time evolved spectrum at time = Age.
   * Also, forego edge bins in order to avoid artefacts.
   */
  Clear2DVector(ParticleSpectrum);
  for (unsigned int j = 1; j < EnergyAxis.size() - 1; j++) {
    e1 = pow(10., EnergyAxis[j]);
    double val = Grid[0][j];
    if(std::isnan(val) || std::isinf(val) || !val)
      continue;
    ParticleSpectrum.push_back(vector<double>());
    ParticleSpectrum[ParticleSpectrum.size() - 1].push_back(e1);
    ParticleSpectrum[ParticleSpectrum.size() - 1].push_back(val);
  }
  /* Important for wrapper function 'ComputeGridInTimeInterval': remove the last
   * vector (Grid[1]) so that Grid has the same shape as in the beginning of
   * this function.
   */
  Grid.pop_back();

  /* for the format of the info writeout */
  if (QUIETMODE == false) {
    cout << endl;
    cout << "    -> DONE!   " << endl;
    cout << endl;
  }

  return;
}

/** slope for MinMod slope limiter method */
double Particles::GetMinModSlope(int i, double deltaX,
                                 vector<vector<double> > *Grid) {
  double a = ((*Grid)[0][i] - (*Grid)[0][i + 1]) / deltaX;
  double b = ((*Grid)[0][i - 1] - (*Grid)[0][i]) / deltaX;
  double sigma = MinMod(a, b);
  return sigma;
}

/** slope for superbee slope limiter method */
double Particles::GetSuperBeeSlope(int i, double deltaX,
                                   vector<vector<double> > *Grid) {
  double a = ((*Grid)[0][i - 1] - (*Grid)[0][i]) / deltaX;
  double b = ((*Grid)[0][i] - (*Grid)[0][i + 1]) / deltaX;

  double sigma1 = MinMod(a, 2. * b);
  double sigma2 = MinMod(2. * a, b);

  double sigma = MaxMod(sigma1, sigma2);
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

/** wrapper function to calculate the grid only in a specified time interval dT
 * = T2-T2
 * This is especially useful for the creation of time-series of spectra and
 * necessary
 * as the 'ComputeGrid' function only stores two spectra (the actual one and
 * that in the
 * time step before) due to memory storage reasons.
 */
void Particles::ComputeGridInTimeInterval(double T1, double T2) {
  if (T1 <= Tmin) {
    cout << "Particles::ComputeGridInTimeInterval T1<internal min time (" << T1
         << "<" << Tmin << "). set it artificially to " << Tmin << "yrs."
         << endl;
    T1 = Tmin * 1.001;
  }

  Clear2DVector(ParticleSpectrum);
  ComputeGrid(grid, energyAxis, T1, T2, timeAxis, yr_to_sec * (T2 - T1) / 100.);
  return;
}

void Particles::CalcSpecSemiAnalyticNoELoss() {
  if (Tmin >= Age) {
    cout << "CalcSpecSemiAnalyticNoELoss: Tmin is larger/equal "
            "than source age... Exiting" << endl;
    return;
  }
  Clear2DVector(ParticleSpectrum);
  double totallum = 0.;
  if(LumConstant)
    totallum = LumConstant*Age;
  else if(LumVector.size()) {
    gsl_interp_accel_reset(accLum);
    if(gsl_spline_eval_integ_e(LumLookup, Tmin, Age, accLum, &totallum))
      totallum = 0.;
  }
  else return;
  LumConstant = totallum*yr_to_sec;
  SetMembers(Age);

  double logstep = (log10(eMax) - log10(Emin)) / ebins;
  for (double e = Emin; e < eMax; e = pow(10., log10(e) + logstep)) {
    double val = SourceSpectrum(e);
    if(std::isnan(val) || std::isinf(val) || !val)
      continue;
    ParticleSpectrum.push_back(vector<double>());
    ParticleSpectrum[ParticleSpectrum.size() - 1].push_back(e);
    ParticleSpectrum[ParticleSpectrum.size() - 1].push_back(val);
  }
  return;
}

void Particles::CalcSpecSemiAnalyticConstELoss() {
  if (Tmin >= Age) {
    cout << "Particles::CalcSpecSemiAnalyticConstELoss: Tmin is larger/equal "
            "than source age... Exiting" << endl;
    return;
  }
  fPointer IntFunc = NULL;
  Clear2DVector(ParticleSpectrum);
  double logstep = (log10(eMax) - log10(Emin)) / ebins;

  /* info writeout. Disable it by using 'ToggleQuietMode()' */
  if (!Type && QUIETMODE == false)
    cout << "** Evolving Electron Spectrum:" << endl;
  else if (Type == 1 && QUIETMODE == false)
    cout << "** Evolving Proton Spectrum:" << endl;
  // steady state solution
  double maxCoolingTime = -100.;
  for (double E = Emin; E < eMax; E *= 1.01) {
    double CoolingTime = E / EnergyLossRate(E);
    if (CoolingTime > maxCoolingTime) maxCoolingTime = CoolingTime;
  }
  maxCoolingTime /= yr_to_sec;
  if (Age > maxCoolingTime) {
    double dummy = 0.;
    IntFunc = &Particles::SourceSpectrumWrapper;
    int tt = 0;
    SetMembers(Age);
    for (double e = Emin; e < eMax; e = pow(10., log10(e) + logstep)) {
      if (QUIETMODE == false)
        cout << "    " << (int)(100. * tt / ebins) << "\% done\r" << std::flush;
      double val = Integrate(IntFunc, &dummy, e, eMax, integratorTolerance) /
                   EnergyLossRate(e);
      if(std::isnan(val) || std::isinf(val) || !val)
        continue;
      ParticleSpectrum.push_back(vector<double>());
      ParticleSpectrum[ParticleSpectrum.size() - 1].push_back(e);
      ParticleSpectrum[ParticleSpectrum.size() - 1].push_back(val);
      tt++;
    }
  }
  // time integration (constant energy losses)
  else {
    Tmin = vETrajectory[0][0];
    IntFunc = &Particles::SemiAnalyticConstELossIntegrand;
    int tt = 0;
    for (double e = Emin; e < eMax; e = pow(10., log10(e) + logstep)) {
      if (QUIETMODE == false)
        cout << "    " << (int)(100. * tt / ebins) << "\% done\r" << std::flush;
      double val = Integrate(IntFunc, &e, log10(Tmin), log10(Age), 5.e-3);
      SetMembers(Age);
      val /= EnergyLossRate(e);
      if(std::isnan(val) || std::isinf(val) || !val)
        continue;
      ParticleSpectrum.push_back(vector<double>());
      ParticleSpectrum[ParticleSpectrum.size() - 1].push_back(e);
      ParticleSpectrum[ParticleSpectrum.size() - 1].push_back(val);
      tt++;
    }
  }
  return;
}


double Particles::SemiAnalyticConstELossIntegrand(double T, void *par) {
  double tdash, E, Enow;
  Enow = *(double *)par;
  if (gsl_spline_eval_e(energyTrajectoryInverse, log10(Enow), accTrInv, &tdash))
    return 0.;
  if (T > tdash) return 0.;
  T = pow(10., T);
  if(Age - T < Tmin) return 0.;
  tdash = pow(10., tdash);
  if (gsl_spline_eval_e(energyTrajectory, log10(tdash - T), accTr, &E))
    return 0.;
  E = pow(10., E);
  if (E > eMax || E <= 0.) return 0.;
  SetMembers(Age - T);
  if (!SACELI_Told) {
    SACELI_Told = T;
    return 0.;
  }
  double dT = (T - SACELI_Told) * yr_to_sec;
  double dlogT = log10(T) - log10(SACELI_Told);
  SACELI_Told = T;
  return EnergyLossRate(E) * SourceSpectrum(E) * dT / dlogT;
}

void Particles::SetType(string type) {
  if (!type.compare("electrons"))
    Type = 0;
  else if (!type.compare("protons"))
    Type = 1;
  else
    cout << "Particles::SetType: What the f***! Specify supported particle "
            "species! " << endl;
  return;
}

double Particles::SourceSpectrumWrapper(double E, void *par) {
  return SourceSpectrum(E);
}

/**
 * Determine the common time boundaries of the provided lookups for the
 * evolution of parameters (Source Luminosity,B-Field, Ambient density etc.).
 * If everything is set to constant values, i.e. in a stationary scenario,
 * these boundaries will be set to {TminInternal,TmaxInternal} = {1.e-10,1.e100}yrs.
 */
void Particles::DetermineLookupTimeBoundaries() {

  // determine starting time and energy
  double lumtmin, emaxtmin, ntmin, btmin, rtmin, vtmin, esctmin, lumtmax,
      emaxtmax, ntmax, btmax, rtmax, vtmax, esctmax;
  lumtmin = emaxtmin = ntmin = btmin = rtmin = vtmin = esctmin = 1.e-10;
  lumtmax = emaxtmax = ntmax = btmax = rtmax = vtmax = esctmax = 1.e100;

  if (LumVector.size()) lumtmin = LumVector[0][0];
  if (eMaxVector.size()) emaxtmin = eMaxVector[0][0];
  if (NVector.size()) ntmin = NVector[0][0];
  if (BVector.size()) btmin = BVector[0][0];
  if (RVector.size()) rtmin = RVector[0][0];
  if (VVector.size()) vtmin = VVector[0][0];
  if (EscapeVector.size()) esctmin = EscapeVector[0][0];

  if (LumVector.size()) lumtmax = LumVector[LumVector.size() - 1][0];
  if (eMaxVector.size()) emaxtmax = eMaxVector[eMaxVector.size() - 1][0];
  if (NVector.size()) ntmax = NVector[NVector.size() - 1][0];
  if (BVector.size()) btmax = BVector[BVector.size() - 1][0];
  if (RVector.size()) rtmax = RVector[RVector.size() - 1][0];
  if (VVector.size()) vtmax = VVector[VVector.size() - 1][0];
  if (EscapeVector.size()) esctmin = EscapeVector[EscapeVector.size() - 1][0];

  double T0, T1, T2, T3, T4, T;

  // first determine highest common lower boundary time of the lookups
  if (!emaxtmin && !ntmin && !btmin && !rtmin && !vtmin)
    T = TminInternal;
  else {
    (rtmin > vtmin) ? T0 = rtmin : T0 = vtmin;
    (ntmin > T0) ? T1 = ntmin : T1 = T0;
    (btmin > T1) ? T2 = btmin : T2 = T1;
    (esctmin > T2) ? T3 = esctmin : T3 = T2;
    (lumtmin > T3) ? T4 = lumtmin : T4 = T3;
    (emaxtmin > T4) ? T = emaxtmin : T = T4;
  }
  TminInternal = T;

  // now do the similar thing to get the lowest common higher lookup boundary
  if (!emaxtmin && !ntmin && !btmin && !rtmin && !vtmin)
    T = TmaxInternal;
  else {
    (rtmax < vtmax) ? T0 = rtmax : T0 = vtmax;
    (ntmax < T0) ? T1 = ntmax : T1 = T0;
    (btmax < T1) ? T2 = btmax : T2 = T1;
    (esctmax < T2) ? T3 = esctmax : T3 = T2;
    (lumtmax < T3) ? T4 = lumtmax : T4 = T3;
    (emaxtmax < T4) ? T = emaxtmax : T = T4;
  }
  TmaxInternal = T;

  return;
}

void Particles::CalculateEnergyTrajectory(double TExt) {
  if (TminInternal < 0.) {
    cout << "Particles::CalculateEnergyTrajectory: Calculate internal Tmin "
            "first by running DetermineLookupStartingTime(). Exiting." << endl;
    return;
  }
  Clear2DVector(vETrajectory);

  double T = TminInternal;
  double E, Edot, dt;
  vETrajectory.clear();
  if (TExt) (TExt > T) ? T = TExt : 1;
  T *= 1.1;
  SetMembers(T);
  E = eMax;

  while (E >= Emin) {
    vETrajectory.push_back(vector<double>());
    vETrajectory[vETrajectory.size() - 1].push_back(T);
    vETrajectory[vETrajectory.size() - 1].push_back(E);
    Edot = EnergyLossRate(E);
    dt = 2.e-2 * E / Edot;
    E -= dt * Edot;
    T += dt / yr_to_sec;
    if(T>TmaxInternal) break;
    SetMembers(T);
  }
  int size = (int)vETrajectory.size();
  double x1[size], y1[size], x2[size], y2[size];
  for (int i = 0; i < size; i++) {
    double xVal = log10(vETrajectory[i][0]);
    double eVal = log10(vETrajectory[i][1]);
    if (std::isinf(xVal) || std::isinf(eVal)) continue;
    if (std::isnan(xVal) || std::isnan(eVal)) continue;
    if (std::isnan(xVal) || std::isnan(eVal)) continue;
    x1[i] = xVal;
    y1[i] = eVal;
    x2[size - 1 - i] = eVal;
    y2[size - 1 - i] = xVal;
  }

  energyTrajectory = gsl_spline_alloc(gsl_interp_linear, size);
  gsl_spline_init(energyTrajectory, x1, y1, size);
  energyTrajectoryInverse = gsl_spline_alloc(gsl_interp_linear, size);
  gsl_spline_init(energyTrajectoryInverse, x2, y2, size);
  return;
}

/**
 * Return a particle SED dN/dE vs E (erg vs TeV)
 */
vector<vector<double> > Particles::GetParticleSED() {
  vector<vector<double> > v;
  for (unsigned int i = 0; i < ParticleSpectrum.size(); i++) {
    double E = ParticleSpectrum[i][0];
    double ETeV = E / TeV_to_erg;
    double N = ParticleSpectrum[i][1];
    if (!N) continue;
    v.push_back(vector<double>());
    v[v.size() - 1].push_back(ETeV);
    v[v.size() - 1].push_back(E * E * N);
  }
  return v;
}

void Particles::Clear2DVector(vector< vector<double> > &v) {
  for (unsigned int i = 0; i < v.size(); i++) v[i].clear();
  v.clear();
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
/**
 * Integration function using the GSL QAG functionality
 *
 */
double Particles::Integrate(fPointer f, double *x, double emin, double emax,
                            double tolerance) {
  double integral, error;
  /* no comment */
  auto ptr = [=](double xx)->double {
    return (this->*f)(xx, (void *)x);
  };
  GSLfuncPart<decltype(ptr)> Fp(ptr);
  gsl_function F = *static_cast<gsl_function *>(&Fp);
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(gslmemory);
  if (gsl_integration_qag(&F, emin, emax, 0, tolerance, gslmemory, 1, w, &integral,
                          &error)) {
    gsl_integration_workspace_free(w);
    return 0.;
  }
  gsl_integration_workspace_free(w);
  return integral;
}
