#include "Astro.h"

Astro::Astro() {

  fUtils = new Utils();

  /* Standard assumption of 4 spiral arms */
  ArmsVector.push_back(1);
  ArmsVector.push_back(2);
  ArmsVector.push_back(3);
  ArmsVector.push_back(4);

  accArm1 = gsl_interp_accel_alloc();
  accArm2 = gsl_interp_accel_alloc();
  accArm3 = gsl_interp_accel_alloc();
  accArm4 = gsl_interp_accel_alloc();

  TaylorCordesArm1 = gsl_spline_alloc(gsl_interp_linear, 7);
  TaylorCordesArm2 = gsl_spline_alloc(gsl_interp_linear, 7);
  TaylorCordesArm3 = gsl_spline_alloc(gsl_interp_linear, 7);
  TaylorCordesArm4 = gsl_spline_alloc(gsl_interp_linear, 7);

  TaylorCordesArm1Inv = gsl_spline_alloc(gsl_interp_linear, 7);
  TaylorCordesArm2Inv = gsl_spline_alloc(gsl_interp_linear, 7);
  TaylorCordesArm3Inv = gsl_spline_alloc(gsl_interp_linear, 7);
  TaylorCordesArm4Inv = gsl_spline_alloc(gsl_interp_linear, 7);

  /*  Spline implementation of the Taylor & Cordes arms */
  x_a1[0]=3.53;   y_a1[0]=164.-90.;
  x_a1[1]=3.76;   y_a1[1]=200.-90.;
  x_a1[2]=4.44;   y_a1[2]=240.-90.;
  x_a1[3]=5.24;   y_a1[3]=280.-90.;
  x_a1[4]=5.36;   y_a1[4]=290.-90.;
  x_a1[5]=5.81;   y_a1[5]=315.-90.;
  x_a1[6]=5.8101; y_a1[6]=330.-90.;

  gsl_spline_init(TaylorCordesArm1, x_a1, y_a1, 7);
  gsl_spline_init(TaylorCordesArm1Inv, y_a1, x_a1, 7);

  x_a2[0]=3.76;  y_a2[0]=63.-90.;
  x_a2[1]=4.56;  y_a2[1]=120.-90.;
  x_a2[2]=4.79;  y_a2[2]=160.-90.;
  x_a2[3]=5.70;  y_a2[3]=200.-90.;
  x_a2[4]=6.49;  y_a2[4]=220.-90.;
  x_a2[5]=7.29;  y_a2[5]=250.-90.;
  x_a2[6]=8.20;  y_a2[6]=288.-90.;

  gsl_spline_init(TaylorCordesArm2, x_a2, y_a2, 7);
  gsl_spline_init(TaylorCordesArm2Inv, y_a2, x_a2, 7);

  x_a3[0]=4.90;  y_a3[0]=52.-90.;
  x_a3[1]=6.27;  y_a3[1]=120.-90.;
  x_a3[2]=6.49;  y_a3[2]=170.-90.;
  x_a3[3]=6.95;  y_a3[3]=180.-90.;
  x_a3[4]=8.20;  y_a3[4]=200.-90.;
  x_a3[5]=8.89;  y_a3[5]=220.-90.;
  x_a3[6]=9.57;  y_a3[6]=252.-90.;

  gsl_spline_init(TaylorCordesArm3, x_a3, y_a3, 7);
  gsl_spline_init(TaylorCordesArm3Inv, y_a3, x_a3, 7);

  x_a4[0]=5.92;  y_a4[0]=20.-90.;
  x_a4[1]=7.06;  y_a4[1]=70.-90.;
  x_a4[2]=7.86;  y_a4[2]=100.-90.;
  x_a4[3]=9.68;  y_a4[3]=160.-90.;
  x_a4[4]=10.37;  y_a4[4]=180.-90.;
  x_a4[5]=11.39;  y_a4[5]=200.-90.;
  x_a4[6]=12.08;  y_a4[6]=223.-90.;

  gsl_spline_init(TaylorCordesArm4, x_a4, y_a4, 7);
  gsl_spline_init(TaylorCordesArm4Inv, y_a4, x_a4, 7);

  TaylorCordesArms.push_back(TaylorCordesArm1);
  TaylorCordesArms.push_back(TaylorCordesArm2);
  TaylorCordesArms.push_back(TaylorCordesArm3);
  TaylorCordesArms.push_back(TaylorCordesArm4);

  TaylorCordesArmsInv.push_back(TaylorCordesArm1Inv);
  TaylorCordesArmsInv.push_back(TaylorCordesArm2Inv);
  TaylorCordesArmsInv.push_back(TaylorCordesArm3Inv);
  TaylorCordesArmsInv.push_back(TaylorCordesArm4Inv);

  /* Vallee Spiral parameters */
  psi0V.push_back(12.);
  psi0V.push_back(12.);
  psi0V.push_back(12.);
  psi0V.push_back(12.);
  theta0V.push_back(0.);
  theta0V.push_back(-90.);
  theta0V.push_back(-180.);
  theta0V.push_back(-270.);
  r0V.push_back(2.72);
  r0V.push_back(2.72);
  r0V.push_back(2.72);
  r0V.push_back(2.72);

  barRadius = 3.6;
  barAngle = 14.;

  armWidth = 0.56;

  rMax = 16.;

  SPIRALARMMODEL = 0;
  SURFACEDENSITYMODEL = 0;
  CENTRALSTRUCTUREMODEL = 0;
  MSBUBBLEMODEL = 0;
  GAMERAALREADYONTHECONSOLE = false;
  WRAPPINGHAPPENED = false;
}

Astro::~Astro() {}


/*    ------        SNR PROGENITOR STUFF       ------    */

/* Dice an initial star mass following the Salpeter law */
double Astro::RandomSalpeterInitialMass() {
  double initialMass = fUtils->PowerLawRandom(2.35,8.,100.,1.)[0];
  return initialMass;
}

/**
 * calculate the radius created by the main sequence wind of a star.
 * Implemented at the moment are the classical Weaver model (1977) and
 * the empirical relation found by Chen et al. 2013
 */
double Astro::MainSequenceBubbleRadius(double mDotMS, double vMSWind, double timeMS, double n, double initialMass) {
  double bubbleRadius = 0.;
  if(!MSBUBBLEMODEL) {
    /* bubble radius
     * Astrophysical Journal, Part 1, vol. 218, Dec. 1, 1977, p. 377-395, Eq. 51
     */
    double L36 = 0.5*mDotMS*(mSol/yr_to_sec)*vMSWind*vMSWind/1.e36;
    bubbleRadius = pow(L36/n,0.2);
    bubbleRadius *= 27.*pow(timeMS/1.e6,0.6);
  }
  else if(MSBUBBLEMODEL==1) {
    /* Chen et al., ApJL, 769:L16 (5pp), 2013 May 20
     * I have my doubts about
     * this because quoted interclump pressures don't
     * fit to the average values I use in my model but
     * rather dense molecular clouds. (i.e. p_5 = 1)
     */
    bubbleRadius = -9.15402+1.21532*initialMass;
    double randomfactor = fUtils->GaussianRandom(0.2*bubbleRadius, 0., 1)[0];
    while(randomfactor<-bubbleRadius) {
      randomfactor = fUtils->GaussianRandom(0.2*bubbleRadius, 0., 1)[0];
    }
    bubbleRadius += randomfactor;
  }
  else std::cout<<"Astro::MainSequenceBubbleRadius: Specify valid model (at the moment either 'Weaver' or 'Chen'! Returning 0 value."<<std::endl;
  bubbleRadius *= pc_to_cm;
  return bubbleRadius;
}


/**
 * Density inside the main-sequence wind blown bubble.
 */
double Astro::CalculateMSBubbleDensity(double mDotMS, double vMSWind, double timeMS, double n, double bubbleRadius) {
  double bubbleDensity = 0.;
  if(!MSBUBBLEMODEL) {
    /* Castor, McCray, Weaver 1975
     *(Astrophysical Journal, vol. 200, Sept. 1, 1975, pt. 2, p. L107-L110),
     * Eq. 8
     */
    double L = (mDotMS/1.e-6)*pow(vMSWind/2.e8,2.);
    bubbleDensity = 0.01*pow(pow(L,6.)*pow(n,19.)*pow(timeMS/1.e6,-22.),1./35.);
  }

  else {
    /* Chevalier et al 1989
     * (ApJ, vol. 344, Sept. 1, 1989, p. 332-340)
     * Eq 2.5
     */
    double totalMassLoss = mDotMS*timeMS;
    bubbleDensity = 2.e-27*totalMassLoss*pow(bubbleRadius/pc_to_cm/20.,-3.)/m_p_g;
  }
  return bubbleDensity;
}

/**
 * Time on the main sequence in years as a function of initial star mass.
 * This is empirical from my own fit to data compiled by
 * Chen et al., ApJL, 769:L16 (5pp), 2013 May 20
 * his time estimates are based on Schaller1992
 */
double Astro::CalculateTimeOnMainSequence(double initialMass) {
  double timeMS = 1.29098e+06*exp(1.11877e+01*pow(initialMass,-0.6));
  return timeMS;
}

/**
 * Main sequence wind mass luminosity as a function of initial star mass.
 * This is empirical from my own fit to data compiled by
 * Chen et al., ApJL, 769:L16 (5pp), 2013 May 20
 */
double Astro::CalculateMSMassLuminosity(double initialMass) {
  double mDotMS = 1.e-29*pow(initialMass,20.);
  mDotMS *= pow(1.+pow(initialMass/11.,1./0.1),(3.3-20.)*0.1);
  return mDotMS;
}

/**
 * Main sequence wind speed as a function of initial star mass.
 * This is empirical from my own fit to data compiled by
 * Chen et al., ApJL, 769:L16 (5pp), 2013 May 20
 */
double Astro::CalculateMSWindSpeed(double initialMass) {
  double vMSWind = 1.e7 + 6.0e6*initialMass;
  return vMSWind;
}

/**
 * Radius of the red giant wind zone (Chevalier2004) //TODO find paper again!
 */
double Astro::CalculateRGWRadius(double pBubble, double mDotRGW, double vRGWind) {
  if(!pBubble) {
    /* Typical bubble pressure (Chevalier) with 30% Gaussian scatter */
    pBubble = 1.36e-16*1.e4;
    double randomfactor = fUtils->GaussianRandom(0.3*pBubble, 0., 1)[0];
    while(randomfactor<=-pBubble) {
      randomfactor = fUtils->GaussianRandom(0.3*pBubble, 0., 1)[0];
    }
    pBubble += randomfactor;
  }
  double RGWRadius = sqrt(0.881*mDotRGW*(mSol/yr_to_sec)*vRGWind/(4.*pi*pBubble));
  return RGWRadius;
}


vector< vector<double> > Astro::CreateDensityProfile(double RGWRadius, double mDotRGWind, double vRGWind, double MSBubbleRadius, double MSBubbleDensity, double ISMDensity, double mEj, double rSWRGW, double rSWMS) {
  vector< vector<double> > dProfile;
  /* calculate the mass that has been swept up by either the RGW or the wind bubble material moving backwards into the RGW */
  double mSWRGW;
  double rTransRGWMS = sqrt(mDotRGWind/(4.*pi*MSBubbleDensity*vRGWind));
  if(RGWRadius>rTransRGWMS) mSWRGW=mDotRGWind*(RGWRadius-rTransRGWMS)/vRGWind/m_p_g;
  else mSWRGW=4.*pi*MSBubbleDensity*(pow(rTransRGWMS,3.)-pow(RGWRadius,3.));
  /* calculate the mass that has been swept by the main sequence wind */
  double mSWMS = 4.*pi*ISMDensity*pow(MSBubbleRadius,3.);
  /* calcutate densities in the shells */
  double nSWRGW = 3.*mSWRGW/(4.*pi*(pow(RGWRadius+rSWRGW,3.)-pow(RGWRadius,3.)));
  double nSWMS = 3.*mSWMS/(4.*pi*(pow(MSBubbleRadius+rSWMS,3.)-pow(MSBubbleRadius,3.)));

  double rMin=1.e-5*pc_to_cm;
  double rMax=100.*pc_to_cm;
  double stepSizeNormal=log(rMax/rMin)/1000.;
  double stepSizeFine=stepSizeNormal/10.;
  double stepSize=stepSizeNormal;
  double n;
  for(double r=rMin;r<rMax;r=pow(10.,log10(r)+stepSize)) {
    /* adapt step size in regions of transition */
    if(r<0.9*RGWRadius) stepSize=stepSizeNormal;
    else if(r>=0.9*RGWRadius && r<1.1*(RGWRadius+rSWRGW)) stepSize=stepSizeFine;
    else if(r>=1.1*(RGWRadius+rSWRGW) && r<0.9*MSBubbleRadius) stepSize=stepSizeNormal;
    else if(r>=0.9*MSBubbleRadius && r<1.1*(MSBubbleRadius+rSWMS)) stepSize=stepSizeFine;
    else stepSize=stepSizeNormal;

    /* calculate density */
    if(r<MSBubbleRadius){
      if(r<RGWRadius) n = mDotRGWind/(4.*pi*vRGWind*r*r)/m_p_g;
      else if(r>=RGWRadius && r<RGWRadius+rSWRGW) n = nSWRGW;
      else if(r>=RGWRadius+rSWRGW && r<MSBubbleRadius) n = MSBubbleDensity;
    }
    else if(r>=MSBubbleRadius && r<MSBubbleRadius+rSWMS) n = nSWMS;
    else n = ISMDensity;
    dProfile.push_back(vector<double>());
    dProfile[dProfile.size()-1].push_back(r);
    dProfile[dProfile.size()-1].push_back(n);
  }
  return dProfile;
}

/*    ------      B-Field stuff      ------ */

/**
 * 2D-Model of the large-scale galactic magnetic field structure
 * by Jaffe et al, 2010.
 * MNRAS, Volume 401, Issue 2, pp. 1013-1028
 * Input x and y coordinate, returns
 * -total B-Field at this position in the plane (B_tot)
 * -large-scale coherent field (B_coh)
 * -irregular field component along the large-scale component (B_ord)
 * -random field component (B_iso)
 * -Direction in the x-y plane of the large scale coherent field (regFieldDirection)
 * Note: it is only in the x-y plane at the moment, thus the z variable doesn't
 * do anything!
 */
void Astro::CalculateBField(double x, double y, double z, double &B_tot, double &B_coh, double &B_ord, double &B_iso, vector<double> &regFieldDirection, vector<double> &isoFieldDirection) {
  regFieldDirection.clear();
  isoFieldDirection.clear();
  double DistanceToClosestArm = 0.;
  int ClosestArm = 0;
  GetDistanceToNearestSpiralArm(x,y,DistanceToClosestArm,ClosestArm);
  double r = sqrt(pow(x,2.)+pow(y,2.));
  double R0 = 20.;
  double R1 = 3.;
  double R2 = 0.5;
  double Rarms = 15.;
  double RmolRing = 5.;
  double Rmin = 3.;
  double C0 = 1.;
  double an = 1.;
  double an1 = 1.26;
  double an2 = 1.04;
  double an3 = 0.61;
  double an4 = 1.65;
  double anMolRing  = 1.;
  switch (ClosestArm){
    case 4:
      an = an4;
      break;
    case 3:
      an = an3;
      break;
    case 2:
      an = an2;
      break;
    case 1:
      an = an1;
      break;
  }
  double B_0 = 1.e-6;
  double d0 = 0.3; //could be replaced by arm-width
  double B_rms = 2.1e-6;
  double ford = 1.9;
  regFieldDirection.push_back(0.);
  regFieldDirection.push_back(0.);
  regFieldDirection.push_back(0.);

  isoFieldDirection.push_back(0.);
  isoFieldDirection.push_back(0.);
  isoFieldDirection.push_back(0.);
  if(r<Rmin||r>R0) {
    B_coh = 0.;
    B_iso = 0.;
    B_ord = 0.;
    B_tot = 0.;

    return;

  }
  double rhoc = C0*exp(-DistanceToClosestArm*DistanceToClosestArm/(d0*d0))+1.;
  if(r>Rarms) rhoc = 1.;
  if(r>Rmin && r<RmolRing) {
    rhoc  =  C0*exp(-pow(r-RmolRing,2.)/(d0*d0))+1.;///< own description, since exact profile is not described in paper
    an = anMolRing;
  }
  /* coherent, large scale B-field along spiral arms */
  B_coh = B_0*(1.-exp(-r*r/(R2*R2)));
  B_coh *= exp(-r*r/(R0*R0))+exp(-r*r*r*r/(R1*R1*R1*R1));
  B_coh *= rhoc*an;
  /* small-scale irregular component of B-field */ //NOTE: This is only random in 2D! If 3D (incl. scatter perpendicular to the plane), the component along the arms would smaller)
  B_iso = rhoc*B_rms;
  /* irregular component along spiral arms */
  B_ord = ford*(rhoc-1.)*B_rms;
  /* dice relative orientation of irregular components to coherent one */
  double phi_iso = 2.*pi*fUtils->Random();
  /* Component along the coherent field lines */
  double B_sumx = B_coh+B_iso*cos(phi_iso)+B_ord*cos(phi_iso);
  /* Component perpendicular to the coherent field lines (only isotropic field) */
  double B_sumy = B_iso*sin(phi_iso);

  B_tot = sqrt(B_sumx*B_sumx+B_sumy*B_sumy);

  double e_B_coh[3];
  double theta = -(180./pi)*atan2(x,y)+90;
//  std::cout<<"(x,y)="<<x<<","<<y<<" "<<theta<<std::endl;
  double xarmclose,yarmclose,rarmclose,xarmfar,yarmfar,rarmfar,x0,y0,x1,y1,rtest;
  PositionOnSpiralArmAngular(theta, ClosestArm, xarmclose, yarmclose);
  rarmclose = sqrt(xarmclose*xarmclose+yarmclose*yarmclose);
  PositionOnSpiralArmAngular(theta+360., ClosestArm, xarmfar, yarmfar);
  rarmfar = sqrt(xarmfar*xarmfar+yarmfar*yarmfar);
  (fabs(r-rarmfar)>fabs(r-rarmclose))?(rtest=rarmclose):
    rtest=rarmfar;
  GalacticPositionXY(0.999*rtest,ClosestArm,x0,y0);
  GalacticPositionXY(1.001*rtest,ClosestArm,x1,y1);
//  std::cout<<r<<" -> ("<<x0<<","<<y0<<")"<<" ("<<x1<<","<<y1<<")"<<std::endl;
  double l = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
  e_B_coh[0] = (x1-x0)/l;
  e_B_coh[1] = (y1-y0)/l;
  e_B_coh[2] = 0.;
  if(r>Rarms) {
    double weight = (r-Rarms)/(R0-Rarms);
    double xrand = 1.-2.*fUtils->Random();
    double yrand = 1.-2.*fUtils->Random();
    double zrand = 1.-2.*fUtils->Random();

    double rrand = sqrt(xrand*xrand+yrand*yrand+zrand*zrand);
    xrand /= rrand;
    yrand /= rrand;
    zrand /= rrand;

    e_B_coh[0] += weight*xrand;
    e_B_coh[1] += weight*yrand;
    e_B_coh[2] += weight*zrand;

    double r_coh = sqrt(e_B_coh[0]*e_B_coh[0]+e_B_coh[1]*e_B_coh[1]+e_B_coh[2]*e_B_coh[2]);

    e_B_coh[0] /= r_coh;
    e_B_coh[1] /= r_coh;
    e_B_coh[2] /= r_coh;
  }
  if(r>RmolRing) {
    regFieldDirection[0] = e_B_coh[0];
    regFieldDirection[1] = e_B_coh[1];
    regFieldDirection[2] = e_B_coh[2];
  }
//  std::cout<<"  -> "<<regFieldDirection[0]<<","<<regFieldDirection[1]<<","<<regFieldDirection[2]<<std::endl;

  double e_B_iso[3];
  e_B_iso[0] = e_B_coh[0];
  e_B_iso[1] = e_B_coh[1];
  e_B_iso[2] = e_B_coh[2];
  RotateCoordinates(e_B_iso[0],e_B_iso[1],e_B_iso[2],0.,0.,phi_iso);
  isoFieldDirection[0] = e_B_iso[0];
  isoFieldDirection[1] = e_B_iso[1];
  isoFieldDirection[2] = e_B_iso[2];

  return;
}

/*  -----     AMBIENT DENSITY STUFF     ----- */

/**
 *
 */
double Astro::HIDensity(double x, double y, double z) {
  double R = sqrt(pow(x,2.)+pow(y,2.));
  double h1density = 0.;
  double norm = CalculateHINorm(R);
  double FWHM = GetHIFWHM(R)/1000.;
  h1density = 0.7*exp(-pow(z/(0.55*FWHM),2.));
  h1density += 0.19*exp(-pow(z/(1.38*FWHM),2.));
  h1density += 0.11*exp(-fabs(z)/(1.75*FWHM));
  h1density *= norm;
  return h1density;
}

/**
 *
 */
double Astro::H2Density(double x, double y, double z) {
  double R = sqrt(pow(x,2.)+pow(y,2.));
  double h2density = 0.;
  double norm = CalculateH2Norm(R);
  double sigma = GetH2FWHM(R)/(2.*sqrt(2.*log(2.)));
  sigma /= 1000.;
  h2density = norm*exp(-pow(z/(sqrt(2.)*sigma),2.));
  return h2density;
}

/**
 *
 */
double Astro::CalculateHINorm(double R) {
  double CD = GetHIColumnDensity(R);
  double FWHM = GetHIFWHM(R);
  double norm = 0.;
  double sigma1 = 0.55*FWHM/sqrt(2.);
  double sigma2 = 1.38*FWHM/sqrt(2.);
  double z0 = 1.75*FWHM;
  double k1=0.7;
  double k2=0.19;
  double k3=0.11;
  norm = CD/(3.08e18*(sqrt(2.*3.1415)*(k1*sigma1+k2*sigma2)+2.*k3*z0));
  return norm;
}

/**
 *
 */
double Astro::CalculateH2Norm(double R) {
  double CD = GetH2ColumnDensity(R);
  double FWHM = GetH2FWHM(R);
  double norm = 0.;
  double sigma = FWHM/(2.*sqrt(2.*log(2.)));
  if(sigma) norm = CD/(3.08e18*sigma*sqrt(2.*3.1415));
  else norm = 0.;
  return norm;
}

/**
 *
 */
double Astro::GetHIColumnDensity(double R) {
  double CD=0.;
  if(R<3.5) CD = 1.e20*6.2*exp(-pow((R-3.5)/0.7,2.)/2.);
  else if(R<13.65) CD = 1.e20*6.2;
  else CD = 1.e20*6.2*exp(-0.28*(R-13.65));
  return CD;
}

/**
 *
 */
double Astro::GetH2ColumnDensity(double R) {
  double CD=0.;
  CD = 1.e20*1.54269e-01*pow(R,8)*exp(-1.69222e+00*R);
  return CD;
}

/**
 *
 */
double Astro::GetHIFWHM(double R) {
  R = R*1000.;
  double FWHM=0.;
  if(R<3500.) FWHM = (230.-165.)*exp(-pow(R-3500.,2.)/(2.*pow(500.,2.)))+165.;
  else if(R<=8500.) FWHM = 230.;
  else FWHM = ((2000.-230.)/(20000.-8500.))*(R-8500.)+230.;
  return FWHM;
}

/**
 *
 */
double Astro::GetH2FWHM(double R) {
  R *= 1000.;
  double FWHM = 2.*sqrt(2.*log(2.))*pow(R/7.81,0.58);
  return FWHM;
}

/**
 *
 */
double Astro::ModulateGasDensityWithSpirals(double n, double x, double y, double z) {
  double r = sqrt(x*x+y*y);
  if(r<=barRadius) return n;
  /* normalisation quantity that conserves the gas mass before and after modulation. Depends on
   * spiral arm width; this is for the gaussian arm width of 560pc (Russeil, D. 2005, A&A, 397, 133-146)
   */
  double n1 = 0.315;
  double DistanceToClosestArm = 0.;
  int ClosestArm = 0;
  GetDistanceToNearestSpiralArm(x,y,DistanceToClosestArm,ClosestArm);
  double rhoc = n1*(3.*exp(-DistanceToClosestArm*DistanceToClosestArm/(2.*armWidth*armWidth))+1.);
  return rhoc*n;
}

// /**
//  * calculate gas column densities from the on-point given by x,y,z
//  * in direction l,b out to radius r
//  */
// double Astro::CalculateGasColumnDensity(vector<double> xyzReference,
//                                         vector<double> GLGB,
//                                         string gascomponent,
//                                         double modulate,
//                                         double range,
//                                         double steps) {
//
//   double GasSpecies = -1.;
//   if(!gascomponent.compare("HI")) GasSpecies = 0.;
//   else if(!gascomponent.compare("H2")) GasSpecies = 1.;
//   else {
//     std::cout<<"Astro::CalculateGasColumnDensity: Specify supported Gas component (currently 'HI' and 'H2')! returning 0 value."<<std::endl;
//     return 0.;
//   }
//
//   TF1 *cdintegrand = new TF1("cd",this,&Astro::nRadial,0.,range,7,"Astro","cdintegrand");
//   cdintegrand->SetParameter(0,GLGB[0]);
//   cdintegrand->SetParameter(1,GLGB[1]);
//   cdintegrand->SetParameter(2,xyzReference[0]);
//   cdintegrand->SetParameter(3,xyzReference[1]);
//   cdintegrand->SetParameter(4,xyzReference[2]);
//   cdintegrand->SetParameter(5,GasSpecies);
//   cdintegrand->SetParameter(6,modulate);
//   double cd = kpc_to_cm*Integrate(cdintegrand,0.,range,steps,false);
//
// //  std::cout<<GasSpecies<<" "<<modulate<<" "<<SPIRALARMMODEL<<" "<<cd<<std::endl;
//   return cd;
// }


/**
 *
 */
double Astro::nRadial(double *x, double *pars) {
  double r = x[0];
  double l = pars[0];
  double b = pars[1];
  double xRef = pars[2];
  double yRef = pars[3];
  double zRef = pars[4];
  int gasspecies = pars[5];
  int modulatewitharms = (int)pars[6];

  double xActual,yActual,zActual;
  GetCartesian(r,l,b,xRef,yRef,zRef,xActual,yActual,zActual);

  double n=0.;
  if(!gasspecies) n = HIDensity(xActual,yActual,zActual);
  else if(gasspecies==1) n = H2Density(xActual,yActual,zActual);
  else {
    std::cout<<"Astro::nRadial: Specify implemented gas species! returning 0.!"<<std::endl;
    return 0.;
  }
  if(modulatewitharms==1) n = ModulateGasDensityWithSpirals(n,xActual,yActual,zActual);
  return n;
}

/*    ------            Galactic Structure stuff              ------    */

/**
 * Galactic coordinates (R,GL,GB) ->Cartesian
 */
void Astro::GetCartesian(double r, double l, double b, double xref, double yref, double zref, double &x, double &y, double &z) {

  double dl, db;
  dl = db = 0.;
  dl = atan2(yref,xref);
  db = atan(zref/sqrt(xref*xref+yref*yref));
  double phi = (pi/180.)*l+dl;
  double theta = pi/2. - (pi/180.)*b+db;
  x = xref - r*sin(theta)*cos(phi);
  y = yref - r*sin(theta)*sin(phi);
  z = zref + r*cos(theta);
  return;
}

/**
 * Cartesian coordinates -> Galactic Coordinates(R,GL,GB)
 */
void Astro::GetGalactic(double x, double y, double z, double xref, double yref, double zref, double &l, double &b) {
  double dx  = x-xref;
  double dy  = y-yref;
  double dz =  z-zref;

  double dr = sqrt(dx*dx+dy*dy+dz*dz);
  double rrref = sqrt(xref*xref+yref*yref);
  double drr = sqrt(dx*dx+dy*dy);
  double prod = -dx*xref - dy*yref;
  l = 180.*acos(prod/(drr*rrref))/pi;
  if(dx<=0.) l*=-1.;
  b = 180.*asin(dz/dr)/pi;
  return;
}

/**
 * rotate coordinate around (0,0,0) by (phi,theta,psi)
 */
void Astro::RotateCoordinates(double &x, double &y, double &z, double phi, double theta, double psi) {

  double xyz[3];
  xyz[0] = x;
  xyz[1] = y;
  xyz[2] = z;

  phi *= pi/180.;// rotation around x
  theta *= pi/180.;// rotation around y
  psi *= pi/180.;// rotation around z

  double rx[3][3];
  double ry[3][3];
  double rz[3][3];

  /* rotation matrices */
  rx[0][0] = 1.; rx[0][1] =   0.    ; rx[0][2] =  0.      ;
  rx[1][0] = 0.; rx[1][1] = cos(phi); rx[1][2] = -sin(phi);
  rx[2][0] = 0.; rx[2][1] = sin(phi); rx[2][2] =  cos(phi);

  ry[0][0] =  cos(theta); ry[0][1] =   0.    ; ry[0][2] = sin(theta);
  ry[1][0] =     0.     ; ry[1][1] =   1.    ; ry[1][2] =     0.    ;
  ry[2][0] = -sin(theta); ry[2][1] =   0.    ; ry[2][2] = cos(theta);

  rz[0][0] = cos(psi); rz[0][1] = -sin(psi) ; rz[0][2] = 0.;
  rz[1][0] = sin(psi); rz[1][1] =  cos(psi) ; rz[1][2] = 0.;
  rz[2][0] =    0.   ; rz[2][1] =    0.     ; rz[2][2] = 1.;

  //TODO: this can be compactified!
  if(phi) {
    double new_coord[3] = {0.,0.,0.};
    for(int j=0;j<3;j++) {
      for(int k=0;k<3;k++) {
         new_coord[j] += rx[j][k]*xyz[k];
      }
    }
    for(int l=0;l<3;l++) xyz[l] = new_coord[l];
  }
  if(theta) {
    double new_coord[3] = {0.,0.,0.};
    for(int j=0;j<3;j++) {
      for(int k=0;k<3;k++) {
        new_coord[j] += ry[j][k]*xyz[k];
      }
    }
    for(int l=0;l<3;l++) xyz[l] = new_coord[l];
  }
  if(psi) {
    double new_coord[3] = {0.,0.,0.};
    for(int j=0;j<3;j++) {
      for(int k=0;k<3;k++) {
        new_coord[j] += rz[j][k]*xyz[k];
      }
    }
    for(int l=0;l<3;l++) xyz[l] = new_coord[l];
  }

  x = xyz[0];
  y = xyz[1];
  z = xyz[2];
  return;
}

/**
 * Switch off an individual arm in Astro spiral model
 */
void Astro::DisableArm(int arm) {
  if(!QUIETMODE) std::cout<<"Disabling Arm no. "<<arm<<"."<<std::endl;
  bool Set = false;
  for(unsigned int i=0;i<ArmsVector.size();i++) {
    if(ArmsVector[i]==arm) {
      ArmsVector[i]=0;
      Set = true;
    }
  }
  if(Set==false) std::cout<<"Arm "<<1<<" already disabled."<<std::endl;
  vector<int> tempvec;
  for(unsigned int i=0;i<ArmsVector.size();i++) {
    if(ArmsVector[i]) tempvec.push_back(ArmsVector[i]);
  }
  ArmsVector = tempvec;
  if(!QUIETMODE) {
    for(unsigned int i=0;i<ArmsVector.size();i++) std::cout<<"Arm["<<i<<"]="<<ArmsVector[i]<<std::endl;
  }
  return;
}

void Astro::GetDistanceToNearestSpiralArm(double x, double y, double &DistanceToClosestArm, int &ClosestArm) {
  double theta,drdr,xx,yy,theta0,psi0fac,r0;
  theta = drdr = xx = yy = theta0 = psi0fac = r0 = 0.;
  double ddr=1.e12;
  int arm=0;

  if(SPIRALARMMODEL==0 && ArmsVector.size()!=4) std::cout<<"GetDistanceToNearestSpiralArm: Number of arms is wrong (n="<<ArmsVector.size()<<"; Vallee model has 4 arms). Returning 0 value."<<std::endl;
  if(SPIRALARMMODEL==1 && ArmsVector.size()!=4) std::cout<<"GetDistanceToNearestSpiralArm: Number of arms is wrong (n="<<ArmsVector.size()<<"; Taylor&Cordes model has 4 arms). Returning 0 value."<<std::endl;

  for(unsigned int i=0;i<ArmsVector.size();i++) {
    if(SPIRALARMMODEL==0) {
      theta0 = (pi/180.)*theta0V[ArmsVector[i]-1];
      psi0fac = pow(tan((pi/180.)*psi0V[ArmsVector[i]-1]),-1.);
      r0 = r0V[ArmsVector[i]-1];
      for(double r=0.00001;r<rMax;r+= .01){
        theta = theta0+psi0fac*(log(r/r0));
        xx = r*cos(theta)-x;
        yy = r*sin(theta)-y;
        drdr = xx*xx + yy*yy;
        if(drdr<ddr) {
          ddr = drdr;
          arm = i;
        }
      }
    }
    else if(SPIRALARMMODEL==1) {
      for(double r=0.00001;r<rMax;r+= .03){
        theta = EvalTaylorCordesArmTheta(r,ArmsVector[i]);
        xx = r*cos(theta)-x;
        yy = r*sin(theta)-y;
        drdr = xx*xx + yy*yy;
        if(drdr<ddr) {
          ddr = drdr;
          arm = i;
        }
      }
    }
    else {
      std::cout<<"GetDistanceToNearestSpiralArm: WTF? Specify proper spiral arm model!"<<std::endl;
      return;
    }
  }
  DistanceToClosestArm = sqrt(ddr);
  ClosestArm = ArmsVector[arm];
//  std::cout<<"Distance = "<<ddr<<" Arm = "<<ClosestArm<<std::endl;
  return;
}

double Astro::GetDistanceToGivenSpiralArm(double x, double y, int arm) {
  double theta,drdr,xx,yy,theta0,psi0fac,r0;
  theta = drdr = xx = yy = theta0 = psi0fac = r0 = 0.;
  double ddr=1.e12;

  if(SPIRALARMMODEL==0 && ArmsVector.size()!=4) std::cout<<"GetDistanceToNearestSpiralArm: Number of arms is wrong (n="<<ArmsVector.size()<<"; Vallee model has 4 arms). Returning 0 value."<<std::endl;
  if(SPIRALARMMODEL==1 && ArmsVector.size()!=4) std::cout<<"GetDistanceToNearestSpiralArm: Number of arms is wrong (n="<<ArmsVector.size()<<"; Taylor&Cordes model has 4 arms). Returning 0 value."<<std::endl;
    if(SPIRALARMMODEL==0) {
      theta0 = (pi/180.)*theta0V[arm-1];
      psi0fac = pow(tan((pi/180.)*psi0V[arm-1]),-1.);
      r0 = r0V[arm-1];
      for(double r=0.00001;r<rMax;r+= .05){
        theta = theta0+psi0fac*(log(r/r0));
        xx = r*cos(theta)-x;
        yy = r*sin(theta)-y;
        drdr = xx*xx + yy*yy;
        if(drdr<ddr) {
          ddr = drdr;
        }
      }
    }
    else if(SPIRALARMMODEL==1) {
      for(double r=0.00001;r<rMax;r+= .03){
        theta = EvalTaylorCordesArmTheta(r,arm);
        xx = r*cos(theta)-x;
        yy = r*sin(theta)-y;
        drdr = xx*xx + yy*yy;
        if(drdr<ddr) {
          ddr = drdr;
        }
      }
    }
    else {
      std::cout<<"GetDistanceToGivenSpiralArm: WTF? Specify proper spiral arm model! Returning 0."<<std::endl;
      return 0.;
    }

  return sqrt(ddr);
}

double Astro::EvalTaylorCordesArmTheta(double r, int arm) {
  double theta = 0.;
  if(!arm) {
    cout << "Astro::EvalTaylorCordesArmTheta: Requested Arm 0! (Arm Index starts a 1, not 0.) "
            "Returning 0. " << endl;
    return 0.;
  }
  gsl_interp_accel *a = NULL;
  double rmin,rmax;
  switch (arm){
    case 1:
      a = accArm1;
      rmin = x_a1[0]; rmax = x_a1[6];
      break;
    case 2:
      a = accArm2;
      rmin = x_a2[0]; rmax = x_a2[6];
      break;
    case 3:
      a = accArm3;
      rmin = x_a3[0]; rmax = x_a3[6];
      break;
    case 4:
      a = accArm4;
      rmin = x_a4[0]; rmax = x_a4[6];
      break;
  }
  if(r<rmin || r>rmax) return 0.;
  if (gsl_spline_eval_e(TaylorCordesArms[arm - 1], r, a, &theta)) {
    cout << "Astro::EvalTaylorCordesArmTheta: Interpolation in arm "
         << arm << " failed. Exiting with 0. return." << endl;
    return 0.;
  }
  if (isnan(theta) || isinf(theta)) return 0.;
  theta *= (pi/180.);
  return theta;
}

double Astro::EvalTaylorCordesArmGalactocentricRadius(double theta, int arm) {
  double r = 0.;
  if(!arm) {
    cout << "Astro::EvalTaylorCordesArmGalactocentricRadius: Requested Arm 0! "
             "(Arm Index starts a 1, not 0.) Returning 0. " << endl;
    return 0.;
  }
  gsl_interp_accel *a = NULL;
  double thetamin,thetamax;
  switch (arm){
    case 1:
      a = accArm1Inv;
      thetamin = y_a1[0]; thetamax = y_a1[6];
      break;
    case 2:
      a = accArm2Inv;
      thetamin = y_a2[0]; thetamax = y_a2[6];
      break;
    case 3:
      a = accArm3Inv;
      thetamin = y_a3[0]; thetamax = y_a3[6];
      break;
    case 4:
      a = accArm4Inv;
      thetamin = y_a4[0]; thetamax = y_a4[6];
      break;
  }
  if(theta<thetamin || theta>thetamax) return 0.;
  if (gsl_spline_eval_e(TaylorCordesArmsInv[arm - 1], theta, a, &r)) {
    cout << "Astro::EvalTaylorCordesArmGalactocentricRadius: Interpolation in arm "
         << arm << " failed. Exiting with 0. return." << endl;
    return 0.;
  }
  if (isnan(r) || isinf(r)) return 0.;
  gsl_interp_accel_free(a);
  return r;
}

void Astro::SetSurfaceDensityModel(string surfacedensitymodel) {
  if(!surfacedensitymodel.compare("CaseBhattacharya")) SURFACEDENSITYMODEL = 0;
  else if(!surfacedensitymodel.compare("IusifovKucuk")) SURFACEDENSITYMODEL = 1;
  else if(!surfacedensitymodel.compare("Linear")) SURFACEDENSITYMODEL = 2;
  else {
    std::cout<<"Astro::SetSurfaceDensityModel: Specify supported surface density model (currently 'CaseBhattacharya', 'IusifovKucuk' and 'Linear')! returning 0 value."<<std::endl;
  }
  return;
}

void Astro::SetCentralStructureModel(string centralstructuremodel) {
  if(!centralstructuremodel.compare("Bar")) CENTRALSTRUCTUREMODEL = 0;
  else if(!centralstructuremodel.compare("Disk")) CENTRALSTRUCTUREMODEL = 1;
  else {
    std::cout<<"Astro::SetCentralStructureModel: Specify central structure model (currently 'Bar' and 'Disk')! returning 0 value."<<std::endl;
  }
  return;
}

void Astro::SetSpiralArmModel(string spiralarmmodel) {
  if(!spiralarmmodel.compare("Vallee")) SPIRALARMMODEL = 0;
  else if(!spiralarmmodel.compare("TaylorCordes")) SPIRALARMMODEL = 1;
  else if(!spiralarmmodel.compare("Flat")) SPIRALARMMODEL = 2;
  else {
    std::cout<<"Astro::SetSpiralArmModel: Specify supported spiral arm model (currently 'Vallee' and 'TaylorCordes')! returning 0 value."<<std::endl;
  }
}

void Astro::SetMainSequenceBubbleModel(string msbubblemodel) {
  if(!msbubblemodel.compare("Weaver")) MSBUBBLEMODEL = 0;
  else if(!msbubblemodel.compare("Chen")) MSBUBBLEMODEL = 1;
  else {
    std::cout<<"Astro::SetMainSequenceBubbleModel: Specify supported main sequence wind bubble model (currently 'Weaver' and 'Chen')! Defaulting to Weaver's model."<<std::endl;
  }
}

vector< vector<double> > Astro::DiceGalacticPositions(int n) {
  vector<double> r = GetRandomGalactocentricRadii(n);
  vector< vector<double> > v;
  double x,y,z;
  for(int i=0;i<n;i++) {
    v.push_back(vector<double>());
    GalacticPositionXY(r[i],GetRandomArm(),x,y);
    RandomGaussianShift(r[i], armWidth, x, y);
    z = GetRandomZValue(scaleHeight);
    v[v.size()-1].push_back(x);
    v[v.size()-1].push_back(y);
    v[v.size()-1].push_back(z);
  }
  return v;
}


/**
 * return a random spiral arm (identifier, i.e. an unsigned integer)
 */
unsigned int Astro::GetRandomArm() {
  unsigned int size = ArmsVector.size();
  double v = fUtils->Random();
  for(unsigned int i=1;i<=size;i++) {
    double x = (double)i/size;
    if(x>v) return ArmsVector[i-1];
  }
  return 0;
}

/**
 * Return x,y coordinates of position along the i-th spiral arm or in the central structure at
 * the Galactocentric radius r
 */
void Astro::GalacticPositionXY(double r, int arm, double &x, double &y) {
  x = y = 0.;
  /* Central Structure */
  if(r<3.6 + fUtils->GaussianRandom(0.36, 0.,1)[0]) {
    if(!CENTRALSTRUCTUREMODEL) Bar(r,x,y);
    else if(CENTRALSTRUCTUREMODEL==1) Disk(r,x,y);
    else {
      std::cout<<"Astro::LogarithmicSpiral: Please specify proper central galactic structure. Options: Bar,Bulge. Exiting"<<std::endl;
    }
    return;
  }
  /* Spiral arms */
  if(!SPIRALARMMODEL) ValleeSpiral(r,arm,x,y);
  else if(SPIRALARMMODEL==1) TaylorCordesSpiral(r,arm,x,y);
  else if(SPIRALARMMODEL==2) Disk(r,x,y);
  else {
    std::cout<<"Astro::LogarithmicSpiral: Please specify proper available spiral arm model. Options: Vallee,TaylorCordes. Returning zero coordinates"<<std::endl;
    x = y = 0.;
  }
  return;
}

void Astro::PositionOnSpiralArmAngular(double theta, int arm, double &x, double &y) {
  if(!SPIRALARMMODEL)  ValleeSpiral(ValleeSpiralAngular(theta,arm),arm,x,y);
  else if(SPIRALARMMODEL==1) TaylorCordesSpiral(TaylorCordesSpiralAngular(theta,arm),arm,x,y);
  else {
    std::cout<<"Astro::LogarithmicSpiral: Please specify proper available spiral arm model. Options: Vallee,TaylorCordes. Returning zero coordinates"<<std::endl;
    x = y = 0.;
  }

//  std::cout<<"theta -> r: "<<theta<<" -> "<<ValleeSpiralAngular(theta,arm)<<" ->  x,y : "<<x<<" "<<y<<std::endl;
  return;
}

vector<double> Astro::GetRandomGalactocentricRadii(int n) {
  double rmin = 0.;
  double rmax = 20.;
  if(SURFACEDENSITYMODEL==2) return fUtils->LinearRandom(0.,rmin,rmax,n);
  vector< vector<double> > v;
  int steps = 200;
  double dr = (rmax-rmin)/steps;
  for(int i=0;i<steps;i++) {
    double r,val;
    r  = rmin+i*dr;

    if(!SURFACEDENSITYMODEL) val = CaseBhattacharyaProfile(r);
    else if(SURFACEDENSITYMODEL==1) val = IusifovKucukProfile(r);
    else break;
    v.push_back(vector<double>());
    v[v.size()-1].push_back(r);
    v[v.size()-1].push_back(val);
  }
  return fUtils->CustomFunctionRandom(v,v[0][0],v[v.size()-1][0],n);
}

double Astro::CaseBhattacharyaProfile(double r) {
  return 2.*pi*r*pow(r/8.5,2.)*exp(-5.1*(r-8.5)/8.5);
}
double Astro::IusifovKucukProfile(double r) {
  return 2.*pi*r*pow((r+0.55)/(8.5+0.55),1.64)*exp(-4.01*(r-8.5)/(8.5+0.55));
}
double Astro::GetRandomZValue(double scaleHeight) {
  return fUtils->SignRandom(1)[0]*fUtils->ExponentialRandom(scaleHeight,0.,5.*scaleHeight,1)[0];
}

/**
 * gaussian shift along concentric circle (retains r-distribution), but
 * looks odd with TaylorCordes spirals
 */
void Astro::RandomTangentialShift(double r, double width, double &x, double &y) {
  double dtheta = 2.*pi*fUtils->GaussianRandom(width, 0.,1)[0]/r;
  double x_r = x*cos(dtheta)-y*sin(dtheta);
  double y_r = x*sin(dtheta)+y*cos(dtheta);
  x = x_r;
  y = y_r;
  return;
}
/**
 * gaussian shift along concentric circle (retains r-distribution), but
 * looks odd with TaylorCordes spirals
 */
void Astro::RandomGaussianShift(double r, double width, double &x, double &y) {
  x += fUtils->GaussianRandom(sqrt(width), 0.,1)[0];
  y += fUtils->GaussianRandom(sqrt(width), 0.,1)[0];
  return;
}

void Astro::ValleeSpiral(double r, int arm, double &x, double &y) {
  if(arm<1) {
    std::cout<<"Astro::ValleeSpiral: Arm numbers start at 1 (here: "<<arm<<"). Returning zero value"<<std::endl;
    return;
  }
  bool ARMNOTTHERE = true;
  for(unsigned int i=0;i<ArmsVector.size();i++) {if(arm==ArmsVector[i]) ARMNOTTHERE = false;}
  if(ARMNOTTHERE==true) {
    std::cout<<"Astro::ValleeSpiral: Arm number "<<arm<<" not there. Returning zero value."<<std::endl;
    return;
  }
  double psi0 = (pi/180.)*psi0V[arm-1];
  double theta0 = (pi/180.)*theta0V[arm-1];
  double r0 = r0V[arm-1];
  double delta_theta = log(r/r0)/tan(psi0);
  double theta = theta0 + delta_theta;
  x = r*cos(theta);
  y = r*sin(theta);
  return;
}

double Astro::ValleeSpiralAngular(double theta, int arm) {
  if(arm<1) {
    std::cout<<"Astro::ValleeSpiral: Arm numbers start at 1 (here: "<<arm<<"). Returning zero value"<<std::endl;
    return 0.;
  }
  bool ARMNOTTHERE = true;
  for(unsigned int i=0;i<ArmsVector.size();i++) {if(arm==ArmsVector[i]) ARMNOTTHERE = false;}
  if(ARMNOTTHERE==true) {
    std::cout<<"Astro::ValleeSpiral: Arm number "<<arm<<" not there. Returning zero value."<<std::endl;
    return 0.;
  }
  double psi0 = (pi/180.)*psi0V[arm-1];
  double theta0 = (pi/180.)*theta0V[arm-1];
  double r0 = r0V[arm-1];
  double delta_theta = (pi/180.)*theta-theta0;
  double r = r0*exp(delta_theta*tan(psi0));
  return r;
}

void Astro::TaylorCordesSpiral(double r, int arm, double &x, double &y) {
  if(arm<1) {
    std::cout<<"Astro::TaylorCordesSpiral: Arm numbers start at 1 (here: "<<arm<<"). Returning zero value"<<std::endl;
    return;
  }
  bool ARMNOTTHERE = true;
  for(unsigned int i=0;i<ArmsVector.size();i++) {if(arm==ArmsVector[i]) ARMNOTTHERE = false;}
  if(ARMNOTTHERE==true) {
    std::cout<<"Astro::TaylorCordesSpiral: Arm number "<<arm<<" not there. Returning zero value."<<std::endl;
    return;
  }
  double theta = EvalTaylorCordesArmTheta(r,arm);
  if(!theta) return;
  x = r*cos(theta);
  y = r*sin(theta);
  return;
}

double Astro::TaylorCordesSpiralAngular(double theta, int arm) {
  if(arm<1) {
    std::cout<<"Astro::TaylorCordesSpiral: Arm numbers start at 1 (here: "<<arm<<"). Returning zero value"<<std::endl;
    return 0.;
  }
  bool ARMNOTTHERE = true;
  for(unsigned int i=0;i<ArmsVector.size();i++) {if(arm==ArmsVector[i]) ARMNOTTHERE = false;}
  if(ARMNOTTHERE==true) {
    std::cout<<"Astro::TaylorCordesSpiral: Arm number "<<arm<<" not there. Returning zero value."<<std::endl;
    return 0.;
  }
  double r = EvalTaylorCordesArmGalactocentricRadius(theta,arm);
  return r;
}

/**
 * Get random xy Coordinates on a circular disk
 */
void Astro::Disk(double r, double &x, double &y) {
  x = y = 0.;
  double theta = 2*pi*fUtils->Random();
  x = r*cos(theta);
  y = r*sin(theta);
  return;
}

/**
 * Get xy Coordinates on a central bar inclined with an angle 'angle' clockwise
 */
void Astro::Bar(double r, double &x, double &y) {
  double angle = 90. - barAngle ; //90 - x because angle is relative to Sun-GC line (y-dir)
  angle *= pi/180.;
  x = y = 0.;
  r *= fUtils->SignRandom(1)[0];
  x = r*cos(angle);
  y = r*sin(angle);
  return;
}

/**
 * radial distribution that gives a flat appearance of the Astro
 */
double Astro::LinearSurfDens(double *x, double *par) {
  double xx = x[0];
  if(xx<=3.) return 0.5*pow(xx,0.5);
  else return 0.87;
}
