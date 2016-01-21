#include "Astro.h"

Astro::Astro() {

  fUtils = new Utils();
  QUIETMODE = 0;
  /* Standard assumption of 4 spiral arms */
  ArmsVector.push_back(1);
  ArmsVector.push_back(2);
  ArmsVector.push_back(3);
  ArmsVector.push_back(4);

  accArm1 = gsl_interp_accel_alloc();
  accArm2 = gsl_interp_accel_alloc();
  accArm3 = gsl_interp_accel_alloc();
  accArm4 = gsl_interp_accel_alloc();

  TaylorCordesArm1 = gsl_spline_alloc(gsl_interp_linear, 6);
  TaylorCordesArm2 = gsl_spline_alloc(gsl_interp_linear, 7);
  TaylorCordesArm3 = gsl_spline_alloc(gsl_interp_linear, 7);
  TaylorCordesArm4 = gsl_spline_alloc(gsl_interp_linear, 7);

  TaylorCordesArm1Inv = gsl_spline_alloc(gsl_interp_linear, 6);
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
  //x_a1[6]=5.8101; y_a1[6]=330.-90.; FIXME this last point breaks interoplation

  gsl_spline_init(TaylorCordesArm1, x_a1, y_a1, 6);
  gsl_spline_init(TaylorCordesArm1Inv, y_a1, x_a1, 6);

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

  rMinInternalBoundary = 1.e-5*pc_to_cm;
  rMaxInternalBoundary = 100.*pc_to_cm;
  tMinInternalBoundary = 0.01*yr_to_sec;
  tMaxInternalBoundary = 1.e7*yr_to_sec;

  thinshellsolutiontmin = 0.;
  thinshellsolutiontmax = 0.;

  atsrad = gsl_interp_accel_alloc();
  atsvel = gsl_interp_accel_alloc();
  aredr = gsl_interp_accel_alloc();
  avedr = gsl_interp_accel_alloc();
  aredf = gsl_interp_accel_alloc();
  avedf = gsl_interp_accel_alloc();
  arstr = gsl_interp_accel_alloc();
  avstr = gsl_interp_accel_alloc();
  arstf = gsl_interp_accel_alloc();
  avstf = gsl_interp_accel_alloc();

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


vector< vector<double> > Astro::CreateDensityProfile(vector<double> pars,
                                                    double rmin, double rmax) {
  vector< vector<double> > dProfile;

  if(pars.size()!=8) {
    cout << "Astro::CreateDensityProfile: wrong number of parameters ("
         << pars.size() << "). Parameter list takes exactly 8 "
            "parameters: " << endl;
    cout << " - [0] red giant wind radius (pc) " << endl;
    cout << " - [1] mass loss rate of progenitor red giant (mSol/yr) " << endl;
    cout << " - [2] red giant wind speed (cm/s) " << endl;
    cout << " - [3] radius of main-seqence wind bubble (MSWB) (pc) " << endl;
    cout << " - [4] density inside MSWB (cm^-3) " << endl;
    cout << " - [5] density of the ISM (cm^-3) " << endl;
    cout << " - [6] width of the shell at RGW->MSWB transition (pc) " << endl;
    cout << " - [7] width of the shell at MSWB->ISM transition (pc) " << endl;
    cout << "Returning empty vector!" <<endl;
    return dProfile;
  }
  if(rmin >= rmax) {
    cout << "Astro::CreateDensityProfile: lower radius boundary larger than "
            "upper boundary (" << rmin << " > " << rmax <<"). Returning "
            "empty vector. " <<endl;
    return dProfile;
  }
  double RGWRadius = pars[0]*pc_to_cm;
  double mDotRGWind = pars[1]*mSol/yr_to_sec;
  double vRGWind = pars[2];
  double MSBubbleRadius = pars[3]*pc_to_cm;
  double MSBubbleDensity = pars[4];
  double ISMDensity = pars[5];
  double rSWRGW = pars[6]*pc_to_cm;
  double rSWMS = pars[7]*pc_to_cm;
  rmin *= pc_to_cm;
  rmax *= pc_to_cm;
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

  double stepSizeNormal=log(rmax/rmin)/1000.;
  double stepSizeFine=stepSizeNormal/10.;
  double stepSize=stepSizeNormal;
  double n = 0.;
  for(double r=rmin;r<rmax;r=pow(10.,log10(r)+stepSize)) {
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
    fUtils->TwoDVectorPushBack(r,n,dProfile);
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

double Astro::HTotalDensity(double x, double y, double z) {
  
  double n = H2Density(x,y,z,false) + HIDensity(x,y,z,false) + HIIDensity(x,y,z,false);
  if(SPIRALARMMODEL != 2) n = ModulateGasDensityWithSpirals(n,x,y,z);
  return n;
}

double Astro::HIDensity(double x, double y, double z, bool MODULATE) {
  double R = sqrt(pow(x,2.)+pow(y,2.));
  double n = 0.;
  if( R < 3.) {
    n = HICMZ(x,y,z);
    n += HIInnerDisk(x,y,z);
  }
  else {
    double norm = CalculateHINorm(R);
    double FWHM = GetHIFWHM(R)/1000.;
    n = 0.7*exp(-pow(z/(0.55*FWHM),2.));
    n += 0.19*exp(-pow(z/(1.38*FWHM),2.));
    n += 0.11*exp(-fabs(z)/(1.75*FWHM));
    n *= norm;
  }
  if(MODULATE && SPIRALARMMODEL != 2) n = ModulateGasDensityWithSpirals(n,x,y,z);
  return n;
}

/**
 *
 */
double Astro::H2Density(double x, double y, double z, bool MODULATE) {
  double R = sqrt(pow(x,2.)+pow(y,2.));
  double n = 0.;
  if( R < 3.) {
    n = H2CMZ(x,y,z);
    n += H2InnerDisk(x,y,z);
    n *= 2.;
  }
  else {
    double norm = CalculateH2Norm(R);
    double sigma = GetH2FWHM(R)/(2.*sqrt(2.*log(2.)));
    sigma /= 1000.;
    n = norm*exp(-pow(z/(sqrt(2.)*sigma),2.));
  }
  if(MODULATE && SPIRALARMMODEL != 2) n = ModulateGasDensityWithSpirals(n,x,y,z);
  return n;
}


/**
 *
 */
double Astro::HIIDensity(double x, double y, double z, bool MODULATE) {
  double Rsquare = x * x + y * y;
  double n = 0.;
  if( Rsquare < 9.) {
    n = WarmHIIDensity(x,y,z);
    n += HotHIIDensity(x,y,z);
    n += VeryHotHIIDensity(x,y,z);
  }
  else n = HIIThickDisk(x,y,z) + HIIThinDisk(x,y,z);

  if(MODULATE && SPIRALARMMODEL != 2) n = ModulateGasDensityWithSpirals(n,x,y,z);
  return n;
}

// TODO: can't implement this solution because infos in paper are not complete!
//double Astro::HIISpiralArms(double x, double y, double z, double fna, double hHa, double wwa, double F) {
//  double Aa = 8.5;
//  double sech2 = 1. / cosh(r - Aa);
//  double f = {}
//  for(i = 0; i < 5; i++) {
//    double sech1 = 1. / cosh( z / hjha[i] );
//    double Ga = f[i] * sech1 * sech1 * sech2 * sech2;  
//    
//  }
//}


double Astro::HIIThickDisk(double x, double y, double z) {

  double n1 = 3.3e-2 / 0.95;
  double H1 = 0.95;
  double A1 = 17.;
  double r = sqrt(x * x + y * y);
  double Rsol = 8.5;
  double sech = 1. / cosh( z / H1 );
  double g1 = 0.;
  if( r < A1 ) {
    g1 = cos(pi * r / (2. * A1)) / cos(pi * Rsol / (2. * A1));
  }
  
  return n1 * g1 * sech * sech;
}

double Astro::HIIThinDisk(double x, double y, double z) {

  double n2 = 9.e-2;
  double H2 = 0.14;
  double A2 = 3.7;
  double A2corr = 1.8;
  double sech = 1. / cosh(z / H2);
  double r = sqrt(x * x + y * y );
  double g2 = exp(- (r - A2) * (r - A2) / (A2corr * A2corr));
  return n2 * g2 * sech * sech;

}


double Astro::WarmHIIDensity(double x, double y, double z) {
  double y3 = -1.e-2;
  double z3 = -2.e-2;
  double L3 = 0.145;
  double H3 = 2.6e-2;
  double L2 = 3.7;
  double H2 = 0.14;
  double L1 = 17.;
  double H1 = 0.95;
  double r = sqrt( x * x + y * y );

  double sech1 = 1. / cosh( z / H1 );
  double sech2 = 1. / cosh( z / H2 );

  double n = exp(-( x * x + (y - y3)*(y - y3) )/( L3 * L3 ))
           * exp(-( z - z3 )*( z - z3 )/( H3 * H3 ))
           + 9.e-3 * exp(-( r - L2 )*( r - L2 )/( 0.25 * L2 * L2 )) *sech2 *sech2;
  if(r<L1) {
    n += 5.e-3 * cos( 0.5 * pi * r / L1 ) * sech1 * sech1;
  }
  n *= 8.;
  return n;
}


double Astro::H_I_2_CMZ(double x, double y, double z, double H_c, double norm) {

  double x_c = -5.e-2;
  double y_c = 5.e-2;
  double theta_c = pi * 70. / 180.;
  double X_c = 0.125;
  double L_c = 0.137;

  double X = ( x - x_c ) * cos( theta_c ) + ( y - y_c ) * sin( theta_c );
  double Y = -( x - x_c ) * sin( theta_c ) + ( y - y_c ) * cos( theta_c );

  double n = exp(-pow( sqrt(X * X + (2.5 * Y) * (2.5 * Y) ) - X_c, 4.) / L_c)
              * exp(- z * z / (H_c * H_c)) * norm;

  return n;
}

double Astro::H2CMZ(double x, double y, double z) {
  double norm = 150.;
  double H_c = 1.8e-2;
  return H_I_2_CMZ(x,y,z,H_c,norm);
}

double Astro::HICMZ(double x, double y, double z) {
  double norm = 8.8;
  double H_c = 5.4e-2;
  return H_I_2_CMZ(x,y,z,H_c,norm);
}

double Astro::H_I_2_InnerDisk(double x, double y, double z, double H_d, double norm) { 

  double X_d = 1.2;
  double L_d = 0.438;

  double alpha = pi * 13.5 / 180.;
  double beta = pi * 20. / 180.;
  double theta_d = pi * 48.5 / 180.;
  
  double sa = sin(alpha);
  double sb = sin(beta);
  double st = sin(theta_d);

  double ca = cos(alpha);
  double cb = cos(beta);
  double ct = cos(theta_d);

  double X = x * cb * ct - y * (sa * sb * ct - ca * st) - z * (ca * sb * ct + sa * st);
  double Y = - x * cb * st + y * (sa * sb * st + ca * ct) + z * (ca * sb * st - sa * ct);
  double Z = x * sb + y * sa * cb + z * ca * cb;

  double n = exp(-pow( sqrt(X * X + (3.1 * Y) * (3.1 * Y) ) - X_d, 4.) / L_d)
           * exp(- Z * Z / (H_d * H_d)) * norm;
 
  return n; 
}


double Astro::H2InnerDisk(double x, double y, double z) {
  double norm = 4.8;
  double H_d = 4.2e-2;
  return H_I_2_InnerDisk(x,y,z,H_d,norm);
}

double Astro::HIInnerDisk(double x, double y, double z) {
  double norm = 0.34;
  double H_d = 0.12;
  return H_I_2_InnerDisk(x,y,z,H_d,norm);
}

double Astro::HotHIIDensity(double x, double y, double z) {

  double r = sqrt( x * x + y * y );
  double n = pow( pow(9.e-3,2./3.) - 1.54e-17 * ( GravitationalPotential( r, z ) - GravitationalPotential(0., 0.) ),1.5);
  return n;

}

double Astro::VeryHotHIIDensity(double x, double y, double z) {

  double alpha = pi * 21. / 180.;
  double eta = y * cos (alpha) + z * sin (alpha);
  double zeta = -y * sin (alpha) + z * cos (alpha);
  double Lvh = 0.162;
  double Hvh = 0.09;
    
  double n = 0.29 * exp( - ( (x * x + eta * eta)/(Lvh * Lvh) + zeta * zeta / (Hvh * Hvh) ) );
  return n;
}
double Astro::GravitationalPotential(double r, double z) {

  double C1 = 8.887;
  double a1 = 6.5;
  double b1 = 0.26;
  double C2 = 3.;
  double a2 = 0.7;
  double C3  = 0.325;
  double a3 = 12.;
  double rh = 210.;
  double k = sqrt(1. + (a3 * a3 + r * r + z * z) / (rh * rh));
  double l = a1 + sqrt(z * z + b1 * b1);

  double phi = C1 / sqrt( r * r + l * l )
             + C2 / (a2 + sqrt(r * r + z * z))
             - C3 * log((k-1)/(k+1));

  phi *= -225.e5 * 225.e5;
  return phi;

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
  vector<double> xyzref;
  xyzref.resize(3);
  xyzref[0] = pars[2];
  xyzref[1] = pars[3];
  xyzref[2] = pars[4];
  int gasspecies = pars[5];
  int modulatewitharms = (int)pars[6];

  vector<double> xyzactual = GetCartesian(r,l,b,xyzref);
  double xActual = xyzactual[0];
  double yActual = xyzactual[1];
  double zActual = xyzactual[2];

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
vector<double> Astro::GetCartesian(double r, double l, double b, vector<double> xyzref) {

  double dl, db;
  vector<double> xyz;
  xyz.resize(3);
  double xref = xyzref[0];
  double yref = xyzref[1];
  double zref = xyzref[2];

  dl = db = 0.;
  dl = atan2(yref,xref);
  db = atan(zref/sqrt(xref*xref+yref*yref));
  double phi = (pi/180.)*l+dl;
  double theta = pi/2. - (pi/180.)*b+db;
  xyz[0] = xref - r*sin(theta)*cos(phi);
  xyz[1] = yref - r*sin(theta)*sin(phi);
  xyz[2] = zref + r*cos(theta);

  return xyz;
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
  bool Set = false;
  for(unsigned int i=0;i<ArmsVector.size();i++) {
    if(ArmsVector[i]==arm) {
      ArmsVector.erase(ArmsVector.begin()+i);
      Set = true;
      break;
    }
  }
  if(Set==false) std::cout<<"Arm "<<arm<<" already disabled."
                            "Nothing to do! "<<std::endl;
 
  return;
}

/**
 * Switch off an individual arm in Astro spiral model
 */
void Astro::EnableArm(int arm) {
  if(!QUIETMODE) std::cout<<"Enabling Arm no. "<<arm<<"."<<std::endl;
  bool AlreadyThere = false;
  for(unsigned int i=0;i<ArmsVector.size();i++) {
    if(ArmsVector[i] == arm) AlreadyThere = true;
  }
  if(AlreadyThere==true)  {
    std::cout<<"Arm "<<arm<<" already there! Nothing to do..."<<std::endl;
    if(!QUIETMODE) {
      for(unsigned int i=0;i<ArmsVector.size();i++) std::cout<<"Arm["<<i<<"]="<<ArmsVector[i]<<std::endl;
    }
    return;
  }
  ArmsVector.push_back(arm);
  sort(ArmsVector.begin(), ArmsVector.end());
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
    cout << "Astro::EvalTaylorCordesArmTheta: Requested Arm 0! (Arm Index "
            "starts a 1, not 0. Returning 0. " << endl;
    return 0.;
  }
  gsl_interp_accel *a = NULL;
  double rmin = 0.;
  double rmax = 0.;
  switch (arm){
    case 1:
      a = accArm1;
      rmin = x_a1[0]; rmax = x_a1[5];
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
  //cout << "-> " << arm << " " << r<<" " << rmin << " " <<rmax << endl;
  theta = fUtils->EvalSpline(r,TaylorCordesArms[arm - 1],a,
                             __func__,__LINE__);
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
  double thetamin = 0.;
  double thetamax = 0.;
  switch (arm){
    case 1:
      a = accArm1Inv;
      thetamin = y_a1[0]; thetamax = y_a1[5];
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
  r = fUtils->EvalSpline(theta,TaylorCordesArmsInv[arm - 1],a,
                             __func__,__LINE__);
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
    GalacticPositionXY(r[i],GetRandomArm(r[i]),x,y);
    if(x && y) {
      z = GetRandomZValue(scaleHeight);
      RandomGaussianShift(r[i], armWidth, x, y);
      v.push_back(vector<double>());
      v[v.size()-1].push_back(x);
      v[v.size()-1].push_back(y);
      v[v.size()-1].push_back(z);
    }
  }
  return v;
}


/**
 * return a random spiral arm (identifier, i.e. an unsigned integer)
 */
unsigned int Astro::GetRandomArm(double r) {
  vector<int> vint = ArmsVector;
    if(SPIRALARMMODEL==1) {

    //QUIETMODE = 1;
    if(r<x_a1[0] || r>x_a1[5]) DisableArm(1);
    if(r<x_a2[0] || r>x_a2[6]) DisableArm(2);
    if(r<x_a3[0] || r>x_a3[6]) DisableArm(3);
    if(r<x_a4[0] || r>x_a4[6]) DisableArm(4);
  }
  unsigned int size = ArmsVector.size();
  double v = fUtils->Random();
  unsigned int i = 1;
  for(;i<=size;i++) {
    //cout << "-- " << i << std::endl;
    double x = (double)i/size;
    if(x>v) break;
  }
  int arm = ArmsVector[i-1];
  ArmsVector = vint;
  return arm;
}

/**
 * Return x,y coordinates of position along the i-th spiral arm or in the central structure at
 * the Galactocentric radius r
 */
void Astro::GalacticPositionXY(double r, int arm, double &x, double &y) {
  x = y = 0.;
  if(!arm) {
    cout << "Astro::GalacticPositionXY: Arm 0 does not exits! Returing 0. "
            "values" << endl;
  }
  /* Central Structure */
  if(r<barRadius) {
    if(!CENTRALSTRUCTUREMODEL) Bar(r,x,y);
    else if(CENTRALSTRUCTUREMODEL==1) Disk(r,x,y);
    else {
      cout << "Astro::GalacticPositionXY: Please specify proper central galactic"
              "structure. Options: Bar,Bulge. Exiting" << endl;
    }
    return;
  }
  /* Spiral arms */
  if(!SPIRALARMMODEL) ValleeSpiral(r,arm,x,y);
  else if(SPIRALARMMODEL==1) TaylorCordesSpiral(r,arm,x,y);
  else if(SPIRALARMMODEL==2) Disk(r,x,y);
  else {
    std::cout<<"Astro::GalacticPositionXY: Please specify proper available spiral arm model. Options: Vallee,TaylorCordes. Returning zero coordinates"<<std::endl;
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
    fUtils->TwoDVectorPushBack(r,val,v);
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
  x += fUtils->GaussianRandom(width, 0.,1)[0];
  y += fUtils->GaussianRandom(width, 0.,1)[0];
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

//**  Here starts the dynamics model collection!  *//


/**
 * Radius and shock speed in the thin shell approximation. Implementation in
 * CalculateThinShellApproximation() member function.
 */
vector<double> Astro::ThinShellRadiusAndSpeed(double t) {
  vector<double> v,vzero;
  vzero.push_back(0.);
  vzero.push_back(0.);
  double R,V;
  t*=yr_to_sec;

  if(!thinshellradius.size() || !thinshellvelocity.size()) {
    cout << "Astro::ThinShellRadiusAndSpeed: You have to run the calculation "
            "first via CalculateThinShellApproximation()! Returning (r=0,v=0)!"
            << std::endl;
    return vzero;
  }
  if(t<=thinshellsolutiontmin || t>=thinshellsolutiontmax) {
    if(!QUIETMODE) {
      cout << "Astro::ThinShellRadiusAndSpeed: time outside solution boundaries "
           << "(t=" << t/yr_to_sec << " outside ["
           << thinshellsolutiontmin/yr_to_sec << ","
           << thinshellsolutiontmax/yr_to_sec << "]). " << endl;
      cout << "Returning (r=0,v=0)!" << std::endl;
    }
    return vzero;
  }
  R = fUtils->EvalSpline(t,thinshellradiuslookup,atsrad,__func__,__LINE__);
  V = fUtils->EvalSpline(t,thinshellvelocitylookup,atsvel,__func__,__LINE__);
  v.push_back(R);
  v.push_back(V);

  return v;
}

/**
 * This is the thin shell approximation for adiabatic shocks.
 * The implementation here is taken from the appendix in
 * Ptuskin&Zirakashvili 2005 (A&A 429, 755–765 (2005)).
 */
void Astro::CalculateThinShellApproximation(vector< vector<double> > dProfile,
                                            double E, double AdiabIndex,
                                            double tmin, double tmax,
                                            int steps) {

  if(!dProfile.size()) {
    cout << "Astro::CalculateThinShellApproximation: density profile empty."
            "Exiting." << endl;
    return;
  }
  if(!E) {
    cout << "Astro::CalculateThinShellApproximation: Blast energy is zero."
            "Exiting." << endl;
    return;
  }
  if(!AdiabIndex) {
    cout << "Astro::CalculateThinShellApproximation: Adiabatic index energy is"
            "zero. Exiting." << endl;
    return;
  }
  if(thinshellradius.size()) {
    fUtils->Clear2DVector(thinshellradius);
    fUtils->Clear2DVector(thinshellvelocity);
    gsl_spline_free(thinshellradiuslookup);
    gsl_spline_free(thinshellvelocitylookup);
    gsl_interp_accel_reset(atsrad);
    gsl_interp_accel_reset(atsvel);
  }
  /* fill a lookup of mass vs. radius, needed in the integral in (A.4.1) */
  vector< vector<double> > massintegrand;
  vector< vector<double> > mass;
  vector< vector<double> > fintegrand;
  vector< vector<double> > fintegral;
  for(unsigned int i=0;i<dProfile.size();i++) {
    double r = dProfile[i][0];
    double n = dProfile[i][1];
    fUtils->TwoDVectorPushBack(r,4.* pi * m_p_g * n * r * r,massintegrand);
  }
  mass = fUtils->IntegratedProfile(massintegrand);

  /* Fill a Graph Velocity vs. Radius first. Eq.   (A.4) */

  /* fintegrand is the integrand in (A.4.I). */
  double r0,r,m,v,t,y;
  double w = 6.*(AdiabIndex-1.)/(AdiabIndex+1.);
  for(unsigned int i=0;i<mass.size();i++) {
    r = mass[i][0];
    m = mass[i][1];
    fUtils->TwoDVectorPushBack(r,pow(r,w-1.)*m,fintegrand);
  }
  fintegral = fUtils->IntegratedProfile(fintegrand);

  gsl_spline *masslookup = fUtils->GSLsplineFromTwoDVector(mass);
  gsl_interp_accel *amass = gsl_interp_accel_alloc();

    /* Finally fill TGraph with V vs R, as well as 1./V vs R for time calculation. */
  vector< vector<double> > velocity,inversevelocity;
  double k=(AdiabIndex+1.)*sqrt(2.*w*E)/2.;
  for(unsigned int i=0;i<fintegral.size();i++) {
    r = fintegral[i][0];
    y = fintegral[i][1];
    m = fUtils->EvalSpline(r,masslookup,amass,__func__,__LINE__);
    v = k*sqrt(y/(m*m*pow(r,w)));
    if(!v) continue;

    fUtils->TwoDVectorPushBack(r,v,velocity);
    fUtils->TwoDVectorPushBack(r,1./v,inversevelocity);
  }

  /* This fills the t vs R vector (A.4.II)*/
  gsl_spline *timelookup = fUtils->GSLsplineFromTwoDVector(
                                    fUtils->IntegratedProfile(inversevelocity)
                                   );
  gsl_interp_accel *atime = gsl_interp_accel_alloc();

  /* but what we want is not t(r) but rather r(t) -> invert the TGraph! */
  for(unsigned int i=1;i<inversevelocity.size();i++) {
    r = inversevelocity[i][0];
    t = fUtils->EvalSpline(r,timelookup,atime,__func__,__LINE__);
    r0 = velocity[i][0];
    v = velocity[i][1];
    if(r!=r0) {
      cout << "Astro::CalculateThinShellApproximation: Radius and Velocity "
                 "calculations do not match! Exiting!" << endl;
      return;
    }
    if(r<rMinInternalBoundary||r>rMaxInternalBoundary) continue;
    if(t<tMinInternalBoundary||t>tMaxInternalBoundary) continue;

    fUtils->TwoDVectorPushBack(t,r,thinshellradius);
    fUtils->TwoDVectorPushBack(t,v,thinshellvelocity);
  }
  thinshellradius = fUtils->SortTwoDVector(thinshellradius,0);
  thinshellvelocity = fUtils->SortTwoDVector(thinshellvelocity,0);
  thinshellsolutiontmin = thinshellradius[0][0];
  thinshellsolutiontmax = thinshellradius[thinshellradius.size()-1][0];

  thinshellradiuslookup = fUtils->GSLsplineFromTwoDVector(thinshellradius);
  thinshellvelocitylookup = fUtils->GSLsplineFromTwoDVector(thinshellvelocity);

  // make the r and v time profiles
  if(steps && tmin && tmax) {
    if(tmin>=tmax) {
      cout << "Astro::CalculateThinShellApproximation: Can't create profiles"
              "because tmin >= tmax (" << tmin <<" >= "<< tmax <<
              ") Exiting!" << endl;
      return;
    }
    double logtmin = log10(tmin);
    double logtmax = log10(tmax);
    double dlogt = (logtmax-logtmin)/steps;
    fUtils->Clear2DVector(forwardshockradiusprofile);
    fUtils->Clear2DVector(forwardshockvelocityprofile);
    vector<double> x;
    double t=0.;
    for(double logt = logtmin ; logt < logtmax ; logt += dlogt) {
      t = pow(10.,logt);
      x = ThinShellRadiusAndSpeed(t);
      fUtils->TwoDVectorPushBack(t,x[0]/pc_to_cm,forwardshockradiusprofile);
      fUtils->TwoDVectorPushBack(t,x[1],forwardshockvelocityprofile);
    }
  }
  return;
}

vector<double> Astro::TrueloveMcKeeRadiusAndSpeed(double t, bool REVERSE) {
  vector<double> v,vzero;
  vzero.push_back(0.);
  vzero.push_back(0.);
  double R,V;
  t*=yr_to_sec;

  if(!r_edforward.size() || !v_edforward.size() || !r_stforward.size() ||
     !v_stforward.size() || !r_edreverse.size() || !v_edreverse.size() ||
     !r_streverse.size() || !v_streverse.size() ||
     truelovemckeeparams.size() != 23) {

     cout << "Astro::TrueloveMcKeeRadiusAndSpeed: You have to run the calculation "
             "first via CalculateThinShellApproximation()! Returning (r=0,v=0)!"
          << std::endl;
  }
  double tmin = 0.;
  double tmax = 0.;
  double ttrans = 0.;
  double sedovtime = truelovemckeeparams[20];
  double coreexittime = truelovemckeeparams[21];
  int index = (int)truelovemckeeparams[22];

  if(REVERSE) (index>5) ? ttrans = coreexittime: ttrans = sedovtime;
  else ttrans = sedovtime;

  if(t<ttrans) {
    if(REVERSE) {
      tmin = truelovemckeetminred;
      tmax = truelovemckeetmaxred;
    }
    else {
      tmin = truelovemckeetminfed;
      tmax = truelovemckeetmaxfed;
    }
  }
  else {
    if(REVERSE) {
      tmin = truelovemckeetminrst;
      tmax = truelovemckeetmaxrst;
    }
    else {
      tmin = truelovemckeetminfst;
      tmax = truelovemckeetmaxfst;
    }
  }
  if(tmin >= tmax) {
    cout << "Astro::TrueloveMcKeeRadiusAndSpeed: Something wrong with time "
            "boundaries "
         << "(tmin=" << tmin/yr_to_sec << " , tmax = "
         << tmax/yr_to_sec << "). Returning (r=0,v=0)!"
         << std::endl;
    return vzero;
  }
  if(t<tmin || t>tmax) {
    if(!QUIETMODE) {
      cout << "Astro::TrueloveMcKeeRadiusAndSpeed: time outside solution "
           << "boundaries (t=" << t/yr_to_sec << " outside ["
           << tmin/yr_to_sec << ","
           << tmax/yr_to_sec << "]). " << endl;
      cout << "Returning (r=0,v=0)!" << std::endl;
    }
    return vzero;
  }

  if(t<ttrans) {
    if(REVERSE) {
      R = fUtils->EvalSpline(t,r_edreverselookup,aredr,__func__,__LINE__);
      V = fUtils->EvalSpline(t,v_edreverselookup,avedr,__func__,__LINE__);
    }
    else {
      R = fUtils->EvalSpline(t,r_edforwardlookup,aredf,__func__,__LINE__);
      V = fUtils->EvalSpline(t,v_edforwardlookup,avedf,__func__,__LINE__);
    }
  }
  else {
    if(REVERSE) {
      R = fUtils->EvalSpline(t,r_streverselookup,arstr,__func__,__LINE__);
      V = fUtils->EvalSpline(t,v_streverselookup,avstr,__func__,__LINE__);
    }
    else {
      R = fUtils->EvalSpline(t,r_stforwardlookup,arstf,__func__,__LINE__);
      V = fUtils->EvalSpline(t,v_stforwardlookup,avstf,__func__,__LINE__);
    }
  }
  v.push_back(R);
  v.push_back(V);
  return v;
}

void Astro::CalculateTrueloveMcKeeSolution(vector<double> pars, double tmin,
                                           double tmax, int steps) {

  if(pars.size()!=4) {
    cout << "Astro::CalculateTrueloveMcKeeSolution: wrong number of parameters ("
         << pars.size() << "). Parameter list takes exactly 4 "
            "parameters: " << endl;
    cout << " - [0] ambient density (cm^-3) " << endl;
    cout << " - [1] SN blast energy (erg) " << endl;
    cout << " - [2] SN ejecta mass (mSol) " << endl;
    cout << " - [3] density index " << endl;
    cout << "Exiting!" <<endl;
    return;
  }

  CalculateTrueLoveMcKeeParams(pars);
  if(!truelovemckeeparams.size()) {
    std::cout<<"Astro::CalculateTrueloveMcKeeSolution: "
               "No paramters given. Exiting."<<std::endl;
    return;
  }


  double charradiustruelovemckee = truelovemckeeparams[17];
  double chartimetruelovemckee = truelovemckeeparams[18];
  double charvelocitytruelovemckee = truelovemckeeparams[19];


  bool SKIPFWS,SKIPRVS;
  double x0,x1,dx,y,yr,t,tr,r,rr,v,vr;
  int index = (int)truelovemckeeparams[22];
  if(index<5) {
    x0 = rMinInternalBoundary/charradiustruelovemckee;
    x1 = rMaxInternalBoundary/charradiustruelovemckee;
  }
  else {
    x0 = tMinInternalBoundary/chartimetruelovemckee;
    x1 = tMaxInternalBoundary/chartimetruelovemckee;
  }
  dx=log10(x1/x0)/steps;
  for(double x=x0;x<x1;x=pow(10.,log10(x)+dx)) {
    SKIPFWS=SKIPRVS=false;
    if(TrueloveMcKeeEjectaDominatedForwardShock(x,y,v)) SKIPFWS=true;
    if(TrueloveMcKeeEjectaDominatedReverseShock(x,yr,vr)) SKIPRVS=true;
    if(index<5) {r=x;t=y;tr=yr,rr=x;}
    else {r=y;t=x;rr=yr;tr=x;}
    t*=chartimetruelovemckee;
    r*=charradiustruelovemckee;
    v*=charvelocitytruelovemckee;
    tr*=chartimetruelovemckee;
    rr*=charradiustruelovemckee;
    vr*=charvelocitytruelovemckee;
    if(SKIPFWS==false) {
      fUtils->TwoDVectorPushBack(t,r,r_edforward);
      fUtils->TwoDVectorPushBack(t,v,v_edforward);
    }
    if(SKIPRVS==false) {
      fUtils->TwoDVectorPushBack(tr,rr,r_edreverse);
      fUtils->TwoDVectorPushBack(tr,vr,v_edreverse);
    }
  }
  x0 = tMinInternalBoundary/chartimetruelovemckee;
  x1 = tMaxInternalBoundary/chartimetruelovemckee;
  dx=log10(x1/x0)/steps;
  for(double x=x0;x<x1;x=pow(10.,log10(x)+dx)) {
    SKIPFWS=SKIPRVS=false;
    if(TrueloveMcKeeSedovTaylorPhaseForwardShock(x,r,v)) SKIPFWS=true;
    if(TrueloveMcKeeSedovTaylorPhaseReverseShock(x,rr,vr)) SKIPRVS=true;
    t=x*chartimetruelovemckee;
    r*=charradiustruelovemckee;
    v*=charvelocitytruelovemckee;
    tr=x*chartimetruelovemckee;
    rr*=charradiustruelovemckee;
    vr*=charvelocitytruelovemckee;
    if(SKIPFWS==false) {
      fUtils->TwoDVectorPushBack(t,r,r_stforward);
      fUtils->TwoDVectorPushBack(t,v,v_stforward);
    }
    if(SKIPRVS==false) {
      fUtils->TwoDVectorPushBack(tr,rr,r_streverse);
      fUtils->TwoDVectorPushBack(tr,vr,v_streverse);
    }
  }

  /* sort TGraphs in time and get minimum and maximum time boundaries */
  r_edforward = fUtils->SortTwoDVector(r_edforward,0);
  v_edforward = fUtils->SortTwoDVector(v_edforward,0);
  r_edreverse = fUtils->SortTwoDVector(r_edreverse,0);
  v_edreverse = fUtils->SortTwoDVector(v_edreverse,0);
  r_stforward = fUtils->SortTwoDVector(r_stforward,0);
  v_stforward = fUtils->SortTwoDVector(v_stforward,0);
  r_streverse = fUtils->SortTwoDVector(r_streverse,0);
  v_streverse = fUtils->SortTwoDVector(v_streverse,0);

  truelovemckeetminfed = r_edforward[0][0];
  truelovemckeetminfst = r_stforward[0][0];
  truelovemckeetmaxfed = r_edforward[r_edforward.size()-1][0];
  truelovemckeetmaxfst = r_stforward[r_stforward.size()-1][0];

  truelovemckeetminred = r_edreverse[0][0];
  truelovemckeetminrst = r_streverse[0][0];
  truelovemckeetmaxred = r_edreverse[r_edreverse.size()-1][0];
  truelovemckeetmaxrst = r_streverse[r_streverse.size()-1][0];

  r_edforwardlookup = fUtils->GSLsplineFromTwoDVector(r_edforward);
  v_edforwardlookup = fUtils->GSLsplineFromTwoDVector(v_edforward);
  r_edreverselookup = fUtils->GSLsplineFromTwoDVector(r_edreverse);
  v_edreverselookup = fUtils->GSLsplineFromTwoDVector(v_edreverse);

  r_stforwardlookup = fUtils->GSLsplineFromTwoDVector(r_stforward);
  v_stforwardlookup = fUtils->GSLsplineFromTwoDVector(v_stforward);
  r_streverselookup = fUtils->GSLsplineFromTwoDVector(r_streverse);
  v_streverselookup = fUtils->GSLsplineFromTwoDVector(v_streverse);

  // make the r and v time profiles
  if(steps && tmin && tmax) {
    if(tmin>=tmax) {
      cout << "Astro::CalculateThinShellApproximation: Can't create profiles"
              "because tmin >= tmax (" << tmin <<" >= "<< tmax <<
              ") Exiting!" << endl;
      return;
    }
    double logtmin = log10(tmin);
    double logtmax = log10(tmax);
    double dlogt = (logtmax-logtmin)/steps;
    fUtils->Clear2DVector(forwardshockradiusprofile);
    fUtils->Clear2DVector(forwardshockvelocityprofile);
    fUtils->Clear2DVector(reverseshockradiusprofile);
    fUtils->Clear2DVector(reverseshockvelocityprofile);
    vector<double> x;
    double t=0.;
    for(double logt = logtmin ; logt < logtmax ; logt += dlogt) {
      t = pow(10.,logt);
      x = TrueloveMcKeeRadiusAndSpeedForward(t);
      fUtils->TwoDVectorPushBack(t,x[0]/pc_to_cm,forwardshockradiusprofile);
      fUtils->TwoDVectorPushBack(t,x[1],forwardshockvelocityprofile);
      x = TrueloveMcKeeRadiusAndSpeedReverse(t);
      fUtils->TwoDVectorPushBack(t,x[0]/pc_to_cm,reverseshockradiusprofile);
      fUtils->TwoDVectorPushBack(t,x[1],reverseshockvelocityprofile);
    }
  }
  return;
}

int Astro::TrueloveMcKeeEjectaDominatedForwardShock(double x,
                                                    double &Y, double &V) {
  double l_ED = truelovemckeeparams[0];
  double phi_ED = truelovemckeeparams[2];
  double phi_EDeff = truelovemckeeparams[3];
  double f_n = truelovemckeeparams[14];
  double alpha = truelovemckeeparams[15];
  int index = (int)truelovemckeeparams[22];

  if(index<5.) {
    Y = 1.-((3.-index)/3.)*sqrt(phi_EDeff/(l_ED*f_n))*pow(x,1.5);
    if(Y<=0.) {Y=0.;V=0.;return 1;}
    Y = pow(Y,-2./(3.-index));
    Y *= sqrt(alpha/2.)*x/l_ED;

    V = 1.-((3.-index)/3.)*sqrt(phi_EDeff/(l_ED*f_n))*pow(x,1.5);
    if(V<=0.) {Y=0.;V=0.;return 1;}
    V = pow(V,(5.-index)/(3.-index));
    V *= sqrt(2./alpha)*l_ED;
    V /= 1.+(index/3.)*sqrt(phi_EDeff/(l_ED*f_n))*pow(x,1.5);
  }
  else {
    Y = 27.*pow(l_ED,index-2.)/(phi_ED*(4.*pi)*index*(index-3.));
    Y*= pow((10./3.)*(index-5.)/(index-3.),(index-3.)/2.);
    if(Y<=0.) {Y=0.;V=0.;return 1;}
    Y = pow(Y,1./index);
    Y *= pow(x,(index-3.)/index);

    V = 27./(4.*pi)*pow(index*(index-3.),-1.)*pow(l_ED,index-2.)/phi_ED;
    V *= pow((10./3.)*(index-5.)/(index-3.),(index-3.)/2.);
    if(V<=0.) {Y=0.;V=0.;return 1;}
    V = pow(V,1./index);
    V *= (index-3.)/index*pow(x,-3./index);
  }
  if(Y<0.) Y=0.;
  return 0;
}

int Astro::TrueloveMcKeeSedovTaylorPhaseForwardShock(double t,
                                                     double &R, double &V) {
  double tSt=truelovemckeeparams[4];
  double rSt=truelovemckeeparams[5];
  double eta=0.4;
  double zeta=2.026;

  double k = pow(rSt,1./eta);
  double l = pow(rSt,2.5);
  double m = sqrt(zeta)*(t-tSt);
  if(-m >= k || -m >= l) {
    R=V=0.;
    return 1;
  }
  R = pow(k+m,eta);
  V = 0.4*sqrt(zeta)*pow(l+m,-0.6);
  if(R<0.) {
    R=V=0.;
    return 1;
  }
  return 0;
}

int Astro::TrueloveMcKeeSedovTaylorPhaseReverseShock(double t, double &R, double &V) {
  double r_rSt = truelovemckeeparams[6];
  double a_rSt = truelovemckeeparams[8];
  double tSt = truelovemckeeparams[4];
  double v_rSt = truelovemckeeparams[7];
  double r_rCore = truelovemckeeparams[11];
  double a_rCore = truelovemckeeparams[13];
  double tCore = truelovemckeeparams[10];
  double v_rCore = truelovemckeeparams[12];
  int index = (int)truelovemckeeparams[22];

  double r_rRel,tRel,a_rRel,v_rRel;
  r_rRel=tRel=a_rRel=v_rRel=0.;
  if(index<3.) {
    tRel = tSt;
    r_rRel = r_rSt;
    a_rRel = a_rSt;
    v_rRel = v_rSt;
  }
  else if(index>5.) {
    tRel = tCore;
    r_rRel = r_rCore;
    a_rRel = a_rCore;
    v_rRel = v_rCore;
  }
  else {
    std::cout<<"Astro::TrueloveMcKeeSedovTaylorPhaseReverseShock: "
               "Whaaat?! I have to implement the reverse shock if index = 4!"
               "Exiting!"<<std::endl;
    R=V=0.;
    return 0;
  }
  R = t*(r_rRel/tRel - a_rRel*(t-tRel)-(v_rRel-a_rRel*tRel)*log(t/tRel));
  V = v_rRel+a_rRel*(t-tRel);
  if(R<0.) {
    R=V=0.;
    return 1;
  }
  return 0;
}

int Astro::TrueloveMcKeeEjectaDominatedReverseShock(double x,
                                                    double &Y, double &V) {
  double l_ED = truelovemckeeparams[0];
  double phi_ED = truelovemckeeparams[2];
  double f_n = truelovemckeeparams[14];
  double alpha = truelovemckeeparams[15];
  int index = (int)truelovemckeeparams[22];

  if(index<5.) {
    Y = 1.-((3.-index)/3.)*sqrt(phi_ED/(l_ED*f_n))*pow(x*l_ED,1.5);
    if(Y<=0.) {Y=0.;V=0.;return 1;}
    Y = pow(Y,-2./(3.-index));
    Y *= sqrt(alpha/2.)*x;

    V = 1.-((3.-index)/3.)*sqrt(phi_ED/f_n)*l_ED*pow(x,1.5);
    if(V<=0.) {Y=0.;V=0.;return 1;}
    V = pow(V,2./(3.-index));
    V *= sqrt(2.*phi_ED/(alpha*f_n))*l_ED*pow(x,1.5);
    V /= 1.+(index/3.)*sqrt(phi_ED/f_n)*l_ED*pow(x,1.5);
  }
  else {
    Y = 27./(4.*pi)*pow(index*(index-3.),-1.)*pow(l_ED,index-2.)/phi_ED;
    Y *= pow((10./3.)*(index-5.)/(index-3.),(index-3.)/2.);
    if(Y<=0.) {Y=0.;V=0.;return 1;}
    Y = pow(Y,1./index);
    V = Y*3./(index*l_ED)*pow(x,-3./index);
    if(V<=0.) {Y=0.;V=0.;return 1;}
    Y *= pow(x,(index-3.)/index)/l_ED;
  }
  if(Y<0.) Y=0.;
  return 0;
}

void Astro::CalculateTrueLoveMcKeeParams(vector<double> pars) {
  double n = pars[0];
  double e = pars[1];
  double mej = pars[2]*mSol;
  double index = pars[3];

  /* Eq (1) */
  double charradiustruelovemckee = pow(mej,1./3.)*pow(n*m_p_g,-1./3.);
  /* Eq (2) */
  double chartimetruelovemckee = pow(e,-1./2.)*pow(mej,5./6.)*pow(n*m_p_g,-1./3.);
  double charvelocitytruelovemckee=charradiustruelovemckee/chartimetruelovemckee;

  /* format: {l_ED,w_core,phi_ED,phi_EDeff,t*_ST,R*_ST,R*_r_ST,v~*_r_ST,a~*_r_ST,t*_ref,t*_core,r*_r_core,v~*_r_core,a~*_r_core,f_n,alpha,eta} */
  if(!index)          truelovemckeeparams = {  1.1,  0., 0.343, 0.0961, 0.495, 0.727, 0.545, 0.585, 0.106,  0.,    0.,    0.,    0.,    0.};
  else if(index== 2.) truelovemckeeparams = {  1.1,  0., 0.343, 0.0947, 0.387, 0.679, 0.503, 0.686,-0.151, 1.6,    0.,    0.,    0.,    0.};
  else if(index== 4.) truelovemckeeparams = {  1.1, 0.1, 0.343, 0.0791, 0.232, 0.587,    0.,    0.,    0.,  0.,   1.7,    0.,    0.,    0.};
  else if(index== 6.) truelovemckeeparams = { 1.39,  0.,  0.39,     0.,  1.04,  1.07,    0.,    0.,    0.,  0., 0.513, 0.541, 0.527, 0.112};
  else if(index== 7.) truelovemckeeparams = { 1.26,  0.,  0.47,     0., 0.732, 0.881,    0.,    0.,    0.,  0., 0.363, 0.469, 0.553, 0.116};
  else if(index== 8.) truelovemckeeparams = { 1.21,  0.,  0.52,     0., 0.605, 0.788,    0.,    0.,    0.,  0., 0.292, 0.413, 0.530, 0.139};
  else if(index== 9.) truelovemckeeparams = { 1.19,  0.,  0.55,     0., 0.523, 0.725,    0.,    0.,    0.,  0., 0.249, 0.371, 0.497, 0.162};
  else if(index==10.) truelovemckeeparams = { 1.17,  0.,  0.57,     0., 0.481, 0.687,    0.,    0.,    0.,  0., 0.220, 0.340, 0.463, 0.192};
  else if(index==12.) truelovemckeeparams = { 1.15,  0.,  0.60,     0., 0.424, 0.636,    0.,    0.,    0.,  0., 0.182, 0.293, 0.403, 0.251};
  else if(index==14.) truelovemckeeparams = { 1.14,  0.,  0.62,     0., 0.389, 0.603,    0.,    0.,    0.,  0., 0.157, 0.259, 0.354, 0.277};
  else {
    cout << "Astro::CalculateTrueLoveMcKeeParams: Density Index (index="
         << index << ") not supported! Returning nothing." << endl;
    truelovemckeeparams={};
    return;
  }

  double f_n,alpha,eta,w_core;
  w_core = truelovemckeeparams[1];
  f_n = (3./(4.*pi))*((1.-index/3.)/(1.-(index/3.)*pow(w_core,3.-index)));
  if(index<3.&&!w_core) alpha = (3.-index)/(5.-index);
  else if(index>5.&&!w_core) alpha = (3./5.)*(index-3.)/(index-5.);
  else alpha = (3.-index)/(5.-index)*(pow(w_core,-(5.-index))-index/5.)
               /(pow(w_core,-(3.-index))-index/3.)*pow(w_core,2.);
  eta = 0.;
  if(index>0.) eta = (index-3.)/index;
  truelovemckeeparams.push_back(f_n);
  truelovemckeeparams.push_back(alpha);
  truelovemckeeparams.push_back(eta);
  truelovemckeeparams.push_back(charradiustruelovemckee);
  truelovemckeeparams.push_back(chartimetruelovemckee);
  truelovemckeeparams.push_back(charvelocitytruelovemckee);
  truelovemckeeparams.push_back(truelovemckeeparams[4]*chartimetruelovemckee);
  truelovemckeeparams.push_back(truelovemckeeparams[10]*chartimetruelovemckee);
  truelovemckeeparams.push_back(index);

  if(truelovemckeeparams.size() != 23) {
    cout << "Astro::CalculateTrueLoveMcKeeParams: Something went wrong in the"
            "Calculation of the necessary parameter for the Truelove & McKee "
            "solution. Exiting!" << endl;
    truelovemckeeparams.clear();
    return;
  }
  if(!QUIETMODE) {
    vector<string> ParNames = {"l_ED","w_core","phi_ED","phi_EDeff","t*_ST",
                               "R*_ST","R*_r_ST","v~*_r_ST","a~*_r_ST","t*_ref",
                               "t*_core","r*_r_core","v~*_r_core","a~*_r_core",
                               "f_n","alpha","eta","charradiustruelovemckee",
                               "chartimetruelovemckee",
                               "charvelocitytruelovemckee","sedovtaylortime",
                               "coreexittime","index"};
    cout << ">> Paramaters for Truelove & McKee's solution (density index = "
         << index << "):" <<std::endl;
    for(unsigned int i=0;i<truelovemckeeparams.size();i++) {
      cout << "   " << ParNames[i] << ": " << truelovemckeeparams[i] << endl;
    }
    cout << endl;
  }
  return;
}

void Astro::CalculateForwardShockInRGWind(vector<double> pars, double tmin,
                                          double tmax, int steps) {

  if(pars.size()!=10) {
    cout << "Astro::CalculateForwardShockInRGWind: wrong number of parameters ("
         << pars.size() << "). Parameter list takes exactly 10 "
            "parameters: " << endl;
    cout << " - [0] red giant wind radius (pc) " << endl;
    cout << " - [1] mass loss rate of progenitor red giant (mSol/yr) " << endl;
    cout << " - [2] red giant wind speed (cm/s) " << endl;
    cout << " - [3] radius of main-seqence wind bubble (MSWB) (pc) " << endl;
    cout << " - [4] density inside MSWB (cm^-3) " << endl;
    cout << " - [5] density of the ISM (cm^-3) " << endl;
    cout << " - [6] width of the shell at RGW->MSWB transition (pc) " << endl;
    cout << " - [7] width of the shell at MSWB->ISM transition (pc) " << endl;
    cout << " - [8] SN ejecta mass (mSol) " << endl;
    cout << " - [9] SN blast energy (erg) " << endl;
    cout << "Returning empty vector!" <<endl;
    return;
  }

  vector<double> pars1 = pars;
  pars1.erase(pars1.end()-2,pars1.end());

  double RGWSedovRadius = pars[8]*mSol*pars[2]/(pars[1]*mSol/yr_to_sec)
                          /pc_to_cm;

  vector< vector<double> > dprofile = CreateDensityProfile(pars1);
  vector< vector<double> > massintegrand;
  vector< vector<double> > mass;
  vector< vector<double> > massinverse;
  for(unsigned int i=0;i<dprofile.size();i++) {
    double r = dprofile[i][0];
    double n = dprofile[i][1];
    fUtils->TwoDVectorPushBack( r, 4.* pi * m_p_g * n * r * r,massintegrand);
  }
  massintegrand = fUtils->SortTwoDVector(massintegrand,0);
  mass = fUtils->IntegratedProfile(massintegrand);
  for(unsigned int i=0;i<mass.size();i++)
    fUtils->TwoDVectorPushBack(mass[i][1],mass[i][0],massinverse);
  massinverse = fUtils->SortTwoDVector(massinverse,0);

  gsl_spline *m = fUtils->GSLsplineFromTwoDVector(massinverse);
  gsl_interp_accel *a = gsl_interp_accel_alloc();
  double RGWSedovRadiusFromProfile = fUtils->EvalSpline(pars[8]*mSol,m,a,
                                                        __func__,__LINE__);
  RGWSedovRadiusFromProfile /= pc_to_cm;
  cout << RGWSedovRadius << " " << RGWSedovRadiusFromProfile << endl;
  double RGVSR = 0.;
  (RGWSedovRadiusFromProfile > RGWSedovRadius) ? RGVSR = RGWSedovRadius:
    RGVSR = RGWSedovRadiusFromProfile;
  vector<double> pars2;
  pars2.push_back(pars[1]);
  pars2.push_back(pars[2]);
  pars2.push_back(pars[8]);
  pars2.push_back(pars[9]);
  pars2.push_back(RGVSR);

  CalculateThinShellApproximation(dprofile,pars[9],1.333,0.,0.,0);

  // make the r and v time profiles
  if(steps && tmin && tmax) {
    if(tmin>=tmax) {
      cout << "Astro::CalculateThinShellApproximation: Can't create profiles"
              "because tmin >= tmax (" << tmin <<" >= "<< tmax <<
              ") Exiting!" << endl;
      return;
    }
    double logtmin = log10(tmin);
    double logtmax = log10(tmax);
    double dlogt = (logtmax-logtmin)/steps;
    fUtils->Clear2DVector(forwardshockradiusprofile);
    fUtils->Clear2DVector(forwardshockvelocityprofile);
    vector<double> x;
    double t=0.;
    for(double logt = logtmin ; logt < logtmax ; logt += dlogt) {
      t = pow(10.,logt);
      x = ForwardShockInRGWind(t,pars2,dprofile);
      fUtils->TwoDVectorPushBack(t,x[0]/pc_to_cm,forwardshockradiusprofile);
      fUtils->TwoDVectorPushBack(t,x[1],forwardshockvelocityprofile);
    }
  }
  return;
}

/**
 * Free expansion phase forward shock dynamics in red giant wind.
 * Eq. 14 from Ptuskin&Zirakashvili 2005 (A&A 429, 755–765 (2005))
 */
 vector<double> Astro::ForwardShockInRGWind(double t, vector<double> pars,
                                            vector< vector<double> > profile){

  double tScaled = t/1.e3;
  double mDotRGWindScaled = pars[0]/1.e-5;
  double vRGWindScaled = pars[1]/1.e6;
  double mEjScaled = pars[2];
  double Escaled = pars[3]/1.e51;
  double RGWSedovRadius = pars[4];
  double k,R,V;
  k = pow(pow(Escaled,3.5)*vRGWindScaled/(mDotRGWindScaled*pow(mEjScaled,2.5)),1./8.);
  R = 7.7*k*pow(tScaled,7./8.);
  V = 6.6e3*k*pow(tScaled,-1./8.);
  R*=1.15;///<this connects the solution to the Thin Shell adiabatic solution
  V*=0.9;///<this connects the solution to the Thin Shell adiabatic solution
  //cout<<t<<" "<<R<<" "<<RGWSedovRadius<<std::endl;
  if(R>RGWSedovRadius) {
    vector<double> v = ThinShellRadiusAndSpeed(t);
    R = v[0]/pc_to_cm;
    V = v[1]/1.e5;
  }
  R *= pc_to_cm;
  V *= 1.e5;
  vector<double> v;
  v.push_back(R);
  v.push_back(V);
  return v;
}
