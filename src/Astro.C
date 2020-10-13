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

  xyzref.push_back(0.); // default origin of cartesian space = {0,0,0}
  xyzref.push_back(0.);
  xyzref.push_back(0.);


}

Astro::~Astro() {}


void Astro::SetGalacticReferencePoint(vector<double> xyz_ref) {
    xyzref = xyz_ref;
    return;
}

void Astro::TranslateToReferencePoint(double &x, double &y, double &z) {
    x += xyzref[0];
    y += xyzref[1];
    z += xyzref[2];
    return;
}


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

int Astro::BFieldGetSpiralArm(vector<double> xyz) {
  double r = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
  double phi = atan2(xyz[1],xyz[0]);
  double rx[] = {15.5, 5.1, 6.3, 7.1, 8.3, 9.8, 11.4, 12.7, 15.5};
  double opening_angle = pi * 11.5 / 180.;

  // Those represent multiple turns on the log spirals. Needed to find the right arm.
  double exp_fac[4];
  exp_fac[0] = exp(((phi-0.5*pi) - 2.*pi) / tan(pi/2. - opening_angle));
  exp_fac[1] = exp(((phi-0.5*pi) ) / tan(pi/2. - opening_angle));
  exp_fac[2] = exp(((phi-0.5*pi) + 2.*pi) / tan(pi/2. - opening_angle));
  exp_fac[3] = exp(((phi-0.5*pi) + 4.*pi) / tan(pi/2. - opening_angle));

  // get the right spiral arm
  double dr = 1e100;
  int spiral_arm = -1;
  if (r>5.) {
    for(int i=0;i<8;i++) {
      for(int j=0;j<4;j++) {
        double r1 =rx[i]*exp_fac[j];
        for(int k=0;k<4;k++) {
          double r2 =rx[i+1]*exp_fac[k];
          if (r1 <= r && r < r2 && fabs(r1-r2) < dr) {
            spiral_arm = i;
            dr = fabs(r1-r2);
          }
        }
      }
    }
  }
  return spiral_arm;
}


/**
 * regular B-Field disk component
 */
vector< double > Astro::BFieldDiskComponent(vector<double> xyz) {

  double r = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
  double z = xyz[2];
  double opening_angle = pi * 11.5 / 180.;

  double r_fac = 5. / r; // factor describing radial field strength decrease along spiral arms
  int spiral_arm = BFieldGetSpiralArm(xyz);

  double b_regular_i[] = {.1, 3.,-.9,-.8,-2.,-4.2,0., 2.7,.1};
  double b_regular_ring = 0.1;
  double h_disk = 0.5;
  double w_disk = 0.27;
  double transition = fUtils->LogisticsFunction(z,h_disk,w_disk);

  vector<double> RegularDisk; // B-field direction in reference coordinate system
  RegularDisk.resize(3,0.);

  double b_regular_disk = (r<5.) ? b_regular_ring: r_fac * b_regular_i[spiral_arm];
  b_regular_disk *= (1.-transition);
  if (r>3. && r < 5.) {
    RegularDisk[1]= b_regular_disk; 
  }
  if (r>=5. && r < 20.) {
    RegularDisk[0]= b_regular_disk * sin(opening_angle); 
    RegularDisk[1]= b_regular_disk * cos(opening_angle); 
  }

  return RotateBFieldComponent(RegularDisk,xyz);
}


/***************************************************************************
 * Regular toroidal halo component. Equation 6 in Jansson & Farrar 2012 (1)
 * https://ui.adsabs.harvard.edu/abs/2012ApJ...757...14J/abstract
 * Input:   - Position in the galaxy (x,y,z) in [kpc]
 * Output:  - B-field vector in [micro Gauss]
 **************************************************************************/
vector<double> Astro::BFieldRegularHaloComponent(vector<double> xyz) {

  double r = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
  double z = xyz[2];

  double b_n = 1.4;
  double b_s = -1.1;
  double r_n = 9.22;
  double r_s = 16.7;
  double w_h = 0.2;
  double z_0_tor = 5.3;
  double h_disk = 0.4;
  double w_disk = 0.27;

  double transition = fUtils->LogisticsFunction(z,h_disk,w_disk);
  double b_tor = exp(-fabs(z) / z_0_tor) * transition;
  b_tor *= (z>0) ? b_n * (1. - fUtils->LogisticsFunction(r,r_n,w_h)):
               b_s * (1. - fUtils->LogisticsFunction(r,r_s,w_h));

  vector<double> RegularHalo; // B-field direction in reference coordinate system
  RegularHalo.resize(3,0.);
  RegularHalo[1] = b_tor;

  return RotateBFieldComponent(RegularHalo,xyz);
}



/**
 * Regular X-shaped (when looking edge-on galaxy) component
 **/
vector<double> Astro::BFieldXComponent(vector<double> xyz) {

  double r = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
  double z = xyz[2];

  //#### X-field component
  double b_x_0 = 4.6;
  double theta_x_0 = 49. * pi / 180.;
  double r_x_c = 4.8;
  double r_x = 2.9;

  double r_p,theta,b_x;
  r_p = r * r_x_c / (r_x_c + fabs(z) / tan(theta_x_0));
  if(r_p < r_x_c) {
    theta = (r!=r_p) ? atan(fabs(z)/(r-r_p)): pi/2.;
    b_x = b_x_0 * exp(-r_p / r_x) * r_p*r_p / (r*r);
  }
  else {
    r_p = r - fabs(z)/tan(theta_x_0);
    theta = theta_x_0;
    b_x = b_x_0 * exp(-r_p / r_x) * r_p / r;//Intensity jumps between both domains ... Intended by the authors?
  }

  if (!z) theta = pi/2.; // field lines to cross the plane in perpendicular manner
  vector<double> RegularX; 
  RegularX.resize(3,0.);
  RegularX[0] = b_x * cos(theta); 
  RegularX[2] = b_x * sin(theta);
  if(z<0) RegularX[0] *= -1.;
  

  return RotateBFieldComponent(RegularX,xyz);
}

/***************************************************************************
 * Random, unordered B-field component from Jansson & Farrar 2012 (2),
 * https://ui.adsabs.harvard.edu/abs/2012ApJ...761L..11J/abstract
 * It returns a random vector.
 * Input:   - vector<double> containing the position
 * Output:  - Random B-field vector
 ***************************************************************************/
vector<double> Astro::BFieldRandomComponent(vector<double> xyz) {
  double r = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
  double z = xyz[2];
  double r_fac = 5. / r; // factor describing radial field strength decrease along spiral arms

  int spiral_arm = BFieldGetSpiralArm(xyz);

  double b_random_i[] = {10.81, 6.96,9.59,6.96,1.96,16.34,37.29,10.35};
  double b_random_int = 7.63;
  double z_0_disk = 0.61;
  double b_0 = 4.68;
  double r_0 = 10.97;
  double z_0 = 2.84;

  //####  random disk component rms strength
  double b_random_disk = (r<5.) ? b_random_int: r_fac * b_random_i[spiral_arm];
  if (r>20.) b_random_disk = 0.;
  b_random_disk *= fUtils->Gaussian1D(z,z_0_disk,0.,1.);
  b_random_disk = fUtils->GaussianRandom(b_random_disk,0,1)[0];

  //#### random halo component rms strength
  double b_random_halo = b_0 * exp(-r / r_0) * exp(-z*z / (2.*z_0*z_0) );
  b_random_halo = fUtils->GaussianRandom(b_random_halo,0,1)[0];

  //#### random sum rms strength
  double b_random = sqrt(b_random_disk*b_random_disk + b_random_halo*b_random_halo);
  b_random = fUtils->GaussianRandom(b_random,0,1)[0];
  
  //#### random direction
  vector<double> Random = fUtils->UniformRandom(-1.,1.,3);
  double norm = sqrt(Random[0]*Random[0]+
                     Random[1]*Random[1]+
                     Random[2]*Random[2]);

  Random[0] *= b_random / norm; 
  Random[1] *= b_random / norm; 
  Random[2] *= b_random / norm; 

  return Random;
}


/***************************************************************************
 * Random, unordered B-field component from Jansson & Farrar 2012 (2),
 * https://ui.adsabs.harvard.edu/abs/2012ApJ...761L..11J/abstract
 * It returns the rms-Strenght
 * Input:   - vector<double> containing the position
 * Output:  - rms strength of the random B-field
 ***************************************************************************/
double Astro::BFieldRandomRmsStrength(vector<double> xyz) {
  double r = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
  double z = xyz[2];
  double r_fac = 5. / r; // factor describing radial field strength decrease along spiral arms

  int spiral_arm = BFieldGetSpiralArm(xyz);

  double b_random_i[] = {10.81, 6.96,9.59,6.96,1.96,16.34,37.29,10.35};
  double b_random_int = 7.63;
  double z_0_disk = 0.61;
  double b_0 = 4.68;
  double r_0 = 10.97;
  double z_0 = 2.84;

  //####  random disk component rms strength
  double b_random_disk = (r<5.) ? b_random_int: r_fac * b_random_i[spiral_arm];
  if (r>20.) b_random_disk = 0.;
  b_random_disk *= fUtils->Gaussian1D(z,z_0_disk,0.,1.);

  //#### random halo component rms strength
  double b_random_halo = b_0 * exp(-r / r_0) * exp(-z*z / (2.*z_0*z_0) );

  //#### random sum rms strength
  double b_random = sqrt(b_random_disk*b_random_disk + b_random_halo*b_random_halo);
  

  return b_random;
}


/***********************************************************************
 * Striated random component (i.e. oriented along regular B-field) from
 * Jansson & Farrar 2012 (1) and (2).
 * Input:   - RegularDisk component
 *          - RegularHalo component
 *          - RegularX component
 * Output:  - Vector of the striated random B-field aligned with the
 *              regular B-field
 ***********************************************************************/
vector<double> Astro::BFieldStriatedComponent(vector<double> RegularDisk, 
                                              vector<double> RegularHalo,
                                              vector<double> RegularX) {

  // sum of all regular components
  vector<double> Regular;
  Regular.resize(3,0.);
  Regular[0] = RegularDisk[0] + RegularHalo[0] + RegularX[0];
  Regular[1] = RegularDisk[1] + RegularHalo[1] + RegularX[1];
  Regular[2] = RegularDisk[2] + RegularHalo[2] + RegularX[2];

  double beta = 1.36;
  double norm = sqrt(Regular[0]*Regular[0]+
                     Regular[1]*Regular[1]+
                     Regular[2]*Regular[2]);

  vector<double> Striated;
  Striated.resize(3);
  double striated_norm = fUtils->GaussianRandom(sqrt(beta) * norm,0,1)[0];

  Striated[0] = striated_norm / norm * Regular[0]; 
  Striated[1] = striated_norm / norm * Regular[1]; 
  Striated[2] = striated_norm / norm * Regular[2]; 

  return Striated;
}

/**
 * Rotate B-field Vector to the right position in the xy plane
 **/
vector< double > Astro::RotateBFieldComponent(vector<double> component, vector<double> xyz) {



  vector<double> plane,reference;
  plane.resize(3,0.);reference.resize(3,0.);
  plane[0] = xyz[0]; plane[1] = xyz[1];
  reference[0] = 1.;
  if(plane[0] < 0. && !plane[1]) {
    vector<double> plane_int;
    plane_int.resize(3,0.);
    plane_int[1] = -plane[0];
     
    component = fUtils->RotateVector(component,reference,plane_int);
    reference = plane_int;
  }

  component = fUtils->RotateVector(component,reference,plane);
  return component;
}


vector< vector<double> > Astro::CalculateBFieldComponents(vector<double> xyz, 
                                                          int component) {
  
  vector<double> reference; 
  reference.resize(3,0.);
  reference[0]=1.; 
    
  vector<double> RegularDisk,RegularHalo,RegularX,Random,Striated;
  RegularDisk.resize(3,0.); RegularHalo.resize(3,0.); RegularX.resize(3,0.);
  Random.resize(3,0.); Striated.resize(3,0.);
  

  // regular field components
  if(component==-1||!component||component==4) RegularDisk = BFieldDiskComponent(xyz);
  if(component==-1||component==1||component==4) RegularHalo = BFieldRegularHaloComponent(xyz);
  if(component==-1||component==2||component==4) RegularX = BFieldXComponent(xyz);
  if(component==-1||component==3) Random = BFieldRandomComponent(xyz);
  if(component==-1||component==4) Striated = BFieldStriatedComponent(RegularDisk,RegularHalo,RegularX);

  vector< vector<double> > components;
  components.push_back(RegularDisk);
  components.push_back(RegularHalo);
  components.push_back(RegularX);
  components.push_back(Striated);
  components.push_back(Random);

  return components;
}


vector< vector<double> > Astro::CalculateBFieldPositions(vector< vector<double> > pos, int component) {

    vector< vector<double> > bvecs;
    for(unsigned int i=0;i<pos.size();i++) {
        bvecs.push_back(CalculateBField(pos[i],component));
    }
    return bvecs;
}

vector<double> Astro::CalculateBField(vector<double> xyz, int component) {

  vector< vector<double> > components = CalculateBFieldComponents(xyz,component);
  vector<double> total;
  total.resize(3,0.);
  
  if(component>-1 && component<(int)components.size()) {
    total[0] += components[component][0];
    total[1] += components[component][1];
    total[2] += components[component][2];
  }
  if(component==-1) {
    for (unsigned int i=0; i<components.size(); i++) {
      total[0] += components[i][0];
      total[1] += components[i][1];
      total[2] += components[i][2];
    }
  }
  if(component>=(int)components.size() || component< -1) {
    cout<<"Astro::CalculateBField: Component "<<component<<" not available. Viable indices are "
        <<0<<" to "<<  components.size() - 1<< ". Returning empty vector."<<endl;
  }

  return total;

}




vector<double> Astro::CalculateBFieldStrengthComposition(vector<double> xyz) {
  double beta = 1.36;
  
  // Calculate the regular field first
  vector<double> RegularDisk,RegularHalo,RegularX, Regular;
  RegularDisk.resize(3,0.); RegularHalo.resize(3,0.); RegularX.resize(3,0.), Regular.resize(3, 0.);
    
  RegularDisk = BFieldDiskComponent(xyz);
  RegularHalo = BFieldRegularHaloComponent(xyz);
  RegularX = BFieldXComponent(xyz);    
 
  for (int i = 0; i < 3; i++) {
      Regular[i] = RegularDisk[i] + RegularHalo[i] + RegularX[i];
  }
  double norm = fUtils->Norm(Regular);
 
  double Random = BFieldRandomRmsStrength(xyz);
  double Striated = sqrt(beta)*norm;
  
  vector<double> temp;
  temp.push_back(norm); temp.push_back(Random); temp.push_back(Striated);
  return temp;
}


double Astro::CalculateBFieldStrength(vector<double> xyz) {
 vector<double> components = CalculateBFieldStrengthComposition(xyz); 
 double result = fUtils->Norm(components);
 return result;
}


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
//void Astro::CalculateBField(double x, double y, double z, double &B_tot, double &B_coh, double &B_ord, double &B_iso, vector<double> &regFieldDirection, vector<double> &isoFieldDirection) {
//  regFieldDirection.clear();
//  isoFieldDirection.clear();
//  double DistanceToClosestArm = 0.;
//  int ClosestArm = 0;
//  GetDistanceToNearestSpiralArm(x,y,DistanceToClosestArm,ClosestArm);
//  double r = sqrt(pow(x,2.)+pow(y,2.));
//  double R0 = 20.;
//  double R1 = 3.;
//  double R2 = 0.5;
//  double Rarms = 15.;
//  double RmolRing = 5.;
//  double Rmin = 3.;
//  double C0 = 1.;
//  double an = 1.;
//  double an1 = 1.26;
//  double an2 = 1.04;
//  double an3 = 0.61;
//  double an4 = 1.65;
//  double anMolRing  = 1.;
//  switch (ClosestArm){
//    case 4:
//      an = an4;
//      break;
//    case 3:
//      an = an3;
//      break;
//    case 2:
//      an = an2;
//      break;
//    case 1:
//      an = an1;
//      break;
//  }
//  double B_0 = 1.e-6;
//  double d0 = 0.3; //could be replaced by arm-width
//  double B_rms = 2.1e-6;
//  double ford = 1.9;
//  regFieldDirection.push_back(0.);
//  regFieldDirection.push_back(0.);
//  regFieldDirection.push_back(0.);

//  isoFieldDirection.push_back(0.);
//  isoFieldDirection.push_back(0.);
//  isoFieldDirection.push_back(0.);
//  if(r<Rmin||r>R0) {
//    B_coh = 0.;
//    B_iso = 0.;
//    B_ord = 0.;
//    B_tot = 0.;

//    return;

//  }
//  double rhoc = C0*exp(-DistanceToClosestArm*DistanceToClosestArm/(d0*d0))+1.;
//  if(r>Rarms) rhoc = 1.;
//  if(r>Rmin && r<RmolRing) {
//    rhoc  =  C0*exp(-pow(r-RmolRing,2.)/(d0*d0))+1.;///< own description, since exact profile is not described in paper
//    an = anMolRing;
//  }
//  /* coherent, large scale B-field along spiral arms */
//  B_coh = B_0*(1.-exp(-r*r/(R2*R2)));
//  B_coh *= exp(-r*r/(R0*R0))+exp(-r*r*r*r/(R1*R1*R1*R1));
//  B_coh *= rhoc*an;
//  /* small-scale irregular component of B-field */ //NOTE: This is only random in 2D! If 3D (incl. scatter perpendicular to the plane), the component along the arms would smaller)
//  B_iso = rhoc*B_rms;
//  /* irregular component along spiral arms */
//  B_ord = ford*(rhoc-1.)*B_rms;
//  /* dice relative orientation of irregular components to coherent one */
//  double phi_iso = 2.*pi*fUtils->Random();
//  /* Component along the coherent field lines */
//  double B_sumx = B_coh+B_iso*cos(phi_iso)+B_ord*cos(phi_iso);
//  /* Component perpendicular to the coherent field lines (only isotropic field) */
//  double B_sumy = B_iso*sin(phi_iso);

//  B_tot = sqrt(B_sumx*B_sumx+B_sumy*B_sumy);

//  double e_B_coh[3];
//  double theta = -(180./pi)*atan2(x,y)+90;
////  std::cout<<"(x,y)="<<x<<","<<y<<" "<<theta<<std::endl;
//  double xarmclose,yarmclose,rarmclose,xarmfar,yarmfar,rarmfar,x0,y0,x1,y1,rtest;
//  PositionOnSpiralArmAngular(theta, ClosestArm, xarmclose, yarmclose);
//  rarmclose = sqrt(xarmclose*xarmclose+yarmclose*yarmclose);
//  PositionOnSpiralArmAngular(theta+360., ClosestArm, xarmfar, yarmfar);
//  rarmfar = sqrt(xarmfar*xarmfar+yarmfar*yarmfar);
//  (fabs(r-rarmfar)>fabs(r-rarmclose))?(rtest=rarmclose):
//    rtest=rarmfar;
//  GalacticPositionXY(0.999*rtest,ClosestArm,x0,y0);
//  GalacticPositionXY(1.001*rtest,ClosestArm,x1,y1);
////  std::cout<<r<<" -> ("<<x0<<","<<y0<<")"<<" ("<<x1<<","<<y1<<")"<<std::endl;
//  double l = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
//  e_B_coh[0] = (x1-x0)/l;
//  e_B_coh[1] = (y1-y0)/l;
//  e_B_coh[2] = 0.;
//  if(r>Rarms) {
//    double weight = (r-Rarms)/(R0-Rarms);
//    double xrand = 1.-2.*fUtils->Random();
//    double yrand = 1.-2.*fUtils->Random();
//    double zrand = 1.-2.*fUtils->Random();

//    double rrand = sqrt(xrand*xrand+yrand*yrand+zrand*zrand);
//    xrand /= rrand;
//    yrand /= rrand;
//    zrand /= rrand;

//    e_B_coh[0] += weight*xrand;
//    e_B_coh[1] += weight*yrand;
//    e_B_coh[2] += weight*zrand;

//    double r_coh = sqrt(e_B_coh[0]*e_B_coh[0]+e_B_coh[1]*e_B_coh[1]+e_B_coh[2]*e_B_coh[2]);

//    e_B_coh[0] /= r_coh;
//    e_B_coh[1] /= r_coh;
//    e_B_coh[2] /= r_coh;
//  }
//  if(r>RmolRing) {
//    regFieldDirection[0] = e_B_coh[0];
//    regFieldDirection[1] = e_B_coh[1];
//    regFieldDirection[2] = e_B_coh[2];
//  }
////  std::cout<<"  -> "<<regFieldDirection[0]<<","<<regFieldDirection[1]<<","<<regFieldDirection[2]<<std::endl;

//  double e_B_iso[3];
//  e_B_iso[0] = e_B_coh[0];
//  e_B_iso[1] = e_B_coh[1];
//  e_B_iso[2] = e_B_coh[2];
//  RotateCoordinates(e_B_iso[0],e_B_iso[1],e_B_iso[2],0.,0.,phi_iso);
//  isoFieldDirection[0] = e_B_iso[0];
//  isoFieldDirection[1] = e_B_iso[1];
//  isoFieldDirection[2] = e_B_iso[2];

//  return;
//}

//double Astro::TotalBField(double x, double y, double z) {
//  double B_tot,B_coh,B_ord,B_iso;
//  vector<double> regFieldDirection, isoFieldDirection;
//  CalculateBField(x, y, z, B_tot, B_coh, B_ord, B_iso, regFieldDirection, isoFieldDirection);
//  std::cout<<x<<" "<<y<<" "<<z<<" "<<B_tot<<" "<<B_coh<<" "<<B_ord<<" "<<B_iso<<std::endl;
//  return B_tot;
//}

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
//  else {
    double norm = CalculateHINorm(R);
    double FWHM = GetHIFWHM(R)/1000.;
    n = 0.7*exp(-pow(z/(0.55*FWHM),2.));
    n += 0.19*exp(-pow(z/(1.38*FWHM),2.));
    n += 0.11*exp(-fabs(z)/(1.75*FWHM));
    n *= norm;
//  }
  if( R < 3.) {
    n += HICMZ(x,y,z);
    n += HIInnerDisk(x,y,z);
  }
  if(MODULATE && SPIRALARMMODEL != 2) n = ModulateGasDensityWithSpirals(n,x,y,z);
  return n;
}

/**
 *
 */
double Astro::H2Density(double x, double y, double z, bool MODULATE) {
//  std::cout<<x<<","<<y<<","<<z<<std::endl;
//  TranslateToReferencePoint(x,y,z);
//  std::cout<<"->"<<x<<","<<y<<","<<z<<std::endl;
//  std::cout<<"__________________"<<std::endl;
  double R = sqrt(pow(x,2.)+pow(y,2.));
  double n = 0.;
//  if( R < 3.) {
//    n = H2CMZ(x,y,z);
//    n += H2InnerDisk(x,y,z);
//    n *= 2.;
//  }
//  else {
    double norm = CalculateH2Norm(R);
    double sigma = GetH2FWHM(R)/(2.*sqrt(2.*log(2.)));
    sigma /= 1000.;
    n = norm*exp(-pow(z/(sqrt(2.)*sigma),2.));
//  }
  
  if( R < 3.) {
    n += 2.* H2CMZ(x,y,z);
    n += 2.*H2InnerDisk(x,y,z);
//    n *= 2.;
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
  return std::isnan(n) ? 0.: n;
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
//double Astro::nRadial(double *x, double *pars) {
//  vector<double> xyzobs;
//  xyzobs.resize(3);
//  vector<double> lbr;
//  lbr.resize(3);
//  xyzobs[0] = pars[2];
//  xyzobs[1] = pars[3];
//  xyzobs[2] = pars[4];
//  lbr[0] = pars[0];
//  lbr[1] = pars[1];
//  lbr[2] = x[0];
//  int gasspecies = pars[5];
//  int modulatewitharms = (int)pars[6];

//  vector<double> xyzactual = GetCartesian(lbr,xyzobs);
//  double xActual = xyzactual[0];
//  double yActual = xyzactual[1];
//  double zActual = xyzactual[2];

//  double n=0.;
//  if(!gasspecies) n = HIDensity(xActual,yActual,zActual);
//  else if(gasspecies==1) n = H2Density(xActual,yActual,zActual);
//  else {
//    std::cout<<"Astro::nRadial: Specify implemented gas species! returning 0.!"<<std::endl;
//    return 0.;
//  }
//  if(modulatewitharms==1) n = ModulateGasDensityWithSpirals(n,xActual,yActual,zActual);
//  return n;
//}




/************************************************************************************
 * Functions for the galactic density model from Ferriere 1998 for
 * the outer galaxy, and Ferriere et al. 2007 for R < 3 kpc.
 *********************************************************************************/


/***********************************************************
 * Function to return the total H density according to
 * Ferriere 1998 and 2997.
 * Input:   - Position in [kpc]
 * Output:  - Number density in [cm^-3]
 * ********************************************************/
double Astro::GetTotalHDensityFerriere(vector< double > xyz){
    return (GetH2DensityFerriere(xyz) + GetCMDensityFerriere(xyz) + GetWNMDensityFerriere(xyz)
            + GetWIMDensityFerriere(xyz) + GetHIMDensityFerriere(xyz));
}

/***********************************************************
 * Function to return the density of molecular hydrogen,
 * equation 2 and 3 in Ferriere 1998.
 * Input:   - Position in [kpc]
 * Output:  - Number density in [cm^-3]
 * ********************************************************/
double Astro::GetH2DensityFerriere(vector< double > xyz){
    double z = xyz[2];
    double R = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
    double R_sun = 8.5;     // position of the sun in [kpc]
    
    double Hm = 81.0* pow(R/R_sun, 0.58) * 1.0e-3;  // in kpc 
    
    double result = 0.58 * exp(-((R-4.5)*(R-4.5) - (R_sun - 4.5)*(R_sun - 4.5))/(2.9*2.9))
                    * pow(R/R_sun, 0.58) * exp(-(z*z/(Hm*Hm)));
    return result;
}

/****************************************************************
 * Function to return the density of the neutral cold Hydrogen,
 * equation 6 in Ferriere 1998.
 * Input:   - Position in [kpc]
 * Output:  - Number density in [cm^-3]
 * **************************************************************/  
double Astro::GetCMDensityFerriere(vector< double > xyz){
    double z = xyz[2]*1.0e3;    // Now in [pc]
    double R = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
    double R_sun = 8.5;
    
    double alpha;
    if ((3.5 <= R) && (R <= R_sun)) alpha = 1.0;
    else alpha = R/R_sun;

    double H1 = 127.0*alpha;
    double H2 = 318.0*alpha;
    double H3 = 403.0*alpha;
    
    double nc = 0.34/(alpha*alpha) * (0.859*exp(-z*z/(H1*H1)) + 0.047*exp(-z*z/(H2*H2)) + 0.094*exp(-abs(z)/H3));
    
    return nc;
}

/***********************************************************************
 * Function to return the density of the warm neutral Hydrogen medium,
 * equation 7 in Ferriere 1998.
 * Input:   - Position in [kpc]
 * Output:  - Number density in [cm^-3]
 * ********************************************************************/  
double Astro::GetWNMDensityFerriere(vector< double > xyz){
    double z = xyz[2]*1.0e3;    // Now in [pc]
    double R = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
    double R_sun = 8.5;
    
    double alpha;
    if ((3.5 <= R) && (R <= R_sun)) alpha = 1.0;
    else alpha = R/R_sun;

    double H1 = 127.0*alpha;
    double H2 = 318.0*alpha;
    double H3 = 403.0*alpha;
    
    double nc = 0.226/alpha* ((1.745 - 1.289/alpha)*exp(-z*z/(H1*H1)) 
                                    + (0.473 - 0.07/alpha)*exp(-z*z/(H2*H2)) 
                                    + (0.283 - 0.142/alpha)*exp(-abs(z)/H3));
    
    return nc;
}

/***********************************************************************
 * Function to return the density of the warm ionized Hydrogen medium,
 * equation 9 in Ferriere 1998.
 * Input:   - Position in [kpc]
 * Output:  - Number density in [cm^-3]
 * ********************************************************************/  
double Astro::GetWIMDensityFerriere(vector< double > xyz){
    double z = xyz[2]*1.0e3;    // Now in [pc]
    double R = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
    double R_sun = 8.5;  // in [kpc]
    
    double result = 0.0237* exp(-(R*R - R_sun*R_sun)/1369.0) * exp(-abs(z)/1.0e3)
                    + 0.0013 * exp(-((R - 4.0)*(R - 4.0) - (R_sun - 4.0)*(R_sun - 4.0))/4.0)
                    * exp(-abs(z)/150.0);
    return result;
}

/***********************************************************************
 * Function to return the density of the hot medium hydrogen,
 * equation 12 and 52 in Ferriere 1998.
 * Input:   - Position in [kpc]
 * Output:  - Number density in [cm^-3]
 * ********************************************************************/  
double Astro::GetHIMDensityFerriere(vector< double > xyz){
    double z = xyz[2];
    double R = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
    double R_sun = 8.5;  // in [kpc]
    
    double Hh = 1.5 * pow(R/R_sun, 1.65);   // Equ. 52 in Ferriere
    
    double result = 4.8e-4*(0.12 * exp(-(R-R_sun)/4.9)
                    + 0.88*exp(-((R - 4.5)*(R - 4.5) - (R_sun - 4.5)*(R_sun - 4.5))/8.41))
                    * pow(R/R_sun, -1.65) * exp(-abs(z)/Hh);
    return result;
}








/*    ------            Galactic Structure stuff              ------    */

/**
 * Galactic coordinates (R,GL,GB) ->Cartesian
 */

/**
 * In there for backwards compability
 */
vector<double> Astro::GetCartesian(double r, double l, double b, vector<double> xyzobs) {
    vector<double> lbr;
    lbr.push_back(l);
    lbr.push_back(b);
    lbr.push_back(r);
    return GetCartesian(lbr,xyzobs);
}
vector<double> Astro::GetCartesian(vector<double> lbr,vector<double> xyzobs) {


  double l = pi * lbr[0] / 180.;
  double b = pi * lbr[1] / 180.;
  double r = lbr[2];

  if(!xyzobs[0] && !xyzobs[1]) {
    cout<<"Astro::GetCartesian: l and b at a observer's position directly above or below "<<
          "the Galactic center (x,y components = 0) are not well defined. Please "<<
          "provide x and/or y values to fix orientation of Galactic coordinates. Returning original vector." << endl;
    return lbr;
  }


  vector<double> v;
  v.push_back(-r*cos(l)*cos(b)); //x component
  v.push_back(r*sin(l)*cos(b)); //y component
  v.push_back(r*sin(b)); //z component

  // unit vector along x. This is the right direction for the xyz components of v
  vector<double> p1;    
  p1.push_back(1.);    
  p1.push_back(0.);    
  p1.push_back(0.);

  vector<double> xyzorigin;    
  xyzorigin.push_back(xyzobs[0]);    
  xyzorigin.push_back(xyzobs[1]);    
  xyzorigin.push_back(0.);
  
  // rotate coordinates so that the (x,0,0)-axis lies along xyzorigin. This defines
  // the orientation of the system, relative to which the coordinates will be transformed. 
  // For instance, the earth would at xyzorigin = (0,8.5,0).
  // This re-defines the coordinate system to point towards the earth's position.

  // put intermediate step if initial rotation is by 180 degree (for technical reasons)
  if(xyzorigin[0]<0 && !xyzorigin[1]) {
    vector<double> xyz_temp;    
    xyz_temp.push_back(xyzorigin[0]);    
    xyz_temp.push_back(fabs(xyzorigin[0]));    
    xyz_temp.push_back(0.);
    v = fUtils->RotateVector(v,p1,xyz_temp);
    p1 = xyz_temp;
  }
  vector<double> v_rel = fUtils->RotateVector(v,p1,xyzorigin);

  // rotate v_rel to the final position in 3D space at the coordinates xyzobs
  vector<double> v_rot = fUtils->RotateVector(v_rel,xyzorigin,xyzobs);
  vector<double> xyz;
  xyz.resize(3);
  xyz[0] = xyzobs[0] + v_rot[0];
  xyz[1] = xyzobs[1] + v_rot[1];
  xyz[2] = xyzobs[2] + v_rot[2];



//  double dl, db;
//  vector<double> xyz;
//  xyz.resize(3);
//  double xobs = xyzobs[0];
//  double yobs = xyzobs[1];
//  double zobs = xyzobs[2];

//  dl = db = 0.;
//  dl = atan2(yobs,xobs);
//  db = atan(zobs/sqrt(xobs*xobs+yobs*yobs));
//  std::cout<<dl/(pi/180.)<<" "<<db/(pi/180.)<<std::endl;
//  double phi = (pi/180.)*l+dl;
//  double theta = pi/2. - (pi/180.)*b+db;
//  theta = pi/2. - (db + (pi/180.)*b);
//  std::cout<<theta/(pi/180.)<<" "<<phi/(pi/180.)<<std::endl;
//  std::cout<<".............."<<std::endl;
//  xyz[0] = xobs - r*sin(theta)*cos(phi);
//  xyz[1] = yobs - r*sin(theta)*sin(phi);
//  xyz[2] = zobs + r*cos(theta);


 // neuer plan : l,b ->dxdydz -> rotation -> dxdydz' -> addition zu vobs
  return xyz;
}

/**
 * Galactic Coordinates(GL,GB,R) coordinates ->  Cartesian for vector of points
 */
vector< vector<double> > Astro::GetCartesianPositions(
                        vector< vector<double> > lbr, vector<double> xyzobs) {
    vector< vector<double> > temp;
    if (!lbr.size()) {
        cout << "Astro::GetCartesianPositions: Input vector of positions is "
                "empty. Returning empty vector." << endl;
        return temp;
    }

    for(unsigned int i=0; i<lbr.size(); i++) {
        vector<double> xyz = GetCartesian(lbr[i], xyzobs);
        temp.push_back(vector<double>());
        temp[temp.size()-1].push_back(xyz[0]);
        temp[temp.size()-1].push_back(xyz[1]);
        temp[temp.size()-1].push_back(xyz[2]);
    }
    return temp;
}
/**
 * In there for backwards compatibility
 */
void Astro::GetGalactic(double x, double y, double z, double xobs, double yobs, double zobs, double &l, double &b) {
    vector<double> xyz;
    xyz.push_back(x);
    xyz.push_back(y);
    xyz.push_back(z);
    vector<double> xyzobs;
    xyzobs.push_back(xobs);
    xyzobs.push_back(yobs);
    xyzobs.push_back(zobs);
    vector<double> gal = GetGalactic(xyz, xyzobs);
    l = gal[0]; b = gal[1];
    return;

}



/**
 * Cartesian coordinates -> Galactic Coordinates(GL,GB,R)
 */
vector<double> Astro::GetGalactic(vector<double> xyz, vector<double> xyzobs) {

  double x = xyz[0];
  double y = xyz[1];
  double z = xyz[2];

  double l,b,r;
  vector<double> lbr;
  double xobs = xyzobs[0];
  double yobs = xyzobs[1];
  double zobs = xyzobs[2];
  lbr.resize(3);

  double dx  = x-xobs;
  double dy  = y-yobs;
  double dz =  z-zobs;

  r = sqrt(dx*dx+dy*dy+dz*dz);
  double rrobs = sqrt(xobs*xobs+yobs*yobs);
  double dr = sqrt(dx*dx+dy*dy);
  double prod = -dx*xobs - dy*yobs;

  double vec_prod = dx*yobs - dy*xobs;
  l = 180.*acos(prod/(dr*rrobs))/pi;
  if(vec_prod<=0.) l*=-1.;
  b = 180.*asin(dz/r)/pi - 180.*asin(zobs/rrobs)/pi;
  lbr[0] = l;
  lbr[1] = b;
  lbr[2] = r;

  return lbr;
}

/**
 * Cartesian coordinates -> Galactic Coordinates(GL,GB,R) for vector of points
 */
vector< vector<double> > Astro::GetGalacticPositions(
                        vector< vector<double> > xyz, vector<double> xyzobs) {
    vector< vector<double> > temp;
    if (!xyz.size()) {
        cout << "Astro::GetGalacticPositions: Input vector of positions is "
                "empty. Returning empty vector." << endl;
        return temp;
    }

    for(unsigned int i=0; i<xyz.size(); i++) {
        vector<double> lbr = GetGalactic(xyz[i], xyzobs);
        temp.push_back(vector<double>());
        temp[temp.size()-1].push_back(lbr[0]);
        temp[temp.size()-1].push_back(lbr[1]);
        temp[temp.size()-1].push_back(lbr[2]);
    }
    return temp;
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
 * Switch on an individual arm in Astro spiral model
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
  else if(!surfacedensitymodel.compare("Flat")) SURFACEDENSITYMODEL = 2;
  else {
    std::cout<<"Astro::SetSurfaceDensityModel: Specify supported surface density model (currently 'CaseBhattacharya', 'IusifovKucuk' and 'Flat')! returning 0 value."<<std::endl;
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


double Astro::ApparentSize(double size, vector<double> lbr) {
    if(lbr.size() != 3) {
        cout << "Astro::ApparentSize: Something's wrong with "
                "position vector. Returning 0." << endl;
        return 0.;
    }
    return 180.*atan(size * 1e-3 / lbr[2]) / pi;
}

vector<double> Astro::GetApparentSizes(vector<double> sizes, vector< vector<double> > lbr) {
    vector<double> angles;
    if(sizes.size() != lbr.size()) {
        cout << "Astro::GetApparentSizes: Size and position vectors are of different "
                "length. Returning empty vector." << endl;
        return angles;
    }

    for(unsigned int i=0;i<sizes.size();i++)
        angles.push_back(ApparentSize(sizes[i],lbr[i]));

    return angles;
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
  if(SURFACEDENSITYMODEL==2) return fUtils->LinearRandom(1,rmin,rmax,n);
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
  return fUtils->CustomFunctionRandom(v,n);
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
      if(!x[0] || !x[1]) continue;
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
      cout << "Astro::CalculateTrueloveMcKeeSolution: Can't create profiles"
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
  if(!index)          { static const double itl[] = {  1.1,  0., 0.343, 0.0961, 0.495, 0.727, 0.545, 0.585, 0.106,  0.,    0.,    0.,    0.,    0.}; truelovemckeeparams = vector<double>(itl, itl + sizeof(itl)/sizeof(double)); }
  else if(index== 2.) { static const double itl[] = {  1.1,  0., 0.343, 0.0947, 0.387, 0.679, 0.503, 0.686,-0.151, 1.6,    0.,    0.,    0.,    0.}; truelovemckeeparams = vector<double>(itl, itl + sizeof(itl)/sizeof(double)); }
  else if(index== 4.) { static const double itl[] = {  1.1, 0.1, 0.343, 0.0791, 0.232, 0.587,    0.,    0.,    0.,  0.,   1.7,    0.,    0.,    0.}; truelovemckeeparams = vector<double>(itl, itl + sizeof(itl)/sizeof(double)); }
  else if(index== 6.) { static const double itl[] = { 1.39,  0.,  0.39,     0.,  1.04,  1.07,    0.,    0.,    0.,  0., 0.513, 0.541, 0.527, 0.112}; truelovemckeeparams = vector<double>(itl, itl + sizeof(itl)/sizeof(double)); }
  else if(index== 7.) { static const double itl[] = { 1.26,  0.,  0.47,     0., 0.732, 0.881,    0.,    0.,    0.,  0., 0.363, 0.469, 0.553, 0.116}; truelovemckeeparams = vector<double>(itl, itl + sizeof(itl)/sizeof(double)); }
  else if(index== 8.) { static const double itl[] = { 1.21,  0.,  0.52,     0., 0.605, 0.788,    0.,    0.,    0.,  0., 0.292, 0.413, 0.530, 0.139}; truelovemckeeparams = vector<double>(itl, itl + sizeof(itl)/sizeof(double)); }
  else if(index== 9.) { static const double itl[] = { 1.19,  0.,  0.55,     0., 0.523, 0.725,    0.,    0.,    0.,  0., 0.249, 0.371, 0.497, 0.162}; truelovemckeeparams = vector<double>(itl, itl + sizeof(itl)/sizeof(double)); }
  else if(index==10.) { static const double itl[] = { 1.17,  0.,  0.57,     0., 0.481, 0.687,    0.,    0.,    0.,  0., 0.220, 0.340, 0.463, 0.192}; truelovemckeeparams = vector<double>(itl, itl + sizeof(itl)/sizeof(double)); }
  else if(index==12.) { static const double itl[] = { 1.15,  0.,  0.60,     0., 0.424, 0.636,    0.,    0.,    0.,  0., 0.182, 0.293, 0.403, 0.251}; truelovemckeeparams = vector<double>(itl, itl + sizeof(itl)/sizeof(double)); }
  else if(index==14.) { static const double itl[] = { 1.14,  0.,  0.62,     0., 0.389, 0.603,    0.,    0.,    0.,  0., 0.157, 0.259, 0.354, 0.277}; truelovemckeeparams = vector<double>(itl, itl + sizeof(itl)/sizeof(double)); }
  else {
    cout << "Astro::CalculateTrueLoveMcKeeParams: Density Index (index="
         << index << ") not supported! Returning nothing." << endl;
    truelovemckeeparams= vector<double>();
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
    string ParNamesArr[] = {"l_ED","w_core","phi_ED","phi_EDeff","t*_ST",
                               "R*_ST","R*_r_ST","v~*_r_ST","a~*_r_ST","t*_ref",
                               "t*_core","r*_r_core","v~*_r_core","a~*_r_core",
                               "f_n","alpha","eta","charradiustruelovemckee",
                               "chartimetruelovemckee",
                               "charvelocitytruelovemckee","sedovtaylortime",
                               "coreexittime","index"};
    vector<string> ParNames( ParNamesArr, ParNamesArr + ( sizeof ( ParNamesArr ) /  sizeof ( std::string ) ) );
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
      cout << "Astro::CalculateForwardShockInRGWind: Can't create profiles"
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

void Astro::CalculateTangChevalierSolution(vector<double> pars, double tmin,
                                           double tmax, int steps) {
    
  fUtils->Clear2DVector(forwardshockradiusprofile);
  fUtils->Clear2DVector(forwardshockvelocityprofile);
  fUtils->Clear2DVector(cdradiusprofile);
  fUtils->Clear2DVector(cdvelocityprofile);

  if(pars.size()!=5) {
    cout << "Astro::CalculateTangChevalierSolution: wrong number of parameters ("
         << pars.size() << "). Parameter list takes exactly 4 "
            "parameters: " << endl;
    cout << " - [0] eta parameter (cm^-3) " << endl;
    cout << " - [1] SN blast energy (erg) " << endl;
    cout << " - [2] SN ejecta mass (mSol) " << endl;
    cout << " - [3] ejecta density index " << endl;
    cout << " - [4] ambient density index " << endl;
    cout << "Exiting!" <<endl;
    return;
  }

  if (pars[3] < 0 || pars[3] == 3 || pars[3] == 5 || pars[3] == 11 || pars[3] ==13 || pars[3]>14){
    cout << "Astro::CalculateTangChevalierSolution: ejecta density index of n = "
         << pars[3] << " not supported. Exiting." <<endl;
    return;
  }

  if (!(pars[4] == 0 || pars[4] == 2)){
    cout << "Astro::CalculateTangChevalierSolution: ambient density index of s = "
         << pars[4] << " not supported. Exiting." <<endl;
    return;
  }

  if (!tmin) { 
    cout << "Astro::CalculateTangChevalierSolution: tmin has to be larger than 0."
            " Exiting." <<endl;
        return;
  }

  /* intercept invalid values here */

  double alphas_i0[] = {1.25,1.19,1.10,0.,0.80,0.,36.1,20.3,14.7,10.4,8.91,0.,7.11,0.,6.23};
  double alphas_i2[] = {0.95,0.91,0.85,0.,0.63,0.,11.3,8.00,6.22,5.16,4.56,0.,3.81,0.,3.40};
  double alphas_cd_i0[] = {0.89,0.94,1.08,0.,1.38,0.,5.54,4.28,3.09,3.01,2.92,0.,2.35,0.,2.20};
  double alphas_cd_i2[] = {1.23,1.16,1.06,0.,0.58,0.,6.04,4.61,3.81,3.35,2.94,0.,2.47,0.,2.21};
  double zetas_b_i0[] = {0.,0.,0.,0.,0.,0.,1.06,1.06,1.08,1.12,1.15,0.,1.21,0.,1.26};
  double zetas_b_i2[] = {0.,0.,0.,0.,0.,0.,0.77,0.83,0.90,0.97,1.03,0.,1.14,0.,1.23};
  double zetas_c_i0[] = {0.,0.,0.,0.,0.,0.,0.84,0.89,0.94,0.98,1.01,0.,1.08,0.,1.13};
  double zetas_c_i2[] = {0.,0.,0.,0.,0.,0.,0.56,0.64,0.71,0.77,0.83,0.,0.93,0.,1.01};
  double w_core_i0[] = {0.,0.,0.,0.,0.05,0.,0.05,0.05,0.05,0.05,0.05,0.,0.05,0.,0.05};
  double w_core_i2[] = {0.,0.,0.,0.,0.1,0.,0.1,0.1,0.1,0.1,0.1,0.,0.1,0.,0.1};
  double b_i0[] = {-0.37,-0.25,-0.06,0.,0.27,0.,-0.1,-0.11,-0.16,-0.13,-0.12,0.,-0.16,0.,-0.16};
  double b_i2[] = {0.50,0.50,0.51,0.,0.47,0.,0.53,0.53,0.53,0.53,0.53,0.,0.52,0.,0.52};
  double c_i0[] = {1.66,1.42,1.05,0.,0.72,0.,1.11,1.09,1.16,1.09,1.07,0.,1.11,0.,1.10};
  double c_i2[] = {0.53,0.54,0.52,0.,0.71,0.,0.47,0.47,0.47,0.47,0.47,0.,0.48,0.,0.48};

  double Esn = pars[1];
  double Mej = pars[2]*mSol;
  double index_e = pars[3];
  double index_a = pars[4];
  double eta = pars[0] * m_p_g;

  double rChar = pow(Mej,1./(3.-index_a)) * pow(eta,-1./(3.-index_a));
  double tChar = pow(Esn,-0.5) * pow(Mej,(5.-index_a)/(2.*(3.-index_a))) 
                   * pow(eta,-1./(3.-index_a));
  double fac = 1.05;
  double w_core = !index_a ? w_core_i0[int(index_e)] : w_core_i2[int(index_e)];
  double t_star=0.;double r_star = 0.;double r = 0.;double v = 0.;
  double r_c_star=0.; double r0=0.;double t0=0.;
  double r_c=0.;double v_c=0.; double r_c0 = 0.;
  double alpha = !index_a ? alphas_i0[int(index_e)] : alphas_i2[int(index_e)];
  double alpha_cd = !index_a ? alphas_cd_i0[int(index_e)] : alphas_cd_i2[int(index_e)];
  double xi = !index_a ? 2.026 : 3./(2.*pi);
  double q_b = !index_a ? 1.1 : 1.19;
  double lambda_c = 0.;
  if (index_e > 3) {
      lambda_c = 2.*pow(w_core,-2.)*((5.-index_e)/(3.-index_e));
      lambda_c *= (pow(w_core,index_e-3.)-index_e/3.) / (pow(w_core,index_e-5.)-index_e/5.);
      lambda_c = sqrt(lambda_c);
  }
  else lambda_c = sqrt(2.*(5.-index_e)/(3.-index_e));
  double lambda_b = index_e<5. ? q_b*lambda_c : 0.;
  double zeta_b = !index_a ? zetas_b_i0[int(index_e)] : zetas_b_i2[int(index_e)];
  double zeta_c = !index_a ? zetas_c_i0[int(index_e)] : zetas_c_i2[int(index_e)];
  double b = !index_a ? b_i0[int(index_e)] : b_i2[int(index_e)];
  double c = !index_a ? c_i0[int(index_e)] : c_i2[int(index_e)];

  cout<<">>Parameters for Tang & Chevalier's solution:"<<endl;
  cout<<">>>Input parameters:"<<endl;
  cout<<"   s = "<<index_a<<endl;
  cout<<"   n = "<<index_e<<endl;
  cout<<"   E_sn = "<<Esn<<" erg"<<endl;
  cout<<"   eta = "<<eta/m_p_g<<"/cm^3"<<endl;
  cout<<"   Mej = "<<Mej/mSol<<" mSol"<<endl;
  cout<<"   w_core = "<<w_core<<endl;
  cout<<">>>Internal model parameters - forward shock:"<<endl;
  cout<<"   alpha = "<<alpha<< endl;
  if (lambda_b) cout<< "   lambda_b = "<<lambda_b<< endl;
  else cout<< "   lambda_b = -"<< endl;
  cout<< "   xi = "<<xi<< endl;
  if (zeta_b) cout<< "   zeta_b = "<<zeta_b<< endl;
  else cout<< "   zeta_b = -"<< endl;
  cout<<">>>Internal model parameters - contact discontinuity:"<<endl;
  cout<<"   alpha = "<<alpha_cd<< endl;
  cout<<"   b = "<<b<<endl;
  cout<<"   c = "<<c<<endl;
  if (zeta_c) cout<< "   zeta_c = "<<zeta_c<< endl;
  else cout<<"   zeta_c = -"<< endl;
  cout<<"   lambda_c = "<<lambda_c<< endl;
  cout<<">>>Characteristic scales:"<<endl;
  cout<<"   t_char = "<<tChar/yr_to_sec<<" yrs"<<endl;
  cout<<"   r_char = "<<rChar/pc_to_cm<<" pc"<<endl;  
  cout<<"Please note that this solution does not include the reverse shock dynamics."<<endl;    
  cout<<"For that, please use the Truelove & McKee model." << endl;

  if (index_e < 5.) {
    for (double t=tmin/fac;t<tmax;t*=fac) {
      t_star = t * yr_to_sec/ tChar;

      r_star = pow(lambda_b*t_star,-2.*alpha);
      r_star += pow(xi*t_star*t_star,-2.*alpha/(5.-index_a));
      r_star = pow(r_star,-1./(2.*alpha));
          
      r_c_star = pow(lambda_c*t_star,-alpha_cd);
      r_c_star += pow(c*pow(t_star,b),-alpha_cd);
      r_c_star = pow(r_c_star,-1./alpha_cd);

      r = r_star*rChar/pc_to_cm;
      r_c = r_c_star*rChar/pc_to_cm;

      if (t0 && r0) v = pc_to_cm * (r-r0) / (t-t0) / yr_to_sec;
      if (t0 && r_c0) v_c = pc_to_cm * (r_c-r_c0) / (t-t0) / yr_to_sec;
      t0 = t;
      r0 = r;
      r_c0 = r_c;
      if(t>=tmin){
        fUtils->TwoDVectorPushBack(t,r,forwardshockradiusprofile);
        fUtils->TwoDVectorPushBack(t,v,forwardshockvelocityprofile);
        fUtils->TwoDVectorPushBack(t,r_c,cdradiusprofile);
        fUtils->TwoDVectorPushBack(t,v_c,cdvelocityprofile);
      }
    }
  }
  else if (index_e > 5.){
    for (double t=tmin/fac;t<tmax;t*=fac) {
      t_star = t * yr_to_sec / tChar;

      r_star = pow(zeta_b*pow(t_star,(index_e-3.)/(index_e-index_a)),-alpha);
      r_star += pow(xi*t_star*t_star,-alpha/(5.-index_a));
      r_star = pow(r_star,-1./alpha);

      r_c_star = pow(zeta_c*pow(t_star,(index_e-3.)/(index_e-index_a)),-alpha_cd);
      r_c_star += pow(c*pow(t_star,b),-alpha_cd);
      r_c_star = pow(r_c_star,-1./alpha_cd);

      r = r_star*rChar/pc_to_cm;
      r_c = r_c_star*rChar/pc_to_cm;

      if (t0 && r0) v = pc_to_cm * (r-r0) / (t-t0) / yr_to_sec;
      if (t0 && r_c0) v_c = pc_to_cm * (r_c-r_c0) / (t-t0) / yr_to_sec;
      t0 = t;
      r0 = r;
      r_c0 = r_c;
      if(t>=tmin){
        fUtils->TwoDVectorPushBack(t,r,forwardshockradiusprofile);
        fUtils->TwoDVectorPushBack(t,v,forwardshockvelocityprofile);
        fUtils->TwoDVectorPushBack(t,r_c,cdradiusprofile);
        fUtils->TwoDVectorPushBack(t,v_c,cdvelocityprofile);
      }
    }
  }  

  return;
}

vector< vector<double> > Astro::LineOfSight(double gl, double gb, vector<double> obs, double r_max, int steps) {
  double dr = r_max / steps;
  vector<double> rvals;
  for (double r = 0; r < r_max; r+= dr) {
    rvals.push_back(r);
  }
  return LineOfSight(gl, gb, obs, rvals);
}


vector< vector<double> > Astro::LineOfSight(double gl, double gb, vector<double> obs, vector<double> rvals) {

  vector< vector<double> > los;
  vector<double> xyz,lbr;
  lbr.resize(3);
  lbr[0] = gl; lbr[1] = gb;
  
  for (unsigned int i = 0; i < rvals.size(); i++) {
    lbr[2] = rvals[i];
    xyz = GetCartesian(lbr,obs);
    fUtils->TwoDVectorPushBack(xyz[0],xyz[1],xyz[2],los);
  }

  return los;
}


