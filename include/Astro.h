#ifndef _ASTRO_
#define _ASTRO_

#include "Utils.h"

/* TeV->erg */
#define TeV_to_erg 1.602
/* parsec to cm */
#define pc_to_cm 3.0857e18
/* kiloparsec to cm */
#define kpc_to_cm 3.0857e21
/* proton mass in g */
#define m_p_g 1.6726e-24
/* year in seconds */
#define yr_to_sec 3.15576e7
/* solar mass */
#define mSol 1.9891e33
/* Thomson cross section */
#define sigma_T 6.6524e-25
/* electron mass in erg */
#define m_e 8.187e-7
/* boltzmann constant (erg/K) */
#define kb 1.380658e-16
/* proton mass in erg */
#define m_p 1.50310854e-3
/* pi0 mass in erg */
#define m_pi 2.1622194e-4
/* parsec to cm */
#define pc_to_cm 3.0857e18
/* year in seconds */
#define yr_to_sec 3.15576e7
/* solar mass */
#define mSol 1.9891e33
/* speed of light cm/s */
#define c_speed 29979245800.
/* elementary charge */
#define el_charge 4.80320427e-10
/* classical electron radius (cm) */
#define eRadius 2.8179e-13
/* planck's constant */
#define hp 6.62606896e-27
/* fine structure constant */
#define fineStructConst 7.2974e-3
/* hour in seconds */
#define h_to_sec 3.6e3
/* pc to lyr */
#define pc_to_lyr 3.26156

using namespace std;

/**
 * @file Astro.C
 * @author Joachim Hahn
 * @short Class that provides utility functions
 * Provides different things, mainly useful for population studies.
 * This includes spiral galaxy models as well as gas distributions,
 * B-Field distributions, coordinate transformations etc.
 * Also includes more technical functions that are helpful in writing
 * GAMERA programs.
 */

class Astro  {
  private:
    Utils *fUtils;
    vector<double> xyzref;
    void TranslateToReferencePoint(double &x, double &y, double &z);
    double x_a1[7];
    double y_a1[7];
    double x_a2[7];
    double y_a2[7];
    double x_a3[7];
    double y_a3[7];
    double x_a4[7];
    double y_a4[7];
    vector<gsl_spline*> TaylorCordesArms,TaylorCordesArmsInv;
    vector<double> psi0V;
    vector<double> theta0V;
    vector<double> r0V;
    vector<int> ArmsVector;
    double barRadius;
    double barAngle;
    double armWidth;
    double scaleHeight;
    bool QUIETMODE;
    double rMax;
    double rMinInternalBoundary;
    double rMaxInternalBoundary;
    double tMinInternalBoundary;
    double tMaxInternalBoundary;
    int SPIRALARMMODEL;
    int SURFACEDENSITYMODEL;
    int CENTRALSTRUCTUREMODEL;
    int MSBUBBLEMODEL;
    gsl_spline *TaylorCordesArm1,*TaylorCordesArm2,*TaylorCordesArm3,
               *TaylorCordesArm4,*TaylorCordesArm1Inv,*TaylorCordesArm2Inv,
               *TaylorCordesArm3Inv,*TaylorCordesArm4Inv,
               *thinshellvelocitylookup,*thinshellradiuslookup,
               *r_edreverselookup,*v_edreverselookup,*r_edforwardlookup,
               *v_edforwardlookup,*r_streverselookup,*v_streverselookup,
               *r_stforwardlookup,*v_stforwardlookup;
    gsl_interp_accel *accArm1,*accArm2,*accArm3,*accArm4,
                     *accArm1Inv,*accArm2Inv,*accArm3Inv,*accArm4Inv,
                     *atsrad,*atsvel,*aredr,*avedr,*aredf,*avedf,*arstr,
                     *avstr,*arstf,*avstf;


    vector< vector<double> > forwardshockradiusprofile;
    vector< vector<double> > forwardshockvelocityprofile;
    vector< vector<double> > reverseshockradiusprofile;
    vector< vector<double> > reverseshockvelocityprofile;
    vector< vector<double> > cdradiusprofile;
    vector< vector<double> > cdvelocityprofile;


    vector< vector<double> > thinshellradius;
    vector< vector<double> > thinshellvelocity;
    double thinshellsolutiontmin;
    double thinshellsolutiontmax;

    double truelovemckeetminfed;
    double truelovemckeetminfst;
    double truelovemckeetmaxfed;
    double truelovemckeetmaxfst;

    double truelovemckeetminred;
    double truelovemckeetminrst;
    double truelovemckeetmaxred;
    double truelovemckeetmaxrst;

    vector<double> truelovemckeeparams;
    vector< vector<double> > r_edforward;
    vector< vector<double> > v_edforward;
    vector< vector<double> > r_edreverse;
    vector< vector<double> > v_edreverse;
    vector< vector<double> > r_stforward;
    vector< vector<double> > v_stforward;
    vector< vector<double> > r_streverse;
    vector< vector<double> > v_streverse;
    vector<double> TrueloveMcKeeRadiusAndSpeed(double t, bool REVERSE);
    int TrueloveMcKeeEjectaDominatedForwardShock(double x,
                                                 double &Y,double &V);
    int TrueloveMcKeeEjectaDominatedReverseShock(double x,
                                                 double &Y, double &V);
    int TrueloveMcKeeSedovTaylorPhaseForwardShock(double t,
                                                  double &R, double &V);
    int TrueloveMcKeeSedovTaylorPhaseReverseShock(double t,
                                                  double &R, double &V);
    void CalculateTrueLoveMcKeeParams(vector<double> pars);
    vector< double > BFieldDiskComponent(vector<double> xyz);
    vector<double> BFieldRegularHaloComponent(vector<double> xyz);
    vector<double> BFieldXComponent(vector<double> xyz);
    vector<double> BFieldRandomComponent(vector<double> xyz);
    double BFieldRandomRmsStrength(vector<double> xyz);
    vector<double> BFieldStriatedComponent(vector<double> RegularDisk, 
                                           vector<double> RegularHalo,
                                           vector<double> RegularX);
    vector< double > RotateBFieldComponent(vector<double> component, 
                                           vector<double> xyz);
    int BFieldGetSpiralArm(vector<double> xyz);
  public:
    Astro();
    ~Astro();
    void SetGalacticReferencePoint(vector<double> xyz_ref);
    void DisableArm(int arm);///< Switch off an individual arm in Galaxy spiral model
    void EnableArm(int arm);///< Switch on an individual arm in Galaxy spiral model
    double RandomSalpeterInitialMass();///< Dice an initial star mass following the Salpeter law
    double MainSequenceBubbleRadius(double mDotMS, double vMSWind, double timeMS, double n, double initialMass);///< calculate the radius created by the main sequence wind of a star
    double CalculateMSBubbleDensity(double mDotMS, double vMSWind, double timeMS, double n, double bubbleRadius);///< Density inside the main-sequence wind blown bubble.
    double CalculateTimeOnMainSequence(double initialMass);///< Time on the main sequence in years as a function of initial star mass
    double CalculateMSMassLuminosity(double initialMass);///< Main sequence wind mass luminosity as a function of initial star mass
    double CalculateMSWindSpeed(double initialMass);///< Main sequence wind speed as a function of initial star mass
    double CalculateRGWRadius(double pBubble, double mDotRGW, double vRGWind);///< Radius of the red giant wind zone
/*    void CalculateBField(double x, double y, double z, double &B_tot, double &B_coh, double &B_ord, double &B_iso, vector<double> &regFieldDirection, vector<double> &isoFieldDirection);///< 2D-Model of the large-scale galactic magnetic field structure*/
    double HTotalDensity(double x, double y, double z);
/*    double TotalBField(double x, double y, double z);*/
    double HIDensity(double x, double y, double z, bool MODULATE = true);///< TODO:COMMENT
    double H2Density(double x, double y, double z, bool MODULATE = true);///< TODO:COMMENT
    double HIIDensity(double x, double y, double z, bool MODULATE = true);///< TODO:COMMENT
    double H_I_2_CMZ(double x, double y, double z, double H_c, double norm);
    double H_I_2_InnerDisk(double x, double y, double z, double H_d, double norm);
    double H2CMZ(double x, double y, double z);
    double HICMZ(double x, double y, double z);
    double H2InnerDisk(double x, double y, double z);
    double HIIThickDisk(double x, double y, double z);
    double HIIThinDisk(double x, double y, double z);
    double HIInnerDisk(double x, double y, double z);
    double WarmHIIDensity(double x, double y, double z);
    double HotHIIDensity(double x, double y, double z);
    double VeryHotHIIDensity(double x, double y, double z);
    double GravitationalPotential(double r, double z);
    double CalculateHINorm(double R);///< TODO:COMMENT
    double CalculateH2Norm(double R);///< TODO:COMMENT
    double GetHIColumnDensity(double R);///< TODO:COMMENT
    double GetH2ColumnDensity(double R);///< TODO:COMMENT
    double GetHIFWHM(double R);///< TODO:COMMENT
    double GetH2FWHM(double R);///< TODO:COMMENT
    double ModulateGasDensityWithSpirals(double n, double x, double y, double z);
    
    // Functions for the milky way hydrogen density according to Ferriere 1998 and 2007
    double GetTotalHDensityFerriere(vector< double > xyz);
    double GetH2DensityFerriere(vector< double > xyz); 
    double GetCMDensityFerriere(vector< double > xyz);
    double GetWNMDensityFerriere(vector< double > xyz);
    double GetWIMDensityFerriere(vector< double > xyz);
    double GetHIMDensityFerriere(vector< double > xyz);
    double GetH2DensityCMZFerriere(vector< double > xyz);
    double GetH2DensityDiskFerriere(vector< double > xyz);
    double GetHIDensityCMZFerriere(vector< double > xyz);
    double GetHIDensityDiskFerriere(vector< double > xyz);
    double GetHIIDensityWIMFerriere(vector< double > xyz);
    double GetHIIDensityHIMFerriere(vector< double > xyz);
    double GetHIIDensityVHIMFerriere(vector < double > xyz);
    double phi_RZ(double r, double z);
    
    
    
    //double CalculateGasColumnDensity(vector<double> xyzReference,vector<double> GLGB,string gascomponent,double modulate,double range,double steps);///< calculate gas column densities
/*    double nRadial(double *x, double *pars);///< TODO:COMMENT*/
    vector<double> GetCartesian(vector<double> lbr, vector<double> xyzobs);///< Galactic coordinates (GL,GB,R) ->Cartesian
    vector<double> GetCartesian(double r, double l, double b, vector<double> xyzobs);
    vector< vector<double> > GetCartesianPositions(vector< vector<double> > lbr, vector<double> xyzobs);
    vector<double> GetGalactic(vector<double> xyz, vector<double> xyzref);///< Cartesian coordinates -> Galactic Coordinates(GL,GB,R)
    void GetGalactic(double x, double y, double z, double xref, double yref, double zref, double &l, double &b);
    vector< vector<double> > GetGalacticPositions(vector< vector<double> > xyz, vector<double> xyzref);
    void RotateCoordinates(double &x, double &y, double &z, double phi, double theta, double psi);///< rotate coordinate around (0,0,0) by (phi,theta,psi)
    void RandomTangentialShift(double r, double width, double &x, double &y);///< gaussian shift along concentric circle (retains r-distribution)
    void RandomGaussianShift(double r, double width, double &x, double &y);
    double ApparentSize(double size, vector<double> lbr);
    vector<double> GetApparentSizes(vector<double> sizes, vector< vector<double> > lbr);
    unsigned int GetRandomArm(double r);///< return a random spiral arm (identifier, i.e. an unsigned integer)
    void Bar(double r, double &x, double &y);///< Get xy Coordinates on a central bar inclined with an angle 'angle' clockwise
    void Disk(double r, double &x, double &y);///< Get random xy Coordinates  on a disk
    double LinearSurfDens(double *x, double *par);///< radial distribution that gives a flat appearance of the galaxy
    void TaylorCordesSpiral(double r, int arm, double &x, double &y);///<TODO: Comment
    double TaylorCordesSpiralAngular(double theta, int arm);///< TODO: Comment
    void ValleeSpiral(double r, int arm, double &x, double &y);///<TODO: Comment
    double ValleeSpiralAngular(double theta, int arm);///<TODO: Comment
    double GetRandomZValue(double scaleHeight);///<TODO: Comment
    vector<double> GetRandomGalactocentricRadii(int n);///<TODO: Comment
    void GalacticPositionXY(double r, int arm, double &x, double &y);///<Return x,y coordinates at given galactocentric radius r and on arm number i
    void PositionOnSpiralArmAngular(double theta, int arm, double &x, double &y); ///< get position on arm in x-y plane given a certain polar angle in the galaxy.
    void GetDistanceToNearestSpiralArm(double x, double y, double &DistanceToClosestArm, int &ClosestArm);
    double GetDistanceToGivenSpiralArm(double x, double y, int arm);
    vector< vector<double> > DiceGalacticPositions(int i);
    double GetRMax() {return rMax;}
    void SetSpiralArmModel(string centralstructuremodel);
    void SetCentralStructureModel(string centralstructuremodel);
    void SetSurfaceDensityModel(string surfacedensitymodel);
    void SetMainSequenceBubbleModel(string msbubblemodel);
    void SetScaleHeight(double SCALEHEIGHT) {scaleHeight=SCALEHEIGHT;}
    void SetArmWidth(double ARMWIDTH) {armWidth=ARMWIDTH;}
    vector< vector<double> > CreateDensityProfile(vector<double> pars,
                                                  double rmin = 1.e-3,
                                                  double rmax = 1.e3);
    double EvalTaylorCordesArmTheta(double r, int arm);
    double EvalTaylorCordesArmGalactocentricRadius(double theta, int arm);
    double CaseBhattacharyaProfile(double r);
    double IusifovKucukProfile(double r);
    void CalculateThinShellApproximation(vector< vector<double> > dProfile,
                                         double E,
                                         double AdiabIndex = 1.333,
                                         double tmin = 1.e-2,
                                         double tmax = 1.e4,
                                         int steps = 200);
    void CalculateTrueloveMcKeeSolution(vector<double> pars,
                                        double tmin = 1.e-2,
                                        double tmax = 1.e4,
                                        int steps = 200);
    void CalculateTangChevalierSolution(vector<double> pars, double tmin=1.e-2,
                                           double tmax=1.e4, int steps=200);
    vector<double> ThinShellRadiusAndSpeed(double t);
    vector<double> TrueloveMcKeeRadiusAndSpeedForward(double t) {
                                  return TrueloveMcKeeRadiusAndSpeed(t,false);}
    vector<double> TrueloveMcKeeRadiusAndSpeedReverse(double t) {
                                  return TrueloveMcKeeRadiusAndSpeed(t,true);}
    vector< vector<double> > GetForwardShockRadiusEvolution() {
                                  return forwardshockradiusprofile;}
    vector< vector<double> > GetForwardShockVelocityEvolution() {
                                  return forwardshockvelocityprofile;}
    vector< vector<double> > GetReverseShockRadiusEvolution() {
                                  return reverseshockradiusprofile;}
    vector< vector<double> > GetReverseShockVelocityEvolution() {
                                  return reverseshockvelocityprofile;}
    vector< vector<double> > GetContactDiscontinuityRadiusEvolution() {
                                  return cdradiusprofile;}
    vector< vector<double> > GetContactDiscontinuityVelocityEvolution() {
                                  return cdvelocityprofile;}
    void CalculateForwardShockInRGWind(vector<double> pars,
                                       double tmin = 1.e-2,
                                       double tmax = 1.e4,
                                       int steps = 200);
    vector<double> ForwardShockInRGWind(double t, vector<double> pars,
                                       vector< vector<double> > profile);

    void SetInterpolationMethod(string intermeth)
       {fUtils->SetInterpolationMethod(intermeth);}
    void ToggleQuietMode() { QUIETMODE = QUIETMODE == true ? false : true; }
    bool GetQuietMode() {return QUIETMODE;}
    vector<double> CalculateBField(vector<double> xyz, int component = -1);
    vector< vector<double> > CalculateBFieldPositions(vector< vector<double> > pos, int component=-1);
    vector< vector<double> > CalculateBFieldComponents(vector<double> xyz, int component = -1);
    vector<double> CalculateBFieldStrengthComposition(vector<double> xyz);
    double CalculateBFieldStrength(vector<double> xyz);
    vector< vector<double> > LineOfSight(double gl, double gb, vector<double> obs, double r_max=50, int steps=200);
    vector< vector<double> > LineOfSight(double gl, double gb, vector<double> obs, vector<double> rvals);
    
};
#endif
