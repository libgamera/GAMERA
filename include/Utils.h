#ifndef _UTILS_
#define _UTILS_

#include <math.h>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <ctime>
#include <algorithm>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_spline.h>
#include "2D_interp/interp2d.h"
#include "2D_interp/interp2d_spline.h"

/* TeV->erg */
#define TeV_to_erg 1.602
/* GeV->erg */
#define GeV_to_erg 1.602e-3
/* eV->erg */
#define eV_to_erg 1.602e-12
/* Thomson cross section */
#define sigma_T 6.6524e-25
/* classical electron radius (cm) */
#define e_radius 2.8179e-13
/* parsec to cm */
#define pc_to_cm 3.0857e18
/* kiloparsec to cm */
#define kpc_to_cm 3.0857e21
/* AU to cm*/
#define AU_to_cm 1.496e13
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
/* electron mass in g */
#define m_e_g = 9.1093837015e-28
/* boltzmann constant (erg/K) */
#define kb 1.380658e-16
/* proton mass in erg */
#define m_p 1.50310854e-3
/* pi0 mass in erg */
#define m_pi 2.1622194e-4
/* parsec to cm */
#define pc_to_cm 3.0857e18
/* well... pi! */
#define pi 3.14159265359
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
/* natural logarithm of 10 */
#define ln10 2.302585

using namespace std;

/**
 * @file Utils.C
 * @author Joachim Hahn
 * @short Class that provides utility functions
 * Provides different things, mainly useful for population studies.
 * This includes spiral galaxy models as well as gas distributions,
 * B-Field distributions, coordinate transformations etc.
 * Also includes more technical functions that are helpful in writing
 * GAMERA programs.
 */

bool sortcriterionfirstcolumn (vector<double> i,vector<double> j);
bool sortcriterionsecondcolumn (vector<double> i,vector<double> j);

class Utils {
 private:
  bool QUIETMODE;
  bool GAMERADESTROYEDTHECONSOLE;
  gsl_interp_type *INTERMETH;
  gsl_rng * r;
  gsl_spline *INTERNAL1DSPLINE;
  interp2d_spline *INTERNAL2DSPLINE;

  gsl_interp_accel *acc;
  gsl_interp_accel *xacc;
  gsl_interp_accel *yacc;

 public:
  Utils(bool DRAWLOGO = false);
  ~Utils();
  int ReadParameterFile(string inputname, vector<string> parameter_names,
                        vector<string> files_names,
                        vector<double> &parameter_values,
                        vector<string> &files);
  void DrawGamera();
  void DrawGappa();
  void WriteOut(vector<vector<double> > sp, string outname);
  void ReadIn(string inname, vector< vector<double> > &sp);
  double Random();
  vector<double>  UniformRandom(double x_min, double x_max,int n);
  vector<double>  LinearRandom(double slope, double x_min, double x_max, int n);
  vector<double>  PowerLawRandom(double index, double x_min, double x_max,
                                 int n);
  vector<double>  GaussianRandom(double width, double offset, int n);
  vector<double>  SignRandom(int n);
  vector<double>  ExponentialRandom(double ind_norm, double x_min,
                                    double x_max, int n);
  vector<double>  CustomFunctionRandom(vector< vector<double> > f, int n,
                                       double xmin = 0., double xmax = 0.);
  vector< vector<double> > CustomFunctionRandom2D(vector< vector<double> > f, 
                                           int n,
                                           double xmin=0., double xmax=0.,
                                           double ymin=0., double ymax=0.);
  double Integrate(vector< vector<double> > f, double xmin=0., double xmax=0.);
  vector< vector< double> > IntegratedProfile(vector< vector<double> > f);
  gsl_spline *GSLsplineFromTwoDVector(vector< vector<double> > v);
  vector< vector<double> > SortTwoDVector(vector< vector<double> > v,
                                                 int column);
  void TwoDVectorPushBack(double x, double y, vector< vector<double> > &v);
  void TwoDVectorPushBack(double x, double y, double z, vector< vector<double> > &v);

  double EvalSpline(double x, gsl_spline *s, gsl_interp_accel *a,
                    const char* t, int l);
  void SetInterpolationMethod(string intermeth);
  void Clear2DVector(vector< vector<double> > &v);
  vector< vector<double> > VectorAxisLogarithm(vector< vector<double> > v,
                                               unsigned int column);

  vector< vector<double> > VectorAxisPow10(vector< vector<double> > v,
                                                    int column);
  interp2d_spline *TwoDsplineFromTwoDVector(vector< vector<double> > v, 
                                            double &xmin, double &xmax,
                                            double &ymin, double &ymax);
  void SetInternal1DSpline(vector< vector<double> > v) {
    gsl_interp_accel_reset(acc);
    if(INTERNAL1DSPLINE!=NULL) gsl_spline_free(INTERNAL1DSPLINE);
    INTERNAL1DSPLINE = GSLsplineFromTwoDVector(v);} 
  void SetInternal2DSpline(vector< vector<double> > v) {
//    std::cout<<" -- - --- --- " <<std::endl;
//    for(unsigned int i=0;i<v.size();i++) {
//      for(unsigned int j=0;j<v[i].size();j++) {
//          std::cout<<" "<<v[i][j];
//      } 
//      std::cout<<std::endl;
//    }
//    std::cout<<" -- - --- --- " <<std::endl;
    gsl_interp_accel_reset(xacc);
    gsl_interp_accel_reset(yacc);
//    std::cout<<INTERNAL2DSPLINE<<std::endl;
//    printf("%p\n", INTERNAL2DSPLINE);
    if(INTERNAL2DSPLINE!=NULL) interp2d_spline_free(INTERNAL2DSPLINE);
//    printf("%p\n", INTERNAL2DSPLINE);
    double xmin,xmax,ymin,ymax;
    INTERNAL2DSPLINE = TwoDsplineFromTwoDVector(v,xmin,xmax,ymin,ymax);}
  double EvalInternal1DSpline(double x) {
    return gsl_spline_eval(INTERNAL1DSPLINE, x, acc);}
  double EvalInternal2DSpline(double x, double y) {
    return interp2d_spline_eval(INTERNAL2DSPLINE,x,y,xacc,yacc);}
  double EnergyContent(vector< vector<double> > f, double emin=0., double emax=0.);
  double EnergyContentFast(vector< vector<double> > f);
  vector< vector<double> > MeshgridToTwoDVector(vector<double> x, vector<double> y, 
                                                     vector< vector<double> > mesh);
/*  void TwoDVectorToMeshgrid(vector< vector<double> > vec, vector<double> &x,*/
/*                                vector<double> &y, vector< vector<double> > &mesh);*/
  vector<double> GetVectorMinMax(vector< vector<double> > v, unsigned int axis);
  vector<unsigned int> GetVectorMinMaxIndices(vector< vector<double> > v, unsigned int axis);
  void ToggleQuietMode() { QUIETMODE = QUIETMODE == true ? false : true; }
  bool GetQuietMode() {return QUIETMODE;}
  vector< vector<double> > RemoveZeroEntries(vector< vector<double> > v);
  vector< vector<double> > TwoDVectorFabs(vector< vector<double> > v);
  vector<double> RotateVector(vector<double> v, vector<double> p1, vector<double> p2);
  double Gaussian1D(double x, double sigma, double mu=0., double norm=0.);
  double LogisticsFunction(double z, double h, double w);
};
#endif
