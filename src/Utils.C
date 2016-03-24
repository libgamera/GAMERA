#include "Utils.h"

Utils::Utils(bool DRAWLOGO) {
  if (DRAWLOGO == true) DrawGamera();
  QUIETMODE = false;
  GAMERADESTROYEDTHECONSOLE = false;
  INTERMETH = (gsl_interp_type*)gsl_interp_cspline;
  r = gsl_rng_alloc (gsl_rng_mt19937);
  gsl_rng_set(r,time(0));

  INTERNAL1DSPLINE = NULL;
  INTERNAL2DSPLINE = NULL;

  acc = gsl_interp_accel_alloc();
  xacc = gsl_interp_accel_alloc();
  yacc = gsl_interp_accel_alloc();
}

Utils::~Utils() {}

/**
 * Generic parameter file reading function. Fill a vector with parameter titles
 * "parameter_names" (e.g. filled with "Distance", "E0") and another one with
 * string titles that are to be involved "strings_names" (e.g. filled with
 * "outputfile", "datafile"). Keep the same nomenklature in the config file
 * (i.e. inputname), here this would be :
 *
 * ####################################
 * # Distance 1.e3                    #
 * # E0 1.e51                         #
 * # outputfile "results/results.dat" #
 * # datafile "data/points.dat"       #
 * ####################################

 * This function will scan through the parameter file and fill the vectors
 "parameter_values"
 * and "strings_values" with the korresponding values provided in the parameter
 file. It keeps the
 * order of the input vectors.
 */
int Utils::ReadParameterFile(string inputname, vector<string> parameter_names,
                             vector<string> strings_names,
                             vector<double> &parameter_values,
                             vector<string> &strings_values) {

  if (!inputname.size()) {
    cout << ">> Please provide parameter file! -> exiting ... " << endl;
    return 1;
  }
  ifstream inlist(inputname.c_str());
  if (inlist.is_open())
    cout << ">> Reading parameters from " << inputname << endl;
  else
    cout << ">> Problems reading parameters from " << inputname << endl;
  parameter_values.resize(parameter_names.size());
  strings_values.resize(strings_names.size());

  string parname;
  string parval;
  while (inlist) {
    if (inlist.eof()) break;
    inlist >> parname >> parval;
    for (unsigned int i = 0; i < parameter_values.size(); i++) {
      if (!parname.compare(parameter_names[i]))
        parameter_values[i] = atof(parval.c_str());
    }
    for (unsigned int i = 0; i < strings_names.size(); i++) {
      if (!parname.compare(strings_names[i]))
        strings_values[i] = parval.c_str();
    }
  }
  inlist.close();
  for (unsigned int i = 0; i < parameter_values.size(); i++) {
    if (parameter_values[i] < 0.) {
      cout << "Utils::ReadParameterFile: Something wrong with the \"" << parameter_names[i]
           << "\" parameter! -> exiting... " << endl;
      return 1;
    }
  }
  for (unsigned int i = 0; i < strings_values.size(); i++) {
    if (!strings_values[i].compare("")) {
      cout << "Utils::ReadParameterFile:  You have to specify the name of \""
           << strings_names[i] << "\"! -> exiting... " << endl;
      return 1;
    }
  }
  cout << ">> PARAMETERS: " << endl;
  for (unsigned int i = 0; i < parameter_values.size(); i++) {
    cout << "    " << parameter_names[i] << ": " << parameter_values[i] << endl;
  }
  for (unsigned int i = 0; i < strings_names.size(); i++) {
    cout << "    " << strings_names[i] << ": " << strings_values[i] << endl;
  }
  cout << endl;
  return 0;
}

void Utils::WriteOut(vector<vector<double> > sp, string outname) {

  ofstream output_file(outname.c_str());
  if (output_file.is_open())
    cout << ">> Writing to " << outname << endl;
  else
    cout << ">> Problem with " << outname << endl;

  for (unsigned int i = 0; i < sp.size(); i++) {
    for (unsigned int j = 0; j < sp[i].size(); j++)
      output_file << sp[i][j] << " ";
    output_file << endl;
  }
  output_file.close();
}

void Utils::ReadIn(string inname, vector< vector<double> > &sp) {
  sp.clear();
  ifstream input_file(inname.c_str());
  if (input_file.is_open())
    cout << ">> Reading from " << inname << endl;
  else
    cout << ">> Problem with " << inname << endl;
  while(1){
    if(input_file.eof()) break;
    double E,N;
    E=N=0.;
    input_file>>E>>N;
    if(E&&N) {
      sp.push_back(vector<double>());
      sp[sp.size()-1].push_back(E);
      sp[sp.size()-1].push_back(N);
    }
  }
  input_file.close();
  return;
}

void Utils::DrawGamera() {
  if (GAMERADESTROYEDTHECONSOLE) return;
  cout << endl;
  cout << endl;
  cout << "                                                  ▄▄▄            " << endl;
  cout << "                                              ▄▓▌▒▀▀▀▀▓▓▄▄▄      " << endl;
  cout << "                                            ▓▀░▄▒▒   ▀    ░▀▄    " << endl;
  cout << "                                ▒▒░░  ░░░▓▒▒▒▀▀▒▒▒░▒░░░ ░▒▒░░░▒  " << endl;
  cout << "                            ▓░░░░      ░▒▒▒░▒▒▒░░▒░░ ▐█░░▌       " << endl;
  cout << "                         ▄▓░ ░░          ░░▒░░░░░░  ▄█▀ ▀  ▄▀    " << endl;
  cout << "                       ▓▓░       ░        ░░░░▒▒▒▀▀▒░░▀▐▓▄▀      " << endl;
  cout << "                     ▓█▀░    ▒░░░░░▄   ░   ░░░▄▓       ▀▒▓▓      " << endl;
  cout << "                    ░▀░░    ▀ ▒▒░░▒░ ▒░░░░▒ ▒▀░░ ░░              " << endl;
  cout << "                    ▒░    ░▐▐▒▒▓▒▐▌░░░▒░▒ ░▒▓░ ░ ░▀░░▓           " << endl;
  cout << "                 ▓▌░░    ▒░░░▒░░▒▓▓░░ ░░░░░░░▒░░ ▀ ░░     ░▌░ ░  " << endl;
  cout << "                ▄▀░░   ▒▒░░░░░░▒▒░░▒▒▄▄▓▓▓▓▓▒▒▄░░ ░░░░░▒░▒░░░░ ░ " << endl;
  cout << "               ▓░░░  ░▓░░░░▒▒░▒▓▓░▓▓▀▓▓▒▓▓▓▓▓▒      ▀            " << endl;
  cout << "              ▄▒░░   ▒░░░▒▒▒░▒▒▒▒▀▀▒▒▀▒▀▒░▀▒▀                    " << endl;
  cout << "            ▒░▒░  ░░▒░░▒▒▒░░▀▒▒░▀░▒▀░▒░▒░▒▓                      " << endl;
  cout << "            ░▒░░ ░░▒▒░▒▒▒░▓▒░░▀▀░░▀░░░░▒▒▀                       " << endl;
  cout << "          ▒░░░░  ░▒░░░▒▒░▒▒▒░▒▒░░▒░░░░▒▓                         " << endl;
  cout << "           ░▒░   ░░░░░░░░░░░░░░▒▒░░ ░░░                          " << endl;
  cout << "          ░░░░   ░░    ░░░░░▒░░▒▒   ░░░░                         " << endl;
  cout << "          ▒░░░ ░░ ▒░▒▒▒▒░░░░░▒░░      ░░▒                        " << endl;
  cout << "           ▒░     ▒▓▒░░░░░▒▒▒▒       ░░▒▒▓                       " << endl;
  cout << " ▓▄▄         ▒░   ▀▒░░░░░░▒░▒    ░  ░  ░░▓                       " << endl;
  cout << "  ░▓▒▄         ▄▒▒▒▒▒░░░▒▒▒▒         ░░░░▒                       " << endl;
  cout << "   ▒░▄▒▄     ▓▓▒▒░ ░░░░░░▒░          ░░░░░                       " << endl;
  cout << "    ░▒▀ ▒▓▒▒▒▓▒░░░░░░░░░▒         ░   ░░░▒                       " << endl;
  cout << "      ░▒▒░░░░▒░░░░░░▒░▀              ░░░░                        " << endl;
  cout << "        ▒░░░▓▒▀░░░▒▀                     ░                       " << endl;
  cout << "           ▓▒░ ░░░                      ░░░░▒▒▒▒                 " << endl;
  cout << "          ▓▓▓░░▓▒▒▒▒▄▒              ▒░      ▒  ░▓▀               " << endl;
  cout << "           ▀░█▄▐█▄ ▐█▄██                                         " << endl;
  cout << "              ▀  ▀                                               " << endl;
  cout << endl;
  cout << "       ▄████  ▄▄▄       ███▄ ▄███▓▓█████  ██▀███   ▄▄▄      " << endl;
  cout << "      ██▒ ▀█▒▒████▄    ▓██▒▀█▀ ██▒▓█   ▀ ▓██ ▒ ██▒▒████▄    " << endl;
  cout << "     ▒██░▄▄▄░▒██  ▀█▄  ▓██    ▓██░▒███   ▓██ ░▄█ ▒▒██  ▀█▄  " << endl;
  cout << "     ░▓█  ██▓░██▄▄▄▄██ ▒██    ▒██ ▒▓█  ▄ ▒██▀▀█▄  ░██▄▄▄▄██ " << endl;
  cout << "     ░▒▓███▀▒ ▓█   ▓██▒▒██▒   ░██▒░▒████▒░██▓ ▒██▒ ▓█   ▓██▒" << endl;
  cout << "      ░▒   ▒  ▒▒   ▓▒█░░ ▒░   ░  ░░░ ▒░ ░░ ▒▓ ░▒▓░ ▒▒   ▓▒█░" << endl;
  cout << "       ░   ░   ▒   ▒▒ ░░  ░      ░ ░ ░  ░  ░▒ ░ ▒░  ▒   ▒▒ ░" << endl;
  cout << "     ░ ░   ░   ░   ▒   ░      ░      ░     ░░   ░   ░   ▒   " << endl;
  cout << "           ░       ░  ░       ░      ░  ░   ░           ░  ░" << endl;
  cout << endl;
  cout << endl;

  GAMERADESTROYEDTHECONSOLE = true;
}


void Utils::DrawGappa() {
  if (GAMERADESTROYEDTHECONSOLE) return;

  cout << endl;
  cout << endl;
  cout << "                                            ▄▄                                  " << endl;
  cout << "                                            ▐▓█▄                                " << endl;
  cout << "                                             ▓▓▒▓                               " << endl;
  cout << "                                           ▀▓▓▓▄▄▀░                             " << endl;
  cout << "                                            ▓▓▒▒░░ ▄░                   ▄▓▓     " << endl;
  cout << "                                            ▓▒▒▀▒▓█▌▒▒             ▄▄▓▓▓▓▓▀     " << endl;
  cout << "                                            ▀▒▒▒░ ░░     ░░░░░░░░▒▓▒▓▒▒▒▒▌      " << endl;
  cout << "  ▄▄             ▄ ▄▄▄▄▄  ▒ ▒▒▀▀▒▀▀  ░░     ░▓▒░░   ░    ░░ ░░░▒▒▓▒▒░▒▀▒▀       " << endl;
  cout << "  ▀▓▓▒░▄▄   ░  ░░░░░░                ░   ░  ▓▓▒░    ░   ░ ░  ░▄▒▒░▒▒▀▓▓         " << endl;
  cout << "   ▀▓▒▒▒▒▒▒▒░░░░                ░ ░  ▒   ░░▓▓▒▒▒        ░ ▒░   ░▀▀░▒▓▓▀         " << endl;
  cout << "    ▐▌▒▒░░░▒░    ░            ░░░       ░▒▌▒░▀░ ▒▄░       ░░ ▀▒░░▒▒▓▓           " << endl;
  cout << "      ▓░▓▀░░░         ▄  ▄▀  ░░░░░░░░ ░ ▓▓▀░  ▐░  ░  ░▒  ░    ▓▒▀▒              " << endl;
  cout << "                  ▄  ▄▒▓ ▄▀ ▒           ▓▓▒▒ ▄▀   ░  ░░      ▒▓                 " << endl;
  cout << "              ▄▄ ▄▄▓▓▀ ▒  ░   ░        ░▐▒▓▓▓▒                                  " << endl;
  cout << "               ▓▓▓▀  ▒                 ░▐▓▐▌░          ░▐                       " << endl;
  cout << "                    ░               ░▒▒▐▓▒░░  ▒      ░ ░                        " << endl;
  cout << "                      ░            ░░░ ▓▓▒    ░░     ░ ░                        " << endl;
  cout << "                                  ▄▒░░░▓▓▒▒░            ▒                       " << endl;
  cout << "                                   ▒▒░▐▓▒▒▒░▒           ▒                       " << endl;
  cout << "                                     ░ ▓▒▓▓░ ▒ ░                                " << endl;
  cout << "                                       ▐▒▒▒▒  ░▒                                " << endl;
  cout << "                                         ░▒▒░  ▒▒▓▒                             " << endl;
  cout << "                                        ▄▓▒▒   ░░ ░                             " << endl;
  cout << "                                      ▐▓▀░▒   ▐▒░                               " << endl;
  cout << "                                     ▄▀▀░     ▀▀                                " << endl;
  cout << "                                    ▓▌▒░                                        " << endl;
  cout << "                                    ▓▓  ▄                                       " << endl;
  cout << "                                    ▓▄   ▀▐▄                                    " << endl;
  cout << endl;
  cout << "                     ▄████  ▄▄▄       ██▓███   ██▓███   ▄▄▄      " << endl;
  cout << "                    ██▒ ▀█▒▒████▄    ▓██░  ██▒▓██░  ██▒▒████▄    " << endl;
  cout << "                   ▒██░▄▄▄░▒██  ▀█▄  ▓██░ ██▓▒▓██░ ██▓▒▒██  ▀█▄  " << endl;
  cout << "                   ░▓█  ██▓░██▄▄▄▄██ ▒██▄█▓▒ ▒▒██▄█▓▒ ▒░██▄▄▄▄██ " << endl;
  cout << "                   ░▒▓███▀▒ ▓█   ▓██▒▒██▒ ░  ░▒██▒ ░  ░ ▓█   ▓██▒" << endl;
  cout << "                    ░▒   ▒  ▒▒   ▓▒█░▒▓▒░ ░  ░▒▓▒░ ░  ░ ▒▒   ▓▒█░" << endl;
  cout << "                     ░   ░   ▒   ▒▒ ░░▒ ░     ░▒ ░       ▒   ▒▒ ░" << endl;
  cout << "                   ░ ░   ░   ░   ▒   ░░       ░░         ░   ▒   " << endl;
  cout << "                         ░       ░  ░                        ░  ░" << endl;
  cout << "                                                                 " << endl;
  cout << endl;
  cout << endl;
  return;
}


/**
 * Return uniform random number in [0,1). GSL Wrapper function.
 */
double Utils::Random() {
  return gsl_rng_uniform(r);
}

/**
* get linearly distributed variates
*/
vector<double> Utils::LinearRandom(double slope, double x_min,
                                   double x_max, int n) {
  vector<double> v;
  for(int i=0;i<n;i++) {
    if(!slope) v.push_back((x_max-x_min)*gsl_rng_uniform(r)+x_min);
    else {
      double rand_min = 0.5*slope*pow(x_min,2.);
      double rand_max = 0.5*slope*pow(x_max,2.);
      double rand = rand_min + (rand_max - rand_min)*gsl_rng_uniform(r);
      v.push_back(sqrt(2.*rand/slope));
    }
  }
  return v;
}

/**
* get power-law distributed variates
*/
vector<double> Utils::PowerLawRandom(double index, double x_min,
                                     double x_max, int n) {
  vector<double> v;
  for(int i=0;i<n;i++) {
    double rand = -100.;
    double k = 0.03;
    while(rand<x_min || rand>x_max) {
      rand=pow(((index-1.)/k)*gsl_rng_uniform(r),-1./(index-1.));
    }
    v.push_back(rand);
  }
  return v;
}

/**
* get gaussian distributed variate using the Box-Muller method
*/
vector<double> Utils::GaussianRandom(double width, double offset, int n) {
  vector<double> v;
  for(int i=0;i<n;i++) {
    double random_gaussian = sqrt(-2.*log(gsl_rng_uniform(r)))
                        * cos(2*pi*gsl_rng_uniform(r));
    v.push_back(offset + width*random_gaussian);
  }
  return v;
}

/**
* return a random sign
*/
vector<double> Utils::SignRandom(int n) {
  vector<double> v;
  for(int i=0;i<n;i++) {
    if(gsl_rng_uniform(r)< 0.5) v.push_back(1.);
    else v.push_back(-1.);
  }
  return v;
}

/**
* return a exponentially distributed variate
*/
vector<double> Utils::ExponentialRandom(double ind_norm, double x_min,
                                        double x_max, int n) {
  vector<double> v;
  ind_norm = -1./ind_norm;
  double rand_min = (1./ind_norm)*(exp(ind_norm*x_min)-1.);
  double rand_max = (1./ind_norm)*(exp(ind_norm*x_max)-1.);
  for(int i=0;i<n;i++) {
    double rand = rand_min + (rand_max-rand_min)*gsl_rng_uniform(r);
    v.push_back(pow(ind_norm,-1.)*log((ind_norm)*rand+1.));
  }
  return v;
}

/**
 * Sample Variate that follows the input 2D-Vector
 */
vector<double> Utils::CustomFunctionRandom(vector< vector<double> > f,
                                           double xmin, double xmax, int n) {

                                             cout<<"1.0"<<std::endl;
  vector<double> v;
  if(!f.size()) {
    cout << "Utils::CustomFunctionRandom: function vector empty! "
            "Exiting!" << endl;
    return v;
  }
  if(!xmin && !xmax) {
    xmin = f[0][0];
    xmax = f[f.size()-1][0];
  }
  if( xmin < f[0][0] || xmax > f[f.size()-1][0] ) {
    cout << "Utils::CustomFunctionRandom: requested sampling range outside "
            "of boundaries req.:(" << xmin << "," << xmax << ") vs. "
            "avail.:(" << f[0][0] << "," << f[f.size()-1][0] << "). "
            "Returning emptyvector." << endl;
    return v;
  }
  cout<<"1.1"<<std::endl;
  int size = (int)f.size();
  double x[size];
  double y[size];
  double ymax = -1.e-100;
  for (unsigned int i=0;i<f.size();i++) {
    x[i] = f[i][0];
    y[i] = f[i][1];
    if(x[i] > xmin && x[i] < xmax && y[i] > ymax) ymax = y[i];
  }
  gsl_spline *lookup = gsl_spline_alloc(gsl_interp_linear, size);
  gsl_spline_init(lookup, x, y, size);
  gsl_interp_accel *a = gsl_interp_accel_alloc();

  cout<<"1.2"<<std::endl;
  for(int i=0;i<n;i++) {
    double x;
    while(1) {
      x = (xmax-xmin)*gsl_rng_uniform(r)+xmin;
      double u = gsl_rng_uniform(r);
      double val = 0.;
      if (gsl_spline_eval_e(lookup, x, a, &val)) {
        cout << "Utils::CustomFunctionRandom: Function interpolation "
        "failed. Exiting!" << endl;
        return v;
      }
      if (std::isnan(val) || std::isinf(val)) {
        cout << "Utils::CustomFunctionRandom: Function interpolation returned "
             << val << "! Exiting! " << endl;
        return v;
      }
      if(u < val/ymax) break;
    }
    v.push_back(x);
  }
  return v;
}

double Utils::EnergyContent(vector< vector<double> > f, double emin, double emax) {
  if(!f.size()) {
    cout << "Utils::EnergyContent: spectrum vector empty! "
    "Exiting & returning 0." << endl;
    return 0.;
  }
  for(unsigned int i=0;i<f.size();i++) f[i][1] = f[i][0] * f[i][1];
  std::cout<<Integrate(f,emin,emax)<<std::endl;
  return Integrate(f,emin,emax);
}

double Utils::EnergyContentFast(vector< vector<double> > f) {
  double sum = 0.;  
  for(unsigned int i=1;i<f.size();i++) 
    sum += f[i][0] * f[i][1] * (f[i][0] - f[i-1][0]);
  return sum;
}

///**
// * get Poissonian variate (just wrapped from TRandom))
// */
// int Utils::PoissonianRandom(double mean) {
//  int rand = randomiser->Poisson(mean);
//  return rand;

//}
/**
 * Integration function using the GSL spline integrator functionality
 *
 */
double Utils::Integrate(vector< vector<double> > f, double xmin, double xmax) {
  if(!f.size()) {
    cout << "Utils::Integrate: function vector empty! "
    "Exiting & returning 0." << endl;
    return 0.;
  }
  if(!xmin && !xmax) {
    xmin = f[0][0];
    xmax = f[f.size()-1][0];
  }
  if( xmin < f[0][0]) xmin = f[0][0];
  if( xmax > f[f.size()-1][0]) xmax = f[f.size()-1][0];

  int size = (int)f.size();
  double x[size];
  double y[size];
  double ymax = -1.e-100;
  for (unsigned int i=0;i<f.size();i++) {
    x[i] = f[i][0];
    y[i] = f[i][1];
    if(x[i] > xmin && x[i] < xmax && y[i] > ymax) ymax = y[i];
  }
  gsl_spline *lookup = gsl_spline_alloc(gsl_interp_linear, size);
  gsl_spline_init(lookup, x, y, size);
  gsl_interp_accel *a = gsl_interp_accel_alloc();
  double integral = 0.;
  int errcode = gsl_spline_eval_integ_e(lookup, xmin, xmax, a, &integral);
  if(errcode) {
    cout << "Utils::Integrate: Someting went wrong in the integration!"
            "Errorcode " << errcode << ". Returning 0. value. " << endl;
    integral = 0.;
  }
  return integral;
}

vector< vector< double> > Utils::IntegratedProfile(vector< vector<double> > f) {
  vector< vector<double> > v;
  double x0 = f[0][0];
  for(unsigned int i=1;i<f.size();i++) {
    double x = f[i][0];
    double y = Integrate(f,x0,x);
    v.push_back(vector<double>());
    v[v.size()-1].push_back(x);
    v[v.size()-1].push_back(y);
  }
  return v;
}

gsl_spline *Utils::GSLsplineFromTwoDVector(vector< vector<double> > v) {
  int size = (int)v.size();
  double x[size];
  double y[size];
  for(int i=0;i<size;i++) {
    x[i] = v[i][0];
    y[i] = v[i][1];
  }
  gsl_spline *s = gsl_spline_alloc((const gsl_interp_type*)INTERMETH, size);
  gsl_spline_init(s, x, y, size);
  return s;
}


interp2d_spline *Utils::TwoDsplineFromTwoDVector(vector< vector<double> > v, 
                                                 double &xmin, double &xmax,
                                                 double &ymin, double &ymax) {

  vector<double> vx;
  vector<double> vy;
  // get x and y -dimensions
//  int xdim = 0;
//  int ydim = 0;
//  double x == v[0][0];
//  while(x == v[0][0]) 
//    x = v[]
//xdim++;
//  ydim = (int)v.size()/xdim;

//  std::cout<<xdim<<" "<<ydim<<std::endl;

  double az[(int)v.size()];
  for(unsigned int i=0;i<v.size();i++) {
      double x = v[i][0];
      double y = v[i][1];
      double z = v[i][2];

      if(!vy.size()) vy.push_back(y);
      else if (y != vy[vy.size()-1]) vy.push_back(y);

      if(vy.size()==1) vx.push_back(x);
      az[i] = z;
  }
  xmin = vx[0];
  xmax = vx[vx.size()-1];
  ymin = vy[0];
  ymax = vy[vy.size()-1];
  unsigned long xs = (unsigned long)vx.size();
  unsigned long ys = (unsigned long)vy.size();
  double ax[(int)xs];
  double ay[(int)ys];
  for(unsigned int i=0;i<vx.size();i++) ax[i] = vx[i];
  for(unsigned int i=0;i<vy.size();i++) ay[i] = vy[i];

//  for(unsigned int i=0;i<vx.size();i++) cout << ax[i] <<endl;
////  cout << " - - - - - - " << endl;
//  for(unsigned int i=0;i<vy.size();i++) cout << ay[i] <<endl;
////  cout << " - - - - - - " << endl;
//  for(unsigned int i=0;i<v.size();i++) cout << az[i] <<endl;
////  cout << " - - - - - - " << endl;
  interp2d_spline *s = interp2d_spline_alloc(interp2d_bilinear, xs, ys);
  interp2d_spline_init (s, ax, ay, az, xs, ys);
  return s;
}



double Utils::EvalSpline(double x, gsl_spline *s, gsl_interp_accel *a,
                         const char* t, int l) {
  double y = 0.;

  if (gsl_spline_eval_e(s, x, a, &y)) {
    cout << t << ",l." << l << ": Interpolation of lookup failed with GSL error"
              "code "
              << GSL_EDOM
              << ". Exiting!" <<endl;
    exit(1);
  }
  if (std::isnan(y) || std::isinf(y)) {
    cout << t << ",l." << l << ": value is " << y << ". Exiting!" <<endl;
    exit(1);
  }
  return y;
}

vector< vector<double> > Utils::SortTwoDVector(vector< vector<double> > v,
                                               int column) {
  if(!v.size()) {
    cout << "Utils::SortTwoDVector: function vector empty! "
    "Exiting!" << endl;
    return v;
  }
  if(column<0 || column>1) {
    cout << "Utils::SortTwoDVector: Sorting column must be 0 or 1! "
    "Exiting!" << endl;
    return v;
  }
  if(!column) sort (v.begin(), v.end(), sortcriterionfirstcolumn);
  else  sort (v.begin(), v.end(), sortcriterionsecondcolumn);

  return v;
}

vector< vector<double> > Utils::VectorAxisLogarithm(vector< vector<double> > v,
                                                    unsigned int column) {
  if(!v.size()) {
    cout << "Utils::VectorAxisLogarithm: Input vector empty."
            "Doing nothing." << endl;
    return v;
  }
  
  if(column >= v[0].size()) {
    cout << "Utils::VectorAxisLogarithm: Input vector has "<< v[0].size() << 
            "columns, but you want to logharithmise column " << column << "!"
            "Doing nothing." << endl;
    return v;
  }
  for(unsigned int i = 0; i < v.size() ; i++) 
      v[i][column] = log10(v[i][column]);
    
  return v;
}

vector< vector<double> > Utils::VectorAxisPow10(vector< vector<double> > v,
                                                    int column) {
  for(unsigned int i = 0; i < v.size() ; i++) {
    double xval = v[i][0];
    double yval = v[i][1];

    if(column == 0) 
      v[i][0] = pow(10.,xval);
    if(column == 1) 
      v[i][1] = pow(10.,yval);
    
  }

  return v;
}

void Utils::TwoDVectorPushBack(double x, double y,
                               vector< vector<double> > &v) {
  v.push_back(vector<double>());
  v[v.size()-1].push_back(x);
  v[v.size()-1].push_back(y);
  return;
}

/**
 * Private method that is used in the sorting of 2D vectors.
 */
bool sortcriterionfirstcolumn (vector<double> i,
                                             vector<double> j) {
  return (i[0]<j[0]);
}
bool sortcriterionsecondcolumn (vector<double> i,
                                              vector<double> j) {
 return (i[1]<j[1]);
}

void Utils::Clear2DVector(vector< vector<double> > &v) {
  for (unsigned int i = 0; i < v.size(); i++) v[i].clear();
  v.clear();
}

void Utils::SetInterpolationMethod(string intermeth) {
  if(!intermeth.compare("CUBICSPLINE")) INTERMETH = (gsl_interp_type*)gsl_interp_cspline;
  else if(!intermeth.compare("LINEAR")) INTERMETH = (gsl_interp_type*)gsl_interp_linear;
  else if(!intermeth.compare("AKIMA")) INTERMETH = (gsl_interp_type*)gsl_interp_akima;
  else if(!intermeth.compare("POLYNOMIAL")) INTERMETH = (gsl_interp_type*)gsl_interp_polynomial;
  else {
    std::cout<<"Utils::SetInterpolationMethod: " << intermeth
             <<" is not a valid interpolation type. Options: " << endl;
    std::cout<<"  'LINEAR','CUBICSPLINE','AKIMA','POLYNOMIAL'" << endl;
  }
  return;
}
