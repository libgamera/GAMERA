#include "Utils.h"

Utils::Utils(bool DRAWLOGO) {
  QUIETMODE = false;
  GAMERADESTROYEDTHECONSOLE = false;
  if (DRAWLOGO == true) DrawGamera();
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
* get n uniform random numbers in [x_min,x_max)
*/
vector<double> Utils::UniformRandom(double x_min, double x_max,int n) {
  vector<double> v;
  for(int i=0;i<n;i++) v.push_back(x_min + (x_max-x_min)*gsl_rng_uniform(r));
  return v;
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
* get power-law distributed variates TODO: Index == 1
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
vector<double> Utils::CustomFunctionRandom(vector< vector<double> > f, int n,
                                           double xmin, double xmax ) {

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

/**
 * Sample Variate that follows the input 2D-Vector of a 2D dimensional function
 */
vector< vector<double> > Utils::CustomFunctionRandom2D(vector< vector<double> > f, int n,
                                           double xmin, double xmax,
                                           double ymin, double ymax ) {

  vector< vector<double> > v;
  if(!f.size()) {
    cout << "Utils::CustomFunctionRandom: function vector empty! "
            "Exiting!" << endl;
    return v;
  }
  if(!xmin && !xmax) {
    xmin = f[0][0];
    xmax = f[f.size()-1][0];
  }
  if(!ymin && !ymax) {
    ymin = f[0][1];
    ymax = f[f.size()-1][1];
  }
  if(  xmin < f[0][0] || xmax > f[f.size()-1][0] ||  
       ymin < f[0][1] || ymax > f[f.size()-1][1] ) {
    cout << "Utils::CustomFunctionRandom2D: requested sampling range outside "
            "of boundaries req.: x-(" << xmin << "," << xmax << ") vs. "
            "avail.:(" << f[0][0] << "," << f[f.size()-1][0] << "),"
            "y-(" << ymin << "," << ymax << ") vs. "
            "avail.:(" << f[0][1] << "," << f[f.size()-1][1] << "). "
            "Returning emptyvector." << endl;
    return v;
  }
  int size = (int)f.size();
  double x[size];
  double y[size];
  double z[size];
  double zmax = -1.e-100;
  for (unsigned int i=0;i<f.size();i++) {
    x[i] = f[i][0];
    y[i] = f[i][1];
    z[i] = f[i][2];
    if(x[i] > xmin && x[i] < xmax &&
       y[i] > ymin && y[i] < ymax && z[i] > zmax) zmax = z[i];
  }
  double x0,x1,y0,y1;
  interp2d_spline *spline2d = TwoDsplineFromTwoDVector(f,x0,x1,y0,y1);
  gsl_interp_accel *xaccsp = gsl_interp_accel_alloc();
  gsl_interp_accel *yaccsp = gsl_interp_accel_alloc();

  for(int i=0;i<n;i++) {
    double x,y;
    while(1) {
      x = (xmax-xmin)*gsl_rng_uniform(r)+xmin;
      y = (ymax-ymin)*gsl_rng_uniform(r)+ymin;
      double u = gsl_rng_uniform(r);
      double val = 0.;
      val = interp2d_spline_eval(spline2d,x,y,xaccsp,yaccsp);
//      if (gsl_spline_eval_e(lookup, x, a, &val)) {
//        cout << "Utils::CustomFunctionRandom: Function interpolation "
//        "failed. Exiting!" << endl;
//        return v;
//      }
      if (std::isnan(val) || std::isinf(val)) {
        cout << "Utils::CustomFunctionRandom: Function interpolation returned "
             << val << "! Exiting! " << endl;
        return v;
      }
      if(u < val/zmax) break;
    }
    TwoDVectorPushBack(x,y,v);
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
  if(f.size() < 3) {
    cout << "Utils::Integrate: function vector has less than 3 bins, "
            "Interpolation not possible! Exiting & returning 0." << endl;
    return 0.;
  }
  if(f[0].size() != 2) {
    cout << "Utils::Integrate: Input is 2D-array. This one has only "
         << f[0].size() 
         << " dimensions. Exiting & returning 0." << endl;
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
  for (unsigned int i=0;i<f.size();i++) {
    x[i] = f[i][0];
    y[i] = f[i][1];
  }
  gsl_spline *lookup = gsl_spline_alloc(gsl_interp_linear, size);
  gsl_spline_init(lookup, x, y, size);
  gsl_interp_accel *a = gsl_interp_accel_alloc();
  double integral = 0.;
  int errcode = gsl_spline_eval_integ_e(lookup, xmin, xmax, a, &integral);
  if(errcode) {
    if (QUIETMODE == false) {
        cout << "Utils::Integrate: Someting went wrong in the integration!"
            "Errorcode " << errcode << ". Returning 0. value. " << endl;}
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
  xmin = vx[0];  xmax = vx[vx.size()-1];
  ymin = vy[0];  ymax = vy[vy.size()-1];
  unsigned long xs = (unsigned long)vx.size();
  unsigned long ys = (unsigned long)vy.size();
  double ax[(int)xs];
  double ay[(int)ys];
  for(unsigned int i=0;i<vx.size();i++) ax[i] = vx[i];
  for(unsigned int i=0;i<vy.size();i++) ay[i] = vy[i];

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
  for(unsigned int i = 0; i < v.size() ; i++) {
      if (v[i][column])
          v[i][column] = log10(v[i][column]);
      else
          v[i][column] = -100.;
  }
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
    if(column == -1) {
      v[i][0] = pow(10.,xval);
      v[i][1] = pow(10.,yval);
    }
  }

  return v;
}

vector<double> Utils::GetVectorMinMax(vector< vector<double> > v, unsigned int axis) {
                
  vector<double> extrema;       
  if(!v.size()) {
    cout<<"Utils::GetVectorMinMax: " 
        <<" vector empty. returning." << endl;
    return extrema;
    
  }                         
  if(axis > v[0].size()-1) {
    cout<<"Utils::GetVectorMinMax: "
             <<" too few columns: Vector has " << 
             v[0].size() << "columns, but you asked for the extrema of the " 
             <<axis;
             if (axis == 1) cout<<"st";
             else if (axis == 2) cout<<"nd";
             else cout<<"th";
             cout<<" column. returning." << endl;
    return extrema;
  }
  double max = -1e100;
  double min = 1e100;
  
  for(unsigned int i=0;i<v.size();i++) {
    if(v[i][axis] < min) min = v[i][axis];
    if(v[i][axis] > max) max = v[i][axis];
  }
  extrema.push_back(min);
  extrema.push_back(max);
  return extrema;
}

vector<unsigned int> Utils::GetVectorMinMaxIndices(vector< vector<double> > v, unsigned int axis) {
                
  vector<unsigned int> extrema_indices;       
  if(!v.size()) {
    cout<<"Utils::GetVectorMinMax: " 
        <<" vector empty. returning." << endl;
    return extrema_indices;
    
  }                         
  if(axis > v[0].size()-1) {
    cout<<"Utils::GetVectorMinMax: "
             <<" too few columns: Vector has " << 
             v[0].size() << "columns, but you asked for the extrema of the " 
             <<axis;
             if (axis == 1) cout<<"st";
             else if (axis == 2) cout<<"nd";
             else cout<<"th";
             cout<<" column. returning." << endl;
    return extrema_indices;
  }
  double max = -1e100;
  double min = 1e100;
  unsigned i_min = 0;
  unsigned i_max = 0;
  for(unsigned int i=0;i<v.size();i++) {
    if(v[i][axis] < min)  {
        min = v[i][axis];
        i_min = i;
    }
    if(v[i][axis] > max) {
        max = v[i][axis];
        i_max = i;
    }
  }
  extrema_indices.push_back(i_min);
  extrema_indices.push_back(i_max);
  return extrema_indices;
}


vector< vector<double> > Utils::TwoDVectorFabs(vector< vector<double> > v) {
  vector< vector<double> > v_out;
  for(unsigned int i = 0; i < v.size() ; i++) {
    TwoDVectorPushBack(v[i][0],fabs(v[i][1]),v_out);  
  }
  return v_out;
}

vector< vector<double> > Utils::RemoveZeroEntries(vector< vector<double> > v) {
  vector< vector<double> > v_out;
  for(unsigned int i = 0; i < v.size() ; i++) {
    if(v[i][1]) TwoDVectorPushBack(v[i][0],v[i][1],v_out);  
  }
  return v_out;
}



vector< vector<double> > Utils::ZipTwoOneDVectors(vector<double> x, 
                                                vector<double> y) {
  vector< vector<double> > zipped;
  if(x.size() != y.size()) {
    cout << "Utils::ZipTwoOneDVectors: 1D vector sizes don't match!"
            "First vector has length "<< x.size() << 
            ", second vector has lenght " << y.size() << "."
            " Returning empty vector." << endl;
  }   
  else {
    for(unsigned int i=0; i<x.size();i++) {
      TwoDVectorPushBack(x[i],y[i],zipped);
    }
  }
  return zipped;
}


void Utils::TwoDVectorPushBack(double x, double y,
                               vector< vector<double> > &v) {
  v.push_back(vector<double>());
  v[v.size()-1].push_back(x);
  v[v.size()-1].push_back(y);
  return;
}

void Utils::TwoDVectorPushBack(double x, double y, double z,
                               vector< vector<double> > &v) {
  v.push_back(vector<double>());
  v[v.size()-1].push_back(x);
  v[v.size()-1].push_back(y);
  v[v.size()-1].push_back(z);
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


vector< vector<double> > Utils::MeshgridToTwoDVector(vector<double> x, vector<double> y, 
                                                     vector< vector<double> > mesh) {
                             
    vector< vector<double > > v;                        
    if(!mesh.size()){
        cout<<"Utils::MeshgridToTwoDVector: mesh has size 0. Returning empty vector. " 
            << endl;
        return v;}
    if(!x.size() || !y.size()) {
        cout<<"Utils::MeshgridToTwoDVector: At least one of the axis vectors has size 0. "
              "Returning empty vector."<< endl;
        return v;}
    if(x.size() != mesh[0].size() && y.size() != mesh.size()) {
        cout<<"Utils::MeshgridToTwoDVector: x and y-axis vector sizes don't match "
              "dimensions of mesh grid. Returning empty vector."<< endl;
        return v;}
    if(x.size() != mesh[0].size()) {
        cout<<"Utils::MeshgridToTwoDVector: x-axis vector size doesn't match x-dimension "
              "of mesh grid. Returning empty vector."<< endl;
        return v;}
    if(y.size() != mesh.size()) {
        cout<<"Utils::MeshgridToTwoDVector: y-axis vector size doesn't match y-dimension "
              "of mesh grid. Returning empty vector."<< endl;
        return v;}
    for(unsigned int i=0; i < mesh.size(); i++) {
         for(unsigned int j=0; j < mesh[i].size(); j++) {
             TwoDVectorPushBack(x[j],y[i],mesh[i][j],v);  
         }
    }
    return v;
}

/**
 * Rotate 3d-space vector v relative to position p1 to position p2. This uses Rodrigues' rotation formula.
 * The result is relative to p2.
 */
vector<double> Utils::RotateVector(vector<double> v, vector<double> p1, vector<double> p2) {

    double p1_norm = sqrt(p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2]);
    double p2_norm = sqrt(p2[0]*p2[0]+p2[1]*p2[1]+p2[2]*p2[2]);

    // dot product to get rotation angle
    double angle = acos((p1[0]*p2[0]+p1[1]*p2[1]+p1[2]*p2[2])/p1_norm/p2_norm);
    if (!angle) return v;
    if (fabs((angle-pi)/pi)<1e-8) {
        cout<<"Utils::RotateVector: Rotation of vector by 180 degrees is not well defined. You should supply an intermediate step in order to fix the orientation of rotation in 3D space. Returning original vector."<< endl;
        return v;
    }
    double sine = sin(angle);
    double cosine_factor = 1.-cos(angle);

    // vector product p1xp2, defining plane of rotation  
    double k_norm = p1_norm * p2_norm * sin(angle);
    double k_0,k_1,k_2;
    k_0 = (p1[1]*p2[2]-p1[2]*p2[1])/k_norm;
    k_1 = (p1[2]*p2[0]-p1[0]*p2[2])/k_norm;
    k_2 = (p1[0]*p2[1]-p1[1]*p2[0])/k_norm;

    // Idendity matrix
    double I_00,I_01,I_02,I_10,I_11,I_12,I_20,I_21,I_22;
    I_00 = 1.; I_01 = 0.; I_02 = 0.;
    I_10 = 0.; I_11 = 1.; I_12 = 0.;
    I_20 = 0.; I_21 = 0.; I_22 = 1.;    
    
    // auxiliary matrix K
    double K_00,K_01,K_02,K_10,K_11,K_12,K_20,K_21,K_22;
    K_00 =   0.; K_01 = -k_2;  K_02 = k_1;
    K_10 =  k_2; K_11 =   0.; K_12 = -k_0;
    K_20 = -k_1; K_21 =  k_0; K_22 =   0.;

    
    // auxiliary matrix K^2
    double Ks_00,Ks_01,Ks_02,Ks_10,Ks_11,Ks_12,Ks_20,Ks_21,Ks_22;
    Ks_00 = K_00*K_00+K_01*K_10+K_02*K_20; 
    Ks_01 = K_00*K_01+K_01*K_11+K_02*K_21; 
    Ks_02 = K_00*K_02+K_01*K_12+K_02*K_22;

    Ks_10 = K_10*K_00+K_11*K_10+K_12*K_20;
    Ks_11 = K_10*K_01+K_11*K_11+K_12*K_21;
    Ks_12 = K_10*K_02+K_11*K_12+K_12*K_22;

    Ks_20 = K_20*K_00+K_21*K_10+K_22*K_20;
    Ks_21 = K_20*K_01+K_21*K_11+K_22*K_21;
    Ks_22 = K_20*K_02+K_21*K_12+K_22*K_22;

    // final rotation matrix R
    double R_00,R_01,R_02,R_10,R_11,R_12,R_20,R_21,R_22;
    R_00 = I_00 + sine*K_00 + cosine_factor*Ks_00;
    R_01 = I_01 + sine*K_01 + cosine_factor*Ks_01;
    R_02 = I_02 + sine*K_02 + cosine_factor*Ks_02;
    R_10 = I_10 + sine*K_10 + cosine_factor*Ks_10;
    R_11 = I_11 + sine*K_11 + cosine_factor*Ks_11;
    R_12 = I_12 + sine*K_12 + cosine_factor*Ks_12;
    R_20 = I_20 + sine*K_20 + cosine_factor*Ks_20;
    R_21 = I_21 + sine*K_21 + cosine_factor*Ks_21;
    R_22 = I_22 + sine*K_22 + cosine_factor*Ks_22;
    
    // multiply vector with rotation matrix to obtain result
    vector<double> v_rot;
    v_rot.push_back(R_00*v[0]+R_01*v[1]+R_02*v[2]);
    v_rot.push_back(R_10*v[0]+R_11*v[1]+R_12*v[2]);
    v_rot.push_back(R_20*v[0]+R_21*v[1]+R_22*v[2]);
                        
    return v_rot;
}
 



//void Utils::TwoDVectorToMeshgrid(vector< vector<double> > vec)//, vector<double> &x,
////                                vector<double> &y, vector< vector<double> > &mesh) {
//    {
//    vector<double> x;
//    vector<double> y;
//    vector< vector<double> > mesh;

//    double yval = vec[0][1];
//    mesh.push_back(vector<double>());
//    y.push_back(vec[0][1]);

//    for(unsigned int i=0;i<vec.size();i++) {
//        if(i==0) x.push_back(vec[i][0]);
//        if(vec[i][1] != yval) {
//            mesh.push_back(vector<double>());
//            yval = vec[i][1];
//            y.push_back(yval);      
//        }
//        mesh[mesh.size()-1].push_back(vec[i][2]);
//    }
//    vector<double> a;
//    a.push_bac
//    return;
//}


void Utils::SetInterpolationMethod(string intermeth) {
  if(!intermeth.compare("CUBICSPLINE")) INTERMETH = (gsl_interp_type*)gsl_interp_cspline;
  else if(!intermeth.compare("LINEAR")) INTERMETH = (gsl_interp_type*)gsl_interp_linear;
  else if(!intermeth.compare("AKIMA")) INTERMETH = (gsl_interp_type*)gsl_interp_akima;
  else if(!intermeth.compare("POLYNOMIAL")) INTERMETH = (gsl_interp_type*)gsl_interp_polynomial;
  else {
    cout<<"Utils::SetInterpolationMethod: " << intermeth
             <<" is not a valid interpolation type. Options: " << endl;
    cout<<"  'LINEAR','CUBICSPLINE','AKIMA','POLYNOMIAL'" << endl;
  }
  return;
}

double Utils::Gaussian1D(double x, double sigma, double mu, double norm) {
    double s = (x - mu) / sigma;
    double gaus = exp(-0.5 * s * s );
    return (norm) ? norm * gaus: gaus / sigma / sqrt(2.*pi);
}

double Utils::LogisticsFunction(double z, double h, double w) {
  double val  = (1. + exp(-2.*(fabs(z)-h)/w));
  return 1./val;
}
