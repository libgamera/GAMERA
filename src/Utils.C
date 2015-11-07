#include "Utils.h"

Utils::Utils(bool DRAWLOGO) {
  if (DRAWLOGO == true) DrawGamera();
  QUIETMODE = false;
  GAMERADESTROYEDTHECONSOLE = false;
  r = gsl_rng_alloc (gsl_rng_mt19937);
  gsl_rng_set(r,time(0));
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
  cout << "                               .==#W,                       "
       << endl;
  cout << "                             ,*''RBEB                       "
       << endl;
  cout << "                              $=gp-'M @                     "
       << endl;
  cout << "                             f  `-^;p$/W=                   "
       << endl;
  cout << "                               ,|5==*$Bp5BK@,               "
       << endl;
  cout << "                           ,=@w=>=@@$ZpB@E8BBBw             "
       << endl;
  cout << "                        ,4MM'^YU@`5Z3@EEEZBM , .E=@@,@p@K   "
       << endl;
  cout << "                ,     ;`  ,,@;E]BBWBWp@5W#UP@P;&;ET@5BBEBBp "
       << endl;
  cout << "              ,*pT-p/@ ',/|@p|@]BP'`!*|`;p5N@EBp@(L`(,^^EBB "
       << endl;
  cout << "              `*=B@/@BB@pp@5@@BM`b,@|@{Ep5E@BB@Bpp^`'''```  "
       << endl;
  cout << "                ''(Z@5BE``'*9y^  `*@bb]pEbE@ZBBU@Lp         "
       << endl;
  cout << "                 ^^']EMP `  `^w@|L|`^'$EEEEEZ@@EPP          "
       << endl;
  cout << "                    $E`E   ` /t u'*L=u@EbEb@UEEH@p          "
       << endl;
  cout << "                    $QL$.  ,,ppp,|uEp^@Z@@@@EBEEZp          "
       << endl;
  cout << "                    ]pp5p(4BEp3$5B@EE@B$@BBB@p5UpM          "
       << endl;
  cout << "                    jE@pBB5p5Bp@5bP;p{53BBBBDBB@            "
       << endl;
  cout << "                     8M `I5|5EEEBU@4b55]BBU5EE(Ey           "
       << endl;
  cout << "                     'b  ',|@EBEBUE6@@$BBB|p,/`$@B          "
       << endl;
  cout << "      ,,              !'``T#p6EBBUEpBp65BBE^$^@uZ6N,        "
       << endl;
  cout << " ,p -``:!pp,         !.,!@y58$EEBBEEEUpEBBM95P/^``$C        "
       << endl;
  cout << "        ^YEP;B,       4/@@E@B]H3BEE@EppBBM 4(*WQ,A!b        "
       << endl;
  cout << "           !' |Bw   , ,/ZE3EBBHEE6@@EEBM   E^Z`|3+@@        "
       << endl;
  cout << "             ' `e!h ,8EbEEEM,B`($EE@E6H    'Zp/ZEE@'        "
       << endl;
  cout << "              `!`.em@@(uB'``,`Y@EEEEEB      j@P5]EBw,       "
       << endl;
  cout << "                ` ,/'^U^ @ @P;ZZKEEB`         Z`u!E^,b@p@w  "
       << endl;
  cout << "           ,,=@uw@P,]EEbppppBBM^             '9Nq@@~(C'9w ' "
       << endl;
  cout << "            ]^ ]B` AM^                                      "
       << endl;
  cout << endl;
  cout << "      ▄▀▀▀▀▄    ▄▀▀█▄   ▄▀▀▄ ▄▀▄  ▄▀▀█▄▄▄▄  ▄▀▀▄▀▀▀▄  ▄▀▀█▄  "
       << endl;
  cout << "     █         ▐ ▄▀ ▀▄ █  █ ▀  █ ▐  ▄▀   ▐ █   █   █ ▐ ▄▀ ▀▄ "
       << endl;
  cout << "     █    ▀▄▄    █▄▄▄█ ▐  █    █   █▄▄▄▄▄  ▐  █▀▀█▀    █▄▄▄█ "
       << endl;
  cout << "     █     █ █  ▄▀   █   █    █    █    ▌   ▄▀    █   ▄▀   █ "
       << endl;
  cout << "     ▐▀▄▄▄▄▀ ▐ █   ▄▀  ▄▀   ▄▀    ▄▀▄▄▄▄   █     █   █   ▄▀  "
       << endl;
  cout << "     ▐         ▐   ▐   █    █     █    ▐   ▐     ▐   ▐   ▐   "
       << endl;
  cout << "                       ▐    ▐     ▐                          "
       << endl;
  cout << endl;
  cout << endl;

  GAMERADESTROYEDTHECONSOLE = true;
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
double Utils::LinearRandom(double slope, double x_min, double x_max) {
 double rand_min = 0.5*slope*pow(x_min,2.);
 double rand_max = 0.5*slope*pow(x_max,2.);
 double rand = rand_min + (rand_max - rand_min)*gsl_rng_uniform(r);
 return sqrt(2.*rand/slope);
}

/**
* get power-law distributed variates
*/
double Utils::PowerLawRandom(double index, double x_min, double x_max) {
 double rand = -100.;
 double k = 0.03;
 while(rand<x_min || rand>x_max) {
   rand=pow(((index-1.)/k)*gsl_rng_uniform(r),-1./(index-1.));
 }
 return rand;
}

/**
* get gaussian distributed variate using the Box-Muller method
*/
double Utils::GaussianRandom(double width, double offset) {
 double random_gaussian = sqrt(-2.*log(gsl_rng_uniform(r)))
                        * cos(2*pi*gsl_rng_uniform(r));
 return offset + width*random_gaussian;
}

/**
* return a random sign
*/
double Utils::SignRandom() {
 if(gsl_rng_uniform(r)< 0.5) return 1.;
 else return -1.;
}

/**
* return a exponentially distributed variate
*/
double Utils::ExponentialRandom(double ind_norm, double x_min, double x_max)
{
 ind_norm = -1./ind_norm;
 double rand_min = (1./ind_norm)*(exp(ind_norm*x_min)-1.);
 double rand_max = (1./ind_norm)*(exp(ind_norm*x_max)-1.);
 double rand = rand_min + (rand_max-rand_min)*gsl_rng_uniform(r);
 return pow(ind_norm,-1.)*log((ind_norm)*rand+1.);
}

///**
// * get Poissonian variate (just wrapped from TRandom))
// */
// int Utils::PoissonianRandom(double mean) {
//  int rand = randomiser->Poisson(mean);
//  return rand;

//}
