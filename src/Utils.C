#include "Utils.h"

Utils::Utils(bool DRAWLOGO) {
  if (DRAWLOGO == true) DrawGamera();
  QUIETMODE = false;
  GAMERADESTROYEDTHECONSOLE = false;
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

/// THIS IS STUFF THAT NEED TO VANISH / BE CHANGED / FREED OF ROOT ///

///**
// * get linearly distributed variates
// */
// double Utils::LinearRandom(double slope, double x_min, double x_max) {
//  double rand_min = 0.5*slope*pow(x_min,2.);
//  double rand_max = 0.5*slope*pow(x_max,2.);
//  double rand = rand_min + (rand_max - rand_min)*randomiser->Rndm(0);
//  return sqrt(2.*rand/slope);
//}

///**
// * get power-law distributed variates
// */
// double Utils::PowerLawRandom(double index, double x_min, double x_max) {
//  double rand = -100.;
//  double k = 0.03;
//  while(rand<x_min || rand>x_max) {
//    rand=pow(((index-1.)/k)*randomiser->Rndm(0),-1./(index-1.));
//  }
//  return rand;
//}

///**
// * get gaussian distributed variate using the Box-Muller method
// */
// double Utils::GaussianRandom(double width, double offset) {
//  double rand1 = randomiser->Rndm(0);
//  double rand2 = randomiser->Rndm(0);
//  double random_gaussian = sqrt(-2.*log(rand1))*cos(2*pi*rand2);
//  return offset + width*random_gaussian;

//}
///**
// * get Poissonian variate (just wrapped from TRandom))
// */
// int Utils::PoissonianRandom(double mean) {
//  int rand = randomiser->Poisson(mean);
//  return rand;

//}

// double Utils::LiMaSignificance(int ON, int OFF, double ALPHA) {
//  double nON,nOFF,alpha,significance,sign,nTOT,p,tON,tOFF;

//  if(!ALPHA) return sqrt((double)ON);

//  double XS = (double)ON - (double)OFF*ALPHA;
//  if (XS > 0) {
//    nON = (double)ON;
//    nOFF = (double)OFF;
//    alpha = ALPHA;
//    sign = 1.;
//  }
//  else {
//    nON = (double)OFF;
//    nOFF = (double)ON;
//    alpha = 1./ALPHA;
//    sign = -1.;
//  }

//  if(!nON)  return sign*sqrt(2.*nOFF*log(1.+alpha));
//  if(!nOFF) return sign*sqrt(2.*nON*log((1.+alpha)/alpha));

//  p = alpha+1.;
//  nTOT = nON+nOFF;
//  tON = nON*log((p/alpha)*(nON/nTOT));
//  tOFF = nOFF*log(p*(nOFF/nTOT));

//  significance = sign*sqrt(2)*sqrt(fabs(tON+tOFF));

//  return significance;
//}

////void Utils::GetRolkeConfidenceIntervals(int ON, int OFF, double ALPHA,
///double CL, double &xsDOWN, double &xsUP) {
////  TRolke *rolke = new TRolke();
////  rolke->SetCL(CL);
////  rolke->SetPoissonBkgKnownEff(ON,OFF,1./ALPHA,1);
////  rolke->GetLimits(xsDOWN,xsUP);
////  delete rolke;
////  return;
////}

///**
// * return a random sign
// */
// double Utils::SignRandom() {
//  double rand = randomiser->Rndm(0);
//  if(rand< 0.5) return 1.;
//  else return -1.;
//}

///**
// * return a exponentially distributed variate
// */
// double Utils::ExponentialRandom(double ind_norm, double x_min, double x_max)
// {
//  ind_norm = -1./ind_norm;
//  double rand_min = (1./ind_norm)*(exp(ind_norm*x_min)-1.);
//  double rand_max = (1./ind_norm)*(exp(ind_norm*x_max)-1.);
//  double rand = rand_min + (rand_max-rand_min)*randomiser->Rndm(0);
//  return pow(ind_norm,-1.)*log((ind_norm)*rand+1.);
//}

// double Utils::IntegrateTGraph(TGraph *graph, double xmin, double xmax, int
// steps, bool logarithmic) {
//  if(xmin>=xmax) {
//    cout<<"Utils::IntegrateTGraph: xmin>=xmax ("<<xmin<<">="<<xmax<<").
// Returning 0!"<<endl;
//    return 0.;
//  }
//  return Integrate(TGraphToTF1(graph),xmin,xmax,steps);
//}

// TF1 *Utils::TGraphToTF1(TGraph *graph) {
//  if(WrappedTGraph!=graph) {
//    WRAPPINGHAPPENED = true;
//    WrappedTGraph=graph;
//    double xmin,xmax,y;
//    graph->GetPoint(0,xmin,y);
//    graph->GetPoint(graph->GetN()-1,xmax,y);
//    TGraphWrapped=new
// TF1("TGraphWrapped",this,&Utils::TGraphWrapper,xmin,xmax,1,"Utils","TGraphWrapped");
//  }
//  return TGraphWrapped;
//}
// double Utils::TGraphWrapper(double *X, double *par) {
////cout<<"->wrap: "<<WrappedTGraph->Eval(X[0])<<endl;
//  return WrappedTGraph->Eval(X[0]);
//}

// TGraph *Utils::MakeTGraphIntegratedProfile(TGraph *graph) {
//  TGraph *graphIntegrated = new TGraph();
//  if(!graph->GetN()) {
// cout<<"Utils::MakeTGraphIntegratedProfile: TGraph you want to integrate is
// empty. Returning empty TGraph."<<endl;
//    return graphIntegrated;
//  }
//  double x0,x,y,u;
//  graph->GetPoint(0,x0,y);
//  u=0.;
//  for(int i=1;i<graph->GetN();i++) {
//    graph->GetPoint(i,x,y);
//    u += IntegrateTGraph(graph,x0,x,10,false);
//    graphIntegrated->SetPoint(graphIntegrated->GetN(),x,u);
//    x0 = x;
//  }
//  return graphIntegrated;
//}

// vector< vector<double> > Utils::TGraphToVector(TGraph* g) {
//  double x,y;
//  vector< vector<double> > v;
//  for(int i=0;i<g->GetN();i++) {
//    g->GetPoint(i,x,y);
//    v.push_back(vector<double> ());
//    v[v.size()-1].push_back(x);
//    v[v.size()-1].push_back(y);
//  }
//  return v;
//}

// TGraph *Utils::VectorToTGraph(vector< vector<double> > v) {
//  TGraph *g = new TGraph();
//  for(unsigned int i=0;i<v.size()-1;i++) {
//    g->SetPoint(g->GetN(),v[i][0],v[i][1]);
//  }
//  double xx,yy;
//  for(int i=0;i<g->GetN();i++) g->GetPoint(i,xx,yy);
//  return g;
//}

// TGraph *Utils::TH1FToTGraph(TH1F *h, int logtolinx) {
//  TGraph *g = new TGraph();
//  double x,y;
//  for(int i=0;i<h->GetNbinsX()+1;i++) {
//    x = h->GetBinCenter(i);
//    y = h->GetBinContent(i);
//    if(!logtolinx) g->SetPoint(g->GetN(),x,y);
//    if(logtolinx==1) g->SetPoint(g->GetN(),pow(10.,x),y);
//  }
//  return g;
//}

// TGraphErrors *Utils::TH1FToTGraphErrors(TH1F *h) {
//  TGraphErrors *g = new TGraphErrors();
//  double x,y,yerr;
//  for(int i=0;i<h->GetNbinsX()+1;i++) {
//    x = h->GetBinCenter(i);
//    y = h->GetBinContent(i);
//    yerr = h->GetBinError(i);
//    g->SetPoint(g->GetN(),x,y);
//    g->SetPointError(g->GetN()-1,0.,yerr);
//  }
//  return g;
//}

// double Utils::GetSpectralIndex(TGraph *g, double Eref) {
//  if(!g->GetN()) {
//    cout<<"Utils::GetSpectralIndex: not spectral points to fit!
// Exiting..."<<endl;
//  }
//  TF1 *fitfunction = new TF1("ff","[0]*(x/[1])^[2]",0.,1.);
//  fitfunction->SetRange(Eref*0.1,Eref*10.);
//  fitfunction->SetParameter(0,g->Eval(Eref));
//  fitfunction->SetParameter(1,Eref);
//  fitfunction->SetParLimits(0,10.,1.);
//  fitfunction->SetParLimits(1,10.,1.);
//  g->Fit(fitfunction,"QWR");
//  return fitfunction->GetParameter(2);
//}

// void Utils::WriteOut(TGraph *sp, string outname) {
//
//  ofstream output_file(outname.c_str());
//  if(output_file.is_open())  cout<<">> Writing to "<<outname<<endl;
//  else  cout<<">> Problem with "<<outname<<endl;
//
//  for(unsigned int i=0;i<sp->GetN();i++) {
//    double E,N;
//    E=N=0.;
//    sp->GetPoint(i,E,N);
////    if(E<=0. || N<=0.) continue;
//    output_file<<E<<" "<<N<<endl;
//  }
//  output_file.close();
//}

// void Utils::WriteOut(TGraphErrors *sp, string outname) {
//
//  ofstream output_file(outname.c_str());
//  if(output_file.is_open())  cout<<">> Writing to "<<outname<<endl;
//  else  cout<<">> Problem with "<<outname<<endl;
//
//  for(unsigned int i=0;i<sp->GetN();i++) {
//    double E,N,dE,dN;
//    E=N=dE=dN=0.;
//    sp->GetPoint(i,E,N);
//    dE = sp->GetErrorX(i);
//    dN = sp->GetErrorY(i);
////    if(E<=0. || N<=0.) continue;
//    output_file<<E<<" "<<dE<<" "<<N<<" "<<dN<<endl;
//  }
//  output_file.close();
//}

// void Utils::WriteOut(TGraphAsymmErrors *sp, string outname) {
//
//  ofstream output_file(outname.c_str());
//  if(output_file.is_open())  cout<<">> Writing to "<<outname<<endl;
//  else  cout<<">> Problem with "<<outname<<endl;
//
//  for(unsigned int i=0;i<sp->GetN();i++) {
//    double E,N,dELow,dEHigh,dNLow,dNHigh;
//    E=N=dELow=dEHigh=dNLow=dNHigh=0.;
//    sp->GetPoint(i,E,N);
//    dELow = sp->GetErrorXlow(i);
//    dEHigh = sp->GetErrorXhigh(i);
//    dNLow = sp->GetErrorYlow(i);
//    dNHigh = sp->GetErrorYhigh(i);
////    if(E<=0. || N<=0.) continue;
//    output_file<<E<<" "<<dELow<<" "<<dEHigh<<" "<<N<<" "<<dNLow<<"
// "<<dNHigh<<endl;
//  }
//  output_file.close();
//}

///**
// * Makes a nicer ROOT style and a cool sky-mappy color palette (black -> blue
// -> yellow -> white).
// */
// void Utils::SetNiceStyle() {

//  gROOT->SetStyle("Plain");
//  gStyle->SetPalette(1);
//  gStyle->SetOptStat(0);
//  gStyle->SetLabelFont(42,"x");
//  gStyle->SetLabelFont(42,"y");
//  gStyle->SetLabelFont(42,"z");
//  gStyle->SetTitleFont(42,"x");
//  gStyle->SetTitleFont(42,"y");
//  gStyle->SetTitleFont(42,"z");
//  int number=5;
//  int nb = 100;
//  double red[5] = {0.,0.,1.,1.,1.};
//  double green[5] = {0.,0.,0.,1.,1.};
//  double blue[5] = {0.,1.,0.,0.,1.};
//  double s[5] = {0.,0.3,0.6,0.9,1.};
//  TColor::CreateGradientColorTable(number,s,red,green,blue,nb);
//  gStyle->SetNumberContours(100);

//}

// TH2D *Utils::GetRightFormatHistogram(TGraph *graph, string xtitle, string
// ytitle, string title, double lowmargin, double highmargin, double leftmargin,
// double rightmargin) {
//  double ymin=1.e50;
//  double ymax=-1.e50;
//  double x,y,xmin,xmax;
//  graph->GetPoint(0,xmin,y);
//  graph->GetPoint(graph->GetN()-1,xmax,y);
//  for(int i=0;i<graph->GetN();i++) {
//    graph->GetPoint(i,x,y);
//    if(y<ymin) ymin = y;
//    if(y>ymax) ymax = y;
//  }
//	ymin *= lowmargin;
//  ymax *= highmargin;
//	xmin *= leftmargin;
//  xmax *= rightmargin;
//  TH2D *h = new TH2D(title.c_str(),"",100,xmin,xmax,100,ymin,ymax);
//  h->GetXaxis()->SetTitle(xtitle.c_str());
//  h->GetYaxis()->SetTitle(ytitle.c_str());
//  return h;
//}
