#include "Radiation.h"
#include "Utils.h"

Utils *fUtils;
Radiation *fRad;
vector<double> parameters;
vector<string> parameters_names;
vector<string> strings;
vector<string> strings_names;
int ReadParameterFile(const char* inputname);

int main(int argc, char **argv){

  fUtils = new Utils();
  fRad = new Radiation();
  const char* parameterfile = argv[1];
  if(ReadParameterFile(parameterfile)) return 0;

  double d = pc_to_cm*parameters[0];
  double n = parameters[1];
  double b = parameters[2];
  double t = parameters[3];
  double e = TeV_to_erg*1.e-12*parameters[4];
  vector< vector <double> >  electronvector,protonvector;
  fUtils->ReadIn(strings[0],electronvector);
  fUtils->ReadIn(strings[1],protonvector);
  /* set radiation object stuff */
  fRad->SetDistance(d);
  fRad->SetBField(b);
  fRad->SetAmbientDensity(n);
  fRad->AddThermalTargetPhotons(2.7,0.25*TeV_to_erg*1.e-12);//CMB is hard-coded!
  fRad->AddThermalTargetPhotons(t,e);
  fRad->SetElectrons(electronvector);
  fRad->SetProtons(protonvector);
  fRad->CalculateDifferentialPhotonSpectrum(300,1.e-20,1.e3);
  string outtag = strings[2];
  fUtils->WriteOut(fRad->GetProtonVector(),outtag+"_ProtonSpectrum.dat");
  fUtils->WriteOut(fRad->GetElectronVector(),outtag+"_ElectronSpectrum.dat");
  fUtils->WriteOut(fRad->GetTotalSED(),outtag+"_TotalSpectrum.dat");
  fUtils->WriteOut(fRad->GetICSED(),outtag+"_ICSpectrum.dat");
  fUtils->WriteOut(fRad->GetBremsstrahlungSED(),outtag+"_BremsSpectrum.dat");
  fUtils->WriteOut(fRad->GetSynchrotronSED(),outtag+"_SynchSpectrum.dat");
  fUtils->WriteOut(fRad->GetPPSED(),outtag+"_PPSpectrum.dat");
  fUtils->WriteOut(fRad->GetTotalTargetPhotonVector(),outtag+"_TargetPhotons.dat");
  return 0;
}

int ReadParameterFile(const char* inputname){
  parameters_names.push_back("Distance");
  parameters_names.push_back("AmbientDensity");
  parameters_names.push_back("BField");
  parameters_names.push_back("TRAD");
  parameters_names.push_back("edensRAD");
  strings_names.push_back("inputspectrumElectrons");
  strings_names.push_back("inputspectrumProtons");
  strings_names.push_back("outfile");
  return fUtils->ReadParameterFile(inputname, parameters_names, strings_names, parameters, strings);
}

