//#include "include/HistManager.h"

void Ana( TString infile = "../Data/Perovskite.root"){
  cout << "Input file: " << infile << endl;
  
  gROOT->LoadMacro("src/Perovskite.cc+");
  gROOT->LoadMacro("src/HistManager.cc+");
  gROOT->LoadMacro("src/AnaMain.cc++");
  
  return;
}
