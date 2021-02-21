#ifndef _HIST_MANAGER_H_
#define _HIST_MANAGER_H_

#include <vector>

#include <TH1.h>
#include <TH2.h>

class HistManager
{
  public:
  HistManager();
  ~HistManager();
  void test();
  int Init();
  int FillSingleDet( std::vector<Int_t> opt, Double_t val );
  int DrawSingleDet( TString mode, Int_t num );
  // single hit hist
  // Det: 0-7, Matching: 0-4, PE: 0-2
  TH1F* h_Dhit_E[8][5][3];
  // coincidence hit hist
  TH1F* h_nhit;
  TH1F* h_Dhit_angle;
  TH1F* h_Dhit_angle_mB12;
  TH1F* h_Dhit_angle_mNO;
  TH2F* h_Dhit_Ecoin[4];
  TH2F* h_Dhit_Ecoin_mB12[4];
  TH2F* h_Dhit_Ecoin_mNO[4];
};

#endif
