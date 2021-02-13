#ifndef _HIST_MANAGER_H_
#define _HIST_MANAGER_H_

#include <vector>

#include <TH1.h>

class HistManager
{
  public:
  HistManager();
  ~HistManager();
  void test();
  int Init();
  int FillSingleDet( std::vector<Int_t> opt, Double_t val );
  int DrawSingleDet();
  // single hit hist
  TH1F* h_Dhit_E[8];
  TH1F* h_Dhit_E_mB1[8];
  TH1F* h_Dhit_E_mB2[8];
  TH1F* h_Dhit_E_mNO[8];
  TH1F* h_Dhit_E_noPE[8];
  TH1F* h_Dhit_E_PE[8];
  TH1F* h_Dhit_E_BB[8];
  // coincidence hit hist
  TH1F* h_Dhit_angle;
  TH1F* h_Dhit_angle_mB12;
  TH1F* h_Dhit_angle_mNO;
};

#endif
