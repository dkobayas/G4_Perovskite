#include <TROOT.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

#include <TH1.h>
#include <TString.h>
#include <THStack.h>

#include "../include/HistManager.h"

HistManager::HistManager(){};

HistManager::~HistManager(){};

void HistManager::test()
{
  std::cout << "HistManager::test" << std::endl;
}

int HistManager::Init()
{
  std::cout << "HistManager::Hist Initialization..." << std::endl;

  // single hit hist
  for ( int i_det=0; i_det<8; i_det++ ) {
    h_Dhit_E[i_det] = new TH1F(Form("h_D%d_E",i_det+1), 
                               Form("h_D%d_E;Energy[MeV];# event",i_det+1),
                               150, 0, 1.5);
    h_Dhit_E_mB1[i_det] = new TH1F(Form("h_D%d_E_mB1",i_det+1), 
                                   Form("h_D%d_E_mB1;Energy[MeV];# event",i_det+1),
                                   150, 0, 1.5);
    h_Dhit_E_mB2[i_det] = new TH1F(Form("h_D%d_E_mB2",i_det+1), 
                                   Form("h_D%d_E_mB2;Energy[MeV];# event",i_det+1),
                                   150, 0, 1.5);
    h_Dhit_E_mNO[i_det] = new TH1F(Form("h_D%d_E_mNO",i_det+1), 
                                   Form("h_D%d_E_mNO;Energy[MeV];# event",i_det+1),
                                   150, 0, 1.5);
    h_Dhit_E_noPE[i_det] = new TH1F(Form("h_D%d_E_noPE",i_det+1), 
                                    Form("h_D%d_E_noPE;Energy[MeV];# event",i_det+1),
                                    150, 0, 1.5);
    h_Dhit_E_PE[i_det] = new TH1F(Form("h_D%d_E_PE",i_det+1), 
                                  Form("h_D%d_E_PE;Energy[MeV];# event",i_det+1),
                                  150, 0, 1.5);
    h_Dhit_E_BB[i_det] = new TH1F(Form("h_D%d_E_BB",i_det+1), 
                                  Form("h_D%d_E_BB;Energy[MeV];# event",i_det+1),
                                  150, 0, 1.5);
  }
  return 0;
}

int HistManager::FillSingleDet( std::vector<Int_t> opt, Double_t val )
{
  /* options
  [0]:i_det 1-8
  [1]:plot edep(0), initxyz(1,2,3)
  [2]:matching no(0), beam1(1), beam2(2), both(3)
  [3]:p.e. no(0), yes(1)
  */
  if( opt.size() != 4 ) return 1;

  h_Dhit_E[opt[0]]->Fill(val);
  // matching
  if( opt[2]==1 ) {
    h_Dhit_E_mB1[opt[0]]->Fill(val);
  } else if( opt[2]==2 ) {
    h_Dhit_E_mB2[opt[0]]->Fill(val);
  } else {
    h_Dhit_E_mNO[opt[0]]->Fill(val);
  }
  // p.e.
  if( opt[3]==0 ) {
    h_Dhit_E_noPE[opt[0]]->Fill(val);
  } else if( opt[3]==1 ) {
    h_Dhit_E_PE[opt[0]]->Fill(val);
  }  
  return 0;
}

int HistManager::DrawSingleDet() 
{
  THStack* hs = new THStack("hs","");
  /*
  h_Dhit_E_mNO[0]->SetFillColor(kRed+2);
  hs->Add(h_Dhit_E_mNO[0]);
  h_Dhit_E_mB2[0]->SetFillColor(kBlue+2);
  hs->Add(h_Dhit_E_mB2[0]);
  h_Dhit_E_mB1[0]->SetFillColor(kGreen+2);
  hs->Add(h_Dhit_E_mB1[0]);
  */
  h_Dhit_E_noPE[0]->SetFillColor(kRed+2);
  hs->Add(h_Dhit_E_noPE[0]);
  h_Dhit_E_PE[0]->SetFillColor(kBlue+2);
  hs->Add(h_Dhit_E_PE[0]);
  hs->Draw();
  return 0;
}
