#include <TROOT.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

#include <TH1.h>
#include <TMath.h>
#include <TString.h>
#include <TCanvas.h>
#include <THStack.h>
#include <TLegend.h>

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

  const Int_t nbin_E = 100;
  const Double_t bin_E[2] = {0.,3.};
  TString s_bm[5] = {"","_mb1","_mb2","_mb12","mbNO"};
  TString s_pe[3] = {"","_pe","_nope"};
  // single hit hist
  for ( int i_det=0; i_det<8; i_det++ ) {
    for ( int i_bm=0; i_bm<5; i_bm++ ) {
      for ( int i_pe=0; i_pe<3; i_pe++ ) {
        h_Dhit_E[i_det][i_bm][i_pe] 
        = new TH1F(Form("h_D%d_E",i_det+1)+s_bm[i_bm]+s_pe[i_pe], 
                   Form("h_D%d_E",i_det+1)+s_bm[i_bm]+s_pe[i_pe]+";Energy[MeV];# event",
                   nbin_E, bin_E[0], bin_E[1]);
      }
    }
  }
  const Int_t nbin_ang = 5;
  const Double_t bin_ang[2] = {TMath::Pi()*(-1./8.), TMath::Pi()*(9./8.)};
  // coincidence hit hist
  h_nhit = new TH1F("h_nhit", "h_nhit;# hits;# entries",
                    9, -0.5, 8.5);
  h_Dhit_angle = new TH1F("h_Dhit_angle", "h_Dhit_angle;angle[deg];# entries",
                          nbin_ang, bin_ang[0], bin_ang[1]);
  h_Dhit_angle_mB12 = new TH1F("h_Dhit_angle_mB12", "h_Dhit_angle_mB12;angle[deg];# entries",
                          nbin_ang, bin_ang[0], bin_ang[1]);
  h_Dhit_angle_mNO = new TH1F("h_Dhit_angle_mNO", "h_Dhit_angle_mNO;angle[deg];# entries",
                          nbin_ang, bin_ang[0], bin_ang[1]);

  for ( int i_co=0; i_co<4; i_co++ ) {
    h_Dhit_Ecoin[i_co] = new TH2F(Form("h_Dhit_Ecoin_%d",(i_co+1)*45),
                                  Form("h_Dhit_Ecoin_%d;Energy[MeV];Energy[MeV]",(i_co+1)*45),
                                  nbin_E, bin_E[0], bin_E[1],
                                  nbin_E, bin_E[0], bin_E[1]);
    h_Dhit_Ecoin_mB12[i_co] = new TH2F(Form("h_Dhit_Ecoin_%d_mB12",(i_co+1)*45),
                                  Form("h_Dhit_Ecoin_%d_mB12;Energy[MeV];Energy[MeV]",(i_co+1)*45),
                                  nbin_E, bin_E[0], bin_E[1],
                                  nbin_E, bin_E[0], bin_E[1]);
    h_Dhit_Ecoin_mNO[i_co] = new TH2F(Form("h_Dhit_Ecoin_%d_mNO",(i_co+1)*45),
                                  Form("h_Dhit_Ecoin_%d_mNO;Energy[MeV];Energy[MeV]",(i_co+1)*45),
                                  nbin_E, bin_E[0], bin_E[1],
                                  nbin_E, bin_E[0], bin_E[1]);
  }
  return 0;
}

int HistManager::FillSingleDet( std::vector<Int_t> opt, Double_t val )
{
  /* options
  [0]:i_det 1-8
  [1]:plot edep(0), initxyz(1,2,3)
  [2]:matching no(4), beam1(1), beam2(2), both(3)
  [3]:p.e. no(2), yes(1)
  */
  if( opt.size() != 4 ) return 1;

  h_Dhit_E[opt[0]][0][0]->Fill(val);
  h_Dhit_E[opt[0]][0][opt[3]]->Fill(val);
  h_Dhit_E[opt[0]][opt[2]][0]->Fill(val);
  h_Dhit_E[opt[0]][opt[2]][opt[3]]->Fill(val);
  
  return 0;
}

int HistManager::DrawSingleDet( TString mode="E_match", Int_t num=0 ) 
{
  THStack* hs = new THStack("hs","");
  TLegend* leg = new TLegend(0.6,0.6,0.9,0.9,"title");
  TCanvas* c_semp[8];
  TCanvas* c_sem[8];
  TCanvas* c_sep[8];
  TString title[2] = {"x title", "y title"};
  ///*
  if( mode=="E_match_pe") {
    if( num > -1 && num < 8 ) {
      c_semp[num] = new TCanvas(Form("c_semp%d",num+1), Form("c_semp%d",num+1), 800,700);
      h_Dhit_E[num][3][2]->SetFillColor(kRed+3);
      hs->Add(h_Dhit_E[num][3][2]);
      h_Dhit_E[num][3][1]->SetFillColor(kRed-3);
      hs->Add(h_Dhit_E[num][3][1]);
      h_Dhit_E[num][4][2]->SetFillColor(kYellow+3);
      hs->Add(h_Dhit_E[num][4][2]);
      h_Dhit_E[num][4][1]->SetFillColor(kYellow-3);
      hs->Add(h_Dhit_E[num][4][1]);
      h_Dhit_E[num][2][2]->SetFillColor(kBlue+3);
      hs->Add(h_Dhit_E[num][2][2]);
      h_Dhit_E[num][2][1]->SetFillColor(kBlue-3);
      hs->Add(h_Dhit_E[num][2][1]);
      h_Dhit_E[num][1][2]->SetFillColor(kGreen+3);
      hs->Add(h_Dhit_E[num][1][2]);
      h_Dhit_E[num][1][1]->SetFillColor(kGreen-3);
      hs->Add(h_Dhit_E[num][1][1]);
      title[0] = h_Dhit_E[num][3][0]->GetXaxis()->GetTitle();
      title[1] = h_Dhit_E[num][3][0]->GetYaxis()->GetTitle();
      leg->SetHeader(Form("Energy Det#%d", num+1));
      leg->AddEntry(h_Dhit_E[num][1][1], "Beam1 & P.E.", "f");
      leg->AddEntry(h_Dhit_E[num][1][2], "Beam1 & Comp.", "f");
      leg->AddEntry(h_Dhit_E[num][2][1], "Beam2 & P.E.", "f");
      leg->AddEntry(h_Dhit_E[num][2][2], "Beam2 Comp.", "f");
      leg->AddEntry(h_Dhit_E[num][4][1], "No & P.E.", "f");
      leg->AddEntry(h_Dhit_E[num][4][2], "No & Comp.", "f");
      leg->AddEntry(h_Dhit_E[num][3][1], "Both & P.E.", "f");
      leg->AddEntry(h_Dhit_E[num][3][2], "Both & Comp.", "f");
    }
  }
  if( mode=="E_match") {
    if( num > -1 && num < 8 ) {
      c_sem[num] = new TCanvas(Form("c_sem%d",num+1), Form("c_sem%d",num+1), 800,700);
      h_Dhit_E[num][3][0]->SetFillColor(kRed+2);
      hs->Add(h_Dhit_E[num][3][0]);
      h_Dhit_E[num][4][0]->SetFillColor(kYellow+2);
      hs->Add(h_Dhit_E[num][4][0]);
      h_Dhit_E[num][2][0]->SetFillColor(kBlue+2);
      hs->Add(h_Dhit_E[num][2][0]);
      h_Dhit_E[num][1][0]->SetFillColor(kGreen+2);
      hs->Add(h_Dhit_E[num][1][0]);
      title[0] = h_Dhit_E[num][3][0]->GetXaxis()->GetTitle();
      title[1] = h_Dhit_E[num][3][0]->GetYaxis()->GetTitle();
      leg->SetHeader(Form("Energy Det#%d", num+1));
      leg->AddEntry(h_Dhit_E[num][1][0], "Beam1 matching", "f");
      leg->AddEntry(h_Dhit_E[num][2][0], "Beam2 matching", "f");
      leg->AddEntry(h_Dhit_E[num][3][0], "Both matching", "f");
      leg->AddEntry(h_Dhit_E[num][4][0], "No matching", "f");
    }
  }
  if( mode=="E_pe") {
    if( num > -1 && num < 8 ) {
      c_sep[num] = new TCanvas(Form("c_sep%d",num+1), Form("c_sep%d",num+1), 800,700);
      h_Dhit_E[num][0][2]->SetFillColor(kRed+2);
      hs->Add(h_Dhit_E[num][0][2]);
      h_Dhit_E[num][0][1]->SetFillColor(kBlue+2);
      hs->Add(h_Dhit_E[num][0][1]);
      title[0] = h_Dhit_E[num][0][1]->GetXaxis()->GetTitle();
      title[1] = h_Dhit_E[num][0][1]->GetYaxis()->GetTitle();
      leg->SetHeader(Form("Energy Det#%d", num+1));
      leg->AddEntry(h_Dhit_E[num][0][2], "Compton", "f");
      leg->AddEntry(h_Dhit_E[num][0][1], "Photo-Electron", "f");
    }
  }
  hs->Draw();
  hs->GetXaxis()->SetTitle( title[0] );
  hs->GetYaxis()->SetTitle( title[1] );
  leg->SetFillStyle(0);
  leg->Draw();
  return 0;
}
