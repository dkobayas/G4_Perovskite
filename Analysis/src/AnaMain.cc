#include <iostream>

#include <TFile.h>
#include <TTree.h>

#include "../include/Perovskite.h"
#include "../include/HistManager.h"

int AnaMain( TString infname = "../Data/Perovskite_Co60.root"){
  //std::cout << "test" << std::endl;
  TTree* intree;
  TFile* infile = new TFile(infname);
  infile->GetObject("Perovskite",intree);
  Perovskite* t_pv = new Perovskite( intree );
  HistManager* hman = new HistManager();
  hman->Init();
  
  const int n_pv = t_pv->fChain->GetEntries();
  std::cout << "infile: " << infname << "(entry:" << n_pv << ")" << std::endl;
  for( int i_pv=0; i_pv<n_pv; i_pv++ ) {
    t_pv->GetEntry(i_pv);
    // variable amtrixnization
    // beam dir
    Double_t beam_x[2] = {t_pv->beam1_dirx, t_pv->beam2_dirx};
    Double_t beam_y[2] = {t_pv->beam1_diry, t_pv->beam2_diry};
    Double_t beam_z[2] = {t_pv->beam1_dirz, t_pv->beam2_dirz};
    // hit position
    Double_t Dhit_x[8] = {t_pv->D1_gam_initX, t_pv->D2_gam_initX,
                          t_pv->D3_gam_initX, t_pv->D4_gam_initX,
                          t_pv->D5_gam_initX, t_pv->D6_gam_initX,
                          t_pv->D7_gam_initX, t_pv->D8_gam_initX};
    Double_t Dhit_y[8] = {t_pv->D1_gam_initY, t_pv->D2_gam_initY,
                          t_pv->D3_gam_initY, t_pv->D4_gam_initY,
                          t_pv->D5_gam_initY, t_pv->D6_gam_initY,
                          t_pv->D7_gam_initY, t_pv->D8_gam_initY};
    Double_t Dhit_z[8] = {t_pv->D1_gam_initZ, t_pv->D2_gam_initZ,
                          t_pv->D3_gam_initZ, t_pv->D4_gam_initZ,
                          t_pv->D5_gam_initZ, t_pv->D6_gam_initZ,
                          t_pv->D7_gam_initZ, t_pv->D8_gam_initZ};
    Double_t Dhit_E[8] = {t_pv->D1_gam_edep, t_pv->D2_gam_edep,
                          t_pv->D3_gam_edep, t_pv->D4_gam_edep,
                          t_pv->D5_gam_edep, t_pv->D6_gam_edep,
                          t_pv->D7_gam_edep, t_pv->D8_gam_edep};
    Int_t Dhit_PE[8] = {t_pv->D1_gam_nphot, t_pv->D2_gam_nphot,
                           t_pv->D3_gam_nphot, t_pv->D4_gam_nphot,
                           t_pv->D5_gam_nphot, t_pv->D6_gam_nphot,
                           t_pv->D7_gam_nphot, t_pv->D8_gam_nphot};

    // hit matching with beam
    Double_t beam_rxz[2];
    for(int i_beam=0; i_beam<2; i_beam++ ) {
      beam_rxz[i_beam] = (beam_z[i_beam]==0)? 1:beam_x[i_beam]/beam_z[i_beam];
    }

    Double_t Dhit_rxz[8];
    Int_t Dhit_beam_sign[8][2];
    Bool_t Dhit_beam_match[8][2];
    std::vector<Int_t> hit_det;
    std::vector<Int_t> hit_det_high;
    for(int i_det=0; i_det<8; i_det++ ) {
      if( Dhit_E[i_det] > 0.05 ) hit_det.push_back(i_det);
      if( Dhit_E[i_det] > 0.1 ) hit_det_high.push_back(i_det);
      Dhit_rxz[i_det] = (Dhit_z[i_det]==0)? -100:Dhit_x[i_det]/Dhit_z[i_det];
      for(int i_beam=0; i_beam<2; i_beam++ ) {
        Dhit_beam_sign[i_det][i_beam] = (Dhit_x[i_det]*beam_x[i_beam]>0)? 1:-1;
        if( Dhit_x[i_det]*beam_x[i_beam]==0 ) {
          Dhit_beam_sign[i_det][i_beam] = (Dhit_z[i_det]*beam_z[i_beam]>0)? 1:-1;
        }
        Dhit_beam_match[i_det][i_beam] = ( Dhit_beam_sign[i_det][i_beam]==1 
                                           && (beam_rxz[i_beam]-Dhit_rxz[i_det])<0.01
                                         );
      }
    }
    
    const Int_t n_hit_det = hit_det.size();
    hman->h_nhit->Fill(n_hit_det);

    for(int i_hit_det=0; i_hit_det<n_hit_det; i_hit_det++ ) {
      //Float_t D1_beam2_x = (t_pv->D1_gam_initX*t_pv->beam2_dirx>0)? 1:-1;
      int i_det = hit_det[i_hit_det];
      //std::cout << "at " << i_pv << ": " << beam_rxz[0] << ", " << Dhit_rxz[0] << std::endl;
      std::vector<Int_t> opt(4,0);
      opt[0] = i_det; // i_det
      if( Dhit_beam_match[i_det][0] && Dhit_beam_match[i_det][1]) opt[2] = 3;
      else if( Dhit_beam_match[i_det][0] ) opt[2] = 1;
      else if( Dhit_beam_match[i_det][1] ) opt[2] = 2;
      else opt[2] = 0;
      if( Dhit_PE[i_det] > 0 ) opt[3] = 1; // p.e. flag
      else opt[3] = 0;
    
      opt[1] = 0; // energy
      hman->FillSingleDet( opt, Dhit_E[i_det]);
      // back to back
      if( i_det>3 ) hman->h_Dhit_E_BB[i_det-4]->Fill(Dhit_E[i_det-4]);
      else hman->h_Dhit_E_BB[i_det+4]->Fill(Dhit_E[i_det+4]);
    }

    const Int_t n_hit_detH = hit_det_high.size();
    if( n_hit_detH==2 ) {
      Int_t comb = abs(hit_det_high[1]-hit_det_high[0])-1;
      Double_t angle[7] = {45., 90., 135., 180., 225., 270., 315. };
      Double_t weight[7] = {7., 6., 5., 4., 3., 2., 1. };
      hman->h_Dhit_angle->Fill(angle[comb],1./weight[comb]);
      if( ( Dhit_beam_match[hit_det_high[0]][0] && Dhit_beam_match[hit_det_high[1]][1] ) 
          || ( Dhit_beam_match[hit_det_high[0]][1] && Dhit_beam_match[hit_det_high[1]][0] )
        ) {
        hman->h_Dhit_angle_mB12->Fill(angle[comb],1./weight[comb]);
      } else {
        hman->h_Dhit_angle_mNO->Fill(angle[comb],1./weight[comb]);
      }
    }
  }

  //hman->DrawSingleDet();
  hman->h_Dhit_angle->Draw();
  hman->h_Dhit_angle_mB12->Draw("same");
  return 0;
}
