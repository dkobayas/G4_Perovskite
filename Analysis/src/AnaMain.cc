#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH2.h>

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
    Double_t beam_cos[2];
    for(int i_beam=0; i_beam<2; i_beam++ ) {
      beam_rxz[i_beam] = (beam_z[i_beam]==0)? 1:beam_x[i_beam]/beam_z[i_beam];
      beam_cos[i_beam] = beam_x[i_beam]/sqrt(pow(beam_x[i_beam],2)+pow(beam_z[i_beam],2));
    }

    Double_t Dhit_rxz[8];
    Double_t Dhit_cos[8];
    Int_t Dhit_beam_sign[8][2];
    Bool_t Dhit_beam_match[8][2];
    std::vector<Int_t> hit_det;
    std::vector<Int_t> hit_det_high;
    for(int i_det=0; i_det<8; i_det++ ) {
      if( Dhit_E[i_det] > 0.1 ) hit_det.push_back(i_det);
      if( Dhit_E[i_det] > 0.5 ) hit_det_high.push_back(i_det);
      Dhit_rxz[i_det] = (Dhit_z[i_det]==0)? -100:Dhit_x[i_det]/Dhit_z[i_det];
      Dhit_cos[i_det] = Dhit_x[i_det]/sqrt(pow(Dhit_x[i_det],2)+pow(Dhit_z[i_det],2));
      for(int i_beam=0; i_beam<2; i_beam++ ) {
        Dhit_beam_sign[i_det][i_beam] = (Dhit_x[i_det]*beam_x[i_beam]>0)? 1:-1;
        if( Dhit_x[i_det]*beam_x[i_beam]==0 ) {
          Dhit_beam_sign[i_det][i_beam] = (Dhit_z[i_det]*beam_z[i_beam]>0)? 1:-1;
        }
        Dhit_beam_match[i_det][i_beam] = ( Dhit_beam_sign[i_det][i_beam]==1 
                                           && fabs(beam_cos[i_beam]-Dhit_cos[i_det])<0.01
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
      // beam matching
      if( Dhit_beam_match[i_det][0] && Dhit_beam_match[i_det][1]) opt[2] = 3;
      else if( Dhit_beam_match[i_det][0] ) opt[2] = 1;
      else if( Dhit_beam_match[i_det][1] ) opt[2] = 2;
      else opt[2] = 4;
      // p.e.
      if( Dhit_PE[i_det] > 0 ) opt[3] = 1; // p.e. flag
      else opt[3] = 2;
    
      opt[1] = 0; // energy
      hman->FillSingleDet( opt, Dhit_E[i_det]);
      // back to back
      //if( i_det>3 ) hman->h_Dhit_E_BB[i_det-4]->Fill(Dhit_E[i_det-4]);
      //else hman->h_Dhit_E_BB[i_det+4]->Fill(Dhit_E[i_det+4]);
    }

    const Int_t n_hit_detH = hit_det_high.size();
    if( n_hit_detH>0 && n_hit_det==2 ) {
      Int_t comb = abs(hit_det[1]-hit_det[0])-1;
      Int_t comb_type = 0;
      if( comb == 0 || comb == 6 ) comb_type = 0;
      else if( comb == 1 || comb == 5 ) comb_type = 1;
      else if( comb == 2 || comb == 4 ) comb_type = 2;
      else if( comb == 3 ) comb_type = 3;
      Double_t angle[7] = {TMath::Pi()*(1./4.), TMath::Pi()*(1./2.), TMath::Pi()*(3./4.), TMath::Pi(), 
                           TMath::Pi()*(3./4.), TMath::Pi()*(1./2.), TMath::Pi()*(1./4.) };
      Double_t weight[7] = {7.*2., 6.*2., 5.*2., 4., 3.*2., 2.*2., 1.*2. };
      hman->h_Dhit_angle->Fill(angle[comb],1./weight[comb]);
      hman->h_Dhit_Ecoin[comb_type]->Fill(Dhit_E[hit_det[0]],Dhit_E[hit_det[1]]);
      if( ( Dhit_beam_match[hit_det[0]][0] && Dhit_beam_match[hit_det[1]][1] ) 
          || ( Dhit_beam_match[hit_det[0]][1] && Dhit_beam_match[hit_det[1]][0] )
        ) {
        hman->h_Dhit_angle_mB12->Fill(angle[comb],1./weight[comb]);
        hman->h_Dhit_Ecoin_mB12[comb_type]->Fill(Dhit_E[hit_det[0]],Dhit_E[hit_det[1]]);
//        cout << beam_cos[0] << ", " << beam_cos[1] << ", " <<  Dhit_cos[hit_det_high[0]] << ", " << Dhit_cos[hit_det_high[1]]<< endl;
      } else {
        hman->h_Dhit_angle_mNO->Fill(angle[comb],1./weight[comb]);
        hman->h_Dhit_Ecoin_mNO[comb_type]->Fill(Dhit_E[hit_det[0]],Dhit_E[hit_det[1]]);
      }
    }
  }

  hman->DrawSingleDet( "E_pe", 0 );
  hman->DrawSingleDet( "E_match", 1 );
  hman->DrawSingleDet( "E_match_pe", 0 );
  //hman->h_Dhit_angle->SetFillColor(kRed+2);
  //hman->h_Dhit_angle->Draw("hist");
  //hman->h_Dhit_angle_mB12->SetFillColor(kBlue+2);
  //hman->h_Dhit_angle_mB12->Draw("hist,same");
  return 0;
}
