//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Feb 13 03:17:15 2021 by ROOT version 6.20/04
// from TTree Perovskite/D1 and D2
// found on file: ../Data/Perovskite.root
//////////////////////////////////////////////////////////

#ifndef Perovskite_h
#define Perovskite_h
//#define Perovskite_cxx

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class Perovskite {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        beam1_energy;
   Double_t        beam1_dirx;
   Double_t        beam1_diry;
   Double_t        beam1_dirz;
   Double_t        beam2_energy;
   Double_t        beam2_dirx;
   Double_t        beam2_diry;
   Double_t        beam2_dirz;
   Double_t        beam1_beam2_angle;
   Double_t        D1_gam_edep;
   Double_t        D1_gam_initX;
   Double_t        D1_gam_initY;
   Double_t        D1_gam_initZ;
   Double_t        D1_gam_diffX;
   Double_t        D1_gam_diffY;
   Double_t        D1_gam_depth;
   Int_t           D1_gam_ncompt;
   Int_t           D1_gam_nphot;
   Double_t        D1_eltrk_l;
   Double_t        D1_eltrk_maxl;
   Double_t        D1_eltrk_maxlz;
   Double_t        D2_gam_edep;
   Double_t        D2_gam_initX;
   Double_t        D2_gam_initY;
   Double_t        D2_gam_initZ;
   Double_t        D2_gam_diffX;
   Double_t        D2_gam_diffY;
   Double_t        D2_gam_depth;
   Int_t           D2_gam_ncompt;
   Int_t           D2_gam_nphot;
   Double_t        D2_eltrk_l;
   Double_t        D2_eltrk_maxl;
   Double_t        D2_eltrk_maxlz;
   Double_t        D3_gam_edep;
   Double_t        D3_gam_initX;
   Double_t        D3_gam_initY;
   Double_t        D3_gam_initZ;
   Double_t        D3_gam_diffX;
   Double_t        D3_gam_diffY;
   Double_t        D3_gam_depth;
   Int_t           D3_gam_ncompt;
   Int_t           D3_gam_nphot;
   Double_t        D3_eltrk_l;
   Double_t        D3_eltrk_maxl;
   Double_t        D3_eltrk_maxlz;
   Double_t        D4_gam_edep;
   Double_t        D4_gam_initX;
   Double_t        D4_gam_initY;
   Double_t        D4_gam_initZ;
   Double_t        D4_gam_diffX;
   Double_t        D4_gam_diffY;
   Double_t        D4_gam_depth;
   Int_t           D4_gam_ncompt;
   Int_t           D4_gam_nphot;
   Double_t        D4_eltrk_l;
   Double_t        D4_eltrk_maxl;
   Double_t        D4_eltrk_maxlz;
   Double_t        D5_gam_edep;
   Double_t        D5_gam_initX;
   Double_t        D5_gam_initY;
   Double_t        D5_gam_initZ;
   Double_t        D5_gam_diffX;
   Double_t        D5_gam_diffY;
   Double_t        D5_gam_depth;
   Int_t           D5_gam_ncompt;
   Int_t           D5_gam_nphot;
   Double_t        D5_eltrk_l;
   Double_t        D5_eltrk_maxl;
   Double_t        D5_eltrk_maxlz;
   Double_t        D6_gam_edep;
   Double_t        D6_gam_initX;
   Double_t        D6_gam_initY;
   Double_t        D6_gam_initZ;
   Double_t        D6_gam_diffX;
   Double_t        D6_gam_diffY;
   Double_t        D6_gam_depth;
   Int_t           D6_gam_ncompt;
   Int_t           D6_gam_nphot;
   Double_t        D6_eltrk_l;
   Double_t        D6_eltrk_maxl;
   Double_t        D6_eltrk_maxlz;
   Double_t        D7_gam_edep;
   Double_t        D7_gam_initX;
   Double_t        D7_gam_initY;
   Double_t        D7_gam_initZ;
   Double_t        D7_gam_diffX;
   Double_t        D7_gam_diffY;
   Double_t        D7_gam_depth;
   Int_t           D7_gam_ncompt;
   Int_t           D7_gam_nphot;
   Double_t        D7_eltrk_l;
   Double_t        D7_eltrk_maxl;
   Double_t        D7_eltrk_maxlz;
   Double_t        D8_gam_edep;
   Double_t        D8_gam_initX;
   Double_t        D8_gam_initY;
   Double_t        D8_gam_initZ;
   Double_t        D8_gam_diffX;
   Double_t        D8_gam_diffY;
   Double_t        D8_gam_depth;
   Int_t           D8_gam_ncompt;
   Int_t           D8_gam_nphot;
   Double_t        D8_eltrk_l;
   Double_t        D8_eltrk_maxl;
   Double_t        D8_eltrk_maxlz;

   // List of branches
   TBranch        *b_beam1_energy;   //!
   TBranch        *b_beam1_dirx;   //!
   TBranch        *b_beam1_diry;   //!
   TBranch        *b_beam1_dirz;   //!
   TBranch        *b_beam2_energy;   //!
   TBranch        *b_beam2_dirx;   //!
   TBranch        *b_beam2_diry;   //!
   TBranch        *b_beam2_dirz;   //!
   TBranch        *b_beam1_beam2_angle;   //!
   TBranch        *b_D1_gam_edep;   //!
   TBranch        *b_D1_gam_initX;   //!
   TBranch        *b_D1_gam_initY;   //!
   TBranch        *b_D1_gam_initZ;   //!
   TBranch        *b_D1_gam_diffX;   //!
   TBranch        *b_D1_gam_diffY;   //!
   TBranch        *b_D1_gam_depth;   //!
   TBranch        *b_D1_gam_ncompt;   //!
   TBranch        *b_D1_gam_nphot;   //!
   TBranch        *b_D1_eltrk_l;   //!
   TBranch        *b_D1_eltrk_maxl;   //!
   TBranch        *b_D1_eltrk_maxlz;   //!
   TBranch        *b_D2_gam_edep;   //!
   TBranch        *b_D2_gam_initX;   //!
   TBranch        *b_D2_gam_initY;   //!
   TBranch        *b_D2_gam_initZ;   //!
   TBranch        *b_D2_gam_diffX;   //!
   TBranch        *b_D2_gam_diffY;   //!
   TBranch        *b_D2_gam_depth;   //!
   TBranch        *b_D2_gam_ncompt;   //!
   TBranch        *b_D2_gam_nphot;   //!
   TBranch        *b_D2_eltrk_l;   //!
   TBranch        *b_D2_eltrk_maxl;   //!
   TBranch        *b_D2_eltrk_maxlz;   //!
   TBranch        *b_D3_gam_edep;   //!
   TBranch        *b_D3_gam_initX;   //!
   TBranch        *b_D3_gam_initY;   //!
   TBranch        *b_D3_gam_initZ;   //!
   TBranch        *b_D3_gam_diffX;   //!
   TBranch        *b_D3_gam_diffY;   //!
   TBranch        *b_D3_gam_depth;   //!
   TBranch        *b_D3_gam_ncompt;   //!
   TBranch        *b_D3_gam_nphot;   //!
   TBranch        *b_D3_eltrk_l;   //!
   TBranch        *b_D3_eltrk_maxl;   //!
   TBranch        *b_D3_eltrk_maxlz;   //!
   TBranch        *b_D4_gam_edep;   //!
   TBranch        *b_D4_gam_initX;   //!
   TBranch        *b_D4_gam_initY;   //!
   TBranch        *b_D4_gam_initZ;   //!
   TBranch        *b_D4_gam_diffX;   //!
   TBranch        *b_D4_gam_diffY;   //!
   TBranch        *b_D4_gam_depth;   //!
   TBranch        *b_D4_gam_ncompt;   //!
   TBranch        *b_D4_gam_nphot;   //!
   TBranch        *b_D4_eltrk_l;   //!
   TBranch        *b_D4_eltrk_maxl;   //!
   TBranch        *b_D4_eltrk_maxlz;   //!
   TBranch        *b_D5_gam_edep;   //!
   TBranch        *b_D5_gam_initX;   //!
   TBranch        *b_D5_gam_initY;   //!
   TBranch        *b_D5_gam_initZ;   //!
   TBranch        *b_D5_gam_diffX;   //!
   TBranch        *b_D5_gam_diffY;   //!
   TBranch        *b_D5_gam_depth;   //!
   TBranch        *b_D5_gam_ncompt;   //!
   TBranch        *b_D5_gam_nphot;   //!
   TBranch        *b_D5_eltrk_l;   //!
   TBranch        *b_D5_eltrk_maxl;   //!
   TBranch        *b_D5_eltrk_maxlz;   //!
   TBranch        *b_D6_gam_edep;   //!
   TBranch        *b_D6_gam_initX;   //!
   TBranch        *b_D6_gam_initY;   //!
   TBranch        *b_D6_gam_initZ;   //!
   TBranch        *b_D6_gam_diffX;   //!
   TBranch        *b_D6_gam_diffY;   //!
   TBranch        *b_D6_gam_depth;   //!
   TBranch        *b_D6_gam_ncompt;   //!
   TBranch        *b_D6_gam_nphot;   //!
   TBranch        *b_D6_eltrk_l;   //!
   TBranch        *b_D6_eltrk_maxl;   //!
   TBranch        *b_D6_eltrk_maxlz;   //!
   TBranch        *b_D7_gam_edep;   //!
   TBranch        *b_D7_gam_initX;   //!
   TBranch        *b_D7_gam_initY;   //!
   TBranch        *b_D7_gam_initZ;   //!
   TBranch        *b_D7_gam_diffX;   //!
   TBranch        *b_D7_gam_diffY;   //!
   TBranch        *b_D7_gam_depth;   //!
   TBranch        *b_D7_gam_ncompt;   //!
   TBranch        *b_D7_gam_nphot;   //!
   TBranch        *b_D7_eltrk_l;   //!
   TBranch        *b_D7_eltrk_maxl;   //!
   TBranch        *b_D7_eltrk_maxlz;   //!
   TBranch        *b_D8_gam_edep;   //!
   TBranch        *b_D8_gam_initX;   //!
   TBranch        *b_D8_gam_initY;   //!
   TBranch        *b_D8_gam_initZ;   //!
   TBranch        *b_D8_gam_diffX;   //!
   TBranch        *b_D8_gam_diffY;   //!
   TBranch        *b_D8_gam_depth;   //!
   TBranch        *b_D8_gam_ncompt;   //!
   TBranch        *b_D8_gam_nphot;   //!
   TBranch        *b_D8_eltrk_l;   //!
   TBranch        *b_D8_eltrk_maxl;   //!
   TBranch        *b_D8_eltrk_maxlz;   //!

   Perovskite(TTree *tree=0);
   virtual ~Perovskite();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Perovskite_cxx
Perovskite::Perovskite(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../Data/Perovskite.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../Data/Perovskite.root");
      }
      f->GetObject("Perovskite",tree);

   }
   Init(tree);
}

Perovskite::~Perovskite()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Perovskite::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Perovskite::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Perovskite::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("beam1_energy", &beam1_energy, &b_beam1_energy);
   fChain->SetBranchAddress("beam1_dirx", &beam1_dirx, &b_beam1_dirx);
   fChain->SetBranchAddress("beam1_diry", &beam1_diry, &b_beam1_diry);
   fChain->SetBranchAddress("beam1_dirz", &beam1_dirz, &b_beam1_dirz);
   fChain->SetBranchAddress("beam2_energy", &beam2_energy, &b_beam2_energy);
   fChain->SetBranchAddress("beam2_dirx", &beam2_dirx, &b_beam2_dirx);
   fChain->SetBranchAddress("beam2_diry", &beam2_diry, &b_beam2_diry);
   fChain->SetBranchAddress("beam2_dirz", &beam2_dirz, &b_beam2_dirz);
   fChain->SetBranchAddress("beam1_beam2_angle", &beam1_beam2_angle, &b_beam1_beam2_angle);
   fChain->SetBranchAddress("D1_gam_edep", &D1_gam_edep, &b_D1_gam_edep);
   fChain->SetBranchAddress("D1_gam_initX", &D1_gam_initX, &b_D1_gam_initX);
   fChain->SetBranchAddress("D1_gam_initY", &D1_gam_initY, &b_D1_gam_initY);
   fChain->SetBranchAddress("D1_gam_initZ", &D1_gam_initZ, &b_D1_gam_initZ);
   fChain->SetBranchAddress("D1_gam_diffX", &D1_gam_diffX, &b_D1_gam_diffX);
   fChain->SetBranchAddress("D1_gam_diffY", &D1_gam_diffY, &b_D1_gam_diffY);
   fChain->SetBranchAddress("D1_gam_depth", &D1_gam_depth, &b_D1_gam_depth);
   fChain->SetBranchAddress("D1_gam_ncompt", &D1_gam_ncompt, &b_D1_gam_ncompt);
   fChain->SetBranchAddress("D1_gam_nphot", &D1_gam_nphot, &b_D1_gam_nphot);
   fChain->SetBranchAddress("D1_eltrk_l", &D1_eltrk_l, &b_D1_eltrk_l);
   fChain->SetBranchAddress("D1_eltrk_maxl", &D1_eltrk_maxl, &b_D1_eltrk_maxl);
   fChain->SetBranchAddress("D1_eltrk_maxlz", &D1_eltrk_maxlz, &b_D1_eltrk_maxlz);
   fChain->SetBranchAddress("D2_gam_edep", &D2_gam_edep, &b_D2_gam_edep);
   fChain->SetBranchAddress("D2_gam_initX", &D2_gam_initX, &b_D2_gam_initX);
   fChain->SetBranchAddress("D2_gam_initY", &D2_gam_initY, &b_D2_gam_initY);
   fChain->SetBranchAddress("D2_gam_initZ", &D2_gam_initZ, &b_D2_gam_initZ);
   fChain->SetBranchAddress("D2_gam_diffX", &D2_gam_diffX, &b_D2_gam_diffX);
   fChain->SetBranchAddress("D2_gam_diffY", &D2_gam_diffY, &b_D2_gam_diffY);
   fChain->SetBranchAddress("D2_gam_depth", &D2_gam_depth, &b_D2_gam_depth);
   fChain->SetBranchAddress("D2_gam_ncompt", &D2_gam_ncompt, &b_D2_gam_ncompt);
   fChain->SetBranchAddress("D2_gam_nphot", &D2_gam_nphot, &b_D2_gam_nphot);
   fChain->SetBranchAddress("D2_eltrk_l", &D2_eltrk_l, &b_D2_eltrk_l);
   fChain->SetBranchAddress("D2_eltrk_maxl", &D2_eltrk_maxl, &b_D2_eltrk_maxl);
   fChain->SetBranchAddress("D2_eltrk_maxlz", &D2_eltrk_maxlz, &b_D2_eltrk_maxlz);
   fChain->SetBranchAddress("D3_gam_edep", &D3_gam_edep, &b_D3_gam_edep);
   fChain->SetBranchAddress("D3_gam_initX", &D3_gam_initX, &b_D3_gam_initX);
   fChain->SetBranchAddress("D3_gam_initY", &D3_gam_initY, &b_D3_gam_initY);
   fChain->SetBranchAddress("D3_gam_initZ", &D3_gam_initZ, &b_D3_gam_initZ);
   fChain->SetBranchAddress("D3_gam_diffX", &D3_gam_diffX, &b_D3_gam_diffX);
   fChain->SetBranchAddress("D3_gam_diffY", &D3_gam_diffY, &b_D3_gam_diffY);
   fChain->SetBranchAddress("D3_gam_depth", &D3_gam_depth, &b_D3_gam_depth);
   fChain->SetBranchAddress("D3_gam_ncompt", &D3_gam_ncompt, &b_D3_gam_ncompt);
   fChain->SetBranchAddress("D3_gam_nphot", &D3_gam_nphot, &b_D3_gam_nphot);
   fChain->SetBranchAddress("D3_eltrk_l", &D3_eltrk_l, &b_D3_eltrk_l);
   fChain->SetBranchAddress("D3_eltrk_maxl", &D3_eltrk_maxl, &b_D3_eltrk_maxl);
   fChain->SetBranchAddress("D3_eltrk_maxlz", &D3_eltrk_maxlz, &b_D3_eltrk_maxlz);
   fChain->SetBranchAddress("D4_gam_edep", &D4_gam_edep, &b_D4_gam_edep);
   fChain->SetBranchAddress("D4_gam_initX", &D4_gam_initX, &b_D4_gam_initX);
   fChain->SetBranchAddress("D4_gam_initY", &D4_gam_initY, &b_D4_gam_initY);
   fChain->SetBranchAddress("D4_gam_initZ", &D4_gam_initZ, &b_D4_gam_initZ);
   fChain->SetBranchAddress("D4_gam_diffX", &D4_gam_diffX, &b_D4_gam_diffX);
   fChain->SetBranchAddress("D4_gam_diffY", &D4_gam_diffY, &b_D4_gam_diffY);
   fChain->SetBranchAddress("D4_gam_depth", &D4_gam_depth, &b_D4_gam_depth);
   fChain->SetBranchAddress("D4_gam_ncompt", &D4_gam_ncompt, &b_D4_gam_ncompt);
   fChain->SetBranchAddress("D4_gam_nphot", &D4_gam_nphot, &b_D4_gam_nphot);
   fChain->SetBranchAddress("D4_eltrk_l", &D4_eltrk_l, &b_D4_eltrk_l);
   fChain->SetBranchAddress("D4_eltrk_maxl", &D4_eltrk_maxl, &b_D4_eltrk_maxl);
   fChain->SetBranchAddress("D4_eltrk_maxlz", &D4_eltrk_maxlz, &b_D4_eltrk_maxlz);
   fChain->SetBranchAddress("D5_gam_edep", &D5_gam_edep, &b_D5_gam_edep);
   fChain->SetBranchAddress("D5_gam_initX", &D5_gam_initX, &b_D5_gam_initX);
   fChain->SetBranchAddress("D5_gam_initY", &D5_gam_initY, &b_D5_gam_initY);
   fChain->SetBranchAddress("D5_gam_initZ", &D5_gam_initZ, &b_D5_gam_initZ);
   fChain->SetBranchAddress("D5_gam_diffX", &D5_gam_diffX, &b_D5_gam_diffX);
   fChain->SetBranchAddress("D5_gam_diffY", &D5_gam_diffY, &b_D5_gam_diffY);
   fChain->SetBranchAddress("D5_gam_depth", &D5_gam_depth, &b_D5_gam_depth);
   fChain->SetBranchAddress("D5_gam_ncompt", &D5_gam_ncompt, &b_D5_gam_ncompt);
   fChain->SetBranchAddress("D5_gam_nphot", &D5_gam_nphot, &b_D5_gam_nphot);
   fChain->SetBranchAddress("D5_eltrk_l", &D5_eltrk_l, &b_D5_eltrk_l);
   fChain->SetBranchAddress("D5_eltrk_maxl", &D5_eltrk_maxl, &b_D5_eltrk_maxl);
   fChain->SetBranchAddress("D5_eltrk_maxlz", &D5_eltrk_maxlz, &b_D5_eltrk_maxlz);
   fChain->SetBranchAddress("D6_gam_edep", &D6_gam_edep, &b_D6_gam_edep);
   fChain->SetBranchAddress("D6_gam_initX", &D6_gam_initX, &b_D6_gam_initX);
   fChain->SetBranchAddress("D6_gam_initY", &D6_gam_initY, &b_D6_gam_initY);
   fChain->SetBranchAddress("D6_gam_initZ", &D6_gam_initZ, &b_D6_gam_initZ);
   fChain->SetBranchAddress("D6_gam_diffX", &D6_gam_diffX, &b_D6_gam_diffX);
   fChain->SetBranchAddress("D6_gam_diffY", &D6_gam_diffY, &b_D6_gam_diffY);
   fChain->SetBranchAddress("D6_gam_depth", &D6_gam_depth, &b_D6_gam_depth);
   fChain->SetBranchAddress("D6_gam_ncompt", &D6_gam_ncompt, &b_D6_gam_ncompt);
   fChain->SetBranchAddress("D6_gam_nphot", &D6_gam_nphot, &b_D6_gam_nphot);
   fChain->SetBranchAddress("D6_eltrk_l", &D6_eltrk_l, &b_D6_eltrk_l);
   fChain->SetBranchAddress("D6_eltrk_maxl", &D6_eltrk_maxl, &b_D6_eltrk_maxl);
   fChain->SetBranchAddress("D6_eltrk_maxlz", &D6_eltrk_maxlz, &b_D6_eltrk_maxlz);
   fChain->SetBranchAddress("D7_gam_edep", &D7_gam_edep, &b_D7_gam_edep);
   fChain->SetBranchAddress("D7_gam_initX", &D7_gam_initX, &b_D7_gam_initX);
   fChain->SetBranchAddress("D7_gam_initY", &D7_gam_initY, &b_D7_gam_initY);
   fChain->SetBranchAddress("D7_gam_initZ", &D7_gam_initZ, &b_D7_gam_initZ);
   fChain->SetBranchAddress("D7_gam_diffX", &D7_gam_diffX, &b_D7_gam_diffX);
   fChain->SetBranchAddress("D7_gam_diffY", &D7_gam_diffY, &b_D7_gam_diffY);
   fChain->SetBranchAddress("D7_gam_depth", &D7_gam_depth, &b_D7_gam_depth);
   fChain->SetBranchAddress("D7_gam_ncompt", &D7_gam_ncompt, &b_D7_gam_ncompt);
   fChain->SetBranchAddress("D7_gam_nphot", &D7_gam_nphot, &b_D7_gam_nphot);
   fChain->SetBranchAddress("D7_eltrk_l", &D7_eltrk_l, &b_D7_eltrk_l);
   fChain->SetBranchAddress("D7_eltrk_maxl", &D7_eltrk_maxl, &b_D7_eltrk_maxl);
   fChain->SetBranchAddress("D7_eltrk_maxlz", &D7_eltrk_maxlz, &b_D7_eltrk_maxlz);
   fChain->SetBranchAddress("D8_gam_edep", &D8_gam_edep, &b_D8_gam_edep);
   fChain->SetBranchAddress("D8_gam_initX", &D8_gam_initX, &b_D8_gam_initX);
   fChain->SetBranchAddress("D8_gam_initY", &D8_gam_initY, &b_D8_gam_initY);
   fChain->SetBranchAddress("D8_gam_initZ", &D8_gam_initZ, &b_D8_gam_initZ);
   fChain->SetBranchAddress("D8_gam_diffX", &D8_gam_diffX, &b_D8_gam_diffX);
   fChain->SetBranchAddress("D8_gam_diffY", &D8_gam_diffY, &b_D8_gam_diffY);
   fChain->SetBranchAddress("D8_gam_depth", &D8_gam_depth, &b_D8_gam_depth);
   fChain->SetBranchAddress("D8_gam_ncompt", &D8_gam_ncompt, &b_D8_gam_ncompt);
   fChain->SetBranchAddress("D8_gam_nphot", &D8_gam_nphot, &b_D8_gam_nphot);
   fChain->SetBranchAddress("D8_eltrk_l", &D8_eltrk_l, &b_D8_eltrk_l);
   fChain->SetBranchAddress("D8_eltrk_maxl", &D8_eltrk_maxl, &b_D8_eltrk_maxl);
   fChain->SetBranchAddress("D8_eltrk_maxlz", &D8_eltrk_maxlz, &b_D8_eltrk_maxlz);
   Notify();
}

Bool_t Perovskite::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Perovskite::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Perovskite::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Perovskite_cxx
