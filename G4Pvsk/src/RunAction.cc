//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "DKAnalysis.hh"
#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "PvskDetectorConstruction.hh"
// #include "Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
: G4UserRunAction(),
  fEdep(0.),
  fEdep2(0.)
{ 
  // add new units for dose
  // 
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;
   
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray); 

  // Register accumulable to the accumulable manager
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->RegisterAccumulable(fEdep);
  accumulableManager->RegisterAccumulable(fEdep2); 


  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);     

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in B4Analysis.hh
  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories 
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
    // Note: merging ntuples is available only with Root output

  // Book histograms, ntuple
  //
  
  // Creating histograms
  // 1D
  //G4double E_range[2] = { 0.*MeV , 1.5*MeV };
  G4double E_range[2] = { 0.*keV , 1500*keV };
  //G4double D_range[2] = { -1.*mm , 10*mm };
  G4double D_range[2] = { -3.*mm , 3*mm };

  analysisManager->CreateH1("ev_pass","ev_pass", 2, -0.5, 1.5);
  analysisManager->CreateH1("beam_angle","beam_angle", 100, 0, 3.14);
  analysisManager->CreateH1("D1_gam_edep","D1_gam_edep", 150, E_range[0], E_range[1]);
  analysisManager->CreateH1("D1_eltrk_l","D1_eltrk_l", 150, 0., 0.3*mm);
  analysisManager->CreateH1("D1_gam_depth","D1_gam_depth", 110, D_range[0], D_range[1]);
  
  analysisManager->CreateH1("D2_gam_edep","D2_gam_edep", 150, E_range[0], E_range[1]);
  analysisManager->CreateH1("D2_eltrk_l","D2_eltrk_l", 150, 0., 0.3*mm);
  analysisManager->CreateH1("D2_gam_depth","D2_gam_depth",  110, -D_range[1], -D_range[0]);

  // 2D
  analysisManager->CreateH2("beam1_x_y","beam1_x_y", 100, -1, 1, 100, -1, 1 );
  analysisManager->CreateH2("D1_gam_i-nit_x_y","D1_gam_initxy", 120, -12*mm, 12*mm, 120, -12*mm, 12*mm);
  analysisManager->CreateH2("D1_gam_diff_x_y","D1_gam_diffxy", 100, -1*mm, 1*mm, 120, -1*mm, 1*mm);
  analysisManager->CreateH2("D1_gam_edep_depth","D1_gam_edep_depth", 150, E_range[0], E_range[1],  110, D_range[0], D_range[1]);
  analysisManager->CreateH2("D1_gam_edep_eltrk_L","D1_gam_edep_eltrk_L", 150, E_range[0], E_range[1], 100, 0, 0.5*mm);
  
  analysisManager->CreateH2("D2_gam_init_x_y","D2_gam_initxy", 120, -12*mm, 12*mm, 120, -12*mm, 12*mm);
  analysisManager->CreateH2("D2_gam_diff_x_y","D2_gam_diffxy", 100, -1*mm, 1*mm, 120, -1*mm, 1*mm);
  analysisManager->CreateH2("D2_gam_edep_depth","D2_gam_edep_depth", 150, E_range[0], E_range[1],  110, -D_range[1], -D_range[0]);
  analysisManager->CreateH2("D2_gam_edep_eltrk_L","D2_gam_edep_eltrk_L", 150, E_range[0], E_range[1], 100, 0, 0.5*mm);
  
  analysisManager->CreateH2("D1_D2_gam_edep","D1_D2_gam_edep", 150, E_range[0], E_range[1], 150, E_range[0], E_range[1]);
  
  // Creating ntuple
  //
  analysisManager->CreateNtuple("Perovskite", "D1 and D2");
  
  analysisManager->CreateNtupleDColumn("beam1_energy");
  analysisManager->CreateNtupleDColumn("beam1_dirx");
  analysisManager->CreateNtupleDColumn("beam1_diry");
  analysisManager->CreateNtupleDColumn("beam1_dirz");
  analysisManager->CreateNtupleDColumn("beam2_energy");
  analysisManager->CreateNtupleDColumn("beam2_dirx");
  analysisManager->CreateNtupleDColumn("beam2_diry");
  analysisManager->CreateNtupleDColumn("beam2_dirz");
  analysisManager->CreateNtupleDColumn("beam1_beam2_angle");
  
  analysisManager->CreateNtupleDColumn("D1_gam_edep");
  analysisManager->CreateNtupleDColumn("D1_gam_initX");
  analysisManager->CreateNtupleDColumn("D1_gam_initY");
  analysisManager->CreateNtupleDColumn("D1_gam_initZ");
  analysisManager->CreateNtupleDColumn("D1_gam_diffX");
  analysisManager->CreateNtupleDColumn("D1_gam_diffY");
  analysisManager->CreateNtupleDColumn("D1_gam_depth");
  analysisManager->CreateNtupleIColumn("D1_gam_ncompt");
  analysisManager->CreateNtupleIColumn("D1_gam_nphot");
  analysisManager->CreateNtupleDColumn("D1_eltrk_l");
  analysisManager->CreateNtupleDColumn("D1_eltrk_maxl");
  analysisManager->CreateNtupleDColumn("D1_eltrk_maxlz");

  analysisManager->CreateNtupleDColumn("D2_gam_edep");
  analysisManager->CreateNtupleDColumn("D2_gam_initX");
  analysisManager->CreateNtupleDColumn("D2_gam_initY");
  analysisManager->CreateNtupleDColumn("D2_gam_initZ");
  analysisManager->CreateNtupleDColumn("D2_gam_diffX");
  analysisManager->CreateNtupleDColumn("D2_gam_diffY");
  analysisManager->CreateNtupleDColumn("D2_gam_depth");
  analysisManager->CreateNtupleIColumn("D2_gam_ncompt");
  analysisManager->CreateNtupleIColumn("D2_gam_nphot");
  analysisManager->CreateNtupleDColumn("D2_eltrk_l");
  analysisManager->CreateNtupleDColumn("D2_eltrk_maxl");
  analysisManager->CreateNtupleDColumn("D2_eltrk_maxlz");
  
  analysisManager->CreateNtupleDColumn("D3_gam_edep");
  analysisManager->CreateNtupleDColumn("D3_gam_initX");
  analysisManager->CreateNtupleDColumn("D3_gam_initY");
  analysisManager->CreateNtupleDColumn("D3_gam_initZ");
  analysisManager->CreateNtupleDColumn("D3_gam_diffX");
  analysisManager->CreateNtupleDColumn("D3_gam_diffY");
  analysisManager->CreateNtupleDColumn("D3_gam_depth");
  analysisManager->CreateNtupleIColumn("D3_gam_ncompt");
  analysisManager->CreateNtupleIColumn("D3_gam_nphot");
  analysisManager->CreateNtupleDColumn("D3_eltrk_l");
  analysisManager->CreateNtupleDColumn("D3_eltrk_maxl");
  analysisManager->CreateNtupleDColumn("D3_eltrk_maxlz");
  
  analysisManager->CreateNtupleDColumn("D4_gam_edep");
  analysisManager->CreateNtupleDColumn("D4_gam_initX");
  analysisManager->CreateNtupleDColumn("D4_gam_initY");
  analysisManager->CreateNtupleDColumn("D4_gam_initZ");
  analysisManager->CreateNtupleDColumn("D4_gam_diffX");
  analysisManager->CreateNtupleDColumn("D4_gam_diffY");
  analysisManager->CreateNtupleDColumn("D4_gam_depth");
  analysisManager->CreateNtupleIColumn("D4_gam_ncompt");
  analysisManager->CreateNtupleIColumn("D4_gam_nphot");
  analysisManager->CreateNtupleDColumn("D4_eltrk_l");
  analysisManager->CreateNtupleDColumn("D4_eltrk_maxl");
  analysisManager->CreateNtupleDColumn("D4_eltrk_maxlz");
  
  analysisManager->CreateNtupleDColumn("D5_gam_edep");
  analysisManager->CreateNtupleDColumn("D5_gam_initX");
  analysisManager->CreateNtupleDColumn("D5_gam_initY");
  analysisManager->CreateNtupleDColumn("D5_gam_initZ");
  analysisManager->CreateNtupleDColumn("D5_gam_diffX");
  analysisManager->CreateNtupleDColumn("D5_gam_diffY");
  analysisManager->CreateNtupleDColumn("D5_gam_depth");
  analysisManager->CreateNtupleIColumn("D5_gam_ncompt");
  analysisManager->CreateNtupleIColumn("D5_gam_nphot");
  analysisManager->CreateNtupleDColumn("D5_eltrk_l");
  analysisManager->CreateNtupleDColumn("D5_eltrk_maxl");
  analysisManager->CreateNtupleDColumn("D5_eltrk_maxlz");
  
  analysisManager->CreateNtupleDColumn("D6_gam_edep");
  analysisManager->CreateNtupleDColumn("D6_gam_initX");
  analysisManager->CreateNtupleDColumn("D6_gam_initY");
  analysisManager->CreateNtupleDColumn("D6_gam_initZ");
  analysisManager->CreateNtupleDColumn("D6_gam_diffX");
  analysisManager->CreateNtupleDColumn("D6_gam_diffY");
  analysisManager->CreateNtupleDColumn("D6_gam_depth");
  analysisManager->CreateNtupleIColumn("D6_gam_ncompt");
  analysisManager->CreateNtupleIColumn("D6_gam_nphot");
  analysisManager->CreateNtupleDColumn("D6_eltrk_l");
  analysisManager->CreateNtupleDColumn("D6_eltrk_maxl");
  analysisManager->CreateNtupleDColumn("D6_eltrk_maxlz");
 
  analysisManager->CreateNtupleDColumn("D7_gam_edep");
  analysisManager->CreateNtupleDColumn("D7_gam_initX");
  analysisManager->CreateNtupleDColumn("D7_gam_initY");
  analysisManager->CreateNtupleDColumn("D7_gam_initZ");
  analysisManager->CreateNtupleDColumn("D7_gam_diffX");
  analysisManager->CreateNtupleDColumn("D7_gam_diffY");
  analysisManager->CreateNtupleDColumn("D7_gam_depth");
  analysisManager->CreateNtupleIColumn("D7_gam_ncompt");
  analysisManager->CreateNtupleIColumn("D7_gam_nphot");
  analysisManager->CreateNtupleDColumn("D7_eltrk_l");
  analysisManager->CreateNtupleDColumn("D7_eltrk_maxl");
  analysisManager->CreateNtupleDColumn("D7_eltrk_maxlz");
  
  analysisManager->CreateNtupleDColumn("D8_gam_edep");
  analysisManager->CreateNtupleDColumn("D8_gam_initX");
  analysisManager->CreateNtupleDColumn("D8_gam_initY");
  analysisManager->CreateNtupleDColumn("D8_gam_initZ");
  analysisManager->CreateNtupleDColumn("D8_gam_diffX");
  analysisManager->CreateNtupleDColumn("D8_gam_diffY");
  analysisManager->CreateNtupleDColumn("D8_gam_depth");
  analysisManager->CreateNtupleIColumn("D8_gam_ncompt");
  analysisManager->CreateNtupleIColumn("D8_gam_nphot");
  analysisManager->CreateNtupleDColumn("D8_eltrk_l");
  analysisManager->CreateNtupleDColumn("D8_eltrk_maxl");
  analysisManager->CreateNtupleDColumn("D8_eltrk_maxlz");
  
  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{ 
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  // reset accumulables to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  //
  G4String fileName = "Perovskite";
  analysisManager->OpenFile(fileName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // Merge accumulables 
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  // Compute dose = total energy deposit in a run and its variance
  //
  G4double edep  = fEdep.GetValue();
  G4double edep2 = fEdep2.GetValue();
  
  G4double rms = edep2 - edep*edep/nofEvents;
  if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;  

  const PvskDetectorConstruction* detectorConstruction
   = static_cast<const PvskDetectorConstruction*>
     (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4double mass = detectorConstruction->GetScoringVolume()->GetMass();
  G4double dose = edep/mass;
  G4double rmsDose = rms/mass;

  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const PrimaryGeneratorAction* generatorAction
   = static_cast<const PrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
  {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }
        
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  // Print
  //  
  if (IsMaster()) {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------";
  }
  else {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------";
  }
/*  
  G4cout
     << G4endl
     << " The run consists of " << nofEvents << " "<< runCondition
     << G4endl
     << " Cumulated dose per run, in scoring volume : " 
     << G4BestUnit(dose,"Dose") << " rms = " << G4BestUnit(rmsDose,"Dose")
     << G4endl
     << "------------------------------------------------------------"
     << G4endl
     << " EAbs : mean = " 
     << G4BestUnit(analysisManager->GetH1(0)->mean(), "Energy") 
     << " rms = " 
     << G4BestUnit(analysisManager->GetH1(0)->rms(),  "Energy") << G4endl;
*/
  // save histograms & ntuple
  //
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::AddEdep(G4double edep)
{
  fEdep  += edep;
  fEdep2 += edep*edep;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

