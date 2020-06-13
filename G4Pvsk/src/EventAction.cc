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
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "PrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"
#include "DKAnalysis.hh"
#include "EventAction.hh"
#include "RunAction.hh"
#include "PvskHit.hh"
#include "PvskSD.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
  fEdep(0.),
  fPerovskiteHCID(-1)
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{    
  fEdep = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PvskHitsCollection* 
EventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  auto hitsCollection 
    = static_cast<PvskHitsCollection*>(
        event->GetHCofThisEvent()->GetHC(hcID));
  
  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID; 
    G4Exception("EventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }         

  return hitsCollection;
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::PrintEventStatistics(
                              G4double absoEdep, G4double absoTrackLength) const
{
  // print event statistics
  G4cout
     << "      - total energy: " 
     << std::setw(7) << G4BestUnit(absoEdep, "Energy") << G4endl
     << "- total track length: " 
     << std::setw(7) << G4BestUnit(absoTrackLength, "Length")
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{   
  const PrimaryGeneratorAction* generatorAction
   = static_cast<const PrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  const G4ParticleGun* particleGun_1 = generatorAction->GetParticleGun(1);
  const G4ParticleGun* particleGun_2 = generatorAction->GetParticleGun(2);
  // accumulate statistics in run action
  // Get hits collections IDs (only once)
  //if ( fPerovskiteHCID == -1 ) {
  //  fPerovskiteHCID 
  //    = G4SDManager::GetSDMpointer()->GetCollectionID("PvskHitsCollection_1");
  //}
  fPerovskiteHCID = G4SDManager::GetSDMpointer()->GetCollectionID("PvskHitsCollection_1");
  // Get hits collections #1
  auto PerovskiteHC_1 = GetHitsCollection(fPerovskiteHCID, event);
  auto PerovskiteHit_1 = (*PerovskiteHC_1)[0];
  
  fPerovskiteHCID = G4SDManager::GetSDMpointer()->GetCollectionID("PvskHitsCollection_2");
  // Get hits collections #2
  auto PerovskiteHC_2 = GetHitsCollection(fPerovskiteHCID, event);
  auto PerovskiteHit_2 = (*PerovskiteHC_2)[0];

  // Get hit with total values
  //auto PerovskiteHit = (*PerovskiteHC)[PerovskiteHC->entries()-1];
 
  // Print per event (modulo n)
  //
  auto eventID = event->GetEventID();
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    G4cout << "---> End of event: " << eventID << G4endl;     
    G4cout << "    - Beam energy: " << particleGun_1->GetParticleEnergy() << ", " 
                                    << particleGun_2->GetParticleEnergy() << G4endl 
           << " - Beam direction: " << particleGun_1->GetParticleMomentumDirection().x() << ", " 
                                    << particleGun_2->GetParticleMomentumDirection().x() 
    << G4endl; 
    G4cout
      << "            Collection: PvskHitsCollection_1" << G4endl
      << "           - injection: " 
      << std::setw(7) << G4BestUnit(PerovskiteHit_1->GetGammaInX(), "Length") << ", "
                      << G4BestUnit(PerovskiteHit_1->GetGammaInY(), "Length") << G4endl
      << "          - final diff: " 
      << std::setw(7) << G4BestUnit(PerovskiteHit_1->GetGammaDiffX(), "Length") << ", "
                      << G4BestUnit(PerovskiteHit_1->GetGammaDiffY(), "Length") << G4endl
      << "     - Gamma Thickness: " 
      << std::setw(7) << G4BestUnit(PerovskiteHit_1->GetGammaThickness(), "Length")
      << "   - (# compt, # phot): " 
      << std::setw(7) << PerovskiteHit_1->GetNCompt() << ", " << PerovskiteHit_1->GetNPhot() << G4endl
      << " - el track max length: "
      << std::setw(7) << G4BestUnit(PerovskiteHit_1->GetEltrkMaxL(), "Length") << ", "
                      << G4BestUnit(PerovskiteHit_1->GetEltrkMaxLz(), "Length")
      << "        - total energy: "
      << std::setw(7) << G4BestUnit(PerovskiteHit_1->GetEdep(), "Energy") << G4endl
      << "  - total track length: "
      << std::setw(7) << G4BestUnit(PerovskiteHit_1->GetTrackLength(), "Length")
      << G4endl
      << "            Collection: PvskHitsCollection_2" << G4endl
      << "           - injection: " 
      << std::setw(7) << G4BestUnit(PerovskiteHit_2->GetGammaInX(), "Length") << ", "
                      << G4BestUnit(PerovskiteHit_2->GetGammaInY(), "Length") << G4endl
      << "          - final diff: " 
      << std::setw(7) << G4BestUnit(PerovskiteHit_2->GetGammaDiffX(), "Length") << ", "
                      << G4BestUnit(PerovskiteHit_2->GetGammaDiffY(), "Length") << G4endl
      << "     - Gamma Thickness: " 
      << std::setw(7) << G4BestUnit(PerovskiteHit_2->GetGammaThickness(), "Length")
      << "   - (# compt, # phot): " 
      << std::setw(7) << PerovskiteHit_2->GetNCompt() << ", " << PerovskiteHit_2->GetNPhot() << G4endl
      << " - el track max length: "
      << std::setw(7) << G4BestUnit(PerovskiteHit_2->GetEltrkMaxL(), "Length") << ", "
                      << G4BestUnit(PerovskiteHit_2->GetEltrkMaxLz(), "Length")
      << "        - total energy: "
      << std::setw(7) << G4BestUnit(PerovskiteHit_2->GetEdep(), "Energy") << G4endl
      << "  - total track length: "
      << std::setw(7) << G4BestUnit(PerovskiteHit_2->GetTrackLength(), "Length")
 
      << G4endl;
  }
  
  fRunAction->AddEdep(fEdep);
  //if( PerovskiteHit_1->GetEdep() > 0 || PerovskiteHit_2->GetEdep() > 0 ) { 
  if( true ) { 

    // get analysis manager
    ///*
    auto analysisManager = G4AnalysisManager::Instance();
  
    // fill histograms
    // 1D
    analysisManager->FillH1(0, PerovskiteHit_1->GetEdep());
    analysisManager->FillH1(1, PerovskiteHit_1->GetTrackLength());
    analysisManager->FillH1(2, PerovskiteHit_1->GetGammaThickness());
    
    analysisManager->FillH1(3, PerovskiteHit_2->GetEdep());
    analysisManager->FillH1(4, PerovskiteHit_2->GetTrackLength());
    analysisManager->FillH1(5, PerovskiteHit_2->GetGammaThickness());
    
    // 2D
    analysisManager->FillH2(0, PerovskiteHit_1->GetGammaInX(), PerovskiteHit_1->GetGammaInY());
    analysisManager->FillH2(1, PerovskiteHit_1->GetGammaDiffX(), PerovskiteHit_1->GetGammaDiffY());
    analysisManager->FillH2(2, PerovskiteHit_1->GetEdep(), PerovskiteHit_1->GetGammaThickness());
    analysisManager->FillH2(3, PerovskiteHit_1->GetEdep(), PerovskiteHit_1->GetEltrkMaxLz());
    
    analysisManager->FillH2(4, PerovskiteHit_2->GetGammaInX(), PerovskiteHit_2->GetGammaInY());
    analysisManager->FillH2(5, PerovskiteHit_2->GetGammaDiffX(), PerovskiteHit_2->GetGammaDiffY());
    analysisManager->FillH2(6, PerovskiteHit_2->GetEdep(), PerovskiteHit_2->GetGammaThickness());
    analysisManager->FillH2(7, PerovskiteHit_2->GetEdep(), PerovskiteHit_2->GetEltrkMaxLz());
    
    analysisManager->FillH2(8, PerovskiteHit_1->GetEdep(), PerovskiteHit_2->GetEdep() );
  
    // fill ntuple
    analysisManager->FillNtupleDColumn(0, PerovskiteHit_1->GetEdep());
    analysisManager->FillNtupleDColumn(1, PerovskiteHit_1->GetGammaInX());
    analysisManager->FillNtupleDColumn(2, PerovskiteHit_1->GetGammaInY());
    analysisManager->FillNtupleDColumn(3, PerovskiteHit_1->GetGammaInZ());
    analysisManager->FillNtupleDColumn(4, PerovskiteHit_1->GetGammaDiffX());
    analysisManager->FillNtupleDColumn(5, PerovskiteHit_1->GetGammaDiffY());
    analysisManager->FillNtupleDColumn(6, PerovskiteHit_1->GetGammaThickness());
    analysisManager->FillNtupleIColumn(7, PerovskiteHit_1->GetNCompt());
    analysisManager->FillNtupleIColumn(8, PerovskiteHit_1->GetNPhot());
    analysisManager->FillNtupleDColumn(9, PerovskiteHit_1->GetTrackLength());
    analysisManager->FillNtupleDColumn(10, PerovskiteHit_1->GetEltrkMaxL());
    analysisManager->FillNtupleDColumn(11, PerovskiteHit_1->GetEltrkMaxLz());
    
    analysisManager->FillNtupleDColumn(12, PerovskiteHit_2->GetEdep());
    analysisManager->FillNtupleDColumn(13, PerovskiteHit_2->GetGammaInX());
    analysisManager->FillNtupleDColumn(14, PerovskiteHit_2->GetGammaInY());
    analysisManager->FillNtupleDColumn(15, PerovskiteHit_2->GetGammaInZ());
    analysisManager->FillNtupleDColumn(16, PerovskiteHit_2->GetGammaDiffX());
    analysisManager->FillNtupleDColumn(17, PerovskiteHit_2->GetGammaDiffY());
    analysisManager->FillNtupleDColumn(18, PerovskiteHit_2->GetGammaThickness());
    analysisManager->FillNtupleIColumn(19, PerovskiteHit_2->GetNCompt());
    analysisManager->FillNtupleIColumn(20, PerovskiteHit_2->GetNPhot());
    analysisManager->FillNtupleDColumn(21, PerovskiteHit_2->GetTrackLength());
    analysisManager->FillNtupleDColumn(22, PerovskiteHit_2->GetEltrkMaxL());
    analysisManager->FillNtupleDColumn(23, PerovskiteHit_2->GetEltrkMaxLz());
    
    analysisManager->FillNtupleDColumn(24, particleGun_1->GetParticleEnergy());
    analysisManager->FillNtupleDColumn(25, particleGun_1->GetParticleMomentumDirection().x());
    analysisManager->FillNtupleDColumn(26, particleGun_1->GetParticleMomentumDirection().y());
    analysisManager->FillNtupleDColumn(27, particleGun_1->GetParticleMomentumDirection().z());
    analysisManager->FillNtupleDColumn(28, particleGun_2->GetParticleEnergy());
    analysisManager->FillNtupleDColumn(29, particleGun_2->GetParticleMomentumDirection().x());
    analysisManager->FillNtupleDColumn(30, particleGun_2->GetParticleMomentumDirection().y());
    analysisManager->FillNtupleDColumn(31, particleGun_2->GetParticleMomentumDirection().z());
    analysisManager->FillNtupleDColumn(32, particleGun_2->GetParticleMomentumDirection().angle(particleGun_1->GetParticleMomentumDirection()));
    
    analysisManager->AddNtupleRow();  
    //*/
  }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
