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
  // Get hits collections #1
  fPerovskiteHCID = G4SDManager::GetSDMpointer()->GetCollectionID("PvskHitsCollection_1");
  auto PerovskiteHC_1 = GetHitsCollection(fPerovskiteHCID, event);
  auto PerovskiteHit_1 = (*PerovskiteHC_1)[0];
  
  // Get hits collections #2
  fPerovskiteHCID = G4SDManager::GetSDMpointer()->GetCollectionID("PvskHitsCollection_2");
  auto PerovskiteHC_2 = GetHitsCollection(fPerovskiteHCID, event);
  auto PerovskiteHit_2 = (*PerovskiteHC_2)[0];

  // Get hits collections #3
  fPerovskiteHCID = G4SDManager::GetSDMpointer()->GetCollectionID("PvskHitsCollection_3");
  auto PerovskiteHC_3 = GetHitsCollection(fPerovskiteHCID, event);
  auto PerovskiteHit_3 = (*PerovskiteHC_3)[0];

  // Get hits collections #4
  fPerovskiteHCID = G4SDManager::GetSDMpointer()->GetCollectionID("PvskHitsCollection_4");
  auto PerovskiteHC_4 = GetHitsCollection(fPerovskiteHCID, event);
  auto PerovskiteHit_4 = (*PerovskiteHC_4)[0];
  
  // Get hits collections #5
  fPerovskiteHCID = G4SDManager::GetSDMpointer()->GetCollectionID("PvskHitsCollection_5");
  auto PerovskiteHC_5 = GetHitsCollection(fPerovskiteHCID, event);
  auto PerovskiteHit_5 = (*PerovskiteHC_5)[0];
  
  // Get hits collections #6
  fPerovskiteHCID = G4SDManager::GetSDMpointer()->GetCollectionID("PvskHitsCollection_6");
  auto PerovskiteHC_6 = GetHitsCollection(fPerovskiteHCID, event);
  auto PerovskiteHit_6 = (*PerovskiteHC_6)[0];

  // Get hits collections #7
  fPerovskiteHCID = G4SDManager::GetSDMpointer()->GetCollectionID("PvskHitsCollection_7");
  auto PerovskiteHC_7 = GetHitsCollection(fPerovskiteHCID, event);
  auto PerovskiteHit_7 = (*PerovskiteHC_7)[0];
  
  // Get hits collections #8
  fPerovskiteHCID = G4SDManager::GetSDMpointer()->GetCollectionID("PvskHitsCollection_8");
  auto PerovskiteHC_8 = GetHitsCollection(fPerovskiteHCID, event);
  auto PerovskiteHit_8 = (*PerovskiteHC_8)[0];
  
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
    
  // get analysis manager
  ///*
  auto analysisManager = G4AnalysisManager::Instance();

  // distribution before pre-selection
  analysisManager->FillH1(0, 0);
  analysisManager->FillH1(1, particleGun_2->GetParticleMomentumDirection().angle(particleGun_1->GetParticleMomentumDirection()));
  analysisManager->FillH2(0, particleGun_1->GetParticleMomentumDirection().x(), particleGun_1->GetParticleMomentumDirection().y());
  if( true ) {
  /*
  if( PerovskiteHit_1->GetEdep() > 0 
      || PerovskiteHit_2->GetEdep() > 0
      || PerovskiteHit_3->GetEdep() > 0
      || PerovskiteHit_4->GetEdep() > 0
      || PerovskiteHit_5->GetEdep() > 0
      || PerovskiteHit_6->GetEdep() > 0
      || PerovskiteHit_7->GetEdep() > 0
      || PerovskiteHit_8->GetEdep() > 0
    ) { 
  */
  
    // fill histograms
    // 1D
    G4int ih1d = 2;
    analysisManager->FillH1(0, 1);ih1d++;
    analysisManager->FillH1(ih1d, PerovskiteHit_1->GetEdep());ih1d++;
    analysisManager->FillH1(ih1d, PerovskiteHit_1->GetTrackLength()); ih1d++;
    analysisManager->FillH1(ih1d, PerovskiteHit_1->GetGammaThickness()); ih1d++;
    
    analysisManager->FillH1(ih1d, PerovskiteHit_2->GetEdep());ih1d++;
    analysisManager->FillH1(ih1d, PerovskiteHit_2->GetTrackLength());ih1d++;
    analysisManager->FillH1(ih1d, PerovskiteHit_2->GetGammaThickness());ih1d++;
    
    // 2D
    G4int ih2d = 1;
    analysisManager->FillH2(ih2d, PerovskiteHit_1->GetGammaInX(), PerovskiteHit_1->GetGammaInY());ih1d++;
    analysisManager->FillH2(ih2d, PerovskiteHit_1->GetGammaDiffX(), PerovskiteHit_1->GetGammaDiffY());ih1d++;
    analysisManager->FillH2(ih2d, PerovskiteHit_1->GetEdep(), PerovskiteHit_1->GetGammaThickness());ih1d++;
    analysisManager->FillH2(ih2d, PerovskiteHit_1->GetEdep(), PerovskiteHit_1->GetEltrkMaxLz());ih1d++;
    
    analysisManager->FillH2(ih2d, PerovskiteHit_2->GetGammaInX(), PerovskiteHit_2->GetGammaInY());ih1d++;
    analysisManager->FillH2(ih2d, PerovskiteHit_2->GetGammaDiffX(), PerovskiteHit_2->GetGammaDiffY());ih1d++;
    analysisManager->FillH2(ih2d, PerovskiteHit_2->GetEdep(), PerovskiteHit_2->GetGammaThickness());ih1d++;
    analysisManager->FillH2(ih2d, PerovskiteHit_2->GetEdep(), PerovskiteHit_2->GetEltrkMaxLz());ih1d++;
    
    analysisManager->FillH2(ih2d, PerovskiteHit_1->GetEdep(), PerovskiteHit_2->GetEdep() );ih1d++;
  
    // fill ntuple
    G4int it = 0;
    analysisManager->FillNtupleDColumn(it, particleGun_1->GetParticleEnergy());it++;
    analysisManager->FillNtupleDColumn(it, particleGun_1->GetParticleMomentumDirection().x());it++;
    analysisManager->FillNtupleDColumn(it, particleGun_1->GetParticleMomentumDirection().y());it++;
    analysisManager->FillNtupleDColumn(it, particleGun_1->GetParticleMomentumDirection().z());it++;
    analysisManager->FillNtupleDColumn(it, particleGun_2->GetParticleEnergy());it++;
    analysisManager->FillNtupleDColumn(it, particleGun_2->GetParticleMomentumDirection().x());it++;
    analysisManager->FillNtupleDColumn(it, particleGun_2->GetParticleMomentumDirection().y());it++;
    analysisManager->FillNtupleDColumn(it, particleGun_2->GetParticleMomentumDirection().z());it++;
    analysisManager->FillNtupleDColumn(it, particleGun_2->GetParticleMomentumDirection().angle(particleGun_1->GetParticleMomentumDirection()));it++;
    
    G4double reso = 0.01;
    //G4double bg = 2.18;
    //reso = sqrt(PerovskiteHit_1->GetEdep()*1e6/bg)/(PerovskiteHit_1->GetEdep()*1e6/bg); 
    analysisManager->FillNtupleDColumn(it, G4RandGauss::shoot(PerovskiteHit_1->GetEdep(),PerovskiteHit_1->GetEdep()*reso));it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_1->GetGammaInX());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_1->GetGammaInY());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_1->GetGammaInZ());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_1->GetGammaDiffX());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_1->GetGammaDiffY());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_1->GetGammaThickness());it++;
    analysisManager->FillNtupleIColumn(it, PerovskiteHit_1->GetNCompt());it++;
    analysisManager->FillNtupleIColumn(it, PerovskiteHit_1->GetNPhot());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_1->GetTrackLength());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_1->GetEltrkMaxL());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_1->GetEltrkMaxLz());it++;
    
    //reso = sqrt(PerovskiteHit_2->GetEdep()*1e6/bg)/(PerovskiteHit_2->GetEdep()*1e6/bg); 
    analysisManager->FillNtupleDColumn(it, G4RandGauss::shoot(PerovskiteHit_2->GetEdep(),PerovskiteHit_2->GetEdep()*reso));it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_2->GetGammaInX());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_2->GetGammaInY());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_2->GetGammaInZ());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_2->GetGammaDiffX());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_2->GetGammaDiffY());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_2->GetGammaThickness());it++;
    analysisManager->FillNtupleIColumn(it, PerovskiteHit_2->GetNCompt());it++;
    analysisManager->FillNtupleIColumn(it, PerovskiteHit_2->GetNPhot());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_2->GetTrackLength());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_2->GetEltrkMaxL());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_2->GetEltrkMaxLz());it++;
    
    //reso = sqrt(PerovskiteHit_3->GetEdep()*1e6/bg)/(PerovskiteHit_3->GetEdep()*1e6/bg); 
    analysisManager->FillNtupleDColumn(it, G4RandGauss::shoot(PerovskiteHit_3->GetEdep(),PerovskiteHit_3->GetEdep()*reso));it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_3->GetGammaInX());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_3->GetGammaInY());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_3->GetGammaInZ());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_3->GetGammaDiffX());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_3->GetGammaDiffY());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_3->GetGammaThickness());it++;
    analysisManager->FillNtupleIColumn(it, PerovskiteHit_3->GetNCompt());it++;
    analysisManager->FillNtupleIColumn(it, PerovskiteHit_3->GetNPhot());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_3->GetTrackLength());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_3->GetEltrkMaxL());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_3->GetEltrkMaxLz());it++;
    
    //reso = sqrt(PerovskiteHit_4->GetEdep()*1e6/bg)/(PerovskiteHit_4->GetEdep()*1e6/bg); 
    analysisManager->FillNtupleDColumn(it, G4RandGauss::shoot(PerovskiteHit_4->GetEdep(),PerovskiteHit_4->GetEdep()*reso));it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_4->GetGammaInX());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_4->GetGammaInY());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_4->GetGammaInZ());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_4->GetGammaDiffX());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_4->GetGammaDiffY());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_4->GetGammaThickness());it++;
    analysisManager->FillNtupleIColumn(it, PerovskiteHit_4->GetNCompt());it++;
    analysisManager->FillNtupleIColumn(it, PerovskiteHit_4->GetNPhot());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_4->GetTrackLength());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_4->GetEltrkMaxL());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_4->GetEltrkMaxLz());it++;
    
    //reso = sqrt(PerovskiteHit_5->GetEdep()*1e6/bg)/(PerovskiteHit_5->GetEdep()*1e6/bg); 
    analysisManager->FillNtupleDColumn(it, G4RandGauss::shoot(PerovskiteHit_5->GetEdep(),PerovskiteHit_5->GetEdep()*reso));it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_5->GetGammaInX());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_5->GetGammaInY());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_5->GetGammaInZ());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_5->GetGammaDiffX());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_5->GetGammaDiffY());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_5->GetGammaThickness());it++;
    analysisManager->FillNtupleIColumn(it, PerovskiteHit_5->GetNCompt());it++;
    analysisManager->FillNtupleIColumn(it, PerovskiteHit_5->GetNPhot());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_5->GetTrackLength());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_5->GetEltrkMaxL());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_5->GetEltrkMaxLz());it++;
    
    //reso = sqrt(PerovskiteHit_6->GetEdep()*1e6/bg)/(PerovskiteHit_6->GetEdep()*1e6/bg); 
    analysisManager->FillNtupleDColumn(it, G4RandGauss::shoot(PerovskiteHit_6->GetEdep(),PerovskiteHit_6->GetEdep()*reso));it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_6->GetGammaInX());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_6->GetGammaInY());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_6->GetGammaInZ());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_6->GetGammaDiffX());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_6->GetGammaDiffY());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_6->GetGammaThickness());it++;
    analysisManager->FillNtupleIColumn(it, PerovskiteHit_6->GetNCompt());it++;
    analysisManager->FillNtupleIColumn(it, PerovskiteHit_6->GetNPhot());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_6->GetTrackLength());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_6->GetEltrkMaxL());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_6->GetEltrkMaxLz());it++;
    
    //reso = sqrt(PerovskiteHit_7->GetEdep()*1e6/bg)/(PerovskiteHit_7->GetEdep()*1e6/bg); 
    analysisManager->FillNtupleDColumn(it, G4RandGauss::shoot(PerovskiteHit_7->GetEdep(),PerovskiteHit_7->GetEdep()*reso));it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_7->GetGammaInX());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_7->GetGammaInY());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_7->GetGammaInZ());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_7->GetGammaDiffX());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_7->GetGammaDiffY());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_7->GetGammaThickness());it++;
    analysisManager->FillNtupleIColumn(it, PerovskiteHit_7->GetNCompt());it++;
    analysisManager->FillNtupleIColumn(it, PerovskiteHit_7->GetNPhot());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_7->GetTrackLength());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_7->GetEltrkMaxL());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_7->GetEltrkMaxLz());it++;
    
    //reso = sqrt(PerovskiteHit_8->GetEdep()*1e6/bg)/(PerovskiteHit_8->GetEdep()*1e6/bg); 
    analysisManager->FillNtupleDColumn(it, G4RandGauss::shoot(PerovskiteHit_8->GetEdep(),PerovskiteHit_8->GetEdep()*reso));it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_8->GetGammaInX());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_8->GetGammaInY());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_8->GetGammaInZ());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_8->GetGammaDiffX());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_8->GetGammaDiffY());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_8->GetGammaThickness());it++;
    analysisManager->FillNtupleIColumn(it, PerovskiteHit_8->GetNCompt());it++;
    analysisManager->FillNtupleIColumn(it, PerovskiteHit_8->GetNPhot());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_8->GetTrackLength());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_8->GetEltrkMaxL());it++;
    analysisManager->FillNtupleDColumn(it, PerovskiteHit_8->GetEltrkMaxLz());it++;
    
    analysisManager->AddNtupleRow();  
    //*/
  }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
