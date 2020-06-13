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
/// \file PvskSD.cc
/// \brief Implementation of the PvskSD class

#include "PvskSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PvskSD::PvskSD(
                            const G4String& name, 
                            const G4String& hitsCollectionName,
                            G4int nofCells)
 : G4VSensitiveDetector(name),
   fHitsCollection(nullptr),
   fNofCells(nofCells)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PvskSD::~PvskSD() 
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PvskSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection 
    = new PvskHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce
  auto hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 

  // Create hits
  // fNofCells for cells + one more for total sums 
  for (G4int i=0; i<fNofCells+1; i++ ) {
    fHitsCollection->insert(new PvskHit());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool PvskSD::ProcessHits(G4Step* step, 
                                     G4TouchableHistory*)
{  
  auto touchable = (step->GetPreStepPoint()->GetTouchable());
    
  // Get calorimeter cell id 
  auto layerNumber = touchable->GetReplicaNumber(1);
  
  // Get hit accounting data for this cell
  auto hit = (*fHitsCollection)[layerNumber];
  if ( ! hit ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hit " << layerNumber; 
    G4Exception("PvskSD::ProcessHits()",
      "MyCode0004", FatalException, msg);
  }         

  // Get hit for total accounting
  auto hitTotal 
    = (*fHitsCollection)[fHitsCollection->entries()-1];

  // energy deposit
  auto edep = step->GetTotalEnergyDeposit();
  
  // step length
  G4double stepLength = 0.;
  
  if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
    stepLength = step->GetStepLength();
    G4StepPoint* prevtx = step->GetPreStepPoint();
    G4StepPoint* postvtx = step->GetPostStepPoint();
    if( prevtx && postvtx ) {
      hit->AddEltrkStep( stepLength, (postvtx->GetPosition()).getZ() - (prevtx->GetPosition()).getZ() );
    }
    //G4cout << "length += " << stepLength << "(" << step->GetTrack()->GetDefinition()->GetPDGCharge() << ")" << G4endl;
  }
  else if ( step->GetTrack()->GetDefinition()->GetParticleName() == "gamma") {
    G4StepPoint* prevtx = step->GetPreStepPoint();
    G4StepPoint* postvtx = step->GetPostStepPoint();
    if ( ( postvtx->GetProcessDefinedStep()->GetProcessName() == "compt" || postvtx->GetProcessDefinedStep()->GetProcessName() == "phot" )
         && prevtx->GetProcessDefinedStep() ) {
      G4String prevtx_proc = prevtx->GetProcessDefinedStep()->GetProcessName();
      G4String postvtx_proc = postvtx->GetProcessDefinedStep()->GetProcessName();
      G4ThreeVector prevtx_pos = prevtx->GetPosition();
      G4ThreeVector postvtx_pos = postvtx->GetPosition();
      //G4cout << prevtx_proc << "->" << postvtx_proc<< G4endl 
      //<< "Lz = " << fabs(postvtx_pos.getZ()-prevtx_pos.getZ()) << G4endl;
      if( prevtx_proc == "Transportation") {
        hit->SetInitialGamma( prevtx_pos.getX(), prevtx_pos.getY(), prevtx_pos.getZ() );
      }
      hit->AddGammaStep( postvtx_pos.getX()-prevtx_pos.getX(),
                              postvtx_pos.getY()-prevtx_pos.getY(),
                              postvtx_pos.getZ()-prevtx_pos.getZ(),
                              postvtx_proc );
    }
  }

  //if ( edep==0. && stepLength == 0. ) return false;      

  
  // Add values
  hit->Add(edep, stepLength);
  hitTotal->Add(edep, stepLength); 
      
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PvskSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) { 
     auto nofHits = fHitsCollection->entries();
     //G4cout
     //  << G4endl 
     //  << "-------->Hits Collection: in this event they are " << nofHits 
     //  << " hits in the tracker chambers: " << G4endl;
     for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
