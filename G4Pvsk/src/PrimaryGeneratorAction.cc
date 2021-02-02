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
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun_1(0), 
  fParticleGun_2(0), 
  fParticleSource(0), 
  fEnvelopeBox(0),
  fCo60_integral(0)
{
  G4int n_particle = 1;
  fParticleGun_1  = new G4ParticleGun(n_particle); // 1333 keV
  fParticleGun_2  = new G4ParticleGun(n_particle); // 1173 keV

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="gamma");

  fParticleGun_1->SetParticleDefinition(particle);
  fParticleGun_1->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun_1->SetParticleEnergy(1.333*MeV);
  //fParticleGun_1->SetParticleEnergy(0.511*MeV);
  
  fParticleGun_2->SetParticleDefinition(particle);
  fParticleGun_2->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun_2->SetParticleEnergy(1.173*MeV);
  //fParticleGun_2->SetParticleEnergy(0.511*MeV);

  for( G4int idiv=0; idiv<nCo60_division; idiv++ ) {
    G4double theta = CLHEP::twopi * (G4double)idiv/(G4double)nCo60_division;
    fCo60_integral += 1. 
                       + 0.1020 * ( 3 * cos(theta) * cos(theta) )/2.
                       + 0.00907 * (35 * pow( cos(theta), 4 ) - 30 * cos(theta) * cos(theta) + 3 )/8.;
    //fCo60_integral += theta; 
    fCo60_dist[idiv] = fCo60_integral;
    G4cout << idiv << "-th bin: " << fCo60_dist[idiv] << "(integral: " << fCo60_integral << ")@" << theta << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun_1;
  delete fParticleGun_2;
  delete fParticleSource;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of ecah event
  //

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get Envelope volume
  // from G4LogicalVolumeStore.
  
  G4double envSizeXY = 0;
  G4double envSizeZ = 0;

  if (!fEnvelopeBox)
  {
    G4LogicalVolume* envLV
      = G4LogicalVolumeStore::GetInstance()->GetVolume("Envelope");
    if ( envLV ) fEnvelopeBox = dynamic_cast<G4Box*>(envLV->GetSolid());
  }

  if ( fEnvelopeBox ) {
    envSizeXY = fEnvelopeBox->GetXHalfLength()*2.;
    envSizeZ = fEnvelopeBox->GetZHalfLength()*2.;
  }  
  else  {
    G4ExceptionDescription msg;
    msg << "Envelope volume of box shape not found.\n"; 
    msg << "Perhaps you have changed geometry.\n";
    msg << "The gun will be place at the center.";
    G4Exception("PrimaryGeneratorAction::GeneratePrimaries()",
     "MyCode0002",JustWarning,msg);
  }

  // Vertex Position
  //G4double size = 0.1; 
  //G4double x0 = size * envSizeXY * (G4UniformRand()-0.5);
  //G4double y0 = size * envSizeXY * (G4UniformRand()-0.5);
  //G4double z0 = -0.5 * envSizeZ;
  
  //G4double r0 = 3 * (G4UniformRand()-0.5);
  //G4double ang0 = 2 * 3.14 * (G4UniformRand()-0.5);
  //G4double x0 = r0 * cos(ang0);
  //G4double y0 = r0 * sin(ang0);
  //G4double z0 = -0.5 * envSizeZ;
  
  G4double x0 = 0;
  G4double y0 = 0;
  G4double z0 = 0;
  
  // Co 60
  // gamma #1
  //G4ThreeVector firstBeam = G4ThreeVector( sin(theta0) * cos(phi0), sin(theta0) * sin(phi0), cos(theta0) );
  G4ThreeVector firstBeam = G4ThreeVector( 0, 0, 1 ); 
  
  // gamma #2
  //G4double theta = 3.14*180./180. + 3.14*1./180. * (G4UniformRand()-0.5);
  //G4double theta = CLHEP::pi*180./180.;
  G4double theta = G4UniformRand();
  for( G4int idiv=0; idiv<nCo60_division; idiv++ ) {
    if( (G4double)fCo60_dist[idiv]/(G4double)fCo60_integral > theta ) {
      theta = CLHEP::twopi * (G4double)idiv/(G4double)nCo60_division;
      break;
    } else if( idiv == nCo60_division-1 ) {
      theta = -1;
      G4cout << "Warning: angle is not assigned for gamma-rays properly, then it is set as 0." << G4endl;
    }
  }
  G4ThreeVector secondBeam = firstBeam;
  secondBeam.rotate( theta, G4ThreeVector(0,1,0) );

  G4double phi0   = CLHEP::twopi * G4UniformRand();
  firstBeam.rotate( phi0, G4ThreeVector(0,1,0) );
  secondBeam.rotate( phi0, G4ThreeVector(0,1,0) );

  G4double phi1   = CLHEP::twopi * G4UniformRand();
  G4double theta1 = CLHEP::pi * (G4UniformRand() - 0.5) * 30./180.;
  G4ThreeVector plane_axis = G4ThreeVector( cos(phi1)*sin(theta1), cos(theta1), sin(phi1)*sin(theta1));
  if( plane_axis != G4ThreeVector(0,1,0) ) {
    G4ThreeVector cross = plane_axis.cross( G4ThreeVector(0,1,0) );
    firstBeam.rotate( theta1, cross );
    secondBeam.rotate( theta1, cross );
  }
  
  //G4double theta1 = CLHEP::pi * (30./180. * G4UniformRand() - 15./180);
  //G4double phi1   = CLHEP::twopi * G4UniformRand();
  //G4ThreeVector cross = ( fabs(cos(theta0))==1. )? firstBeam.cross( G4ThreeVector(1,0,0) ):firstBeam.cross( G4ThreeVector(0,0,1) );
  //if( fabs(cos(theta))==1. ) secondBeam.rotate( theta, cross );
  //else (secondBeam.rotate( theta, cross )).rotate( phi, firstBeam );
  //G4ThreeVector secondBeam = G4ThreeVector(0,0,-1);
  //G4cout << "angle: " << secondBeam.angle( firstBeam ) << G4endl;

  fParticleGun_1->SetParticleMomentumDirection( firstBeam );
  fParticleGun_1->SetParticleEnergy(1.333*MeV);
  fParticleGun_1->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  fParticleGun_1->GeneratePrimaryVertex(anEvent);
  fParticleGun_2->SetParticleMomentumDirection( secondBeam );
  fParticleGun_2->SetParticleEnergy(1.173*MeV);
  fParticleGun_2->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  fParticleGun_2->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

