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
/// \file PvskHit.cc
/// \brief Implementation of the PvskHit class

#include "PvskHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

G4ThreadLocal G4Allocator<PvskHit>* PvskHitAllocator = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PvskHit::PvskHit()
 : G4VHit(),
   fEdep(0.),
   fTrackLength(0.),
   fEltrkMaxL(-99999.),
   fEltrkMaxLz(-99999.),
   fGammaThickness(-99999.),
   fGammaInX(-99999.),
   fGammaInY(-99999.),
   fGammaInZ(-99999.),
   fGammaDiffX(-99999.),
   fGammaDiffY(-99999.),
   fnCompt(0),
   fnPhot(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PvskHit::~PvskHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PvskHit::PvskHit(const PvskHit& right)
  : G4VHit()
{
  fEdep           = right.fEdep;
  fTrackLength    = right.fTrackLength;
  fEltrkMaxL      = right.fEltrkMaxL;
  fEltrkMaxLz     = right.fEltrkMaxLz;
  fGammaInX       = right.fGammaInX;
  fGammaInY       = right.fGammaInY;
  fGammaInZ       = right.fGammaInZ;
  fGammaThickness = right.fGammaThickness;
  fnCompt         = right.fnCompt;
  fnPhot          = right.fnPhot;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const PvskHit& PvskHit::operator=(const PvskHit& right)
{
  fEdep           = right.fEdep;
  fTrackLength    = right.fTrackLength;
  fEltrkMaxL      = right.fEltrkMaxL;
  fEltrkMaxLz     = right.fEltrkMaxLz;
  fGammaInX       = right.fGammaInX;
  fGammaInY       = right.fGammaInY;
  fGammaInZ       = right.fGammaInZ;
  fGammaThickness = right.fGammaThickness;
  fnCompt         = right.fnCompt;
  fnPhot          = right.fnPhot;

  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool PvskHit::operator==(const PvskHit& right) const
{
  return ( this == &right ) ? true : false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PvskHit::Print()
{
  G4cout
     << "Edep: " 
     << std::setw(7) << G4BestUnit(fEdep,"Energy")
     << " track length: " 
     << std::setw(7) << G4BestUnit( fTrackLength,"Length")
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
