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
/// \file PvskHit.hh
/// \brief Definition of the PvskHit class

#ifndef PvskHit_h
#define PvskHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"

/// Calorimeter hit class
///
/// It defines data members to store the the energy deposit and track lengths
/// of charged particles in a selected volume:
/// - fEdep, fTrackLength

class PvskHit : public G4VHit
{
  public:
    PvskHit();
    PvskHit(const PvskHit&);
    virtual ~PvskHit();

    // operators
    const PvskHit& operator=(const PvskHit&);
    G4bool operator==(const PvskHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw() {}
    virtual void Print();

    // methods to handle data
    void Add(G4double de, G4double dl);
    void SetInitialGamma(G4double x, G4double y, G4double z);
    void AddGammaStep(G4double x, G4double y, G4double z, G4String proc);
    void AddEltrkStep(G4double l, G4double lz);

    // get methods
    G4double GetEdep() const;
    G4double GetTrackLength() const;
    G4double GetEltrkMaxL() const;
    G4double GetEltrkMaxLz() const;
    G4double GetGammaThickness() const; 
    G4double GetGammaInX() const;
    G4double GetGammaInY() const;
    G4double GetGammaInZ() const;
    G4double GetGammaDiffX() const;
    G4double GetGammaDiffY() const;
    G4int GetNCompt() const;
    G4int GetNPhot() const;
      
  private:
    G4double fEdep;        ///< Energy deposit in the sensitive volume
    G4double fTrackLength; ///< Track length in the  sensitive volume
    G4double fEltrkMaxL, fEltrkMaxLz; ///< Track length in the  sensitive volume
    G4double fGammaThickness; 
    G4double fGammaInX, fGammaInY, fGammaInZ; 
    G4double fGammaDiffX, fGammaDiffY; 
    G4int fnCompt, fnPhot;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using PvskHitsCollection = G4THitsCollection<PvskHit>;

extern G4ThreadLocal G4Allocator<PvskHit>* PvskHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* PvskHit::operator new(size_t)
{
  if (!PvskHitAllocator) {
    PvskHitAllocator = new G4Allocator<PvskHit>;
  }
  void *hit;
  hit = (void *) PvskHitAllocator->MallocSingle();
  return hit;
}

inline void PvskHit::operator delete(void *hit)
{
  if (!PvskHitAllocator) {
    PvskHitAllocator = new G4Allocator<PvskHit>;
  }
  PvskHitAllocator->FreeSingle((PvskHit*) hit);
}

inline void PvskHit::Add(G4double de, G4double dl) {
  fEdep += de; 
  fTrackLength += dl;
}

inline void PvskHit::SetInitialGamma(G4double x, G4double y, G4double z) {
  fGammaInX = x; 
  fGammaInY = y;
  fGammaInZ = z;
  fGammaDiffX = 0; 
  fGammaDiffY = 0;
  fGammaThickness = 0;
}

inline void PvskHit::AddEltrkStep(G4double l, G4double lz) {
  if( fEltrkMaxL < l ) fEltrkMaxL = l;
  if( fEltrkMaxLz < lz ) fEltrkMaxLz = lz;
}

inline void PvskHit::AddGammaStep(G4double x, G4double y, G4double z, G4String proc) {
  fGammaDiffX += x; 
  fGammaDiffY += y;
  fGammaThickness += z;
  if( proc=="compt" ) fnCompt++;
  else if( proc=="phot" ) fnPhot++;
  //G4cout << "dx: " << x << ", dy: " << y << ", dz(thickness): " << z << G4endl
  //<< "total dx: " << fGammaDiffX << ", dy: " << fGammaDiffY << ", dz(thickness): " << fGammaThickness << G4endl
  //<< "# compt: " << fnCompt << ", # phot: " << fnPhot << G4endl;
}

inline G4double PvskHit::GetEdep() const { 
  return fEdep; 
}

inline G4double PvskHit::GetTrackLength() const { 
  return fTrackLength; 
}

inline G4double PvskHit::GetEltrkMaxL() const { 
  return fEltrkMaxL; 
}

inline G4double PvskHit::GetEltrkMaxLz() const { 
  return fEltrkMaxLz; 
}

inline G4double PvskHit::GetGammaInX() const { 
  return fGammaInX; 
}

inline G4double PvskHit::GetGammaInY() const { 
  return fGammaInY; 
}

inline G4double PvskHit::GetGammaInZ() const { 
  return fGammaInZ; 
}

inline G4double PvskHit::GetGammaDiffX() const { 
  return fGammaDiffX; 
}

inline G4double PvskHit::GetGammaDiffY() const { 
  return fGammaDiffY; 
}

inline G4double PvskHit::GetGammaThickness() const { 
  return fGammaThickness; 
}

inline G4int PvskHit::GetNCompt() const { 
  return fnCompt; 
}

inline G4int PvskHit::GetNPhot() const { 
  return fnPhot; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
