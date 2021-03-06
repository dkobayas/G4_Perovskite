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
/// \file PvskDetectorConstruction.cc
/// \brief Implementation of the PvskDetectorConstruction class

#include "PvskDetectorConstruction.hh"

#include "PvskSD.hh"
#include "G4SDManager.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PvskDetectorConstruction::PvskDetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PvskDetectorConstruction::~PvskDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* PvskDetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Envelope parameters
  //
  G4double env_sizeXY = 30*cm, env_sizeZ = 30*cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  //     
  // Envelope
  //  
  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
 
  //==============================================================
  // MAPbI3 shared parameters
  //==============================================================
  //     
  // MAPBI3 = CH3NH3PbI3(C1H6N1Pb1I3) 
  // only 0.62[kg/mol] for information... So, calculated using molar volume of BaTiO3( perovskite construction )
  // Then, 620[g/mol]/38[cm3/mol] = 16.3[g/cm3]
  //  
  G4double a = 12.01*g/mole;
  G4Element* elC  = new G4Element("Carbon", "C" , 6., a);
  a = 1.01*g/mole;
  G4Element* elH  = new G4Element("Hydrogen", "H" , 1., a);
  a = 14.01*g/mole;
  G4Element* elN  = new G4Element("Nitrogen", "N" , 7., a);
  a = 207.2*g/mole;
  G4Element* elPb = new G4Element("Lead"    ,"Pb", 82., a);
  a = 132.9*g/mole;
  G4Element* elCs = new G4Element("Cs"    ,"Cs", 55., a);
  a = 126.9*g/mole;
  G4Element* elI  = new G4Element("I", "I" , 53., a);
  a = 79.9*g/mole;
  G4Element* elBr  = new G4Element("Br", "Br" , 35., a);
  //G4double density = 16.3*g/cm3;
  G4double density = 4*g/cm3;
  G4Material* MAPbI3 = new G4Material("MAPbI3", density, 5);
  G4Material* CsPbBr3 = new G4Material("CsPbBr3", density, 3);
  //
  // 1 single cristal
  //
  MAPbI3->AddElement(elC , 1);
  MAPbI3->AddElement(elH , 6);
  MAPbI3->AddElement(elN , 1);
  MAPbI3->AddElement(elPb , 1);
  MAPbI3->AddElement(elI , 3);
  G4double MAPbI3_dxa = 6*mm, MAPbI3_dxb = 6*mm;
  G4double MAPbI3_dya = 6*mm, MAPbI3_dyb = 6*mm;
  G4double MAPbI3_dz  = 2*mm;      
  G4Trd* solidShapeMAPbI3 =    
    new G4Trd("ShapeMAPbI3",                      //its name
              0.5*MAPbI3_dxa, 0.5*MAPbI3_dxb, 
              0.5*MAPbI3_dya, 0.5*MAPbI3_dyb, 0.5*MAPbI3_dz); //its size
  CsPbBr3->AddElement(elCs, 1);
  CsPbBr3->AddElement(elPb, 1);
  CsPbBr3->AddElement(elBr, 3);
  G4double CsPbBr3_dxa = 6*mm, CsPbBr3_dxb = 6*mm;
  G4double CsPbBr3_dya = 6*mm, CsPbBr3_dyb = 6*mm;
  G4double CsPbBr3_dz  = 2*mm;      
  G4Trd* solidShapeCsPbBr3 =    
    new G4Trd("ShapeCsPbBr3",                      //its name
              0.5*CsPbBr3_dxa, 0.5*CsPbBr3_dxb, 
              0.5*CsPbBr3_dya, 0.5*CsPbBr3_dyb, 0.5*CsPbBr3_dz); //its size
  
  // Detector #1
  G4ThreeVector pos1 = G4ThreeVector(0*cm, 0*cm, 3*cm);
  G4LogicalVolume* logicShapeMAPbI3_1 =                         
    new G4LogicalVolume(solidShapeMAPbI3,         //its solid
                        MAPbI3,          //its material
                        "ShapeMAPbI3_1");           //its name
  new G4PVPlacement(0,                       //no rotation
                    pos1,                    //at position
                    logicShapeMAPbI3_1,             //its logical volume
                    "ShapeMAPbI3_1",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  
  fScoringVolume = logicShapeMAPbI3_1;
  
  // Detector #2
  G4ThreeVector pos2 = G4ThreeVector(0*cm, 0*cm, -3*cm);
  G4LogicalVolume* logicShapeCsPbBr3_2 =                         
    new G4LogicalVolume(solidShapeCsPbBr3,         //its solid
                        CsPbBr3,          //its material
                        "ShapeCsPbBr3_2");           //its name
  new G4PVPlacement(0,                       //no rotation
                    pos2,                    //at position
                    logicShapeCsPbBr3_2,             //its logical volume
                    "ShapeCsPbBr3_2",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  
  //fScoringVolume = logicShapeMAPbI3_1;

  //
  //always return the physical World
  //
  return physWorld;
}

void PvskDetectorConstruction::ConstructSDandField()
{
  // G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // 
  // Sensitive detectors
  //
  auto perovskiteSD_1
    = new PvskSD("PvskSD_1", "PvskHitsCollection_1", 1);
  G4SDManager::GetSDMpointer()->AddNewDetector(perovskiteSD_1);
  SetSensitiveDetector("ShapeMAPbI3_1",perovskiteSD_1);
  auto perovskiteSD_2
    = new PvskSD("PvskSD_2", "PvskHitsCollection_2", 1);
  G4SDManager::GetSDMpointer()->AddNewDetector(perovskiteSD_2);
  SetSensitiveDetector("ShapeCsPbBr3_2",perovskiteSD_2);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
