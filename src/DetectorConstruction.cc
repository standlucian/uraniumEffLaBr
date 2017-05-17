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
// $Id: B1DetectorConstruction.cc 75117 2013-10-28 09:38:37Z gcosmo $
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "G4Material.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSTrackLength.hh"
#include "G4PhysicalConstants.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{  
  std::vector<G4int> natoms;
  std::vector<G4String> elements;

  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Envelope parameters
  //

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;
  //     
  // World
  //
  G4double world_sizeXY = 33.5*cm;
  G4double world_sizeZ  = 55*cm;
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
  // Shape 1
  //  
  G4Material* pb_mat     = nist->FindOrBuildMaterial("G4_Pb");

  elements.push_back("C");     natoms.push_back(5);
  elements.push_back("H");     natoms.push_back(8);
  elements.push_back("O");     natoms.push_back(2);

  G4Material* PMMA_mat = nist->ConstructNewMaterial("PMMA", elements, natoms, 1.190*g/cm3);

  elements.clear();
  natoms.clear();

  elements.push_back("C");     natoms.push_back(12);
  elements.push_back("H");     natoms.push_back(5);
  elements.push_back("Cl");     natoms.push_back(5);

  G4Material* pcb_mat = nist->ConstructNewMaterial("pcb", elements, natoms, 1.374*g/cm3);

  elements.clear();
  natoms.clear();

  elements.push_back("Cd");     natoms.push_back(1);
  elements.push_back("Zn");     natoms.push_back(1);
  elements.push_back("Te");     natoms.push_back(1);

  G4Material* czt_mat = nist->ConstructNewMaterial("czt", elements, natoms, 5.56*g/cm3);

  elements.clear();
  natoms.clear();

  elements.push_back("U");     natoms.push_back(3);
  elements.push_back("O");     natoms.push_back(8);

  G4Material* uo_mat = nist->ConstructNewMaterial("uo", elements, natoms, 3.29*g/cm3);

  G4Material* pvc_mat     = nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");

  G4Material* al_mat     = nist->FindOrBuildMaterial("G4_Al");

   // Outside Shielding

 G4VSolid* ShieldingCyl1 = new G4Tubs ( "ShieldingCyl1", 0*mm, 335/2.*mm, 505/2.*mm, 0.*deg, 360.*deg );

 G4VSolid* ShieldingCyl2 = new G4Tubs ( "ShieldingCyl2", 0*mm, 205/2.*mm, 396/2.*mm, 0.*deg, 360.*deg );

 G4VSolid* ShieldingCyl3 = new G4Tubs ( "ShieldingCyl3", 0*mm, 85/2.*mm, 53/2.*mm, 0.*deg, 360.*deg );

 G4VSolid* Shield_subst1 = new G4SubtractionSolid ("Shield_subst1", ShieldingCyl1, ShieldingCyl2, 0, G4ThreeVector(0.,0.,-2.5*mm));

 G4VSolid* Shield_subst2 = new G4SubtractionSolid ("Shield_subst2", Shield_subst1, ShieldingCyl3, 0, G4ThreeVector(0.,0.,-226.5*mm));

 G4LogicalVolume* logicShield  = new G4LogicalVolume (Shield_subst2, pb_mat, "logicShield");

 new G4PVPlacement(0, G4ThreeVector(0.,0.,-22.5*mm), logicShield, "logicShield", logicWorld, false, 0, checkOverlaps);

  // Sample holder

 G4VSolid* SampleHolder = new G4Tubs ( "SampleHolder", 0*mm, 85/2.*mm, 5/2.*mm, 0.*deg, 360.*deg );

 G4LogicalVolume* logicHolder  = new G4LogicalVolume (SampleHolder, PMMA_mat, "logicHolder");

 new G4PVPlacement(0, G4ThreeVector(0.,0.,-5.5*mm), logicHolder, "logicHolder", logicWorld, false, 0, checkOverlaps);

 // Sample spacer

 G4VSolid* SampleSpacer1 = new G4Tubs ( "SampleSpacer1", 41.6/2.*mm, 58.2/2.*mm, 40/2.*mm, 0.*deg, 360.*deg );

 G4VSolid* SampleSpacer2 = new G4Tubs ( "SampleSpacer2", 0*mm, 49.9/2.*mm, 28.15/2.*mm, 0.*deg, 360.*deg );

 G4VSolid* Spacer_subst = new G4SubtractionSolid ("Spacer_subst", SampleSpacer1, SampleSpacer2, 0, G4ThreeVector(0.,0.,5.925*mm));

 G4LogicalVolume* logicSpacer  = new G4LogicalVolume (Spacer_subst, PMMA_mat, "logicSpacer");

 new G4PVPlacement(0, G4ThreeVector(0.,0.,-28*mm), logicSpacer, "logicSpacer", logicWorld, false, 0, checkOverlaps);

 // Detector Holder

 G4VSolid* DetectorHolder1 = new G4Tubs ( "DetectorHolder1", 24/2.*mm, 62/2.*mm, 147/2.*mm, 0.*deg, 360.*deg );

 G4LogicalVolume* logicDetHolder  = new G4LogicalVolume (DetectorHolder1, pvc_mat, "logicDetHolder");

 new G4PVPlacement(0, G4ThreeVector(0.,0.,-201.5*mm), logicDetHolder, "logicDetHolder", logicWorld, false, 0, checkOverlaps);

 G4VSolid* DetectorHolder2 = new G4Tubs ( "DetectorHolder2", 24/2.*mm, 41.4/2.*mm, 29.4/2.*mm, 0.*deg, 360.*deg );

 G4LogicalVolume* logicDetHolder2  = new G4LogicalVolume (DetectorHolder2, pvc_mat, "logicDetHolder2");

  new G4PVPlacement(0, G4ThreeVector(0.,0.,-62.7*mm), logicDetHolder2, "logicDetHolder2", logicWorld, false, 0, checkOverlaps);

 G4VSolid* DetectorHolder3 = new G4Cons ( "DetectorHolder3", 0*mm, 62/2.*mm, 0*mm, 45.5/2.*mm, 50.6/2.*mm,0.*deg, 360.*deg);

 G4VSolid* DetectorHolder4 = new G4Tubs ( "DetectorHolder4", 0*mm, 24/2.*mm, 51.5/2.*mm, 0.*deg, 360.*deg );

 G4VSolid* DetHolder_subst = new G4SubtractionSolid ("DetHolder_subst", DetectorHolder3, DetectorHolder4, 0, G4ThreeVector(0.,0.,0.*mm));

 G4LogicalVolume* logicDetHolder3  = new G4LogicalVolume (DetHolder_subst, pvc_mat, "logicDetHolder3");

 new G4PVPlacement(0, G4ThreeVector(0.,0.,-102.7*mm), logicDetHolder3, "logicDetHolder3", logicWorld, false, 0, checkOverlaps);

 // Detector casing

 G4VSolid* DetectorCasing1 = new G4Tubs ( "DetectorCasing1", 21/2.*mm, 24/2.*mm, 121.1/2.*mm, 0.*deg, 360.*deg );

 G4LogicalVolume* logicDetCase1  = new G4LogicalVolume (DetectorCasing1, al_mat, "logicDetCase1");

 new G4PVPlacement(0, G4ThreeVector(0.,0.,-84.25*mm), logicDetCase1, "logicDetCase1", logicWorld, false, 0, checkOverlaps);

 G4VSolid* DetectorCasing2 = new G4Tubs ( "DetectorCasing2", 0*mm, 24/2.*mm, 1.2/2.*mm, 0.*deg, 360.*deg );

 G4LogicalVolume* logicDetCase2  = new G4LogicalVolume (DetectorCasing2, al_mat, "logicDetCase2");

 new G4PVPlacement(0, G4ThreeVector(0.,0.,-23.1*mm), logicDetCase2, "logicDetCase2", logicWorld, false, 0, checkOverlaps);

 G4VSolid* DetectorCasing3 = new G4Tubs ( "DetectorCasing3", 0*mm, 21/2.*mm, 4.1/2.*mm, 0.*deg, 360.*deg );

 G4LogicalVolume* logicDetCase3  = new G4LogicalVolume (DetectorCasing3, al_mat, "logicDetCase3");

 new G4PVPlacement(0, G4ThreeVector(0.,0.,-35.45*mm), logicDetCase3, "logicDetCase3", logicWorld, false, 0, checkOverlaps);

 G4VSolid* DetectorCasing4 = new G4Tubs ( "DetectorCasing4", 0*mm, 2.3/2.*mm, 106.2/2.*mm, 0.*deg, 360.*deg );

 G4LogicalVolume* logicDetCase4 = new G4LogicalVolume (DetectorCasing4, al_mat, "logicDetCase4");

 new G4PVPlacement(0, G4ThreeVector(0.,0.,-90.6*mm), logicDetCase4, "logicDetCase4", logicWorld, false, 0, checkOverlaps);

 G4VSolid* DetectorCasing5 = new G4Tubs ( "DetectorCasing5", 0*mm, 21/2.*mm, 0.8/2.*mm, 0.*deg, 360.*deg );

 G4LogicalVolume* logicDetCase5 = new G4LogicalVolume (DetectorCasing5, pcb_mat, "logicDetCase5");

 new G4PVPlacement(0, G4ThreeVector(0.,0.,-33*mm), logicDetCase5, "logicDetCase5", logicWorld, false, 0, checkOverlaps);

 //Detector volume

 G4VSolid* Detector = new G4Box ( "Detector", 10/2.*mm, 10/2.*mm, 5/2.*mm );


 //Active detector volume

 G4VSolid* DetectorActive = new G4Box ( "DetectorActive", 8.875/2.*mm, 8.875/2.*mm, 3.25/2.*mm );

 G4LogicalVolume* logicDetector = new G4LogicalVolume (DetectorActive, czt_mat, "logicDetector");

 new G4PVPlacement(0, G4ThreeVector(0.,0.,-29.35*mm), logicDetector, "logicDetector", logicWorld, false, 0, checkOverlaps);

 G4VSolid* DeadLayer = new G4SubtractionSolid ("DeadLayer", Detector, DetectorActive, 0, G4ThreeVector(0.,0.,0.75*mm));

 G4LogicalVolume* logicDead = new G4LogicalVolume (DeadLayer, czt_mat, "logicDead");

 new G4PVPlacement(0, G4ThreeVector(0.,0.,-30.1*mm), logicDead, "logicDead", logicWorld, false, 0, checkOverlaps);

 G4VSolid* SourceCasing1 = new G4Tubs ( "SourceCasing1", 0.*mm, 80/2.*mm, 89/2.*mm, 0.*deg, 360.*deg );

 G4VSolid* Source = new G4Tubs ( "Source", 0.*mm, 70/2.*mm, 20.8/2.*mm, 0.*deg, 360.*deg );

 G4VSolid* SourceCasing2 = new G4SubtractionSolid ("SourceCasing2", SourceCasing1, Source, 0, G4ThreeVector(0.,0.,-32.1*mm));

 G4LogicalVolume* logicSourceCasing  = new G4LogicalVolume (SourceCasing2, al_mat, "SourceCasing2");

 new G4PVPlacement(0, G4ThreeVector(0.,0.,43.5*mm), logicSourceCasing, "logicSourceCasing", logicWorld, false, 0, checkOverlaps);

 G4LogicalVolume* logicSource  = new G4LogicalVolume (Source, uo_mat, "Source");

 new G4PVPlacement(0, G4ThreeVector(0.,0.,11.4*mm), logicSource, "logicSource", logicWorld, false, 0, checkOverlaps);

 fScoringVolume = logicDetector;
                                 
  return physWorld;
}

void DetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
  // 
  // Scorers
  //

  // declare Absorber as a MultiFunctionalDetector scorer
  //  
  G4MultiFunctionalDetector* absDetector
    = new G4MultiFunctionalDetector("Absorber");

  G4VPrimitiveScorer* primitive;
  primitive = new G4PSEnergyDeposit("Edep");
  absDetector->RegisterPrimitive(primitive);


  SetSensitiveDetector("logicDetector",absDetector, true);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
