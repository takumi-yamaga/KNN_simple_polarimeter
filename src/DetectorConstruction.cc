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
/// \copied from B5DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "DriftChamberSD.hh"
#include "SolenoidMagneticField.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"
#include "G4UserLimits.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal SolenoidMagneticField* DetectorConstruction::magnetic_field_ = 0;
G4ThreadLocal G4FieldManager* DetectorConstruction::field_manager_ = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(), 
  dcin_wireplane_logical_(nullptr), dcout_wireplane_logical_(nullptr)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::~DetectorConstruction()
{
  //delete fMessenger;

  for (auto visAttributes: fVisAttributes) {
    delete visAttributes;
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Construct materials
  ConstructMaterials();
  auto air = G4Material::GetMaterial("G4_AIR");
  auto vacuum = G4Material::GetMaterial("G4_Galactic");
  auto argonGas = G4Material::GetMaterial("G4_Ar");
  auto scintillator = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  auto csI = G4Material::GetMaterial("G4_CESIUM_IODIDE");
  auto lead = G4Material::GetMaterial("G4_Pb");
  auto carbon = G4Material::GetMaterial("G4_C");
  auto hydrogen = new G4Material("hydrogen", 1., 1.01*g/mole, 1.*g/cm3);

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  // geometries --------------------------------------------------------------

  // experimental hall (world volume)
  auto worldSolid 
    = new G4Box("worldBox",10.*m,3.*m,10.*m);
  auto worldLogical
    = new G4LogicalVolume(worldSolid,vacuum,"worldLogical");
  auto worldPhysical
    = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0,
        false,0,checkOverlaps);

  // magnetic field
  auto magnetic_size_x = 1.*m;
  auto magnetic_size_y = 1.*m;
  auto magnetic_size_z = 3.*m;
  auto magnetic_solid
    = new G4Box("dcin_Box",magnetic_size_x/2.,magnetic_size_y/2.,magnetic_size_z/2.);
  magnetic_logical_ = new G4LogicalVolume(magnetic_solid,air,"magnetic_logical");
  new G4PVPlacement(0,G4ThreeVector(),magnetic_logical_,"magnetic_physical",
      worldLogical,false,0,checkOverlaps);
  // set step limit in magnetic field  
  G4UserLimits* magnetic_userlimits = new G4UserLimits(magnetic_size_x);
  magnetic_logical_->SetUserLimits(magnetic_userlimits);

  // target 
  auto target_size_x = 50.*cm;
  auto target_size_y = 50.*cm;
  auto target_thickness = 30.*mm; 
  auto target_position_z = target_thickness/2. + 5.*mm;
  auto targetSolid 
    = new G4Box("targetBox",target_size_x/2.,target_size_y/2.,target_thickness/2.);
  auto targetLogical
    = new G4LogicalVolume(targetSolid,vacuum,"targetLogical");
  auto targetPhysical
    = new G4PVPlacement(0,G4ThreeVector(0.,0.,target_position_z),targetLogical,"targetPhysical",
        magnetic_logical_,false,0,checkOverlaps);

  // drift chamber (in)
  auto dc_size_x = target_size_x;
  auto dc_size_y = target_size_y;
  auto dc_thickness = 1.*mm;
  auto dcin_position = G4ThreeVector(0.,0.,target_position_z-(target_thickness/2.+dc_thickness/2.+kSpace));
  auto dcin_Solid 
    = new G4Box("dcin_Box",dc_size_x/2.,dc_size_y/2.,dc_thickness/2.);
  auto dcin_Logical
    = new G4LogicalVolume(dcin_Solid,air,"dcin_Logical");
  auto dcin_Physical
    = new G4PVPlacement(0,dcin_position,dcin_Logical,"dcin_Physical",
        magnetic_logical_,false,0,checkOverlaps);
  // wireplane
  auto dc_wireplane_thickness = 1.*nm;
  auto dcin_wireplane_solid
    = new G4Box("dcin_wireplane_box", dc_size_x/2., dc_size_y/2., dc_wireplane_thickness/2.);
  dcin_wireplane_logical_
    = new G4LogicalVolume(dcin_wireplane_solid,argonGas,"dcin_wireplane_logical");
  auto dcin_wireplane_physical
    = new G4PVPlacement(0,G4ThreeVector(),dcin_wireplane_logical_,"dcin_wireplane_physical",
        dcin_Logical, false,0,checkOverlaps);

  // drift chamber (out)
  auto dcout_position = G4ThreeVector(0.,0.,target_position_z+(target_thickness/2.+dc_thickness/2.+kSpace));
  auto dcout_Solid 
    = new G4Box("dcout_Box",dc_size_x/2.,dc_size_y/2.,dc_thickness/2.);
  auto dcout_Logical
    = new G4LogicalVolume(dcout_Solid,air,"dcout_Logical");
  auto dcout_Physical
    = new G4PVPlacement(0,dcout_position,dcout_Logical,"dcout_Physical", 
        magnetic_logical_, false,0,checkOverlaps);
  // wireplane
  auto dcout_wireplane_solid
    = new G4Box("dcout_wireplane_box", dc_size_x/2., dc_size_y/2., dc_wireplane_thickness/2.);
  dcout_wireplane_logical_
    = new G4LogicalVolume(dcout_wireplane_solid,argonGas,"dcout_wireplane_logical");
  auto dcout_wireplane_physical
    = new G4PVPlacement(0,G4ThreeVector(),dcout_wireplane_logical_,"dcout_wireplane_physical",
        dcout_Logical, false,0,checkOverlaps);


  // visualization attributes ------------------------------------------------

  auto visAttributes = new G4VisAttributes(G4Colour::White());
  visAttributes->SetVisibility(false);
  worldLogical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  magnetic_logical_->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  visAttributes = new G4VisAttributes(G4Colour::Blue());
  targetLogical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  visAttributes = new G4VisAttributes(G4Colour::Yellow());
  dcin_Logical->SetVisAttributes(visAttributes);
  dcout_Logical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  visAttributes = new G4VisAttributes(G4Colour::Green());
  dcin_wireplane_logical_->SetVisAttributes(visAttributes);
  dcout_wireplane_logical_->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);


  // return the world physical volume ----------------------------------------

  return worldPhysical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  auto sdManager = G4SDManager::GetSDMpointer();
  G4String SDname;

  // sensitive detectors -----------------------------------------------------
  auto dcin = new DriftChamberSD(SDname="/dcin");
  sdManager->AddNewDetector(dcin);
  dcin_wireplane_logical_->SetSensitiveDetector(dcin);

  auto dcout = new DriftChamberSD(SDname="/dcout");
  sdManager->AddNewDetector(dcout);
  dcout_wireplane_logical_->SetSensitiveDetector(dcout);

  // magnetic field ----------------------------------------------------------
  magnetic_field_ = new SolenoidMagneticField();
  field_manager_ = new G4FieldManager();
  field_manager_->SetDetectorField(magnetic_field_);
  field_manager_->CreateChordFinder(magnetic_field_);
  G4bool force_to_all_daughters = true; // if true, all daughters have the same field
  magnetic_logical_->SetFieldManager(field_manager_, force_to_all_daughters);
  // -------------------------------------------------------------------------
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructMaterials()
{
  auto nistManager = G4NistManager::Instance();

  // Air 
  nistManager->FindOrBuildMaterial("G4_AIR");

  // Argon gas
  nistManager->FindOrBuildMaterial("G4_Ar");

  // Scintillator
  // (PolyVinylToluene, C_9H_10)
  nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

  // CsI
  nistManager->FindOrBuildMaterial("G4_CESIUM_IODIDE");

  // Lead
  nistManager->FindOrBuildMaterial("G4_Pb");

  // Carbon
  nistManager->FindOrBuildMaterial("G4_C");

  // Vacuum "Galactic"
  nistManager->FindOrBuildMaterial("G4_Galactic");


  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
