/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "DriftChamberSD.hh"
#include "SolenoidMagneticField.hh"
#include "Constants.hh"

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

#include <TMath.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal SolenoidMagneticField* DetectorConstruction::magnetic_field_ = 0;
G4ThreadLocal G4FieldManager* DetectorConstruction::field_manager_ = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(), 
  cdc_logical_(nullptr), chcbarrel_logical_(nullptr), ncbarrel_layer1_logical_(nullptr), tracker_layer1_logical_(nullptr), tracker_layer2_logical_(nullptr)
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
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // Construct materials ////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  
  ConstructMaterials();
  auto air = G4Material::GetMaterial("G4_AIR");
  auto vacuum = G4Material::GetMaterial("G4_Galactic");
  auto argonGas = G4Material::GetMaterial("G4_Ar");
  auto scintillator = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  auto csI = G4Material::GetMaterial("G4_CESIUM_IODIDE");
  auto lead = G4Material::GetMaterial("G4_Pb");
  auto carbon = G4Material::GetMaterial("G4_C");

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////




  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // Option t oswitch on/off checking of volumes overlaps ///////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  G4bool checkOverlaps = true;
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////




  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // geometries /////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  // experimental hall ==================================================================================
  auto world_size_x = 5.*m;
  auto world_size_y = 5.*m;
  auto world_size_z = 5.*m;
  // -----
  auto world_solid
    = new G4Box("world_solid",world_size_x/2.,world_size_y/2.,world_size_z/2.);
  // -----
  auto world_logical
    = new G4LogicalVolume(world_solid,vacuum,"world_logical");
  // -----
  auto world_physical
    = new G4PVPlacement(0,G4ThreeVector(),world_logical,"world_physical",0,
        false,0,checkOverlaps);
  // ====================================================================================================


  // magnetic field =====================================================================================
  auto magnetic_size_x = 4.*m;
  auto magnetic_size_y = 4.*m;
  auto magnetic_size_z = 4.*m;
  // -----
  auto magnetic_solid
    = new G4Box("magnetic_solid",magnetic_size_x/2.,magnetic_size_y/2.,magnetic_size_z/2.);
  // -----
  magnetic_logical_ 
    = new G4LogicalVolume(magnetic_solid,vacuum,"magnetic_logical");
  // -----
  new G4PVPlacement(0,G4ThreeVector(),magnetic_logical_,"magnetic_physical",
      world_logical,false,0,checkOverlaps);
  // -----
  G4UserLimits* magnetic_userlimits = new G4UserLimits(magnetic_size_x);
  magnetic_logical_->SetUserLimits(magnetic_userlimits);
  // ====================================================================================================


  // simplest CDC (only a virtual wire plane at most ouside part) =======================================
  auto cdc_size_r = 530.*mm;
  auto cdc_size_thickness = 1.*nm;
  auto cdc_size_z = 2570.*mm;
  // -----
  auto cdc_solid
    = new G4Tubs("cdc_solid",cdc_size_r-cdc_size_thickness/2.,cdc_size_r+cdc_size_thickness/2.,cdc_size_z/2.,0.,2.*TMath::Pi());
  // -----
  cdc_logical_
    = new G4LogicalVolume(cdc_solid,vacuum,"cdc_logical");
  // -----
  auto cdc_physical
    = new G4PVPlacement(0,G4ThreeVector(),cdc_logical_,"cdc_logical",
        magnetic_logical_,false,0,checkOverlaps);
  // ====================================================================================================


  // CHCbarrel ==========================================================================================
  auto chcbarrel_size_r = 535.*mm;
  auto chcbarrel_size_thickness = 3.*mm;
  auto chcbarrel_size_z = 2570.*mm;
  // -----
  auto chcbarrel_solid
    = new G4Tubs("chcbarrel_solid",chcbarrel_size_r-(chcbarrel_size_thickness-kSpace)/2.,chcbarrel_size_r+(chcbarrel_size_thickness-kSpace)/2.,chcbarrel_size_z/2.,0.,2.*TMath::Pi());
  // -----
  chcbarrel_logical_
    = new G4LogicalVolume(chcbarrel_solid,scintillator,"chcbarrel_logical");
  // -----
  auto chcbarrel_physical
    = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),chcbarrel_logical_,"chcbarrel_logical",
        magnetic_logical_,false,0,checkOverlaps);
  // ====================================================================================================


  // NCbarrel_layer1 ====================================================================================
  auto ncbarrel_layer1_size_thickness = 30.*mm;
  auto ncbarrel_layer1_size_r = chcbarrel_size_r + chcbarrel_size_thickness/2. + ncbarrel_layer1_size_thickness/2.;
  auto ncbarrel_layer1_size_z = 2570.*mm;
  // -----
  auto ncbarrel_layer1_solid
    = new G4Tubs("ncbarrel_layer1_solid",ncbarrel_layer1_size_r-(ncbarrel_layer1_size_thickness-kSpace)/2.,ncbarrel_layer1_size_r+(ncbarrel_layer1_size_thickness-kSpace)/2.,ncbarrel_layer1_size_z/2.,0.,2.*TMath::Pi());
  // -----
  ncbarrel_layer1_logical_
    = new G4LogicalVolume(ncbarrel_layer1_solid,scintillator,"ncbarrel_layer1_logical");
  // -----
  auto ncbarrel_layer1_physical
    = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),ncbarrel_layer1_logical_,"ncbarrel_layer1_logical",
        magnetic_logical_,false,0,checkOverlaps);
  // ====================================================================================================


  // Tracker_layer1 =====================================================================================
  auto tracker_layer1_size_thickness = 5.*mm;
  auto tracker_layer1_size_r = ncbarrel_layer1_size_r + ncbarrel_layer1_size_thickness/2. + tracker_layer1_size_thickness/2.;
  auto tracker_layer1_size_z = 2570.*mm;
  // -----
  auto tracker_layer1_solid
    = new G4Tubs("tracker_layer1_solid",tracker_layer1_size_r-(tracker_layer1_size_thickness-kSpace)/2.,tracker_layer1_size_r+(tracker_layer1_size_thickness-kSpace)/2.,tracker_layer1_size_z/2.,0.,2.*TMath::Pi());
  // -----
  tracker_layer1_logical_
    = new G4LogicalVolume(tracker_layer1_solid,scintillator,"tracker_layer1_logical");
  // -----
  auto tracker_layer1_physical
    = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),tracker_layer1_logical_,"tracker_layer1_logical",
        magnetic_logical_,false,0,checkOverlaps);
  // ====================================================================================================


  // NCbarrel_layer2 ====================================================================================
  auto ncbarrel_layer2_size_thickness = 30.*mm;
  auto ncbarrel_layer2_size_r = tracker_layer1_size_r + tracker_layer1_size_thickness/2. + ncbarrel_layer2_size_thickness/2.;
  auto ncbarrel_layer2_size_z = 2570.*mm;
  // -----
  auto ncbarrel_layer2_solid
    = new G4Tubs("ncbarrel_layer2_solid",ncbarrel_layer2_size_r-(ncbarrel_layer2_size_thickness-kSpace)/2.,ncbarrel_layer2_size_r+(ncbarrel_layer2_size_thickness-kSpace)/2.,ncbarrel_layer2_size_z/2.,0.,2.*TMath::Pi());
  // -----
  ncbarrel_layer2_logical_
    = new G4LogicalVolume(ncbarrel_layer2_solid,scintillator,"ncbarrel_layer2_logical");
  // -----
  auto ncbarrel_layer2_physical
    = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),ncbarrel_layer2_logical_,"ncbarrel_layer2_logical",
        magnetic_logical_,false,0,checkOverlaps);
  // ====================================================================================================


  // Tracker_layer2 =====================================================================================
  auto tracker_layer2_size_thickness = 5.*mm;
  auto tracker_layer2_size_r = ncbarrel_layer2_size_r + ncbarrel_layer2_size_thickness/2. + tracker_layer2_size_thickness/2.;
  auto tracker_layer2_size_z = 2570.*mm;
  // -----
  auto tracker_layer2_solid
    = new G4Tubs("tracker_layer2_solid",tracker_layer2_size_r-(tracker_layer2_size_thickness-kSpace)/2.,tracker_layer2_size_r+(tracker_layer2_size_thickness-kSpace)/2.,tracker_layer2_size_z/2.,0.,2.*TMath::Pi());
  // -----
  tracker_layer2_logical_
    = new G4LogicalVolume(tracker_layer2_solid,scintillator,"tracker_layer2_logical");
  // -----
  auto tracker_layer2_physical
    = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),tracker_layer2_logical_,"tracker_layer2_logical",
        magnetic_logical_,false,0,checkOverlaps);
  // ====================================================================================================

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////





  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // visualization attributes ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  G4VisAttributes*  visAttributes;
  // -----
  visAttributes = new G4VisAttributes(Colors::Transparent());
  world_logical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  // -----
  visAttributes = new G4VisAttributes(Colors::Transparent());
  magnetic_logical_->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  // -----
  visAttributes = new G4VisAttributes(Colors::WirePlane());
  cdc_logical_->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  // -----
  visAttributes = new G4VisAttributes(Colors::Scintillator());
  chcbarrel_logical_->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  // -----
  visAttributes = new G4VisAttributes(Colors::Scintillator());
  ncbarrel_layer1_logical_->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  // -----
  visAttributes = new G4VisAttributes(Colors::Scintillator());
  ncbarrel_layer2_logical_->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  // -----
  visAttributes = new G4VisAttributes(Colors::WirePlane());
  tracker_layer1_logical_->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  // -----
  visAttributes = new G4VisAttributes(Colors::WirePlane());
  tracker_layer2_logical_->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////




  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // return the world physical volume ///////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  return world_physical;

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  auto sdManager = G4SDManager::GetSDMpointer();
  G4String SDname;

  // sensitive detectors -----------------------------------------------------
  auto cdc = new DriftChamberSD(SDname="/cdc");
  sdManager->AddNewDetector(cdc);
  cdc_logical_->SetSensitiveDetector(cdc);
  // -----
  auto tracker_layer1 = new DriftChamberSD(SDname="/tracker_layer1");
  sdManager->AddNewDetector(tracker_layer1);
  tracker_layer1_logical_->SetSensitiveDetector(tracker_layer1);
  // -----
  auto tracker_layer2 = new DriftChamberSD(SDname="/tracker_layer2");
  sdManager->AddNewDetector(tracker_layer2);
  tracker_layer2_logical_->SetSensitiveDetector(tracker_layer2);

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
