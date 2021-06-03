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
  cdc_logical_(nullptr),
  chcbarrel_logical_(nullptr),
  ncbarrel_layer1_logical_(nullptr), ncbarrel_layer2_logical_(nullptr), ncbarrel_layer3_logical_(nullptr), ncbarrel_layer4_logical_(nullptr), ncbarrel_layer5_logical_(nullptr),
  tracker_layer1_logical_(nullptr), tracker_layer2_logical_(nullptr),
  tracker_layer3_logical_(nullptr), tracker_layer4_logical_(nullptr), tracker_layer5_logical_(nullptr), tracker_layer6_logical_(nullptr),
  tracker_layer7_logical_(nullptr), tracker_layer8_logical_(nullptr), tracker_layer9_logical_(nullptr), tracker_layer10_logical_(nullptr)
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


  // ====================================================================================================
  auto ncbarrel_thickness = 50.*mm;
  auto carbon_thickness   = 50.*mm;
  auto tracker_thickness  =  5.*mm;
  auto space_thickness    = 50.*mm;
  // ====================================================================================================


  // NCbarrel_layer1 ====================================================================================
  auto ncbarrel_layer1_size_thickness = ncbarrel_thickness;
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
  auto tracker_layer1_size_thickness = tracker_thickness;
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


  // Space ==============================================================================================
  auto space_layer1to2_size_thickness = space_thickness;
  // ====================================================================================================


  // Tracker_layer2 =====================================================================================
  auto tracker_layer2_size_thickness = tracker_thickness;
  auto tracker_layer2_size_r = tracker_layer1_size_r + tracker_layer1_size_thickness/2. + space_layer1to2_size_thickness + tracker_layer2_size_thickness/2.;
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


  // NCbarrel_layer2 ====================================================================================
  auto ncbarrel_layer2_size_thickness = carbon_thickness;
  auto ncbarrel_layer2_size_r = tracker_layer2_size_r + tracker_layer2_size_thickness/2. + ncbarrel_layer2_size_thickness/2.;
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


  // Tracker_layer3 =====================================================================================
  auto tracker_layer3_size_thickness = tracker_thickness;
  auto tracker_layer3_size_r = ncbarrel_layer2_size_r + ncbarrel_layer2_size_thickness/2. + tracker_layer3_size_thickness/2.;
  auto tracker_layer3_size_z = 2570.*mm;
  // -----
  auto tracker_layer3_solid
    = new G4Tubs("tracker_layer3_solid",tracker_layer3_size_r-(tracker_layer3_size_thickness-kSpace)/2.,tracker_layer3_size_r+(tracker_layer3_size_thickness-kSpace)/2.,tracker_layer3_size_z/2.,0.,2.*TMath::Pi());
  // -----
  tracker_layer3_logical_
    = new G4LogicalVolume(tracker_layer3_solid,scintillator,"tracker_layer3_logical");
  // -----
  auto tracker_layer3_physical
    = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),tracker_layer3_logical_,"tracker_layer3_logical",
        magnetic_logical_,false,0,checkOverlaps);
  // ====================================================================================================


  // Space ==============================================================================================
  auto space_layer3to4_size_thickness = space_thickness;
  // ====================================================================================================


  // Tracker_layer4 =====================================================================================
  auto tracker_layer4_size_thickness = tracker_thickness;
  auto tracker_layer4_size_r = tracker_layer3_size_r + tracker_layer3_size_thickness/2. + space_layer3to4_size_thickness + tracker_layer4_size_thickness/2.;
  auto tracker_layer4_size_z = 2570.*mm;
  // -----
  auto tracker_layer4_solid
    = new G4Tubs("tracker_layer4_solid",tracker_layer4_size_r-(tracker_layer4_size_thickness-kSpace)/2.,tracker_layer4_size_r+(tracker_layer4_size_thickness-kSpace)/2.,tracker_layer4_size_z/2.,0.,2.*TMath::Pi());
  // -----
  tracker_layer4_logical_
    = new G4LogicalVolume(tracker_layer4_solid,scintillator,"tracker_layer4_logical");
  // -----
  auto tracker_layer4_physical
    = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),tracker_layer4_logical_,"tracker_layer4_logical",
        magnetic_logical_,false,0,checkOverlaps);
  // ====================================================================================================


  // NCbarrel_layer3 ====================================================================================
  auto ncbarrel_layer3_size_thickness = carbon_thickness;
  auto ncbarrel_layer3_size_r = tracker_layer4_size_r + tracker_layer4_size_thickness/2. + ncbarrel_layer3_size_thickness/2.;
  auto ncbarrel_layer3_size_z = 2570.*mm;
  // -----
  auto ncbarrel_layer3_solid
    = new G4Tubs("ncbarrel_layer3_solid",ncbarrel_layer3_size_r-(ncbarrel_layer3_size_thickness-kSpace)/2.,ncbarrel_layer3_size_r+(ncbarrel_layer3_size_thickness-kSpace)/2.,ncbarrel_layer3_size_z/2.,0.,2.*TMath::Pi());
  // -----
  ncbarrel_layer3_logical_
    = new G4LogicalVolume(ncbarrel_layer3_solid,scintillator,"ncbarrel_layer3_logical");
  // -----
  auto ncbarrel_layer3_physical
    = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),ncbarrel_layer3_logical_,"ncbarrel_layer3_logical",
        magnetic_logical_,false,0,checkOverlaps);
  // ====================================================================================================


  // Tracker_layer5 =====================================================================================
  auto tracker_layer5_size_thickness = tracker_thickness;
  auto tracker_layer5_size_r = ncbarrel_layer3_size_r + ncbarrel_layer3_size_thickness/2. + tracker_layer5_size_thickness/2.;
  auto tracker_layer5_size_z = 2570.*mm;
  // -----
  auto tracker_layer5_solid
    = new G4Tubs("tracker_layer5_solid",tracker_layer5_size_r-(tracker_layer5_size_thickness-kSpace)/2.,tracker_layer5_size_r+(tracker_layer5_size_thickness-kSpace)/2.,tracker_layer5_size_z/2.,0.,2.*TMath::Pi());
  // -----
  tracker_layer5_logical_
    = new G4LogicalVolume(tracker_layer5_solid,scintillator,"tracker_layer5_logical");
  // -----
  auto tracker_layer5_physical
    = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),tracker_layer5_logical_,"tracker_layer5_logical",
        magnetic_logical_,false,0,checkOverlaps);
  // ====================================================================================================


  // Space ==============================================================================================
  auto space_layer5to6_size_thickness = space_thickness;
  // ====================================================================================================


  // Tracker_layer6 =====================================================================================
  auto tracker_layer6_size_thickness = tracker_thickness;
  auto tracker_layer6_size_r = tracker_layer5_size_r + tracker_layer5_size_thickness/2. + space_layer5to6_size_thickness + tracker_layer6_size_thickness/2.;
  auto tracker_layer6_size_z = 2570.*mm;
  // -----
  auto tracker_layer6_solid
    = new G4Tubs("tracker_layer6_solid",tracker_layer6_size_r-(tracker_layer6_size_thickness-kSpace)/2.,tracker_layer6_size_r+(tracker_layer6_size_thickness-kSpace)/2.,tracker_layer6_size_z/2.,0.,2.*TMath::Pi());
  // -----
  tracker_layer6_logical_
    = new G4LogicalVolume(tracker_layer6_solid,scintillator,"tracker_layer6_logical");
  // -----
  auto tracker_layer6_physical
    = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),tracker_layer6_logical_,"tracker_layer6_logical",
        magnetic_logical_,false,0,checkOverlaps);
  // ====================================================================================================


  // NCbarrel_layer4 ====================================================================================
  auto ncbarrel_layer4_size_thickness = carbon_thickness;
  auto ncbarrel_layer4_size_r = tracker_layer6_size_r + tracker_layer6_size_thickness/2. + ncbarrel_layer4_size_thickness/2.;
  auto ncbarrel_layer4_size_z = 2570.*mm;
  // -----
  auto ncbarrel_layer4_solid
    = new G4Tubs("ncbarrel_layer4_solid",ncbarrel_layer4_size_r-(ncbarrel_layer4_size_thickness-kSpace)/2.,ncbarrel_layer4_size_r+(ncbarrel_layer4_size_thickness-kSpace)/2.,ncbarrel_layer4_size_z/2.,0.,2.*TMath::Pi());
  // -----
  ncbarrel_layer4_logical_
    = new G4LogicalVolume(ncbarrel_layer4_solid,scintillator,"ncbarrel_layer4_logical");
  // -----
  auto ncbarrel_layer4_physical
    = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),ncbarrel_layer4_logical_,"ncbarrel_layer4_logical",
        magnetic_logical_,false,0,checkOverlaps);
  // ====================================================================================================


  // Tracker_layer7 =====================================================================================
  auto tracker_layer7_size_thickness = tracker_thickness;
  auto tracker_layer7_size_r = ncbarrel_layer4_size_r + ncbarrel_layer4_size_thickness/2. + tracker_layer7_size_thickness/2.;
  auto tracker_layer7_size_z = 2570.*mm;
  // -----
  auto tracker_layer7_solid
    = new G4Tubs("tracker_layer7_solid",tracker_layer7_size_r-(tracker_layer7_size_thickness-kSpace)/2.,tracker_layer7_size_r+(tracker_layer7_size_thickness-kSpace)/2.,tracker_layer7_size_z/2.,0.,2.*TMath::Pi());
  // -----
  tracker_layer7_logical_
    = new G4LogicalVolume(tracker_layer7_solid,scintillator,"tracker_layer7_logical");
  // -----
  auto tracker_layer7_physical
    = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),tracker_layer7_logical_,"tracker_layer7_logical",
        magnetic_logical_,false,0,checkOverlaps);
  // ====================================================================================================


  // Space ==============================================================================================
  auto space_layer7to8_size_thickness = space_thickness;
  // ====================================================================================================


  // Tracker_layer8 =====================================================================================
  auto tracker_layer8_size_thickness = tracker_thickness;
  auto tracker_layer8_size_r = tracker_layer7_size_r + tracker_layer7_size_thickness/2. + space_layer7to8_size_thickness + tracker_layer8_size_thickness/2.;
  auto tracker_layer8_size_z = 2570.*mm;
  // -----
  auto tracker_layer8_solid
    = new G4Tubs("tracker_layer8_solid",tracker_layer8_size_r-(tracker_layer8_size_thickness-kSpace)/2.,tracker_layer8_size_r+(tracker_layer8_size_thickness-kSpace)/2.,tracker_layer8_size_z/2.,0.,2.*TMath::Pi());
  // -----
  tracker_layer8_logical_
    = new G4LogicalVolume(tracker_layer8_solid,scintillator,"tracker_layer8_logical");
  // -----
  auto tracker_layer8_physical
    = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),tracker_layer8_logical_,"tracker_layer8_logical",
        magnetic_logical_,false,0,checkOverlaps);
  // ====================================================================================================


  // NCbarrel_layer5 ====================================================================================
  auto ncbarrel_layer5_size_thickness = carbon_thickness;
  auto ncbarrel_layer5_size_r = tracker_layer8_size_r + tracker_layer8_size_thickness/2. + ncbarrel_layer5_size_thickness/2.;
  auto ncbarrel_layer5_size_z = 2570.*mm;
  // -----
  auto ncbarrel_layer5_solid
    = new G4Tubs("ncbarrel_layer5_solid",ncbarrel_layer5_size_r-(ncbarrel_layer5_size_thickness-kSpace)/2.,ncbarrel_layer5_size_r+(ncbarrel_layer5_size_thickness-kSpace)/2.,ncbarrel_layer5_size_z/2.,0.,2.*TMath::Pi());
  // -----
  ncbarrel_layer5_logical_
    = new G4LogicalVolume(ncbarrel_layer5_solid,scintillator,"ncbarrel_layer5_logical");
  // -----
  auto ncbarrel_layer5_physical
    = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),ncbarrel_layer5_logical_,"ncbarrel_layer5_logical",
        magnetic_logical_,false,0,checkOverlaps);
  // ====================================================================================================


  // Tracker_layer9 =====================================================================================
  auto tracker_layer9_size_thickness = tracker_thickness;
  auto tracker_layer9_size_r = ncbarrel_layer5_size_r + ncbarrel_layer5_size_thickness/2. + tracker_layer9_size_thickness/2.;
  auto tracker_layer9_size_z = 2570.*mm;
  // -----
  auto tracker_layer9_solid
    = new G4Tubs("tracker_layer9_solid",tracker_layer9_size_r-(tracker_layer9_size_thickness-kSpace)/2.,tracker_layer9_size_r+(tracker_layer9_size_thickness-kSpace)/2.,tracker_layer9_size_z/2.,0.,2.*TMath::Pi());
  // -----
  tracker_layer9_logical_
    = new G4LogicalVolume(tracker_layer9_solid,scintillator,"tracker_layer9_logical");
  // -----
  auto tracker_layer9_physical
    = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),tracker_layer9_logical_,"tracker_layer9_logical",
        magnetic_logical_,false,0,checkOverlaps);
  // ====================================================================================================


  // Space ==============================================================================================
  auto space_layer9to10_size_thickness = space_thickness;
  // ====================================================================================================


  // Tracker_layer10 =====================================================================================
  auto tracker_layer10_size_thickness = tracker_thickness;
  auto tracker_layer10_size_r = tracker_layer9_size_r + tracker_layer9_size_thickness/2. + space_layer9to10_size_thickness + tracker_layer10_size_thickness/2.;
  auto tracker_layer10_size_z = 2570.*mm;
  // -----
  auto tracker_layer10_solid
    = new G4Tubs("tracker_layer10_solid",tracker_layer10_size_r-(tracker_layer10_size_thickness-kSpace)/2.,tracker_layer10_size_r+(tracker_layer10_size_thickness-kSpace)/2.,tracker_layer10_size_z/2.,0.,2.*TMath::Pi());
  // -----
  tracker_layer10_logical_
    = new G4LogicalVolume(tracker_layer10_solid,scintillator,"tracker_layer10_logical");
  // -----
  auto tracker_layer10_physical
    = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),tracker_layer10_logical_,"tracker_layer10_logical",
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
  visAttributes = new G4VisAttributes(Colors::Scintillator());
  ncbarrel_layer3_logical_->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  // -----
  visAttributes = new G4VisAttributes(Colors::Scintillator());
  ncbarrel_layer4_logical_->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  // -----
  visAttributes = new G4VisAttributes(Colors::Scintillator());
  ncbarrel_layer5_logical_->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  // -----
  visAttributes = new G4VisAttributes(Colors::WirePlane());
  tracker_layer1_logical_->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  // -----
  visAttributes = new G4VisAttributes(Colors::WirePlane());
  tracker_layer2_logical_->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  // -----
  visAttributes = new G4VisAttributes(Colors::WirePlane());
  tracker_layer3_logical_->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  // -----
  visAttributes = new G4VisAttributes(Colors::WirePlane());
  tracker_layer4_logical_->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  // -----
  visAttributes = new G4VisAttributes(Colors::WirePlane());
  tracker_layer5_logical_->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  // -----
  visAttributes = new G4VisAttributes(Colors::WirePlane());
  tracker_layer6_logical_->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  // -----
  visAttributes = new G4VisAttributes(Colors::WirePlane());
  tracker_layer7_logical_->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  // -----
  visAttributes = new G4VisAttributes(Colors::WirePlane());
  tracker_layer8_logical_->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  // -----
  visAttributes = new G4VisAttributes(Colors::WirePlane());
  tracker_layer9_logical_->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  // -----
  visAttributes = new G4VisAttributes(Colors::WirePlane());
  tracker_layer10_logical_->SetVisAttributes(visAttributes);
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
  // -----
  auto tracker_layer3 = new DriftChamberSD(SDname="/tracker_layer3");
  sdManager->AddNewDetector(tracker_layer3);
  tracker_layer3_logical_->SetSensitiveDetector(tracker_layer3);
  // -----
  auto tracker_layer4 = new DriftChamberSD(SDname="/tracker_layer4");
  sdManager->AddNewDetector(tracker_layer4);
  tracker_layer4_logical_->SetSensitiveDetector(tracker_layer4);
  // -----
  auto tracker_layer5 = new DriftChamberSD(SDname="/tracker_layer5");
  sdManager->AddNewDetector(tracker_layer5);
  tracker_layer5_logical_->SetSensitiveDetector(tracker_layer5);
  // -----
  auto tracker_layer6 = new DriftChamberSD(SDname="/tracker_layer6");
  sdManager->AddNewDetector(tracker_layer6);
  tracker_layer6_logical_->SetSensitiveDetector(tracker_layer6);
  // -----
  auto tracker_layer7 = new DriftChamberSD(SDname="/tracker_layer7");
  sdManager->AddNewDetector(tracker_layer7);
  tracker_layer7_logical_->SetSensitiveDetector(tracker_layer7);
  // -----
  auto tracker_layer8 = new DriftChamberSD(SDname="/tracker_layer8");
  sdManager->AddNewDetector(tracker_layer8);
  tracker_layer8_logical_->SetSensitiveDetector(tracker_layer8);
  // -----
  auto tracker_layer9 = new DriftChamberSD(SDname="/tracker_layer9");
  sdManager->AddNewDetector(tracker_layer9);
  tracker_layer9_logical_->SetSensitiveDetector(tracker_layer9);
  // -----
  auto tracker_layer10 = new DriftChamberSD(SDname="/tracker_layer10");
  sdManager->AddNewDetector(tracker_layer10);
  tracker_layer10_logical_->SetSensitiveDetector(tracker_layer10);

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
