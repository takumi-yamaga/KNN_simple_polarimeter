/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "TMath.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),     
  particlegun_(nullptr), messenger_(nullptr), 
  proton_(nullptr),
  momentum_(600.*MeV),
  randomize_primary_(false)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particlegun_;
  delete messenger_;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::DefineCommands()
{
  // Define /KNN_simple_polarimeter/generator command directory using generic messenger class
  messenger_ 
    = new G4GenericMessenger(this, 
        "/KNN_simple_polarimeter/generator/", 
        "Primary generator control");

  // momentum command
  auto& momentumCmd
    = messenger_->DeclarePropertyWithUnit("momentum", "GeV", momentum_, 
        "Mean momentum of primaries.");
  momentumCmd.SetParameterName("p", true);
  momentumCmd.SetRange("p>=0.");                                
  momentumCmd.SetDefaultValue("1.");

  // randomizePrimary command
  auto& randomCmd
    = messenger_->DeclareProperty("randomizePrimary", randomize_primary_);
  G4String guidance
    = "Boolean flag for randomizing primary particle types.\n";   
  guidance
    += "In case this flag is false, you can select the primary particle\n";
  guidance += "  with /gun/particle command.";                               
  randomCmd.SetGuidance(guidance);
  randomCmd.SetParameterName("flg", true);
  randomCmd.SetDefaultValue("true");
}

//..oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
