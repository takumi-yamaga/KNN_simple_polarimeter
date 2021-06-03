/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4GenericMessenger.hh"
#include "Randomize.hh"

#include "TMath.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),     
  kpp_generator_(nullptr),
  messenger_(nullptr) 
{
  // for beam
  beam_particle_name_ = "kaon-";
  auto particle_table = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* beam_particle = particle_table->FindParticle(beam_particle_name_);
  auto beam_mass = beam_particle->GetPDGMass();
  auto beam_momentum = 1.0*GeV;
  auto beam_energy = sqrt(beam_mass*beam_mass+beam_momentum*beam_momentum);
  lv_beam_ = TLorentzVector(0.,0.,beam_momentum,beam_energy);

  // for scattered
  total_primaries_=3;
  particle_names_.push_back("lambda");
  particle_names_.push_back("proton");
  particle_names_.push_back("neutron");

  // construction of particle guns
  ConstructParticleGun();

  // create KppGenerator
  kpp_generator_ = new KppGenerator(1,-1);

  // define commands for this class
  DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete kpp_generator_;
  delete messenger_;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  GenerateKpp(event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::ConstructParticleGun() 
{
  if(particleguns_.size()){
    for(auto* particlegun : particleguns_){
      delete particlegun;
    }
    particleguns_.clear();
    particleguns_.shrink_to_fit();
  }
  if(particle_masses_.size()){
    particle_masses_.clear();
    particle_masses_.shrink_to_fit();
  }
  // construct particlegun
  auto particle_table = G4ParticleTable::GetParticleTable();
  for(auto particle_name : particle_names_){
    G4ParticleDefinition* particle_definition = particle_table->FindParticle(particle_name);
    auto particle_mass = particle_definition->GetPDGMass();
    auto* particlegun = new G4ParticleGun();
    particlegun->SetParticleDefinition(particle_definition);
    particlegun->SetParticleEnergy(particle_mass);
    particlegun->SetParticleMomentum(G4ThreeVector());
    particlegun->SetParticlePosition(G4ThreeVector());
    particleguns_.push_back(particlegun);
    particle_masses_.push_back(particle_mass);
  }
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
  //auto& momentumCmd
  //  = messenger_->DeclarePropertyWithUnit("momentum", "GeV", momentum_, 
  //      "Mean momentum of primaries.");
  //momentumCmd.SetParameterName("p", true);
  //momentumCmd.SetRange("p>=0.");                                
  //momentumCmd.SetDefaultValue("1.");

  // randomizePrimary command
  //auto& randomCmd
  //  = messenger_->DeclareProperty("randomizePrimary", randomize_primary_);
  //G4String guidance
  //  = "Boolean flag for randomizing primary particle types.\n";   
  //guidance
  //  += "In case this flag is false, you can select the primary particle\n";
  //guidance += "  with /gun/particle command.";                               
  //randomCmd.SetGuidance(guidance);
  //randomCmd.SetParameterName("flg", true);
  //randomCmd.SetDefaultValue("true");
}

//..oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
