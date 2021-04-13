/// \brief Definition of the PrimaryGeneratorAction class

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "KppGenerator.hh"
#include "Constants.hh"

#include "TVector3.h"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

class G4GenericMessenger;
class G4Event;
class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();
    virtual ~PrimaryGeneratorAction();
    
    void ConstructParticleGun();

    virtual void GeneratePrimaries(G4Event*);
    inline void GenerateKpp(G4Event*);
    
  private:
    void DefineCommands();

    G4GenericMessenger* messenger_;

    KppGenerator* kpp_generator_;

    G4int total_primaries_;
    std::vector <G4ParticleGun*> particleguns_;
    std::vector <G4String> particle_names_;
    std::vector <G4double> particle_masses_;

    G4String beam_particle_name_;
    TLorentzVector lv_beam_;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void PrimaryGeneratorAction::GenerateKpp(G4Event* event){

  // generate particles -------------------------
  kpp_generator_->Generate(lv_beam_);

  // for lambda
  auto* particlegun_lambda = particleguns_[0];
  TLorentzVector lv_lambda = kpp_generator_->LvLambda();
  TVector3 vec_lambda_spin_direction = kpp_generator_->VecLambdaSpin();
  G4ThreeVector g4vec_lambda_momentum_direction = ConvertVector(lv_lambda.Vect().Unit());
  G4ThreeVector g4vec_lambda_spin_direction = ConvertVector(vec_lambda_spin_direction);
  particlegun_lambda->SetParticleEnergy((lv_lambda.E()-lv_lambda.M())*GeV);
  particlegun_lambda->SetParticleMomentumDirection(g4vec_lambda_momentum_direction);
  particlegun_lambda->SetParticlePolarization(g4vec_lambda_spin_direction);
  particlegun_lambda->GeneratePrimaryVertex(event);

  // for proton
  auto* particlegun_proton = particleguns_[1];
  TLorentzVector lv_proton = kpp_generator_->LvProton();
  TVector3 vec_proton_spin_direction = kpp_generator_->VecProtonSpin();
  G4ThreeVector g4vec_proton_momentum_direction = ConvertVector(lv_proton.Vect().Unit());
  G4ThreeVector g4vec_proton_spin_direction = ConvertVector(vec_proton_spin_direction);
  particlegun_proton->SetParticleEnergy((lv_proton.E()-lv_proton.M())*GeV);
  particlegun_proton->SetParticleMomentumDirection(g4vec_proton_momentum_direction);
  particlegun_proton->SetParticlePolarization(g4vec_proton_spin_direction);
  particlegun_proton->GeneratePrimaryVertex(event);

  // --------------------------------------------

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif
