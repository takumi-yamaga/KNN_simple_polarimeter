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

  // for proton
  auto* particlegun_proton = particleguns_[1];
  TLorentzVector lv_proton = kpp_generator_->LvProton();
  G4ThreeVector vec_proton_momentum_direction = ConvertVector(lv_proton.Vect().Unit());
  particlegun_proton->SetParticleEnergy((lv_proton.E()-lv_proton.M())*GeV);
  particlegun_proton->SetParticleMomentumDirection(vec_proton_momentum_direction);
  particlegun_proton->GeneratePrimaryVertex(event);
  // --------------------------------------------

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif
