// SteppingAction.cc

#include "G4ios.hh"
#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SystemOfUnits.hh"

#include "TrackInformation.hh"
#include "SteppingAction.hh"
#include "LambdaDecay.hh"
#include "ProtonScatteringWithSpin.hh"
#include "Constants.hh"

#include "TLorentzVector.h"

SteppingAction::SteppingAction()
{ 
  lambda_decay_ = new LambdaDecay();
  proton_scattering_with_spin_ = new ProtonScatteringWithSpin();
}

SteppingAction::~SteppingAction()
{ 
  delete lambda_decay_;
  delete proton_scattering_with_spin_;
}

void SteppingAction::UserSteppingAction(const G4Step* step)
{ 
  G4SteppingManager* stepping_manager = fpSteppingManager;
  G4Track* track = step->GetTrack();
  TrackInformation* track_information = (TrackInformation*)track->GetUserInformation();

  // ====================================================================================================
  // for lambda-decay
  // ====================================================================================================
  if(track->GetDefinition()->GetParticleName()=="lambda"){
    G4ThreeVector g4vec_lambda_momentum = track->GetMomentum();
    double lambda_energy = track->GetTotalEnergy();
    TLorentzVector lv_lambda(ConvertVector(g4vec_lambda_momentum),lambda_energy);
    const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();
    if(process){
      if(process->GetProcessName()=="Decay"){
        G4TrackVector* secondaries = stepping_manager -> GetfSecondary();
        int number_of_decayproducts = (int)(*secondaries).size();
        if(number_of_decayproducts==2){
          G4Track* track_decayproduct1 = (*secondaries)[0];
          G4Track* track_decayproduct2 = (*secondaries)[1];
          G4String decayproduct1_name = track_decayproduct1->GetDefinition()->GetParticleName(); 
          if(decayproduct1_name=="proton"||decayproduct1_name=="pi-"){
            G4ThreeVector g4vec_lambda_spin_direction = track_information->GetCurrentSpinDirection();
            TVector3 vec_lambda_spin_direction = ConvertVector(g4vec_lambda_spin_direction);

            lambda_decay_->Generate(lv_lambda,vec_lambda_spin_direction);
            TLorentzVector lv_proton = lambda_decay_->LvProton();
            TLorentzVector lv_pim = lambda_decay_->LvPim();
            G4ThreeVector g4vec_proton_momentum_direction = ConvertVector(lv_proton.Vect().Unit());
            TVector3 vec_proton_spin_direction = lambda_decay_->VecProtonSpin();
            G4ThreeVector g4vec_proton_spin_direction = ConvertVector(vec_proton_spin_direction);
            double proton_kinetic_energy = lv_proton.E() - lv_proton.M();
            G4ThreeVector g4vec_pim_momentum_direction = ConvertVector(lv_pim.Vect().Unit());
            double pim_kinetic_energy = lv_pim.E() - lv_pim.M();
            if(decayproduct1_name=="proton"){
              track_decayproduct1->SetMomentumDirection(g4vec_proton_momentum_direction);
              track_decayproduct1->SetKineticEnergy(proton_kinetic_energy);
              track_decayproduct1->SetPolarization(g4vec_proton_spin_direction);
              track_decayproduct2->SetMomentumDirection(g4vec_pim_momentum_direction);
              track_decayproduct2->SetKineticEnergy(pim_kinetic_energy);
            }
            else{
              track_decayproduct1->SetMomentumDirection(g4vec_pim_momentum_direction);
              track_decayproduct1->SetKineticEnergy(pim_kinetic_energy);
              track_decayproduct2->SetMomentumDirection(g4vec_proton_momentum_direction);
              track_decayproduct2->SetKineticEnergy(proton_kinetic_energy);
              track_decayproduct2->SetPolarization(g4vec_proton_spin_direction);
            }
          }
        }
      }
    }
  }
  // ====================================================================================================



  // ====================================================================================================
  // for proton scattering asymmetry
  // ====================================================================================================
  if(track->GetDefinition()->GetParticleName()=="proton"){

    // changing momentum direction
    auto pre_steppoint = step->GetPreStepPoint();
    auto post_steppoint = step->GetPostStepPoint();
    if(pre_steppoint&&post_steppoint){
      const G4VProcess* post_process = post_steppoint->GetProcessDefinedStep();
      if(post_process){
        auto post_physical_volume = post_steppoint->GetPhysicalVolume();
        if(post_physical_volume){
          auto post_physical_volume_name = post_physical_volume->GetName();
          if(post_physical_volume_name.contains("nc")){
            G4String post_process_name = post_process->GetProcessName();

            // ==================== hadElastic ====================
            if(post_process_name=="hadElastic"){ // primary proton will be scattered.
              TVector3 vec_proton_momentum_initial = ConvertVector(pre_steppoint->GetMomentum());
              TVector3 vec_proton_momentum_final_pre = ConvertVector(post_steppoint->GetMomentum());

              G4ThreeVector g4vec_proton_spin_direction = track_information->GetCurrentSpinDirection();
              TVector3 vec_proton_spin_direction = ConvertVector(g4vec_proton_spin_direction);
              Double_t theta_scattering = acos(vec_proton_momentum_final_pre.Unit().Dot(vec_proton_momentum_initial.Unit()));
              proton_scattering_with_spin_->Generate(vec_proton_momentum_initial,vec_proton_momentum_final_pre,vec_proton_spin_direction,theta_scattering);
              TVector3 vec_proton_momentum_final = proton_scattering_with_spin_->VecProtonMomentumFinal();

              G4ThreeVector momentum_direction_after_scattering = ConvertVector(vec_proton_momentum_final.Unit());
              track->SetMomentumDirection(momentum_direction_after_scattering);
              step->GetPostStepPoint()->SetMomentumDirection(momentum_direction_after_scattering);
            }
            // ==================== hadElastic ====================

            // ==================== proton Inelastic ====================
            else if(post_process_name=="protonInelastic"){ // primary proton will be vanished and secondaries to be created.
              TVector3 vec_proton_momentum_initial = ConvertVector(pre_steppoint->GetMomentum());
              TVector3 vec_proton_momentum_final_pre(0.,0.,0.);

              G4TrackVector* secondaries = stepping_manager -> GetfSecondary();
              G4Track* secondary_proton = nullptr;
              int number_of_secondaries = (int)(*secondaries).size();
              bool is_there_proton = false;
              for(int i_secondary=0; i_secondary<number_of_secondaries; ++i_secondary){
                // checking secondary protons and pickup the highest energy proton
                auto secondary = (*secondaries)[i_secondary];
                if(secondary->GetDefinition()->GetParticleName() == "proton"){
                  if(vec_proton_momentum_final_pre.Mag()<secondary->GetMomentum().mag()){
                    vec_proton_momentum_final_pre = ConvertVector(secondary->GetMomentum());
                    secondary_proton = secondary;
                    is_there_proton = true;
                  }
                }
              }
              if(is_there_proton){
                G4ThreeVector g4vec_proton_spin_direction = track_information->GetCurrentSpinDirection();
                TVector3 vec_proton_spin_direction = ConvertVector(g4vec_proton_spin_direction);
                Double_t theta_scattering = acos(vec_proton_momentum_final_pre.Unit().Dot(vec_proton_momentum_initial.Unit()));
                proton_scattering_with_spin_->Generate(vec_proton_momentum_initial,vec_proton_momentum_final_pre,vec_proton_spin_direction,theta_scattering);
                TVector3 vec_proton_momentum_final = proton_scattering_with_spin_->VecProtonMomentumFinal();
                G4ThreeVector momentum_direction_after_scattering = ConvertVector(vec_proton_momentum_final.Unit());
                secondary_proton->SetMomentumDirection(momentum_direction_after_scattering);
              }
            }
            // ==================== proton Inelastic ====================

          }
        }
      }
    }
  }
  //// -----

  // ====================================================================================================


}

