// src/StackingAction.cc

#include "StackingAction.hh"
#include "TrackInformation.hh"

#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4SystemOfUnits.hh"    
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
StackingAction::StackingAction() : G4UserStackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
StackingAction::~StackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
G4ClassificationOfNewTrack StackingAction::ClassifyNewTrack(const G4Track* track)
{
  G4ClassificationOfNewTrack classification = fUrgent;

  if(track->GetParentID() == 0){
    TrackInformation* track_information = new TrackInformation(track);
    track->SetUserInformation(track_information);
  }
  else{
    TrackInformation* track_information = (TrackInformation*)track->GetUserInformation();
    track_information->PushTrackInformation(track);
  }

  return classification;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
void StackingAction::NewStage() 
  // called when urgent stack is empty
  // processing for wait stack
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......     
void StackingAction::PrepareNewEvent()
{}
