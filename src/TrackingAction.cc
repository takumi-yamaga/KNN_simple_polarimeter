/// \file src/TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
//

#include "TrackingAction.hh"
#include "TrackInformation.hh"

#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4Trajectory.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
TrackingAction::TrackingAction()
  :G4UserTrackingAction()
{;}

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
  fpTrackingManager->SetStoreTrajectory(true);
  fpTrackingManager->SetTrajectory(new G4Trajectory(track));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
void TrackingAction::PostUserTrackingAction(const G4Track* track)
{

  G4TrackVector* daughter_tracks = fpTrackingManager->GimmeSecondaries();
  if(daughter_tracks)
  {
    TrackInformation* current_track_information = (TrackInformation*)(track->GetUserInformation());
    G4int total_daughters = (G4int)( daughter_tracks->size() );
    if(total_daughters){
      for(G4int i_daughter=0; i_daughter<total_daughters; ++i_daughter){
        G4Track* daughter_track = (*daughter_tracks)[i_daughter];
        TrackInformation* daughter_track_information = new TrackInformation(current_track_information);
        daughter_track->SetUserInformation(daughter_track_information);
      }
    }
  }
}


