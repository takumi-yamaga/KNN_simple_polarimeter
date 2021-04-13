/// \brief Implementation of the EventAction class

#include "EventAction.hh"
#include "DriftChamberHit.hh"
#include "Constants.hh"
#include "Analysis.hh"
#include "RootFileManager.hh"

#include "G4Event.hh"
#include "G4Trajectory.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

using std::array;
using std::vector;


namespace {

// Utility function which finds a hit collection with the given Id
// and print warnings if not found 
G4VHitsCollection* GetHC(const G4Event* event, G4int collId) {
  auto hce = event->GetHCofThisEvent();
  if (!hce) {
      G4ExceptionDescription msg;
      msg << "No hits collection of this event found." << G4endl; 
      G4Exception("EventAction::EndOfEventAction()",
                  "Code001", JustWarning, msg);
      return nullptr;
  }

  auto hc = hce->GetHC(collId);
  if ( ! hc) {
    G4ExceptionDescription msg;
    msg << "Hits collection " << collId << " of this event not found." << G4endl; 
    G4Exception("EventAction::EndOfEventAction()",
                "Code001", JustWarning, msg);
  }
  return hc;  
}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
: G4UserEventAction(), 
  driftchamber_hitscollection_ids_{{ -1, -1, -1 }}
{
  // set printing per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
  // Find hit collections and histogram Ids by names (just once)
  // and save them in the data members of this class

  // driftchamber
  if (driftchamber_hitscollection_ids_[0] == -1) {
    auto sd_manager = G4SDManager::GetSDMpointer();

    for (auto i_driftchamber = 0; i_driftchamber < Driftchamber::kTotalNumber; ++i_driftchamber) {
      auto collection_name = Driftchamber::kDetectorNames[i_driftchamber];
      collection_name.append("/driftchamber_hitscollection");
      driftchamber_hitscollection_ids_[i_driftchamber]
        = sd_manager->GetCollectionID(collection_name);
    }
  }

}     

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::EndOfEventAction(const G4Event* event)
{
  // get instance of RootFileManager
  RootFileManager& rootfile_manager = RootFileManager::Instance();

  // check event id
  int event_id = event->GetEventID();
  rootfile_manager.FillNtupleIColumn("tree_event","event_id",event_id);

  // primary vertex
  int number_of_primaries = event->GetNumberOfPrimaryVertex();
  rootfile_manager.FillNtupleIColumn("tree_event","number_of_primaries",number_of_primaries);
  for(G4int i_primary=0; i_primary<number_of_primaries; ++i_primary){
    auto primary_vertex = event->GetPrimaryVertex(i_primary);
    auto primary = primary_vertex->GetPrimary();
    if(primary_vertex&&primary){
      rootfile_manager.FillNtupleIColumn("tree_primary","event_id",event_id);
      // -----
      int particle_id = primary->GetPDGcode();
      rootfile_manager.FillNtupleIColumn("tree_primary","particle_id",particle_id);
      // -----
      G4ThreeVector position = primary_vertex->GetPosition();
      rootfile_manager.FillNtupleFColumn("tree_primary","position_x",(float)position.x());
      rootfile_manager.FillNtupleFColumn("tree_primary","position_y",(float)position.y());
      rootfile_manager.FillNtupleFColumn("tree_primary","position_z",(float)position.z());
      // -----
      G4ThreeVector momentum = primary->GetMomentum();
      rootfile_manager.FillNtupleFColumn("tree_primary","momentum_x",(float)momentum.x());
      rootfile_manager.FillNtupleFColumn("tree_primary","momentum_y",(float)momentum.y());
      rootfile_manager.FillNtupleFColumn("tree_primary","momentum_z",(float)momentum.z());
      // -----
      G4ThreeVector spin = primary->GetPolarization();
      rootfile_manager.FillNtupleFColumn("tree_primary","spin_x",(float)spin.x());
      rootfile_manager.FillNtupleFColumn("tree_primary","spin_y",(float)spin.y());
      rootfile_manager.FillNtupleFColumn("tree_primary","spin_z",(float)spin.z());
      // -----
      rootfile_manager.AddNtupleRow("tree_primary");
    }
  }

  // trajectory
  auto trajectories = event->GetTrajectoryContainer();
  int number_of_trajectories = 0;
  if(trajectories){
    int total_trajectories = (int)trajectories->size();
    for(int i_trajectory=0; i_trajectory<total_trajectories; ++i_trajectory){
      G4Trajectory* trajectory = (G4Trajectory*)(*trajectories)[i_trajectory];
      if(trajectory){
        rootfile_manager.FillNtupleIColumn("tree_trajectory","event_id",event_id);
        // -----
        int track_id = trajectory->GetTrackID();
        rootfile_manager.FillNtupleIColumn("tree_trajectory","track_id",track_id);
        // -----
        int parent_id = trajectory->GetParentID();
        rootfile_manager.FillNtupleIColumn("tree_trajectory","parent_id",parent_id);
        // -----
        int particle_id = trajectory->GetPDGEncoding();
        rootfile_manager.FillNtupleIColumn("tree_trajectory","particle_id",particle_id);
        // -----
        G4ThreeVector initial_position = trajectory->GetPoint(0)->GetPosition();
        rootfile_manager.FillNtupleFColumn("tree_trajectory","initial_position_x",(float)initial_position.x());
        rootfile_manager.FillNtupleFColumn("tree_trajectory","initial_position_y",(float)initial_position.y());
        rootfile_manager.FillNtupleFColumn("tree_trajectory","initial_position_z",(float)initial_position.z());
        // -----
        G4ThreeVector initial_momentum = trajectory->GetInitialMomentum();
        rootfile_manager.FillNtupleFColumn("tree_trajectory","initial_momentum_x",(float)initial_momentum.x());
        rootfile_manager.FillNtupleFColumn("tree_trajectory","initial_momentum_y",(float)initial_momentum.y());
        rootfile_manager.FillNtupleFColumn("tree_trajectory","initial_momentum_z",(float)initial_momentum.z());
        // -----
        if(particle_id==3122||particle_id==2212||particle_id==-211){
          rootfile_manager.AddNtupleRow("tree_trajectory");
          number_of_trajectories++;
        }
      }
    }
  }
  rootfile_manager.FillNtupleIColumn("tree_event","number_of_trajectories",number_of_trajectories);


  // ======================================================
  // Driftchamber =========================================
  // ======================================================

  // cdc
  auto cdc_hitscollection = GetHC(event,driftchamber_hitscollection_ids_[Driftchamber::kCDCId]);
  auto number_of_hits_in_cdc = (int)cdc_hitscollection->GetSize();
  rootfile_manager.FillNtupleIColumn("tree_event","number_of_hits_in_cdc",number_of_hits_in_cdc);
  for(G4int i_hit=0; i_hit<number_of_hits_in_cdc; ++i_hit){
    DriftChamberHit* hit = (DriftChamberHit*)cdc_hitscollection->GetHit(i_hit);
    if(hit){
      rootfile_manager.FillNtupleIColumn("tree_cdc_hit","event_id",event_id);
      // -----
      int track_id = hit->GetTrackID();
      rootfile_manager.FillNtupleIColumn("tree_cdc_hit","track_id",track_id);
      // -----
      int parent_id = hit->GetParentID();
      rootfile_manager.FillNtupleIColumn("tree_cdc_hit","parent_id",parent_id);
      // -----
      int particle_id = hit->GetParticleID();
      rootfile_manager.FillNtupleIColumn("tree_cdc_hit","particle_id",particle_id);
      // -----
      int layer_id = hit->GetLayerID();
      rootfile_manager.FillNtupleIColumn("tree_cdc_hit","layer_id",layer_id);
      // -----
      float hit_time = (float)hit->GetHitTime();
      rootfile_manager.FillNtupleFColumn("tree_cdc_hit","hit_time",hit_time);
      // -----
      G4ThreeVector hit_position = hit->GetGlobalPosition();
      rootfile_manager.FillNtupleFColumn("tree_cdc_hit","hit_position_x",(float)hit_position.x());
      rootfile_manager.FillNtupleFColumn("tree_cdc_hit","hit_position_y",(float)hit_position.y());
      rootfile_manager.FillNtupleFColumn("tree_cdc_hit","hit_position_z",(float)hit_position.z());
      // -----
      G4ThreeVector momentum = hit->GetMomentum();
      rootfile_manager.FillNtupleFColumn("tree_cdc_hit","momentum_x",(float)momentum.x());
      rootfile_manager.FillNtupleFColumn("tree_cdc_hit","momentum_y",(float)momentum.y());
      rootfile_manager.FillNtupleFColumn("tree_cdc_hit","momentum_z",(float)momentum.z());
      // -----
      int is_asymmetric_scattering = (int)hit->IsAsymmetricScattering();
      rootfile_manager.FillNtupleIColumn("tree_cdc_hit","is_asymmetric_scattering",is_asymmetric_scattering);
      // -----
      rootfile_manager.AddNtupleRow("tree_cdc_hit");
    }
  }
  // -----

  // tracker_layer1
  auto tracker_layer1_hitscollection = GetHC(event,driftchamber_hitscollection_ids_[Driftchamber::kTrackerLayer1Id]);
  auto number_of_hits_in_tracker_layer1 = (int)tracker_layer1_hitscollection->GetSize();
  rootfile_manager.FillNtupleIColumn("tree_event","number_of_hits_in_tracker_layer1",number_of_hits_in_tracker_layer1);
  for(G4int i_hit=0; i_hit<number_of_hits_in_tracker_layer1; ++i_hit){
    DriftChamberHit* hit = (DriftChamberHit*)tracker_layer1_hitscollection->GetHit(i_hit);
    if(hit){
      rootfile_manager.FillNtupleIColumn("tree_tracker_layer1_hit","event_id",event_id);
      // -----
      int track_id = hit->GetTrackID();
      rootfile_manager.FillNtupleIColumn("tree_tracker_layer1_hit","track_id",track_id);
      // -----
      int parent_id = hit->GetParentID();
      rootfile_manager.FillNtupleIColumn("tree_tracker_layer1_hit","parent_id",parent_id);
      // -----
      int particle_id = hit->GetParticleID();
      rootfile_manager.FillNtupleIColumn("tree_tracker_layer1_hit","particle_id",particle_id);
      // -----
      int layer_id = hit->GetLayerID();
      rootfile_manager.FillNtupleIColumn("tree_tracker_layer1_hit","layer_id",layer_id);
      // -----
      float hit_time = (float)hit->GetHitTime();
      rootfile_manager.FillNtupleFColumn("tree_tracker_layer1_hit","hit_time",hit_time);
      // -----
      G4ThreeVector hit_position = hit->GetGlobalPosition();
      rootfile_manager.FillNtupleFColumn("tree_tracker_layer1_hit","hit_position_x",(float)hit_position.x());
      rootfile_manager.FillNtupleFColumn("tree_tracker_layer1_hit","hit_position_y",(float)hit_position.y());
      rootfile_manager.FillNtupleFColumn("tree_tracker_layer1_hit","hit_position_z",(float)hit_position.z());
      // -----
      G4ThreeVector momentum = hit->GetMomentum();
      rootfile_manager.FillNtupleFColumn("tree_tracker_layer1_hit","momentum_x",(float)momentum.x());
      rootfile_manager.FillNtupleFColumn("tree_tracker_layer1_hit","momentum_y",(float)momentum.y());
      rootfile_manager.FillNtupleFColumn("tree_tracker_layer1_hit","momentum_z",(float)momentum.z());
      // -----
      int is_asymmetric_scattering = (int)hit->IsAsymmetricScattering();
      rootfile_manager.FillNtupleIColumn("tree_tracker_layer1_hit","is_asymmetric_scattering",is_asymmetric_scattering);
      // -----
      rootfile_manager.AddNtupleRow("tree_tracker_layer1_hit");
    }
  }
  // -----

  // tracker_layer2
  auto tracker_layer2_hitscollection = GetHC(event,driftchamber_hitscollection_ids_[Driftchamber::kTrackerLayer2Id]);
  auto number_of_hits_in_tracker_layer2 = (int)tracker_layer2_hitscollection->GetSize();
  rootfile_manager.FillNtupleIColumn("tree_event","number_of_hits_in_tracker_layer2",number_of_hits_in_tracker_layer2);
  for(G4int i_hit=0; i_hit<number_of_hits_in_tracker_layer2; ++i_hit){
    DriftChamberHit* hit = (DriftChamberHit*)tracker_layer2_hitscollection->GetHit(i_hit);
    if(hit){
      rootfile_manager.FillNtupleIColumn("tree_tracker_layer2_hit","event_id",event_id);
      // -----
      int track_id = hit->GetTrackID();
      rootfile_manager.FillNtupleIColumn("tree_tracker_layer2_hit","track_id",track_id);
      // -----
      int parent_id = hit->GetParentID();
      rootfile_manager.FillNtupleIColumn("tree_tracker_layer2_hit","parent_id",parent_id);
      // -----
      int particle_id = hit->GetParticleID();
      rootfile_manager.FillNtupleIColumn("tree_tracker_layer2_hit","particle_id",particle_id);
      // -----
      int layer_id = hit->GetLayerID();
      rootfile_manager.FillNtupleIColumn("tree_tracker_layer2_hit","layer_id",layer_id);
      // -----
      float hit_time = (float)hit->GetHitTime();
      rootfile_manager.FillNtupleFColumn("tree_tracker_layer2_hit","hit_time",hit_time);
      // -----
      G4ThreeVector hit_position = hit->GetGlobalPosition();
      rootfile_manager.FillNtupleFColumn("tree_tracker_layer2_hit","hit_position_x",(float)hit_position.x());
      rootfile_manager.FillNtupleFColumn("tree_tracker_layer2_hit","hit_position_y",(float)hit_position.y());
      rootfile_manager.FillNtupleFColumn("tree_tracker_layer2_hit","hit_position_z",(float)hit_position.z());
      // -----
      G4ThreeVector momentum = hit->GetMomentum();
      rootfile_manager.FillNtupleFColumn("tree_tracker_layer2_hit","momentum_x",(float)momentum.x());
      rootfile_manager.FillNtupleFColumn("tree_tracker_layer2_hit","momentum_y",(float)momentum.y());
      rootfile_manager.FillNtupleFColumn("tree_tracker_layer2_hit","momentum_z",(float)momentum.z());
      // -----
      int is_asymmetric_scattering = (int)hit->IsAsymmetricScattering();
      rootfile_manager.FillNtupleIColumn("tree_tracker_layer2_hit","is_asymmetric_scattering",is_asymmetric_scattering);
      // -----
      rootfile_manager.AddNtupleRow("tree_tracker_layer2_hit");
    }
  }
  // -----

  rootfile_manager.AddNtupleRow("tree_event");

  // set printing per each event
  if(event->GetEventID()){
    G4int print_progress = (G4int)log10(event->GetEventID());
    G4RunManager::GetRunManager()->SetPrintProgress(pow(10,print_progress));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
