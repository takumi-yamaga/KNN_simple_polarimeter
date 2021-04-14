/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "Analysis.hh"
#include "RootFileManager.hh"

#include "time.h"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
 : G4UserRunAction()
{ 
  RootFileManager& rootfile_manager = RootFileManager::Instance();

  // ====================================================================================================
  rootfile_manager.CreateNtuple("tree_event", "Tree for event information");
  // -----
  rootfile_manager.CreateNtupleIColumn("event_id");
  rootfile_manager.CreateNtupleIColumn("number_of_primaries");
  rootfile_manager.CreateNtupleIColumn("number_of_trajectories");
  rootfile_manager.CreateNtupleIColumn("number_of_hits_in_cdc");
  rootfile_manager.CreateNtupleIColumn("number_of_hits_in_tracker_layer1");
  rootfile_manager.CreateNtupleIColumn("number_of_hits_in_tracker_layer2");
  // -----
  rootfile_manager.FinishNtuple();
  // ====================================================================================================


  // ====================================================================================================
  rootfile_manager.CreateNtuple("tree_primary", "Tree for primary vertex information");
  // -----
  rootfile_manager.CreateNtupleIColumn("event_id");
  rootfile_manager.CreateNtupleIColumn("particle_id");
  rootfile_manager.CreateNtupleFColumn("position_x");
  rootfile_manager.CreateNtupleFColumn("position_y");
  rootfile_manager.CreateNtupleFColumn("position_z");
  rootfile_manager.CreateNtupleFColumn("momentum_x");
  rootfile_manager.CreateNtupleFColumn("momentum_y");
  rootfile_manager.CreateNtupleFColumn("momentum_z");
  rootfile_manager.CreateNtupleFColumn("spin_x");
  rootfile_manager.CreateNtupleFColumn("spin_y");
  rootfile_manager.CreateNtupleFColumn("spin_z");
  // -----
  rootfile_manager.FinishNtuple();
  // ====================================================================================================


  // ====================================================================================================
  rootfile_manager.CreateNtuple("tree_trajectory", "Tree for trajectory information");
  // -----
  rootfile_manager.CreateNtupleIColumn("event_id");
  rootfile_manager.CreateNtupleIColumn("track_id");
  rootfile_manager.CreateNtupleIColumn("parent_id");
  rootfile_manager.CreateNtupleIColumn("particle_id");
  rootfile_manager.CreateNtupleFColumn("initial_position_x");
  rootfile_manager.CreateNtupleFColumn("initial_position_y");
  rootfile_manager.CreateNtupleFColumn("initial_position_z");
  rootfile_manager.CreateNtupleFColumn("initial_momentum_x");
  rootfile_manager.CreateNtupleFColumn("initial_momentum_y");
  rootfile_manager.CreateNtupleFColumn("initial_momentum_z");
  // -----
  rootfile_manager.FinishNtuple();
  // ====================================================================================================


  // ====================================================================================================
  rootfile_manager.CreateNtuple("tree_cdc_hit", "Tree for CDC hit information");
  // -----
  rootfile_manager.CreateNtupleIColumn("event_id");
  rootfile_manager.CreateNtupleIColumn("track_id");
  rootfile_manager.CreateNtupleIColumn("parent_id");
  rootfile_manager.CreateNtupleIColumn("particle_id");
  rootfile_manager.CreateNtupleIColumn("layer_id");
  rootfile_manager.CreateNtupleFColumn("hit_time");
  rootfile_manager.CreateNtupleFColumn("hit_position_x");
  rootfile_manager.CreateNtupleFColumn("hit_position_y");
  rootfile_manager.CreateNtupleFColumn("hit_position_z");
  rootfile_manager.CreateNtupleFColumn("momentum_x");
  rootfile_manager.CreateNtupleFColumn("momentum_y");
  rootfile_manager.CreateNtupleFColumn("momentum_z");
  rootfile_manager.CreateNtupleIColumn("is_asymmetric_scattering");
  // -----
  rootfile_manager.FinishNtuple();
  // ====================================================================================================


  // ====================================================================================================
  rootfile_manager.CreateNtuple("tree_tracker_layer1_hit", "Tree for TrackerLayer1 hit information");
  // -----
  rootfile_manager.CreateNtupleIColumn("event_id");
  rootfile_manager.CreateNtupleIColumn("track_id");
  rootfile_manager.CreateNtupleIColumn("parent_id");
  rootfile_manager.CreateNtupleIColumn("particle_id");
  rootfile_manager.CreateNtupleIColumn("layer_id");
  rootfile_manager.CreateNtupleFColumn("hit_time");
  rootfile_manager.CreateNtupleFColumn("hit_position_x");
  rootfile_manager.CreateNtupleFColumn("hit_position_y");
  rootfile_manager.CreateNtupleFColumn("hit_position_z");
  rootfile_manager.CreateNtupleFColumn("momentum_x");
  rootfile_manager.CreateNtupleFColumn("momentum_y");
  rootfile_manager.CreateNtupleFColumn("momentum_z");
  rootfile_manager.CreateNtupleIColumn("is_asymmetric_scattering");
  // -----
  rootfile_manager.FinishNtuple();
  // ====================================================================================================


  // ====================================================================================================
  rootfile_manager.CreateNtuple("tree_tracker_layer2_hit", "Tree for TrackerLayer2 hit information");
  // -----
  rootfile_manager.CreateNtupleIColumn("event_id");
  rootfile_manager.CreateNtupleIColumn("track_id");
  rootfile_manager.CreateNtupleIColumn("parent_id");
  rootfile_manager.CreateNtupleIColumn("particle_id");
  rootfile_manager.CreateNtupleIColumn("layer_id");
  rootfile_manager.CreateNtupleFColumn("hit_time");
  rootfile_manager.CreateNtupleFColumn("hit_position_x");
  rootfile_manager.CreateNtupleFColumn("hit_position_y");
  rootfile_manager.CreateNtupleFColumn("hit_position_z");
  rootfile_manager.CreateNtupleFColumn("momentum_x");
  rootfile_manager.CreateNtupleFColumn("momentum_y");
  rootfile_manager.CreateNtupleFColumn("momentum_z");
  rootfile_manager.CreateNtupleIColumn("is_asymmetric_scattering");
  // -----
  rootfile_manager.FinishNtuple();
  // ====================================================================================================

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
  G4long random_seed  = time(NULL);
  G4int random_luxury = 5;
  CLHEP::HepRandom::setTheSeed(random_seed,random_luxury);

  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  G4RunManager::GetRunManager()->SetRandomNumberStoreDir("./rndm/");

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file 
  // The default file name is set in RunAction::RunAction(),
  // it can be overwritten in a macro
  analysisManager->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // save histograms & ntuple
  //
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
