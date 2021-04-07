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

  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Default settings
  //analysisManager->SetNtupleMerging(true); // for multi threading
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFileName("mc_out");

  // Creating 1D histograms
  //analysisManager // H1-ID = 0
  //  ->CreateH1("dcin_numhit","dcin : number of hits", 10, 0., 10.);
  
  // Creating 2D histograms
  //analysisManager  // H2-ID = 0                                              
  //  ->CreateH2("dcin_hitposition_xy","dcin : hit position on x-y plane;x;y",
  //             50, -100., 100, 50, -100., 100.); 

  // Creating tree_event
  // ====================================================================================================
  rootfile_manager.CreateNtuple("tree_event", "Tree for event information");
  // -----
  rootfile_manager.CreateNtupleIColumn("event_id");
  rootfile_manager.CreateNtupleIColumn("number_of_primaries");
  rootfile_manager.CreateNtupleIColumn("number_of_hits_in_cdc");
  rootfile_manager.CreateNtupleIColumn("number_of_hits_in_tracker_layer1");
  // -----
  rootfile_manager.FinishNtuple();
  // ====================================================================================================


  // ====================================================================================================
  rootfile_manager.CreateNtuple("tree_cdc_hit", "Tree for CDC hit information");
  // -----
  rootfile_manager.CreateNtupleIColumn("event_id");
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
