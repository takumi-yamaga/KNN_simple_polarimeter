/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "Analysis.hh"

#include "time.h"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
 : G4UserRunAction()
{ 
  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Default settings
  analysisManager->SetNtupleMerging(true);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFileName("hodoscope");

  // Creating 1D histograms
  //analysisManager // H1-ID = 0
  //  ->CreateH1("dcin_numhit","dcin : number of hits", 10, 0., 10.);
  
  // Creating 2D histograms
  //analysisManager  // H2-ID = 0                                              
  //  ->CreateH2("dcin_hitposition_xy","dcin : hit position on x-y plane;x;y",
  //             50, -100., 100, 50, -100., 100.); 

  // Creating tree
  analysisManager->CreateNtuple("tree_event", "Tree for event information");

  analysisManager->CreateNtupleIColumn("event_id");

  analysisManager->CreateNtupleIColumn("dcout_nhit");       // column Id = 7
  analysisManager->CreateNtupleFColumn("dcout_position_x"); // column Id = 8
  analysisManager->CreateNtupleFColumn("dcout_position_y"); // column Id = 9
  analysisManager->CreateNtupleFColumn("dcout_position_z"); // column Id =10
  analysisManager->CreateNtupleFColumn("dcout_momentum_x"); // column Id =11
  analysisManager->CreateNtupleFColumn("dcout_momentum_y"); // column Id =12
  analysisManager->CreateNtupleFColumn("dcout_momentum_z"); // column Id =13

  analysisManager->FinishNtuple();
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
