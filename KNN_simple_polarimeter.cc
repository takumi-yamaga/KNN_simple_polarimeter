//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file KNN_simple_polarimeter.cc
/// \brief Main program of KNN_simple_polarimeter

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "Analysis.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  std::string physics_list_name = "QGSP_BERT_HP";
  std::string out_file_name = "root/mc_out";
  std::string macro_name = "macro/run.mac";

  std::istringstream iss;
  std::cout<<"argc : "<<argc<<std::endl;
  for (int i = 0 ; i < argc ; i++) {
    std::string arg = argv[i];
    std::cout<<"argv["<<i<<"] : "<<argv[i]<<std::endl;
    iss.str("");
    iss.clear();
    if (arg.substr(0, 10) == "--physics=") {
      iss.str(arg.substr(10));
      iss >> physics_list_name;
    }
    else if (arg.substr(0, 10) == "--outfile=") {
      iss.str(arg.substr(10));
      iss >> out_file_name;
    }
    else if (arg.substr(0, 8) == "--macro=") {
      iss.str(arg.substr(8));
      iss >> macro_name;
    }
  }
  std::cout<<"#############################################################"<<std::endl;
  std::cout<<"--- optional input files [argument] ---"<<std::endl;
  std::cout<<"OutFileName       [--outfile=] = "<<out_file_name<<std::endl;
  std::cout<<"PhysicsListName   [--physics=] = "<<physics_list_name<<std::endl;
  std::cout<<"MacroFileName     [--macro=]   = "<<macro_name<<std::endl;
  std::cout<<"#############################################################"<<std::endl;
    

  // Detect interactive mode (if no arguments) and define UI session
  //
  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }

  // Construct the default run manager
  //
#ifdef G4MULTITHREADED
  auto runManager = new G4MTRunManager;
#else
  auto runManager = new G4RunManager;
#endif

  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFileName(out_file_name.data());

  // Mandatory user initialization classes
  runManager->SetUserInitialization(new DetectorConstruction);

  auto physicslist = new PhysicsList();
  physicslist->AddPhysicsList(physics_list_name.data());
  runManager->SetUserInitialization(physicslist);

  // User action initialization
  runManager->SetUserInitialization(new ActionInitialization());

  // Visualization manager construction
  auto visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  auto UImanager = G4UImanager::GetUIpointer();

  if ( !ui ) {
    // execute an argument macro file if exist
    G4String command = "/control/execute ";
    G4String fileName = macro_name.data();
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    UImanager->ApplyCommand("/control/execute macro/init_vis.mac");
    if (ui->IsGUI()) {
      UImanager->ApplyCommand("/control/execute macro/gui.mac");
    }     
    // start interactive session
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !

  delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
