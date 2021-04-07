// RootFileManager.cc

#include "Analysis.hh"
#include "RootFileManager.hh"

// -----
int RootFileManager::CreateNtuple(const char* ntuple_name, const char* ntuple_info){
  auto analysis_manager = G4AnalysisManager::Instance();
  int ntuple_id = analysis_manager->CreateNtuple(ntuple_name,ntuple_info);
  ++number_of_ntuples_;
  ntuple_ids_.push_back(ntuple_id);
  ntuple_names_.push_back(std::string(ntuple_name));

  number_of_columns_ = 0;
  columns_ids_.push_back(std::vector<int>());
  columns_names_.push_back(std::vector<std::string>());

  return ntuple_id;
}
// -----
void RootFileManager::FinishNtuple(){
  auto analysis_manager = G4AnalysisManager::Instance();
  analysis_manager->FinishNtuple();

  numbers_of_columns_.push_back(number_of_columns_);
  number_of_columns_ = 0;

  return;
}
// -----
int RootFileManager::CreateNtupleIColumn(const char* column_name){
  auto analysis_manager = G4AnalysisManager::Instance();
  int column_id = analysis_manager->CreateNtupleIColumn(column_name);
  ++number_of_columns_;
  columns_ids_[number_of_ntuples_-1].push_back(column_id);
  columns_names_[number_of_ntuples_-1].push_back(std::string(column_name));

  return column_id;
}
// -----
int RootFileManager::CreateNtupleFColumn(const char* column_name){
  auto analysis_manager = G4AnalysisManager::Instance();
  int column_id = analysis_manager->CreateNtupleFColumn(column_name);
  ++number_of_columns_;
  columns_ids_[number_of_ntuples_-1].push_back(column_id);
  columns_names_[number_of_ntuples_-1].push_back(std::string(column_name));

  return column_id;
}
// -----
int RootFileManager::CreateNtupleDColumn(const char* column_name){
  auto analysis_manager = G4AnalysisManager::Instance();
  int column_id = analysis_manager->CreateNtupleDColumn(column_name);
  ++number_of_columns_;
  columns_ids_[number_of_ntuples_-1].push_back(column_id);
  columns_names_[number_of_ntuples_-1].push_back(std::string(column_name));

  return column_id;
}
// -----
void RootFileManager::FillNtupleIColumn(const char* ntuple_name, const char* column_name, int val){
  int ntuple_id = FindNtupleId(ntuple_name);
  int column_id = FindColumnId(ntuple_name, column_name);

  auto analysis_manager = G4AnalysisManager::Instance();
  analysis_manager->FillNtupleIColumn(ntuple_id,column_id,val);

  return;
}
// -----
void RootFileManager::FillNtupleFColumn(const char* ntuple_name, const char* column_name, float val){
  int ntuple_id = FindNtupleId(ntuple_name);
  int column_id = FindColumnId(ntuple_name, column_name);

  auto analysis_manager = G4AnalysisManager::Instance();
  analysis_manager->FillNtupleFColumn(ntuple_id,column_id,val);

  return;
}
// -----
void RootFileManager::FillNtupleDColumn(const char* ntuple_name, const char* column_name, double val){
  int ntuple_id = FindNtupleId(ntuple_name);
  int column_id = FindColumnId(ntuple_name, column_name);

  auto analysis_manager = G4AnalysisManager::Instance();
  analysis_manager->FillNtupleDColumn(ntuple_id,column_id,val);

  return;
}
// -----
void RootFileManager::AddNtupleRow(const char* ntuple_name){
  int ntuple_id = FindNtupleId(ntuple_name);

  auto analysis_manager = G4AnalysisManager::Instance();
  analysis_manager->AddNtupleRow(ntuple_id);

  return;
}
