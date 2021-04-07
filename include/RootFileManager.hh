// RootFileManager.hh

#ifndef RootFileManager_hh
#define RootFileManager_hh

#include<vector>

class RootFileManager final
{
  private:
    RootFileManager() = default;
    ~RootFileManager() = default;

  public:
    RootFileManager(const RootFileManager&) = delete;
    RootFileManager& operator=(const RootFileManager&) = delete;
    RootFileManager(const RootFileManager&&) = delete;
    RootFileManager& operator=(const RootFileManager&&) = delete;

    static RootFileManager& Instance(){
      static RootFileManager instance;
      return instance;
    }

  private:
    int number_of_ntuples_;
    std::vector<int> ntuple_ids_;
    std::vector<std::string> ntuple_names_;

    int number_of_columns_;
    std::vector<int> numbers_of_columns_;
    std::vector< std::vector<int> > columns_ids_;
    std::vector< std::vector<std::string> > columns_names_;

  public:
    int CreateNtuple(const char*, const char*);
    void FinishNtuple();
    int CreateNtupleIColumn(const char*);
    int CreateNtupleFColumn(const char*);
    int CreateNtupleDColumn(const char*);
    void FillNtupleIColumn(const char*, const char*,int);
    void FillNtupleFColumn(const char*, const char*,float);
    void FillNtupleDColumn(const char*, const char*,double);
    void AddNtupleRow(const char*);

  private:
    inline int FindNtupleId(const char*);
    inline int FindColumnId(const char*, const char*);

  public:
    inline int GetNumberOfNtuples();
    inline std::vector<int> GetNtupleIds();
    inline int GetNtupleId(int);
    inline std::vector<std::string> GetNtupleNames();
    inline std::string GetNtupleName(int);

    inline std::vector<int> GetNumbersOfColumns();
    inline int GetNumberOfColumns(int);
    inline std::vector< std::vector<int> > GetColumnsIds();
    inline std::vector<int> GetColumnIds(int);
    inline int GetColumnId(int,int);
    inline std::vector< std::vector<std::string> > GetColumnsNames();
    inline std::vector<std::string> GetColumnNames(int);
    inline std::string GetColumnName(int,int);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Inline functions
// -----
inline int RootFileManager::GetNumberOfNtuples(){
  return number_of_ntuples_;
}
// -----
inline std::vector<int> RootFileManager::GetNtupleIds(){
  return ntuple_ids_;
}
// -----
inline int RootFileManager::GetNtupleId(int id){
  if(id<number_of_ntuples_){
    return ntuple_ids_[id];
  }
  else{
    return -999;
  }
}
// -----
inline std::vector<std::string> RootFileManager::GetNtupleNames(){
  return ntuple_names_;
}
// -----
inline std::string RootFileManager::GetNtupleName(int id){
  if(id<number_of_ntuples_){
    return ntuple_names_[id];
  }
  else{
    return std::string("default");
  }
}
// -----
inline std::vector<int> RootFileManager::GetNumbersOfColumns(){
  return numbers_of_columns_;
}
// -----
inline int RootFileManager::GetNumberOfColumns(int id){
  if(id<number_of_ntuples_){
    return numbers_of_columns_[id];
  }
  else{
    return -999;
  }
}
// -----
inline std::vector< std::vector<int> > RootFileManager::GetColumnsIds(){
  return columns_ids_;
}
// -----
inline std::vector<int> RootFileManager::GetColumnIds(int id){
  if(id<number_of_ntuples_){
    return columns_ids_[id];
  }
  else{
    std::vector<int> vect;
    return vect;
  }
}
// -----
inline int RootFileManager::GetColumnId(int ntuple_id, int column_id){
  if(ntuple_id<number_of_ntuples_){
    if(column_id<numbers_of_columns_[ntuple_id]){
      return columns_ids_[ntuple_id][column_id];
    }
    else{
      return -999;
    }
  }
  else{
    return -999;
  }
}
// -----
inline std::vector< std::vector<std::string> > RootFileManager::GetColumnsNames(){
  return columns_names_;
}
// -----
inline std::vector<std::string> RootFileManager::GetColumnNames(int id){
  if(id<number_of_ntuples_){
    return columns_names_[id];
  }
  else{
    std::vector<std::string> vect;
    return vect;
  }
}
// -----
inline std::string RootFileManager::GetColumnName(int ntuple_id, int column_id){
  if(ntuple_id<number_of_ntuples_){
    if(column_id<numbers_of_columns_[ntuple_id]){
      return columns_names_[ntuple_id][column_id];
    }
    else{
      return std::string("default");
    }
  }
  else{
    return std::string("default");
  }
}
// -----
inline int RootFileManager::FindNtupleId(const char* ntuple_name){
  size_t index = 0;
  for(std::string name : ntuple_names_){
    if(name==ntuple_name){
      return ntuple_ids_[index];
    }
    ++index;
  }
  return -999;
}
// -----
inline int RootFileManager::FindColumnId(const char* ntuple_name, const char* column_name){
  int ntuple_id = FindNtupleId(ntuple_name);
  size_t index = 0;
  for(std::string name : columns_names_[ntuple_id]){
    if(name==column_name){
      return columns_ids_[ntuple_id][index];
    }
    ++index;
  }
  return -999;
}

#endif
