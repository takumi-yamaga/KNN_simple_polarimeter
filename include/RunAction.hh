/// \brief Definition of the RunAction class

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

#include <vector>

class G4Run;

/// Run action class

class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    virtual ~RunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

  private:
    std::vector<int> ntuple_ids_; 
    std::vector<std::string> ntuple_names_; 

    std::vector< std::vector<int> > columns_ids_;
    std::vector< std::vector<std::string> > columns_names_; 

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
