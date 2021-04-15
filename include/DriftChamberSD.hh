/// \brief Definition of the DriftChamberSD class

#ifndef DriftChamberSD_h
#define DriftChamberSD_h 1

#include "G4VSensitiveDetector.hh"

#include "DriftChamberHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

/// Drift chamber sensitive detector

class DriftChamberSD : public G4VSensitiveDetector
{
  public:
    DriftChamberSD(G4String name);
    virtual ~DriftChamberSD();
    
    virtual void Initialize(G4HCofThisEvent*HCE);
    virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
    
  private:
    DriftChamberHitsCollection* fHitsCollection;
    G4int fHCID;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
