// SteppingSction.hh

#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "globals.hh"

class G4Step;
class LambdaDecay;
class ProtonScatteringWithSpin;

class SteppingAction: public G4UserSteppingAction
{
  public:
    SteppingAction();
    ~SteppingAction();

    void UserSteppingAction(const G4Step* aStep);

  private:
    LambdaDecay* lambda_decay_;
    ProtonScatteringWithSpin* proton_scattering_with_spin_;

};
#endif
