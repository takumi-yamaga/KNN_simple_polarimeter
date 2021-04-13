/// \file include/StackingAction.hh
/// \brief Definition of the StackingAction class

#ifndef StackingAction_h
#define StackingAction_h 1

#include "globals.hh"
#include "G4UserStackingAction.hh"

class G4Track;
class G4VHitsCollection;

class StackingAction : public G4UserStackingAction
{
public:
  StackingAction();
  virtual ~StackingAction();

public:
  virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* track);
  virtual void NewStage();
  virtual void PrepareNewEvent();
};

#endif
