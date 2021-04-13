/// \file include/TrackInformation.hh
/// \brief Definition of the TrackInformation class

#ifndef TrackInformation_h
#define TrackInformation_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Allocator.hh"
#include "G4VUserTrackInformation.hh"

#include "vector"

class TrackInformation : public G4VUserTrackInformation 
{
  public:
    TrackInformation();
    TrackInformation(const G4Track* track);
    TrackInformation(const TrackInformation* track_information);
    virtual ~TrackInformation();

    inline void *operator new(size_t);
    inline void operator delete(void *track_information);

    TrackInformation& operator =(const TrackInformation& right);

    virtual void Print() const;

    // setter -----------------------------------------------------------------
    inline void SetCurrentSpinDirection(const G4ThreeVector vect){ current_spin_direction_ = vect; }
    inline void PushTrackInformation(const G4Track*);
    inline void PushTrackInformation(const TrackInformation*);
    // getter -----------------------------------------------------------------
    inline G4int GetTrackID() const {return track_id_;}
    inline G4int GetParentID() const {return parent_id_;}
    inline G4String GetParticleName() const {return particle_name_;}
    inline G4ThreeVector GetInitialMomentum() const {return initial_momentum_;}
    inline G4ThreeVector GetInitialSpinDirection() const {return initial_spin_direction_;}
    inline G4ThreeVector GetCurrentSpinDirection() const {return current_spin_direction_;}
    inline G4ThreeVector GetInitialPosition() const {return initial_position_;}
    inline std::vector <G4int> GetParentTrackIDs() const {return parent_track_ids_;}
    inline std::vector <G4int> GetParentParentIDs() const {return parent_parent_ids_;}
    inline std::vector <G4String> GetParentParticleNames() const {return parent_particle_names_;}
    inline std::vector <G4ThreeVector> GetParentInitialMomenta() const {return parent_initial_momenta_;}
    inline std::vector <G4ThreeVector> GetParentInitialPositions() const {return parent_initial_positions_;}
    // ------------------------------------------------------------------------

  private:
    G4int track_id_;
    G4int parent_id_;
    G4String particle_name_;
    G4ThreeVector initial_momentum_;
    G4ThreeVector initial_position_;
    G4ThreeVector initial_spin_direction_;
    G4ThreeVector current_spin_direction_;
    std::vector <G4int> parent_track_ids_;
    std::vector <G4int> parent_parent_ids_;
    std::vector <G4String> parent_particle_names_;
    std::vector <G4ThreeVector> parent_initial_momenta_;
    std::vector <G4ThreeVector> parent_initial_positions_;
    std::vector <G4ThreeVector> parent_initial_spin_directions_;
    std::vector <G4ThreeVector> parent_current_spin_directions_;
};

extern G4ThreadLocal
G4Allocator<TrackInformation> * aTrackInformationAllocator;

inline void* TrackInformation::operator new(size_t)
{
  if(!aTrackInformationAllocator)
    aTrackInformationAllocator = new G4Allocator<TrackInformation>;
  return (void*)aTrackInformationAllocator->MallocSingle();
}

inline void TrackInformation::operator delete(void *track_information)
{ aTrackInformationAllocator->FreeSingle((TrackInformation*)track_information);}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline void TrackInformation::PushTrackInformation(const G4Track* track)
{
  // push current track_information as the newest parent_information
  parent_track_ids_.push_back(this->track_id_);
  parent_parent_ids_.push_back(this->parent_id_);
  parent_particle_names_.push_back(this->particle_name_);
  parent_initial_momenta_.push_back(this->initial_momentum_);
  parent_initial_positions_.push_back(this->initial_position_);
  parent_initial_spin_directions_.push_back(this->initial_spin_direction_);
  parent_current_spin_directions_.push_back(this->current_spin_direction_);

  // set current track_information
  track_id_ = track->GetTrackID();
  parent_id_ = track->GetParentID();
  particle_name_ = track->GetParticleDefinition()->GetParticleName();
  initial_momentum_ = track->GetMomentum();
  initial_position_ = track->GetPosition();
  initial_spin_direction_ = track->GetPolarization();
  current_spin_direction_ = initial_spin_direction_;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline void TrackInformation::PushTrackInformation(const TrackInformation* track_information)
{
  // push current track_information as the newest parent_information
  parent_track_ids_.push_back(this->track_id_);
  parent_parent_ids_.push_back(this->parent_id_);
  parent_particle_names_.push_back(this->particle_name_);
  parent_initial_momenta_.push_back(this->initial_momentum_);
  parent_initial_positions_.push_back(this->initial_position_);
  parent_initial_spin_directions_.push_back(this->initial_spin_direction_);
  parent_current_spin_directions_.push_back(this->current_spin_direction_);

  // set current track_information
  track_id_ = track_information->track_id_;
  parent_id_ = track_information->parent_id_;
  particle_name_ = track_information->particle_name_;
  initial_momentum_ = track_information->initial_momentum_;
  initial_position_ = track_information->initial_position_;
  initial_spin_direction_ = track_information->initial_spin_direction_;
}

#endif

