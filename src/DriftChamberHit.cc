/// \brief Implementation of the DriftChamberHit class

#include "DriftChamberHit.hh"
#include "Constants.hh"

#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4Allocator<DriftChamberHit>* DriftChamberHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DriftChamberHit::DriftChamberHit()
: G4VHit(), 
  track_id_(-1), parent_id_(-1), particle_id_(-1), layer_id_(-1), hit_time_(0.), local_position_(0), global_position_(0), momentum_(0), polarization_(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DriftChamberHit::DriftChamberHit(G4int layer_id)
: G4VHit(), 
  track_id_(layer_id),parent_id_(layer_id),particle_id_(layer_id),layer_id_(layer_id), hit_time_(0.), local_position_(0), global_position_(0), momentum_(0), polarization_(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DriftChamberHit::~DriftChamberHit()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DriftChamberHit::DriftChamberHit(const DriftChamberHit &right)
: G4VHit(),
  track_id_(right.track_id_),
  parent_id_(right.parent_id_),
  particle_id_(right.particle_id_),
  layer_id_(right.layer_id_),
  hit_time_(right.hit_time_),
  local_position_(right.local_position_),
  global_position_(right.global_position_),
  momentum_(right.momentum_),
  polarization_(right.polarization_)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const DriftChamberHit& DriftChamberHit::operator=(const DriftChamberHit &right)
{
  track_id_ = right.track_id_;
  parent_id_ = right.parent_id_;
  particle_id_ = right.particle_id_;
  layer_id_ = right.layer_id_;
  hit_time_ = right.hit_time_;
  local_position_ = right.local_position_;
  global_position_ = right.global_position_;
  momentum_ = right.momentum_;
  polarization_ = right.polarization_;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool DriftChamberHit::operator==(const DriftChamberHit &/*right*/) const
{
  return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DriftChamberHit::Draw()
{
  auto visManager = G4VVisManager::GetConcreteInstance();
  if (! visManager) return;

  G4Circle circle(global_position_);
  circle.SetScreenSize(10);
  circle.SetFillStyle(G4Circle::filled);
  G4VisAttributes attribs(Colors::DriftchamberHit());
  circle.SetVisAttributes(attribs);
  visManager->Draw(circle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const std::map<G4String,G4AttDef>* DriftChamberHit::GetAttDefs() const
{
  G4bool isNew;
  auto store = G4AttDefStore::GetInstance("DriftChamberHit",isNew);

  if (isNew) {
      (*store)["HitType"] 
        = G4AttDef("HitType","Hit Type","Physics","","G4String");
      
      (*store)["ID"] 
        = G4AttDef("ID","ID","Physics","","G4int");
      
      (*store)["Time"] 
        = G4AttDef("Time","Time","Physics","G4BestUnit","G4double");
      
      (*store)["Pos"] 
        = G4AttDef("Pos", "Position", "Physics","G4BestUnit","G4ThreeVector");
  }
  
  return store;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<G4AttValue>* DriftChamberHit::CreateAttValues() const
{
  auto values = new std::vector<G4AttValue>;
  
  values
    ->push_back(G4AttValue("HitType","DriftChamberHit",""));
  values
    ->push_back(G4AttValue("ID",G4UIcommand::ConvertToString(layer_id_),""));
  values
    ->push_back(G4AttValue("Time",G4BestUnit(global_position_,"Time"),""));
  values
    ->push_back(G4AttValue("Position",G4BestUnit(global_position_,"Length"),""));
  
  return values;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DriftChamberHit::Print()
{
  //G4cout << "  Layer[" << fLayerID << "] : time " << fTime/ns
  //<< " (nsec) --- local (x,y) " << fLocalPos.x()
  //<< ", " << fLocalPos.y() << G4endl;
  G4cout << "polarization : (" << polarization_.x() << ", " << polarization_.y() << ", " << polarization_.z() << ")" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
