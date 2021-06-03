/// \brief Definition of constants.

#ifndef Constants_h
#define Constants_h 1

#include <array>
#include "globals.hh"
#include "G4Colour.hh"
#include "G4ThreeVector.hh"

#include <TVector3.h>

using std::array;

// driftchambers
namespace Driftchamber{
  constexpr G4int kTotalNumber = 11;
  const array<G4String, kTotalNumber> kDetectorNames
      = {{ "cdc",
          "tracker_layer1", 
          "tracker_layer2",
          "tracker_layer3",
          "tracker_layer4",
          "tracker_layer5",
          "tracker_layer6",
          "tracker_layer7",
          "tracker_layer8",
          "tracker_layer9",
          "tracker_layer10"
      }};
  const G4int kCDCId = 0;
  const G4int kTrackerLayer1Id  = 1;
  const G4int kTrackerLayer2Id  = 2;
  const G4int kTrackerLayer3Id  = 3;
  const G4int kTrackerLayer4Id  = 4;
  const G4int kTrackerLayer5Id  = 5;
  const G4int kTrackerLayer6Id  = 6;
  const G4int kTrackerLayer7Id  = 7;
  const G4int kTrackerLayer8Id  = 8;
  const G4int kTrackerLayer9Id  = 9;
  const G4int kTrackerLayer10Id = 10;
}

// colors
namespace Colors{
    // G4Colour(red, green, blue, alpha)
    // alpha = 1. - transparency

    inline G4Colour Scintillator(){ return G4Colour(0.2,1.0,1.0,0.2); }
    inline G4Colour ScintillatorHasHit(){ return G4Colour(1.0,0.0,0.0,0.2); }
    inline G4Colour DriftchamberHit(){ return G4Colour(1.0,0.0,0.0,0.2); }
    inline G4Colour Reflector(){ return G4Colour(0.7882,0.7922,0.7922,0.1); }
    inline G4Colour LightShield(){ return G4Colour(0.1,0.1,0.1,0.1); }
    inline G4Colour LightGuide(){ return G4Colour(0.0,0.0,1.0,0.5); }
    inline G4Colour Transparent(){ return G4Colour(0.0,0.0,0.0,0.0); }
    inline G4Colour PMTWindow(){ return G4Colour(139./255.,69./255.,19./255.,0.8); }
    inline G4Colour PMTInner(){ return G4Colour(0.05,0.05,0.9,0.1); }
    inline G4Colour PMTOuter(){ return G4Colour(55./255.,55./255.,55./255.,1.0); }
    inline G4Colour Magnetic(){ return G4Colour(0.3,0.3,0.3,0.2); }
    inline G4Colour WirePlane(){ return G4Colour(240./255.,145./255.,153./255.,0.2); }

    inline G4Colour Hit(){ return G4Colour(1.0,0.0,0.0,1.0); }
}

// utility functions

inline G4ThreeVector ConvertVector(TVector3 vector){
    return G4ThreeVector(vector.x(),vector.y(),vector.z());
}
inline TVector3 ConvertVector(G4ThreeVector vector){
    return TVector3(vector.x(),vector.y(),vector.z());
}
#endif
