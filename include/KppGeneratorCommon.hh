// KppGeneratorCommon.hh

#ifndef KppGeneratorCommon_hh
#define KppGeneratorCommon_hh

#include "G4SystemOfUnits.hh"

namespace KppGeneratorCommon{

const Double_t kPiMass = 0.13957*GeV;
const Double_t kPi0Mass= 0.1349766*GeV;
const Double_t kProtonMass= 0.938272*GeV;
const Double_t kNeutronMass = 0.939565*GeV;
const Double_t kDeutronMass = 1.87561*GeV;
const Double_t kLambdaMass = 1.115683*GeV;
const Double_t kSigma0Mass= 1.192642*GeV;
const Double_t kSigmaPlusMass= 1.18937*GeV;
const Double_t kSigmaMinusMass= 1.197449*GeV;
const Double_t kKaonMass = 0.4936*GeV;
const Double_t kKaon0Mass = 0.497614*GeV;
const Double_t kThreeHeMass = 2.80839*GeV;
const Double_t kLightVelocity = 29.97; // cm/ns

const Double_t kMaxWeightTowBody = 1.0;
const Double_t kMaxWeightThreeBody = 0.5;
const Double_t kMaxRatioDecay = 3.5;
}

#endif
