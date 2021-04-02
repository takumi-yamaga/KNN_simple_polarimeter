// KppGenerator.hh

#ifndef KppGenerator_hh
#define KppGenerator_hh

#include <iostream>
#include <iomanip>

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGenPhaseSpace.h"
#include "TRandom3.h"
#include "TF1.h"

#include "KppGeneratorCommon.hh"

class KppGenerator
{
  public:
    KppGenerator();
    KppGenerator(Int_t kpp_spin, Int_t kpp_parity);
    ~KppGenerator();

  public:
    void Generate(TLorentzVector& lv_beam);
    void Print();

    TLorentzVector LvKpp(){ return lv_kpp_; }
    TLorentzVector LvLambda(){ return lv_lambda_; }
    TLorentzVector LvProton(){ return lv_proton_; }
    TLorentzVector LvNeutron(){ return lv_neutron_; }

    Int_t KppSpin(){ return kpp_spin_; }
    void SetKppSpin(Int_t val){ kpp_spin_ = val; }

    Int_t KppParity(){ return kpp_parity_; }
    void SetKppParity(Int_t val){ kpp_parity_ = val; }

    Double_t KppMass(){ return kpp_mass_; }
    inline void SetKppMass(Double_t val){
      kpp_mass_ = val; 
      if(f_kpp_mass_){
        f_kpp_mass_->SetParameter(0,kpp_mass_);
      }
    }

    Double_t KppWidth(){ return kpp_width_; }
    inline void SetKppWidth(Double_t val){ 
      kpp_width_ = val;
      if(f_kpp_mass_){
        f_kpp_mass_->SetParameter(1,kpp_width_);
      }
    }

    Double_t KppFormFactor(){ return kpp_formfactor_; }
    inline void SetKppFormFactor(Double_t val){
      kpp_formfactor_ = val; 
      if(f_kpp_momtrans_){
        f_kpp_momtrans_->SetParameter(0,kpp_formfactor_);
      }
    }

    Int_t MLDecay(){ return ml_decay_; }
    Int_t MLSpin(){ return ml_spin_; }
    TVector3 VecReference(){ return vec_reference_direction_; }
    TVector3 VecLambdaSpin(){ return vec_lambda_spin_direction_; }
    TVector3 VecProtonSpin(){ return vec_proton_spin_direction_; }

  private:
    void KppZeroMinusDecay();
    void KppOneMinusDecay();
    void KppZeroPlusDecay();
    void KppOnePlusDecay();
    TRandom3* random_;

    TF1* f_kpp_mass_;
    TF1* f_kpp_momtrans_;
    
    TLorentzVector lv_target_;
    TGenPhaseSpace generator_;
    Double_t cm_masses_[3] = { KppGeneratorCommon::kLambdaMass, KppGeneratorCommon::kProtonMass, KppGeneratorCommon::kNeutronMass};

    TLorentzVector lv_kpp_;
    TLorentzVector lv_lambda_;
    TLorentzVector lv_proton_;
    TLorentzVector lv_neutron_;

    Int_t kpp_spin_;
    Int_t kpp_parity_;
    Int_t ml_decay_;
    Int_t ml_spin_;
    Double_t kpp_mass_;
    Double_t kpp_width_;
    Double_t kpp_formfactor_;
    TVector3 vec_reference_direction_;
    TVector3 vec_lambda_spin_direction_;
    TVector3 vec_proton_spin_direction_;
};

#endif
