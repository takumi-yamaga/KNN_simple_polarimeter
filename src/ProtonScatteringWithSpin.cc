// ProtonScatteringWithSpin.cc

#include "ProtonScatteringWithSpin.hh"

ProtonScatteringWithSpin::ProtonScatteringWithSpin()
{
  vec_proton_momentum_initial_ = TVector3(0.,0.,0.);
  vec_proton_momentum_final_ = TVector3(0.,0.,0.);
  vec_proton_spin_direction_ = TVector3(0.,0.,0.);

  random_ = new TRandom3(0);
}

ProtonScatteringWithSpin::~ProtonScatteringWithSpin()
{
  delete random_;
}

void ProtonScatteringWithSpin::Generate(TVector3& vec_proton_momentum_initial, TVector3& vec_proton_momentum_final_pre, TVector3& vec_proton_spin_direction, Double_t theta_scattering)
{
  // input values
  theta_scattering_ = theta_scattering;
  //theta_scattering_ = 10./180.*TMath::Pi();
  vec_proton_momentum_initial_ = vec_proton_momentum_initial;
  vec_proton_spin_direction_ = vec_proton_spin_direction;

  if(vec_proton_spin_direction_.Mag()==0.){
    Double_t theta_spin = random_->Uniform(0.,TMath::Pi());
    Double_t phi_spin = random_->Uniform(-TMath::Pi(),TMath::Pi());
    vec_proton_spin_direction_ = TVector3(sin(theta_spin)*cos(phi_spin),sin(theta_spin)*sin(phi_spin),cos(theta_spin));
  }

  // effective spin component (transversal component)
  Double_t effective_spin = vec_proton_spin_direction_.Cross(vec_proton_momentum_initial_.Unit()).Mag();

  // x, y unit vectors
  vec_x_ = vec_proton_spin_direction_.Cross(vec_proton_momentum_initial_.Unit()).Unit();
  vec_y_ = vec_proton_momentum_initial_.Unit().Cross(vec_x_).Unit();

  // analyzing power
  analyzing_power_ = val_A(vec_proton_momentum_initial_.Mag(),theta_scattering_);
  //analyzing_power_ = 0.4; // set to be constant for the present.

  // scattering angle phi
  Double_t phi_scattering = -999.;
  while(true){
    phi_scattering = random_->Uniform(-TMath::Pi(),TMath::Pi());
    Double_t val_ratio = 1. + analyzing_power_ * effective_spin * cos(phi_scattering);
    Double_t val_random = random_->Uniform(0.,1.+analyzing_power_*effective_spin);
    if(val_ratio>val_random){
      break;
    }
  }
  vec_proton_momentum_final_ = vec_proton_momentum_initial_;
  vec_proton_momentum_final_.Rotate(theta_scattering_,vec_y_);
  vec_proton_momentum_final_.Rotate(phi_scattering,vec_proton_momentum_initial_.Unit());
  vec_proton_momentum_final_.SetMag(vec_proton_momentum_final_pre.Mag());

  return;

}

void ProtonScatteringWithSpin::Print()
{
  std::cout << "ProtonScatteringWithSpin::Print()------------------ " << std::endl;
  std::cout << "analysing_power : " << analyzing_power_ << std::endl;
  std::cout << std::scientific << std::setprecision(2) 
    << "vec_proton_momentum_initial : " << vec_proton_momentum_initial_[0] << ", " << vec_proton_momentum_initial_[1] << ", " << vec_proton_momentum_initial_[2] << std::endl;
  std::cout << std::scientific << std::setprecision(2) 
    << "vec_proton_momentum_final : " << vec_proton_momentum_final_[0] << ", " << vec_proton_momentum_final_[1] << ", " << vec_proton_momentum_final_[2] << std::endl;
  std::cout << std::scientific << std::setprecision(2) 
    << "vec_proton_spin_direction : " << vec_proton_spin_direction_[0] << ", " << vec_proton_spin_direction_[1] << ", " << vec_proton_spin_direction_[2] << std::endl;
  std::cout << std::scientific << std::setprecision(2) 
    << "scattering angle (theta) : " << TMath::ACos(vec_proton_momentum_initial_.Unit().Dot(vec_proton_momentum_final_.Unit()))/TMath::Pi()*180. << std::endl;
}
