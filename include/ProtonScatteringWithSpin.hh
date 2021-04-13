// ProtonScatteringWithSpin.hh

#ifndef ProtonScatteringWithSpin_hh
#define ProtonScatteringWithSpin_hh

#include <iostream>
#include <iomanip>

#include "TRandom3.h"
#include "TVector3.h"

class ProtonScatteringWithSpin
{
  public:
    ProtonScatteringWithSpin();
    ~ProtonScatteringWithSpin();

  public:
    void Generate(TVector3& vec_proton_momentum_initial, TVector3& vec_proton_momentum_final_pre, TVector3& vec_proton_spin_direction, Double_t theta_scattering);
    void Print();

    TVector3 VecProtonMomentumInitial() { return vec_proton_momentum_initial_; }
    TVector3 VecProtonMomentumFinal() { return vec_proton_momentum_final_; }
    TVector3 VecProtonSpin() { return vec_proton_spin_direction_; }
    TVector3 VecX(){ return vec_x_; }
    TVector3 VecY(){ return vec_y_; }
    Double_t AnalyzingPower() { return analyzing_power_; }

  private:
    TVector3 vec_proton_momentum_initial_;
    TVector3 vec_proton_momentum_final_;
    TVector3 vec_proton_spin_direction_;
    TVector3 vec_x_;
    TVector3 vec_y_;

    Double_t theta_scattering_;
    Double_t analyzing_power_;

  private:
    inline Double_t val_A(Double_t, Double_t);
    inline Double_t val_pol(Double_t, Double_t*);

    TRandom3* random_;

};

inline Double_t ProtonScatteringWithSpin::val_pol(Double_t p, Double_t* par){
  // function for alpha, beta, gamma, and delta

  // p : Momentum in GeV/c

  // par[0] : polynomial_parameter_0
  // par[1] : polynomial_parameter_1
  // par[2] : polynomial_parameter_2
  // par[3] : polynomial_parameter_3
  // par[4] : polynomial_parameter_4
  // par[5] : p_central

  Double_t X = p - par[5];
  //std::cout << p << " : " << X << std::endl;

  Double_t val = 0.;
  val += par[0];
  val += pow(X,1) * par[1];
  val += pow(X,2) * par[2];
  val += pow(X,3) * par[3];
  val += pow(X,4) * par[4];

  return val;
}

inline Double_t ProtonScatteringWithSpin::val_A(Double_t p, Double_t theta){
  // inputs
  // p : mometum (MeV/c)
  // theta : scattering angle (rad.)
  p /= 1000.; // MeV/c -> GeV/c

  //parameters from NIM201 Low energy fitting
  Double_t par[21];

  // par[0] : p_central
  par[0] = 0.7;

  // for alpha(p)
  // par[1] : alpha_0
  // par[2] : alpha_1
  // par[3] : alpha_2
  // par[4] : alpha_3
  // par[5] : alpha_4
  par[1] = 5.4771;
  par[2] = -4.2906;
  par[3] = -25.379;
  par[4] = 121.15;
  par[5] = -108.37;

  // for beta(p)
  // par[6]  : beta_0
  // par[7]  : beta_1
  // par[8]  : beta_2
  // par[9]  : beta_3
  // par[10] : beta_4
  par[6]  = -10.475;
  par[7]  = -40.170;
  par[8]  = 525.84;
  par[9]  = -899.29;
  par[10] = 1616.6;

  // for gamma(p)
  // par[11] : gamma_0
  // par[12] : gamma_1
  // par[13] : gamma_2
  // par[14] : gamma_3
  // par[15] : gamma_4
  par[11] = 1052.3;
  par[12] = 628.51;
  par[13] = -13215.0;
  par[14] = 19083.0;
  par[15] = -5485.5;

  // for delta(p)
  // par[16] : delta_0
  // par[17] : delta_1
  // par[18] : delta_2
  // par[19] : delta_3
  // par[20] : delta_4
  par[16] = 0.;
  par[17] = 0.; 
  par[18] = 0.; 
  par[19] = 0.; 
  par[20] = 0.; 

  // r
  Double_t r = p * sin(theta);

  // alpha(T)
  Double_t par_alpha[6] = {par[1],par[2],par[3],par[4],par[5],par[0]};
  Double_t val_alpha = val_pol(p,par_alpha);

  // beta(T)
  Double_t par_beta[6] = {par[6],par[7],par[8],par[9],par[10],par[0]};
  Double_t val_beta = val_pol(p,par_beta);

  // gamma(T)
  Double_t par_gamma[6] = {par[11],par[12],par[13],par[14],par[15],par[0]};
  Double_t val_gamma = val_pol(p,par_gamma);

  // delta(T)
  Double_t par_delta[6] = {par[16],par[17],par[18],par[19],par[20],par[0]};
  Double_t val_delta = val_pol(p,par_delta);

  // Ac(theta,T)
  Double_t val = val_alpha*r / (1. + val_beta*r*r + val_gamma*r*r*r*r);
  val += val_delta*p*sin(5.*theta);

  //std::cout << val << std::endl;

  if(val<0.){ 
    val = 0.;
  }

  return val;
  
}

#endif
