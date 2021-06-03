// analysis_functions.cc

#include"analysis_functions.hh"
#include"resolution_functions.hh"
#include"DrawTools.h"
#include"MyConstants.h"

#include"TF1.h"
#include"TLorentzVector.h"
#include"TLatex.h"

// -----
TTree* tree_event;
Int_t event_id;
Int_t number_of_primaries;
Int_t number_of_trajectories;
Int_t number_of_hits_in_cdc;
Int_t number_of_hits_in_tracker[10];
// -----
TTree* tree_primary;
Int_t primary_event_id;
Int_t primary_particle_id;
Float_t primary_position_x;
Float_t primary_position_y;
Float_t primary_position_z;
Float_t primary_momentum_x;
Float_t primary_momentum_y;
Float_t primary_momentum_z;
Float_t primary_spin_x;
Float_t primary_spin_y;
Float_t primary_spin_z;
// -----
TTree* tree_trajectory;
Int_t trajectory_event_id;
Int_t trajectory_track_id;
Int_t trajectory_parent_id;
Int_t trajectory_particle_id;
Float_t trajectory_initial_position_x;
Float_t trajectory_initial_position_y;
Float_t trajectory_initial_position_z;
Float_t trajectory_initial_momentum_x;
Float_t trajectory_initial_momentum_y;
Float_t trajectory_initial_momentum_z;
// -----
TTree* tree_cdc_hit;
Int_t cdc_event_id;
Int_t cdc_track_id;
Int_t cdc_parent_id;
Int_t cdc_particle_id;
Int_t cdc_layer_id;
Float_t cdc_hit_time;
Float_t cdc_hit_position_x;
Float_t cdc_hit_position_y;
Float_t cdc_hit_position_z;
Float_t cdc_momentum_x;
Float_t cdc_momentum_y;
Float_t cdc_momentum_z;
Int_t cdc_is_asymmetric_scattering;
// -----
TTree* tree_tracker_hit[10];
Int_t tracker_event_id[10];
Int_t tracker_track_id[10];
Int_t tracker_parent_id[10];
Int_t tracker_particle_id[10];
Int_t tracker_layer_id[10];
Float_t tracker_hit_time[10];
Float_t tracker_hit_position_x[10];
Float_t tracker_hit_position_y[10];
Float_t tracker_hit_position_z[10];
Float_t tracker_momentum_x[10];
Float_t tracker_momentum_y[10];
Float_t tracker_momentum_z[10];
Int_t tracker_is_asymmetric_scattering[10];
// -----

void DrawHistograms(TFile* outfile, std::string pdf_name){
  return;
  TCanvas* c_title = new TCanvas("c_title","c_title",1500,1000);
  c_title->cd();
  c_title->Divide(3,2);

  TCanvas* canvas = new TCanvas("canvas","canvas",1500,1000);
  canvas->cd();
  canvas->Divide(3,2);

  // open pdf file
  canvas->Print(std::string(pdf_name + "[").data());

  TH1F* hist_1d = nullptr;
  TH1F* hist_1d_clone = nullptr;
  TH2F* hist_2d = nullptr;

  TLatex* text = new TLatex();
  text->SetTextFont(132);
  text->SetTextAlign(22);
  text->SetTextSize(0.1);

  TString sel_name[12] = {
    "",
    "_charged",
    "_charged_mom",
    "_charged_mom_cdc",
    "_charged_mom_cdc_proton_both",
    "_charged_mom_cdc_proton_both_asym",
    "_charged_mom_cdc_proton_both_sel",
    "_charged_mom_cdc_proton_both_sel_asym",
    "_charged_mom_cdc_proton_from_lambda_both",
    "_charged_mom_cdc_proton_from_lambda_both_asym",
    "_charged_mom_cdc_proton_from_lambda_both_sel",
    "_charged_mom_cdc_proton_from_lambda_both_sel_asym"
  };

  for(int i_sel=0; i_sel<12; ++i_sel){
    std::cout << sel_name[i_sel].Data() << std::endl;
    c_title->Clear();
    c_title->cd();
    text->DrawLatex(0.5,0.5,Form("%s",sel_name[i_sel].Data()));
    c_title->Print(pdf_name.data());

    // page-1
    // lambda_momentum
    canvas->cd(1);
    hist_1d = (TH1F*)outfile->Get(Form("lambda_momentum%s",sel_name[i_sel].Data()));
    Draw(hist_1d);

    // proton_momentum
    canvas->cd(2);
    hist_1d = (TH1F*)outfile->Get(Form("proton_momentum%s",sel_name[i_sel].Data()));
    Draw(hist_1d);

    // proton_from_lambda_momentum
    canvas->cd(4);
    hist_1d = (TH1F*)outfile->Get(Form("proton_from_lambda_momentum%s",sel_name[i_sel].Data()));
    Draw(hist_1d);

    // pim_from_lambda_momentum
    canvas->cd(5);
    hist_1d = (TH1F*)outfile->Get(Form("pim_from_lambda_momentum%s",sel_name[i_sel].Data()));
    Draw(hist_1d);

    canvas->Print(pdf_name.data());
    for(int i_pad=1; i_pad<=6; i_pad++){
      canvas->cd(i_pad);
      gPad->SetLogx(0);
      gPad->SetLogy(0);
      gPad->SetLogz(0);
      gPad->Clear();
    }
    // -----

    // page-2
    // lp_mass
    canvas->cd(1);
    hist_1d = (TH1F*)outfile->Get(Form("lp_mass%s",sel_name[i_sel].Data()));
    Draw(hist_1d);

    // lp_momentum
    canvas->cd(2);
    hist_1d = (TH1F*)outfile->Get(Form("lp_momentum%s",sel_name[i_sel].Data()));
    Draw(hist_1d);

    // lp_mass_vs_lp_momentum
    canvas->cd(3);
    hist_2d = (TH2F*)outfile->Get(Form("lp_mass_vs_lp_momentum%s",sel_name[i_sel].Data()));
    Draw(hist_2d,"col");

    // ppim_mass
    canvas->cd(4);
    hist_1d = (TH1F*)outfile->Get(Form("ppim_mass%s",sel_name[i_sel].Data()));
    Draw(hist_1d);

    // ppim_momentum
    canvas->cd(5);
    hist_1d = (TH1F*)outfile->Get(Form("ppim_momentum%s",sel_name[i_sel].Data()));
    Draw(hist_1d);

    canvas->Print(pdf_name.data());
    for(int i_pad=1; i_pad<=6; i_pad++){
      canvas->cd(i_pad)->Clear();
    }
    // -----

    // page-3
    // proton_momentum_at_cdc
    canvas->cd(1);
    hist_1d = (TH1F*)outfile->Get(Form("proton_momentum_at_cdc%s",sel_name[i_sel].Data()));
    Draw(hist_1d);

    // proton_momentum_at_tracker_layer1
    canvas->cd(2);
    hist_1d = (TH1F*)outfile->Get(Form("proton_momentum_at_tracker_layer1%s",sel_name[i_sel].Data()));
    Draw(hist_1d);

    // proton_momentum_at_tracker_layer2
    canvas->cd(3);
    hist_1d = (TH1F*)outfile->Get(Form("proton_momentum_at_tracker_layer2%s",sel_name[i_sel].Data()));
    Draw(hist_1d);

    // proton_momentum_difference_between_cdc_and_tracker_layer1
    canvas->cd(4);
    hist_1d = (TH1F*)outfile->Get(Form("proton_momentum_difference_between_cdc_and_tracker_layer1%s",sel_name[i_sel].Data()));
    Draw(hist_1d);

    // scattering_angle_theta_mc
    canvas->cd(5);
    gPad->SetLogy();
    hist_1d = (TH1F*)outfile->Get(Form("proton_scattering_angle_theta_mc%s",sel_name[i_sel].Data()));
    Draw(hist_1d);

    canvas->Print(pdf_name.data());
    for(int i_pad=1; i_pad<=6; i_pad++){
      canvas->cd(i_pad);
      gPad->SetLogx(0);
      gPad->SetLogy(0);
      gPad->SetLogz(0);
      gPad->Clear();
    }
    // -----
    //
    TF1* f_phi = new TF1("f_phi","1 + [0]*cos(x)",-TMath::Pi(),TMath::Pi());

    // page-4
    // scattering_angle_theta
    canvas->cd(1);
    gPad->SetLogy();
    hist_1d = (TH1F*)outfile->Get(Form("proton_scattering_angle_theta%s",sel_name[i_sel].Data()));
    Draw(hist_1d);

    // phi_of_spins
    canvas->cd(2);
    hist_1d = (TH1F*)outfile->Get(Form("proton_phi_of_spins%s",sel_name[i_sel].Data()));
    Draw(hist_1d);
    hist_1d->SetMinimum(0);

    // phi_of_spins (scaled & fitting)
    canvas->cd(3);
    hist_1d_clone = new TH1F();
    hist_1d_clone = (TH1F*)hist_1d->Clone();
    Draw(hist_1d_clone,"e");
    hist_1d_clone->Scale(1./(hist_1d->Integral()/hist_1d->GetNbinsX()));
    hist_1d_clone->SetMinimum(0.8);
    hist_1d_clone->SetMaximum(1.2);
    if(hist_1d_clone->GetEntries()){
      hist_1d_clone->Fit("f_phi");
    }

    canvas->Print(pdf_name.data());
    for(int i_pad=1; i_pad<=6; i_pad++){
      canvas->cd(i_pad);
      gPad->SetLogx(0);
      gPad->SetLogy(0);
      gPad->SetLogz(0);
      gPad->Clear();
    }
    // -----

    // page-5
    // proton_from_lambda_momentum_at_cdc
    canvas->cd(1);
    hist_1d = (TH1F*)outfile->Get(Form("proton_from_lambda_momentum_at_cdc%s",sel_name[i_sel].Data()));
    Draw(hist_1d);

    // proton_from_lambda_momentum_at_tracker_layer1
    canvas->cd(2);
    hist_1d = (TH1F*)outfile->Get(Form("proton_from_lambda_momentum_at_tracker_layer1%s",sel_name[i_sel].Data()));
    Draw(hist_1d);

    // proton_from_lambda_momentum_at_tracker_layer2
    canvas->cd(3);
    hist_1d = (TH1F*)outfile->Get(Form("proton_from_lambda_momentum_at_tracker_layer2%s",sel_name[i_sel].Data()));
    Draw(hist_1d);

    // proton_from_lambda_momentum_difference_between_cdc_and_tracker_layer1
    canvas->cd(4);
    hist_1d = (TH1F*)outfile->Get(Form("proton_from_lambda_momentum_difference_between_cdc_and_tracker_layer1%s",sel_name[i_sel].Data()));
    Draw(hist_1d);

    // scattering_angle_theta_mc
    canvas->cd(5);
    gPad->SetLogy();
    hist_1d = (TH1F*)outfile->Get(Form("proton_from_lambda_scattering_angle_theta_mc%s",sel_name[i_sel].Data()));
    Draw(hist_1d);

    canvas->Print(pdf_name.data());
    for(int i_pad=1; i_pad<=6; i_pad++){
      canvas->cd(i_pad);
      gPad->SetLogx(0);
      gPad->SetLogy(0);
      gPad->SetLogz(0);
      gPad->Clear();
    }
    // -----

    // page-6
    // scattering_angle_theta
    canvas->cd(1);
    gPad->SetLogy();
    hist_1d = (TH1F*)outfile->Get(Form("proton_from_lambda_scattering_angle_theta%s",sel_name[i_sel].Data()));
    Draw(hist_1d);

    // phi_of_spins
    canvas->cd(2);
    hist_1d = (TH1F*)outfile->Get(Form("proton_from_lambda_phi_of_spins%s",sel_name[i_sel].Data()));
    Draw(hist_1d);
    hist_1d->SetMinimum(0);

    // phi_of_spins (scaled & fitting)
    canvas->cd(3);
    hist_1d_clone = new TH1F();
    hist_1d_clone = (TH1F*)hist_1d->Clone();
    Draw(hist_1d_clone,"e");
    hist_1d_clone->Scale(1./(hist_1d->Integral()/hist_1d->GetNbinsX()));
    hist_1d_clone->SetMinimum(0.8);
    hist_1d_clone->SetMaximum(1.2);
    if(hist_1d_clone->GetEntries()){
      hist_1d_clone->Fit("f_phi");
    }

    canvas->Print(pdf_name.data());
    for(int i_pad=1; i_pad<=6; i_pad++){
      canvas->cd(i_pad);
      gPad->SetLogx(0);
      gPad->SetLogy(0);
      gPad->SetLogz(0);
      gPad->Clear();
    }
    // -----

  }

  // close pdf file
  canvas->Print(std::string(pdf_name + "]").data());
}

void CreateHistograms(TFile* outfile){
  outfile->cd();

  TString sel_name[84] = {
    /*  0 */ "",
    /*  1 */ "_charged",
    /*  2 */ "_charged_mom",
    /*  3 */ "_charged_mom_cdc",

    /*  4 */ "_charged_mom_cdc_proton_pair1",
    /*  5 */ "_charged_mom_cdc_proton_pair1_asym",
    /*  6 */ "_charged_mom_cdc_proton_pair1_sel",
    /*  7 */ "_charged_mom_cdc_proton_pair1_sel_asym",
    /*  8 */ "_charged_mom_cdc_proton_pair1_sel_mqsel",
    /*  9 */ "_charged_mom_cdc_proton_pair1_selmeas",
    /* 10 */ "_charged_mom_cdc_proton_pair1_selmeas_asym",
    /* 11 */ "_charged_mom_cdc_proton_pair1_selmeas_mqsel",
    /* 12 */ "_charged_mom_cdc_proton_pair2",
    /* 13 */ "_charged_mom_cdc_proton_pair2_asym",
    /* 14 */ "_charged_mom_cdc_proton_pair2_sel",
    /* 15 */ "_charged_mom_cdc_proton_pair2_sel_asym",
    /* 16 */ "_charged_mom_cdc_proton_pair2_sel_mqsel",
    /* 17 */ "_charged_mom_cdc_proton_pair2_selmeas",
    /* 18 */ "_charged_mom_cdc_proton_pair2_selmeas_asym",
    /* 19 */ "_charged_mom_cdc_proton_pair2_selmeas_mqsel",
    /* 20 */ "_charged_mom_cdc_proton_pair3",
    /* 21 */ "_charged_mom_cdc_proton_pair3_asym",
    /* 22 */ "_charged_mom_cdc_proton_pair3_sel",
    /* 23 */ "_charged_mom_cdc_proton_pair3_sel_asym",
    /* 24 */ "_charged_mom_cdc_proton_pair3_sel_mqsel",
    /* 25 */ "_charged_mom_cdc_proton_pair3_selmeas",
    /* 26 */ "_charged_mom_cdc_proton_pair3_selmeas_asym",
    /* 27 */ "_charged_mom_cdc_proton_pair3_selmeas_mqsel",
    /* 28 */ "_charged_mom_cdc_proton_pair4",
    /* 29 */ "_charged_mom_cdc_proton_pair4_asym",
    /* 30 */ "_charged_mom_cdc_proton_pair4_sel",
    /* 31 */ "_charged_mom_cdc_proton_pair4_sel_asym",
    /* 32 */ "_charged_mom_cdc_proton_pair4_sel_mqsel",
    /* 33 */ "_charged_mom_cdc_proton_pair4_selmeas",
    /* 34 */ "_charged_mom_cdc_proton_pair4_selmeas_asym",
    /* 35 */ "_charged_mom_cdc_proton_pair4_selmeas_mqsel",
    /* 36 */ "_charged_mom_cdc_proton_pair5",
    /* 37 */ "_charged_mom_cdc_proton_pair5_asym",
    /* 38 */ "_charged_mom_cdc_proton_pair5_sel",
    /* 39 */ "_charged_mom_cdc_proton_pair5_sel_asym",
    /* 40 */ "_charged_mom_cdc_proton_pair5_sel_mqsel",
    /* 41 */ "_charged_mom_cdc_proton_pair5_selmeas",
    /* 42 */ "_charged_mom_cdc_proton_pair5_selmeas_asym",
    /* 43 */ "_charged_mom_cdc_proton_pair5_selmeas_mqsel",

    /* 44 */ "_charged_mom_cdc_proton_from_lambda_pair1",
    /* 45 */ "_charged_mom_cdc_proton_from_lambda_pair1_asym",
    /* 46 */ "_charged_mom_cdc_proton_from_lambda_pair1_sel",
    /* 47 */ "_charged_mom_cdc_proton_from_lambda_pair1_sel_asym",
    /* 48 */ "_charged_mom_cdc_proton_from_lambda_pair1_sel_mqsel",
    /* 49 */ "_charged_mom_cdc_proton_from_lambda_pair1_selmeas",
    /* 50 */ "_charged_mom_cdc_proton_from_lambda_pair1_selmeas_asym",
    /* 51 */ "_charged_mom_cdc_proton_from_lambda_pair1_selmeas_mqsel",
    /* 52 */ "_charged_mom_cdc_proton_from_lambda_pair2",
    /* 53 */ "_charged_mom_cdc_proton_from_lambda_pair2_asym",
    /* 54 */ "_charged_mom_cdc_proton_from_lambda_pair2_sel",
    /* 55 */ "_charged_mom_cdc_proton_from_lambda_pair2_sel_asym",
    /* 56 */ "_charged_mom_cdc_proton_from_lambda_pair2_sel_mqsel",
    /* 57 */ "_charged_mom_cdc_proton_from_lambda_pair2_selmeas",
    /* 58 */ "_charged_mom_cdc_proton_from_lambda_pair2_selmeas_asym",
    /* 59 */ "_charged_mom_cdc_proton_from_lambda_pair2_selmeas_mqsel",
    /* 60 */ "_charged_mom_cdc_proton_from_lambda_pair3",
    /* 61 */ "_charged_mom_cdc_proton_from_lambda_pair3_asym",
    /* 62 */ "_charged_mom_cdc_proton_from_lambda_pair3_sel",
    /* 63 */ "_charged_mom_cdc_proton_from_lambda_pair3_sel_asym",
    /* 64 */ "_charged_mom_cdc_proton_from_lambda_pair3_sel_mqsel",
    /* 65 */ "_charged_mom_cdc_proton_from_lambda_pair3_selmeas",
    /* 66 */ "_charged_mom_cdc_proton_from_lambda_pair3_selmeas_asym",
    /* 67 */ "_charged_mom_cdc_proton_from_lambda_pair3_selmeas_mqsel",
    /* 68 */ "_charged_mom_cdc_proton_from_lambda_pair4",
    /* 69 */ "_charged_mom_cdc_proton_from_lambda_pair4_asym",
    /* 70 */ "_charged_mom_cdc_proton_from_lambda_pair4_sel",
    /* 71 */ "_charged_mom_cdc_proton_from_lambda_pair4_sel_asym",
    /* 72 */ "_charged_mom_cdc_proton_from_lambda_pair4_sel_mqsel",
    /* 73 */ "_charged_mom_cdc_proton_from_lambda_pair4_selmeas",
    /* 74 */ "_charged_mom_cdc_proton_from_lambda_pair4_selmeas_asym",
    /* 75 */ "_charged_mom_cdc_proton_from_lambda_pair4_selmeas_mqsel",
    /* 76 */ "_charged_mom_cdc_proton_from_lambda_pair5",
    /* 77 */ "_charged_mom_cdc_proton_from_lambda_pair5_asym",
    /* 78 */ "_charged_mom_cdc_proton_from_lambda_pair5_sel",
    /* 79 */ "_charged_mom_cdc_proton_from_lambda_pair5_sel_asym",
    /* 80 */ "_charged_mom_cdc_proton_from_lambda_pair5_sel_mqsel",
    /* 81 */ "_charged_mom_cdc_proton_from_lambda_pair5_selmeas",
    /* 82 */ "_charged_mom_cdc_proton_from_lambda_pair5_selmeas_asym",
    /* 83 */ "_charged_mom_cdc_proton_from_lambda_pair5_selmeas_mqsel",
  };

  for(int i_sel=0; i_sel<84; ++i_sel){
    // generated
    new TH1F(Form("lambda_momentum%s",sel_name[i_sel].Data()),"lambda momentum;lambda momentum (GeV/c);counts",200,0.,2.);
    new TH1F(Form("proton_momentum%s",sel_name[i_sel].Data()),"proton momentum;proton momentum (GeV/c);counts",200,0.,2.);
    new TH1F(Form("proton_from_lambda_momentum%s",sel_name[i_sel].Data()),"proton from lambda momentum;proton from #Lambda momentum (GeV/c);counts",200,0.,2.);
    new TH1F(Form("pim_from_lambda_momentum%s",sel_name[i_sel].Data()),"pim from lambda momentum;#pi^{#minus} from #Lambda momentum (GeV/c);counts",200,0.,2.);
    new TH1F(Form("lp_mass%s",sel_name[i_sel].Data()),"lp mass;#Lambdap mass (GeV/c^{2});counts",100,2.,3.);
    new TH1F(Form("lp_momentum%s",sel_name[i_sel].Data()),"lp momentum;#Lambdap momentum (GeV/c);counts",200,0.,2.);
    new TH1F(Form("ppim_mass%s",sel_name[i_sel].Data()),"ppim mass;p#pi^{#minus} mass (GeV/c^{2});counts",100,1.,2.);
    new TH1F(Form("ppim_momentum%s",sel_name[i_sel].Data()),"ppim momentum;p#pi^{#minus} momentum (GeV/c);counts",200,0.,2.);

    new TH2F(Form("lp_mass_vs_lp_momentum%s",sel_name[i_sel].Data()),"lp mass vs lp momentum;#Lambdap mass (GeV/c^{2});#Lambdap momentum;counts",100,2.,3.,100,0.,1.);


    // momentum at detector
    new TH1F(Form("proton_momentum_at_cdc%s",sel_name[i_sel].Data()),"proton momentum at cdc;proton momentum (GeV/c);counts",200,0.,2.);
    new TH1F(Form("proton_from_lambda_momentum_at_cdc%s",sel_name[i_sel].Data()),"proton_from_lambda momentum at cdc;proton_from_lambda momentum (GeV/c);counts",200,0.,2.);
    for(int layer=1; layer<=10; ++layer){
      new TH1F(Form("proton_momentum_at_tracker_layer%d%s",layer,sel_name[i_sel].Data()),Form("proton momentum at tracker_layer%d;proton momentum (GeV/c);counts",layer),200,0.,2.);
      new TH1F(Form("proton_from_lambda_momentum_at_tracker_layer%d%s",layer,sel_name[i_sel].Data()),Form("proton_from_lambda momentum at tracker_layer%d;proton_from_lambda momentum (GeV/c);counts",layer),200,0.,2.);
    }


    // momenutm difference
    for(int pair=0; pair<5; ++pair){ 
      if(pair==0){
        new TH1F(Form("proton_momentum_difference_between_cdc_and_tracker_layer1%s",sel_name[i_sel].Data()),"proton momentum difference between cdc and tracker layer1;proton momentum difference(GeV/c);counts",200,0.,0.5);
        new TH1F(Form("proton_from_lambda_momentum_difference_between_cdc_and_tracker_layer1%s",sel_name[i_sel].Data()),"proton_from_lambda momentum difference between cdc and tracker layer1;proton_from_lambda momentum difference(GeV/c);counts",200,0.,0.5);
      }
      else{
        int first_layer  = pair*2;
        int second_layer = first_layer+1;
        new TH1F(Form("proton_momentum_difference_between_tracker_layer%d_and_layer%d%s",first_layer,second_layer,sel_name[i_sel].Data()),Form("proton momentum difference between tracker layer%d and layer%d;proton momentum difference(GeV/c);counts",first_layer,second_layer),200,0.,0.5);
        new TH1F(Form("proton_from_lambda_momentum_difference_between_tracker_layer%d_and_layer%d%s",first_layer,second_layer,sel_name[i_sel].Data()),Form("proton_from_lambda momentum difference between tracker layer%d and layer%d;proton_from_lambda momentum difference(GeV/c);counts",first_layer,second_layer),200,0.,0.5);
      }
    }


    // scattering angle and phi
    for(int pair=1; pair<=5; ++pair){
      new TH1F(Form("proton_scattering_angle_theta_mc_pair%d%s",pair,sel_name[i_sel].Data()),"scattering angle theta mc value;#theta (deg.);counts",200,0.,50.);
      new TH1F(Form("proton_scattering_angle_theta_pair%d%s",pair,sel_name[i_sel].Data()),"scattering angle theta;#theta (deg.);counts",200,0.,50.);
      new TH1F(Form("proton_scattering_angle_theta_meas_pair%d%s",pair,sel_name[i_sel].Data()),"scattering angle theta;#theta (deg.);counts",200,0.,50.);
      new TH1F(Form("proton_phi_of_spins_pair%d%s",pair,sel_name[i_sel].Data()),"phi between spins;#phi (rad.);counts",10,-1.,1.);
      new TH1F(Form("proton_phi_of_spins_meas_pair%d%s",pair,sel_name[i_sel].Data()),"phi between spins;#phi (rad.);counts",10,-1.,1.);
      new TH2F(Form("proton_phi_of_spins_pair%d_vs_scattering_angle_theta_pair%d%s",pair,pair,sel_name[i_sel].Data()),"phi between spins;#phi (rad.);counts",10,-1.,1.,200,0.,50.);
      new TH2F(Form("proton_phi_of_spins_meas_pair%d_vs_scattering_angle_theta_meas_pair%d%s",pair,pair,sel_name[i_sel].Data()),"phi between spins;#phi (rad.);counts",10,-1.,1.,200,0.,50.);

      new TH1F(Form("proton_from_lambda_scattering_angle_theta_mc_pair%d%s",pair,sel_name[i_sel].Data()),"scattering angle theta mc value;#theta (deg.);counts",200,0.,50.);
      new TH1F(Form("proton_from_lambda_scattering_angle_theta_pair%d%s",pair,sel_name[i_sel].Data()),"scattering angle theta;#theta (deg.);counts",200,0.,50.);
      new TH1F(Form("proton_from_lambda_scattering_angle_theta_meas_pair%d%s",pair,sel_name[i_sel].Data()),"scattering angle theta;#theta (deg.);counts",200,0.,50.);
      new TH1F(Form("proton_from_lambda_phi_of_spins_pair%d%s",pair,sel_name[i_sel].Data()),"phi between spins;#phi (rad.);counts",10,-1.,1.);
      new TH1F(Form("proton_from_lambda_phi_of_spins_meas_pair%d%s",pair,sel_name[i_sel].Data()),"phi between spins;#phi (rad.);counts",10,-1.,1.);
      new TH2F(Form("proton_from_lambda_phi_of_spins_pair%d_vs_scattering_angle_theta_pair%d%s",pair,pair,sel_name[i_sel].Data()),"phi between spins;#phi (rad.);counts",10,-1.,1.,200,0.,50.);
      new TH2F(Form("proton_from_lambda_phi_of_spins_meas_pair%d_vs_scattering_angle_theta_meas_pair%d%s",pair,pair,sel_name[i_sel].Data()),"phi between spins;#phi (rad.);counts",10,-1.,1.,200,0.,50.);
    }

    // resolution
    new TH1F(Form("proton_diff_phi_direction_at_cdc%s",sel_name[i_sel].Data()),"proton diff of phi;#phi (deg.);counts",200,-1.,1.);
    new TH1F(Form("proton_diff_theta_direction_at_cdc%s",sel_name[i_sel].Data()),"proton diff of theta;#theta (deg.);counts",200,-5.,5.);
    new TH1F(Form("proton_diff_phirho_position_at_cdc%s",sel_name[i_sel].Data()),"proton diff of phirho;#phi#rho (cm);counts",200,-0.2,0.2);
    new TH1F(Form("proton_diff_z_position_at_cdc%s",sel_name[i_sel].Data()),"proton diff of z;z (cm);counts",200,-2.,2.);
    new TH2F(Form("proton_diff_phi_direction_at_cdc_vs_diff_phirho_position_at_cdc%s",sel_name[i_sel].Data()),"proton diff of phi vs phirho;#phi (deg.);#phi#rho (cm);Counts",200,-1.,1.,200,-0.2,0.2);
    new TH2F(Form("proton_diff_theta_direction_at_cdc_vs_diff_z_position_at_cdc%s",sel_name[i_sel].Data()),"proton diff of theta vs z;#theta (deg.);z (cm);Counts",200,-5.,5.,200,-2.,2.);
    // -----
    for(int layer=1; layer<=10; layer++){
      new TH1F(Form("proton_diff_phirho_position_at_tracker_layer%d%s",layer,sel_name[i_sel].Data()),"proton diff of phirho;#phi#rho (cm);counts",200,-1.,1.);
      new TH1F(Form("proton_diff_z_position_at_tracker_layer%d%s",layer,sel_name[i_sel].Data()),"proton diff of z;z (cm);counts",200,-1.,1.);
    }
  }
}

void InitializeTrees(TFile* infile){
  infile->cd();

  tree_event = (TTree*)infile->Get("tree_event");     
  tree_event->SetBranchAddress("event_id",&event_id);
  tree_event->SetBranchAddress("number_of_primaries",&number_of_primaries);
  tree_event->SetBranchAddress("number_of_trajectories",&number_of_trajectories);
  tree_event->SetBranchAddress("number_of_hits_in_cdc",&number_of_hits_in_cdc);
  tree_event->SetBranchAddress("number_of_hits_in_tracker_layer1",&number_of_hits_in_tracker[0]);
  tree_event->SetBranchAddress("number_of_hits_in_tracker_layer2",&number_of_hits_in_tracker[1]);
  tree_event->SetBranchAddress("number_of_hits_in_tracker_layer3",&number_of_hits_in_tracker[2]);
  tree_event->SetBranchAddress("number_of_hits_in_tracker_layer4",&number_of_hits_in_tracker[3]);
  tree_event->SetBranchAddress("number_of_hits_in_tracker_layer5",&number_of_hits_in_tracker[4]);
  tree_event->SetBranchAddress("number_of_hits_in_tracker_layer6",&number_of_hits_in_tracker[5]);
  tree_event->SetBranchAddress("number_of_hits_in_tracker_layer7",&number_of_hits_in_tracker[6]);
  tree_event->SetBranchAddress("number_of_hits_in_tracker_layer8",&number_of_hits_in_tracker[7]);
  tree_event->SetBranchAddress("number_of_hits_in_tracker_layer9",&number_of_hits_in_tracker[8]);
  tree_event->SetBranchAddress("number_of_hits_in_tracker_layer10",&number_of_hits_in_tracker[9]);
  // -----
  tree_primary = (TTree*)infile->Get("tree_primary");     
  tree_primary->SetBranchAddress("event_id",&primary_event_id);
  tree_primary->SetBranchAddress("particle_id",&primary_particle_id);
  tree_primary->SetBranchAddress("position_x",&primary_position_x);
  tree_primary->SetBranchAddress("position_y",&primary_position_y);
  tree_primary->SetBranchAddress("position_z",&primary_position_z);
  tree_primary->SetBranchAddress("momentum_x",&primary_momentum_x);
  tree_primary->SetBranchAddress("momentum_y",&primary_momentum_y);
  tree_primary->SetBranchAddress("momentum_z",&primary_momentum_z);
  tree_primary->SetBranchAddress("spin_x",&primary_spin_x);
  tree_primary->SetBranchAddress("spin_y",&primary_spin_y);
  tree_primary->SetBranchAddress("spin_z",&primary_spin_z);
  // -----
  tree_trajectory = (TTree*)infile->Get("tree_trajectory");     
  tree_trajectory->SetBranchAddress("event_id",&trajectory_event_id);
  tree_trajectory->SetBranchAddress("track_id",&trajectory_track_id);
  tree_trajectory->SetBranchAddress("parent_id",&trajectory_parent_id);
  tree_trajectory->SetBranchAddress("particle_id",&trajectory_particle_id);
  tree_trajectory->SetBranchAddress("initial_position_x",&trajectory_initial_position_x);
  tree_trajectory->SetBranchAddress("initial_position_y",&trajectory_initial_position_y);
  tree_trajectory->SetBranchAddress("initial_position_z",&trajectory_initial_position_z);
  tree_trajectory->SetBranchAddress("initial_momentum_x",&trajectory_initial_momentum_x);
  tree_trajectory->SetBranchAddress("initial_momentum_y",&trajectory_initial_momentum_y);
  tree_trajectory->SetBranchAddress("initial_momentum_z",&trajectory_initial_momentum_z);
  // -----
  tree_cdc_hit = (TTree*)infile->Get("tree_cdc_hit");     
  tree_cdc_hit->SetBranchAddress("event_id",&cdc_event_id);
  tree_cdc_hit->SetBranchAddress("track_id",&cdc_track_id);
  tree_cdc_hit->SetBranchAddress("parent_id",&cdc_parent_id);
  tree_cdc_hit->SetBranchAddress("particle_id",&cdc_particle_id);
  tree_cdc_hit->SetBranchAddress("layer_id",&cdc_layer_id);
  tree_cdc_hit->SetBranchAddress("hit_time",&cdc_hit_time);
  tree_cdc_hit->SetBranchAddress("hit_position_x",&cdc_hit_position_x);
  tree_cdc_hit->SetBranchAddress("hit_position_y",&cdc_hit_position_y);
  tree_cdc_hit->SetBranchAddress("hit_position_z",&cdc_hit_position_z);
  tree_cdc_hit->SetBranchAddress("momentum_x",&cdc_momentum_x);
  tree_cdc_hit->SetBranchAddress("momentum_y",&cdc_momentum_y);
  tree_cdc_hit->SetBranchAddress("momentum_z",&cdc_momentum_z);
  tree_cdc_hit->SetBranchAddress("is_asymmetric_scattering",&cdc_is_asymmetric_scattering);
  // -----
  for(int layer=0; layer<10; ++layer){
    tree_tracker_hit[layer] = (TTree*)infile->Get(Form("tree_tracker_layer%d_hit",layer+1));     
    tree_tracker_hit[layer]->SetBranchAddress("event_id",&tracker_event_id[layer]);
    tree_tracker_hit[layer]->SetBranchAddress("track_id",&tracker_track_id[layer]);
    tree_tracker_hit[layer]->SetBranchAddress("parent_id",&tracker_parent_id[layer]);
    tree_tracker_hit[layer]->SetBranchAddress("particle_id",&tracker_particle_id[layer]);
    tree_tracker_hit[layer]->SetBranchAddress("layer_id",&tracker_layer_id[layer]);
    tree_tracker_hit[layer]->SetBranchAddress("hit_time",&tracker_hit_time[layer]);
    tree_tracker_hit[layer]->SetBranchAddress("hit_position_x",&tracker_hit_position_x[layer]);
    tree_tracker_hit[layer]->SetBranchAddress("hit_position_y",&tracker_hit_position_y[layer]);
    tree_tracker_hit[layer]->SetBranchAddress("hit_position_z",&tracker_hit_position_z[layer]);
    tree_tracker_hit[layer]->SetBranchAddress("momentum_x",&tracker_momentum_x[layer]);
    tree_tracker_hit[layer]->SetBranchAddress("momentum_y",&tracker_momentum_y[layer]);
    tree_tracker_hit[layer]->SetBranchAddress("momentum_z",&tracker_momentum_z[layer]);
    tree_tracker_hit[layer]->SetBranchAddress("is_asymmetric_scattering",&tracker_is_asymmetric_scattering[layer]);
  }
  // -----
}

void Analysis(TFile* outfile){
  Long64_t total_entries = tree_event->GetEntries();
  Int_t i_entry_primary = 0;
  Int_t i_entry_trajectory = 0;
  Int_t i_entry_cdc = 0;
  Int_t i_entry_tracker[10] = {0};

  for(Long64_t i_entry=0; i_entry<total_entries; ++i_entry){
    // begin of event loop

    if((i_entry+1)%(Int_t)pow(10, (Int_t)log10(i_entry+1)) ==0 ){
      std::cout << ">> " << i_entry+1 << std::endl;
    }
    tree_event->GetEntry(i_entry);


    // ====================================================================================================
    // Read trees
    // ====================================================================================================

    // primary
    std::vector<Int_t> primary_particle_ids;
    std::vector<TVector3> vec_primary_momenta;
    std::vector<TVector3> vec_primary_spins;
    for(int i_primary=0; i_primary<number_of_primaries; ++i_primary){
      tree_primary->GetEntry(i_entry_primary);
      ++i_entry_primary;
      if(event_id != primary_event_id){
        std::cout << "event mismatch!!!" << std::endl;
        return;
      }
      primary_particle_ids.push_back(primary_particle_id);
      TVector3 vec_primary_momentum((Double_t)primary_momentum_x,(Double_t)primary_momentum_y,(Double_t)primary_momentum_z);
      TVector3 vec_primary_spin((Double_t)primary_spin_x,(Double_t)primary_spin_y,(Double_t)primary_spin_z);
      vec_primary_momenta.push_back(vec_primary_momentum);
      vec_primary_spins.push_back(vec_primary_spin);
    }

    // trajectory
    std::vector<Int_t> trajectory_track_ids;
    std::vector<Int_t> trajectory_parent_ids;
    std::vector<Int_t> trajectory_particle_ids;
    std::vector<TVector3> vec_trajectory_initial_momenta;
    std::vector<TVector3> vec_trajectory_initial_positions;
    for(int i_trajectory=0; i_trajectory<number_of_trajectories; ++i_trajectory){
      tree_trajectory->GetEntry(i_entry_trajectory);
      ++i_entry_trajectory;
      if(event_id != trajectory_event_id){
        std::cout << "event mismatch!!!" << std::endl;
        return;
      }
      trajectory_track_ids.push_back(trajectory_track_id);
      trajectory_parent_ids.push_back(trajectory_parent_id);
      trajectory_particle_ids.push_back(trajectory_particle_id);
      TVector3 vec_trajectory_initial_position((Double_t)trajectory_initial_position_x,(Double_t)trajectory_initial_position_y,(Double_t)trajectory_initial_position_z);
      TVector3 vec_trajectory_initial_momentum((Double_t)trajectory_initial_momentum_x,(Double_t)trajectory_initial_momentum_y,(Double_t)trajectory_initial_momentum_z);
      vec_trajectory_initial_positions.push_back(vec_trajectory_initial_position);
      vec_trajectory_initial_momenta.push_back(vec_trajectory_initial_momentum);
    }

    // cdc
    std::vector<Int_t> cdc_parent_ids;
    std::vector<Int_t> cdc_particle_ids;
    std::vector<TVector3> vec_cdc_positions;
    std::vector<TVector3> vec_cdc_momenta;
    std::vector<Int_t> cdc_are_asymmetric_scattering;
    for(int i_cdc=0; i_cdc<number_of_hits_in_cdc; ++i_cdc){
      tree_cdc_hit->GetEntry(i_entry_cdc);
      ++i_entry_cdc;
      if(event_id != cdc_event_id){
        std::cout << "event mismatch!!!" << std::endl;
        return;
      }
      cdc_parent_ids.push_back(cdc_parent_id);
      cdc_particle_ids.push_back(cdc_particle_id);
      TVector3 vec_cdc_position((Double_t)cdc_hit_position_x,(Double_t)cdc_hit_position_y,(Double_t)cdc_hit_position_z);
      TVector3 vec_cdc_momentum((Double_t)cdc_momentum_x,(Double_t)cdc_momentum_y,(Double_t)cdc_momentum_z);
      vec_cdc_positions.push_back(vec_cdc_position);
      vec_cdc_momenta.push_back(vec_cdc_momentum);
      cdc_are_asymmetric_scattering.push_back(cdc_is_asymmetric_scattering);
    }

    // tracker
    std::vector< std::vector<Int_t> > tracker_parent_ids(10);
    std::vector< std::vector<Int_t> > tracker_particle_ids(10);
    std::vector< std::vector<TVector3> > vec_tracker_positions(10);
    std::vector< std::vector<TVector3> > vec_tracker_momenta(10);
    std::vector< std::vector<Int_t> > tracker_are_asymmetric_scattering(10);
    for(int layer=0; layer<10; ++layer){
      for(int i_tracker=0; i_tracker<number_of_hits_in_tracker[layer]; ++i_tracker){
        tree_tracker_hit[layer]->GetEntry(i_entry_tracker[layer]);
        ++i_entry_tracker[layer];
        if(event_id != tracker_event_id[layer]){
          std::cout << "event mismatch!!!" << std::endl;
          return;
        }
        tracker_parent_ids[layer].push_back(tracker_parent_id[layer]);
        tracker_particle_ids[layer].push_back(tracker_particle_id[layer]);
        TVector3 vec_tracker_position((Double_t)tracker_hit_position_x[layer],(Double_t)tracker_hit_position_y[layer],(Double_t)tracker_hit_position_z[layer]);
        TVector3 vec_tracker_momentum((Double_t)tracker_momentum_x[layer],(Double_t)tracker_momentum_y[layer],(Double_t)tracker_momentum_z[layer]);
        vec_tracker_positions[layer].push_back(vec_tracker_position);
        vec_tracker_momenta[layer].push_back(vec_tracker_momentum);
        tracker_are_asymmetric_scattering[layer].push_back(tracker_is_asymmetric_scattering[layer]);
      }
    }

    // ====================================================================================================
    // ====================================================================================================





    // ====================================================================================================
    // Analysis
    // ====================================================================================================

    // primary
    TVector3 vec_lambda_spin_direction_mc(0.,0.,0.);
    TVector3 vec_proton_spin_direction_mc(0.,0.,0.);
    for(int i_primary=0; i_primary<number_of_primaries; ++i_primary){
      if(primary_particle_ids[i_primary]==3122){
        vec_lambda_spin_direction_mc = vec_primary_spins[i_primary];
      }
      if(primary_particle_ids[i_primary]==2212){
        vec_proton_spin_direction_mc = vec_primary_spins[i_primary];
      }
    }

    // trajectory
    Int_t lambda_id = -999;
    Int_t proton_id = -999;
    TLorentzVector lvec_lambda_mc(0.,0.,0.,0.);
    TLorentzVector lvec_proton_mc(0.,0.,0.,0.);
    for(int i_trajectory=0; i_trajectory<number_of_trajectories; ++i_trajectory){
      if(trajectory_parent_ids[i_trajectory]==0&&trajectory_particle_ids[i_trajectory]==3122){
        lambda_id = trajectory_track_ids[i_trajectory];
        lvec_lambda_mc.SetVectM(vec_trajectory_initial_momenta[i_trajectory], mass::lambda);
      }
      if(trajectory_parent_ids[i_trajectory]==0&&trajectory_particle_ids[i_trajectory]==2212){
        proton_id = trajectory_track_ids[i_trajectory];
        lvec_proton_mc.SetVectM(vec_trajectory_initial_momenta[i_trajectory], mass::proton);
      }
    }
    TLorentzVector lvec_lp_mc = lvec_proton_mc + lvec_lambda_mc;

    Int_t proton_from_lambda_id = -999;
    TLorentzVector lvec_proton_from_lambda_mc(0.,0.,0.,0.);
    TLorentzVector lvec_pim_from_lambda_mc(0.,0.,0.,0.);
    bool is_lambda_charged_decay = false;
    for(int i_trajectory=0; i_trajectory<number_of_trajectories; ++i_trajectory){
      if(trajectory_parent_ids[i_trajectory]==lambda_id&&trajectory_particle_ids[i_trajectory]==2212){
        lvec_proton_from_lambda_mc.SetVectM(vec_trajectory_initial_momenta[i_trajectory], mass::proton);
        proton_from_lambda_id = trajectory_track_ids[i_trajectory];
        is_lambda_charged_decay = true;
      }
      if(trajectory_parent_ids[i_trajectory]==lambda_id&&trajectory_particle_ids[i_trajectory]==-211){
        lvec_pim_from_lambda_mc.SetVectM(vec_trajectory_initial_momenta[i_trajectory], mass::pi_plus);
      }
    }
    TLorentzVector lvec_ppim_mc = lvec_proton_from_lambda_mc + lvec_pim_from_lambda_mc;

    TLorentzVector clvec_proton_from_lambda_lambda_rest_mc = lvec_proton_from_lambda_mc;
    clvec_proton_from_lambda_lambda_rest_mc.Boost(-lvec_lambda_mc.BoostVector());

    // checking cdc (two protons and pim to be detected by cdc)
    bool is_proton_detected_by_cdc = false;
    bool is_proton_from_lambda_detected_by_cdc = false;
    bool is_pim_from_lambda_detected_by_cdc = false;
    for(int i_hit=0; i_hit<number_of_hits_in_cdc; ++i_hit){
      if(cdc_parent_ids[i_hit]==0&&cdc_particle_ids[i_hit]==2212){
        is_proton_detected_by_cdc = true;
      }
      if(cdc_parent_ids[i_hit]==lambda_id&&cdc_particle_ids[i_hit]==2212){
        is_proton_from_lambda_detected_by_cdc = true;
      }
      if(cdc_parent_ids[i_hit]==lambda_id&&cdc_particle_ids[i_hit]==-211){
        is_pim_from_lambda_detected_by_cdc = true;
      }
    }

    // checking tracker (proton to be detected)
    bool is_proton_detected_by_tracker[10] = {false, false, false, false, false, false, false, false, false, false};
    bool is_proton_asymmetric_scattering[10] = {false, false, false, false, false, false, false, false, false, false};
    bool is_proton_from_lambda_detected_by_tracker[10] = {false, false, false, false, false, false, false, false, false, false};
    bool is_proton_from_lambda_asymmetric_scattering[10] = {false, false, false, false, false, false, false, false, false, false};
    for(int layer=0; layer<10; ++layer){
      for(int i_hit=0; i_hit<number_of_hits_in_tracker[layer]; ++i_hit){
        if(tracker_parent_ids[layer][i_hit]==0 || tracker_parent_ids[layer][i_hit]==proton_id){
          if(tracker_particle_ids[layer][i_hit]==2212){
            is_proton_detected_by_tracker[layer] = true;
            if(tracker_are_asymmetric_scattering[layer][i_hit]){
              is_proton_asymmetric_scattering[layer] = true;
            }
          }
        }
        else if(tracker_parent_ids[layer][i_hit]==lambda_id || tracker_parent_ids[layer][i_hit]==proton_from_lambda_id){
          if(tracker_particle_ids[layer][i_hit]==2212){
            is_proton_from_lambda_detected_by_tracker[layer] = true;
            if(tracker_are_asymmetric_scattering[layer][i_hit]){
              is_proton_from_lambda_asymmetric_scattering[layer] = true;
            }
          }
        }
      }
    }

    // checking cdc, tracker1, and tracker2 (to measure proton scattering angle)
    TVector3 vec_proton_position_at_cdc(0.,0.,0.);
    TVector3 vec_proton_momentum_at_cdc(0.,0.,0.);
    TVector3 vec_proton_position_at_tracker[10];
    TVector3 vec_proton_momentum_at_tracker[10];

    TVector3 vec_proton_from_lambda_position_at_cdc(0.,0.,0.);
    TVector3 vec_proton_from_lambda_momentum_at_cdc(0.,0.,0.);
    TVector3 vec_proton_from_lambda_position_at_tracker[10];
    TVector3 vec_proton_from_lambda_momentum_at_tracker[10];

    // position & momentum at cdc
    for(int i_hit=0; i_hit<number_of_hits_in_cdc; ++i_hit){
      if(cdc_parent_ids[i_hit]==0&&cdc_particle_ids[i_hit]==2212){
        vec_proton_momentum_at_cdc = vec_cdc_momenta[i_hit];
        vec_proton_position_at_cdc = vec_cdc_positions[i_hit];
      }
      else if(cdc_parent_ids[i_hit]==lambda_id&&cdc_particle_ids[i_hit]==2212){
        vec_proton_from_lambda_momentum_at_cdc = vec_cdc_momenta[i_hit];
        vec_proton_from_lambda_position_at_cdc = vec_cdc_positions[i_hit];
      }
    }

    // position & momentum at tracker
    for(int layer=0; layer<10; ++layer){
      if(is_proton_detected_by_tracker[layer]){
        for(int i_hit=0; i_hit<number_of_hits_in_tracker[layer]; ++i_hit){
          if(tracker_parent_ids[layer][i_hit]==0&&tracker_particle_ids[layer][i_hit]==2212){
            vec_proton_momentum_at_tracker[layer] = vec_tracker_momenta[layer][i_hit];
            vec_proton_position_at_tracker[layer] = vec_tracker_positions[layer][i_hit];
          }
        }
        if(vec_proton_position_at_tracker[layer].Mag()==0.){ // proton would be vanished.
          for(int i_hit=0; i_hit<number_of_hits_in_tracker[layer]; ++i_hit){
            if(tracker_parent_ids[layer][i_hit]==proton_id&&tracker_particle_ids[layer][i_hit]==2212){
              vec_proton_momentum_at_tracker[layer] = vec_tracker_momenta[layer][i_hit];
              vec_proton_position_at_tracker[layer] = vec_tracker_positions[layer][i_hit];
            }
          }
        }
      }

      if(is_proton_from_lambda_detected_by_tracker[layer]){
        for(int i_hit=0; i_hit<number_of_hits_in_tracker[layer]; ++i_hit){
          if(tracker_parent_ids[layer][i_hit]==lambda_id&&tracker_particle_ids[layer][i_hit]==2212){
            vec_proton_from_lambda_momentum_at_tracker[layer] = vec_tracker_momenta[layer][i_hit];
            vec_proton_from_lambda_position_at_tracker[layer] = vec_tracker_positions[layer][i_hit];
          }
        }
        if(vec_proton_from_lambda_position_at_tracker[layer].Mag()==0.){ // proton would be vanished.
          for(int i_hit=0; i_hit<number_of_hits_in_tracker[layer]; ++i_hit){
            if(tracker_parent_ids[layer][i_hit]==proton_from_lambda_id&&tracker_particle_ids[layer][i_hit]==2212){
              vec_proton_from_lambda_momentum_at_tracker[layer] = vec_tracker_momenta[layer][i_hit];
              vec_proton_from_lambda_position_at_tracker[layer] = vec_tracker_positions[layer][i_hit];
            }
          }
        }
      }
    }

    // obtained directions
    TVector3 vec_proton_direction_at_cdc = vec_proton_momentum_at_cdc.Unit();
    TVector3 vec_proton_from_lambda_direction_at_cdc = vec_proton_from_lambda_momentum_at_cdc.Unit();
    TVector3 vec_proton_direction_at_tracker[10];
    TVector3 vec_proton_from_lambda_direction_at_tracker[10];
    for(int layer=0; layer<10; ++layer){
      vec_proton_direction_at_tracker[layer] = vec_proton_momentum_at_tracker[layer].Unit();
      vec_proton_from_lambda_direction_at_tracker[layer] = vec_proton_from_lambda_momentum_at_tracker[layer].Unit();
    }

    // resolution
    TVector3 vec_proton_direction_at_cdc_meas = vec_proton_direction_at_cdc;
    TVector3 vec_proton_position_at_cdc_meas = vec_proton_position_at_cdc;
    CDCResolution(vec_proton_position_at_cdc_meas,vec_proton_direction_at_cdc_meas);
    // ----- 
    TVector3 vec_proton_from_lambda_direction_at_cdc_meas = vec_proton_from_lambda_direction_at_cdc;
    TVector3 vec_proton_from_lambda_position_at_cdc_meas = vec_proton_from_lambda_position_at_cdc;
    CDCResolution(vec_proton_from_lambda_position_at_cdc_meas,vec_proton_from_lambda_direction_at_cdc_meas);
    // ----- 
    TVector3 vec_proton_position_at_tracker_meas[10];
    TVector3 vec_proton_from_lambda_position_at_tracker_meas[10];
    for(int layer=0; layer<10; ++layer){ 
      vec_proton_position_at_tracker_meas[layer] = vec_proton_position_at_tracker[layer];
      TrackerResolution(vec_proton_position_at_tracker_meas[layer]);
      vec_proton_from_lambda_position_at_tracker_meas[layer] = vec_proton_from_lambda_position_at_tracker[layer];
      TrackerResolution(vec_proton_from_lambda_position_at_tracker_meas[layer]);
    }

    // differencde between true and measured
    Double_t proton_diff_phi_direction_at_cdc = vec_proton_direction_at_cdc.Phi() - vec_proton_direction_at_cdc_meas.Phi();
    Double_t proton_diff_theta_direction_at_cdc = vec_proton_direction_at_cdc.Theta() - vec_proton_direction_at_cdc_meas.Theta();
    Double_t proton_diff_phirho_position_at_cdc = (vec_proton_position_at_cdc.Phi() - vec_proton_position_at_cdc_meas.Phi()) * vec_proton_position_at_cdc.Perp();
    Double_t proton_diff_z_position_at_cdc = -(vec_proton_position_at_cdc.Z() - vec_proton_position_at_cdc_meas.Z());
    // ----- 
    Double_t proton_diff_phirho_position_at_tracker[10];
    Double_t proton_diff_z_position_at_tracker[10];
    for(int layer=0; layer<10; ++layer){
      proton_diff_phirho_position_at_tracker[layer] = (vec_proton_position_at_tracker[layer].Phi() - vec_proton_position_at_tracker_meas[layer].Phi()) * vec_proton_position_at_tracker[layer].Perp();
      proton_diff_z_position_at_tracker[layer] = -(vec_proton_position_at_tracker[layer].Z() - vec_proton_position_at_tracker_meas[layer].Z());
    }
    // ----- 


    // obtained direction after scattering
    TVector3 vec_proton_direction_between_trackers[5];
    TVector3 vec_proton_from_lambda_direction_between_trackers[5];
    TVector3 vec_proton_direction_between_trackers_meas[5];
    TVector3 vec_proton_from_lambda_direction_between_trackers_meas[5];
    for(int pair=0; pair<5; ++pair){
      int first_layer  = pair*2;
      int second_layer = first_layer+1;
      vec_proton_direction_between_trackers[pair] = (vec_proton_position_at_tracker[second_layer] - vec_proton_position_at_tracker[first_layer]).Unit();
      vec_proton_from_lambda_direction_between_trackers[pair] = (vec_proton_from_lambda_position_at_tracker[second_layer] - vec_proton_from_lambda_position_at_tracker[first_layer]).Unit();
      // -----
      vec_proton_direction_between_trackers_meas[pair] = (vec_proton_position_at_tracker_meas[second_layer] - vec_proton_position_at_tracker_meas[first_layer]).Unit();
      vec_proton_from_lambda_direction_between_trackers_meas[pair] = (vec_proton_from_lambda_position_at_tracker_meas[second_layer] - vec_proton_from_lambda_position_at_tracker_meas[first_layer]).Unit();
    }
    // -----


    // expected spin direction of lambda
    TVector3 vec_lambda_expected_spin_direction = clvec_proton_from_lambda_lambda_rest_mc.Vect().Unit();
    TVector3 vec_lambda_expected_spin_direction_perp_to_proton[5];
    TVector3 vec_lambda_expected_spin_direction_perp_to_proton_from_lambda[5];
    TVector3 vec_lambda_expected_spin_direction_perp_to_proton_meas[5];
    TVector3 vec_lambda_expected_spin_direction_perp_to_proton_from_lambda_meas[5];
    for(int pair=0; pair<5; ++pair){
      if(pair==0){
        vec_lambda_expected_spin_direction_perp_to_proton[pair] = (vec_lambda_expected_spin_direction - (vec_lambda_expected_spin_direction.Dot(vec_proton_direction_at_cdc))*vec_proton_direction_at_cdc).Unit();
        vec_lambda_expected_spin_direction_perp_to_proton_from_lambda[pair] = (vec_lambda_expected_spin_direction - (vec_lambda_expected_spin_direction.Dot(vec_proton_from_lambda_direction_at_cdc))*vec_proton_from_lambda_direction_at_cdc).Unit();
        vec_lambda_expected_spin_direction_perp_to_proton_meas[pair] = (vec_lambda_expected_spin_direction - (vec_lambda_expected_spin_direction.Dot(vec_proton_direction_at_cdc_meas))*vec_proton_direction_at_cdc_meas).Unit();
        vec_lambda_expected_spin_direction_perp_to_proton_from_lambda_meas[pair] = (vec_lambda_expected_spin_direction - (vec_lambda_expected_spin_direction.Dot(vec_proton_from_lambda_direction_at_cdc_meas))*vec_proton_from_lambda_direction_at_cdc_meas).Unit();
      }
      else{
        vec_lambda_expected_spin_direction_perp_to_proton[pair] = (vec_lambda_expected_spin_direction - (vec_lambda_expected_spin_direction.Dot(vec_proton_direction_between_trackers[pair-1]))*vec_proton_direction_between_trackers[pair-1]).Unit();
        vec_lambda_expected_spin_direction_perp_to_proton_from_lambda[pair] = (vec_lambda_expected_spin_direction - (vec_lambda_expected_spin_direction.Dot(vec_proton_from_lambda_direction_between_trackers[pair-1]))*vec_proton_from_lambda_direction_between_trackers[pair-1]).Unit();
        vec_lambda_expected_spin_direction_perp_to_proton_meas[pair] = (vec_lambda_expected_spin_direction - (vec_lambda_expected_spin_direction.Dot(vec_proton_direction_between_trackers_meas[pair-1]))*vec_proton_direction_between_trackers_meas[pair-1]).Unit();
        vec_lambda_expected_spin_direction_perp_to_proton_from_lambda_meas[pair] = (vec_lambda_expected_spin_direction - (vec_lambda_expected_spin_direction.Dot(vec_proton_from_lambda_direction_between_trackers_meas[pair-1]))*vec_proton_from_lambda_direction_between_trackers_meas[pair-1]).Unit();
      }
    }
    // -----


    // scattering angle and momentum difference (mc value)
    Double_t proton_scattering_angle_theta_mc[5];
    Double_t proton_from_lambda_scattering_angle_theta_mc[5];
    Double_t proton_momentum_difference[5];
    Double_t proton_from_lambda_momentum_difference[5];
    for(int pair=0; pair<5; pair++){
      if(pair==0){
        proton_scattering_angle_theta_mc[pair] = acos(vec_proton_direction_at_cdc.Dot(vec_proton_direction_at_tracker[pair])) /TMath::Pi()*180.;
        proton_from_lambda_scattering_angle_theta_mc[pair] = acos(vec_proton_from_lambda_direction_at_cdc.Dot(vec_proton_from_lambda_direction_at_tracker[pair])) /TMath::Pi()*180.;
        proton_momentum_difference[pair] = fabs(vec_proton_momentum_at_cdc.Mag() - vec_proton_momentum_at_tracker[pair].Mag());
        proton_from_lambda_momentum_difference[pair] = fabs(vec_proton_from_lambda_momentum_at_cdc.Mag() - vec_proton_from_lambda_momentum_at_tracker[pair].Mag());
      }
      else{
        int first_layer = pair*2-1; // 1, 3, 5
        int second_layer = first_layer+1; // 2, 4, 6
        proton_scattering_angle_theta_mc[pair] = acos(vec_proton_direction_at_tracker[first_layer].Dot(vec_proton_direction_at_tracker[second_layer])) /TMath::Pi()*180.;
        proton_from_lambda_scattering_angle_theta_mc[pair] = acos(vec_proton_from_lambda_direction_at_tracker[first_layer].Dot(vec_proton_from_lambda_direction_at_tracker[second_layer])) /TMath::Pi()*180.;
        proton_momentum_difference[pair] = fabs(vec_proton_momentum_at_tracker[first_layer].Mag() - vec_proton_momentum_at_tracker[second_layer].Mag());
        proton_from_lambda_momentum_difference[pair] = fabs(vec_proton_from_lambda_momentum_at_cdc.Mag() - vec_proton_from_lambda_momentum_at_tracker[second_layer].Mag());
      }
    }
    // -----


    // scattering angle with two trackers
    Double_t proton_scattering_angle_theta[5];
    Double_t proton_from_lambda_scattering_angle_theta[5];
    Double_t proton_scattering_angle_theta_meas[5];
    Double_t proton_from_lambda_scattering_angle_theta_meas[5];
    for(int pair=0; pair<5; ++pair){
      if(pair==0){
        proton_scattering_angle_theta[pair] = acos(vec_proton_direction_at_cdc.Dot(vec_proton_direction_between_trackers[pair])) /TMath::Pi()*180.;
        proton_from_lambda_scattering_angle_theta[pair] = acos(vec_proton_from_lambda_direction_at_cdc.Dot(vec_proton_from_lambda_direction_between_trackers[pair])) /TMath::Pi()*180.;
        proton_scattering_angle_theta_meas[pair] = acos(vec_proton_direction_at_cdc_meas.Dot(vec_proton_direction_between_trackers_meas[pair])) /TMath::Pi()*180.;
        proton_from_lambda_scattering_angle_theta_meas[pair] = acos(vec_proton_from_lambda_direction_at_cdc_meas.Dot(vec_proton_from_lambda_direction_between_trackers_meas[pair])) /TMath::Pi()*180.;
      }
      else{
        proton_scattering_angle_theta[pair] = acos(vec_proton_direction_between_trackers[pair-1].Dot(vec_proton_direction_between_trackers[pair])) /TMath::Pi()*180.;
        proton_from_lambda_scattering_angle_theta[pair] = acos(vec_proton_from_lambda_direction_between_trackers[pair-1].Dot(vec_proton_from_lambda_direction_between_trackers[pair])) /TMath::Pi()*180.;
        proton_scattering_angle_theta_meas[pair] = acos(vec_proton_direction_between_trackers_meas[pair-1].Dot(vec_proton_direction_between_trackers_meas[pair])) /TMath::Pi()*180.;
        proton_from_lambda_scattering_angle_theta_meas[pair] = acos(vec_proton_from_lambda_direction_between_trackers_meas[pair-1].Dot(vec_proton_from_lambda_direction_between_trackers_meas[pair])) /TMath::Pi()*180.;
      }
    }
    // -----


    // expected spin direction of proton with measured scattering direction
    TVector3 vec_proton_expected_spin_direction[5];
    TVector3 vec_proton_from_lambda_expected_spin_direction[5];
    TVector3 vec_proton_expected_spin_direction_meas[5];
    TVector3 vec_proton_from_lambda_expected_spin_direction_meas[5];
    for(int pair=0; pair<5; ++pair){
      if(pair==0){
        vec_proton_expected_spin_direction[pair] = vec_proton_direction_at_cdc.Cross(vec_proton_direction_between_trackers[pair]).Unit();
        vec_proton_from_lambda_expected_spin_direction[pair] = vec_proton_from_lambda_direction_at_cdc.Cross(vec_proton_from_lambda_direction_between_trackers[pair]).Unit();
        vec_proton_expected_spin_direction_meas[pair] = vec_proton_direction_at_cdc_meas.Cross(vec_proton_direction_between_trackers_meas[pair]).Unit();
        vec_proton_from_lambda_expected_spin_direction_meas[pair] = vec_proton_from_lambda_direction_at_cdc_meas.Cross(vec_proton_from_lambda_direction_between_trackers_meas[pair]).Unit();
      }
      else{
        vec_proton_expected_spin_direction[pair] = vec_proton_direction_between_trackers[pair-1].Cross(vec_proton_direction_between_trackers[pair]).Unit();
        vec_proton_from_lambda_expected_spin_direction[pair] = vec_proton_from_lambda_direction_between_trackers[pair-1].Cross(vec_proton_from_lambda_direction_between_trackers[pair]).Unit();
        vec_proton_expected_spin_direction_meas[pair] = vec_proton_direction_between_trackers_meas[pair-1].Cross(vec_proton_direction_between_trackers_meas[pair]).Unit();
        vec_proton_from_lambda_expected_spin_direction_meas[pair] = vec_proton_from_lambda_direction_between_trackers_meas[pair-1].Cross(vec_proton_from_lambda_direction_between_trackers_meas[pair]).Unit();
      }
    }
    // -----


    // phi angle between lambda and proton spins with two trackers
    Double_t proton_phi_of_spins[5];
    Double_t proton_from_lambda_phi_of_spins[5];
    Double_t proton_phi_of_spins_meas[5];
    Double_t proton_from_lambda_phi_of_spins_meas[5];
    for(int pair=0; pair<5; ++pair){
      proton_phi_of_spins[pair] = acos(vec_lambda_expected_spin_direction_perp_to_proton[pair].Dot(vec_proton_expected_spin_direction[pair]));
      proton_from_lambda_phi_of_spins[pair] = acos(vec_lambda_expected_spin_direction_perp_to_proton_from_lambda[pair].Dot(vec_proton_from_lambda_expected_spin_direction[pair]));
      proton_phi_of_spins_meas[pair] = acos(vec_lambda_expected_spin_direction_perp_to_proton_meas[pair].Dot(vec_proton_expected_spin_direction_meas[pair]));
      proton_from_lambda_phi_of_spins_meas[pair] = acos(vec_lambda_expected_spin_direction_perp_to_proton_from_lambda_meas[pair].Dot(vec_proton_from_lambda_expected_spin_direction_meas[pair]));
      Double_t proton_phi_sign;
      Double_t proton_from_lambda_phi_sign;
      Double_t proton_phi_sign_meas;
      Double_t proton_from_lambda_phi_sign_meas;
      if(pair==0){
        proton_phi_sign = vec_lambda_expected_spin_direction_perp_to_proton[pair].Dot(vec_proton_direction_at_cdc.Cross(vec_proton_expected_spin_direction[pair]));
        proton_from_lambda_phi_sign = vec_lambda_expected_spin_direction_perp_to_proton_from_lambda[pair].Dot(vec_proton_from_lambda_direction_at_cdc.Cross(vec_proton_from_lambda_expected_spin_direction[pair]));
        proton_phi_sign_meas = vec_lambda_expected_spin_direction_perp_to_proton_meas[pair].Dot(vec_proton_direction_at_cdc_meas.Cross(vec_proton_expected_spin_direction_meas[pair]));
        proton_from_lambda_phi_sign_meas = vec_lambda_expected_spin_direction_perp_to_proton_from_lambda_meas[pair].Dot(vec_proton_from_lambda_direction_at_cdc_meas.Cross(vec_proton_from_lambda_expected_spin_direction_meas[pair]));
      }
      else{
        proton_phi_sign = vec_lambda_expected_spin_direction_perp_to_proton[pair].Dot(vec_proton_direction_between_trackers[pair-1].Cross(vec_proton_expected_spin_direction[pair]));
        proton_from_lambda_phi_sign = vec_lambda_expected_spin_direction_perp_to_proton_from_lambda[pair].Dot(vec_proton_from_lambda_direction_between_trackers[pair-1].Cross(vec_proton_from_lambda_expected_spin_direction[pair]));
        proton_phi_sign_meas = vec_lambda_expected_spin_direction_perp_to_proton_meas[pair].Dot(vec_proton_direction_between_trackers_meas[pair-1].Cross(vec_proton_expected_spin_direction_meas[pair]));
        proton_from_lambda_phi_sign_meas = vec_lambda_expected_spin_direction_perp_to_proton_from_lambda_meas[pair].Dot(vec_proton_from_lambda_direction_between_trackers_meas[pair-1].Cross(vec_proton_from_lambda_expected_spin_direction_meas[pair]));
      }

      proton_phi_of_spins[pair] *= proton_phi_sign/fabs(proton_phi_sign);
      proton_from_lambda_phi_of_spins[pair] *= proton_from_lambda_phi_sign/fabs(proton_from_lambda_phi_sign);
      proton_phi_of_spins_meas[pair] *= proton_phi_sign_meas/fabs(proton_phi_sign_meas);
      proton_from_lambda_phi_of_spins_meas[pair] *= proton_from_lambda_phi_sign_meas/fabs(proton_from_lambda_phi_sign_meas);
    }

    // ====================================================================================================
    // event selections
    // ====================================================================================================
    bool are_all_detected_by_cdc = false;
    if(is_proton_detected_by_cdc&&is_proton_from_lambda_detected_by_cdc&&is_pim_from_lambda_detected_by_cdc){
      are_all_detected_by_cdc = true;
    }
    bool are_momenta_over_threshold = false;
    Double_t proton_momentum_ll = 0.1;
    Double_t pim_momentum_ll = 0.5;
    if(lvec_proton_mc.P()>proton_momentum_ll&&lvec_proton_from_lambda_mc.P()>proton_momentum_ll&&lvec_pim_from_lambda_mc.P()>pim_momentum_ll){
      are_momenta_over_threshold = true;
    }

    bool is_proton_detected_by_both_trackers[5] = { false, false, false, false, false };
    bool is_proton_asymmetric_scattering_in_nc[5] = { false, false, false, false, false };
    bool is_proton_from_lambda_detected_by_both_trackers[5] = { false, false, false, false, false };
    bool is_proton_from_lambda_asymmetric_scattering_in_nc[5] = { false, false, false, false, false };
    Double_t theta_sel_ll = 6.;
    Double_t theta_sel_ul = 30.;
    Double_t theta_meas_sel_ll = 8.;
    Double_t theta_meas_sel_ul = 30.;
    bool is_proton_scattering_angle_in_region[5] = { false, false, false, false, false };
    bool is_proton_scattering_angle_meas_in_region[5] = { false, false, false, false, false };
    bool is_proton_from_lambda_scattering_angle_in_region[5] = { false, false, false, false, false };
    bool is_proton_from_lambda_scattering_angle_meas_in_region[5] = { false, false, false, false, false };
    for(int pair=0; pair<5; ++pair){
      int first_layer  = pair*2;
      int second_layer = first_layer+1;
      if(is_proton_detected_by_tracker[first_layer]&&is_proton_detected_by_tracker[second_layer]){
        is_proton_detected_by_both_trackers[pair] = true;
      }
      if(is_proton_asymmetric_scattering[first_layer]){
        is_proton_asymmetric_scattering_in_nc[pair] = true;
      }
      if(is_proton_from_lambda_detected_by_tracker[first_layer]&&is_proton_from_lambda_detected_by_tracker[second_layer]){
        is_proton_from_lambda_detected_by_both_trackers[pair] = true;
      }
      if(is_proton_from_lambda_asymmetric_scattering[first_layer]){
        is_proton_from_lambda_asymmetric_scattering_in_nc[pair] = true;
      }
      if(theta_sel_ll<proton_scattering_angle_theta[pair]&&proton_scattering_angle_theta[pair]<theta_sel_ul){
        is_proton_scattering_angle_in_region[pair] = true;
      }
      if(theta_meas_sel_ll<proton_scattering_angle_theta_meas[pair]&&proton_scattering_angle_theta_meas[pair]<theta_meas_sel_ul){
        is_proton_scattering_angle_meas_in_region[pair] = true;
      }
      if(theta_sel_ll<proton_from_lambda_scattering_angle_theta[pair]&&proton_from_lambda_scattering_angle_theta[pair]<theta_sel_ul){
        is_proton_from_lambda_scattering_angle_in_region[pair] = true;
      }
      if(theta_meas_sel_ll<proton_from_lambda_scattering_angle_theta_meas[pair]&&proton_from_lambda_scattering_angle_theta_meas[pair]<theta_meas_sel_ul){
        is_proton_from_lambda_scattering_angle_meas_in_region[pair] = true;
      }
    }

    // for lp_mass & lp_momentum selection
    bool is_lp_mass_in_kpp_region = false;
    if(2.1<lvec_lp_mc.M()/1000.&&lvec_lp_mc.M()/1000.<2.6){
      is_lp_mass_in_kpp_region = true;
    }
    bool is_lp_momentum_in_kpp_region = false;
    if(0.25<lvec_lp_mc.P()&&lvec_lp_mc.P()/1000.<0.8){
      is_lp_momentum_in_kpp_region = true;
    }

    // ====================================================================================================

    TString sel_name[84] = {
      /*  0 */ "",
      /*  1 */ "_charged",
      /*  2 */ "_charged_mom",
      /*  3 */ "_charged_mom_cdc",

      /*  4 */ "_charged_mom_cdc_proton_pair1",
      /*  5 */ "_charged_mom_cdc_proton_pair1_asym",
      /*  6 */ "_charged_mom_cdc_proton_pair1_sel",
      /*  7 */ "_charged_mom_cdc_proton_pair1_sel_asym",
      /*  8 */ "_charged_mom_cdc_proton_pair1_sel_mqsel",
      /*  9 */ "_charged_mom_cdc_proton_pair1_selmeas",
      /* 10 */ "_charged_mom_cdc_proton_pair1_selmeas_asym",
      /* 11 */ "_charged_mom_cdc_proton_pair1_selmeas_mqsel",
      /* 12 */ "_charged_mom_cdc_proton_pair2",
      /* 13 */ "_charged_mom_cdc_proton_pair2_asym",
      /* 14 */ "_charged_mom_cdc_proton_pair2_sel",
      /* 15 */ "_charged_mom_cdc_proton_pair2_sel_asym",
      /* 16 */ "_charged_mom_cdc_proton_pair2_sel_mqsel",
      /* 17 */ "_charged_mom_cdc_proton_pair2_selmeas",
      /* 18 */ "_charged_mom_cdc_proton_pair2_selmeas_asym",
      /* 19 */ "_charged_mom_cdc_proton_pair2_selmeas_mqsel",
      /* 20 */ "_charged_mom_cdc_proton_pair3",
      /* 21 */ "_charged_mom_cdc_proton_pair3_asym",
      /* 22 */ "_charged_mom_cdc_proton_pair3_sel",
      /* 23 */ "_charged_mom_cdc_proton_pair3_sel_asym",
      /* 24 */ "_charged_mom_cdc_proton_pair3_sel_mqsel",
      /* 25 */ "_charged_mom_cdc_proton_pair3_selmeas",
      /* 26 */ "_charged_mom_cdc_proton_pair3_selmeas_asym",
      /* 27 */ "_charged_mom_cdc_proton_pair3_selmeas_mqsel",
      /* 28 */ "_charged_mom_cdc_proton_pair4",
      /* 29 */ "_charged_mom_cdc_proton_pair4_asym",
      /* 30 */ "_charged_mom_cdc_proton_pair4_sel",
      /* 31 */ "_charged_mom_cdc_proton_pair4_sel_asym",
      /* 32 */ "_charged_mom_cdc_proton_pair4_sel_mqsel",
      /* 33 */ "_charged_mom_cdc_proton_pair4_selmeas",
      /* 34 */ "_charged_mom_cdc_proton_pair4_selmeas_asym",
      /* 35 */ "_charged_mom_cdc_proton_pair4_selmeas_mqsel",
      /* 36 */ "_charged_mom_cdc_proton_pair5",
      /* 37 */ "_charged_mom_cdc_proton_pair5_asym",
      /* 38 */ "_charged_mom_cdc_proton_pair5_sel",
      /* 39 */ "_charged_mom_cdc_proton_pair5_sel_asym",
      /* 40 */ "_charged_mom_cdc_proton_pair5_sel_mqsel",
      /* 41 */ "_charged_mom_cdc_proton_pair5_selmeas",
      /* 42 */ "_charged_mom_cdc_proton_pair5_selmeas_asym",
      /* 43 */ "_charged_mom_cdc_proton_pair5_selmeas_mqsel",

      /* 44 */ "_charged_mom_cdc_proton_from_lambda_pair1",
      /* 45 */ "_charged_mom_cdc_proton_from_lambda_pair1_asym",
      /* 46 */ "_charged_mom_cdc_proton_from_lambda_pair1_sel",
      /* 47 */ "_charged_mom_cdc_proton_from_lambda_pair1_sel_asym",
      /* 48 */ "_charged_mom_cdc_proton_from_lambda_pair1_sel_mqsel",
      /* 49 */ "_charged_mom_cdc_proton_from_lambda_pair1_selmeas",
      /* 50 */ "_charged_mom_cdc_proton_from_lambda_pair1_selmeas_asym",
      /* 51 */ "_charged_mom_cdc_proton_from_lambda_pair1_selmeas_mqsel",
      /* 52 */ "_charged_mom_cdc_proton_from_lambda_pair2",
      /* 53 */ "_charged_mom_cdc_proton_from_lambda_pair2_asym",
      /* 54 */ "_charged_mom_cdc_proton_from_lambda_pair2_sel",
      /* 55 */ "_charged_mom_cdc_proton_from_lambda_pair2_sel_asym",
      /* 56 */ "_charged_mom_cdc_proton_from_lambda_pair2_sel_mqsel",
      /* 57 */ "_charged_mom_cdc_proton_from_lambda_pair2_selmeas",
      /* 58 */ "_charged_mom_cdc_proton_from_lambda_pair2_selmeas_asym",
      /* 59 */ "_charged_mom_cdc_proton_from_lambda_pair2_selmeas_mqsel",
      /* 60 */ "_charged_mom_cdc_proton_from_lambda_pair3",
      /* 61 */ "_charged_mom_cdc_proton_from_lambda_pair3_asym",
      /* 62 */ "_charged_mom_cdc_proton_from_lambda_pair3_sel",
      /* 63 */ "_charged_mom_cdc_proton_from_lambda_pair3_sel_asym",
      /* 64 */ "_charged_mom_cdc_proton_from_lambda_pair3_sel_mqsel",
      /* 65 */ "_charged_mom_cdc_proton_from_lambda_pair3_selmeas",
      /* 66 */ "_charged_mom_cdc_proton_from_lambda_pair3_selmeas_asym",
      /* 67 */ "_charged_mom_cdc_proton_from_lambda_pair3_selmeas_mqsel",
      /* 68 */ "_charged_mom_cdc_proton_from_lambda_pair4",
      /* 69 */ "_charged_mom_cdc_proton_from_lambda_pair4_asym",
      /* 70 */ "_charged_mom_cdc_proton_from_lambda_pair4_sel",
      /* 71 */ "_charged_mom_cdc_proton_from_lambda_pair4_sel_asym",
      /* 72 */ "_charged_mom_cdc_proton_from_lambda_pair4_sel_mqsel",
      /* 73 */ "_charged_mom_cdc_proton_from_lambda_pair4_selmeas",
      /* 74 */ "_charged_mom_cdc_proton_from_lambda_pair4_selmeas_asym",
      /* 75 */ "_charged_mom_cdc_proton_from_lambda_pair4_selmeas_mqsel",
      /* 76 */ "_charged_mom_cdc_proton_from_lambda_pair5",
      /* 77 */ "_charged_mom_cdc_proton_from_lambda_pair5_asym",
      /* 78 */ "_charged_mom_cdc_proton_from_lambda_pair5_sel",
      /* 79 */ "_charged_mom_cdc_proton_from_lambda_pair5_sel_asym",
      /* 80 */ "_charged_mom_cdc_proton_from_lambda_pair5_sel_mqsel",
      /* 81 */ "_charged_mom_cdc_proton_from_lambda_pair5_selmeas",
      /* 82 */ "_charged_mom_cdc_proton_from_lambda_pair5_selmeas_asym",
      /* 83 */ "_charged_mom_cdc_proton_from_lambda_pair5_selmeas_mqsel",
    };
    bool sel_flag[84];
    sel_flag[0] = true;
    sel_flag[1] = is_lambda_charged_decay;
    sel_flag[2] = is_lambda_charged_decay & are_momenta_over_threshold;
    sel_flag[3] = is_lambda_charged_decay & are_momenta_over_threshold & are_all_detected_by_cdc;
    for(int pair=0; pair<5; ++pair){
      sel_flag[4 + pair*8 + 0] = sel_flag[3] & is_proton_detected_by_both_trackers[pair];
      sel_flag[4 + pair*8 + 1] = sel_flag[4 + pair*8 + 0] & is_proton_asymmetric_scattering_in_nc[pair];
      sel_flag[4 + pair*8 + 2] = sel_flag[4 + pair*8 + 0] & is_proton_scattering_angle_in_region[pair];
      sel_flag[4 + pair*8 + 3] = sel_flag[4 + pair*8 + 0] & is_proton_scattering_angle_in_region[pair] & is_proton_asymmetric_scattering_in_nc[pair];
      sel_flag[4 + pair*8 + 4] = sel_flag[4 + pair*8 + 0] & is_proton_scattering_angle_in_region[pair] & is_lp_mass_in_kpp_region & is_lp_momentum_in_kpp_region;
      sel_flag[4 + pair*8 + 5] = sel_flag[4 + pair*8 + 0] & is_proton_scattering_angle_meas_in_region[pair];
      sel_flag[4 + pair*8 + 6] = sel_flag[4 + pair*8 + 0] & is_proton_scattering_angle_meas_in_region[pair] & is_proton_asymmetric_scattering_in_nc[pair];
      sel_flag[4 + pair*8 + 7] = sel_flag[4 + pair*8 + 0] & is_proton_scattering_angle_meas_in_region[pair] & is_lp_mass_in_kpp_region & is_lp_momentum_in_kpp_region;
      if(pair!=0){
        sel_flag[4 + pair*8 + 1] = sel_flag[4 + pair*8 + 1] & !is_proton_asymmetric_scattering_in_nc[pair-1];
        sel_flag[4 + pair*8 + 3] = sel_flag[4 + pair*8 + 3] & !is_proton_asymmetric_scattering_in_nc[pair-1];
        sel_flag[4 + pair*8 + 6] = sel_flag[4 + pair*8 + 6] & !is_proton_asymmetric_scattering_in_nc[pair-1];
      }
    }
    for(int pair=0; pair<5; ++pair){
      sel_flag[44 + pair*8 + 0] = sel_flag[3] & is_proton_from_lambda_detected_by_both_trackers[pair];
      sel_flag[44 + pair*8 + 1] = sel_flag[44 + pair*8 + 0] & is_proton_from_lambda_asymmetric_scattering_in_nc[pair];
      sel_flag[44 + pair*8 + 2] = sel_flag[44 + pair*8 + 0] & is_proton_from_lambda_scattering_angle_in_region[pair];
      sel_flag[44 + pair*8 + 3] = sel_flag[44 + pair*8 + 0] & is_proton_from_lambda_scattering_angle_in_region[pair] & is_proton_from_lambda_asymmetric_scattering_in_nc[pair];
      sel_flag[44 + pair*8 + 4] = sel_flag[44 + pair*8 + 0] & is_proton_from_lambda_scattering_angle_in_region[pair] & is_lp_mass_in_kpp_region & is_lp_momentum_in_kpp_region;
      sel_flag[44 + pair*8 + 5] = sel_flag[44 + pair*8 + 0] & is_proton_from_lambda_scattering_angle_meas_in_region[pair];
      sel_flag[44 + pair*8 + 6] = sel_flag[44 + pair*8 + 0] & is_proton_from_lambda_scattering_angle_meas_in_region[pair] & is_proton_from_lambda_asymmetric_scattering_in_nc[pair];
      sel_flag[44 + pair*8 + 7] = sel_flag[44 + pair*8 + 0] & is_proton_from_lambda_scattering_angle_meas_in_region[pair] & is_lp_mass_in_kpp_region & is_lp_momentum_in_kpp_region;
      if(pair!=0){
        sel_flag[44+ pair*8 + 1] = sel_flag[44+ pair*8 + 1] & !is_proton_from_lambda_asymmetric_scattering_in_nc[pair-1];
        sel_flag[44+ pair*8 + 3] = sel_flag[44+ pair*8 + 3] & !is_proton_from_lambda_asymmetric_scattering_in_nc[pair-1];
        sel_flag[44+ pair*8 + 6] = sel_flag[44+ pair*8 + 6] & !is_proton_from_lambda_asymmetric_scattering_in_nc[pair-1];
      }
    }


    // Fill histograms
    for(int i_sel = 0; i_sel<84; ++i_sel){
      if(sel_flag[i_sel]){
        ((TH1F*)outfile->Get(Form("lambda_momentum%s",sel_name[i_sel].Data())))->Fill(lvec_lambda_mc.P()/1000.);
        ((TH1F*)outfile->Get(Form("proton_momentum%s",sel_name[i_sel].Data())))->Fill(lvec_proton_mc.P()/1001.);
        ((TH1F*)outfile->Get(Form("proton_from_lambda_momentum%s",sel_name[i_sel].Data())))->Fill(lvec_proton_from_lambda_mc.P()/1000.);
        ((TH1F*)outfile->Get(Form("pim_from_lambda_momentum%s",sel_name[i_sel].Data())))->Fill(lvec_pim_from_lambda_mc.P()/1000.);
        ((TH1F*)outfile->Get(Form("lp_mass%s",sel_name[i_sel].Data())))->Fill(lvec_lp_mc.M()/1000.);
        ((TH1F*)outfile->Get(Form("lp_momentum%s",sel_name[i_sel].Data())))->Fill(lvec_lp_mc.P()/1000.);

        ((TH2F*)outfile->Get(Form("lp_mass_vs_lp_momentum%s",sel_name[i_sel].Data())))->Fill(lvec_lp_mc.M()/1000.,lvec_lp_mc.P()/1000.);

        ((TH1F*)outfile->Get(Form("ppim_mass%s",sel_name[i_sel].Data())))->Fill(lvec_ppim_mc.M()/1000.);
        ((TH1F*)outfile->Get(Form("ppim_momentum%s",sel_name[i_sel].Data())))->Fill(lvec_ppim_mc.P()/1000.);

        ((TH1F*)outfile->Get(Form("proton_momentum_at_cdc%s",sel_name[i_sel].Data())))->Fill(vec_proton_momentum_at_cdc.Mag()/1000.);
        ((TH1F*)outfile->Get(Form("proton_from_lambda_momentum_at_cdc%s",sel_name[i_sel].Data())))->Fill(vec_proton_from_lambda_momentum_at_cdc.Mag()/1000.);

        for(int layer=0; layer<10; ++layer){
          ((TH1F*)outfile->Get(Form("proton_momentum_at_tracker_layer%d%s",layer+1,sel_name[i_sel].Data())))->Fill(vec_proton_momentum_at_tracker[layer].Mag()/1000.);
          ((TH1F*)outfile->Get(Form("proton_from_lambda_momentum_at_tracker_layer%d%s",layer+1,sel_name[i_sel].Data())))->Fill(vec_proton_from_lambda_momentum_at_tracker[layer].Mag()/1000.);
        }

        for(int pair=0; pair<5; ++pair){
          if(pair==0){
            ((TH1F*)outfile->Get(Form("proton_momentum_difference_between_cdc_and_tracker_layer1%s",sel_name[i_sel].Data())))->Fill(proton_momentum_difference[pair]/1000.);
            ((TH1F*)outfile->Get(Form("proton_from_lambda_momentum_difference_between_cdc_and_tracker_layer1%s",sel_name[i_sel].Data())))->Fill(proton_from_lambda_momentum_difference[pair]/1000.);
          }
          else{
            int first_layer = pair*2;
            int second_layer = first_layer+1;
            ((TH1F*)outfile->Get(Form("proton_momentum_difference_between_tracker_layer%d_and_layer%d%s",first_layer,second_layer,sel_name[i_sel].Data())))->Fill(proton_momentum_difference[pair]/1000.);
            ((TH1F*)outfile->Get(Form("proton_from_lambda_momentum_difference_between_tracker_layer%d_and_layer%d%s",first_layer,second_layer,sel_name[i_sel].Data())))->Fill(proton_from_lambda_momentum_difference[pair]/1000.);
          }
        }

        for(int pair=0; pair<5; ++pair){
          ((TH1F*)outfile->Get(Form("proton_scattering_angle_theta_mc_pair%d%s",pair+1,sel_name[i_sel].Data())))->Fill(proton_scattering_angle_theta_mc[pair]);
          ((TH1F*)outfile->Get(Form("proton_scattering_angle_theta_pair%d%s",pair+1,sel_name[i_sel].Data())))->Fill(proton_scattering_angle_theta[pair]);
          ((TH1F*)outfile->Get(Form("proton_scattering_angle_theta_meas_pair%d%s",pair+1,sel_name[i_sel].Data())))->Fill(proton_scattering_angle_theta_meas[pair]);
          ((TH1F*)outfile->Get(Form("proton_phi_of_spins_pair%d%s",pair+1,sel_name[i_sel].Data())))->Fill(proton_phi_of_spins[pair]/TMath::Pi());
          ((TH1F*)outfile->Get(Form("proton_phi_of_spins_meas_pair%d%s",pair+1,sel_name[i_sel].Data())))->Fill(proton_phi_of_spins_meas[pair]/TMath::Pi());
          ((TH2F*)outfile->Get(Form("proton_phi_of_spins_pair%d_vs_scattering_angle_theta_pair%d%s",pair+1,pair+1,sel_name[i_sel].Data())))->Fill(proton_phi_of_spins[pair]/TMath::Pi(),proton_scattering_angle_theta[pair]);
          ((TH2F*)outfile->Get(Form("proton_phi_of_spins_meas_pair%d_vs_scattering_angle_theta_meas_pair%d%s",pair+1,pair+1,sel_name[i_sel].Data())))->Fill(proton_phi_of_spins_meas[pair]/TMath::Pi(),proton_scattering_angle_theta_meas[pair]);

          ((TH1F*)outfile->Get(Form("proton_from_lambda_scattering_angle_theta_mc_pair%d%s",pair+1,sel_name[i_sel].Data())))->Fill(proton_from_lambda_scattering_angle_theta_mc[pair]);
          ((TH1F*)outfile->Get(Form("proton_from_lambda_scattering_angle_theta_pair%d%s",pair+1,sel_name[i_sel].Data())))->Fill(proton_from_lambda_scattering_angle_theta[pair]);
          ((TH1F*)outfile->Get(Form("proton_from_lambda_scattering_angle_theta_meas_pair%d%s",pair+1,sel_name[i_sel].Data())))->Fill(proton_from_lambda_scattering_angle_theta_meas[pair]);
          ((TH1F*)outfile->Get(Form("proton_from_lambda_phi_of_spins_pair%d%s",pair+1,sel_name[i_sel].Data())))->Fill(proton_from_lambda_phi_of_spins[pair]/TMath::Pi());
          ((TH1F*)outfile->Get(Form("proton_from_lambda_phi_of_spins_meas_pair%d%s",pair+1,sel_name[i_sel].Data())))->Fill(proton_from_lambda_phi_of_spins_meas[pair]/TMath::Pi());
          ((TH2F*)outfile->Get(Form("proton_from_lambda_phi_of_spins_pair%d_vs_scattering_angle_theta_pair%d%s",pair+1,pair+1,sel_name[i_sel].Data())))->Fill(proton_from_lambda_phi_of_spins[pair]/TMath::Pi(),proton_from_lambda_scattering_angle_theta[pair]);
          ((TH2F*)outfile->Get(Form("proton_from_lambda_phi_of_spins_meas_pair%d_vs_scattering_angle_theta_meas_pair%d%s",pair+1,pair+1,sel_name[i_sel].Data())))->Fill(proton_from_lambda_phi_of_spins_meas[pair]/TMath::Pi(),proton_from_lambda_scattering_angle_theta_meas[pair]);
        }


        ((TH1F*)outfile->Get(Form("proton_diff_phi_direction_at_cdc%s",sel_name[i_sel].Data())))->Fill(proton_diff_phi_direction_at_cdc/TMath::Pi()*180.);
        ((TH1F*)outfile->Get(Form("proton_diff_theta_direction_at_cdc%s",sel_name[i_sel].Data())))->Fill(proton_diff_theta_direction_at_cdc/TMath::Pi()*180.);
        ((TH1F*)outfile->Get(Form("proton_diff_phirho_position_at_cdc%s",sel_name[i_sel].Data())))->Fill(proton_diff_phirho_position_at_cdc);
        ((TH1F*)outfile->Get(Form("proton_diff_z_position_at_cdc%s",sel_name[i_sel].Data())))->Fill(proton_diff_z_position_at_cdc);
        ((TH2F*)outfile->Get(Form("proton_diff_theta_direction_at_cdc_vs_diff_z_position_at_cdc%s",sel_name[i_sel].Data())))->Fill(proton_diff_theta_direction_at_cdc/TMath::Pi()*180.,proton_diff_z_position_at_cdc);
        ((TH2F*)outfile->Get(Form("proton_diff_phi_direction_at_cdc_vs_diff_phirho_position_at_cdc%s",sel_name[i_sel].Data())))->Fill(proton_diff_phi_direction_at_cdc/TMath::Pi()*180.,proton_diff_phirho_position_at_cdc);
        for(int layer=0; layer<10; ++layer){
          ((TH1F*)outfile->Get(Form("proton_diff_phirho_position_at_tracker_layer%d%s",layer+1,sel_name[i_sel].Data())))->Fill(proton_diff_phirho_position_at_tracker[layer]);
          ((TH1F*)outfile->Get(Form("proton_diff_z_position_at_tracker_layer%d%s",layer+1,sel_name[i_sel].Data())))->Fill(proton_diff_z_position_at_tracker[layer]);
        }
      }
    }

  } // end of event loop
}
