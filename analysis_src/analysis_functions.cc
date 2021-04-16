// analysis_functions.cc

#include"analysis_functions.hh"
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
Int_t number_of_hits_in_tracker_layer1;
Int_t number_of_hits_in_tracker_layer2;
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
TTree* tree_tracker_layer1_hit;
Int_t tracker_layer1_event_id;
Int_t tracker_layer1_track_id;
Int_t tracker_layer1_parent_id;
Int_t tracker_layer1_particle_id;
Int_t tracker_layer1_layer_id;
Float_t tracker_layer1_hit_time;
Float_t tracker_layer1_hit_position_x;
Float_t tracker_layer1_hit_position_y;
Float_t tracker_layer1_hit_position_z;
Float_t tracker_layer1_momentum_x;
Float_t tracker_layer1_momentum_y;
Float_t tracker_layer1_momentum_z;
Int_t tracker_layer1_is_asymmetric_scattering;
// -----
TTree* tree_tracker_layer2_hit;
Int_t tracker_layer2_event_id;
Int_t tracker_layer2_track_id;
Int_t tracker_layer2_parent_id;
Int_t tracker_layer2_particle_id;
Int_t tracker_layer2_layer_id;
Float_t tracker_layer2_hit_time;
Float_t tracker_layer2_hit_position_x;
Float_t tracker_layer2_hit_position_y;
Float_t tracker_layer2_hit_position_z;
Float_t tracker_layer2_momentum_x;
Float_t tracker_layer2_momentum_y;
Float_t tracker_layer2_momentum_z;
Int_t tracker_layer2_is_asymmetric_scattering;
// -----

void DrawHistograms(TFile* outfile, std::string pdf_name){
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

  TString sel_name[10] = {
    "",
    "_charged",
    "_charged_mom",
    "_charged_mom_cdc",
    "_charged_mom_cdc_both",
    "_charged_mom_cdc_both_asym",
    "_charged_mom_cdc_one",
    "_charged_mom_cdc_one_asym",
    "_charged_mom_cdc_lay1",
    "_charged_mom_cdc_lay1_asym"
  };

  for(int i_sel=0; i_sel<10; ++i_sel){
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

    // page-1
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
    hist_1d = (TH1F*)outfile->Get(Form("scattering_angle_theta_mc%s",sel_name[i_sel].Data()));
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

    // page-4
    TF1* f_phi = new TF1("f_phi","1 + [0]*cos(x)",-TMath::Pi(),TMath::Pi());
    // scattering_angle_theta
    canvas->cd(1);
    gPad->SetLogy();
    hist_1d = (TH1F*)outfile->Get(Form("scattering_angle_theta%s",sel_name[i_sel].Data()));
    Draw(hist_1d);

    // phi_of_spins
    canvas->cd(2);
    hist_1d = (TH1F*)outfile->Get(Form("phi_of_spins%s",sel_name[i_sel].Data()));
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

    // scattering_angle_theta
    canvas->cd(4);
    gPad->SetLogy();
    hist_1d = (TH1F*)outfile->Get(Form("scattering_angle_theta_rough%s",sel_name[i_sel].Data()));
    Draw(hist_1d);

    // phi_of_spins
    canvas->cd(5);
    hist_1d = (TH1F*)outfile->Get(Form("phi_of_spins_rough%s",sel_name[i_sel].Data()));
    Draw(hist_1d);
    hist_1d->SetMinimum(0);

    // phi_of_spins (scaled & fitting)
    canvas->cd(6);
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

  TString sel_name[10] = {
    "",
    "_charged",
    "_charged_mom",
    "_charged_mom_cdc",
    "_charged_mom_cdc_both",
    "_charged_mom_cdc_both_asym",
    "_charged_mom_cdc_one",
    "_charged_mom_cdc_one_asym",
    "_charged_mom_cdc_lay1",
    "_charged_mom_cdc_lay1_asym"
  };

  for(int i_sel=0; i_sel<10; ++i_sel){
    // 1D histograms
    new TH1F(Form("lambda_momentum%s",sel_name[i_sel].Data()),"lambda momentum;lambda momentum (GeV/c);counts",200,0.,2.);
    new TH1F(Form("proton_momentum%s",sel_name[i_sel].Data()),"proton momentum;proton momentum (GeV/c);counts",200,0.,2.);
    new TH1F(Form("proton_from_lambda_momentum%s",sel_name[i_sel].Data()),"proton from lambda momentum;proton from #Lambda momentum (GeV/c);counts",200,0.,2.);
    new TH1F(Form("pim_from_lambda_momentum%s",sel_name[i_sel].Data()),"pim from lambda momentum;#pi^{#minus} from #Lambda momentum (GeV/c);counts",200,0.,2.);
    new TH1F(Form("lp_mass%s",sel_name[i_sel].Data()),"lp mass;#Lambdap mass (GeV/c^{2});counts",100,2.,3.);
    new TH1F(Form("lp_momentum%s",sel_name[i_sel].Data()),"lp momentum;#Lambdap momentum (GeV/c);counts",200,0.,2.);
    new TH1F(Form("ppim_mass%s",sel_name[i_sel].Data()),"ppim mass;p#pi^{#minus} mass (GeV/c^{2});counts",100,1.,2.);
    new TH1F(Form("ppim_momentum%s",sel_name[i_sel].Data()),"ppim momentum;p#pi^{#minus} momentum (GeV/c);counts",200,0.,2.);

    new TH1F(Form("proton_momentum_at_cdc%s",sel_name[i_sel].Data()),"proton momentum at cdc;proton momentum (GeV/c);counts",200,0.,2.);
    new TH1F(Form("proton_momentum_at_tracker_layer1%s",sel_name[i_sel].Data()),"proton momentum at tracker_layer1;proton momentum (GeV/c);counts",200,0.,2.);
    new TH1F(Form("proton_momentum_at_tracker_layer2%s",sel_name[i_sel].Data()),"proton momentum at tracker_layer2;proton momentum (GeV/c);counts",200,0.,2.);
    new TH1F(Form("proton_momentum_difference_between_cdc_and_tracker_layer1%s",sel_name[i_sel].Data()),"proton momentum difference between cdc and tracker layer1;proton momentum difference(GeV/c);counts",200,0.,0.5);

    new TH1F(Form("scattering_angle_theta_mc%s",sel_name[i_sel].Data()),"scattering angle theta mc value;#theta (deg.);counts",200,0.,50.);

    new TH1F(Form("scattering_angle_theta%s",sel_name[i_sel].Data()),"scattering angle theta;#theta (deg.);counts",200,0.,50.);
    new TH1F(Form("phi_of_spins%s",sel_name[i_sel].Data()),"phi between spins;#phi (rad.);counts",10,-TMath::Pi(),TMath::Pi());

    new TH1F(Form("scattering_angle_theta_rough%s",sel_name[i_sel].Data()),"scattering angle theta;#theta (deg.);counts",200,0.,50.);
    new TH1F(Form("phi_of_spins_rough%s",sel_name[i_sel].Data()),"phi between spins;#phi (rad.);counts",10,-TMath::Pi(),TMath::Pi());

    // 2D histograms
    new TH2F(Form("lp_mass_vs_lp_momentum%s",sel_name[i_sel].Data()),"lp mass vs lp momentum;#Lambdap mass (GeV/c^{2});#Lambdap momentum;counts",100,2.,3.,100,0.,1.);
  }
}

void InitializeTrees(TFile* infile){
  infile->cd();

  tree_event = (TTree*)infile->Get("tree_event");     
  tree_event->SetBranchAddress("event_id",&event_id);
  tree_event->SetBranchAddress("number_of_primaries",&number_of_primaries);
  tree_event->SetBranchAddress("number_of_trajectories",&number_of_trajectories);
  tree_event->SetBranchAddress("number_of_hits_in_cdc",&number_of_hits_in_cdc);
  tree_event->SetBranchAddress("number_of_hits_in_tracker_layer1",&number_of_hits_in_tracker_layer1);
  tree_event->SetBranchAddress("number_of_hits_in_tracker_layer2",&number_of_hits_in_tracker_layer2);
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
  tree_tracker_layer1_hit = (TTree*)infile->Get("tree_tracker_layer1_hit");     
  tree_tracker_layer1_hit->SetBranchAddress("event_id",&tracker_layer1_event_id);
  tree_tracker_layer1_hit->SetBranchAddress("track_id",&tracker_layer1_track_id);
  tree_tracker_layer1_hit->SetBranchAddress("parent_id",&tracker_layer1_parent_id);
  tree_tracker_layer1_hit->SetBranchAddress("particle_id",&tracker_layer1_particle_id);
  tree_tracker_layer1_hit->SetBranchAddress("layer_id",&tracker_layer1_layer_id);
  tree_tracker_layer1_hit->SetBranchAddress("hit_time",&tracker_layer1_hit_time);
  tree_tracker_layer1_hit->SetBranchAddress("hit_position_x",&tracker_layer1_hit_position_x);
  tree_tracker_layer1_hit->SetBranchAddress("hit_position_y",&tracker_layer1_hit_position_y);
  tree_tracker_layer1_hit->SetBranchAddress("hit_position_z",&tracker_layer1_hit_position_z);
  tree_tracker_layer1_hit->SetBranchAddress("momentum_x",&tracker_layer1_momentum_x);
  tree_tracker_layer1_hit->SetBranchAddress("momentum_y",&tracker_layer1_momentum_y);
  tree_tracker_layer1_hit->SetBranchAddress("momentum_z",&tracker_layer1_momentum_z);
  tree_tracker_layer1_hit->SetBranchAddress("is_asymmetric_scattering",&tracker_layer1_is_asymmetric_scattering);
  // -----
  tree_tracker_layer2_hit = (TTree*)infile->Get("tree_tracker_layer2_hit");     
  tree_tracker_layer2_hit->SetBranchAddress("event_id",&tracker_layer2_event_id);
  tree_tracker_layer2_hit->SetBranchAddress("track_id",&tracker_layer2_track_id);
  tree_tracker_layer2_hit->SetBranchAddress("parent_id",&tracker_layer2_parent_id);
  tree_tracker_layer2_hit->SetBranchAddress("particle_id",&tracker_layer2_particle_id);
  tree_tracker_layer2_hit->SetBranchAddress("layer_id",&tracker_layer2_layer_id);
  tree_tracker_layer2_hit->SetBranchAddress("hit_time",&tracker_layer2_hit_time);
  tree_tracker_layer2_hit->SetBranchAddress("hit_position_x",&tracker_layer2_hit_position_x);
  tree_tracker_layer2_hit->SetBranchAddress("hit_position_y",&tracker_layer2_hit_position_y);
  tree_tracker_layer2_hit->SetBranchAddress("hit_position_z",&tracker_layer2_hit_position_z);
  tree_tracker_layer2_hit->SetBranchAddress("momentum_x",&tracker_layer2_momentum_x);
  tree_tracker_layer2_hit->SetBranchAddress("momentum_y",&tracker_layer2_momentum_y);
  tree_tracker_layer2_hit->SetBranchAddress("momentum_z",&tracker_layer2_momentum_z);
  tree_tracker_layer2_hit->SetBranchAddress("is_asymmetric_scattering",&tracker_layer2_is_asymmetric_scattering);
}

void Analysis(TFile* outfile){
  Long64_t total_entries = tree_event->GetEntries();
  Int_t i_entry_primary = 0;
  Int_t i_entry_trajectory = 0;
  Int_t i_entry_cdc = 0;
  Int_t i_entry_tracker_layer1 = 0;
  Int_t i_entry_tracker_layer2 = 0;
  Int_t total_entries_cdc = tree_cdc_hit->GetEntries();
  Int_t total_entries_tracker_layer1 = tree_tracker_layer1_hit->GetEntries();
  Int_t total_entries_tracker_layer2 = tree_tracker_layer2_hit->GetEntries();

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

    // tracker_layer1
    std::vector<Int_t> tracker_layer1_parent_ids;
    std::vector<Int_t> tracker_layer1_particle_ids;
    std::vector<TVector3> vec_tracker_layer1_positions;
    std::vector<TVector3> vec_tracker_layer1_momenta;
    std::vector<Int_t> tracker_layer1_are_asymmetric_scattering;
    for(int i_tracker_layer1=0; i_tracker_layer1<number_of_hits_in_tracker_layer1; ++i_tracker_layer1){
      tree_tracker_layer1_hit->GetEntry(i_entry_tracker_layer1);
      ++i_entry_tracker_layer1;
      if(event_id != tracker_layer1_event_id){
        std::cout << "event mismatch!!!" << std::endl;
        return;
      }
      tracker_layer1_parent_ids.push_back(tracker_layer1_parent_id);
      tracker_layer1_particle_ids.push_back(tracker_layer1_particle_id);
      TVector3 vec_tracker_layer1_position((Double_t)tracker_layer1_hit_position_x,(Double_t)tracker_layer1_hit_position_y,(Double_t)tracker_layer1_hit_position_z);
      TVector3 vec_tracker_layer1_momentum((Double_t)tracker_layer1_momentum_x,(Double_t)tracker_layer1_momentum_y,(Double_t)tracker_layer1_momentum_z);
      vec_tracker_layer1_positions.push_back(vec_tracker_layer1_position);
      vec_tracker_layer1_momenta.push_back(vec_tracker_layer1_momentum);
      tracker_layer1_are_asymmetric_scattering.push_back(tracker_layer1_is_asymmetric_scattering);
    }

    // tracker_layer2
    std::vector<Int_t> tracker_layer2_parent_ids;
    std::vector<Int_t> tracker_layer2_particle_ids;
    std::vector<TVector3> vec_tracker_layer2_positions;
    std::vector<TVector3> vec_tracker_layer2_momenta;
    std::vector<Int_t> tracker_layer2_are_asymmetric_scattering;
    for(int i_tracker_layer2=0; i_tracker_layer2<number_of_hits_in_tracker_layer2; ++i_tracker_layer2){
      tree_tracker_layer2_hit->GetEntry(i_entry_tracker_layer2);
      ++i_entry_tracker_layer2;
      if(event_id != tracker_layer2_event_id){
        std::cout << "event mismatch!!!" << std::endl;
        return;
      }
      tracker_layer2_parent_ids.push_back(tracker_layer2_parent_id);
      tracker_layer2_particle_ids.push_back(tracker_layer2_particle_id);
      TVector3 vec_tracker_layer2_position((Double_t)tracker_layer2_hit_position_x,(Double_t)tracker_layer2_hit_position_y,(Double_t)tracker_layer2_hit_position_z);
      TVector3 vec_tracker_layer2_momentum((Double_t)tracker_layer2_momentum_x,(Double_t)tracker_layer2_momentum_y,(Double_t)tracker_layer2_momentum_z);
      vec_tracker_layer2_positions.push_back(vec_tracker_layer2_position);
      vec_tracker_layer2_momenta.push_back(vec_tracker_layer2_momentum);
      tracker_layer2_are_asymmetric_scattering.push_back(tracker_layer2_is_asymmetric_scattering);
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

    TLorentzVector lvec_proton_from_lambda_mc(0.,0.,0.,0.);
    TLorentzVector lvec_pim_from_lambda_mc(0.,0.,0.,0.);
    bool is_lambda_charged_decay = false;
    for(int i_trajectory=0; i_trajectory<number_of_trajectories; ++i_trajectory){
      if(trajectory_parent_ids[i_trajectory]==lambda_id&&trajectory_particle_ids[i_trajectory]==2212){
        lvec_proton_from_lambda_mc.SetVectM(vec_trajectory_initial_momenta[i_trajectory], mass::proton);
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
    for(int i_cdc=0; i_cdc<number_of_hits_in_cdc; ++i_cdc){
      if(cdc_parent_ids[i_cdc]==0&&cdc_particle_ids[i_cdc]==2212){
        is_proton_detected_by_cdc = true;
      }
      if(cdc_parent_ids[i_cdc]==lambda_id&&cdc_particle_ids[i_cdc]==2212){
        is_proton_from_lambda_detected_by_cdc = true;
      }
      if(cdc_parent_ids[i_cdc]==lambda_id&&cdc_particle_ids[i_cdc]==-211){
        is_pim_from_lambda_detected_by_cdc = true;
      }
    }

    // checking tracker_layer1 & layer2 (proton to be detected)
    bool is_proton_detected_by_tracker_layer1 = false;
    bool is_proton_detected_by_tracker_layer2 = false;
    bool is_asymmetric_scattering_layer1 = false;
    bool is_asymmetric_scattering_layer2 = false;
    for(int i_tracker_layer1=0; i_tracker_layer1<number_of_hits_in_tracker_layer1; ++i_tracker_layer1){
      if(tracker_layer1_parent_ids[i_tracker_layer1]==0 || tracker_layer1_parent_ids[i_tracker_layer1]==proton_id){
        if(tracker_layer1_particle_ids[i_tracker_layer1]==2212){
          is_proton_detected_by_tracker_layer1 = true;
          if(tracker_layer1_are_asymmetric_scattering[i_tracker_layer1]){
            is_asymmetric_scattering_layer1 = true;
          }
        }
      }
    }
    for(int i_tracker_layer2=0; i_tracker_layer2<number_of_hits_in_tracker_layer2; ++i_tracker_layer2){
      if(tracker_layer2_parent_ids[i_tracker_layer2]==0 || tracker_layer2_parent_ids[i_tracker_layer2]==proton_id){
        if(tracker_layer2_particle_ids[i_tracker_layer2]==2212){
          is_proton_detected_by_tracker_layer2 = true;
        }
        if(tracker_layer2_are_asymmetric_scattering[i_tracker_layer2]){
          is_asymmetric_scattering_layer2 = true;
        }
      }
    }

    // checking cdc, tracker1, and tracker2 (to measure proton scattering angle)
    TVector3 vec_proton_position_at_cdc(0.,0.,0.);
    TVector3 vec_proton_momentum_at_cdc(0.,0.,0.);
    TVector3 vec_proton_position_at_tracker_layer1(0.,0.,0.);
    TVector3 vec_proton_momentum_at_tracker_layer1(0.,0.,0.);
    TVector3 vec_proton_position_at_tracker_layer2(0.,0.,0.);
    TVector3 vec_proton_momentum_at_tracker_layer2(0.,0.,0.);

    // position & momentum at cdc
    for(int i_cdc=0; i_cdc<number_of_hits_in_cdc; ++i_cdc){
      if(cdc_parent_ids[i_cdc]==0&&cdc_particle_ids[i_cdc]==2212){
        vec_proton_momentum_at_cdc = vec_cdc_momenta[i_cdc];
        vec_proton_position_at_cdc = vec_cdc_positions[i_cdc];
      }
    }

    // position & momentum at tracker layer1
    if(is_proton_detected_by_tracker_layer1){
      for(int i_tracker_layer1=0; i_tracker_layer1<number_of_hits_in_tracker_layer1; ++i_tracker_layer1){
        if(tracker_layer1_parent_ids[i_tracker_layer1]==0&&tracker_layer1_particle_ids[i_tracker_layer1]==2212){
          vec_proton_momentum_at_tracker_layer1 = vec_tracker_layer1_momenta[i_tracker_layer1];
          vec_proton_position_at_tracker_layer1 = vec_tracker_layer1_positions[i_tracker_layer1];
        }
      }
      if(vec_proton_position_at_tracker_layer1.Mag()==0.){
        for(int i_tracker_layer1=0; i_tracker_layer1<number_of_hits_in_tracker_layer1; ++i_tracker_layer1){
          if(tracker_layer1_parent_ids[i_tracker_layer1]==proton_id&&tracker_layer1_particle_ids[i_tracker_layer1]==2212){
            vec_proton_momentum_at_tracker_layer1 = vec_tracker_layer1_momenta[i_tracker_layer1];
            vec_proton_position_at_tracker_layer1 = vec_tracker_layer1_positions[i_tracker_layer1];
          }
        }
      }
    }

    // position & momentum at tracker layer2
    if(is_proton_detected_by_tracker_layer2){
      for(int i_tracker_layer2=0; i_tracker_layer2<number_of_hits_in_tracker_layer2; ++i_tracker_layer2){
        if(tracker_layer2_parent_ids[i_tracker_layer2]==0&&tracker_layer2_particle_ids[i_tracker_layer2]==2212){
          vec_proton_momentum_at_tracker_layer2 = vec_tracker_layer2_momenta[i_tracker_layer2];
          vec_proton_position_at_tracker_layer2 = vec_tracker_layer2_positions[i_tracker_layer2];
        }
      }
      if(vec_proton_position_at_tracker_layer2.Mag()==0.){
        for(int i_tracker_layer2=0; i_tracker_layer2<number_of_hits_in_tracker_layer2; ++i_tracker_layer2){
          if(tracker_layer2_parent_ids[i_tracker_layer2]==proton_id&&tracker_layer2_particle_ids[i_tracker_layer2]==2212){
            vec_proton_momentum_at_tracker_layer2 = vec_tracker_layer2_momenta[i_tracker_layer2];
            vec_proton_position_at_tracker_layer2 = vec_tracker_layer2_positions[i_tracker_layer2];
          }
        }
      }
    }

    Double_t ncbarrel_radius = 535. + 3./2. + 50./2.; /* mm */
    TVector3 vec_proton_position_at_ncbarrel = ncbarrel_radius * vec_proton_position_at_cdc.Unit();

    TVector3 vec_proton_direction_at_cdc = vec_proton_momentum_at_cdc.Unit();
    TVector3 vec_proton_direction_at_tracker_layer1 = vec_proton_momentum_at_tracker_layer1.Unit();
    TVector3 vec_proton_direction_at_tracker_layer2 = vec_proton_momentum_at_tracker_layer2.Unit();

    TVector3 vec_proton_direction_after_scattering = (vec_proton_position_at_tracker_layer2 - vec_proton_position_at_tracker_layer1).Unit();
    TVector3 vec_proton_direction_after_scattering_rough = (vec_proton_position_at_tracker_layer1 - vec_proton_position_at_ncbarrel).Unit();

    // expected spin direction of lambda
    TVector3 vec_lambda_expected_spin_direction = clvec_proton_from_lambda_lambda_rest_mc.Vect().Unit();
    TVector3 vec_lambda_expectec_spin_direction_perp = (vec_lambda_expected_spin_direction - (vec_lambda_expected_spin_direction.Dot(vec_proton_direction_at_cdc))*vec_proton_direction_at_cdc).Unit();

    // scattering angle mc value
    Double_t scattering_angle_theta_mc = acos(vec_proton_direction_at_cdc.Dot(vec_proton_direction_at_tracker_layer1)) /TMath::Pi()*180.;

    // momentum difference between cdc and tracker layer
    Double_t proton_momentum_difference_between_cdc_and_tracker_layer1 = fabs(vec_proton_momentum_at_tracker_layer1.Mag() - vec_proton_momentum_at_cdc.Mag());

    // scattering angle with two trackers
    Double_t scattering_angle_theta = acos(vec_proton_direction_at_cdc.Dot(vec_proton_direction_after_scattering)) /TMath::Pi()*180.;

    // expected spin direction of proton with two trackers
    TVector3 vec_proton_expected_spin_direction = vec_proton_direction_at_cdc.Cross(vec_proton_direction_after_scattering).Unit();

    // phi angle between lambda and proton spins with two trackers
    Double_t phi_of_spins = acos(vec_lambda_expectec_spin_direction_perp.Dot(vec_proton_expected_spin_direction));
    Double_t phi_sign = vec_lambda_expected_spin_direction.Dot(vec_proton_direction_at_cdc.Cross(vec_proton_expected_spin_direction));
    phi_of_spins *= phi_sign/fabs(phi_sign);

    // scattering angle with one tracker
    Double_t scattering_angle_theta_rough = acos(vec_proton_direction_at_cdc.Dot(vec_proton_direction_after_scattering_rough)) /TMath::Pi()*180.;

    // expected spin direction of proton with one tracker
    TVector3 vec_proton_expected_spin_direction_rough = vec_proton_direction_at_cdc.Cross(vec_proton_direction_after_scattering_rough).Unit();

    // phi angle between lambda and proton spins with one tracker
    Double_t phi_of_spins_rough = acos(vec_lambda_expectec_spin_direction_perp.Dot(vec_proton_expected_spin_direction_rough));
    Double_t phi_sign_rough = vec_lambda_expected_spin_direction.Dot(vec_proton_direction_at_cdc.Cross(vec_proton_expected_spin_direction_rough));
    phi_of_spins_rough *= phi_sign_rough/fabs(phi_sign_rough);


    // ====================================================================================================
    // event selections
    // ====================================================================================================
    bool are_all_detected_by_cdc = false;
    if(is_proton_detected_by_cdc&&is_proton_from_lambda_detected_by_cdc&&is_pim_from_lambda_detected_by_cdc){
      are_all_detected_by_cdc = true;
    }
    bool are_momenta_over_threshold = false;
    Double_t proton_momentum_ll = 100.;
    Double_t pim_momentum_ll = 100.;
    if(lvec_proton_mc.P()>proton_momentum_ll&&lvec_proton_from_lambda_mc.P()>proton_momentum_ll&&lvec_pim_from_lambda_mc.P()>pim_momentum_ll){
      are_momenta_over_threshold = true;
    }

    bool is_tracker_layer1_fired = false;
    if(is_proton_detected_by_tracker_layer1){
      is_tracker_layer1_fired = true;
    }

    bool are_both_trackers_fired = false;
    if(is_proton_detected_by_tracker_layer1&&is_proton_detected_by_tracker_layer2){
      are_both_trackers_fired = true;
    }

    bool is_asymmetric_scattering = false;
    if(is_asymmetric_scattering_layer1){
      is_asymmetric_scattering = true;
    }

    bool is_scattering_angle_in_region = false;
    Double_t theta_sel_ll = 6.;
    Double_t theta_sel_ul = 14.;
    if(theta_sel_ll<scattering_angle_theta&&scattering_angle_theta<theta_sel_ul){
      is_scattering_angle_in_region = true;
    }

    // ====================================================================================================


    TString sel_name[10] = {
      "",
      "_charged",
      "_charged_mom",
      "_charged_mom_cdc",
      "_charged_mom_cdc_both",
      "_charged_mom_cdc_both_asym",
      "_charged_mom_cdc_one",
      "_charged_mom_cdc_one_asym",
      "_charged_mom_cdc_lay1",
      "_charged_mom_cdc_lay1_asym"
    };
    bool sel_flag[10];
    sel_flag[0] = true;
    sel_flag[1] = is_lambda_charged_decay;
    sel_flag[2] = is_lambda_charged_decay & are_momenta_over_threshold;
    sel_flag[3] = is_lambda_charged_decay & are_momenta_over_threshold & are_all_detected_by_cdc;
    sel_flag[4] = is_lambda_charged_decay & are_momenta_over_threshold & are_all_detected_by_cdc  & are_both_trackers_fired;
    sel_flag[5] = is_lambda_charged_decay & are_momenta_over_threshold & are_all_detected_by_cdc  & are_both_trackers_fired & is_asymmetric_scattering;
    sel_flag[6] = is_lambda_charged_decay & are_momenta_over_threshold & are_all_detected_by_cdc  & is_tracker_layer1_fired & !are_both_trackers_fired;
    sel_flag[7] = is_lambda_charged_decay & are_momenta_over_threshold & are_all_detected_by_cdc  & is_tracker_layer1_fired & !are_both_trackers_fired & is_asymmetric_scattering;
    sel_flag[8] = is_lambda_charged_decay & are_momenta_over_threshold & are_all_detected_by_cdc  & is_proton_detected_by_tracker_layer1;
    sel_flag[9] = is_lambda_charged_decay & are_momenta_over_threshold & are_all_detected_by_cdc  & is_proton_detected_by_tracker_layer1 & is_asymmetric_scattering;

    // Fill histograms
    for(int i_sel = 0; i_sel<10; ++i_sel){
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
        ((TH1F*)outfile->Get(Form("proton_momentum_at_tracker_layer1%s",sel_name[i_sel].Data())))->Fill(vec_proton_momentum_at_tracker_layer1.Mag()/1000.);
        ((TH1F*)outfile->Get(Form("proton_momentum_at_tracker_layer2%s",sel_name[i_sel].Data())))->Fill(vec_proton_momentum_at_tracker_layer2.Mag()/1000.);

        ((TH1F*)outfile->Get(Form("proton_momentum_difference_between_cdc_and_tracker_layer1%s",sel_name[i_sel].Data())))->Fill(proton_momentum_difference_between_cdc_and_tracker_layer1/1000.);

        ((TH1F*)outfile->Get(Form("scattering_angle_theta_mc%s",sel_name[i_sel].Data())))->Fill(scattering_angle_theta_mc);

        ((TH1F*)outfile->Get(Form("scattering_angle_theta%s",sel_name[i_sel].Data())))->Fill(scattering_angle_theta);
        ((TH1F*)outfile->Get(Form("scattering_angle_theta_rough%s",sel_name[i_sel].Data())))->Fill(scattering_angle_theta_rough);
        ((TH1F*)outfile->Get(Form("phi_of_spins%s",sel_name[i_sel].Data())))->Fill(phi_of_spins);
        ((TH1F*)outfile->Get(Form("phi_of_spins_rough%s",sel_name[i_sel].Data())))->Fill(phi_of_spins_rough);
      }
    }


    // end of event loop
  }

}
