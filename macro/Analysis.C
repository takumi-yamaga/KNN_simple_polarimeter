// Analysis.C

#include"DrawTools.h"
#include"MyConstants.h"

TFile* input_file = new TFile("mc_out.root","read");
TFile* output_file = nullptr;

TString pdf_name = "analysis.pdf";

void CreateHistograms();
void DrawHistograms();

void Analysis(){
  gROOT->SetBatch(true);

  // initialization ---------------------------------------
  TTree* tree_event = (TTree*)input_file->Get("tree_event");     
  Int_t event_id;
  Int_t number_of_primaries;
  Int_t number_of_trajectories;
  Int_t number_of_hits_in_cdc;
  Int_t number_of_hits_in_tracker_layer1;
  Int_t number_of_hits_in_tracker_layer2;
  tree_event->SetBranchAddress("event_id",&event_id);
  tree_event->SetBranchAddress("number_of_primaries",&number_of_primaries);
  tree_event->SetBranchAddress("number_of_trajectories",&number_of_trajectories);
  tree_event->SetBranchAddress("number_of_hits_in_cdc",&number_of_hits_in_cdc);
  tree_event->SetBranchAddress("number_of_hits_in_tracker_layer1",&number_of_hits_in_tracker_layer1);
  tree_event->SetBranchAddress("number_of_hits_in_tracker_layer2",&number_of_hits_in_tracker_layer2);
  // -----
  TTree* tree_primary = (TTree*)input_file->Get("tree_primary");     
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
  TTree* tree_trajectory = (TTree*)input_file->Get("tree_trajectory");     
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
  TTree* tree_cdc_hit = (TTree*)input_file->Get("tree_cdc_hit");     
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
  TTree* tree_tracker_layer1_hit = (TTree*)input_file->Get("tree_tracker_layer1_hit");     
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
  TTree* tree_tracker_layer2_hit = (TTree*)input_file->Get("tree_tracker_layer2_hit");     
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
  // ------------------------------------------------------
  output_file = new TFile("ana_out.root","recreate");
  CreateHistograms();
  // ------------------------------------------------------


  // analysis ---------------------------------------------
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
    Double_t lambda_spin_direction_theta = vec_lambda_spin_direction_mc.Theta();
    ((TH1F*)output_file->Get("lambda_spin_direction_theta"))->Fill(lambda_spin_direction_theta);
    Double_t lambda_spin_direction_phi = vec_lambda_spin_direction_mc.Phi();
    ((TH1F*)output_file->Get("lambda_spin_direction_phi"))->Fill(lambda_spin_direction_phi);

    // trajectory
    Int_t lambda_id = -999;
    TLorentzVector lvec_lambda_mc(0.,0.,0.,0.);
    TLorentzVector lvec_proton_mc(0.,0.,0.,0.);
    for(int i_trajectory=0; i_trajectory<number_of_trajectories; ++i_trajectory){
      if(trajectory_parent_ids[i_trajectory]==0&&trajectory_particle_ids[i_trajectory]==3122){
        lambda_id = trajectory_track_ids[i_trajectory];
        lvec_lambda_mc.SetVectM(vec_trajectory_initial_momenta[i_trajectory], mass::lambda);
      }
      if(trajectory_parent_ids[i_trajectory]==0&&trajectory_particle_ids[i_trajectory]==2212){
        lvec_proton_mc.SetVectM(vec_trajectory_initial_momenta[i_trajectory], mass::proton);
      }
    }

    TLorentzVector lvec_proton_from_lambda_mc(0.,0.,0.,0.);
    TLorentzVector lvec_pim_from_lambda_mc(0.,0.,0.,0.);
    for(int i_trajectory=0; i_trajectory<number_of_trajectories; ++i_trajectory){
      if(trajectory_parent_ids[i_trajectory]==lambda_id&&trajectory_particle_ids[i_trajectory]==2212){
        lvec_proton_from_lambda_mc.SetVectM(vec_trajectory_initial_momenta[i_trajectory], mass::proton);
      }
      if(trajectory_parent_ids[i_trajectory]==lambda_id&&trajectory_particle_ids[i_trajectory]==-211){
        lvec_pim_from_lambda_mc.SetVectM(vec_trajectory_initial_momenta[i_trajectory], mass::pi_plus);
      }
    }
    TLorentzVector clvec_proton_from_lambda_lambda_rest_mc = lvec_proton_from_lambda_mc;
    clvec_proton_from_lambda_lambda_rest_mc.Boost(-lvec_lambda_mc.BoostVector());
    //std::cout << clvec_proton_from_lambda_lambda_rest_mc.Px() 
    //  << ", " << clvec_proton_from_lambda_lambda_rest_mc.Py()
    //  << ", " << clvec_proton_from_lambda_lambda_rest_mc.Pz()
    //  << ", " << clvec_proton_from_lambda_lambda_rest_mc.E()
    //  << std::endl;

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
    if(!is_proton_detected_by_cdc||!is_proton_from_lambda_detected_by_cdc||!is_pim_from_lambda_detected_by_cdc){
      continue;
    }

    // checking tracker_layer1 & layer2 (proton to be detected)
    bool is_proton_detected_by_tracker_layer1 = false;
    bool is_proton_detected_by_tracker_layer2 = false;
    bool is_proton_asymmetric_scattering = false;
    for(int i_tracker_layer1=0; i_tracker_layer1<number_of_hits_in_tracker_layer1; ++i_tracker_layer1){
      if(tracker_layer1_parent_ids[i_tracker_layer1]==0&&tracker_layer1_particle_ids[i_tracker_layer1]==2212){
        is_proton_detected_by_tracker_layer1 = true;
        if(tracker_layer1_are_asymmetric_scattering[i_tracker_layer1]){
          is_proton_asymmetric_scattering = true;
        }
      }
    }
    for(int i_tracker_layer2=0; i_tracker_layer2<number_of_hits_in_tracker_layer2; ++i_tracker_layer2){
      if(tracker_layer2_parent_ids[i_tracker_layer2]==0&&tracker_layer2_particle_ids[i_tracker_layer2]==2212){
        is_proton_detected_by_tracker_layer2 = true;
      }
    }
    if(!is_proton_detected_by_tracker_layer1||!is_proton_detected_by_tracker_layer2){
      continue;
    }

    // checking cdc, tracker1, and tracker2 (to measure proton scattering angle)
    TVector3 vec_proton_direction_at_cdc(0.,0.,0.);
    TVector3 vec_proton_position_at_cdc(0.,0.,0.);
    TVector3 vec_proton_position_at_tracker_layer1(0.,0.,0.);
    TVector3 vec_proton_position_at_tracker_layer2(0.,0.,0.);
    for(int i_cdc=0; i_cdc<number_of_hits_in_cdc; ++i_cdc){
      if(cdc_parent_ids[i_cdc]==0&&cdc_particle_ids[i_cdc]==2212){
        vec_proton_direction_at_cdc = vec_cdc_momenta[i_cdc].Unit();
        vec_proton_position_at_cdc = vec_cdc_positions[i_cdc];
      }
    }
    for(int i_tracker_layer1=0; i_tracker_layer1<number_of_hits_in_tracker_layer1; ++i_tracker_layer1){
      if(tracker_layer1_parent_ids[i_tracker_layer1]==0&&tracker_layer1_particle_ids[i_tracker_layer1]==2212){
        vec_proton_position_at_tracker_layer1 = vec_tracker_layer1_positions[i_tracker_layer1];
      }
    }
    for(int i_tracker_layer2=0; i_tracker_layer2<number_of_hits_in_tracker_layer2; ++i_tracker_layer2){
      if(tracker_layer2_parent_ids[i_tracker_layer2]==0&&tracker_layer2_particle_ids[i_tracker_layer2]==2212){
        vec_proton_position_at_tracker_layer2 = vec_tracker_layer2_positions[i_tracker_layer2];
      }
    }
    TVector3 vec_proton_direction_afeter_scattering = (vec_proton_position_at_tracker_layer2 - vec_proton_position_at_tracker_layer1).Unit();

    // expected spin direction of lambda
    TVector3 vec_lambda_expected_spin_direction = clvec_proton_from_lambda_lambda_rest_mc.Vect().Unit();
    TVector3 vec_lambda_expectec_spin_direction_perp = (vec_lambda_expected_spin_direction - (vec_lambda_expected_spin_direction.Dot(vec_proton_direction_at_cdc))*vec_proton_direction_at_cdc).Unit();

    Double_t difference_of_lambda_spin = acos(vec_lambda_spin_direction_mc.Dot(vec_lambda_expected_spin_direction));
    ((TH1F*)output_file->Get("difference_of_lambda_spin"))->Fill(difference_of_lambda_spin);


    // -----
    Double_t scattering_angle_theta = acos(vec_proton_direction_at_cdc.Dot(vec_proton_direction_afeter_scattering)) /TMath::Pi()*180.;
    ((TH1F*)output_file->Get("scattering_angle_theta"))->Fill(scattering_angle_theta);

    // expected spin direction of proton
    TVector3 vec_proton_expected_spin_direction = vec_proton_direction_at_cdc.Cross(vec_proton_direction_afeter_scattering).Unit();

    // phi angle between lambda and proton spins 
    Double_t phi_of_spins = acos(vec_lambda_expectec_spin_direction_perp.Dot(vec_proton_expected_spin_direction));
    Double_t phi_sign = vec_lambda_expected_spin_direction.Dot(vec_proton_direction_at_cdc.Cross(vec_proton_expected_spin_direction));
    phi_of_spins *= phi_sign/fabs(phi_sign);
    ((TH1F*)output_file->Get("phi_of_spins"))->Fill(phi_of_spins);

    // event selection
    Double_t theta_sel_ll = 6.;
    Double_t theta_sel_ul = 14.;
    bool is_theta_selected = false;
    if(theta_sel_ll<scattering_angle_theta&&scattering_angle_theta<theta_sel_ul){
      is_theta_selected = true;
    }

    // after selectino
    //if(!is_theta_selected){
    //  continue;
    //}
    if(!is_proton_asymmetric_scattering){
      continue;
    }
    ((TH1F*)output_file->Get("scattering_angle_theta_sel"))->Fill(scattering_angle_theta);
    ((TH1F*)output_file->Get("phi_of_spins_sel"))->Fill(phi_of_spins);
    ((TH1F*)output_file->Get("difference_of_lambda_spin_sel"))->Fill(difference_of_lambda_spin);

  }
  // ------------------------------------------------------

  DrawHistograms();

  output_file->Write();
  output_file->Close();
}

void CreateHistograms(){
  output_file->cd();

  TString cut[2] = {"", "_sel"};

  for(int icut=0; icut<2; ++icut){

    // 1D histograms
    new TH1F(Form("lambda_spin_direction_theta%s",cut[icut].Data()),"lambda spin direction theta;#theta (rad.);counts",200,0.,TMath::Pi());
    new TH1F(Form("lambda_spin_direction_phi%s",cut[icut].Data()),"lambda spin direction phi;#phi (rad.);counts",200,-TMath::Pi(),TMath::Pi());
    new TH1F(Form("scattering_angle_theta%s",cut[icut].Data()),"scattering angle theta;#theta (deg.);counts",200,0.,50.);
    new TH1F(Form("difference_of_lambda_spin%s",cut[icut].Data()),"difference of measured & generated lambda spins;#theta (rad.);counts",100,0.,TMath::Pi());
    new TH1F(Form("phi_of_spins%s",cut[icut].Data()),"phi between spins;#phi (rad.);counts",10,-TMath::Pi(),TMath::Pi());
  }

}

void DrawHistograms(){
  TCanvas* canvas = new TCanvas("canvas","canvas",500,500);
  canvas->cd();

  // open pdf file
  canvas->Print(pdf_name + "[");

  TH1F* hist_1d = nullptr;
  TH2F* hist_2d = nullptr;

  // --------------
  //gPad->SetLogy();
  // --------------

  // lambda_spin_direction_theta
  hist_1d = (TH1F*)output_file->Get("lambda_spin_direction_theta");
  Draw(hist_1d);
  canvas->Print(pdf_name);

  // lambda_spin_direction_phi
  hist_1d = (TH1F*)output_file->Get("lambda_spin_direction_phi");
  Draw(hist_1d);
  canvas->Print(pdf_name);

  // scattering_angle_theta
  gPad->SetLogy();
  hist_1d = (TH1F*)output_file->Get("scattering_angle_theta");
  Draw(hist_1d);
  canvas->Print(pdf_name);
  gPad->SetLogy(0);

  // difference_of_lambda_spin
  hist_1d = (TH1F*)output_file->Get("difference_of_lambda_spin");
  Draw(hist_1d);
  canvas->Print(pdf_name);

  // phi_of_spins
  hist_1d = (TH1F*)output_file->Get("phi_of_spins");
  Draw(hist_1d);
  hist_1d->SetMinimum(0);
  canvas->Print(pdf_name);

  // -----
  canvas->Clear();
  canvas->Print(pdf_name);
  // -----

  // scattering_angle_theta
  gPad->SetLogy();
  hist_1d = (TH1F*)output_file->Get("scattering_angle_theta_sel");
  Draw(hist_1d);
  canvas->Print(pdf_name);
  gPad->SetLogy(0);

  // difference_of_lambda_spin
  hist_1d = (TH1F*)output_file->Get("difference_of_lambda_spin_sel");
  Draw(hist_1d);
  canvas->Print(pdf_name);

  // phi_of_spins
  hist_1d = (TH1F*)output_file->Get("phi_of_spins_sel");
  Draw(hist_1d);
  hist_1d->SetMinimum(0);
  canvas->Print(pdf_name);

  // phi_of_spins
  hist_1d = (TH1F*)output_file->Get("phi_of_spins_sel");
  Draw(hist_1d);
  hist_1d->Scale(1./(hist_1d->Integral()/hist_1d->GetNbinsX()));
  hist_1d->SetMinimum(0);
  canvas->Print(pdf_name);

  // close pdf file
  canvas->Print(pdf_name + "]");
}
