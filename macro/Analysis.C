// Analysis.C

#include"DrawTools.h"

TFile* input_file = new TFile("root/KNN_simple_polarimeter_scinti_50mm.root");
TFile* output_file = nullptr;

TString pdf_name = "fig/KNN_simple_polarimeter_scinti_50mm_ana.pdf";

void CreateHistograms();
void DrawHistograms();

void Analysis(){
  gROOT->SetBatch(true);


  // initialization ---------------------------------------
  TTree* tree = (TTree*)input_file->Get("EventTree");     
  Int_t dcin_nhit;
  Float_t dcin_position_x, dcin_position_y, dcin_position_z;
  Float_t dcin_momentum_x, dcin_momentum_y, dcin_momentum_z;
  Int_t dcout_nhit;
  Float_t dcout_position_x, dcout_position_y, dcout_position_z;
  Float_t dcout_momentum_x, dcout_momentum_y, dcout_momentum_z;
  Float_t generate_momentum_x, generate_momentum_y, generate_momentum_z;
  Float_t normal_x, normal_y, normal_z;
  Float_t reference_x, reference_y, reference_z;

  tree->SetBranchAddress("dcin_nhit",&dcin_nhit);
  tree->SetBranchAddress("dcin_position_x",&dcin_position_x);
  tree->SetBranchAddress("dcin_position_y",&dcin_position_y);
  tree->SetBranchAddress("dcin_position_z",&dcin_position_z);
  tree->SetBranchAddress("dcin_momentum_x",&dcin_momentum_x);
  tree->SetBranchAddress("dcin_momentum_y",&dcin_momentum_y);
  tree->SetBranchAddress("dcin_momentum_z",&dcin_momentum_z);

  tree->SetBranchAddress("dcout_nhit",&dcout_nhit);
  tree->SetBranchAddress("dcout_position_x",&dcout_position_x);
  tree->SetBranchAddress("dcout_position_y",&dcout_position_y);
  tree->SetBranchAddress("dcout_position_z",&dcout_position_z);
  tree->SetBranchAddress("dcout_momentum_x",&dcout_momentum_x);
  tree->SetBranchAddress("dcout_momentum_y",&dcout_momentum_y);
  tree->SetBranchAddress("dcout_momentum_z",&dcout_momentum_z);

  tree->SetBranchAddress("generate_momentum_x",&generate_momentum_x);
  tree->SetBranchAddress("generate_momentum_y",&generate_momentum_y);
  tree->SetBranchAddress("generate_momentum_z",&generate_momentum_z);
  tree->SetBranchAddress("normal_x",&normal_x);
  tree->SetBranchAddress("normal_y",&normal_y);
  tree->SetBranchAddress("normal_z",&normal_z);
  tree->SetBranchAddress("reference_x",&reference_x);
  tree->SetBranchAddress("reference_y",&reference_y);
  tree->SetBranchAddress("reference_z",&reference_z);

  output_file = new TFile("root/KNN_simple_polarimeter_scinti_50mm_ana.root","RECREATE");
  CreateHistograms();
  // ------------------------------------------------------


  // analysis ---------------------------------------------
  Long64_t total_entries = tree->GetEntries();
  for(Long64_t i_entry=0; i_entry<total_entries; ++i_entry){
    if((i_entry+1)%(Int_t)pow(10, (Int_t)log10(i_entry+1)) ==0 ){
      std::cout << ">> " << i_entry+1 << std::endl;
    }
    tree->GetEntry(i_entry);

    TVector3 vec_dcin_position((Double_t)dcin_position_x,(Double_t)dcin_position_y,(Double_t)dcin_position_z);
    TVector3 vec_dcin_momentum((Double_t)dcin_momentum_x,(Double_t)dcin_momentum_y,(Double_t)dcin_momentum_z);
    TVector3 vec_dcout_position((Double_t)dcout_position_x,(Double_t)dcout_position_y,(Double_t)dcout_position_z);
    TVector3 vec_dcout_momentum((Double_t)dcout_momentum_x,(Double_t)dcout_momentum_y,(Double_t)dcout_momentum_z);
    TVector3 vec_generate_momentum((Double_t)generate_momentum_x,(Double_t)generate_momentum_y,(Double_t)generate_momentum_z);
    TVector3 vec_normal((Double_t)normal_x,(Double_t)normal_y,(Double_t)normal_z);
    TVector3 vec_reference((Double_t)reference_x,(Double_t)reference_y,(Double_t)reference_z);
    const Double_t proton_mass = 0.938272;
    Double_t kinetic_energy = sqrt( pow(proton_mass,2.) + pow(vec_generate_momentum.Mag(),2.) ) - proton_mass;

    bool is_dcin_fired = false;
    if(dcin_nhit){
      is_dcin_fired = true;
    }
    bool is_dcout_fired = false;
    if(dcout_nhit){
      is_dcout_fired = true;
    }

    // angles
    Double_t momentum_cos_theta_to_normal = vec_generate_momentum.Dot(vec_normal) / vec_generate_momentum.Mag() / vec_normal.Mag();
    Double_t scattering_angle = acos(vec_dcin_momentum.Unit().Dot(vec_dcout_momentum.Unit())) * 180./TMath::Pi();

    ((TH1F*)output_file->Get("momentum"))->Fill(vec_generate_momentum.Mag());
    ((TH1F*)output_file->Get("momentum_cos_theta_to_normal"))->Fill(momentum_cos_theta_to_normal);
    ((TH1F*)output_file->Get("scattering_angle"))->Fill(scattering_angle);
    ((TH1F*)output_file->Get("kinetic_energy"))->Fill(kinetic_energy);

    if(is_dcin_fired && is_dcout_fired && 6.<scattering_angle && scattering_angle<30.){
      ((TH1F*)output_file->Get("momentum_sel"))->Fill(vec_generate_momentum.Mag());
      ((TH1F*)output_file->Get("momentum_cos_theta_to_normal_sel"))->Fill(momentum_cos_theta_to_normal);
      ((TH1F*)output_file->Get("scattering_angle_sel"))->Fill(scattering_angle);
      ((TH1F*)output_file->Get("kinetic_energy_sel"))->Fill(kinetic_energy);
    }
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
    new TH1F(Form("momentum%s",cut[icut].Data()),"momentum (generate);#font[12]{p} (MeV/#font[12]{c});counts",200,0.,1.);
    new TH1F(Form("momentum_cos_theta_to_normal%s",cut[icut].Data()),"momentum (generate);#font[12]{p} (MeV/#font[12]{c});counts",200,-1.,1.);
    new TH1F(Form("scattering_angle%s",cut[icut].Data()),"scattering angle;#theta (rad.);counts",400,0.,40.);
    new TH1F(Form("kinetic_energy%s",cut[icut].Data()),"kinetic energy;T (MeV);counts",100,0.,1.);

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
  gPad->SetLogy();
  // --------------

  // momentum
  hist_1d = (TH1F*)output_file->Get("momentum");
  Draw(hist_1d);
  canvas->Print(pdf_name);

  // momentum_cos_theta_to_normal
  hist_1d = (TH1F*)output_file->Get("momentum_cos_theta_to_normal");
  Draw(hist_1d);
  canvas->Print(pdf_name);

  // scattering_angle
  hist_1d = (TH1F*)output_file->Get("scattering_angle");
  Draw(hist_1d);
  canvas->Print(pdf_name);

  // kinetic_energy
  hist_1d = (TH1F*)output_file->Get("kinetic_energy");
  Draw(hist_1d);
  canvas->Print(pdf_name);

  // close pdf file
  canvas->Print(pdf_name + "]");
}
