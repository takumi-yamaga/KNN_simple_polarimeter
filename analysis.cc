// analysis.cc

#include"analysis_functions.hh"

#include<string>
#include<iostream>
#include<sstream>

#include"TFile.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  std::string in_file_name = "root/mc_out.root";
  std::string out_file_name ="root/ana_out.root";
  std::string out_pdf_name = "fig/ana_out.pdf";

  std::istringstream iss;
  std::cout<<"argc : "<<argc<<std::endl;
  for (int i = 0 ; i < argc ; i++) {
    std::string arg = argv[i];
    std::cout<<"argv["<<i<<"] : "<<argv[i]<<std::endl;
    iss.str("");
    iss.clear();
    if (arg.substr(0, 9) == "--infile=") {
      iss.str(arg.substr(9));
      iss >> in_file_name;
    }
    else if (arg.substr(0, 10) == "--outfile=") {
      iss.str(arg.substr(10));
      iss >> out_file_name;
    }
    else if (arg.substr(0, 9) == "--outpdf=") {
      iss.str(arg.substr(9));
      iss >> out_pdf_name;
    }
  }
  std::cout<<"#############################################################"<<std::endl;
  std::cout<<"--- optional input files [argument] ---"<<std::endl;
  std::cout<<"InFileName        [--infile=]  = "<<in_file_name<<std::endl;
  std::cout<<"OutFileName       [--outfile=] = "<<out_file_name<<std::endl;
  std::cout<<"OutPDFName        [--outpdf=] = "<<out_pdf_name<<std::endl;
  std::cout<<"#############################################################"<<std::endl;


  TFile* infile = new TFile(in_file_name.data(),"read");
  InitializeTrees(infile);

  TFile* outfile = new TFile(out_file_name.data(),"recreate");
  CreateHistograms(outfile);

  Analysis(outfile);
  outfile->Write();

  DrawHistograms(outfile, out_pdf_name);

  outfile->Close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
