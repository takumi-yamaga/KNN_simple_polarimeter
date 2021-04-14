// analysis_functions.hh

#ifndef analysis_functions_hh
#define analysis_functions_hh

#include"TFile.h"
#include"TTree.h"

extern void InitializeTrees(TFile*);
extern void CreateHistograms(TFile*);
extern void DrawHistograms(TFile*, std::string);
extern void Analysis(TFile*);

#endif
