// resolution_functions.hh

#ifndef resolution_functions_hh
#define resolution_functions_hh

#include"TFile.h"
#include"TTree.h"
#include"TVector3.h"

extern void CDCResolution(TVector3& pos, TVector3& mom);
extern void TrackerResolution(TVector3& pos);

#endif
