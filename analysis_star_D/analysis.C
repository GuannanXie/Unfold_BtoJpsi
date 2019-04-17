#include <iostream>
#include <fstream>
#include <cstdlib>
#include "sys/types.h"
#include "dirent.h"

#include "math.h"
#include "string.h"
#include <iomanip>
#include <stdio.h>

#ifndef __CINT__
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TSystem.h"

#include "TH3.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TTree.h"
#include "TString.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TUnixSystem.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#endif

#include "bt.h"


// void analysis(TString input = "WILL BE REPlACED"){
void analysis(TString input = "/star/u/xgn1992/Jpsi/BtoJpsi/sim/jobs_mine/Bdecays/output/output/BtoAll_1.root"){
        std::cout<<"starting analysis..."<<std::endl;
	bt *evt = new bt(input);
	evt->Loop();
}
