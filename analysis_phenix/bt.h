//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct  4 12:24:02 2017 by ROOT version 5.34/30
// from TTree bt/bt
// found on file: ../output/BtoAll_0_1.root
// last edited by Ziyue Zhang
//////////////////////////////////////////////////////////

#ifndef bt_h
#define bt_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TF2.h>
#include <TH1.h>
#include <TH2.h>
#include <TMatrixD.h>
#include <TString.h>
#include <TRandom.h>
#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
//#include "src/RooUnfoldResponse.h"
//#include "src/RooUnfoldBayes.h"
//#include "src/RooUnfoldSvd.h"
//#include "src/RooUnfoldTUnfold.h"
#endif

using namespace std;

TH1D* pthist_Btoall;

TH1D* Bscale_e;
TH1D* Bscale_D0;
TH1D* Bscale_Jpsi;
//Btoepecific daughter hist, used for scale

TH1D* hist_efromB_selected_rebin;
TH1D* hist_D0fromB;
TH1D* hist_JpsifromB;

TH2D* Mresponse_B_e_rebin;
TH2D* Mresponse_e_rebin_B;
TH2D* Mresponse_D0_B;
TH2D* Mresponse_Jpsi_B;


TH2D* decay_Btoe_selected_rebin;
TH2D* decay_BtoD0;
TH2D* decay_BtoJpsi;

// Header file for the classes stored in the TTree if any.
// Fixed size dimensions of array or collections stored in the TTree if any.

class bt {
   public :
      TTree          *fChain;   //!pointer to the analyzed TTree or TChain
      Int_t           fCurrent; //!current Tree number in a TChain

      // Declaration of leaf types
      Int_t           np;
      Float_t         px[2000];   //[np]
      Float_t         py[2000];   //[np]
      Float_t         pz[2000];   //[np]
      Float_t         pt[2000];   //[np]
      Float_t         eta[2000];   //[np]
      Float_t         phi[2000];   //[np]
      Float_t         m[2000];   //[np]
      Float_t         y[2000];   //[np]
      Int_t           status[2000];   //[np]
      Int_t           pdgid[2000];   //[np]
      Short_t         motherId1[2000];   //[np]
      Short_t         motherId2[2000];   //[np]
      Short_t         dauId1[2000];   //[np]
      Short_t         dauId2[2000];   //[np]
      Short_t         bMotherId[2000];   //[np]

      // List of branches
      TBranch        *b_np;   //!
      TBranch        *b_px;   //!
      TBranch        *b_py;   //!
      TBranch        *b_pz;   //!
      TBranch        *b_pt;   //!
      TBranch        *b_eta;   //!
      TBranch        *b_phi;   //!
      TBranch        *b_m;   //!
      TBranch        *b_y;   //!
      TBranch        *b_status;   //!
      TBranch        *b_pdgid;   //!
      TBranch        *b_motherId1;   //!
      TBranch        *b_motherId2;   //!
      TBranch        *b_dauId1;   //!
      TBranch        *b_dauId2;   //!
      TBranch        *b_bMotherId;   //!

      bt(TString input, TTree *tree=0);
      virtual ~bt();
      virtual Int_t    Cut(Long64_t entry);
      virtual Int_t    GetEntry(Long64_t entry);
      virtual Long64_t LoadTree(Long64_t entry);
      virtual void     Init(TTree *tree);
      virtual void     Loop();
      virtual Bool_t   Notify();
      virtual void     Show(Long64_t entry = -1);

		TString         fOutputFileName;
      void            CreateHist();
      void            SaveHist();
};

#endif

bt::bt(TString input, TTree *tree) : fChain(0) 
{
   // if parameter tree is not specified (or zero), connect the file
   // used to generate this class and read the Tree.
   fOutputFileName = "BtoAll.hist.root";
   //fOutputFileName.ReplaceAll(".root",".hist.root");
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(input);
      if (!f || !f->IsOpen()) {
         f = new TFile(input);
      }
      f->GetObject("bt",tree);
   }
   Init(tree);
}

bt::~bt()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}


Int_t bt::GetEntry(Long64_t entry)
{
   // Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t bt::LoadTree(Long64_t entry)
{
   // Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void bt::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("np", &np, &b_np);
   fChain->SetBranchAddress("px", px, &b_px);
   fChain->SetBranchAddress("py", py, &b_py);
   fChain->SetBranchAddress("pz", pz, &b_pz);
   fChain->SetBranchAddress("pt", pt, &b_pt);
   fChain->SetBranchAddress("eta", eta, &b_eta);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   fChain->SetBranchAddress("m", m, &b_m);
   fChain->SetBranchAddress("y", y, &b_y);
   fChain->SetBranchAddress("status", status, &b_status);
   fChain->SetBranchAddress("pdgid", pdgid, &b_pdgid);
   fChain->SetBranchAddress("motherId1", motherId1, &b_motherId1);
   fChain->SetBranchAddress("motherId2", motherId2, &b_motherId2);
   fChain->SetBranchAddress("dauId1", dauId1, &b_dauId1);
   fChain->SetBranchAddress("dauId2", dauId2, &b_dauId2);
   fChain->SetBranchAddress("bMotherId", bMotherId, &b_bMotherId);
   Notify();
}

Bool_t bt::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void bt::Show(Long64_t entry)
{
   // Print contents of entry.
   // If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t bt::Cut(Long64_t entry)
{
   // This function may be called from Loop.
   // returns  1 if entry is accepted.
   // returns -1 otherwise.
   int Binfo[3][2000]={0};
   for(int i=0;i<np;++i)
   {
      if(bMotherId[i]<0) continue;
      double ptTmp=sqrt(px[i]*px[i]+py[i]*py[i]);
      int pId = bMotherId[i];
      if(fabs(pdgid[i])==11&&status[i]>0&&fabs(y[i])<0.5&&ptTmp>1.&&ptTmp<10.)
         Binfo[0][pId]=Binfo[0][pId]+1;
      if(fabs(pdgid[i])==421&&fabs(y[i])<0.5&&ptTmp<8.)
         Binfo[1][pId]=Binfo[1][pId]+1;
      if(fabs(pdgid[i])==443&&fabs(y[i])<0.5&&ptTmp<8.)
         Binfo[2][pId]=Binfo[2][pId]+1;
   }//Get B counting multiplicity
   for(int i=0;i<np;++i)
   {
      double ptTmp=sqrt(px[i]*px[i]+py[i]*py[i]);
      if(fabs(pdgid[i])>500&&fabs(pdgid[i])<600)
      {
         if(ptTmp<29.5)
         {
            pthist_Btoall->Fill(ptTmp);
         }else{
            pthist_Btoall->Fill(29.75);
         }
      }//Get B real histo
      if(bMotherId[i]<0) continue;
      int pId = bMotherId[i];   
      double parentpt = sqrt(px[pId]*px[pId]+py[pId]*py[pId]);
      if(fabs(pdgid[i])==11&&status[i]>0&&fabs(y[i])<0.5&&ptTmp>1.&&ptTmp<10.)//e
      {
         hist_efromB_selected_rebin->Fill(ptTmp);//Measured
         if(parentpt<29.5){
            Mresponse_e_rebin_B->Fill(parentpt,ptTmp);
            Mresponse_B_e_rebin->Fill(ptTmp,parentpt);
            decay_Btoe_selected_rebin->Fill(parentpt,ptTmp);
            Bscale_e->Fill(parentpt,1./Binfo[0][pId]);//Get B to e real
         }else{
            Mresponse_e_rebin_B->Fill(29.75,ptTmp);
            Mresponse_B_e_rebin->Fill(ptTmp,29.75);
            decay_Btoe_selected_rebin->Fill(29.75,ptTmp);
            Bscale_e->Fill(29.75,1./Binfo[0][pId]);
         }
      }else{
         if(parentpt<29.5){
            decay_Btoe_selected_rebin->Fill(parentpt,ptTmp+10.);
         }else{
            decay_Btoe_selected_rebin->Fill(29.75,ptTmp+10.);
         }
      }



      if(fabs(pdgid[i])==421&&fabs(y[i])<0.5&&ptTmp<8.)//D0 histo
      {
         hist_D0fromB->Fill(ptTmp);
         if(parentpt<29.5){
            Mresponse_D0_B->Fill(parentpt,ptTmp);
            decay_BtoD0->Fill(parentpt,ptTmp);
            Bscale_D0->Fill(parentpt,1./Binfo[1][pId]);
         }else{
            Mresponse_D0_B->Fill(29.75,ptTmp);
            decay_BtoD0->Fill(29.75,ptTmp);
            Bscale_D0->Fill(29.75,1./Binfo[1][pId]);
         }
      }else{
         if(parentpt<29.5){
            decay_BtoD0->Fill(parentpt,ptTmp+8.);
         }else{
            decay_BtoD0->Fill(29.75,ptTmp+8.);
         }
      }
      

      
      if(fabs(pdgid[i])==443&&fabs(y[i])<0.5&&ptTmp<8.)//Jpsi histo
      {
         hist_JpsifromB->Fill(ptTmp);
         if(parentpt<29.5){
            Mresponse_Jpsi_B->Fill(parentpt,ptTmp);
            decay_BtoJpsi->Fill(parentpt,ptTmp);
            Bscale_Jpsi->Fill(parentpt,1./Binfo[2][pId]);
         }else{
            Mresponse_Jpsi_B->Fill(29.75,ptTmp);
            decay_BtoJpsi->Fill(29.75,ptTmp);
            Bscale_Jpsi->Fill(29.75,1./Binfo[2][pId]);
         }
      }else{
         if(parentpt<29.5){
            decay_BtoJpsi->Fill(parentpt,ptTmp+8.);
         }else{
            decay_BtoJpsi->Fill(29.75,ptTmp+8.);
         }
      }
   }   
   return 1;
}

void bt::CreateHist()
{

   //mother binning
   const int nPtBins = 60;
   double mPtMin = 0.;
   double mPtMax = 30.;
   
   //input binning
   // const int nPtBinsInput = 8;
   // const double binEdgesInput[nPtBinsInput+1] = {3,3.5,4,4.5,5,5.5,6.,8.,10.};
   const int nPtBinsInput = 21;//phenix b->e binning
   const double binEdgesInput[nPtBinsInput+1] = {1.0, 1.2, 1.4, 1.6, 1.8, 
                                                 2.0, 2.2, 2.4, 2.6, 2.8,
                                                 3.0, 3.2, 3.4, 3.6, 3.8,
                                                 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0}; 
   
   
   //decay binning
   const int dPtBins = 4;
   const double binEdgesD0[dPtBins+1]={0.,2.,3.,5.,8.};
   const double binEdgesJpsi[dPtBins+1]={0.,2.,3.,5.,8.};
   
   pthist_Btoall = new TH1D("pthist_Btoall","pthist_Btoall",nPtBins,mPtMin,mPtMax);
   Bscale_e = new TH1D("Bscale_e","Bscale_e",nPtBins,mPtMin,mPtMax);
   Bscale_D0 = new TH1D("Bscale_D0","Bscale_D0",nPtBins,mPtMin,mPtMax);
   Bscale_Jpsi = new TH1D("Bscale_Jpsi","Bscale_Jpsi",nPtBins,mPtMin,mPtMax);

   hist_efromB_selected_rebin = new TH1D("hist_efromB_selected_rebin","hist_efromB_selected_rebin",nPtBinsInput,binEdgesInput);
   hist_D0fromB = new TH1D("hist_D0fromB","hist_D0fromB",dPtBins,binEdgesD0);
   hist_JpsifromB = new TH1D("hist_JpsifromB","hist_JpsifromB",dPtBins,binEdgesJpsi);
   
   Mresponse_B_e_rebin = new TH2D("Mresponse_B_e_rebin","Mresponse_B_e_rebin",nPtBinsInput,binEdgesInput,nPtBins,mPtMin,mPtMax);
   Mresponse_e_rebin_B = new TH2D("Mresponse_e_rebin_B","Mresponse_e_rebin_B",nPtBins,mPtMin,mPtMax,nPtBinsInput,binEdgesInput);
   Mresponse_D0_B = new TH2D("Mresponse_D0_B","Mresponse_D0_B",nPtBins,mPtMin,mPtMax,dPtBins,binEdgesD0);
   Mresponse_Jpsi_B = new TH2D("Mresponse_Jpsi_B","Mresponse_Jpsi_B",nPtBins,mPtMin,mPtMax,dPtBins,binEdgesJpsi);
   
   decay_Btoe_selected_rebin = new TH2D("decay_Btoe_selected_rebin","decay_Btoe_selected_rebin",nPtBins,mPtMin,mPtMax,nPtBinsInput,binEdgesInput);
   decay_BtoD0 = new TH2D("decay_BtoD0","decay_BtoD0",nPtBins,mPtMin,mPtMax,dPtBins,binEdgesD0);
   decay_BtoJpsi = new TH2D("decay_BtoJpsi","decay_BtoJpsi",nPtBins,mPtMin,mPtMax,dPtBins,binEdgesJpsi);
}

void bt::SaveHist()
{
   TFile *output=new TFile(fOutputFileName,"recreate");
   output->cd();
   
   pthist_Btoall->Write();
   Bscale_e->Write();
   Bscale_D0->Write();
   Bscale_Jpsi->Write();
   
   hist_efromB_selected_rebin->Write();
   hist_D0fromB->Write();
   hist_JpsifromB->Write();
   
   Mresponse_B_e_rebin->Write();
   Mresponse_e_rebin_B->Write();
   Mresponse_D0_B->Write();
   Mresponse_Jpsi_B->Write();
   
   decay_Btoe_selected_rebin->Write();
   decay_BtoD0->Write();
   decay_BtoJpsi->Write();

   output->Close();
}

void bt::Loop()
{
   //   In a ROOT session, you can do:
   //      Root > .L bt.C
   //      Root > bt t
   //      Root > t.GetEntry(12); // Fill t data members with entry number 12
   //      Root > t.Show();       // Show values of entry 12
   //      Root > t.Show(16);     // Read and show values of entry 16
   //      Root > t.Loop();       // Loop on all entries
   //

   //     This is the loop skeleton where:
   //    jentry is the global entry number in the chain
   //    ientry is the entry number in the current Tree
   //  Note that the argument to GetEntry must be:
   //    jentry for TChain::GetEntry
   //    ientry for TTree::GetEntry and TBranch::GetEntry
   //
   //       To read only selected branches, Insert statements like:
   // METHOD1:
   //    fChain->SetBranchStatus("*",0);  // disable all branches
   //    fChain->SetBranchStatus("branchname",1);  // activate branchname
   // METHOD2: replace line
   //    fChain->GetEntry(jentry);       //read all branches
   //by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   CreateHist();
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (Cut(ientry) < 0) continue;
   }

   SaveHist();
}
