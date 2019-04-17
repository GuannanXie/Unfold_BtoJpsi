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
#endif

using namespace std;

TH1D* pthist_Dtoall;
TH1D* pthist_D0toall;

TH1D* Dscale_e;
TH1D* D0scale_e;

TH1D* hist_efromD_selected_rebin;

TH2D* Mresponse_D_e_rebin;
TH2D* Mresponse_D0_e_rebin;
TH2D* Mresponse_e_rebin_D;
TH2D* Mresponse_e_rebin_D0;

TH2D* decay_Dtoe_selected_rebin;
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
   fOutputFileName = "DtoAll.hist.root";
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
   int Dinfo[2][2000]={0};
   for(int i=0;i<np;++i)
   {
      if(bMotherId[i]<0) continue;
      double ptTmp=sqrt(px[i]*px[i]+py[i]*py[i]);
      int pId = bMotherId[i];
      if(fabs(pdgid[pId])==443) continue;
      if(fabs(pdgid[i])==11&&status[i]>0&&fabs(eta[i])<0.35&&ptTmp>1.&&ptTmp<9.)
      {
         Dinfo[0][pId]=Dinfo[0][pId]+1;
         if(fabs(pdgid[pId])==421&&fabs(y[pId])<1.)
            Dinfo[1][pId]=Dinfo[1][pId]+1;
      }
   }//Get D counting multiplicity
   for(int i=0;i<np;++i)
   {
      double ptTmp=sqrt(px[i]*px[i]+py[i]*py[i]);
      if(fabs(pdgid[i])>400&&fabs(pdgid[i])<500&&fabs(pdgid[i])!=443)
      {
         if(ptTmp<14.5)
         {
            pthist_Dtoall->Fill(ptTmp);
         }else{
            pthist_Dtoall->Fill(14.75);
         }
      }//Get D real histo
      if(fabs(pdgid[i])==421&&fabs(y[i])<1.)
      {
         if(ptTmp<14.5)
         {
            pthist_D0toall->Fill(ptTmp);
         }else{
            pthist_D0toall->Fill(14.75);
         }
      }//Get D0 real histo
      
      if(bMotherId[i]<0) continue;
      int pId = bMotherId[i];   
      if(fabs(pdgid[pId])==443) continue;
      double parentpt = sqrt(px[pId]*px[pId]+py[pId]*py[pId]);
      if(fabs(pdgid[i])==11&&status[i]>0&&fabs(eta[i])<0.35&&ptTmp>1.&&ptTmp<9.)//e
      {
         hist_efromD_selected_rebin->Fill(ptTmp);//Measured
         if(parentpt<14.5){
            Mresponse_D_e_rebin->Fill(ptTmp,parentpt);
            Mresponse_e_rebin_D->Fill(parentpt,ptTmp);
            decay_Dtoe_selected_rebin->Fill(parentpt,ptTmp);
            Dscale_e->Fill(parentpt,1./Dinfo[0][pId]);//Get D to e real
            if(fabs(pdgid[pId])==421&&fabs(y[pId])<1.){
               Mresponse_D0_e_rebin->Fill(ptTmp,parentpt);
               Mresponse_e_rebin_D0->Fill(parentpt,ptTmp);
               D0scale_e->Fill(parentpt,1./Dinfo[1][pId]);
            }
         }else{
            Mresponse_D_e_rebin->Fill(ptTmp,14.75);
            Mresponse_e_rebin_D->Fill(14.75,ptTmp);
            decay_Dtoe_selected_rebin->Fill(14.75,ptTmp);
            Dscale_e->Fill(14.75,1./Dinfo[0][pId]);
            if(fabs(pdgid[pId])==421&&fabs(y[pId])<1.){
               Mresponse_D0_e_rebin->Fill(ptTmp,14.75);
               Mresponse_e_rebin_D0->Fill(14.75,ptTmp);
               D0scale_e->Fill(14.75,1./Dinfo[1][pId]);
            }
         }
      }else{
         if(parentpt<14.5){
            decay_Dtoe_selected_rebin->Fill(parentpt,ptTmp+9.);
         }else{
            decay_Dtoe_selected_rebin->Fill(14.75,ptTmp+9.);
         }
      }
   }   
   return 1;
}

void bt::CreateHist()
{

   //mother binning
   const int nPtBins = 30;
   double mPtMin = 0.;
   double mPtMax = 15.;
   
   //input binning
   // const int nPtBinsInput = 8;
   // const double binEdgesInput[nPtBinsInput+1] = {3,3.5,4,4.5,5,5.5,6.,8.,10.};
   const int nPtBinsInput = 21;//phenix b->e binning
   const double binEdgesInput[nPtBinsInput+1] = {1.0, 1.2, 1.4, 1.6, 1.8, 
                                                 2.0, 2.2, 2.4, 2.6, 2.8,
                                                 3.0, 3.2, 3.4, 3.6, 3.8,
                                                 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0}; 
   
   
   pthist_Dtoall = new TH1D("pthist_Dtoall","pthist_Dtoall",nPtBins,mPtMin,mPtMax);
   pthist_D0toall = new TH1D("pthist_D0toall","pthist_D0toall",nPtBins,mPtMin,mPtMax);
   
   Dscale_e = new TH1D("Dscale_e","Dscale_e",nPtBins,mPtMin,mPtMax);
   D0scale_e = new TH1D("D0scale_e","D0scale_e",nPtBins,mPtMin,mPtMax);

   hist_efromD_selected_rebin = new TH1D("hist_efromD_selected_rebin","hist_efromD_selected_rebin",nPtBinsInput,binEdgesInput);
   
   Mresponse_D_e_rebin = new TH2D("Mresponse_D_e_rebin","Mresponse_D_e_rebin",nPtBinsInput,binEdgesInput,nPtBins,mPtMin,mPtMax);
   Mresponse_D0_e_rebin = new TH2D("Mresponse_D0_e_rebin","Mresponse_D0_e_rebin",nPtBinsInput,binEdgesInput,nPtBins,mPtMin,mPtMax);
   
   Mresponse_e_rebin_D = new TH2D("Mresponse_e_rebin_D","Mresponse_e_rebin_D",nPtBins,mPtMin,mPtMax,nPtBinsInput,binEdgesInput);
   Mresponse_e_rebin_D0 = new TH2D("Mresponse_e_rebin_D0","Mresponse_e_rebin_D0",nPtBins,mPtMin,mPtMax,nPtBinsInput,binEdgesInput);
   
   
   decay_Dtoe_selected_rebin = new TH2D("decay_Dtoe_selected_rebin","decay_Dtoe_selected_rebin",nPtBins,mPtMin,mPtMax,nPtBinsInput,binEdgesInput);
}

void bt::SaveHist()
{
   TFile *output=new TFile(fOutputFileName,"recreate");
   output->cd();
   
   pthist_Dtoall->Write();
   pthist_D0toall->Write();

   Dscale_e->Write();
   D0scale_e->Write();

   hist_efromD_selected_rebin->Write();

   Mresponse_D_e_rebin->Write();
   Mresponse_D0_e_rebin->Write();
   
   Mresponse_e_rebin_D->Write();
   Mresponse_e_rebin_D0->Write();

   decay_Dtoe_selected_rebin->Write();
   
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
