#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;
#include "TRandom.h"
#include "TH1D.h"
#include "TFile.h"

#include "TMath.h"
#include "TH2D.h"
#include "TMatrixD.h"

#include "src/RooUnfoldResponse.h"
#include "src/RooUnfoldInvert.h"
#include "src/RooUnfoldBayes.h"
#include "src/RooUnfoldSvd.h"
#include "src/RooUnfoldTUnfold.h"
#endif

void calcInvYield(TH1 *h){
    for(int i=0;i<h->GetNbinsX();++i){
        double x = h->GetBinCenter(i+1);
        double ex = h->GetBinWidth(i+1);
        double y = h->GetBinContent(i+1);
        double ey = h->GetBinError(i+1);
        h->SetBinContent(i+1,y/x/ex);
        h->SetBinError(i+1,ey/x/ex);
        // h->SetBinContent(i+1,y/x/ex/x);
        // h->SetBinError(i+1,ey/x/ex/x);
    }
}

TH1D* calculatescale(TH1D* histreal_mom, TH2D* response_sontomom)
{
    for(int i=0;i<histreal_mom->GetNbinsX();i++)
    {
        double temp= 0.;
        for(int j=0;j<response_sontomom->GetNbinsY();j++)
            temp = temp + response_sontomom->GetBinContent(i+1,j+1);//now temp is the total number including multiplicity;
        if(histreal_mom->GetBinContent(i+1)!=0.)
        {
            double ratio = temp/histreal_mom->GetBinContent(i+1);
            histreal_mom->SetBinError(i+1,TMath::Sqrt((1.+ratio)/histreal_mom->GetBinContent(i+1)));
            histreal_mom->SetBinContent(i+1,ratio);
        }else{
            histreal_mom->SetBinContent(i+1,1.);
            histreal_mom->SetBinError(i+1,0.);
        }
    }
    return(histreal_mom);
}

TH2D* equivalentdecaymatrix(TH2D*decaymatrix, TH1D* scalef)//calculate the decaymatrix*smear matrix
{
    for(int i=0; i< decaymatrix->GetNbinsX();i++)
    {
        for(int j=0;j<decaymatrix->GetNbinsY()+1;j++)
            decaymatrix->SetBinContent(i+1,j+1,decaymatrix->GetBinContent(i+1,j+1)*scalef->GetBinContent(i+1));
    }
    return(decaymatrix);
}

TH2D *DecayMatrixNormalization(TH2D* decaymatrix, TH1D* histmom)
{
    int NoC=decaymatrix->GetNbinsX();
    int NoR=decaymatrix->GetNbinsY();
    for(int i=0;i<NoC;i++)
    {   
        if(histmom->GetBinContent(i+1)!=0.)
        {
            for(int j=0;j<NoR+1;j++)
                decaymatrix->SetBinContent(i+1,j+1,decaymatrix->GetBinContent(i+1,j+1)/histmom->GetBinContent(i+1));
        }
    }
    return(decaymatrix);
}

TH1D *Decay(TH2D* decaymatrix, TH1D* histmom, TMatrixD* histmomcov, TH1D* histson)
{
    double err[1000]={0.};
    TMatrixD* temp= new TMatrixD(decaymatrix->GetNbinsY()+1,decaymatrix->GetNbinsX());
    for(int i=0;i<decaymatrix->GetNbinsY()+1;i++)
    {
        for(int j=0;j<decaymatrix->GetNbinsX();j++)
            temp(i,j)=decaymatrix->GetBinContent(j+1,i+1);
    }
    for(int i=0;i<decaymatrix->GetNbinsY()+1;i++)
    {
        for(int j=0;j<1000;j++)
            err[j]=0.;
        double sum=0.;
        for(int j=0;j<decaymatrix->GetNbinsX();j++)
        {
            //sum=sum+decaymatrix->GetBinContent(j+1,i+1) * histmom->GetBinContent(j+1);
            sum=sum+temp(i,j) * histmom->GetBinContent(j+1);
            for(int k=0;k<decaymatrix->GetNbinsX();k++)
            {
                double ASAT=temp(i,j)*(histmomcov(j,k)+histmomcov(k,j))/2.*temp(i,k);
                if(ASAT!=0.)
                {
                    int idx =(int)TMath::Log10(fabs(ASAT));
                    err[idx+500]=err[idx+500]+ASAT;
                }else{
                    err[500]=err[500]+ASAT;
                }
            }
        }
        double errsquare=0.;
        for(int j=0;j<1000;j++)
            errsquare=errsquare+err[j];
        histson->SetBinContent(i+1,sum);
        histson->SetBinError(i+1,TMath::Sqrt(errsquare));
    }
    return(histson);
}

TH1D*Bscale(TH1D*Bhist, TH1D*scalef)//from real to smeared
{
    for(int i=0;i<Bhist->GetNbinsX();i++)
    {
        Bhist->SetBinContent(i+1,Bhist->GetBinContent(i+1)*scalef->GetBinContent(i+1));
        Bhist->SetBinError(i+1,Bhist->GetBinError(i+1)*scalef->GetBinContent(i+1));
    }
    return(Bhist);
}

void unfold_decay()
{
    // TFile* fin = new TFile("out/DataInput.root");
    // TH1D* dataBtoe=(TH1D*)fin->Get("gCountsBtoe_highpt");
    TFile* fin = new TFile("out/DataInput_phenix.root");//update to phenix input
    TH1D* dataBtoe=(TH1D*)fin->Get("gCountsBtoe_phenix"); //Side product: Btoe pt:
    TH1D* gPhenixB=(TH1D*)fin->Get("gYieldB_phenix"); //Side product: B pt:


    TFile *fsim = new TFile("rootfiles/BtoAll.phenix.hist.root");//new covMatrix root file with Jpsi decay opened
    TH1D* pthist_Btoall = (TH1D*)fsim->Get("pthist_Btoall");
    TH1D* hist_efromB_selected_rebin = (TH1D*)fsim->Get("hist_efromB_selected_rebin");
    TH1D* hist_D0fromB = (TH1D*)fsim->Get("hist_D0fromB");
    TH1D* hist_JpsifromB = (TH1D*)fsim->Get("hist_JpsifromB");
    
    //Now it's just the Btoe real value, we are going to convert it into scale factor
    TH1D* Bscale_e =(TH1D*)fsim->Get("Bscale_e");
    TH1D* Bscale_D0 =(TH1D*)fsim->Get("Bscale_D0");
    TH1D* Bscale_Jpsi =(TH1D*)fsim->Get("Bscale_Jpsi");
    
    TH2D* Mresponse_B_e_rebin = (TH2D*)fsim->Get("Mresponse_B_e_rebin");
    TH2D* Mresponse_e_rebin_B = (TH2D*)fsim->Get("Mresponse_e_rebin_B");
    TH2D* Mresponse_D0_B = (TH2D*)fsim->Get("Mresponse_D0_B");
    TH2D* Mresponse_Jpsi_B = (TH2D*)fsim->Get("Mresponse_Jpsi_B");
    
    //calculate scale factor
    Bscale_e = calculatescale(Bscale_e, Mresponse_e_rebin_B);
    Bscale_D0 = calculatescale(Bscale_D0, Mresponse_D0_B);
    Bscale_Jpsi = calculatescale(Bscale_Jpsi, Mresponse_Jpsi_B);
    
    //calculate Xovere scale factor
    TH1D* Bscale_eovere=(TH1D*)Bscale_e->Clone("Bscale_eovere");
    Bscale_eovere->Reset("ICES");
    for(int i=0;i<Bscale_eovere->GetNbinsX();i++)
        Bscale_eovere->SetBinContent(i+1,1./Bscale_e->GetBinContent(i+1));
    
    TH1D* Bscale_D0overe=(TH1D*)Bscale_D0->Clone("Bscale_D0overe");
    Bscale_D0overe->Reset("ICES");
    for(int i=0;i<Bscale_D0overe->GetNbinsX();i++)
        Bscale_D0overe->SetBinContent(i+1,1./Bscale_e->GetBinContent(i+1));
        
    TH1D* Bscale_Jpsiovere=(TH1D*)Bscale_Jpsi->Clone("Bscale_Jpsiovere");
    Bscale_Jpsiovere->Reset("ICES");
    for(int i=0;i<Bscale_Jpsiovere->GetNbinsX();i++)
        Bscale_Jpsiovere->SetBinContent(i+1,1./Bscale_e->GetBinContent(i+1));
    
    //create container for smeared B
    TH1D* B_smeared_e = (TH1D*)pthist_Btoall->Clone("B_smeared_e");
    TH1D* B_smeared_D0 = (TH1D*)pthist_Btoall->Clone("B_smeared_D0");
    TH1D* B_smeared_Jpsi = (TH1D*)pthist_Btoall->Clone("B_smeared_Jpsi");
    
    //smear the B spectra (this is for RooUnfoldResponse constructing)
    B_smeared_e = Bscale(B_smeared_e,Bscale_e);
    B_smeared_D0 = Bscale(B_smeared_D0,Bscale_D0);
    B_smeared_Jpsi = Bscale(B_smeared_Jpsi,Bscale_Jpsi);
    
    //construct the correct RooUnfoldResponse, but the unfolded B is the smeared B, in the end we have to make the output antiscale; for fold we need to input the smeared/scaled B spectrum
    RooUnfoldResponse* response_B_e_rebin = new RooUnfoldResponse(hist_efromB_selected_rebin,B_smeared_e,Mresponse_B_e_rebin,"response_B_e_rebin","response_B_e_rebin");
    RooUnfoldResponse* response_e_rebin_B = new RooUnfoldResponse(B_smeared_e,hist_efromB_selected_rebin,Mresponse_e_rebin_B,"response_e_rebin_B","response_e_rebin_B");
    RooUnfoldResponse* response_D0_B = new RooUnfoldResponse(B_smeared_D0,hist_D0fromB,Mresponse_D0_B,"response_D0_B","response_D0_B");
    RooUnfoldResponse* response_Jpsi_B = new RooUnfoldResponse(B_smeared_Jpsi,hist_JpsifromB,Mresponse_Jpsi_B,"response_Jpsi_B","response_Jpsi_B");

    TH2D* decay_Btoe_selected_rebin = (TH2D*)fsim->Get("decay_Btoe_selected_rebin");
    TH2D* decay_BtoD0 = (TH2D*)fsim->Get("decay_BtoD0");
    TH2D* decay_BtoJpsi = (TH2D*)fsim->Get("decay_BtoJpsi");
    
    TH2D* EQdecaymatrixe = new TH2D("EQdecaymatrixe","EQdecaymatrixe",pthist_Btoall->GetNbinsX(),pthist_Btoall->GetBinLowEdge(1),pthist_Btoall->GetBinLowEdge(pthist_Btoall->GetNbinsX()+1),pthist_Btoall->GetNbinsX(),pthist_Btoall->GetBinLowEdge(1),pthist_Btoall->GetBinLowEdge(pthist_Btoall->GetNbinsX()+1));
    
    TH2D* EQdecaymatrixD0 = EQdecaymatrixe->Clone("EQdecaymatrixD0");
    TH2D* EQdecaymatrixJpsi = EQdecaymatrixe->Clone("EQdecaymatrixJpsi");
    
    for(int i =0; i<pthist_Btoall->GetNbinsX();i++)
    {
        EQdecaymatrixe->SetBinContent(i+1,i+1,1./Bscale_e->GetBinContent(i+1));
        EQdecaymatrixD0->SetBinContent(i+1,i+1,Bscale_D0->GetBinContent(i+1)/Bscale_e->GetBinContent(i+1));
        EQdecaymatrixJpsi->SetBinContent(i+1,i+1,Bscale_Jpsi->GetBinContent(i+1)/Bscale_e->GetBinContent(i+1));
    }
        
    //decay matrix normalization
    decay_Btoe_selected_rebin=(TH2D*)DecayMatrixNormalization(decay_Btoe_selected_rebin,pthist_Btoall);
    decay_BtoD0=(TH2D*)DecayMatrixNormalization(decay_BtoD0,pthist_Btoall);
    decay_BtoJpsi=(TH2D*)DecayMatrixNormalization(decay_BtoJpsi,pthist_Btoall);
    
    //calculate the equivalent decay matrix for D0 and Jpsi
    decay_Btoe_selected_rebin=equivalentdecaymatrix(decay_Btoe_selected_rebin,Bscale_eovere);
    decay_BtoD0=equivalentdecaymatrix(decay_BtoD0,Bscale_D0overe);
    decay_BtoJpsi=equivalentdecaymatrix(decay_BtoJpsi,Bscale_Jpsiovere);
    
    //container for the reinversed
    TH1D* reinversed_efromB_selected_rebin = (TH1D*)hist_efromB_selected_rebin->Clone("reinversed_efromB_selected_rebin");
    TH1D* reinversed_D0fromB = (TH1D*)hist_D0fromB->Clone("reinversed_D0fromB");
    TH1D* reinversed_JpsifromB = (TH1D*)hist_JpsifromB->Clone("reinversed_JpsifromB");
    
    //container for the decayed
    TH1D* decayed_efromB_selected_rebin = (TH1D*)hist_efromB_selected_rebin->Clone("decayed_efromB_selected_rebin");
    TH1D* decayed_D0fromB = (TH1D*)hist_D0fromB->Clone("decayed_D0fromB");
    TH1D* decayed_JpsifromB = (TH1D*)hist_JpsifromB->Clone("decayed_JpsifromB");
    
    //reset those containers
    reinversed_efromB_selected_rebin->Reset("ICES");
    reinversed_D0fromB->Reset("ICES");
    reinversed_JpsifromB->Reset("ICES");
    
    decayed_efromB_selected_rebin->Reset("ICES");
    decayed_D0fromB->Reset("ICES");
    decayed_JpsifromB->Reset("ICES");
    
    int nIter = 25;//100;//6;
    //unfold==================================================================================================
    //unfold->inverse for consistency check
    RooUnfoldBayes inverse_BtoAll(response_B_e_rebin,hist_efromB_selected_rebin,nIter);
    TH1D* inversedBtoAll_smear_e = (TH1D*)((TH1D*)(inverse_BtoAll.Hreco(2)))->Clone("inversedBtoAll_smear_e");
    TMatrixD*inversedBtoAll_smear_ecov = (TMatrixD*)(inverse_BtoAll.Ereco(2));
    
    //unsmear
    TH1D* inversedBtoAll=(TH1D*)inversedBtoAll_smear_e->Clone("inversedBtoAll");
    inversedBtoAll->Reset("ICES");
    inversedBtoAll=Decay(EQdecaymatrixe,inversedBtoAll_smear_e,inversedBtoAll_smear_ecov,inversedBtoAll);
    
    //smear to the other daughters
    TH1D* inversedBtoAll_smear_D0 = inversedBtoAll_smear_e->Clone("inversedBtoAll_smear_D0");
    inversedBtoAll_smear_D0->Reset("ICES");
    inversedBtoAll_smear_D0 = Decay(EQdecaymatrixD0,inversedBtoAll_smear_e,inversedBtoAll_smear_ecov,inversedBtoAll_smear_D0);
    
    TH1D* inversedBtoAll_smear_Jpsi = inversedBtoAll_smear_e->Clone("inversedBtoAll_smear_Jpsi");
    inversedBtoAll_smear_Jpsi->Reset("ICES");
    inversedBtoAll_smear_Jpsi = Decay(EQdecaymatrixJpsi,inversedBtoAll_smear_e,inversedBtoAll_smear_ecov,inversedBtoAll_smear_Jpsi);
    
    //unfold for data
    RooUnfoldBayes unfold_BtoAll(response_B_e_rebin,dataBtoe,nIter);
    TH1D* unfoldedBtoAll_smear_e = (TH1D*)((TH1D*)unfold_BtoAll.Hreco(2))->Clone("unfoldedBtoAll_smear_e");
    TMatrixD*unfoldedBtoAll_smear_ecov = (TMatrixD*)(unfold_BtoAll.Ereco(2));
    
    //unsmear
    TH1D* unfoldedBtoAll=(TH1D*)unfoldedBtoAll_smear_e->Clone("unfoldedBtoAll");
    unfoldedBtoAll->Reset("ICES");
    unfoldedBtoAll=Decay(EQdecaymatrixe,unfoldedBtoAll_smear_e,unfoldedBtoAll_smear_ecov,unfoldedBtoAll);
    
    //smear to the other daughters
    TH1D* unfoldedBtoAll_smear_D0 = unfoldedBtoAll_smear_e->Clone("unfoldedBtoAll_smear_D0");
    unfoldedBtoAll_smear_D0->Reset("ICES");
    unfoldedBtoAll_smear_D0 = Decay(EQdecaymatrixD0,unfoldedBtoAll_smear_e,unfoldedBtoAll_smear_ecov,unfoldedBtoAll_smear_D0);
    
    TH1D* unfoldedBtoAll_smear_Jpsi = unfoldedBtoAll_smear_e->Clone("unfoldedBtoAll_smear_Jpsi");
    unfoldedBtoAll_smear_Jpsi->Reset("ICES");
    unfoldedBtoAll_smear_Jpsi = Decay(EQdecaymatrixJpsi,unfoldedBtoAll_smear_e,unfoldedBtoAll_smear_ecov,unfoldedBtoAll_smear_Jpsi);
    
    //reinverse for consistency check
    reinversed_efromB_selected_rebin = (TH1D*)Decay(decay_Btoe_selected_rebin,inversedBtoAll_smear_e,inversedBtoAll_smear_ecov,reinversed_efromB_selected_rebin);
    reinversed_D0fromB = (TH1D*)Decay(decay_BtoD0,inversedBtoAll_smear_e,inversedBtoAll_smear_ecov,reinversed_D0fromB);
    reinversed_JpsifromB = (TH1D*)Decay(decay_BtoJpsi,inversedBtoAll_smear_e,inversedBtoAll_smear_ecov,reinversed_JpsifromB);
    
    //decay the data==========================================================================================
    decayed_efromB_selected_rebin = (TH1D*)Decay(decay_Btoe_selected_rebin,unfoldedBtoAll_smear_e,unfoldedBtoAll_smear_ecov,decayed_efromB_selected_rebin);
    decayed_D0fromB = (TH1D*)Decay(decay_BtoD0,unfoldedBtoAll_smear_e,unfoldedBtoAll_smear_ecov,decayed_D0fromB);
    decayed_JpsifromB = (TH1D*)Decay(decay_BtoJpsi,unfoldedBtoAll_smear_e,unfoldedBtoAll_smear_ecov,decayed_JpsifromB);
    
    //fold the SIMU to check==================================================================================
    RooUnfoldBayes reinverse_fold_efromB(response_e_rebin_B,inversedBtoAll_smear_e,1);
    TH1D* reinversed_folded_efromB = (TH1D*)((TH1D*)(reinverse_fold_efromB.Hreco(2)))->Clone("reinversed_folded_efromB");
    
    RooUnfoldBayes reinverse_fold_D0fromB(response_D0_B,inversedBtoAll_smear_D0,1);
    TH1D* reinversed_folded_D0fromB = (TH1D*)((TH1D*)(reinverse_fold_D0fromB.Hreco(2)))->Clone("reinversed_folded_D0fromB");
    
    RooUnfoldBayes reinverse_fold_JpsifromB(response_Jpsi_B,inversedBtoAll_smear_Jpsi,1);
    TH1D* reinversed_folded_JpsifromB = (TH1D*)((TH1D*)(reinverse_fold_JpsifromB.Hreco(2)))->Clone("reinversed_folded_JpsifromB");
    
    //fold the Data to check the decay matrix mult method=====================================================
    RooUnfoldBayes decay_fold_efromB(response_e_rebin_B,unfoldedBtoAll_smear_e,1);
    TH1D* decayed_folded_efromB = (TH1D*)((TH1D*)(decay_fold_efromB.Hreco(2)))->Clone("decayed_folded_efromB");
    
    RooUnfoldBayes decay_fold_D0fromB(response_D0_B,unfoldedBtoAll_smear_D0,1);
    TH1D* decayed_folded_D0fromB = (TH1D*)((TH1D*)(decay_fold_D0fromB.Hreco(2)))->Clone("decayed_folded_D0fromB");
    
    RooUnfoldBayes decay_fold_JpsifromB(response_Jpsi_B,unfoldedBtoAll_smear_Jpsi,1);
    TH1D* decayed_folded_JpsifromB = (TH1D*)((TH1D*)(decay_fold_JpsifromB.Hreco(2)))->Clone("decayed_folded_JpsifromB");

    //===========
    //cv5
    calcInvYield(decayed_D0fromB);
    calcInvYield(decayed_folded_D0fromB);
    //cv5-2
    calcInvYield(decayed_JpsifromB);
    calcInvYield(decayed_folded_JpsifromB);
    calcInvYield(unfoldedBtoAll);
    //cv1
    calcInvYield(pthist_Btoall);
    calcInvYield(inversedBtoAll);
    //cv2
    calcInvYield(hist_efromB_selected_rebin);
    calcInvYield(reinversed_efromB_selected_rebin);
    calcInvYield(reinversed_folded_efromB);
    //cv3
    calcInvYield(hist_D0fromB);
    calcInvYield(reinversed_D0fromB);
    calcInvYield(reinversed_folded_D0fromB);
    calcInvYield(hist_JpsifromB);
    calcInvYield(reinversed_JpsifromB);
    calcInvYield(reinversed_folded_JpsifromB);
    //cv4
    calcInvYield(dataBtoe);
    calcInvYield(decayed_efromB_selected_rebin);
    calcInvYield(decayed_folded_efromB);
    //
    //===========
    //canvas set A: simu-self check
    //c1: Step 1: unfolded B and direct statistics of B
    TCanvas*c1=new TCanvas("c1","c1",1000,1000);
    gPad->SetLogy();

    pthist_Btoall->SetLineStyle(1);
    pthist_Btoall->SetMarkerStyle(1);
    pthist_Btoall->SetLineColor(1);
    pthist_Btoall->SetMarkerColor(1);
    pthist_Btoall->SetLineWidth(2);
    pthist_Btoall->SetMarkerSize(2);
    pthist_Btoall->Draw("E1");
    
    inversedBtoAll->SetLineStyle(2);
    inversedBtoAll->SetMarkerStyle(2);
    inversedBtoAll->SetLineColor(2);
    inversedBtoAll->SetMarkerColor(2);
    inversedBtoAll->SetLineWidth(1);
    inversedBtoAll->SetMarkerSize(1);
    inversedBtoAll->Draw("E1 SAME");

    TLegend* leg1 = new TLegend(0.4,0.7,0.92,0.92);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.035);
    leg1->AddEntry(pthist_Btoall,"SIMU B pT","l");
    leg1->AddEntry(inversedBtoAll,"Unfolded SIMU B pT (from e pT)","l");
    leg1->Draw();
    // c1->SaveAs("inversedBtoAll.eps");
    c1->SaveAs("inversedBtoAll.png");

    //c2: Step 2: folded/decayed e pT by using unfolded SIMU B pT
    TCanvas*c2=new TCanvas("c2","c2",1000,1000);
    gPad->SetLogy();

    hist_efromB_selected_rebin->SetLineStyle(1);
    hist_efromB_selected_rebin->SetMarkerStyle(1);
    hist_efromB_selected_rebin->SetLineColor(1);
    hist_efromB_selected_rebin->SetMarkerColor(1);
    hist_efromB_selected_rebin->SetLineWidth(3);
    hist_efromB_selected_rebin->SetMarkerSize(3);
    hist_efromB_selected_rebin->Draw("E1");
    
    reinversed_efromB_selected_rebin->SetLineStyle(2);
    reinversed_efromB_selected_rebin->SetMarkerStyle(2);
    reinversed_efromB_selected_rebin->SetLineColor(2);
    reinversed_efromB_selected_rebin->SetMarkerColor(2);
    reinversed_efromB_selected_rebin->SetLineWidth(2);
    reinversed_efromB_selected_rebin->SetMarkerSize(2);
    reinversed_efromB_selected_rebin->Draw("E1 SAME");
    
    reinversed_folded_efromB->SetLineStyle(3);
    reinversed_folded_efromB->SetMarkerStyle(3);
    reinversed_folded_efromB->SetLineColor(3);
    reinversed_folded_efromB->SetMarkerColor(3);
    reinversed_folded_efromB->SetLineWidth(1);
    reinversed_folded_efromB->SetMarkerSize(1);
    reinversed_folded_efromB->Draw("E1 SAME");
    
    TLegend* leg2 = new TLegend(0.4,0.7,0.92,0.92);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.035);
    leg2->AddEntry(hist_efromB_selected_rebin,"SIMU e from B pT","l");
    leg2->AddEntry(reinversed_efromB_selected_rebin,"SIMU e: Unfold---SIMU e pT--->Matrix Mult","l");
    leg2->AddEntry(reinversed_folded_efromB,"SIMU e: Unfold---SIMU e pT--->Fold","l");
    leg2->Draw();
    // c2->SaveAs("reinversed_folded_efromB.eps");
    c2->SaveAs("reinversed_folded_efromB.png");
    
    //Step 3: folded/decayed D0/Jpsi pT by using unfolded SIMU B pT
    TCanvas*c3=new TCanvas("c3","c3",1200,600);
    c3->Divide(2,1,0.0001,0.0001);
    //c3_1: D0
    c3->cd(1);
    gPad->SetLogy();
    
    hist_D0fromB->SetLineStyle(1);
    hist_D0fromB->SetMarkerStyle(1);
    hist_D0fromB->SetLineColor(1);
    hist_D0fromB->SetMarkerColor(1);
    hist_D0fromB->SetLineWidth(3);
    hist_D0fromB->SetMarkerSize(3);
    hist_D0fromB->Draw("E1");
    
    reinversed_D0fromB->SetLineStyle(2);
    reinversed_D0fromB->SetMarkerStyle(2);
    reinversed_D0fromB->SetLineColor(2);
    reinversed_D0fromB->SetMarkerColor(2);
    reinversed_D0fromB->SetLineWidth(2);
    reinversed_D0fromB->SetMarkerSize(2);
    reinversed_D0fromB->Draw("E1 SAME");
    
    reinversed_folded_D0fromB->SetLineStyle(3);
    reinversed_folded_D0fromB->SetMarkerStyle(3);
    reinversed_folded_D0fromB->SetLineColor(3);
    reinversed_folded_D0fromB->SetMarkerColor(3);
    reinversed_folded_D0fromB->SetLineWidth(1);
    reinversed_folded_D0fromB->SetMarkerSize(1);
    reinversed_folded_D0fromB->Draw("E1 SAME");
    
    TLegend* leg3_1 = new TLegend(0.4,0.7,0.92,0.92);
    leg3_1->SetBorderSize(0);
    leg3_1->SetTextSize(0.035);
    leg3_1->AddEntry(hist_D0fromB,"SIMU D0 from B pT","l");
    leg3_1->AddEntry(reinversed_D0fromB,"SIMU D0: Unfold---SIMU e pT--->Matrix Mult","l");
    leg3_1->AddEntry(reinversed_folded_D0fromB,"SIMU D0: Unfold---SIMU e PT--->Fold","l");
    leg3_1->Draw();
    
    //c3_2: Jpsi
    c3->cd(2);
    gPad->SetLogy();
    
    hist_JpsifromB->SetLineStyle(1);
    hist_JpsifromB->SetMarkerStyle(1);
    hist_JpsifromB->SetLineColor(1);
    hist_JpsifromB->SetMarkerColor(1);
    hist_JpsifromB->SetLineWidth(3);
    hist_JpsifromB->SetMarkerSize(3);
    hist_JpsifromB->Draw("E1");
    
    reinversed_JpsifromB->SetLineStyle(2);
    reinversed_JpsifromB->SetMarkerStyle(2);
    reinversed_JpsifromB->SetLineColor(2);
    reinversed_JpsifromB->SetMarkerColor(2);
    reinversed_JpsifromB->SetLineWidth(2);
    reinversed_JpsifromB->SetMarkerSize(2);
    reinversed_JpsifromB->Draw("E1 SAME");
    
    reinversed_folded_JpsifromB->SetLineStyle(3);
    reinversed_folded_JpsifromB->SetMarkerStyle(3);
    reinversed_folded_JpsifromB->SetLineColor(3);
    reinversed_folded_JpsifromB->SetMarkerColor(3);
    reinversed_folded_JpsifromB->SetLineWidth(1);
    reinversed_folded_JpsifromB->SetMarkerSize(1);
    reinversed_folded_JpsifromB->Draw("E1 SAME");
    
    TLegend* leg3_2 = new TLegend(0.4,0.7,0.92,0.92);
    leg3_2->SetBorderSize(0);
    leg3_2->SetTextSize(0.035);
    leg3_2->AddEntry(hist_JpsifromB,"SIMU Jpsi from B pT","l");
    leg3_2->AddEntry(reinversed_JpsifromB,"SIMU Jpsi: Unfold---SIMU e pT--->Matrix Mult","l");
    leg3_2->AddEntry(reinversed_folded_JpsifromB,"SIMU Jpsi: Unfold---SIMU e PT--->Fold","l");
    leg3_2->Draw();
    // c3->SaveAs("reinversed_folded_DandJpsifromB.eps");
    c3->SaveAs("reinversed_folded_DandJpsifromB.png");
    
    //Canvas set B: data self check via e pT
    TCanvas*c4=new TCanvas("c4","c4",1000,1000);
    gPad->SetLogy();

    dataBtoe->GetXaxis()->SetTitle("p_{T}(GeV/c)");
    dataBtoe->GetYaxis()->SetTitle("Ed^{3}#sigma/dp^{3}(mbGeV^{-2}c^{3})");
    dataBtoe->GetXaxis()->SetLabelFont(42);
    dataBtoe->GetYaxis()->SetLabelFont(42);
    dataBtoe->GetXaxis()->SetTitleFont(42);
    dataBtoe->GetYaxis()->SetTitleFont(42);
    dataBtoe->SetLineStyle(1);
    dataBtoe->SetMarkerStyle(1);
    dataBtoe->SetLineColor(1);
    dataBtoe->SetMarkerColor(1);
    dataBtoe->SetLineWidth(3);
    dataBtoe->SetMarkerSize(3);
    dataBtoe->Draw("E1");
    
    decayed_efromB_selected_rebin->SetLineStyle(2);
    decayed_efromB_selected_rebin->SetMarkerStyle(2);
    decayed_efromB_selected_rebin->SetLineColor(2);
    decayed_efromB_selected_rebin->SetMarkerColor(2);
    decayed_efromB_selected_rebin->SetLineWidth(2);
    decayed_efromB_selected_rebin->SetMarkerSize(2);
    decayed_efromB_selected_rebin->Draw("E1 SAME");
    
    decayed_folded_efromB->SetLineStyle(3);
    decayed_folded_efromB->SetMarkerStyle(3);
    decayed_folded_efromB->SetLineColor(3);
    decayed_folded_efromB->SetMarkerColor(3);
    decayed_folded_efromB->SetLineWidth(1);
    decayed_folded_efromB->SetMarkerSize(1);
    decayed_folded_efromB->Draw("E1 SAME");
    
    TLegend* leg4 = new TLegend(0.4,0.7,0.92,0.92);
    leg4->SetBorderSize(0);
    leg4->SetTextSize(0.035);
    leg4->AddEntry(dataBtoe,"Data e from B pT","l");
    leg4->AddEntry(decayed_efromB_selected_rebin,"Data e: Unfold---Data e pT--->Matrix Mult","l");
    leg4->AddEntry(decayed_folded_efromB,"Data e: Unfold---Data e pT--->Fold","l");
    leg4->Draw();
    // c4->SaveAs("decayed_folded_efromB.eps");
    c4->SaveAs("decayed_folded_efromB.png");
    
    //Canvas set C: result for D0 and Jpsi
    TCanvas*c5=new TCanvas("c5","c5",1200,600);
    c5->Divide(2,1,0.0001,0.0001);
    //c5_1: D0
    c5->cd(1);
    gPad->SetLogy();
    
    decayed_D0fromB->SetLineStyle(4);
    decayed_D0fromB->SetMarkerStyle(4);
    decayed_D0fromB->SetLineColor(4);
    decayed_D0fromB->SetMarkerColor(4);
    decayed_D0fromB->SetLineWidth(2);
    decayed_D0fromB->SetMarkerSize(2);
    decayed_D0fromB->Draw("E1");
    
    decayed_folded_D0fromB->SetLineStyle(6);
    decayed_folded_D0fromB->SetMarkerStyle(6);
    decayed_folded_D0fromB->SetLineColor(6);
    decayed_folded_D0fromB->SetMarkerColor(6);
    decayed_folded_D0fromB->SetLineWidth(1);
    decayed_folded_D0fromB->SetMarkerSize(1);
    decayed_folded_D0fromB->Draw("E1 SAME");
    
    TLegend* leg5_1 = new TLegend(0.4,0.7,0.92,0.92);
    leg5_1->SetBorderSize(0);
    leg5_1->SetTextSize(0.035);
    leg5_1->AddEntry(decayed_D0fromB,"D0 from B (with Matrix Mult)","l");
    leg5_1->AddEntry(decayed_folded_D0fromB,"D0 from B (with fold)","l");    
    leg5_1->Draw();
    
    //c5_2: Jpsi
    c5->cd(2);
    gPad->SetLogy();
    
    decayed_JpsifromB->SetLineStyle(4);
    decayed_JpsifromB->SetMarkerStyle(4);
    decayed_JpsifromB->SetLineColor(4);
    decayed_JpsifromB->SetMarkerColor(4);
    decayed_JpsifromB->SetLineWidth(2);
    decayed_JpsifromB->SetMarkerSize(2);
    decayed_JpsifromB->Draw("E1");
    
    decayed_folded_JpsifromB->SetLineStyle(6);
    decayed_folded_JpsifromB->SetMarkerStyle(6);
    decayed_folded_JpsifromB->SetLineColor(6);
    decayed_folded_JpsifromB->SetMarkerColor(6);
    decayed_folded_JpsifromB->SetLineWidth(1);
    decayed_folded_JpsifromB->SetMarkerSize(1);
    decayed_folded_JpsifromB->Draw("E1 SAME");
    
    TLegend* leg5_2 = new TLegend(0.4,0.7,0.92,0.92);
    leg5_2->SetBorderSize(0);
    leg5_2->SetTextSize(0.035);
    leg5_2->AddEntry(decayed_JpsifromB,"Jpsi from B (with Matrix Mult)","l");
    leg5_2->AddEntry(decayed_folded_JpsifromB,"Jpsi from B (with fold)","l");
    leg5_2->Draw();
    // c5->SaveAs("decayed_folded_DandJpsifromB.eps");
    c5->SaveAs("decayed_folded_DandJpsifromB.png");
    

    TCanvas*c6=new TCanvas("c6","c6",800,600);
    gPad->SetLogy();
    
    unfoldedBtoAll->SetLineStyle(4);
    unfoldedBtoAll->SetMarkerStyle(4);
    unfoldedBtoAll->SetLineColor(4);
    unfoldedBtoAll->SetMarkerColor(4);
    unfoldedBtoAll->SetLineWidth(1);
    unfoldedBtoAll->SetMarkerSize(1);
    unfoldedBtoAll->GetXaxis()->SetRangeUser(0,20);
    unfoldedBtoAll->SetMaximum(unfoldedBtoAll->GetMaximum()*10);
    unfoldedBtoAll->Draw("E1");
    unfoldedBtoAll->GetXaxis()->SetLabelFont(42);
    unfoldedBtoAll->GetYaxis()->SetLabelFont(42);
    unfoldedBtoAll->GetXaxis()->SetTitleFont(42);
    unfoldedBtoAll->GetYaxis()->SetTitleFont(42);
    unfoldedBtoAll->GetXaxis()->SetTitle("B p_{T}(GeV/c)");
    unfoldedBtoAll->GetYaxis()->SetTitle("(1/2#pip_{T})d#sigma/dp_{T}[mb GeV/c)^{-2}]");

    gPhenixB->SetLineStyle(2);
    gPhenixB->SetMarkerStyle(33);
    gPhenixB->SetLineColor(8);
    gPhenixB->SetMarkerColor(8);
    gPhenixB->SetLineWidth(1);
    gPhenixB->SetMarkerSize(1.3);
    gPhenixB->Draw("psame");
    
    TLegend* leg6 = new TLegend(0.4,0.7,0.92,0.92);
    leg6->SetBorderSize(0);
    leg6->SetTextSize(0.035);
    leg6->AddEntry(unfoldedBtoAll,"Unfolded B pT (from STAR e pT)","l");
    leg6->AddEntry(gPhenixB,"Unfolded B pT (phenix)","l");
    leg6->AddEntry(gPhenixB,"arXiv:1901.08405v1","");
    leg6->Draw();
    
    // c6->SaveAs("unfoldedBtoAll.eps");
    c6->SaveAs("unfoldedBtoAll.png");

    TFile* out= new TFile("etoBtoD0andJpsi.root","recreate");
    out->cd();
    c1->Write();
    c2->Write();
    c3->Write();
    c4->Write();
    c5->Write();
    c6->Write();
    decayed_D0fromB->Write();
    decayed_folded_D0fromB->Write();
    decayed_JpsifromB->Write();
    decayed_folded_JpsifromB->Write();
    unfoldedBtoAll->Write();
    out->Close();
}
