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

double calcRelativeDiff(TH1 *h1, TH1 *h2){//h1 base line
  double relativeDiff = 1.;
  int NoEmptyBin = 0.;
    for(int i=0;i<h1->GetNbinsX();++i){
        double x1 = h1->GetBinCenter(i+1);
        double y1 = h1->GetBinContent(i+1);
        double ey1 = h1->GetBinError(i+1);
        double x2 = h2->GetBinCenter(i+1);
        double y2 = h2->GetBinContent(i+1);
        double ey2 = h2->GetBinError(i+1);
        if(y2==0) {
            NoEmptyBin++;
            continue;
        }
        if(fabs(x1-x2)>1e-2) return -1;
        double tmp = fabs((y2-y1)/y1);
        relativeDiff= relativeDiff * tmp;
    }
    return (TMath::Power(relativeDiff,1./(h1->GetNbinsX()-NoEmptyBin)));
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
    TFile*fin0 = new TFile("out/d0_pp_PRD_update.root");
    TH1D* dataD0= (TH1D*)fin0->Get("pp_Dnew_err");
    
    
    TFile* fin = new TFile("out/DataInput_star.root");//update to phenix input
    TGraphErrors* gdataDtoe=(TGraphErrors*)fin->Get("gCountsCtoe_star"); //this is TGraph,needs to be converted into hist
    
    TFile* finphe =new TFile("out/DataInput_phenix.root");
    TH1D* gPhenixD=(TH1D*)finphe->Get("gYieldC_phenix"); //Phenix group unfold B result

    TFile *fsim = new TFile("rootfiles/DtoAll.star.hist.root");
    TH1D* pthist_Dtoall = (TH1D*)fsim->Get("pthist_Dtoall");
    TH1D* pthist_D0toall = (TH1D*)fsim->Get("pthist_D0toall");
    
    TH1D* hist_efromD_selected_rebin = (TH1D*)fsim->Get("hist_efromD_selected_rebin");
    
    TH1D* dataDtoe = hist_efromD_selected_rebin->Clone("dataDtoe");
    dataDtoe->Reset("ICES");
    double *ey;
    ey = gdataDtoe->GetEY();
    for(int i=0;i<dataDtoe->GetNbinsX();i++)
    {
        double x,y;
        gdataDtoe->GetPoint(i,x,y);
        dataDtoe->SetBinContent(i+1,y);
        dataDtoe->SetBinError(i+1,ey[i]);
    }
    
    //Now it's just the Btoe real value, we are going to convert it into scale factor
    
    TH1D* Dscale_e =(TH1D*)fsim->Get("Dscale_e");
    TH1D* D0scale_e =(TH1D*)fsim->Get("D0scale_e");

    
    TH2D* Mresponse_D_e_rebin = (TH2D*)fsim->Get("Mresponse_D_e_rebin");
    TH2D* Mresponse_e_rebin_D = (TH2D*)fsim->Get("Mresponse_e_rebin_D");
    TH2D* Mresponse_D0_e_rebin = (TH2D*)fsim->Get("Mresponse_D0_e_rebin");
    TH2D* Mresponse_e_rebin_D0 = (TH2D*)fsim->Get("Mresponse_e_rebin_D0");
    
    //calculate scale factor

    Dscale_e = calculatescale(Dscale_e, Mresponse_e_rebin_D);
    D0scale_e = calculatescale(D0scale_e, Mresponse_e_rebin_D0);
    //calculate Xovere scale factor
    TH1D* Dscale_eovere=(TH1D*)Dscale_e->Clone("Dscale_eovere");
    
    Dscale_eovere->Reset("ICES");
    for(int i=0;i<Dscale_eovere->GetNbinsX();i++){
        Dscale_eovere->SetBinContent(i+1,1./Dscale_e->GetBinContent(i+1));
    }
    
    //create container for smeared D
    TH1D* D_smeared_e = (TH1D*)pthist_Dtoall->Clone("D_smeared_e");
    TH1D* D0_smeared_e = (TH1D*)pthist_D0toall->Clone("D0_smeared_e");
    
    //smear the B spectra (this is for RooUnfoldResponse constructing)
    D_smeared_e = Bscale(D_smeared_e,Dscale_e);
    D0_smeared_e = Bscale(D0_smeared_e,D0scale_e);
    
    //construct the correct RooUnfoldResponse, but the unfolded B is the smeared B, in the end we have to make the output antiscale; for fold we need to input the smeared/scaled B spectrum
    RooUnfoldResponse* response_D_e_rebin = new RooUnfoldResponse(hist_efromD_selected_rebin,D_smeared_e,Mresponse_D_e_rebin,"response_D_e_rebin","response_D_e_rebin");
    RooUnfoldResponse* response_D0_e_rebin = new RooUnfoldResponse(hist_efromD_selected_rebin,D0_smeared_e,Mresponse_D0_e_rebin,"response_D0_e_rebin","response_D0_e_rebin");
    RooUnfoldResponse* response_e_rebin_D0 = new RooUnfoldResponse(D0_smeared_e,hist_efromD_selected_rebin,Mresponse_e_rebin_D0,"response_e_rebin_D0","response_e_rebin_D0");

    TH2D* decay_Dtoe_selected_rebin = (TH2D*)fsim->Get("decay_Dtoe_selected_rebin");
    
    TH2D* EQdecaymatrixe = new TH2D("EQdecaymatrixe","EQdecaymatrixe",pthist_Dtoall->GetNbinsX(),pthist_Dtoall->GetBinLowEdge(1),pthist_Dtoall->GetBinLowEdge(pthist_Dtoall->GetNbinsX()+1),pthist_Dtoall->GetNbinsX(),pthist_Dtoall->GetBinLowEdge(1),pthist_Dtoall->GetBinLowEdge(pthist_Dtoall->GetNbinsX()+1));
    
    for(int i =0; i<pthist_Dtoall->GetNbinsX();i++)
        EQdecaymatrixe->SetBinContent(i+1,i+1,1./Dscale_e->GetBinContent(i+1));
        
        
        
    TH2D* EQdecaymatrixf = new TH2D("EQdecaymatrixf","EQdecaymatrixf",pthist_D0toall->GetNbinsX(),pthist_D0toall->GetBinLowEdge(1),pthist_D0toall->GetBinLowEdge(pthist_D0toall->GetNbinsX()+1),pthist_D0toall->GetNbinsX(),pthist_D0toall->GetBinLowEdge(1),pthist_D0toall->GetBinLowEdge(pthist_D0toall->GetNbinsX()+1));
    
    for(int i =0; i<pthist_D0toall->GetNbinsX();i++)
        EQdecaymatrixf->SetBinContent(i+1,i+1,1./D0scale_e->GetBinContent(i+1));
    
    //decay matrix normalization
    decay_Dtoe_selected_rebin=(TH2D*)DecayMatrixNormalization(decay_Dtoe_selected_rebin,pthist_Dtoall);
    
    //calculate the equivalent decay matrix for D0 and Jpsi
    decay_Dtoe_selected_rebin=equivalentdecaymatrix(decay_Dtoe_selected_rebin,Dscale_eovere);
    
    //container for the reinversed
    TH1D* reinversed_efromD_selected_rebin = (TH1D*)hist_efromD_selected_rebin->Clone("reinversed_efromD_selected_rebin");
    TH1D* reinversed_folded_efromD_viaD0 = (TH1D*)hist_efromD_selected_rebin->Clone("reinversed_folded_efromD_viaD0");
    
    //container for the decayed
    TH1D* decayed_efromD_selected_rebin = (TH1D*)hist_efromD_selected_rebin->Clone("decayed_efromD_selected_rebin");
    TH1D* decayed_folded_efromD_viaD0 = (TH1D*)hist_efromD_selected_rebin->Clone("decayed_folded_efromD_viaD0");
    
    //reset those containers
    reinversed_efromD_selected_rebin->Reset("ICES");
    reinversed_folded_efromD_viaD0->Reset("ICES");
    
    decayed_efromD_selected_rebin->Reset("ICES");
    decayed_folded_efromD_viaD0->Reset("ICES");
    
    int nIter = 215;//100;//6;
    
    
    
    double diff = 1.;
    while(diff>0.01)
    {
    //unfold==================================================================================================
    //unfold->inverse for consistency check
    //Dtoe -> D
    RooUnfoldBayes inverse_DtoAll(response_D_e_rebin,hist_efromD_selected_rebin,1);
    TH1D* inversedDtoAll_smear_e = (TH1D*)((TH1D*)(inverse_DtoAll.Hreco(2)))->Clone("inversedDtoAll_smear_e");
    TMatrixD*inversedDtoAll_smear_ecov = (TMatrixD*)(inverse_DtoAll.Ereco(2));
    
    //unsmear
    TH1D* inversedDtoAll=(TH1D*)inversedDtoAll_smear_e->Clone("inversedDtoAll");
    inversedDtoAll->Reset("ICES");
    inversedDtoAll=Decay(EQdecaymatrixe,inversedDtoAll_smear_e,inversedDtoAll_smear_ecov,inversedDtoAll);
    
    //unfold for data
    //unfold to D
    RooUnfoldBayes unfold_DtoAll(response_D_e_rebin,dataDtoe,nIter);
    TH1D* unfoldedDtoAll_smear_e = (TH1D*)((TH1D*)unfold_DtoAll.Hreco(2))->Clone("unfoldedDtoAll_smear_e");
    TMatrixD*unfoldedDtoAll_smear_ecov = (TMatrixD*)(unfold_DtoAll.Ereco(2));
    
    //unsmear
    TH1D* unfoldedDtoAll=(TH1D*)unfoldedDtoAll_smear_e->Clone("unfoldedDtoAll");
    unfoldedDtoAll->Reset("ICES");
    unfoldedDtoAll=Decay(EQdecaymatrixe,unfoldedDtoAll_smear_e,unfoldedDtoAll_smear_ecov,unfoldedDtoAll);
    
    //reinverse for consistency check 
    reinversed_efromD_selected_rebin = (TH1D*)Decay(decay_Dtoe_selected_rebin,inversedDtoAll_smear_e,inversedDtoAll_smear_ecov,reinversed_efromD_selected_rebin);
    
    //decay the data to check convergence========================================================================================
    decayed_efromD_selected_rebin = (TH1D*)Decay(decay_Dtoe_selected_rebin,unfoldedDtoAll_smear_e,unfoldedDtoAll_smear_ecov,decayed_efromD_selected_rebin);
    
    diff = calcRelativeDiff(dataDtoe, decayed_efromD_selected_rebin);
    nIter++;
    }
    nIter--;
    cout << "nIter = " << nIter << " ; diff = " << diff << "for D" << endl;
    
    //independent section for D0
    nIter = 1;
    diff = 1.;
    double olddiff = 1.;
    while(diff>0.01)
    {
    //unfold==================================================================================================
    //unfold->inverse for consistency check
    //Dtoe ->D0
    RooUnfoldBayes inverse_D0toAll(response_D0_e_rebin,hist_efromD_selected_rebin,1);
    TH1D* inversedD0toAll_smear_e = (TH1D*)((TH1D*)(inverse_D0toAll.Hreco(2)))->Clone("inversedD0toAll_smear_e");
    TMatrixD*inversedD0toAll_smear_ecov = (TMatrixD*)(inverse_D0toAll.Ereco(2));
    
    //unsmear
    TH1D* inversedD0toAll=(TH1D*)inversedD0toAll_smear_e->Clone("inversedD0toAll");
    inversedD0toAll->Reset("ICES");
    inversedD0toAll=Decay(EQdecaymatrixf,inversedD0toAll_smear_e,inversedD0toAll_smear_ecov,inversedD0toAll);
    
    //unfold for data
    //unfold to D0
    RooUnfoldBayes unfold_D0toAll(response_D0_e_rebin,dataDtoe,nIter);
    TH1D* unfoldedD0toAll_smear_e = (TH1D*)((TH1D*)unfold_D0toAll.Hreco(2))->Clone("unfoldedD0toAll_smear_e");
    TMatrixD*unfoldedD0toAll_smear_ecov = (TMatrixD*)(unfold_D0toAll.Ereco(2));
    
    //unsmear
    TH1D* unfoldedD0toAll=(TH1D*)unfoldedD0toAll_smear_e->Clone("unfoldedD0toAll");
    unfoldedD0toAll->Reset("ICES");
    unfoldedD0toAll=Decay(EQdecaymatrixf,unfoldedD0toAll_smear_e,unfoldedD0toAll_smear_ecov,unfoldedD0toAll);
    
    //reinverse the for consistency check
    RooUnfoldBayes reinverse_fold_efromD_viaD0(response_e_rebin_D0,inversedD0toAll_smear_e,1);
    reinversed_folded_efromD_viaD0 = (TH1D*)((TH1D*)reinverse_fold_efromD_viaD0.Hreco(2))->Clone("reinversed_folded_efromD_viaD0");
    
    //fold the unfolded D0 to check convergence
    RooUnfoldBayes decay_fold_efromD_viaD0(response_e_rebin_D0,unfoldedD0toAll_smear_e,1);
    decayed_folded_efromD_viaD0 = (TH1D*)((TH1D*)decay_fold_efromD_viaD0.Hreco(2))->Clone("decayed_folded_efromD_viaD0");
    
    diff = calcRelativeDiff(dataDtoe, decayed_folded_efromD_viaD0);
    if(diff>=olddiff)
        break;
    olddiff=diff;    
    nIter++;
    }
    
    cout << "nIter = " << nIter << " ; diff = " << olddiff << endl;
    
    
    {
    //unfold==================================================================================================
    //unfold->inverse for consistency check
    //Dtoe ->D0
    RooUnfoldBayes inverse_D0toAll(response_D0_e_rebin,hist_efromD_selected_rebin,1);
    TH1D* inversedD0toAll_smear_e = (TH1D*)((TH1D*)(inverse_D0toAll.Hreco(2)))->Clone("inversedD0toAll_smear_e");
    TMatrixD*inversedD0toAll_smear_ecov = (TMatrixD*)(inverse_D0toAll.Ereco(2));
    
    //unsmear
    TH1D* inversedD0toAll=(TH1D*)inversedD0toAll_smear_e->Clone("inversedD0toAll");
    inversedD0toAll->Reset("ICES");
    inversedD0toAll=Decay(EQdecaymatrixf,inversedD0toAll_smear_e,inversedD0toAll_smear_ecov,inversedD0toAll);
    
    //unfold for data
    //unfold to D0
    RooUnfoldBayes unfold_D0toAll(response_D0_e_rebin,dataDtoe,nIter);
    TH1D* unfoldedD0toAll_smear_e = (TH1D*)((TH1D*)unfold_D0toAll.Hreco(2))->Clone("unfoldedD0toAll_smear_e");
    TMatrixD*unfoldedD0toAll_smear_ecov = (TMatrixD*)(unfold_D0toAll.Ereco(2));
    
    //unsmear
    TH1D* unfoldedD0toAll=(TH1D*)unfoldedD0toAll_smear_e->Clone("unfoldedD0toAll");
    unfoldedD0toAll->Reset("ICES");
    unfoldedD0toAll=Decay(EQdecaymatrixf,unfoldedD0toAll_smear_e,unfoldedD0toAll_smear_ecov,unfoldedD0toAll);
    
    //reinverse the for consistency check
    RooUnfoldBayes reinverse_fold_efromD_viaD0(response_e_rebin_D0,inversedD0toAll_smear_e,1);
    reinversed_folded_efromD_viaD0 = (TH1D*)((TH1D*)reinverse_fold_efromD_viaD0.Hreco(2))->Clone("reinversed_folded_efromD_viaD0");
    
    //fold the unfolded D0 to check convergence
    RooUnfoldBayes decay_fold_efromD_viaD0(response_e_rebin_D0,unfoldedD0toAll_smear_e,1);
    decayed_folded_efromD_viaD0 = (TH1D*)((TH1D*)decay_fold_efromD_viaD0.Hreco(2))->Clone("decayed_folded_efromD_viaD0");
    }

    
    
    
    //===========
    //cv1
    calcInvYield(pthist_Dtoall);
    calcInvYield(inversedDtoAll);
    
    calcInvYield(pthist_D0toall);
    calcInvYield(inversedD0toAll);
    
    //cv2
    calcInvYield(hist_efromD_selected_rebin);
    calcInvYield(reinversed_efromD_selected_rebin);
    calcInvYield(reinversed_folded_efromD_viaD0);
    
    //cv3
    calcInvYield(dataDtoe);
    calcInvYield(decayed_efromD_selected_rebin);
    calcInvYield(decayed_folded_efromD_viaD0);
    
    //cv4_1
    calcInvYield(unfoldedDtoAll);
    
    //cv4_2
    calcInvYield(unfoldedD0toAll);


    //
    //===========
    //canvas set A: simu-self check
    //c1: Step 1: unfolded B and direct statistics of D
    TCanvas*c1=new TCanvas("c1","c1",1200,600);
    c1->Divide(2,1,0.0001,0.0001);
    //c1_1 D
    c1->cd(1);
    gPad->SetLogy();

    pthist_Dtoall->SetLineStyle(1);
    pthist_Dtoall->SetMarkerStyle(1);
    pthist_Dtoall->SetLineColor(1);
    pthist_Dtoall->SetMarkerColor(1);
    pthist_Dtoall->SetLineWidth(2);
    pthist_Dtoall->SetMarkerSize(2);
    pthist_Dtoall->Draw("E1");
    
    inversedDtoAll->SetLineStyle(2);
    inversedDtoAll->SetMarkerStyle(2);
    inversedDtoAll->SetLineColor(2);
    inversedDtoAll->SetMarkerColor(2);
    inversedDtoAll->SetLineWidth(1);
    inversedDtoAll->SetMarkerSize(1);
    inversedDtoAll->Draw("E1 SAME");

    TLegend* leg1_1 = new TLegend(0.4,0.7,0.92,0.92);
    leg1_1->SetBorderSize(0);
    leg1_1->SetTextSize(0.035);
    leg1_1->AddEntry(pthist_Dtoall,"SIMU D pT","l");
    leg1_1->AddEntry(inversedDtoAll,"Unfolded SIMU D pT (from e pT)","l");
    leg1_1->Draw();
    
    //c1_2 D0
    c1->cd(2);
    gPad->SetLogy();

    pthist_D0toall->SetLineStyle(1);
    pthist_D0toall->SetMarkerStyle(1);
    pthist_D0toall->SetLineColor(1);
    pthist_D0toall->SetMarkerColor(1);
    pthist_D0toall->SetLineWidth(2);
    pthist_D0toall->SetMarkerSize(2);
    pthist_D0toall->Draw("E1");
    
    inversedD0toAll->SetLineStyle(2);
    inversedD0toAll->SetMarkerStyle(2);
    inversedD0toAll->SetLineColor(2);
    inversedD0toAll->SetMarkerColor(2);
    inversedD0toAll->SetLineWidth(1);
    inversedD0toAll->SetMarkerSize(1);
    inversedD0toAll->Draw("E1 SAME");

    TLegend* leg1_2 = new TLegend(0.4,0.7,0.92,0.92);
    leg1_2->SetBorderSize(0);
    leg1_2->SetTextSize(0.035);
    leg1_2->AddEntry(pthist_Dtoall,"SIMU D0 pT","l");
    leg1_2->AddEntry(inversedDtoAll,"Unfolded SIMU D0 pT (from e pT)","l");
    leg1_2->Draw();
    
    // c1->SaveAs("inversedBtoAll.eps");
    //c1->SaveAs("inversedBtoAll.png");

    //c2: Step 2: folded/decayed e pT by using unfolded SIMU D pT
    TCanvas*c2=new TCanvas("c2","c2",1200,600);
    c2->Divide(2,1,0.0001,0.0001);
    //c2_1: D
    c2->cd(1);
    gPad->SetLogy();

    hist_efromD_selected_rebin->SetLineStyle(1);
    hist_efromD_selected_rebin->SetMarkerStyle(1);
    hist_efromD_selected_rebin->SetLineColor(1);
    hist_efromD_selected_rebin->SetMarkerColor(1);
    hist_efromD_selected_rebin->SetLineWidth(3);
    hist_efromD_selected_rebin->SetMarkerSize(3);
    hist_efromD_selected_rebin->Draw("E1");
    
    reinversed_efromD_selected_rebin->SetLineStyle(2);
    reinversed_efromD_selected_rebin->SetMarkerStyle(2);
    reinversed_efromD_selected_rebin->SetLineColor(2);
    reinversed_efromD_selected_rebin->SetMarkerColor(2);
    reinversed_efromD_selected_rebin->SetLineWidth(2);
    reinversed_efromD_selected_rebin->SetMarkerSize(2);
    reinversed_efromD_selected_rebin->Draw("E1 SAME");
    
    TLegend* leg2_1 = new TLegend(0.4,0.7,0.92,0.92);
    leg2_1->SetBorderSize(0);
    leg2_1->SetTextSize(0.035);
    leg2_1->AddEntry(hist_efromD_selected_rebin,"SIMU e from B pT","l");
    leg2_1->AddEntry(reinversed_efromD_selected_rebin,"SIMU e: Unfold---SIMU e pT--->Matrix Mult","l");
    leg2_1->Draw();
    
    //c2_2: D0
    c2->cd(2);
    gPad->SetLogy();

    hist_efromD_selected_rebin->SetLineStyle(1);
    hist_efromD_selected_rebin->SetMarkerStyle(1);
    hist_efromD_selected_rebin->SetLineColor(1);
    hist_efromD_selected_rebin->SetMarkerColor(1);
    hist_efromD_selected_rebin->SetLineWidth(3);
    hist_efromD_selected_rebin->SetMarkerSize(3);
    hist_efromD_selected_rebin->Draw("E1");
    
    reinversed_folded_efromD_viaD0->SetLineStyle(2);
    reinversed_folded_efromD_viaD0->SetMarkerStyle(2);
    reinversed_folded_efromD_viaD0->SetLineColor(2);
    reinversed_folded_efromD_viaD0->SetMarkerColor(2);
    reinversed_folded_efromD_viaD0->SetLineWidth(2);
    reinversed_folded_efromD_viaD0->SetMarkerSize(2);
    reinversed_folded_efromD_viaD0->Draw("E1 SAME");
    
    TLegend* leg2_2 = new TLegend(0.4,0.7,0.92,0.92);
    leg2_2->SetBorderSize(0);
    leg2_2->SetTextSize(0.035);
    leg2_2->AddEntry(hist_efromD_selected_rebin,"SIMU e from B pT","l");
    leg2_2->AddEntry(reinversed_folded_efromD_viaD0,"(SIMU e unfolded D0) folded e","l");
    leg2_2->Draw();
    // c2->SaveAs("reinversed_folded_efromB.eps");
    //c2->SaveAs("reinversed_folded_efromB.png");
    
    //Canvas set B: data self check via e pT
    TCanvas*c3=new TCanvas("c3","c3",1200,600);
    c3->Divide(2,1,0.0001,0.0001);
    //c3_1: D
    c3->cd(1);
    gPad->SetLogy();

    dataDtoe->GetXaxis()->SetTitle("p_{T}(GeV/c)");
    dataDtoe->GetYaxis()->SetTitle("Ed^{3}#sigma/dp^{3}(mbGeV^{-2}c^{3})");
    dataDtoe->GetXaxis()->SetLabelFont(42);
    dataDtoe->GetYaxis()->SetLabelFont(42);
    dataDtoe->GetXaxis()->SetTitleFont(42);
    dataDtoe->GetYaxis()->SetTitleFont(42);
    dataDtoe->SetLineStyle(1);
    dataDtoe->SetMarkerStyle(1);
    dataDtoe->SetLineColor(1);
    dataDtoe->SetMarkerColor(1);
    dataDtoe->SetLineWidth(3);
    dataDtoe->SetMarkerSize(3);
    dataDtoe->Draw("E1");
    
    decayed_efromD_selected_rebin->SetLineStyle(2);
    decayed_efromD_selected_rebin->SetMarkerStyle(2);
    decayed_efromD_selected_rebin->SetLineColor(2);
    decayed_efromD_selected_rebin->SetMarkerColor(2);
    decayed_efromD_selected_rebin->SetLineWidth(2);
    decayed_efromD_selected_rebin->SetMarkerSize(2);
    decayed_efromD_selected_rebin->Draw("E1 SAME");
    
    TLegend* leg3_1 = new TLegend(0.4,0.7,0.92,0.92);
    leg3_1->SetBorderSize(0);
    leg3_1->SetTextSize(0.035);
    leg3_1->AddEntry(dataDtoe,"Data e from D pT","l");
    leg3_1->AddEntry(decayed_efromD_selected_rebin,"Data e: Unfold---Data e pT--->Matrix Mult","l");
    leg3_1->Draw();
    //c3_2: D0
    c3->cd(2);
    gPad->SetLogy();

    dataDtoe->GetXaxis()->SetTitle("p_{T}(GeV/c)");
    dataDtoe->GetYaxis()->SetTitle("Ed^{3}#sigma/dp^{3}(mbGeV^{-2}c^{3})");
    dataDtoe->GetXaxis()->SetLabelFont(42);
    dataDtoe->GetYaxis()->SetLabelFont(42);
    dataDtoe->GetXaxis()->SetTitleFont(42);
    dataDtoe->GetYaxis()->SetTitleFont(42);
    dataDtoe->SetLineStyle(1);
    dataDtoe->SetMarkerStyle(1);
    dataDtoe->SetLineColor(1);
    dataDtoe->SetMarkerColor(1);
    dataDtoe->SetLineWidth(3);
    dataDtoe->SetMarkerSize(3);
    dataDtoe->Draw("E1");
    
    decayed_folded_efromD_viaD0->SetLineStyle(2);
    decayed_folded_efromD_viaD0->SetMarkerStyle(2);
    decayed_folded_efromD_viaD0->SetLineColor(2);
    decayed_folded_efromD_viaD0->SetMarkerColor(2);
    decayed_folded_efromD_viaD0->SetLineWidth(2);
    decayed_folded_efromD_viaD0->SetMarkerSize(2);
    decayed_folded_efromD_viaD0->Draw("E1 SAME");
    
    TLegend* leg3_2 = new TLegend(0.4,0.7,0.92,0.92);
    leg3_2->SetBorderSize(0);
    leg3_2->SetTextSize(0.035);
    leg3_2->AddEntry(dataDtoe,"Data e from D pT","l");
    leg3_2->AddEntry(decayed_folded_efromD_viaD0,"Unfolded D0 folded e from D","l");
    leg3_2->Draw();
    // c3->SaveAs("decayed_folded_efromB.eps");
    //c3->SaveAs("decayed_folded_efromB.png");
    
    //Canvas set C: unfolded D and unfolded D0 results
    TCanvas*c4=new TCanvas("c4","c4",1200,600);
    c4->Divide(2,1,0.0001,0.0001);
    //c4_1:D
    c4->cd(1);
    gPad->SetLogy();
    
    unfoldedDtoAll->SetLineStyle(4);
    unfoldedDtoAll->SetMarkerStyle(4);
    unfoldedDtoAll->SetLineColor(4);
    unfoldedDtoAll->SetMarkerColor(4);
    unfoldedDtoAll->SetLineWidth(1);
    unfoldedDtoAll->SetMarkerSize(1);
    unfoldedDtoAll->GetXaxis()->SetRangeUser(0,20);
    //unfoldedDtoAll->SetMaximum(unfoldedDtoAll->GetMaximum()*10);
    unfoldedDtoAll->Draw("E1");
    unfoldedDtoAll->GetXaxis()->SetLabelFont(42);
    unfoldedDtoAll->GetYaxis()->SetLabelFont(42);
    unfoldedDtoAll->GetXaxis()->SetTitleFont(42);
    unfoldedDtoAll->GetYaxis()->SetTitleFont(42);
    unfoldedDtoAll->GetXaxis()->SetTitle("D p_{T}(GeV/c)");
    unfoldedDtoAll->GetYaxis()->SetTitle("(1/2#pip_{T})d#sigma/dp_{T}[mb GeV/c)^{-2}]");

    gPhenixD->SetLineStyle(2);
    gPhenixD->SetMarkerStyle(33);
    gPhenixD->SetLineColor(8);
    gPhenixD->SetMarkerColor(8);
    gPhenixD->SetLineWidth(1);
    gPhenixD->SetMarkerSize(1.3);
    gPhenixD->Draw("psame");
    
    TLegend* leg4_1 = new TLegend(0.4,0.7,0.92,0.92);
    leg4_1->SetBorderSize(0);
    leg4_1->SetTextSize(0.035);
    leg4_1->AddEntry(unfoldedDtoAll,"Unfolded D pT","l");
    leg4_1->AddEntry(gPhenixD,"Unfolded D pT (phenix group)","l");
    leg4_1->AddEntry(gPhenixD,"arXiv:1901.08405v1","");
    leg4_1->Draw();
    
    //c4_2:D0
    c4->cd(2);
    gPad->SetLogy();
    
    unfoldedD0toAll->Scale(0.5);//divided by dY=1-(-1)=2
    unfoldedD0toAll->SetLineStyle(4);
    unfoldedD0toAll->SetMarkerStyle(4);
    unfoldedD0toAll->SetLineColor(4);
    unfoldedD0toAll->SetMarkerColor(4);
    unfoldedD0toAll->SetLineWidth(1);
    unfoldedD0toAll->SetMarkerSize(1);
    unfoldedD0toAll->GetXaxis()->SetRangeUser(0,20);
    //unfoldedD0toAll->SetMaximum(unfoldedD0toAll->GetMaximum()*10);
    unfoldedD0toAll->Draw("E1");
    unfoldedD0toAll->GetXaxis()->SetLabelFont(42);
    unfoldedD0toAll->GetYaxis()->SetLabelFont(42);
    unfoldedD0toAll->GetXaxis()->SetTitleFont(42);
    unfoldedD0toAll->GetYaxis()->SetTitleFont(42);
    unfoldedD0toAll->GetXaxis()->SetTitle("D0 p_{T}(GeV/c)");
    unfoldedD0toAll->GetYaxis()->SetTitle("(1/2#pip_{T})d#sigma/dp_{T}[mb GeV/c)^{-2}]");

    dataD0->SetLineStyle(2);
    dataD0->SetMarkerStyle(33);
    dataD0->SetLineColor(8);
    dataD0->SetMarkerColor(8);
    dataD0->SetLineWidth(1);
    dataD0->SetMarkerSize(1.3);
    dataD0->Draw("psame");
    
    TLegend* leg4_2 = new TLegend(0.4,0.7,0.92,0.92);
    leg4_2->SetBorderSize(0);
    leg4_2->SetTextSize(0.035);
    leg4_2->AddEntry(unfoldedDtoAll,"Unfolded D0 pT","l");
    leg4_2->AddEntry(dataD0,"Data D0 pT (phenix group)","l");
    leg4_2->Draw();
    
    // c4->SaveAs("unfoldedBtoAll.eps");
    //c4->SaveAs("unfoldedBtoAll.png");

    TFile* out= new TFile("etoDandD0.phenix.root","recreate");
    out->cd();
    c1->Write();
    c2->Write();
    c3->Write();
    c4->Write();
    unfoldedDtoAll->Write();
    unfoldedD0toAll->Write();
    out->Close();
}
