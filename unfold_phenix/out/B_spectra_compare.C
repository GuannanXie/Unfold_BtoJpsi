#include "math.h"
void B_spectra_compare()
{
    //gROOT->LoadMacro("../rootinit.C");
    gROOT->LoadMacro("../funcUtil.C");
    
    
    TFile*fstar = new TFile("etoBtoD0andJpsi.root");
    //TH1D*Bstar= (TH1D*)fstar->Get("unfoldedBtoAll");
    TH1D*Bstar= (TH1D*)fstar->Get("unfoldedBtoAll");
    TH1D*decayed_JpsifromB= (TH1D*)fstar->Get("decayed_JpsifromB");
    decayed_JpsifromB->Scale(0.0594/30.);
    
    for(int i=0;i<decayed_JpsifromB->GetNbinsX();i++)
    {
        decayed_JpsifromB->SetBinError(i+1,decayed_JpsifromB->GetBinError(i+1)/(decayed_JpsifromB->GetBinLowEdge(i+2)-decayed_JpsifromB->GetBinLowEdge(i+1))/decayed_JpsifromB->GetBinCenter(i+1));
        decayed_JpsifromB->SetBinContent(i+1,decayed_JpsifromB->GetBinContent(i+1)/(decayed_JpsifromB->GetBinLowEdge(i+2)-decayed_JpsifromB->GetBinLowEdge(i+1))/decayed_JpsifromB->GetBinCenter(i+1));
    }

    
    TFile*fphnx = new TFile("DataInput.root");
    TGraph* Bphnx = (TGraph*)fphnx->Get("PHENIXgraph");
    
    //calculate integral
    const int NoB = 60;
    double min= -0.25;
    double max = 29.75;
    TH1D*itgN=new TH1D("itgN","itgN",NoB,min,max);

    for(int i=0;i<NoB;i++)
    {
        double temp = 0.;
        for(int j=0 ; j<=i ;j++)
            temp=temp+Bstar->GetBinContent(j+1);
        itgN->SetBinContent(i+1,temp);
    }
    
    //fit the integral
    TCanvas*c1=new TCanvas("c1","c1",1000,1000);
    gPad->SetLogx();
    
    itgN->SetLineStyle(1);
    itgN->SetMarkerStyle(1);
    itgN->SetLineColor(1);
    itgN->SetMarkerColor(1);
    itgN->SetLineWidth(2);
    itgN->SetMarkerSize(2);
    itgN->Draw();
    
    TF1* myfun = new TF1("myfun","[0]-[0]/(1+[1]*pow(x,[2])+[3]*pow(x,[4]))",0.,29.5);
    myfun->SetParameters(2.06329e-04,1.09404e-01,1.56152e+00,1.97658e-03,4.56939e+00);
    itgN->Fit(myfun,"R S N");
    
    myfun->SetLineStyle(2);
    myfun->SetMarkerStyle(2);
    myfun->SetLineColor(2);
    myfun->SetMarkerColor(2);
    myfun->SetLineWidth(1);
    myfun->SetMarkerSize(1);
    myfun->Draw("SAME");
    
    //calculate the derivative and divided by pT
    TH1D*cBstar= new TH1D("cBstar","cBstar",NoB,min,max);
    for(int i=1;i<40;i++)
        cBstar->SetBinContent(i+1,myfun->Derivative(0.5*i)/0.5/i);
    
    
    TCanvas*c2=new TCanvas("c2","c2",1000,1000);
    gPad->SetLogy();
    cBstar->SetLineStyle(3);
    cBstar->SetMarkerStyle(3);
    cBstar->SetLineColor(3);
    cBstar->SetMarkerColor(3);
    cBstar->SetLineWidth(1);
    cBstar->SetMarkerSize(1);
    cBstar->Draw("P");
    
    Bphnx->SetLineStyle(2);
    Bphnx->SetMarkerStyle(2);
    Bphnx->SetLineColor(2);
    Bphnx->SetMarkerColor(2);
    Bphnx->SetLineWidth(1);
    Bphnx->SetMarkerSize(1);
    Bphnx->Draw("P SAME");
    
    TLegend* leg2 = new TLegend(0.4,0.7,0.92,0.92);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.035);
    leg2->AddEntry(cBstar,"STAR","l");
    leg2->AddEntry(Bphnx,"PHENIX","l");
    leg2->Draw();
    
    TFile*fnoll= new TFile("../rootfiles/fonll_Btojpsi_pp200.root");
    TGraphAsymmErrors *gBJpsiFonll = (TGraphAsymmErrors*)fnoll->Get("gBJpsiFonll");
    scaleGraph(gBJpsiFonll,1./30./1e6);
    TGraphAsymmErrors *gBJpsipp200 = (TGraphAsymmErrors*)fnoll->Get("gBJpsipp200");
    
    
    
    
    
    
    


    
/*    TF1* myfun = new TF1("myfun","[0]-[0]/(1+[1]*pow(x,[2])+[3]*pow(x,[4]))",0.,29.5);
    myfun->SetParameters(2.06329e-04,1.09404e-01,1.56152e+00,1.97658e-03,4.56939e+00);
    itgN->Fit(myfun,"R S N");
    
    myfun->SetLineStyle(2);
    myfun->SetMarkerStyle(2);
    myfun->SetLineColor(2);
    myfun->SetMarkerColor(2);
    myfun->SetLineWidth(1);
    myfun->SetMarkerSize(1);
    myfun->Draw("SAME");*/



























    
    
    
    
    TCanvas * c3 = new TCanvas("c3","c3",1000,800);
    gPad->SetLogy();
    
    gBJpsiFonll->SetFillColor(2);
    gBJpsiFonll->GetXaxis()->SetRangeUser(0.,10.);
    gBJpsiFonll->Draw("A4 0");
    
    gBJpsipp200->SetLineStyle(1);
    gBJpsipp200->SetMarkerStyle(1);
    gBJpsipp200->SetLineColor(1);
    gBJpsipp200->SetMarkerColor(1);
    gBJpsipp200->SetLineWidth(1);
    gBJpsipp200->SetMarkerSize(1);
    gBJpsipp200->Draw("P E1 SAME");
    

    decayed_JpsifromB->SetLineStyle(4);
    decayed_JpsifromB->SetMarkerStyle(4);
    decayed_JpsifromB->SetLineColor(4);
    decayed_JpsifromB->SetMarkerColor(4);
    decayed_JpsifromB->SetLineWidth(4);
    decayed_JpsifromB->SetMarkerSize(4);
    decayed_JpsifromB->Draw("E1 SAME");
    
    
    
}
