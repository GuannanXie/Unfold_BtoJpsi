//#include "math.h"
void Create()
{
    const int NoB=8;
    const double EoB[NoB+1]={3.,3.5,4.,4.5,5.,5.5,6.,8.,10.};
    double y[NoB]={1.556987e-007,7.859455e-008,4.223772e-008,2.429661e-008,1.030599e-008,6.863500e-009,1.033828e-009,3.149538e-010};
    double ey[NoB]={3.728824e-008,1.680373e-008,9.208601e-009,5.522727e-009,2.576287e-009,1.371894e-009,1.986033e-010,9.865579e-011};
    double x[NoB]={3.25,3.75,4.25,4.75,5.25,5.75,7.,9.};
    double ex[NoB]={0.};
    
    TH1D*gcrosssection=new TH1D("gcrosssection","gcrosssection",NoB,EoB);
    for(int i=0;i<NoB;i++)
    {
        gcrosssection->SetBinContent(i+1,y[i]);
        gcrosssection->SetBinError(i+1,ey[i]);
    }
    TF1* myf = new TF1("myf","[0]*exp( -x*[1] )",3.,10.);
    myf->SetNpx(1000);
    
    TFitResultPtr test= gcrosssection->Fit("myf","S N R");
    gcrosssection->Print("all");

    
    double CL = 0.995;
    TH1D *hint = new TH1D("hint", Form("Fitted with %3.3f conf.band",CL), 100, 3, 10);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint,CL);

    TH1D *hint_pt = (TH1D*)hint->Clone("hint_pt");
    for(int i=0;i<hint_pt->GetNbinsX();i++)
    {
      hint_pt->SetBinContent(i,hint_pt->GetBinContent(i)*hint_pt->GetBinCenter(i));
      hint_pt->SetBinError(i,hint_pt->GetBinError(i)*hint_pt->GetBinCenter(i));
    }
    
    double paratemp[2];
    myf->GetParameters(paratemp);
    
    const double para[2]={paratemp[0],paratemp[1]};
    
    TCanvas*c1= new TCanvas("c1","c1",1600,800);
    gPad->SetLogy();

    hint->SetStats(kFALSE);
    hint->SetFillColor(4);
    hint->Draw("e3");
    
    gcrosssection->SetLineStyle(1);
    gcrosssection->SetMarkerStyle(1);
    gcrosssection->SetLineColor(1);
    gcrosssection->SetMarkerColor(1);
    gcrosssection->SetLineWidth(2);
    gcrosssection->SetMarkerSize(2);
    gcrosssection->Draw("E1 same");

    myf->SetLineStyle(2);
    myf->SetMarkerStyle(2);
    myf->SetLineColor(2);
    myf->SetMarkerColor(2);
    myf->SetLineWidth(1);
    myf->SetMarkerSize(1);
    myf->Draw("E1 SAME");


    TLegend* leg1 = new TLegend(0.4,0.7,0.92,0.92);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.035);
    leg1->AddEntry(gcrosssection,"e From B E*d^{3}#sigma/dpT^{3} ","l");
    leg1->AddEntry(myf,"exponetial fit","l");
    leg1->Draw();
    
    TF1* myf_pT = new TF1("myf_pT","x*[0]*exp( -x*[1] )",3.,10.);
    myf_pT->SetNpx(1000);
    
    TH1D*gCountsBtoe_highpt=new TH1D("gCountsBtoe_highpt","gCountsBtoe_highpt",NoB,EoB);
    for(int i=0;i<NoB;i++)
    {
        gCountsBtoe_highpt->SetBinContent(i+1,myf_pT->Integral(EoB[i],EoB[i+1],para));
        // gCountsBtoe_highpt->SetBinError(i+1,myf_pT->IntegralError(EoB[i],EoB[i+1],para,NULL,1.E-10));
        gCountsBtoe_highpt->SetBinError(i+1,myf_pT->IntegralError(EoB[i],EoB[i+1],para,test->GetCovarianceMatrix()->GetMatrixArray(),1.E-2));
    }

    gCountsBtoe_highpt->Print("all");

      
    double integral[NoB];
    double integralErr[NoB];
    for(int i=0;i<NoB;i++)
    {
      int binx1 = hint_pt->GetXaxis()->FindBin(EoB[i]+1e-6);
      int binx2 = hint_pt->GetXaxis()->FindBin(EoB[i+1]-1e-6);
      integral[i] = hint_pt->IntegralAndError(binx1, binx2, integralErr[i], "width");
      cout << " " <<  EoB[i] << " - "<< EoB[i+1] << "  " << integral[i] << " , " << integralErr[i] << endl;
    }
    
    TCanvas*c2= new TCanvas("c2","c2",1600,800);
    gPad->SetLogy();
    gCountsBtoe_highpt->SetLineStyle(4);
    gCountsBtoe_highpt->SetMarkerStyle(4);
    gCountsBtoe_highpt->SetLineColor(4);
    gCountsBtoe_highpt->SetMarkerColor(4);
    gCountsBtoe_highpt->SetLineWidth(1);
    gCountsBtoe_highpt->SetMarkerSize(1);
    gCountsBtoe_highpt->Draw("E1");

    
    TLegend* leg2 = new TLegend(0.4,0.7,0.92,0.92);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.035);
    leg2->AddEntry(gCountsBtoe_highpt,"e From B yields","l");
    leg2->Draw();
    
    double x2[15]={1.25658,1.75109,2.25829,2.75534,3.25239,
                  3.75452,4.25537,4.74989,5.50307,6.50605,
                  7.50016,8.50314,10.5066,13.5066,17.5033};
    double y2[15]={0.00010547,8.79463e-5,7.01873e-5,5.37789e-5,3.99356e-5,
                  2.85169e-5,2.01726e-5,1.39823e-5,7.95573e-6,3.65818e-6,
                  1.66115e-6,7.63824e-7,1.9735e-7,2.4962e-8,2.75941e-9};
    TGraph* PHENIXgraph = new TGraph(15,x2,y2);
    PHENIXgraph->SetName("PHENIXgraph");
    
    TFile*out=new TFile("out/DataInput.root","recreate");
    out->cd();
    gCountsBtoe_highpt->Write();
    PHENIXgraph->Write();
    out->Close();
}
