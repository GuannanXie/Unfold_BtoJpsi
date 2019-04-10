
void calcChi2(TH1F *hData, TH1F *hMc, double &chi2, int &ndf){

   chi2 = 0;
   ndf = 0;
   for(int i=0;i<hData->GetNbinsX();i++){
      double xLo = hData->GetBinCenter(i+1)-hData->GetBinWidth(i+1)*0.5;
      double xHi = hData->GetBinCenter(i+1)+hData->GetBinWidth(i+1)*0.5;
      double y = hData->GetBinContent(i+1);
      double ey = hData->GetBinError(i+1);

      //int b1 = hMc->FindBin(xLo);
      //int b2 = hMc->FindBin(xHi)-1;
      //double yMc = hMc->Integral(b1,b2);
      
      double yMc = hMc->GetBinContent(i+1);
      double eyMc = hMc->GetBinError(i+1);

      if(ey>0&&fabs(yMc)>1e-2&&fabs(y)>1e-2){ 
         cout<<"yMc = "<<yMc<<" y = "<<y<<" ey = "<<ey<<" chi2 = "<<pow((y-yMc)/ey,2.)<<endl;
         chi2 += pow((y-yMc)/eyMc,2.);
         //chi2 += pow((y-yMc)/ey,2.);
         ndf++;
      }
   }
}

void eff_err(double m, double n, double* eff, double* err)
{
   *eff = 0;
   *err = 0;

   if (m<0 || n<=0) return;
   if (m>n) return;

   *eff = m/n;
   *err = sqrt((m+1)/(n+2)*((m+2)/(n+3)-(m+1)/(n+2)));

   return;
}

void eff_err(TH1* hNum, TH1* hDen, TH1* hEff)
{
   if(hNum->GetNbinsX()!=hDen->GetNbinsX()){ 
      cout<<"bins do not match!"<<endl;
      return;
   }
   for(int i=0;i<hNum->GetNbinsX();i++){
      double m = hNum->GetBinContent(i+1);
      double n = hDen->GetBinContent(i+1);
      if (m<0 || n<=0) continue;
      if (m>n) continue;

      double eff = m/n;
      double err = sqrt((m+1)/(n+2)*((m+2)/(n+3)-(m+1)/(n+2)));
      hEff->SetBinContent(i+1,eff);
      hEff->SetBinError(i+1,err);
   }
   return;
}

void eff_err2D(TH2 *hNum, TH2 *hDen, int nBins, float *range, TH1 **hNum1D,  TH1 **hDen1D, TH1 **hEff, int axis = 1, const char *tag = ""){

   for(int i=0;i<nBins;i++){
      int b1 = 0, b2 = 0;
      if(axis==1){
         b1 = hDen->GetYaxis()->FindBin(range[i]);
         b2 = hDen->GetYaxis()->FindBin(range[i+1])-1;
         hDen1D[i] = (TH1*)hDen->ProjectionX(Form("hDen1D_%s_%d",tag,i),b1,b2);
         hNum1D[i] = (TH1*)hNum->ProjectionX(Form("hNum1D_%s_%d",tag,i),b1,b2);
      }
      if(axis==2){
         b1 = hDen->GetXaxis()->FindBin(range[i]);
         b2 = hDen->GetXaxis()->FindBin(range[i+1])-1;
         hDen1D[i] = (TH1*)hDen->ProjectionY(Form("hDen1D_%s_%d",tag,i),b1,b2);
         hNum1D[i] = (TH1*)hNum->ProjectionY(Form("hNum1D_%s_%d",tag,i),b1,b2);
      }

      hEff[i] = (TH1*)hNum1D[i]->Clone(Form("hEff_%s_%d",tag,i));
      hEff[i]->SetName(Form("hEff_%s_%d",tag,i));
      hEff[i]->Reset();
      eff_err(hNum1D[i],hDen1D[i],hEff[i]);
   }
}

//3D efficiency, one axis integrated. 
void eff_err3D_integral(TH3 *hNum, TH3 *hDen, int nBins, float *range,  int intaxis = 3, int bin1, int bin2, TH1 **hNum1D, TH1 **hDen1D, TH1 **hEff, int axis = 1, const char *tag = ""){

   for(int i=0;i<nBins;i++){
      //cout<<"eff_err3D_integral: i = "<<i<<endl;
      int b1 = 0, b2 = 0;
      int b3 = 0, b4 = 0;
      if(axis==1){
         if(intaxis==2){
            b1 = bin1;
            b2 = bin2;
            b3 = hDen->GetZaxis()->FindBin(range[i]);
            b4 = hDen->GetZaxis()->FindBin(range[i+1])-1;
         }else if(intaxis==3){
            b1 = hDen->GetYaxis()->FindBin(range[i]);
            b2 = hDen->GetYaxis()->FindBin(range[i+1])-1;
            b3 = bin1;
            b4 = bin2;
         }else{
            cout<<"integrated axis can't be same as projection axis!"<<endl;
            return;
         }
         hDen1D[i] = (TH1*)hDen->ProjectionX(Form("hDen1D_%s_%d",tag,i),b1,b2,b3,b4);
         hNum1D[i] = (TH1*)hNum->ProjectionX(Form("hNum1D_%s_%d",tag,i),b1,b2,b3,b4);
      }
      if(axis==2){
         if(intaxis==1){
            b1 = bin1;
            b2 = bin2;
            b3 = hDen->GetZaxis()->FindBin(range[i]);
            b4 = hDen->GetZaxis()->FindBin(range[i+1])-1;
         }else if(intaxis==3){
            b1 = hDen->GetXaxis()->FindBin(range[i]);
            b2 = hDen->GetXaxis()->FindBin(range[i+1])-1;
            b3 = bin1;
            b4 = bin2;
         }else{
            cout<<"integrated axis can't be same as projection axis!"<<endl;
            return;
         }
         hDen1D[i] = (TH1*)hDen->ProjectionY(Form("hDen1D_%s_%d",tag,i),b1,b2,b3,b4);
         hNum1D[i] = (TH1*)hNum->ProjectionY(Form("hNum1D_%s_%d",tag,i),b1,b2,b3,b4);
      }
      if(axis==3){
         if(intaxis==1){
            b1 = bin1;
            b2 = bin2;
            b3 = hDen->GetYaxis()->FindBin(range[i]);
            b4 = hDen->GetYaxis()->FindBin(range[i+1])-1;
         }else if(intaxis==2){
            b1 = hDen->GetXaxis()->FindBin(range[i]);
            b2 = hDen->GetXaxis()->FindBin(range[i+1])-1;
            b3 = bin1;
            b4 = bin2;
         }else{
            cout<<"integrated axis can't be same as projection axis!"<<endl;
            return;
         }
         hDen1D[i] = (TH1*)hDen->ProjectionZ(Form("hDen1D_%s_%d",tag,i),b1,b2,b3,b4);
         hNum1D[i] = (TH1*)hNum->ProjectionZ(Form("hNum1D_%s_%d",tag,i),b1,b2,b3,b4);
      }
      if(axis<1||axis>3){ cout<<"wrong axis setting!"<<endl; return;}
      hEff[i] = (TH1*)hNum1D[i]->Clone(Form("hEff_%s_%d",tag,i));
      hEff[i]->SetName(Form("hEff_%s_%d",tag,i));
      hEff[i]->Reset();
      eff_err(hNum1D[i],hDen1D[i],hEff[i]);
   }
}

//3D efficiency, project to one bin. 
void eff_err3D(const TH3 *hNum, const TH3 *hDen, const float low1, const float high1, const float low2, const float high2, TH1 &*hNum1D, TH1 *&hDen1D, TH1 *&hEff, int axis = 1, const char *tag = ""){

      int b1 = 0, b2 = 0;
      int b3 = 0, b4 = 0;
      if(axis==1){
         b1 = hDen->GetYaxis()->FindBin(low1);
         b2 = hDen->GetYaxis()->FindBin(high1)-1;
         b3 = hDen->GetZaxis()->FindBin(low2);
         b4 = hDen->GetZaxis()->FindBin(high2)-1;
         hDen1D = (TH1*)hDen->ProjectionX(Form("hDen1D_%s",tag),b1,b2,b3,b4);
         hNum1D = (TH1*)hNum->ProjectionX(Form("hNum1D_%s",tag),b1,b2,b3,b4);
      }
      if(axis==2){
         b1 = hDen->GetXaxis()->FindBin(low1);
         b2 = hDen->GetXaxis()->FindBin(high1)-1;
         b3 = hDen->GetZaxis()->FindBin(low2);
         b4 = hDen->GetZaxis()->FindBin(high2)-1;
         hDen1D = (TH1*)hDen->ProjectionY(Form("hDen1D_%s",tag),b1,b2,b3,b4);
         hNum1D = (TH1*)hNum->ProjectionY(Form("hNum1D_%s",tag),b1,b2,b3,b4);
      }
      if(axis==3){
         b1 = hDen->GetXaxis()->FindBin(low1);
         b2 = hDen->GetXaxis()->FindBin(high1)-1;
         b3 = hDen->GetYaxis()->FindBin(low2);
         b4 = hDen->GetYaxis()->FindBin(high2)-1;
         hDen1D = (TH1*)hDen->ProjectionZ(Form("hDen1D_%s_%d",tag),b1,b2,b3,b4);
         hNum1D = (TH1*)hNum->ProjectionZ(Form("hNum1D_%s_%d",tag),b1,b2,b3,b4);
      }
      if(axis<1||axis>3){ cout<<"wrong axis setting!"<<endl; return;}
      cout<<"b1 = "<<b1<<" b2 = "<<b2<<" b3 = "<<b3<<" b4 = "<<b4<<endl;
      hEff = (TH1*)hNum1D->Clone(Form("hEff_%s",tag));
      hEff->SetName(Form("hEff_%s",tag));
      hEff->Reset();
      eff_err(hNum1D,hDen1D,hEff);
}

void rebinHisto(TH1 *h, TH1 *hreb, bool corrBW = false){

   int bin = h->GetNbinsX();
   int newbin = hreb->GetNbinsX();
   if(newbin>bin){
      cout<<" failed to rebin histogram, bins are more than original histogram!"<<endl;
      return;
   }

   //double min = h->GetXaxis()->GetXmin();
   //double max = h->GetXaxis()->GetXmax();
   //double dm  = (max-min)/bin;

   double yi = 0;
   double syi = 0;
   double dw = h->GetBinWidth(1);
   for(int i=0;i<newbin;i++){
      double mrcen = hreb->GetBinCenter(i+1);
      double mrwid = hreb->GetBinWidth(i+1);
      double mrlow = mrcen-mrwid/2.+dw/2.;
      double mrup  = mrcen+mrwid/2.+dw/2.;

      int binlow = h->FindBin(mrlow);
      int binup  = h->FindBin(mrup)-1;
      int nb    = binup-binlow+1;
      //cout<<"mrcen = "<<mrcen<<" mrlow:"<<mrlow<<" mrup:"<<mrup<<endl;
      //cout<<"binlow:"<<binlow<<" binup:"<<binup<<" nb:"<<nb<<endl;

      double ey1 = 0.;
      double se1 = 0.;
      double y1  = h->Integral(binlow,binup);
      for(int j=binlow;j<=binup;j++){
         se1+= pow(h->GetBinError(j),2.);
      }
      ey1 = sqrt(se1);

      if(corrBW){
         hreb->SetBinContent(i+1,y1/nb);
         hreb->SetBinError(i+1,ey1/nb);
      }else{
         hreb->SetBinContent(i+1,y1);
         hreb->SetBinError(i+1,ey1);
      }
   }	
   return;
}

void rebinHisto(TH2 *h, TH2 *hreb){

   int nBinsX = h->GetNbinsX();
   int nBinsY = h->GetNbinsY();

   int nNewBinsX = hreb->GetNbinsX();
   int nNewBinsY = hreb->GetNbinsY();
   if(nNewBinsX>nBinsX||nNewBinsY>nBinsY){
      cout<<" failed to rebin histogram, bins are more than original histogram!"<<endl;
      cout<<" old binx:"<<nBinsX<<" biny:"<<nBinsY<<endl;
      cout<<" new binx:"<<nNewBinsX<<" NewBiny:"<<nNewBinsY<<endl;
      return;
   }

   for(int i=0;i<nNewBinsX;i++){
      for(int j=0;j<nNewBinsY;j++){
         double xi = hreb->GetXaxis()->GetBinCenter(i+1);
         double wdxi = hreb->GetXaxis()->GetBinWidth(i+1);
         double yi = hreb->GetYaxis()->GetBinCenter(j+1);
         double wdyi = hreb->GetYaxis()->GetBinWidth(j+1);

         double xlowi = xi-wdxi/2.+1e-10;
         double xupi = xi+wdxi/2.-1e-10;

         double ylowi = yi-wdyi/2.+1e-10;
         double yupi = yi+wdyi/2.-1e-10;

         int binxlow = h->GetXaxis()->FindBin(xlowi);
         int binxup  = h->GetXaxis()->FindBin(xupi);
         double xi0 = h->GetXaxis()->GetBinCenter(binxlow);
         double wdxi0 = h->GetXaxis()->GetBinWidth(binxlow);
         if(xi0-wdxi0/2.>xlowi) binxlow-=1;
         else if(xi0+wdxi0/2.<=xlowi) binxlow+=1;
         else binxlow = binxlow; 
         double xi1 = h->GetXaxis()->GetBinCenter(binxup);
         double wdxi1 = h->GetXaxis()->GetBinWidth(binxup);
         if(xi1-wdxi1/2.>xupi) binxup-=1;
         else if(xi1+wdxi1/2.<=xupi) binxup+=1;
         else binxup = binxup; 

         int binylow = h->GetYaxis()->FindBin(ylowi);
         int binyup  = h->GetYaxis()->FindBin(yupi);
         double yi0 = h->GetYaxis()->GetBinCenter(binylow);
         double wdyi0 = h->GetYaxis()->GetBinWidth(binylow);
         if(yi0-wdyi0/2.>ylowi) binylow-=1;
         else if(yi0+wdyi0/2.<=ylowi) binylow+=1;
         else binylow = binylow; 
         double yi1 = h->GetYaxis()->GetBinCenter(binyup);
         double wdyi1 = h->GetYaxis()->GetBinWidth(binyup);
         if(yi1-wdyi1/2.>yupi) binyup-=1;
         else if(yi1+wdyi1/2.<=yupi) binyup+=1;
         else binyup = binyup; 
         if(binxlow>binxup||binylow>binyup){
            cout<<" Error: Invalid integral bin xlow="<<binxlow<<" xup="<<binxup<<" ylow="<<binylow<<" yup="<<binyup<<endl;
            return;
         }
         double n = h->Integral(binxlow,binxup,binylow,binyup);
         double ey = 0., sey = 0;
         for(int ib=binxlow;ib<=binxup;ib++){
            for(int jb=binylow;jb<=binyup;jb++){
               sey+=pow(h->GetBinError(ib,jb),2.);
            }
         }
         double ey = sqrt(sey);
         hreb->SetBinContent(i+1,j+1,n);
         hreb->SetBinError(i+1,j+1,ey);
      }
   }

   return;
}

void drawOverSizeErr(TGraphErrors *gr, float scale=0.85){

   double errX   = gStyle->GetErrorX();
   double endErr = gStyle->GetEndErrorSize();
   int first = gr->GetXaxis()->GetFirst();
   int last  = gr->GetXaxis()->GetLast();
   //int n = gr->GetN();
   //for(int i=0;i<n;i++)
   for(int i=first-1;i<last;i++){

      double x,y,ex,ey;
      gr->GetPoint(i,x,y);
      ex=gr->GetErrorX(i);
      ey=gr->GetErrorY(i);

      Style_t sty = gr->GetMarkerStyle();
      Size_t  siz = gr->GetMarkerSize();
      Color_t col = gr->GetMarkerColor();
      Color_t lCol = gr->GetLineColor();
      Float_t wid = gr->GetLineWidth()/100.;
      TMarker *tm = new TMarker(x,y,sty);
      tm->SetMarkerColor(col);
      tm->SetMarkerSize(siz);
      tm->SetMarkerStyle(sty);

      TArrow *arX,*arY;
      if(ey>=y){
         if(endErr>0){
            arY = new TArrow(x,y+ey,x,y-scale*y,wid,"|---------|>");
         }else{
            arY = new TArrow(x,y+ey,x,y-scale*y,wid,"---------|>");
         }
         arY->SetLineColor(lCol);
         arY->SetFillColor(lCol);
         arY->DrawClone();
      }else{
         if(endErr>0){
            arY = new TArrow(x,y+ey,x,y-ey,wid,"|---------|");
         }else{
            arY = new TArrow(x,y+ey,x,y-ey,wid,"---------");
         }
         arY->SetLineColor(lCol);
         arY->SetFillColor(lCol);
         arY->DrawClone();
      }
      if(errX>0){
         if(endErr>0){
            arX = new TArrow(x-ex,y,x+ex,y,wid,"|---------|");
         }else{
            arX = new TArrow(x-ex,y,x+ex,y,wid,"---------");
         }
         arX->SetLineColor(lCol);
         arX->SetFillColor(lCol);
         arX->DrawClone();
      }

      tm->DrawClone();
      arY->Delete();
      tm->Delete();
      if(arX) arX->Delete();
   }
}

void drawOverSizeErr(TGraphErrors *gr, float scale=0.85, double xmin, double xmax){

   double errX   = gStyle->GetErrorX();
   double endErr = gStyle->GetEndErrorSize();
   int first = gr->GetXaxis()->GetFirst();
   int last  = gr->GetXaxis()->GetLast();
   //int n = gr->GetN();
   //for(int i=0;i<n;i++)
   for(int i=first-1;i<last;i++){

      double x,y,ex,ey;
      gr->GetPoint(i,x,y);
      ex=gr->GetErrorX(i);
      ey=gr->GetErrorY(i);
      if(x<xmin||x>xmax) continue;

      Style_t sty = gr->GetMarkerStyle();
      Size_t  siz = gr->GetMarkerSize();
      Color_t col = gr->GetMarkerColor();
      Color_t lCol = gr->GetLineColor();
      Float_t wid = gr->GetLineWidth()/100.;
      TMarker *tm = new TMarker(x,y,sty);
      tm->SetMarkerColor(col);
      tm->SetMarkerSize(siz);
      tm->SetMarkerStyle(sty);

      TArrow *arX,*arY;
      if(ey>=y){
         if(endErr>0){
            arY = new TArrow(x,y+ey,x,y-scale*y,wid,"|---------|>");
         }else{
            arY = new TArrow(x,y+ey,x,y-scale*y,wid,"---------|>");
         }
         arY->SetLineColor(lCol);
         arY->SetFillColor(lCol);
         arY->DrawClone();
      }else{
         if(endErr>0){
            arY = new TArrow(x,y+ey,x,y-ey,wid,"|---------|");
         }else{
            arY = new TArrow(x,y+ey,x,y-ey,wid,"---------");
         }
         arY->SetLineColor(lCol);
         arY->SetFillColor(lCol);
         arY->DrawClone();
      }
      if(errX>0){
         double xlow = x-ex;
         double xhigh = x+ex;
         if(xlow<xmin) xlow = xmin;
         if(xhigh>xmax) xhigh = xmax;
         if(endErr>0){
            arX = new TArrow(xlow,y,xhigh,y,wid,"|---------|");
         }else{
            arX = new TArrow(xlow,y,xhigh,y,wid,"---------");
         }
         arX->SetLineColor(lCol);
         arX->SetFillColor(lCol);
         arX->DrawClone();
      }

      tm->DrawClone();
      arY->Delete();
      tm->Delete();
      if(arX) arX->Delete();
   }
}

void drawOverSizeErr(TH1 *h, float scale=0.85, TCanvas *c1, int iPad){
   double endErr = gStyle->GetEndErrorSize();
   double errX   = gStyle->GetErrorX();
   int first = h->GetXaxis()->GetFirst();
   int last  = h->GetXaxis()->GetLast();
   for(int i=first-1;i<last;i++){

      double x,y,ex,ey;
      x =h->GetBinCenter(i+1);
      y =h->GetBinContent(i+1);
      ex=h->GetBinWidth(i+1)/2.;
      ey=h->GetBinError(i+1);

      Style_t sty = h->GetMarkerStyle();
      Size_t  siz = h->GetMarkerSize();
      Color_t col = h->GetMarkerColor();
      Color_t lCol = h->GetLineColor();
      Float_t wid = h->GetLineWidth()/100.;
      TMarker *tm = new TMarker(x,y,sty);
      tm->SetMarkerColor(col);
      tm->SetMarkerSize(siz);
      tm->SetMarkerStyle(sty);

      c1->cd(iPad);
      TArrow *arX,*arY;
      if(ey>=0.99*y){
         if(endErr>0){
            arY = new TArrow(x,y+ey,x,y-scale*y,wid,"|---------|>");
         }else{
            arY = new TArrow(x,y+ey,x,y-scale*y,wid,"---------|>");
         }
         arY->SetLineColor(lCol);
         arY->SetFillColor(lCol);
         arY->DrawClone();
      }else{
         if(endErr>0){
            arY = new TArrow(x,y+ey,x,y-ey,wid,"|---------|");
         }else{
            arY = new TArrow(x,y+ey,x,y-ey,wid,"---------");
         }
         arY->SetLineColor(lCol);
         arY->SetFillColor(lCol);
         arY->DrawClone();
      }
      if(errX>0){
         if(endErr>0){
            arX = new TArrow(x-ex,y,x+ex,y,wid,"|---------|");
         }else{
            arX = new TArrow(x-ex,y,x+ex,y,wid,"---------");
         }
         arX->SetLineColor(lCol);
         arX->SetFillColor(lCol);
         arX->DrawClone();
      }
      tm->DrawClone();
      arY->Delete();
      if(arX) arX->Delete();
      tm->Delete();
   }
}

void drawOverSizeErr(TH1 *h, float scale=0.85, float xmin, float xmax){
   double endErr = gStyle->GetEndErrorSize();
   double errX   = gStyle->GetErrorX();
   int first = h->GetXaxis()->GetFirst();
   int last  = h->GetXaxis()->GetLast();
   for(int i=first-1;i<last;i++){

      double x,y,ex,ey;
      x =h->GetBinCenter(i+1);
      y =h->GetBinContent(i+1);
      ex=h->GetBinWidth(i+1)/2.;
      ey=h->GetBinError(i+1);
      if(x<xmin||x>xmax) continue;

      Style_t sty = h->GetMarkerStyle();
      Size_t  siz = h->GetMarkerSize();
      Color_t col = h->GetMarkerColor();
      Color_t lCol = h->GetLineColor();
      Float_t wid = h->GetLineWidth()/100.;
      TMarker *tm = new TMarker(x,y,sty);
      tm->SetMarkerColor(col);
      tm->SetMarkerSize(siz);
      tm->SetMarkerStyle(sty);

      TArrow *arX,*arY;
      if(ey>=0.99*y){
         if(endErr>0){
            arY = new TArrow(x,y+ey,x,y-scale*y,wid,"|---------|>");
         }else{
            arY = new TArrow(x,y+ey,x,y-scale*y,wid,"---------|>");
         }
         arY->SetLineColor(lCol);
         arY->SetFillColor(lCol);
         arY->DrawClone();
      }else{
         if(endErr>0){
            arY = new TArrow(x,y+ey,x,y-ey,wid,"|---------|");
         }else{
            arY = new TArrow(x,y+ey,x,y-ey,wid,"---------");
         }
         arY->SetLineColor(lCol);
         arY->SetFillColor(lCol);
         arY->DrawClone();
      }
      if(errX>0){
         if(endErr>0){
            arX = new TArrow(x-ex,y,x+ex,y,wid,"|---------|");
         }else{
            arX = new TArrow(x-ex,y,x+ex,y,wid,"---------");
         }
         arX->SetLineColor(lCol);
         arX->SetFillColor(lCol);
         arX->DrawClone();
      }
      tm->DrawClone();
      arY->Delete();
      if(arX) arX->Delete();
      tm->Delete();
   }
}

void scaleGraph(TGraphAsymmErrors *gr, double factor){

   int np = gr->GetN();
   for(int i=0;i<np;i++){
      double x,y,exl,exh,eyl,eyh;
      gr->GetPoint(i,x,y);
      exl = gr->GetErrorXlow(i);
      exh = gr->GetErrorXhigh(i);
      eyl = gr->GetErrorYlow(i);
      eyh = gr->GetErrorYhigh(i);

      y = y*factor;
      eyl = eyl*factor;
      eyh = eyh*factor;
      gr->SetPoint(i,x,y);
      gr->SetPointEYlow(i,eyl);
      gr->SetPointEYhigh(i,eyh);
   }
   return;
}



void scaleGraph(TGraphErrors *gr, double factor){

   int np = gr->GetN();
   for(int i=0;i<np;i++){
      double x,y,ex,ey;
      gr->GetPoint(i,x,y);
      ex = gr->GetErrorX(i);
      ey = gr->GetErrorY(i);

      y = y*factor;
      ey = ey*factor;
      gr->SetPoint(i,x,y);
      gr->SetPointError(i,ex,ey);
   }
   return;
}

void scaleGraph(TGraph *gr, double factor){

   int np = gr->GetN();
   for(int i=0;i<np;i++){
      double x,y;
      gr->GetPoint(i,x,y);
      y = y*factor;
      gr->SetPoint(i,x,y);
   }
   return;
}


double evalHist(TH1 *h, double x){

   if(!h) return 0;
   int bin = h->FindBin(x);
   int binlw = -1;
   int binup = -1;
   double xlw=0.,xup=0.;
   if(h->GetBinCenter(bin)<x){
      binlw = bin-1;
      binup = bin;
   }else{
      binlw = bin;
      binup = bin+1;
   }
   if(binlw<0||binup<0||binup>h->GetNbinsX()) return 0;

   xlw = h->GetBinCenter(binlw);
   xup = h->GetBinCenter(binup);

   double a1 = h->GetBinContent(binlw);
   double a2 = h->GetBinContent(binup);
   double p0 = (a1-a2)/(xlw-xup);
   double p1 = (a2*xlw-a1*xup)/(xlw-xup);

   return p0*x+p1;
}

void formatHist(TH1 *h, int mark=20, int markCol=1, float size=1, int col=1, int sty=1, float wid=2){
   h->SetMarkerStyle(mark);	
   h->SetMarkerColor(markCol);	
   h->SetMarkerSize(size);	
   h->SetLineColor(col);	
   h->SetLineStyle(sty);	
   h->SetLineWidth(wid);	
}

void formatHist(TProfile *h, int mark=20, int markCol=1, float size=1, int col=1, int sty=1, float wid=2){
   h->SetMarkerStyle(mark);	
   h->SetMarkerColor(markCol);	
   h->SetMarkerSize(size);	
   h->SetLineColor(col);	
   h->SetLineStyle(sty);	
   h->SetLineWidth(wid);	
}

void formatGraph(TGraph *g, int mark=20, int markCol=1, float size=1, int col=1, int sty=1, float wid=2){
   g->SetMarkerStyle(mark);	
   g->SetMarkerColor(markCol);	
   g->SetMarkerSize(size);	
   g->SetLineColor(col);	
   g->SetLineStyle(sty);	
   g->SetLineWidth(wid);	
}

void formatGraph(TGraphErrors *g, int mark=20, int markCol=1, float size=1, int col=1, int sty=1, float wid=2){
   g->SetMarkerStyle(mark);	
   g->SetMarkerColor(markCol);	
   g->SetMarkerSize(size);	
   g->SetLineColor(col);	
   g->SetLineStyle(sty);	
   g->SetLineWidth(wid);	
}


void formatLegend(TLegend *l, int border=0, int col=0, int txtsize=0.04, int txtcol=1){
   l->SetBorderSize(border);
   l->SetFillColor(col);
   l->SetTextSize(txtsize);
   l->SetTextColor(txtcol);
}

void formatPad(TPad *pad,float left, float right, float top, float bottom){
   pad->SetFillColor(10);
   pad->SetBorderMode(0);
   pad->SetBorderSize(0);
   pad->SetFrameFillColor(10);
   pad->SetFrameBorderMode(0);
   pad->SetFrameBorderSize(0);
   pad->SetLeftMargin(left);
   pad->SetRightMargin(right);
   pad->SetTopMargin(top);
   pad->SetBottomMargin(bottom);
}

void exportTGraphErrors(TGraphErrors *g, const char *outname){

   ofstream outf(outname);
   int nPoints = g->GetN();
   if(nPoints<=0) return;
   double *x = new double [nPoints];
   double *y = new double [nPoints];
   double *ex = new double [nPoints];
   double *ey = new double [nPoints];
   for(int i=0;i<nPoints;i++){
      g->GetPoint(i,x[i],y[i]); 
      ex[i] = g->GetErrorX(i);
      ey[i] = g->GetErrorY(i);
   }
   for(int i=0;i<nPoints;i++){
      outf<<x[i]<<" \t "<<y[i]<<" \t "<<ex[i]<<" \t "<<ey[i]<<endl;
   }
   out.close();
   delete [] x;
   delete [] y;
   delete [] ex;
   delete [] ey;
}

void importTGraphErrors(TGraphErrors *g, const char *fname){

   ifstream inf(fname);
   int nPoints = 0;
   if(inf.fail()){
      cout<<"Error: No input file!"<<endl;
      return;
   }
   string tmp;
   while(getline(inf,tmp,'\n')){
      nPoints++;
   }
   double *x = new double [nPoints];
   double *y = new double [nPoints];
   double *ex = new double [nPoints];
   double *ey = new double [nPoints];
   for(int i=0;i<nPoints;i++){
      inf>>x[i]>>" \t ">>y[i]>>" \t ">>ex[i]>>" \t ">>ey[i]>>endl;
   }
   g = new TGraphErrors(nPoints,x,y,ex,ey);
   inf.close();
   delete [] x;
   delete [] y;
   delete [] ex;
   delete [] ey;
}

double integralGraph(TGraph *g, double xmin, double xmax){
   double sum = 0.;
   for(int i=0;i<g->GetN();i++){
      double x,y;
      g->GetPoint(i,x,y);
      if(x>xmin && x<xmax) sum+=y;
   }
   return sum;
}

double integralGraph(TGraphErrors *g, double xmin, double xmax){
   double sum = 0.;
   for(int i=0;i<g->GetN();i++){
      double x,y;
      g->GetPoint(i,x,y);
      if(x>xmin && x<xmax) sum+=y;
   }
   return sum;
}

TGraphErrors *rebinGraph(TGraphErrors *g, int nRebins){

   int N = g->GetN();
   if(N%nRebins !=0 ){ cout<<" warning: rebin "<<N<<" by ngroup="<<nRebins<<". Last few bins will be dropped!"<<endl;}
   double *x = new double [N/nRebins];
   double *y = new double [N/nRebins];
   int nNewBins = 0;
   double *xtmp = new double [nRebins];
   double *ytmp = new double [nRebins];
   int ireb = 0;
   for(int i=0;i<N/nRebins;i++){
      x[i] = 0;
      y[i] = 0;
   }
   for(int i=0;i<nRebins;i++){
      xtmp[i] = 0;
      ytmp[i] = 0;
   }
   for(int i=0;i<N-N%nRebins;i++){
      double xi=0,yi=0;
      g->GetPoint(i,xi,yi);
      xtmp[ireb] = xi;
      y[nNewBins] += yi;
      ireb++;
      if(ireb%nRebins==0){
         x[nNewBins] = (xtmp[ireb-1]+xtmp[0])/2.;
         nNewBins++;
         ireb=0;
      }
   }
   TGraphErrors *gnew = new TGraphErrors(nNewBins,x,y,0,0);

   //delete [] x;
   //delete [] y;
   //delete [] xtmp;
   //delete [] ytmp;
   return gnew;
}

TGraphErrors *rebinAvgGraph(TGraphErrors *g, int nRebins){

   int N = g->GetN();
   if(N%nRebins !=0 ){ cout<<" warning: rebin "<<N<<" by ngroup="<<nRebins<<". Last few bins will be dropped!"<<endl;}
   double *x = new double [N/nRebins];
   double *y = new double [N/nRebins];
   int nNewBins = 0;
   double *xtmp = new double [nRebins];
   double *ytmp = new double [nRebins];
   int ireb = 0;
   for(int i=0;i<N/nRebins;i++){
      x[i] = 0;
      y[i] = 0;
   }
   for(int i=0;i<nRebins;i++){
      xtmp[i] = 0;
      ytmp[i] = 0;
   }
   for(int i=0;i<N-N%nRebins;i++){
      double xi=0,yi=0;
      g->GetPoint(i,xi,yi);
      xtmp[ireb] = xi;
      y[nNewBins] += yi;
      ireb++;
      if(ireb%nRebins==0){
         x[nNewBins] = (xtmp[ireb-1]+xtmp[0])/2.;
         y[nNewBins] = y[nNewBins]/nRebins;
         nNewBins++;
         ireb=0;
      }
   }
   TGraphErrors *gnew = new TGraphErrors(nNewBins,x,y,0,0);

   //delete [] x;
   //delete [] y;
   //delete [] xtmp;
   //delete [] ytmp;
   return gnew;
}

void printGraph(TGraph *g){

   for(int i=0;i<g->GetN();i++){
      double xi,yi;
      g->GetPoint(i,xi,yi);
      cout<<"point "<<i<<" xi="<<xi<<" yi="<<yi<<endl;
   }
}

void printGraphErrors(TGraphErrors *g){

   for(int i=0;i<g->GetN();i++){
      double xi,yi;
      g->GetPoint(i,xi,yi);
      cout<<"point "<<i<<" xi="<<xi<<" yi="<<yi<<" xerr = "<<g->GetErrorX(i)<<" yerr = "<<g->GetErrorY(i)<<endl;
   }
}


Double_t squareWave(Double_t *x, Double_t *par)
{
   Double_t y = (TMath::Abs(x[0]-par[0])-par[2])/par[3];
   return par[1]/(1.+TMath::Exp(y))+par[4];
}

TH1F *graphToHistogram(TGraphErrors *g){
   int n = g->GetN();
   double *x = new double [n];
   double *y = new double [n];
   double *ex = new double [n];
   double *ey = new double [n];
   double *edge = new double [n+1];
   for(int i=0;i<n;i++){
      g->GetPoint(i,x[i],y[i]);
      ex[i] = g->GetErrorX(i);
      ey[i] = g->GetErrorY(i);
      edge[i] = x[i]-ex[i];
      cout<<"x = "<<x[i]<<" dx = "<<ex[i]<<" edge = "<<edge[i]<<endl;
   }
   edge[n] = x[n-1]+ex[n-1];

   TH1F *h = new TH1F(Form("h%s",g->GetName()),Form("h%s",g->GetTitle()),n,edge);
   for(int i=0;i<n-1;i++){
      h->SetBinContent(i+1,y[i]);
      h->SetBinError(i+1,ey[i]);
   }
   delete [] x;
   delete [] y;
   delete [] ex;
   delete [] ey;
   delete [] edge;
   return h;
}

void histToGraph(TH1F *h1, TGraphErrors &*g1){
   g1 = new TGraphErrors(h1->GetNbinsX());
   for(int i=0;i<h1->GetNbinsX();i++){
      double x = h1->GetBinCenter(i+1);
      double ex = h1->GetBinWidth(i+1)/2.;
      double y = h1->GetBinContent(i+1);
      double ey = h1->GetBinError(i+1);
      g1->SetPoint(i,x,y);
      g1->SetPointError(i,ex,ey);
   }
}

void histRatio(TH1F *h1, TH1F *h2, TH1F &*hr){
   hr = (TH1F*)h1->Clone(Form("hR_%s",h1->GetName()));
   if(!hr->GetDefaultSumw2()) hr->Sumw2();
   hr->Divide(h2);
}

void compare1D(TH1 *h1, TH1 *h2, char *tag1, char *tag2, TCanvas *c1, TPDF *pdf){
   h1->SetMarkerStyle(24);
   c1->cd();
   h1->SetAxisRange(1,h1->GetMaximum()*10,"Y");
   h1->Draw("p");
   h2->SetLineColor(2);
   h2->Draw("histsame");
   TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
   leg->SetFillColor(10);
   leg->SetBorderSize(0);
   leg->AddEntry(h1,tag1,"p");
   leg->AddEntry(h2,tag2,"l");
   leg->Draw();
   c1->Update();
}

void compare2D(TH2D *h1, TH2D *h2, TCanvas *c1, TPDF *pdf){
   c1->cd();
   TH2D *r = (TH2D*)h1->Clone();
   r->Divide(h2);
   gStyle->SetPaintTextFormat("1.2g");
   r->Draw("colz");
   r->Draw("textsame");
   TLatex *tt1 = new TLatex(0.22,0.87,Form("#frac{%s}{%s}",h1->GetName(),h2->GetName()));
   tt1->SetNDC(1);
   tt1->SetTextSize(0.04);
   tt1->Draw();
   c1->Update();
}

void extractSignal2D(TH2 *hUS2D, TH2 *hLS2D, int axis, int bin1, int bin2, int i, TH1 &*hus, TH1 &*hls, TH1 &*hs){

   if(axis==1){ //X
      hus = (TH1F*)hUS2D->ProjectionX(Form("%s_%d",hUS2D->GetName(),i),bin1,bin2);
      hls = (TH1F*)hLS2D->ProjectionX(Form("%s_%d",hLS2D->GetName(),i),bin1,bin2);
   }
   if(axis==2){
      hus = (TH1F*)hUS2D->ProjectionY(Form("%s_%d",hUS2D->GetName(),i),bin1,bin2);
      hls = (TH1F*)hLS2D->ProjectionY(Form("%s_%d",hLS2D->GetName(),i),bin1,bin2);
   }
   hs = (TH1F*)hus->Clone(Form("hSig_%s_%d",hUS2D->GetName(),i));
   if(!hs->GetDefaultSumw2()) hs->Sumw2();
   hs->Add(hls,-1.);
}

void extractSignal3D(TH3 *hUS3D, TH3 *hLS3D, int axis, int bin11, int bin12, int bin21, int bin22, int i, TH1F &*hus, TH1F &*hls, TH1F &*hs){

   if(axis==1){ //X
      hus = (TH1F*)hUS3D->ProjectionX(Form("%s_%d",hUS3D->GetName(),i),bin11,bin12,bin21,bin22);
      hls = (TH1F*)hLS3D->ProjectionX(Form("%s_%d",hLS3D->GetName(),i),bin11,bin12,bin21,bin22);
   }
   if(axis==2){ //Y
      //hUS3D->GetXaxis()->SetRange(bin11,bin12);
      //hUS3D->GetZaxis()->SetRange(bin21,bin22);
      //hLS3D->GetXaxis()->SetRange(bin11,bin12);
      //hLS3D->GetZaxis()->SetRange(bin21,bin22);
      hus = (TH1F*)hUS3D->ProjectionY(Form("%s_%d",hUS3D->GetName(),i),bin11,bin12,bin21,bin22);
      hls = (TH1F*)hLS3D->ProjectionY(Form("%s_%d",hLS3D->GetName(),i),bin11,bin12,bin21,bin22);
   }
   if(axis==3){ //Z
      hus = (TH1F*)hUS3D->ProjectionZ(Form("%s_%d",hUS3D->GetName(),i),bin11,bin12,bin21,bin22);
      hls = (TH1F*)hLS3D->ProjectionZ(Form("%s_%d",hLS3D->GetName(),i),bin11,bin12,bin21,bin22);
   }
   hs = (TH1F*)hus->Clone(Form("hSig_%s_%d",hUS3D->GetName(),i));
   if(!hs->GetDefaultSumw2()) hs->Sumw2();
   hs->Add(hls,-1.);
}

void extractSignalnD(THnSparseF &*hUS, THnSparseF &*hLS, int axis, int nOtherAxis, int *otherAxis, double *min, double *max, const char *tag, TH1F &*hus, TH1F &*hls, TH1F &*hs){
   for(int i=0;i<nOtherAxis;i++){
      //cout<<"otherAxis = "<<otherAxis[i]<<endl;
      hUS->GetAxis(otherAxis[i])->SetRange(min[otherAxis[i]],max[otherAxis[i]]);
      hLS->GetAxis(otherAxis[i])->SetRange(min[otherAxis[i]],max[otherAxis[i]]);
   }

   hus = (TH1F*)hUS->Projection(axis)->Clone(Form("hus_%s",tag));
   hls = (TH1F*)hLS->Projection(axis)->Clone(Form("hls_%s",tag));
   if(!hus->GetDefaultSumw2()) hus->Sumw2();
   if(!hls->GetDefaultSumw2()) hls->Sumw2();
   hs = (TH1F*)hus->Clone(Form("hsig_%s",tag));
   if(!hs->GetDefaultSumw2()) hs->Sumw2();
   hs->Add(hls,-1.);
}

void extractSignalnDTo2D(THnSparseF *hUS, THnSparseF *hLS, int axis1, int axis2, int nOtherAxis, int *otherAxis, double *min, double *max, const char *tag, TH2F &*hus, TH2F &*hls, TH2F &*hs){
   for(int i=0;i<nOtherAxis;i++){
      cout<<"otherAxis = "<<otherAxis[i]<<" min = "<<min[otherAxis[i]]<<" max = "<<max[otherAxis[i]]<<endl;
      hUS->GetAxis(otherAxis[i])->SetRange(min[otherAxis[i]],max[otherAxis[i]]);
      hLS->GetAxis(otherAxis[i])->SetRange(min[otherAxis[i]],max[otherAxis[i]]);
   }
      
   hus = (TH2F*)hUS->Projection(axis1,axis2)->Clone(Form("hus2D_%s",tag));
   hls = (TH2F*)hLS->Projection(axis1,axis2)->Clone(Form("hls2D_%s",tag));
   if(!hus->GetDefaultSumw2()) hus->Sumw2();
   if(!hls->GetDefaultSumw2()) hls->Sumw2();
   hs = (TH2F*)hus->Clone(Form("hsig2D_%s",tag));
   if(!hs->GetDefaultSumw2()) hs->Sumw2();
   hs->Add(hls,-1.);
}

void extractSignalnDTo3D(THnSparseF *hUS, THnSparseF *hLS, int axis1, int axis2, int axis3, int nOtherAxis, int *otherAxis, double *min, double *max, const char *tag, TH3F &*hus, TH3F &*hls, TH3F &*hs){
   for(int i=0;i<nOtherAxis;i++){
      cout<<"otherAxis = "<<otherAxis[i]<<" min = "<<min[otherAxis[i]]<<" max = "<<max[otherAxis[i]]<<endl;
      hUS->GetAxis(otherAxis[i])->SetRange(min[otherAxis[i]],max[otherAxis[i]]);
      hLS->GetAxis(otherAxis[i])->SetRange(min[otherAxis[i]],max[otherAxis[i]]);
   }
      
   hus = (TH3F*)hUS->Projection(axis1,axis2,axis3)->Clone(Form("hus3D_%s",tag));
   hls = (TH3F*)hLS->Projection(axis1,axis2,axis3)->Clone(Form("hls3D_%s",tag));
   if(!hus->GetDefaultSumw2()) hus->Sumw2();
   if(!hls->GetDefaultSumw2()) hls->Sumw2();
   hs = (TH3F*)hus->Clone(Form("hsig3D_%s",tag));
   if(!hs->GetDefaultSumw2()) hs->Sumw2();
   hs->Add(hls,-1.);
}


void drawFit(TH1 **h1, TF1 *f1, double *par, double *parErr, TCanvas *c1, TPDF *pdf, int nBins){

   for(int i=0;i<nBins;i++){
      c1->cd(i+1);
      gPad->SetLogy(0);
      f1->SetRange(-0.01,0.01);
      if(h1[i]->GetEntries()>10){
         h1[i]->Fit(f1,"R");
         par[i] = f1->GetParameter(2);
         parErr[i] = f1->GetParError(2);
      }else{
         par[i] = -99;
         parErr[i] = 0;
      }

      h1[i]->SetMarkerSize(0.3);
      h1[i]->SetMarkerStyle(20);
      h1[i]->Draw("");
      TLatex *tt1 = new TLatex();
      tt1->SetNDC(1);
      tt1->DrawLatex(0.2,0.9,Form("nTrksBin %d",i));
   }
   c1->cd();
   pdf->On();
   c1->Update();
   pdf->NewPage();
   pdf->Off();
}

void drawGraph(TGraphErrors *g1, TCanvas *c, TPDF *pdf, TH1 *hFrame){
   c->cd();
   hFrame->Draw();
   g1->Draw("pe1same");
   pdf->On();
   c->Update();
   pdf->NewPage();
   pdf->Off();
}

void drawGraph(TGraphErrors *g1, TCanvas *c, TPDF *pdf, TH1 *hFrame, TF1 *f1){
   c->cd();
   g1->Fit(f1,"R");
   hFrame->Draw();
   g1->Draw("pe1same");
   pdf->On();
   c->Update();
   pdf->NewPage();
   pdf->Off();

}

void effCorrection(TH1F *hRaw, TH1F *hEff, TH1F &*hCor, char *name){
   hCor = (TH1F*)hRaw->Clone(name);
   hCor->Reset();
   for(int j=0;j<hRaw->GetNbinsX();j++){
      double x = hRaw->GetBinCenter(j+1);
      double ex = hRaw->GetBinWidth(j+1);
      double y = hRaw->GetBinContent(j+1);
      double ey = hRaw->GetBinError(j+1);
      int bin = hEff->FindBin(x);
      double exEff = hEff->GetBinWidth(bin);
      if(exEff != ex) cout<<"unequal bin width! "<<ex<<" : "<<exEff<<endl;
      double eff = hEff->GetBinContent(bin);
      double yCor = 0.;
      double eyCor = 0.;
      if(eff>0){ 
         yCor = y/eff;
         eyCor = ey/eff;
      }
      hCor->SetBinContent(j+1,yCor);
      hCor->SetBinError(j+1,eyCor);
   }
}

void getRatioOfHistToGraph(TH1F *h, TH1F *r, TGraphErrors *g){
   if(h->GetNbinsX()!=r->GetNbinsX()){ cout<<"unequal bins! exit!"<<endl; return;}
   for(int i=0;i<h->GetNbinsX();i++){
      double x = h->GetBinCenter(i+1);
      double xg = 0, yg = 0;
      //double dx = 99999;
      //for(int j=0;j<g->GetN();j++){
      //  double xtmp=0, ytmp=0; 
      //  g->GetPoint(j,xtmp,ytmp);
      //  if(fabs(xtmp-x)<dx){
      //     dx = fabs(xtmp-x);
      //     yg = ytmp;
      //     xg = xtmp;
      //  }
      //}
      yg = g->Eval(x);
      double y = h->GetBinContent(i+1);
      double ey = h->GetBinError(i+1);
      if(yg!=0){
         r->SetBinContent(i+1,y/yg);
         r->SetBinError(i+1,ey/yg);
      }else{
         r->SetBinContent(i+1,0);
         r->SetBinError(i+1,0);
      }
   }
}

TGraph *getRatioToGraph(TGraph *g1, TGraphErrors *gref, int flag = 1){
   TGraph *gr;
   if(flag==1) gr = new TGraph(g1->GetN()); 
   else  gr = new TGraph(gref->GetN());

   if(flag==1){
      for(int i=0;i<g1->GetN();i++){
         double x=0,y=0;
         double ex=0,ey=0;
         double ref = 1;
         g1->GetPoint(i,x,y);
         ref = gref->Eval(x);
         if(ref!=0){
            gr->SetPoint(i,x,y/ref);
         }else{
            gr->SetPoint(i,x,0);
         }
      }
   }else{
      for(int i=0;i<gref->GetN();i++){
         double x=0,y=0;
         double ex=0,ey=0;
         double ref = 1;
         gref->GetPoint(i,x,ref);
         y = g1->Eval(x);

         double *xg = g1->GetX();
         int bin = 0;
         for(int j=0;j<g1->GetN()-1;j++){
            if(xg[j]<x&&x<xg[j+1]){
               bin = j; break;
            }
         }
         if(y!=0){
            gr->SetPoint(i,x,y/ref);
         }else{
            gr->SetPoint(i,x,0);
         }
      }
   }

   return gr;
}


TGraphErrors *getRatioToGraph(TGraphErrors *g1, TGraphErrors *gref, int flag = 1){
   TGraphErrors *gr;
   if(flag==1) gr = new TGraphErrors(g1->GetN()); 
   else  gr = new TGraphErrors(gref->GetN());

   if(flag==1){
      for(int i=0;i<g1->GetN();i++){
         double x=0,y=0;
         double ex=0,ey=0;
         double ref = 1;
         g1->GetPoint(i,x,y);
         ref = gref->Eval(x);
         if(ref!=0){
            gr->SetPoint(i,x,y/ref);
            gr->SetPointError(i,g1->GetErrorX(i),g1->GetErrorY(i)/ref);
         }else{
            gr->SetPoint(i,x,0);
            gr->SetPointError(i,g1->GetErrorX(i),0);
         }
      }
   }else{
      for(int i=0;i<gref->GetN();i++){
         double x=0,y=0;
         double ex=0,ey=0;
         double ref = 1;
         gref->GetPoint(i,x,ref);
         y = g1->Eval(x);

         double *xg = g1->GetX();
         int bin = 0;
         for(int j=0;j<g1->GetN()-1;j++){
            if(xg[j]<x&&x<xg[j+1]){
               bin = j; break;
            }
         }
         if(y!=0){
            gr->SetPoint(i,x,y/ref);
            gr->SetPointError(i,gref->GetErrorX(i),g1->GetErrorY(bin)/ref);
         }else{
            gr->SetPoint(i,x,0);
            gr->SetPointError(i,gref->GetErrorX(i),0);
         }
      }
   }

   return gr;
}

TGraphAsymmErrors *getRatioToGraph(TGraphAsymmErrors *g1, TGraphErrors *gref){
   TGraphAsymmErrors *gr = new TGraphAsymmErrors(g1->GetN()); 
   for(int i=0;i<g1->GetN();i++){
      double x=0,y=0;
      double exl=0,exh=0,eyl=0,eyh=0;
      g1->GetPoint(i,x,y);
      double ref = gref->Eval(x);
      exl = g1->GetErrorXlow(i);
      exh = g1->GetErrorXhigh(i);
      eyl = g1->GetErrorYlow(i);
      eyh = g1->GetErrorYhigh(i);
      if(ref!=0){
         gr->SetPoint(i,x,y/ref);
         gr->SetPointError(i,exl,exh,eyl/ref,eyh/ref);
      }else{
         gr->SetPoint(i,x,0);
         gr->SetPointError(i,exl,exh,0,0);
      }
   }
   return gr;
}

TGraphAsymmErrors *getRatioToGraphRefBin(TGraph *g1, TGraphErrors *gref){
   TGraphAsymmErrors *gr = new TGraphAsymmErrors(gref->GetN()); 
   for(int i=0;i<gref->GetN();i++){
      double x=0,y=0;
      double ref=0;
      gref->GetPoint(i,x,ref);
      double y = g1->Eval(x);
      double ex = gref->GetErrorX(i);
      double *xg = g1->GetX();
      int bin = 0;
      for(int j=0;j<g1->GetN()-1;j++){
         if(xg[j]<=x&&x<xg[j+1]){
            bin = j; break;
         }
      }
      if(ref!=0){
         gr->SetPoint(i,x,y/ref);
      }else{
         gr->SetPoint(i,x,0);
      }
   }
   return gr;
}

TGraphAsymmErrors *getRatioToGraphRefBin(TGraphAsymmErrors *g1, TGraphErrors *gref){
   TGraphAsymmErrors *gr = new TGraphAsymmErrors(gref->GetN()); 
   for(int i=0;i<gref->GetN();i++){
      double x=0,y=0;
      double exl=0,exh=0,eyl=0,eyh=0;
      double ref=0;
      gref->GetPoint(i,x,ref);
      double y = g1->Eval(x);
      double ex = gref->GetErrorX(i);
      double *xg = g1->GetX();
      int bin = -1;
      for(int j=0;j<g1->GetN()-1;j++){
         if(xg[j]<=x&&x<xg[j+1]){
            bin = j; break;
         }
      }
      if(bin>=0){
         exl = g1->GetErrorXlow(bin);
         exh = g1->GetErrorXhigh(bin);
         eyl = g1->GetErrorYlow(bin);
         eyh = g1->GetErrorYhigh(bin);
         if(ref!=0){
            gr->SetPoint(i,x,y/ref);
            gr->SetPointError(i,ex,ex,eyl/ref,eyh/ref);
         }else{
            gr->SetPoint(i,x,-1);
            gr->SetPointError(i,ex,ex,0,0);
         }
      }else{
         gr->SetPoint(i,x,-1);
         gr->SetPointError(i,ex,ex,0,0);
      }
   }
   return gr;
}


//-----------------------------//
//        fit functions        //
//-----------------------------//

Double_t tofEff(Double_t *x, Double_t *par)
{
   return par[0]*exp(-pow(par[1]/x[0],par[2]))+par[3]*exp(-pow((x[0]-par[4])/par[5],2));
}

Double_t tpcEff(Double_t *x, Double_t *par)
{
   return par[0]*exp(-pow(par[1]/x[0],par[2]));
   //return 0.8304*exp(-pow(0.07207/x[0],1.893))+0.07667*exp(-pow((x[0]-0.3166)/0.5872,2));
}
Double_t tpcEffRc(Double_t *x, Double_t *par)
{
   return par[0]*exp(-pow(par[1]/x[0],par[2]))+par[3]*exp(-pow((x[0]-par[4])/par[5],2));
}

Double_t nSigEEff(Double_t *x, Double_t *par)
{
   return par[0]/(exp((x[0]-par[1])/par[2])+par[3])+par[4];
}

//------------------------------------------------------ //bremsstrahlung tail
Double_t CrystalBall(Double_t *x, Double_t *par)
{
   Double_t N = par[0];
   Double_t mu = par[1];
   Double_t s = par[2];
   Double_t n = par[3];
   Double_t alpha1 = par[4];

   Double_t A = TMath::Power(n/fabs(alpha1), n) * TMath::Exp(-alpha1*alpha1/2.);
   Double_t B = n/fabs(alpha1) - fabs(alpha1);
   Double_t norm = (x[0]-mu)/s;

   if(norm > -alpha1) {
      return N * TMath::Exp(-0.5*norm*norm);
   } else {
      return N * A * TMath::Power(B-norm, -n);
   }
}


//------------------------------------------------------
Double_t CrystalBall2(Double_t *x, Double_t *par)
{
   Double_t N = par[0];
   Double_t mu = par[1];
   Double_t s = par[2];
   Double_t n = par[3];
   Double_t alpha2 = par[4];
   Double_t m = par[5];
   Double_t beta = par[6];

   Double_t A = TMath::Power(n/fabs(alpha2), n) * TMath::Exp(-alpha2*alpha2/2.);
   Double_t B = n/fabs(alpha2) - fabs(alpha2);

   Double_t C = TMath::Power(m/fabs(beta), m) * TMath::Exp(-beta*beta/2.);
   Double_t D = m/fabs(beta) - fabs(beta);

   Double_t norm = (x[0]-mu)/s;

   if(norm < -alpha2) {
      return N * A * TMath::Power(B-norm, -n);
   } else if(norm < beta) {
      return N * TMath::Exp(-0.5*norm*norm);
   } else {
      return N * C * TMath::Power(D+norm, -m);
   }
}

//------------------------------------------------------
Double_t momRes(Double_t *x, Double_t *par)
{
   double pt = x[0];
   const Double_t a = par[0];
   const Double_t b = par[1];
   return TMath::Sqrt(a*a*pt*pt+b*b);
}

Double_t ellipsoid(Double_t *x Double_t *par){
   return par[0] + sqrt(par[1]+x[0]*x[0]) + par[2]*x[1]*x[1]/(1.+sqrt(1.-par[3]*x[1]*x[1]));
}
