TH1D *binShiftCorrection(TH1D *h){
  
  TF1 *fpow = new TF1("fpow","[0]*pow(x,[1])",70,200);
  fpow->SetLineWidth(1);

  string hName = h->GetName();
  //cout<<" bin shift correcting:  "<<hName<<endl;
  string hTempName = hName+"Orig";

  TH1D *hOrig=h->Clone(hTempName.c_str());
  
  for(int iter=0;iter<4;iter++){
    
    //cout<<" iteration # "<<iter<<endl;
    h->Fit(fpow,"","N",70,250);
    for(int i=0;i<h->GetNbinsX();i++){
      float meanVal = fpow->Integral(h->GetBinLowEdge(i+1),h->GetBinLowEdge(i+1)+h->GetBinWidth(i+1))/h->GetBinWidth(i+1);
      float centVal = fpow->Eval(h->GetBinCenter(i+1));
      float binShiftCorr = centVal/meanVal;
      //cout<<" i "<<i<<" corr "<<binShiftCorr<<endl;

      float val =     hOrig->GetBinContent(i+1);
      float err =     hOrig->GetBinError(i+1);

      h->SetBinContent(i+1,val*binShiftCorr);
      h->SetBinError(i+1,err*binShiftCorr);

    }
  }

  delete fpow;
  delete hOrig;
  return h;
}

TH1D *binWidthCorrection(TH1D *h){
  for(int i=0;i<h->GetNbinsX();i++){
   float val =     h->GetBinContent(i+1);
   float err =     h->GetBinError(i+1);
   float width =   h->GetBinWidth(i+1);
   h->SetBinContent(i+1,val/width);
   h->SetBinError(i+1,err/width);
 }
 return h;
}


TH1D *extractSys(TH1D *h, TH1D *hDef, string nameIndex){

  string sysName  = "hSys"+nameIndex;
  TH1D *hSys = (TH1D*)hDef->Clone(sysName.c_str());
  hSys->Reset();

  for(int i=0;i<hDef->GetNbinsX();i++){

    float defVal = hDef->GetBinContent(i+1);
    float binCenter = hDef->GetBinCenter(i+1);
    if(binCenter < h->GetBinLowEdge(1)) continue;
    
    float binContent = h->GetBinContent(h->FindBin(binCenter));
    
    hSys->SetBinContent(i+1,defVal*(binContent+1.));
    /*
    cout<<" default Val "<<defVal<<endl;
    cout<<" binContent "<<binContent<<endl;
    cout<<" binCenter "<<binCenter<<endl;
    cout<<(binContent+1.)*defVal<<endl;
    */

  }

  return hSys;
}

void drawRAA(int writeRaw=0, int grabUnfolded=0, int doBjets=0, int fixedBins=0){

  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gROOT->ForceStyle(1);
  
  TFile *fPbPb, *fpp, *fUnfolded;
  TTree *tPbPb, *tpp;
  TH1D *hPbPb, *hpp;

  if(!grabUnfolded){
    
    fPbPb = new TFile("../histos/PbPbdata_pt30by3_jpHICalibRepass_withDup_PU_jet556580.root");
    fpp = new TFile("../histos/ppdata_ppReco_ak3PF_jetTrig_noIPupperCut.root");
    
    tPbPb = (TTree *)fPbPb->Get("nt");
    tpp = (TTree *)fpp->Get("nt");


    const int nBins = 7;
    double ptBin[nBins+1] = {60,70,80,90,110,130,170,250};
    
    if(fixedBins){
      hPbPb = new TH1D("hPbPb","PbPb Spectrum; p_{T} (GeV/c); d#sigma/dp_{T} [mb (GeV/c)];",21,40,250);
      hpp = new TH1D("hpp","pp Spectrum; p_{T} (GeV/c); d#sigma/dp_{T} [mb (GeV/c)];",21,40,250);
    }
    else{
      hPbPb = new TH1D("hPbPb","PbPb Spectrum; p_{T} (GeV/c); d#sigma/dp_{T} [mb (GeV/c)];",nBins,ptBin);
      hpp = new TH1D("hpp","pp Spectrum; p_{T} [GeV/c]; Counts",nBins,ptBin);
    }
    hPbPb->Sumw2();
    hpp->Sumw2();
    
    cout<<" drawing trees "<<endl;
    tPbPb->Draw("jtptA>>hPbPb","weight");
    tpp->Draw("jtpt>>hpp","weight","same");
    cout<<" done "<<endl;
  }    
  else{
    cout<<" grabbing unfolded spectra "<<endl;
    TFile *fUnfolded = NULL;    
    if(doBjets) fUnfolded = new TFile("pbpb_Unfo_akPu3PF_jtpt60_bJets_v4.root");
    else{
      if(fixedBins)fUnfolded = new TFile("~/Work/bTagging/bJetRAA/pbpb_Unfo_akPu3PF_jtpt60_Inc_v4_fixedBins.root");
      else fUnfolded = new TFile("pbpb_Unfo_akPu3PF_jtpt60_Inc_v4_anaBins.root");
    }
    hPbPb = (TH1D*)fUnfolded->Get("hReco_cent0");
    hpp = (TH1D*)fUnfolded->Get("hReco_cent1");        
    //hPbPb = (TH1D*)fUnfolded->Get("UnfoldedBinByBin_cent0");
    //hpp = (TH1D*)fUnfolded->Get("UnfoldedBinByBin_cent1");        
  }
  
 if(writeRaw&&!grabUnfolded&&!doBjets){
   cout<<" writing raw files "<<endl;
  TFile *fOut=NULL;
   if(fixedBins) fOut=new TFile("rawIncSpectra_MB_fixedBinSize_v4.root","recreate");
   else fOut=new TFile("rawIncSpectra_MB_varBinSize_v4.root","recreate");
   hPbPb->Write();
   hpp->Write();
   return;
 }
  

  
  hPbPb->SetLineColor(2);
  
  
  hPbPb->GetXaxis()->SetNdivisions(505);
  
  hPbPb->SetMarkerStyle(8);
  hpp->SetMarkerStyle(8);
  hpp->SetMarkerColor(2);
  
  
  
  
  TCanvas *c1=new TCanvas("c1","c1",600,600);
  
 
  cout<<" applying scale factors "<<endl;
  hPbPb->Scale(1./5.67);  // taa
  //hPbPb->Scale(1./25.9);  // taa
  //hPbPb->Scale(1./0.17);  // taa
  float nMB = 1.09079385600000000e+09;
  //float nMB = 1.09079385600000000e+09*0.05;
  //float nMB = 1.09079385600000000e+09*0.2;
  hPbPb->Scale(1./nMB); //nMB
  hpp->Scale(1./5.33e9);  // luminosity

  


  cout<<" bin width correction "<<endl;
  // divide out the bin-width   
  hPbPb = binWidthCorrection(hPbPb);
  hpp = binWidthCorrection(hpp);

  cout<<" jet55 trig eff correction "<<endl;

 // correct trig eff:
  if(!grabUnfolded){ 
    float trigEff[7]={
      0.966051,
      0.99591,
      0.999665,
      1.,
      1.,
      1.,
      1.};

    float trigEffSmallBins[21]={
      0.523891,
      0.814003,
      0.966051,
      0.99591	,
      0.999665,
      1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};
    
    for(int i=0;i<hPbPb->GetNbinsX();i++){
      float val = hPbPb->GetBinContent(i+1);
      float err = hPbPb->GetBinError(i+1);

      float trigEffBin = 0.;
      if(fixedBins) trigEffBin = trigEffSmallBins[i];
      else trigEffBin = trigEff[i];
      float newVal = val/trigEffBin;
      float newErr = err/trigEffBin;
            
      hPbPb->SetBinContent(i+1,newVal);
      hPbPb->SetBinError(i+1,newErr);
      
    }
  }

  cout<<" bin shift correction "<<endl;
  hPbPb = binShiftCorrection(hPbPb);
  hpp = binShiftCorrection(hpp);


 
  hPbPb->Draw();
  hpp->Draw("same");
  

  TLegend *leg=new TLegend(0.6,0.6,0.9,0.9);
  leg->AddEntry(hPbPb,"PbPb / T_{AA}","p");
  leg->AddEntry(hpp,"pp","p");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();


  TCanvas *c2=new TCanvas("c2","c2",600,600);
  TH1D *hRatio = hPbPb->Clone("hPbPb");
  hRatio->Reset();
  hRatio->SetYTitle("R_{AA}");
  //hRatio->Divide(hPbPb,hpp,1.,1.,"B");
  hRatio->Divide(hPbPb,hpp);

  hRatio->GetYaxis()->SetRangeUser(0.,1.);

 

 
 
  hRatio->GetXaxis()->SetRangeUser(70,250);
  hRatio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hRatio->Draw();
  /*
  cout<<" numerator error "<<hPbPb->GetBinError(17)<<"/"<<hPbPb->GetBinContent(17)<<" = "<<hPbPb->GetBinError(17)/hPbPb->GetBinContent(17)<<endl;
  cout<<" denominator error "<<hpp->GetBinError(17)<<"/"<<hpp->GetBinContent(17)<<" = "<<hpp->GetBinError(17)/hpp->GetBinContent(17)<<endl;
  cout<<" ratio error "<<hRatio->GetBinError(17)<<"/"<<hRatio->GetBinContent(17)<<" = "<<hRatio->GetBinError(17)/hRatio->GetBinContent(17)<<endl;
  */

  cout<<" writing output file "<<endl;

  string outFileName;
  if(!doBjets){ 
    if(grabUnfolded){
      if(fixedBins)outFileName = "unfIncJetRAA_fixedBins.root";
      else outFileName = "unfIncJetRAA_anaBins.root";
    }
    else {
      if(fixedBins)outFileName = "rawIncJetRAA_fixedBins.root";
      else outFileName = "rawIncJetRAA_anaBins.root";
    }
  }
  else{
    if(grabUnfolded){
      outFileName = "unfBJetRAA.root";
    }
    else{
      cout<<" not writing output "<<endl;
      return;
    }
  }

  // Add systematics
/*
  TCanvas *c4=new TCanvas("c4","c4",600,600);

  if(doBjets){
    
    TFile *fSysPbPb = new TFile("../systematics_components_PbPb.root");
    TFile *fSysPP  = new TFile("../systematics_components_pp.root");

    TH1D *hFracLTJP_PbPb  = fSysPbPb->Get("hFracLTJP");
    TH1D *hFracWorkingPointUp_PbPb  = fSysPbPb->Get("hFracWorkingPointUp");
    TH1D *hFracWorkingPointDown_PbPb  = fSysPbPb->Get("hFracWorkingPointDown");
    TH1D *hFracDataDriven_PbPb  = fSysPbPb->Get("hFracDataDriven");
    TH1D *hFracCharm_PbPb  = fSysPbPb->Get("hFracCharm");

    TH1D *hFracLTJP_PP  = fSysPP->Get("hFracLTJP");
    TH1D *hFracWorkingPointUp_PP  = fSysPP->Get("hFracWorkingPointUp");
    TH1D *hFracWorkingPointDown_PP  = fSysPP->Get("hFracWorkingPointDown");
    TH1D *hFracDataDriven_PP  = fSysPP->Get("hFracDataDriven");
    TH1D *hFracCharm_PP  = fSysPP->Get("hFracCharm");

    TH1D *hSysPbPb1 =extractSys(hFracLTJP_PbPb, hPbPb, "PbPb1"); 
    TH1D *hSysPbPb2a =extractSys(hFracWorkingPointUp_PbPb, hPbPb, "PbPb2a"); 
    TH1D *hSysPbPb2b =extractSys(hFracWorkingPointDown_PbPb, hPbPb, "PbPb2b"); 
    TH1D *hSysPbPb3 =extractSys(hFracDataDriven_PbPb, hPbPb, "PbPb3"); 
    TH1D *hSysPbPb4 =extractSys(hFracCharm_PbPb, hPbPb, "PbPb5"); 

    TH1D *hSysPP1 =extractSys(hFracLTJP_PP, hpp, "PP1"); 
    TH1D *hSysPP2a =extractSys(hFracWorkingPointUp_PP, hpp, "PP2a"); 
    TH1D *hSysPP2b =extractSys(hFracWorkingPointDown_PP, hpp, "PP2b"); 
    TH1D *hSysPP3 =extractSys(hFracDataDriven_PP, hpp, "PP3"); 
    TH1D *hSysPP4 =extractSys(hFracCharm_PP, hpp, "PP4"); 


    TH1D *hRatio1 = hSysPbPb1->Clone("hRatio1");
    TH1D *hRatio2a = hSysPbPb2a->Clone("hRatio2a");
    TH1D *hRatio2b = hSysPbPb2b->Clone("hRatio2b");
    TH1D *hRatio3 = hSysPbPb3->Clone("hRatio3");
    TH1D *hRatio4 = hSysPbPb4->Clone("hRatio4");

    hRatio1->Divide(hSysPbPb1,hSysPP1);
    hRatio2a->Divide(hSysPbPb2a,hSysPP2a);
    hRatio2b->Divide(hSysPbPb2b,hSysPP2b);
    hRatio3->Divide(hSysPbPb3,hSysPP3);
    hRatio4->Divide(hSysPbPb4,hSysPP4);

    hRatio1->SetLineColor(2);
    hRatio2a->SetLineColor(3);
    hRatio2b->SetLineColor(4);
    hRatio3->SetLineColor(6);
    hRatio4->SetLineColor(8);

    hRatio->Draw();
    hRatio1->Draw("same");
    hRatio2a->Draw("same");
    hRatio2b->Draw("same");
    hRatio3->Draw("same");
    hRatio4->Draw("same");

    float xVal[7], yVal[7], xErr[7], yErr[7];

    int j=0;

    for(int i=0;i<hRatio->GetNbinsX();i++){
      float defVal = hRatio->GetBinContent(i+1);
      
      if(hRatio->GetBinLowEdge(i+1)<60) continue;
      
      float sys1 = fabs(defVal-hRatio1->GetBinContent(i+1));
      float sys2a = fabs(defVal-hRatio2a->GetBinContent(i+1));
      float sys2b = fabs(defVal-hRatio2b->GetBinContent(i+1));
      float sys3 = fabs(defVal-hRatio3->GetBinContent(i+1));
      float sys4 = fabs(defVal-hRatio4->GetBinContent(i+1));
      float sys2 = TMath::Max(sys2a, sys2b);
   
      float sysTot = sqrt(sys1*sys1 + sys2*sys2 + sys3*sys3 + sys4*sys4);
      cout<<" sysTot "<<sysTot<<endl;

      xVal[j]=hRatio->GetBinCenter(i+1);
      yVal[j]=defVal;
      yErr[j]=sysTot;
      xErr[j]=hRatio->GetBinWidth(i+1)/2.;
      j++;
    }

    TGraphErrors *gSys=new TGraphErrors(7,xVal,yVal,xErr,yErr);
    gSys->SetFillColor(5);
    gSys->Draw("2");
    gSys->SetName("gSys");
    hRatio->Draw("same");
    
  }
*/

  TFile *fout=new TFile(outFileName.c_str(),"recreate");
  hRatio->SetName("hRatio");
  hPbPb->Write();
  hpp->Write();
  hRatio->Write();
//  if(doBjets) gSys->Write();
  fout->Close();    


}
