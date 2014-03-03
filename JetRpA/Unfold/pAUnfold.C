#include <iostream>
#include <stdio.h>

#include <TRandom2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>

#include "utilities.h"
#include "bayesianUnfold.h"
#include "prior.h"

using namespace std;


//==============================================================================
// Unfolding Ying Lu 08 07 11
// Update Yen-Jie Lee 06.22.12
//==============================================================================

void pAUnfold(int algo= 3,bool useSpectraFromFile=0, bool useMatrixFromFile=0, int doToy = 0, int isMC = 0,char *spectraFileName = (char*)"ppb_spectra_akPu3PF.root",double recoJetPtCut = 30.,double trackMaxPtCut = 0, int nBayesianIter = 4, int doBjets=0, int doTrigEffCorr=0) // algo 2 =akpu2 ; 3 =akpu3 ; 4 =akpu4 ;1 = icpu5
{
  
  gStyle->SetErrorX(0.);
  gStyle->SetPaintTextFormat("3.2f");
  gStyle->SetOptLogz(1);
  gStyle->SetPadRightMargin(0.13);	
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);


  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  const float ppsigma=70.;  //units in mb
  const float nMB = 2.60259184640000000e+10 ; 
  const float Ncoll = 6.9 ; // for inclusive pPb 
//  const float Ncoll = 7.5 ; // for pPb in 0-90% 
  const float etamin = -1.0 ;
  const float etamax = 1.0 ;
  TString coll = "PPb" ; 
  const bool SavePlot=kFALSE;  
  // input files
  char *fileNamePP_mc = NULL;
  char *fileNamePbPb_mc = NULL;
  char *fileNamePP_data = NULL;
  char *fileNamePbPb_data = NULL;
  

  // pp file needs replacing
  if(doBjets)fileNamePP_data = (char*)"/net/hidsk0001/d00/scratch/maoyx/Btag/Unfold/bJetRAA/NewFormatV5_bFractionMCTemplate_pppp1_SSVHEat2.0FixCL0_bin_0_40_eta_0_2.root";
//  else fileNamePP_data = (char*)"/afs/cern.ch/work/y/ymao/analysis/AsymmetryPA/Unfold/histos/PPbJetTrigHFsumEta4Bin1PYTHIAGenLevelAkPu3PFJetSpectra.root" ;
//  else fileNamePP_data = (char*)"/afs/cern.ch/work/y/ymao/analysis/AsymmetryPA/Unfold/histos/PPbJetTrigHFsumEta4Bin1PYTHIAAkPu3PFJetSpectra.root" ;
  else fileNamePP_data = (char*)Form("/afs/cern.ch/work/y/ymao/analysis/AsymmetryPA/Unfold/histos/ppMCProd16_ppReco_%s_QCDjetTrig_JetPt0noIPupperCut.root",algoName[algo]) ;
//  else fileNamePP_data = (char*)"/afs/cern.ch/work/y/ymao/analysis/AsymmetryPA/Unfold/histos/ppMC_ppReco_akPu3PF_QCDjetTrig_noIPupperCut.root" ;
  //if(doBjets)fileNamePbPb_data = (char*)"~/Work/bTagging/outputTowardsFinal/AltBinningV6_bFractionMCTemplate_ppPbPb1_SSVHEat2.0FixCL0_bin_0_40_eta_0_2.root";
  if(doBjets)fileNamePbPb_data = (char*)"/net/hidsk0001/d00/scratch/maoyx/Btag/Unfold/bJetRAA/bFractionMCTemplate_ppPbPb1_SSVHEat2.0FixCL0_bin_0_40_eta_0_2_binomErrors_jet55_wideBin_v2.root";
  //else fileNamePbPb_data = (char*)Form("/afs/cern.ch/work/y/ymao/analysis/AsymmetryPA/Unfold/histos/%sCombinedJetTrig%sJetNoResidualTrkEffHIN12017v5TrkCorr2DCutAllHistHFsumEta4Bin1.root",coll.Data(), algoName[algo]) ;
  else fileNamePbPb_data = (char*)Form("/afs/cern.ch/work/y/ymao/analysis/AsymmetryPA/Unfold/histos/%sCombinedJetTrig%sJetTrkEffHIN12017v5TrkCorr2DCutAllHistHFsumEta4Bin1.root",coll.Data(), algoName[algo]) ;
//  else fileNamePbPb_data = (char*)"/afs/cern.ch/work/y/ymao/analysis/AsymmetryPA/Unfold/histos/PPbJetTrigHF90CentMBAllPYTHIAAkPu3PFJetNotNormalized.root" ; 
//   else fileNamePbPb_data = (char*)"/afs/cern.ch/work/y/ymao/analysis/AsymmetryPA/Unfold/histos/pPbdata_ppReco_akPu3PF_jetTrig40_60_100_PtCut0noIPupperCut.root" ;
  if(doBjets) fileNamePP_mc = (char*)"/net/hidsk0001/d00/scratch/kjung/histos/ppMC_ppReco_ak3PF_BjetTrig_noIPupperCut.root";
  else fileNamePP_mc = (char*)Form("/afs/cern.ch/work/y/ymao/analysis/AsymmetryPA/Unfold/histos/ppMCProd16_ppReco_%s_QCDjetTrig_JetPt0noIPupperCut.root", algoName[algo]) ;
  if(doBjets)fileNamePbPb_mc = (char*) "/net/hisrv0001/home/mnguyen/scratch/bTaggingOutput/ntuples/PbPbBMC_pt30by3_ipHICalibCentWeight_noTrig.root";
  else {
    if(coll=="PPb")fileNamePbPb_mc = (char*) Form("/afs/cern.ch/work/y/ymao/analysis/AsymmetryPA/Unfold/histos/pPbMCProd16_ppReco_%s_QCDjetTrig_JetPt0noIPupperCut.root", algoName[algo]) ;
//  else fileNamePbPb_mc = (char*) Form("/afs/cern.ch/work/y/ymao/analysis/AsymmetryPA/Unfold/histos/PbpMCProd24_ppReco_%s_QCDjetTrig_JetPt0noIPupperCut.root", algoName[algo]) ;
  else fileNamePbPb_mc = (char*) Form("/afs/cern.ch/work/y/ymao/analysis/AsymmetryPA/Unfold/histos/pPbMCProd16_ppReco_%s_QCDjetTrig_JetPt0noIPupperCut.root", algoName[algo]) ;
    }


  // grab ntuples
//  TFile *infPbPb_mc = new TFile(fileNamePbPb_mc);
//  TFile *infPP_mc = new TFile(fileNamePP_mc);
  

  string bJetString = "Inc";
  if(doBjets) bJetString = "bJets";

  // Output file
  TFile *pbpb_Unfo;
  if (isMC) pbpb_Unfo = new TFile(Form("/afs/cern.ch/work/y/ymao/analysis/AsymmetryPA/Unfold/histos/%s_UnfoPriorGen_%sNoResidual_MC_jtpt%.0f_EtaBin%.f_%.f_%s_v2.root",coll.Data(), algoName[algo],recoJetPtCut,etamin*10, etamax*10, bJetString.c_str()),"RECREATE");
  else pbpb_Unfo  = new TFile(Form("/afs/cern.ch/work/y/ymao/analysis/AsymmetryPA/Unfold/histos/%s_UnfoPriorGen_%s_jtpt%.0f_EtaBin%.f_%.f_%s_v5.root",coll.Data(), algoName[algo],recoJetPtCut,etamin*10, etamax*10, bJetString.c_str()),"RECREATE");

// Histograms used by RooUnfold
  UnfoldingHistos *uhist[nbins_cent+1];
		
  // Initialize Histograms   
	
  for (int i=0;i<=nbins_cent;i++) uhist[i] = new UnfoldingHistos(i);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
	
  // Initialize reweighting functions
  /*
  TCut dataSelection;
  TCut dataSelectionPP;
  TCut TriggerSelectionPP;
  TCut TriggerSelectionPbPb80;

  if(doBjets)dataSelection = "weight*(abs(refparton_flavorForB)==5&&abs(jteta)<2)";
  else dataSelection = "weight*(abs(jteta)<2)";
  */

  if (isMC) cout<<"This is a MC closure test"<<endl;
  else cout<< "This is a data analysis"<<endl;    		     
	     	
  // Setup jet data branches, basically the jet tree branches are assigned to this object when we loop over the events
	
  JetDataPbPb *dataPbPb   = new JetDataPbPb(fileNamePbPb_mc,(char*)"nt"); // PbPb data	
  JetDataPP *dataPP = new JetDataPP(fileNamePP_mc,(char*)"nt");	// pp data
	
  TFile *fSpectra(0);		
  if (useSpectraFromFile||useMatrixFromFile){
    fSpectra = new TFile(spectraFileName,"read");
  }
  
  // Come back to the output file dir
  pbpb_Unfo->cd();

  cout <<"MC = " << isMC <<endl;	
  
  TH1F *hCent = new TH1F("hCent","",nbins_cent,boundaries_cent);
 

  // if you change the binning you have to change these, too
  // inclusive trig eff
  /*
    float trigEffInc[6]={0.777265,
			 0.95765,
			 0.998357,
			 0.999941,
			 1.,
			 1.};
  */

    // b-jet trig eff
    /*
    float trigEffbJet[6]={0.660276,
		      0.908997,
		      0.980793,
		      0.998767,
		      0.999442,
		      1.};
    */



 
		
  // Fill pPb MC
   double etashift = 0. ;
   if(coll=="PbP") etashift = 0.465 ;
   else etashift = 0.465 ;   
  if (!useMatrixFromFile) {
    for (Long64_t jentry2=0; jentry2<dataPbPb->tJet->GetEntries();jentry2++) {
      dataPbPb->tJet->GetEntry(jentry2);
      
      // change when we switch to centrality binning
      int cBin = 0;
            
      if ( dataPbPb->refpt  < 0. ) continue;
      if ( (dataPbPb->jteta+etashift)  > etamax || (dataPbPb->jteta+etashift) < etamin ) continue;
    //  if ( fabs(dataPbPb->jteta+0.465)  > 1. ) continue;
      if ( dataPbPb->refpt<0) dataPbPb->refpt=0;
      if (doBjets && fabs(dataPbPb->refparton_flavorForB)!=5) continue;
      //if (doBjets&& dataPbPb->discr_ssvHighEff<2) continue;
      if (doBjets && dataPbPb->jtptB < recoJetPtCut) continue;
    //  if (!doBjets && dataPbPb->jtpt < recoJetPtCut) continue;
      //if (!doTrigEffCorr && dataPbPb->isTrig <1) continue;
     // if ( dataPbPb->isTrig <1) continue;
      
      //if(!doBjets)if(dataPbPb->refpt < 50 && dataPbPb->jtptA>120) continue;
      //if(doBjets)if(dataPbPb->refpt < 50 && dataPbPb->jtptB>120) continue;

      //for centrality selection using HF sum energy with |eta|>4
   //   if((dataPbPb->hiHFplusEta4+dataPbPb->hiHFminusEta4)<2.87) continue ;

      float weight = dataPbPb->weight;

      if(doTrigEffCorr){
	for(int iBin = 0; iBin<nbins_rec; iBin++){
	  float myJtPt = 0.;
	  if(doBjets) myJtPt = dataPbPb->jtptB;
	  else myJtPt = dataPbPb->jtpt;
	  if(myJtPt > boundaries_rec[iBin] && myJtPt < boundaries_rec[iBin+1]){
	    if(doBjets) weight/= trigEffbJet[iBin];
	    else weight/= trigEffInc[iBin];
	  }							  
	}
      }
     if(isMC) {
      if (jentry2 % 2 == 1) {
	if(doBjets)uhist[cBin]-> hMatrix->Fill(dataPbPb->refpt,dataPbPb->jtptB,weight);
	else uhist[cBin]-> hMatrix->Fill(dataPbPb->refpt,dataPbPb->jtpt,weight);
      }	  
      if (jentry2 % 2 == 0) {
	uhist[cBin]-> hGen->Fill(dataPbPb->refpt,weight);   
	//if(doBjets)uhist[cBin]-> hMeas->Fill(dataPbPb->jtptB,weight);  	 
	//else 
          uhist[cBin]-> hMeas->Fill(dataPbPb->jtpt,weight);  	 
	//uhist[cBin]-> hMeasJECSys->Fill(dataPbPb->jtpt*(1.+0.02/nbins_cent*(nbins_cent-i)),weight); 
	// FIXME!!!!!!  i is supposed to be a loop over centrality !!!!
	if(doBjets)uhist[cBin]-> hMeasJECSys->Fill(dataPbPb->jtptB*(1.+0.02/nbins_cent*(nbins_cent-0)),weight); 
	else uhist[cBin]-> hMeasJECSys->Fill(dataPbPb->jtpt*(1.+0.02/nbins_cent*(nbins_cent-0)),weight); 
      }
   }
  else {
      uhist[cBin]-> hMatrix->Fill(dataPbPb->refpt,dataPbPb->jtpt,weight);
      uhist[cBin]-> hGen->Fill(dataPbPb->refpt,weight);
      uhist[cBin]-> hMeas->Fill(dataPbPb->jtpt,weight);
     uhist[cBin]-> hMeasJECSys->Fill(dataPbPb->jtpt*(1.+0.02/nbins_cent*(nbins_cent-0)),weight);
      } 
    } //PbPb entries loop 

    //pp will just fill the last index of the centrality array

    // fill pp MC
    for (Long64_t jentry2=0; jentry2<dataPP->tJet->GetEntries();jentry2++) {
      dataPP->tJet->GetEntry(jentry2);
      
      if ( dataPP->refpt<0) continue;
      if ( (dataPP->jteta+0.465)  > etamax || (dataPP->jteta+0.465) < etamin ) continue;
   //    if ( fabs(dataPP->jteta+0.465)  > 1. ) continue; 
     if ( dataPP->refpt<0) dataPP->refpt=0;
      if ( doBjets && fabs(dataPP->refparton_flavorForB)!=5) continue;
      //if ( doBjets && dataPP->discr_ssvHighEff<2) continue;
    //  if ( dataPP->jtpt < recoJetPtCut) continue;
      
      if(isMC){ 
      if (jentry2 % 2 == 1) {
	uhist[nbins_cent]-> hMatrix->Fill(dataPP->refpt,dataPP->jtpt,dataPP->weight);
      }	  
      if (jentry2 % 2 == 0) {
	uhist[nbins_cent]-> hGen->Fill(dataPP->refpt,dataPP->weight);   
//	uhist[nbins_cent]-> hGen->Fill(dataPP->refpt);   
	uhist[nbins_cent]-> hMeas->Fill(dataPP->jtpt,dataPP->weight); 
      }        
   }
  else  {
      uhist[nbins_cent]-> hMatrix->Fill(dataPP->refpt,dataPP->jtpt,dataPP->weight);
      uhist[nbins_cent]-> hGen->Fill(dataPP->refpt,dataPP->weight);      
   //   uhist[nbins_cent]-> hGen->Fill(dataPP->refpt);      
      uhist[nbins_cent]-> hMeas->Fill(dataPP->jtpt,dataPP->weight);  
     }   
    } //pp entries loop
  } //Fill matrix for pp and PbPb
  /*
  for (int i=0;i<=nbins_cent;i++){
    for (int x=1;x<=uhist[i]->hMatrix->GetNbinsX();x++) {
	float binContent = uhist[i]->hGen->GetBinContent(x);
	float binError = uhist[i]->hGen->GetBinError(x);
	float binWidth =  uhist[i]->hGen->GetXaxis()->GetBinWidth(x);
	uhist[i]->hGen->SetBinContent(x,binContent/binWidth);
	uhist[i]->hGen->SetBinError(x,binError/binWidth);
    }
  }

<<<<<<< HEAD
  for (int i=0;i<=nbins_cent;i++){
    for (int x=1;x<=uhist[i]->hMatrix->GetNbinsX();x++) {
	float binContent = uhist[i]->hMeas->GetBinContent(x);
	float binError = uhist[i]->hMeas->GetBinError(x);
	float binWidth =  uhist[i]->hMeas->GetXaxis()->GetBinWidth(x);
	uhist[i]->hMeas->SetBinContent(x,binContent/binWidth);
	uhist[i]->hMeas->SetBinError(x,binError/binWidth);
    }
  }

  for (int i=0;i<=nbins_cent;i++){
    for (int x=1;x<=uhist[i]->hMatrix->GetNbinsX();x++) {
      for (int y=1;y<=uhist[i]->hMatrix->GetNbinsY();y++) {
	float binContent = uhist[i]->hMatrix->GetBinContent(x,y);
	float binError = uhist[i]->hMatrix->GetBinError(x,y);
	float binWidth2 =  uhist[i]->hMatrix->GetXaxis()->GetBinWidth(x)*uhist[i]->hMatrix->GetYaxis()->GetBinWidth(y);
	uhist[i]->hMatrix->SetBinContent(x,y,binContent/binWidth2);
	uhist[i]->hMatrix->SetBinError(x,y,binError/binWidth2);
      }      
    }
  }	
  */


  if (isMC==0) {
    // Use measured histogram from Matt & Kurt's file
	   
    // PbPb file:

    TFile *infMatt = new TFile(fileNamePbPb_data);
    TH1F *hMattPbPb = NULL;
    TH1F *hRebinPbPb = NULL;
    TH1F *hTagEffPbPb = NULL;

    if(doBjets){
      hMattPbPb = (TH1F*) infMatt->Get("hRawBData");
      hTagEffPbPb = (TH1F*) infMatt->Get("hBEfficiencyMC");
    }
//    else hMattPbPb = (TH1F*) infMatt->Get(Form("DataJetInEtaBin%.f_%.f",etamin*10, etamax*10));
       else {
        if(TMath::Abs(etamin)==1.)
        //   hMattPbPb = (TH1F*) infMatt->Get("hjtpt");
           hMattPbPb = (TH1F*) infMatt->Get("jetpt_0-100%");
        else
        //  hMattPbPb = (TH1F*) infMatt->Get(Form("jetptEtaBin%.f_%.f", etamin*10, etamax*10));
          hMattPbPb = (TH1F*) infMatt->Get(Form("jetptEtaBin%.f_%.f_Cen0-100%%", etamin*10, etamax*10));
        }
//       hMattPbPb->Scale(1./nMB);
   //  hMattPbPb->Scale(1./(etamax-etamin));    
   //  hMattPbPb->Scale(1./2.);    
  //  else hMattPbPb = (TH1F*) infMatt->Get("hjtpt");
  //  divideBinWidth(hMattPbPb);
           
    // Need to match the binning carefully, please double check whenever you change the binning
     hRebinPbPb = (TH1F*)hMattPbPb->Rebin(nbins_rec, Form("DataJetPtEtaBin%.f_%.f",etamin*10, etamax*10),boundaries_rec);
 /*   for (int i=1;i<=hMattPbPb->GetNbinsX();i++)
      {
     	float binContent = hMattPbPb->GetBinContent(i);  
	float binError = hMattPbPb->GetBinError(i); 

	if(doBjets){
	  float tagEff =hTagEffPbPb->GetBinContent(i);
	  float tagEffErr =     hTagEffPbPb->GetBinError(i);   
	  
	  if(tagEff>0){
	    // careful of the order here!
	    binError=binContent/tagEff *sqrt(tagEffErr*tagEffErr/tagEff/tagEff + binError*binError/binContent/binContent);
	    binContent /= tagEff;
	  }
	  else cout<<" TAGEFF = 0"<<endl;	  
	}

	float binCenter = hMattPbPb->GetBinCenter(i);  
	if(binCenter - hMattPbPb->GetBinWidth(i)/2.  < recoJetPtCut) continue;
	
	int ibin=0;
	
	if(doTrigEffCorr){
	  float trigEff = 0;
	  if(doBjets) trigEff = trigEffbJet[i-1];
	  else  trigEff = trigEffInc[i-1];

	  cout<<" binCenter = "<<binCenter<<" trigEff "<<trigEff<<endl;

	  if(trigEff>0){
	    // careful of the order here!
	    binContent /= trigEff;
	    binError /= trigEff;
	  }
	  else cout<<" TRIGEFF = 0"<<endl;	  
	}


          for(ibin=1;ibin<=uhist[0]->hMeas->GetNbinsX();ibin++){
          float testLowEdge = uhist[0]->hMeas->GetBinLowEdge(ibin);
          float testBinWidth = uhist[0]->hMeas->GetBinWidth(ibin);
          cout <<"ibin =" << ibin << " bin edge =" <<testLowEdge << "  wdith =" <<testBinWidth <<endl ;
        if(binCenter>testLowEdge && binCenter < testLowEdge+testBinWidth) break;
        }   

	uhist[0]->hMeas->SetBinContent(ibin,binContent);  
	uhist[0]->hMeas->SetBinError(ibin,binError);  
      }
*/
       //fixed by Yaxian
       for(int ibin=1;ibin<=uhist[0]->hMeas->GetNbinsX();ibin++){
        uhist[0]->hMeas->SetBinContent(ibin,0);
        uhist[0]->hMeas->SetBinError(ibin,0); 
        float binCenter = uhist[0]->hMeas->GetBinCenter(ibin);
        float testBinWidth = uhist[0]->hMeas->GetBinWidth(ibin);
        float testLowEdge = uhist[0]->hMeas->GetBinLowEdge(ibin);
        float binContent =  hRebinPbPb->GetBinContent(hRebinPbPb->FindBin(binCenter));
         float binError = hRebinPbPb->GetBinError(hRebinPbPb->FindBin(binCenter));   

      //  cout <<"get bins " << testLowEdge<<" content =" <<binContent<<endl ;
        uhist[0]->hMeas->SetBinContent(ibin,binContent);  
        uhist[0]->hMeas->SetBinError(ibin,binError);  
       }
 
   TFile *infMattPP = new TFile(fileNamePP_data);
    TH1F *hMattPP = NULL;
    TH1F *hTagEffPP = NULL;
/*    if(doBjets){
      hMattPP = (TH1F*) infMattPP->Get("hRawBData");
      hTagEffPP = (TH1F*) infMattPP->Get("hBEfficiencyMC");
    }
    else hMattPP = (TH1F*) infMattPP->Get(Form("PYTHIAJetInEtaBin%.f_%.f_Cen0-100%%",etamin*10, etamax*10));
*/
 /*   else {
        if(TMath::Abs(etamin)==1.)
           hMattPP = (TH1F*) infMattPP->Get("hjtpt");
        else 
          hMattPP = (TH1F*) infMattPP->Get(Form("jetptEtaBin%.f_%.f", etamin*10, etamax*10));
        }
 */   //divideBinWidth(hMattPP);	   
  /*    for (int i=1;i<=hMattPP->GetNbinsX();i++)
      {
      	float binContent = hMattPP->GetBinContent(i);  
	float binError = hMattPP->GetBinError(i);  
	
	if(doBjets){
	  float tagEff =hTagEffPP->GetBinContent(i);
	  float tagEffErr =     hTagEffPP->GetBinError(i);   
	  if(tagEff>0){
	    // careful of the order here!
	    binError=binContent/tagEff *sqrt(tagEffErr*tagEffErr/tagEff/tagEff + binError*binError/binContent/binContent);
	    binContent /= tagEff;
	  }
	  else cout<<" TAGEFF = 0"<<endl;	  
	}
	
     	float binCenter = hMattPP->GetBinCenter(i);  
	if(binCenter - hMattPP->GetBinWidth(i)/2.  < recoJetPtCut) continue;

	int ibin=0;

	for(ibin=1;ibin<=uhist[nbins_cent]->hMeas->GetNbinsX();ibin++){
	  float testLowEdge = uhist[nbins_cent]->hMeas->GetBinLowEdge(ibin);
	  float testBinWidth = uhist[nbins_cent]->hMeas->GetBinWidth(ibin);
                   cout <<"PP ibin =" << ibin << " bin edge =" <<testLowEdge << "  wdith =" <<testBinWidth <<endl ;
	  if(binCenter>testLowEdge && binCenter < testLowEdge+testBinWidth) break;

	}
    
        uhist[nbins_cent]->hMeas->SetBinContent(ibin,binContent);
        uhist[nbins_cent]->hMeas->SetBinError(ibin,binError);
      }
*/
   
/*    // pp file: MC sample, not necessary for unfolding, using Gen infor directly, comment out this part moment
     //fixed by Yaxian
       for(int ibin=1;ibin<=uhist[nbins_cent]->hMeas->GetNbinsX();ibin++){
          uhist[nbins_cent]->hMeas->SetBinContent(ibin,0);
          uhist[nbins_cent]->hMeas->SetBinError(ibin,0);
         float binCenter = uhist[nbins_cent]->hMeas->GetBinCenter(ibin);
         float testBinWidth = uhist[nbins_cent]->hMeas->GetBinWidth(ibin);
         float testLowEdge = uhist[nbins_cent]->hMeas->GetBinLowEdge(ibin);
         float binContent =  hMattPP->GetBinContent(hMattPP->FindBin(binCenter));
         float binError = hMattPP->GetBinError(hMattPP->FindBin(binCenter));
        uhist[nbins_cent]->hMeas->SetBinContent(ibin,binContent);
        uhist[nbins_cent]->hMeas->SetBinError(ibin,binError);
       }
*/
  }



   //=========================================Response Matrix========================================================= 

  cout <<"Response Matrix..."<<endl;
	
  TCanvas * cMatrix = new TCanvas("cMatrix","cMatrix",800,400);
  cMatrix->Divide(2,1);
// cMatrix->cd(1);
  TH2F * hResponse[nbins_cent+1];
  for (int i=0;i<=nbins_cent;i++){
    cMatrix->cd(i+1);
    if (!useMatrixFromFile) {
      TF1 *f = new TF1("f","[0]*pow(x+[2],[1])");
      f->SetParameters(1e10,-8.8,40);
      for (int y=1;y<=uhist[i]->hMatrix->GetNbinsY();y++) {
	double sum=0;
	for (int x=1;x<=uhist[i]->hMatrix->GetNbinsX();x++) {
	  if (uhist[i]->hMatrix->GetBinContent(x,y)<=1*uhist[i]->hMatrix->GetBinError(x,y)) {
	    uhist[i]->hMatrix->SetBinContent(x,y,0);
	    uhist[i]->hMatrix->SetBinError(x,y,0);
	  }
	  sum+=uhist[i]->hMatrix->GetBinContent(x,y);
	}
				
	for (int x=1;x<=uhist[i]->hMatrix->GetNbinsX();x++) {	   
	  double ratio = 1;
	  uhist[i]->hMatrix->SetBinContent(x,y,uhist[i]->hMatrix->GetBinContent(x,y)*ratio);
	  uhist[i]->hMatrix->SetBinError(x,y,uhist[i]->hMatrix->GetBinError(x,y)*ratio);
	}
      }
    }
    uhist[i]->hResponse = (TH2F*)uhist[i]->hMatrix->Clone(Form("hResponse_cent%d",i));


    for (int y=1;y<=uhist[i]->hResponse->GetNbinsY();y++) {
      double sum=0;
      for (int x=1;x<=uhist[i]->hResponse->GetNbinsX();x++) {
	if (uhist[i]->hResponse->GetBinContent(x,y)<=1*uhist[i]->hResponse->GetBinError(x,y)) {
	  uhist[i]->hResponse->SetBinContent(x,y,0);
	  uhist[i]->hResponse->SetBinError(x,y,0);
	}
	sum+=uhist[i]->hResponse->GetBinContent(x,y);
      }
			
      for (int x=1;x<=uhist[i]->hResponse->GetNbinsX();x++) {  	
	if (sum==0) continue;
	double ratio =1.;
	if (uhist[i]->hMeas->GetBinContent(y)==0) ratio = 1e-100/sum;
         else ratio = uhist[i]->hMeas->GetBinContent(y)/sum;
     //  ratio = 1./sum ;
    //   if (hProj->GetBinContent(y)==0) ratio = 1e-100/sum;
    //     else ratio = hProj->GetBinContent(y)/sum;
//       uhist[i]->hResponse->SetBinContent(x,y,uhist[i]->hResponse->GetBinContent(x,y)*ratio);
//        uhist[i]->hResponse->SetBinError(x,y,uhist[i]->hResponse->GetBinError(x,y)*ratio);
      }
    }
		
    
    TH1F *hProj = (TH1F*)uhist[i]->hResponse->ProjectionY(Form("hProj_cent%d",i));
    for (int y=1;y<=uhist[i]->hResponse->GetNbinsY();y++) {
      for (int x=1;x<=uhist[i]->hResponse->GetNbinsX();x++) {  	
	double sum=hProj->GetBinContent(y);
//	cout <<y<<" "<<x<<" "<<sum<<endl;
	if (sum==0) continue;
	double ratio ;
	if (uhist[i]->hMeas->GetBinContent(y)==0) ratio = 1e-100/sum;
         else ratio = uhist[i]->hMeas->GetBinContent(y)/sum;
        uhist[i]->hResponse->SetBinContent(x,y,uhist[i]->hResponse->GetBinContent(x,y)*ratio);
	uhist[i]->hResponse->SetBinError(x,y,uhist[i]->hResponse->GetBinError(x,y)*ratio);
      }
    }

    uhist[i]->hResponseNorm = (TH2F*)uhist[i]->hMatrix->Clone(Form("hResponseNorm_cent%d",i));
    for (int x=1;x<=uhist[i]->hResponseNorm->GetNbinsX();x++) {
      double sum=0;
      for (int y=1;y<=uhist[i]->hResponseNorm->GetNbinsY();y++) {
	if (uhist[i]->hResponseNorm->GetBinContent(x,y)<=0*uhist[i]->hResponseNorm->GetBinError(x,y)) {
	  uhist[i]->hResponseNorm->SetBinContent(x,y,0);
	  uhist[i]->hResponseNorm->SetBinError(x,y,0);
	}
	sum+=uhist[i]->hResponseNorm->GetBinContent(x,y);
      }
			
      for (int y=1;y<=uhist[i]->hResponseNorm->GetNbinsY();y++) {  	
	if (sum==0) continue;
	double ratio = 1./sum;
	uhist[i]->hResponseNorm->SetBinContent(x,y,uhist[i]->hResponseNorm->GetBinContent(x,y)*ratio);
	uhist[i]->hResponseNorm->SetBinError(x,y,uhist[i]->hResponseNorm->GetBinError(x,y)*ratio);
      }
    }

    if (!useMatrixFromFile) uhist[i]->hMatrixFit = uhist[i]->hMatrix;
    uhist[i]->hMatrixFit->SetName(Form("hMatrixFit_cent%d",i));
/*  } 

  TCanvas * cMatrix = new TCanvas("cMatrix","Matrix",800,400);
  cMatrix->Divide(2,1);
//  cMatrix->cd(1);
    for (int i=0;i<=nbins_cent;i++){
*/    cMatrix->cd(i+1); 
 //     cMatrix->cd(i+1);		
    uhist[i]->hResponse->SetTitleOffset(1.4, "Y");
    uhist[i]->hResponse->SetTitleOffset(1.2, "X");
  //  if(isMC)uhist[i]->hResponse->SetMinimum(1.e-8);
  //  else
      uhist[i]->hResponse->SetMinimum(1.e-10);
    uhist[i]->hResponse->DrawCopy("colz");
/*    uhist[i]->hMatrix->SetTitleOffset(1.4, "Y");
    uhist[i]->hMatrix->SetTitleOffset(1.2, "X");
    uhist[i]->hMatrix->Draw("colz");
		
   cMatrix->cd(i+1)->Modified();
   cMatrix->cd(i+1);
*/  }

  // cMatrix->SetSelected(cMatrix);
  cMatrix->Update();
	
  cout << "==================================== UNFOLD ===================================" << endl;
	
  //char chmet[100]; 
	
  // ======================= Reconstructed pp and PbPb spectra =========================================================
	
  TH1F *hRecoBW[nbins_cent+1], *hRecoBinByBinBW[nbins_cent+1], *hMeasBW[nbins_cent+1], *hGenBW[nbins_cent+1];
 TH1F *hReproduced[nbins_cent+1]; 

	
  pbpb_Unfo->cd();

  for (int i=0;i<=nbins_cent;i++) {
    // Do Bin-by-bin
    TH1F *hBinByBinCorRaw = (TH1F*)uhist[i]->hResponse->ProjectionY(); 
     hBinByBinCorRaw->Sumw2();
    TH1F *hMCGen           = (TH1F*)uhist[i]->hResponse->ProjectionX(); // gen
    hMCGen->Sumw2();
    hBinByBinCorRaw->Divide(hMCGen);
    TF1 *f = new TF1("f","[0]+[1]*x");
    hBinByBinCorRaw->Fit("f","LL ","",recoJetPtCut,850);
    TH1F* hBinByBinCor = (TH1F*)hBinByBinCorRaw->Clone();//functionHist(f,hBinByBinCorRaw,Form("hBinByBinCor_cent%d",i));
    uhist[i]->hRecoBinByBin = (TH1F*) uhist[i]->hMeas->Clone(Form("hRecoBinByBin_cent%d",i));
    uhist[i]->hRecoBinByBin->Divide(hBinByBinCor);
  
    TH1F *hPrior;//=(TH1F*) functionHist(fPow,uhist[i]->hMeas,Form("hPrior_cent%d",i));
    if(isMC) 
       hPrior=(TH1F*)hMCGen->Clone(Form("hPrior_cent%d", i));
    else {
    //  if(i==nbins_cent) 
         hPrior=(TH1F*)hMCGen->Clone(Form("hPrior_cent%d", i));
   //  else 
   //   hPrior=(TH1F*)uhist[i]->hRecoBinByBin->Clone(Form("hPrior_cent%d", i));
   //    hPrior=(TH1F*)uhist[i]->hMeas->Clone(Form("hPrior_cent%d", i));
     }
     removeZero(hPrior);
		
    // Do unfolding
    //if (isMC) uhist[i]->hMeas = (TH1F*)uhist[i]->hMatrix->ProjectionY()->Clone(Form("hMeas_cent%d",i));
    prior myPrior(uhist[i]->hMatrixFit,hPrior, 0);
  //  prior myPrior(uhist[i]->hMatrixFit,uhist[i]->hMeas,0);
   // myPrior.unfold(uhist[i]->hMeas,1);
    myPrior.unfold(uhist[i]->hMeas,nBayesianIter);
//    hPrior = (TH1F*)uhist[i]->hGen->Clone("hPrior");
//    hPrior = (TH1F*)uhist[i]->hMeas->Clone(Form("hPrior_cent%d",i));
  //  TH1F *hReweighted = (TH1F*)(TH1F*)uhist[i]->hResponse->ProjectionY(Form("hReweighted_cent%d",i));
		
    bayesianUnfold myUnfoldingJECSys(uhist[i]->hMatrixFit,hPrior,0);
    myUnfoldingJECSys.unfold(uhist[i]->hMeasJECSys,nBayesianIter);
    bayesianUnfold myUnfoldingSmearSys(uhist[i]->hMatrixFit,hPrior,0);
    myUnfoldingSmearSys.unfold(uhist[i]->hMeasSmearSys,nBayesianIter);
    bayesianUnfold myUnfolding(uhist[i]->hMatrixFit,myPrior.hPrior,0);
    myUnfolding.unfold(uhist[i]->hMeas,nBayesianIter);
    cout <<"Unfolding bin "<<i<<endl;

    delete hBinByBinCorRaw;
    delete hMCGen;

    // Iteration Systematics
    for (int j=2;j<=10;j++)
      {

	bayesianUnfold myUnfoldingSys(uhist[i]->hMatrixFit,hPrior,0);
	myUnfoldingSys.unfold(uhist[i]->hMeas,j);
	uhist[i]->hRecoIterSys[j]  = (TH1F*) myUnfoldingSys.hPrior->Clone(Form("hRecoRAA_IterSys%d_cent%d",j,i));
      }
    
    
    uhist[i]->hReco         = (TH1F*) uhist[i]->hRecoIterSys[nBayesianIter]->Clone(Form("Unfolded_cent%i",i));
    uhist[i]->hRecoJECSys   = (TH1F*) myUnfoldingJECSys.hPrior->Clone(Form("UnfoldedJeCSys_cent%i",i));
    uhist[i]->hRecoSmearSys   = (TH1F*) myUnfoldingSmearSys.hPrior->Clone(Form("UnfoldedSmearSys_cent%i",i));
    uhist[i]->hRecoBinByBin->SetName(Form("UnfoldedBinByBin_cent%i",i));
    
    if (doToy) {
      TCanvas *cToy = new TCanvas("cToy","toy",600,600);
      cToy->cd();
      int nExp=1000;
      TH1F *hTmp[nbins_truth+1];
      TH1F *hTmp2[nbins_truth+1];
      for (int j=1;j<=nbins_truth;j++) {
	hTmp[j] = new TH1F(Form("hTmp%d",j),"",200,0,10.+uhist[i]->hReco->GetBinContent(j)*2);
	hTmp2[j] = new TH1F(Form("hTmp2%d",j),"",200,0,10.+uhist[i]->hRecoBinByBin->GetBinContent(j)*2);
      }
      for (int exp =0; exp<nExp; exp++) {
	TH1F *hToy = (TH1F*)uhist[i]->hMeas->Clone();   
	TH2F *hMatrixToy = (TH2F*)uhist[i]->hMatrixFit->Clone();
	hToy->SetName("hToy");
	if (exp%100==0) cout <<"Pseudo-experiment "<<exp<<endl;
	for (int j=1;j<=hToy->GetNbinsX();j++) {
	  double value = gRandom->Poisson(uhist[i]->hMeas->GetBinContent(j));
	  hToy->SetBinContent(j,value);
	}
				
	for (int j=1;j<=hMatrixToy->GetNbinsX();j++) {
	  for (int k=1;k<=hMatrixToy->GetNbinsY();k++) {
	    double value = gRandom->Gaus(uhist[i]->hMatrixFit->GetBinContent(j,k),uhist[i]->hMatrixFit->GetBinError(j,k));
	    hMatrixToy->SetBinContent(j,k,value);
	  }
	}

	prior myPriorToy(hMatrixToy,hToy,0.0);
	myPriorToy.unfold(hToy,1);
	bayesianUnfold myUnfoldingToy(hMatrixToy,myPriorToy.hPrior,0.0);
	myUnfoldingToy.unfold(hToy,nBayesianIter);
	TH1F *hRecoTmp = (TH1F*) myUnfoldingToy.hPrior->Clone();
				
	for (int j=1;j<=hRecoTmp->GetNbinsX();j++) {
	  hTmp[j]->Fill(hRecoTmp->GetBinContent(j));
	}
	delete hToy;
	delete hRecoTmp;
	delete hMatrixToy;
      }
      TF1 *fGaus = new TF1("fGaus","[0]*TMath::Gaus(x,[1],[2])");
      for (int j=1;j<=nbins_truth;j++)
	{

	  f->SetParameters(hTmp[j]->GetMaximum(),hTmp[j]->GetMean(),hTmp[j]->GetRMS());
				
	  if (hTmp[j]->GetMean()>0) {
	    hTmp[j]->Fit(fGaus,"LL Q ");
	    hTmp[j]->Fit(fGaus,"LL Q ");
	    uhist[i]->hReco->SetBinError(j,f->GetParameter(2));
	  }	       
	  f->SetParameters(hTmp2[j]->GetMaximum(),hTmp2[j]->GetMean(),hTmp2[j]->GetRMS());
	  if (hTmp2[j]->GetMean()>0) {
	    hTmp2[j]->Fit(fGaus,"LL Q ");
	    hTmp2[j]->Fit(fGaus,"LL Q ");
	    uhist[i]->hRecoBinByBin->SetBinError(j,f->GetParameter(2));
	  }	       
	  delete hTmp[j];
	  delete hTmp2[j];
	}
    }  //finish the doToy part

 /*  //do binwidth normalization
   divideBinWidth(uhist[i]->hMeas);
   divideBinWidth(uhist[i]->hReco);
   divideBinWidth(uhist[i]->hRecoBinByBin);
   divideBinWidth(uhist[i]->hGen);
   */ 
    uhist[i]->hMeas->SetMarkerStyle(20);
    uhist[i]->hMeas->SetMarkerColor(2);
    uhist[i]->hReco->SetMarkerStyle(25);
    uhist[i]->hReco->SetName(Form("hReco_cent%d",i));
    
    uhist[i]->hReco->SetTitle("Baysian Unfolded");
    uhist[i]->hRecoBinByBin->SetTitle("Bin-by-bin Unfolded");
    uhist[i]->hReco->SetXTitle("p_{T} (GeV/c)");    
    uhist[i]->hReco->SetYTitle("Counts");    
    uhist[i]->hReco->GetXaxis()->SetNdivisions(505);
    uhist[i]->hReco ->GetYaxis()->SetTitleOffset(1.3);
    uhist[i]->hReco->GetXaxis()->SetTitleOffset(1.2);
    //uhist[i]->hReco->Draw("");    
    uhist[i]->hReco->SetAxisRange(recoJetPtCut,600,"X");
	    uhist[i]->hReco->SetAxisRange(recoJetPtCut,600, "X");
	    uhist[i]->hReco->GetXaxis()->SetRangeUser(recoJetPtCut,600);
	    uhist[i]->hReco->GetXaxis()->SetNdivisions(505);
	   // uhist[i]->hReco->DrawCopy("");   
	     
	    uhist[i]->hGen->SetLineWidth(2);
	    uhist[i]->hGen->SetLineColor(2);
	   // if(isMC)uhist[i]->hGen->Draw("hist same");
	   // uhist[i]->hReco->Draw("same");    
	    uhist[i]->hRecoBinByBin->SetMarkerStyle(28);
	    uhist[i]->hRecoBinByBin->SetMarkerColor(4);
	    uhist[i]->hRecoBinByBin->SetLineColor(4);
	  //  uhist[i]->hRecoBinByBin->DrawCopy("same");    

	    hReproduced[i] = (TH1F*)myUnfolding.hReproduced->Clone(Form("hReproduced_cent%d",i));
	    hReproduced[i]->SetMarkerColor(1);
	    hReproduced[i]->SetMarkerStyle(24);
	    //uhist[i]->hMeas->Draw("same");    

	    hRecoBW[i] = (TH1F*)uhist[i]->hReco->Clone(Form("hReco%d",i));
	    hRecoBinByBinBW[i] = (TH1F*)uhist[i]->hRecoBinByBin->Clone(Form("hRecoBinByBin%d",i));
	    hMeasBW[i] = (TH1F*)uhist[i]->hMeas->Clone(Form("hMeas%d",i));
	  //  if(isMC)
               hGenBW[i] = (TH1F*)uhist[i]->hGen->Clone(Form("hGen%d",i));

	    divideBinWidth(hRecoBW[i]);    
	  //  if(isMC)
                 divideBinWidth(hGenBW[i]); 
           //   if(i==0) 
           //      hGenBW[i]->Scale(1./(etamax-etamin));  
	    divideBinWidth(hRecoBinByBinBW[i]);    
	    divideBinWidth(hMeasBW[i]);    
	    divideBinWidth(hReproduced[i]);
	    hRecoBW[i]->SetTitle("Baysian Unfolded");
	    hRecoBinByBinBW[i]->SetTitle("Bin-by-bin Unfolded");
   }

  TCanvas * cPbPb = new TCanvas("cPbPb","Comparison",1200,600);
  cPbPb->Divide(2,1);
  cPbPb->cd(1);
    for (int i=0;i<=nbins_cent;i++) {
           cPbPb->cd(i+1);
           cPbPb->cd(i+1)->SetLogy();   
	    hMeasBW[i]->SetAxisRange(recoJetPtCut,600, "X");
	    hMeasBW[i]->GetXaxis()->SetRangeUser(recoJetPtCut,600);
	    hMeasBW[i]->GetYaxis()->SetRangeUser(1.e-10, 1.e-3);
            hMeasBW[i]->GetXaxis()->SetNdivisions(505);
	    hMeasBW[i]->DrawCopy("");
	  //  hRecoBW[i]->DrawCopy();
	  //  if(isMC)
               hGenBW[i]->DrawCopy("hist,same");
	    hRecoBinByBinBW[i]->DrawCopy("same");
	    hRecoBW[i]->DrawCopy("same");
	    


	    TLegend *leg = new TLegend(0.5,0.5,0.9,0.9);
	    leg->SetBorderSize(0);
	    leg->SetFillStyle(0);
	/*    if(i==0)leg->AddEntry(uhist[i]->hMeas,"PPb","");
	    if(i==1)leg->AddEntry(uhist[i]->hMeas,"PP",""); 
	   leg->AddEntry(uhist[i]->hMeas,"Measured","pl");
	    leg->AddEntry(uhist[i]->hReco,"Bayesian unfolded","pl");
	    leg->AddEntry(uhist[i]->hRecoBinByBin,"Bin-by-bin unfolded","pl");
	    if(isMC)leg->AddEntry(uhist[i]->hGen,"Generator level truth","l");
	 */   if(i==0)leg->AddEntry(hRecoBW[i],"PPb","");
	    if(i==1)leg->AddEntry(hRecoBW[i],"PP","");
	   leg->AddEntry(hMeasBW[i],"Measured","pl");
	    leg->AddEntry(hRecoBW[i],"Bayesian unfolded","pl");
	    leg->AddEntry(hRecoBinByBinBW[i],"Bin-by-bin unfolded","pl");
	   // if(isMC)
               leg->AddEntry(hGenBW[i],"Generator level truth","l");
	    leg->Draw("same");
       }
     	    cPbPb->Update();

    TCanvas * cPbPbMeas = new TCanvas("cPbPbMeas","Measurement",1200,600);
    cPbPbMeas->Divide(2,1);
    cPbPbMeas->cd(1);
         for (int i=0;i<=nbins_cent;i++) {
            cPbPbMeas->cd(i+1);
	    cPbPbMeas->cd(i+1)->SetLogy();   
	  //  uhist[i]->hMeas->SetAxisRange(30,300,"X");
	    hMeasBW[i]->SetAxisRange(recoJetPtCut,600,"X");
	    hMeasBW[i]->GetXaxis()->SetRangeUser(recoJetPtCut,600);
	    hMeasBW[i]->GetXaxis()->SetNdivisions(505);
            // divideBinWidth(hMeasBW[i]);
	  //   divideBinWidth(hReproduced);
	  //  uhist[i]->hMeas->Draw();
	    hMeasBW[i]->DrawCopy();
	    hReproduced[i]->DrawCopy("same");

	    TLegend *leg2 = new TLegend(0.5,0.5,0.85,0.9);
	    leg2->SetBorderSize(0);
	    leg2->SetFillStyle(0);
	  // if(isMC) leg2->AddEntry(hReproduced,"PYTHIA+HIJING","");
	  //  else leg2->AddEntry(hReproduced,"pPb","");
	    if(i==0)leg2->AddEntry(hMeasBW[i],"PPb","");
	    if(i==1)leg2->AddEntry(hMeasBW[i],"PP","");
	    leg2->AddEntry(hMeasBW[i],"Measured","pl");
	    leg2->AddEntry(hReproduced[i],"Reproduced","pl");

	    leg2->Draw("same");
	  }	     
	    cPbPbMeas->Update();

	 // ======================= Unfolding closure in MC =========================================================
	 TCanvas * cRatio ;
	    TH1F * hReco[nbins_cent+1];
	    TH1F * hRecoBinByBin[nbins_cent+1];
	    TH1F * hMeas[nbins_cent+1];
	    TH1F * hGenScale[nbins_cent+1];
	    TLegend *leg[nbins_cent+1];
	  TLine *line = new TLine(recoJetPtCut,1,600,1);
		line->SetLineStyle(2);
		line->SetLineWidth(2);

	    for (int i=0;i<=nbins_cent;i++) {
		hReco[i]          = (TH1F*)uhist[i]->hReco->Clone(Form("hRecoScaled_Cent%d", i));
		hRecoBinByBin[i]          = (TH1F*)uhist[i]->hRecoBinByBin->Clone(Form("hRecoBinByBinScaled_Cent%d", i));
		hMeas[i]          = (TH1F*)uhist[i]->hMeas->Clone(Form("hMeasScaled_Cent%d", i));
	//	if(isMC) 
                    hGenScale[i]          = (TH1F*)uhist[i]->hGen->Clone(Form("hGenScaled_Cent%d", i));
                divideBinWidth(hReco[i]);
               // if(isMC)
                    divideBinWidth(hGenScale[i]);
                divideBinWidth(hRecoBinByBin[i]);
                divideBinWidth(hMeas[i]);
	    }
	  if(isMC){
	  cRatio = new TCanvas("cRatio","Ratio",1200,600);
	    cRatio->Divide(2,1);
	      for (int i=0;i<=nbins_cent;i++) {
            /*    divideBinWidth(hReco[i]);
                divideBinWidth(hGen[i]);
                divideBinWidth(hRecoBinByBin[i]);
                divideBinWidth(hMeas[i]);
*/
		  hMeas[i]->Divide(hGenScale[i]);
		  hRecoBinByBin[i]->Divide(hGenScale[i]);
		  hReco[i]->Divide(hGenScale[i]);
		  cRatio->cd(i+1);

			//hRecoPP->SetAxisRange(90,300,"X");
			hReco[i]->SetAxisRange(0.5,1.5,"Y");
		//	hReco[i]->SetAxisRange(20.,600.,"X");
			hReco[i]->GetXaxis()->SetRangeUser(recoJetPtCut,600.);
			hReco[i]->GetXaxis()->SetNdivisions(505);
			hReco[i]->SetMarkerStyle(25);
			hReco[i] ->SetLineColor(1);
			hReco[i] ->SetMarkerColor(1);
			hMeas[i]->SetMarkerStyle(20);
			hMeas[i]->SetLineColor(2);
			hMeas[i]->SetMarkerColor(2);
			hRecoBinByBin[i]->SetMarkerStyle(28);
			hRecoBinByBin[i]->SetLineColor(4);
			hRecoBinByBin[i]->SetMarkerColor(4);

			makeHistTitle(hReco[i],"",Form("Jet p_{T} (GeV/c)"),Form("Reco / Truth"));
			hReco[i]->GetYaxis()->SetTitleOffset(1.3);
	                hReco[i]->GetXaxis()->SetNdivisions(505);
         		hReco[i]->GetXaxis()->SetTitleOffset(1.2);
			hReco[i]->DrawCopy("");
			hRecoBinByBin[i]->Draw("same");
			hMeas[i]->DrawCopy("same");
			line->Draw("same");
			leg[i] = myLegend(0.52,0.65,0.85,0.9);
			 if(i==0)leg[i]->AddEntry(hReco[i],"PYTHIA+HIJING","");
			 if(i==1)leg[i]->AddEntry(hReco[i],"PYTHIA","");
			leg[i]->AddEntry(hReco[i],"Bayesian","pl");
			leg[i]->AddEntry(hRecoBinByBin[i],"Bin-by-bin","pl");
			leg[i]->AddEntry(hMeas[i],"no unfolding","pl");
			leg[i]->Draw("same");
			putCMSPrel(0.2,0.83,0.04);
			drawText("Anti-k_{T} Particle Flow Jets   R = 0.3",0.2,0.23,20);
			drawText("CMS Simulation",0.5,0.4,20);
			drawText(Form("%.1f< #eta_{CM} <%.1f ",etamin, etamax),0.5,0.31,20);
	      }
	  }
	  else {
	/*      hMeas[nbins_cent]->Scale(1./5.33e9);
	      hRecoBinByBin[nbins_cent]->Scale(1./5.33e9);
	      hReco[nbins_cent]->Scale(1./5.33e9);
	      float nMB = 1.09079385600000000e+09;
	    
               divideBinWidth(hReco[nbins_cent]);
               divideBinWidth(hRecoBinByBin[nbins_cent]);
               divideBinWidth(hMeas[nbins_cent]);
         */     hMeas[nbins_cent]->Scale(1./ppsigma);
	      hRecoBinByBin[nbins_cent]->Scale(1./ppsigma);
	      hReco[nbins_cent]->Scale(1./ppsigma);
              hGenScale[nbins_cent]->Scale(1./ppsigma); 

       /*       hMeas[nbins_cent]->Scale(2.) ; //normalize to per eta when using the MC file directly
              hRecoBinByBin[nbins_cent]->Scale(2.);
              hReco[nbins_cent]->Scale(2.);
      */
	     cRatio = new TCanvas("cRatio","Ratio",600,600);
	      cRatio->cd(1);

	    
	      for (int i=0;i<nbins_cent;i++) {
	    /*      hMeas[i]            ->Scale(1./TAA[i]/nMB);
		  hRecoBinByBin[i]            ->Scale(1./TAA[i]/nMB);
		  hReco[i]            ->Scale(1./TAA[i]/nMB);

                  hMeas[i]            ->Scale(1./nMB);
                  hRecoBinByBin[i]            ->Scale(1./nMB);
                  hReco[i]            ->Scale(1./nMB);            
*/		  hMeas[i]            ->Scale(1./Ncoll);
		  hRecoBinByBin[i]            ->Scale(1./Ncoll);
		  hReco[i]            ->Scale(1./Ncoll);
/* 
		  hMeas[i]->Divide(hMeas[nbins_cent]);
		  hRecoBinByBin[i]->Divide(hRecoBinByBin[nbins_cent]);
		  hReco[i]->Divide(hReco[nbins_cent]);
*/
                  hMeas[i]->Divide(hGenScale[nbins_cent]);
                  hRecoBinByBin[i]->Divide(hGenScale[nbins_cent]);
                  hReco[i]->Divide(hGenScale[nbins_cent]);

		  hMeas[i]->SetAxisRange(0,3.,"Y");
		  hMeas[i]->SetAxisRange(recoJetPtCut,600,"X");
		  hMeas[i]->GetXaxis()->SetRangeUser(recoJetPtCut,600);
		  hMeas[i]->GetXaxis()->SetNdivisions(505);
		  hReco[i]->SetMarkerStyle(25);
		  hReco[i] ->SetLineColor(1);
		  hReco[i] ->SetMarkerColor(1);
		  hMeas[i]->SetMarkerStyle(20);
		  hMeas[i]->SetLineColor(2);
		  hMeas[i]->SetMarkerColor(2);
		  hRecoBinByBin[i]->SetMarkerStyle(28);
		  hRecoBinByBin[i]->SetLineColor(4);
		  hRecoBinByBin[i]->SetMarkerColor(4);
		 // if(i==0){
		  makeHistTitle(hMeas[i],"",Form("Jet p_{T} (GeV/c)"),Form("Jet R_{pA}"));
		  hMeas[i]->GetYaxis()->SetTitleOffset(1.2);
		  hMeas[i]->GetXaxis()->SetTitleOffset(1.2);
		  hMeas[i]->Draw("");
		      leg[i] = myLegend(0.62,0.75,0.85,0.9);
		      if(i==0)leg[i]->AddEntry(hReco[i],"PPb","");
		      if(i==1)leg[i]->AddEntry(hReco[i],"PP",""); 
		      leg[i]->AddEntry(hReco[i],"Bayesian","pl");
		      leg[i]->AddEntry(hRecoBinByBin[i],"Bin-by-bin","pl");
		      leg[i]->AddEntry(hMeas[i],"no unfolding","pl");
		      leg[i]->Draw("same");
		//  hReco[i]->Draw("same");
		  hRecoBinByBin[i]->Draw("same");
		  hReco[i]->Draw("same");
	      //    }
	       }
		 line->Draw("same");

		  putCMSPrel(0.2,0.83,0.04);
		 // drawText(Form("#intL dt = %.f #mub^{-1}",150.),0.2,0.78,20);
		  drawText("Anti-k_{T} PF Jets R = 0.3",0.2,0.33,20);
		  drawText(Form("%.1f< #eta_{CM} <%.1f ",etamin, etamax),0.2,0.27,20);

	  }
	    cRatio->Update();
	 
	  pbpb_Unfo->Write();

	  SysData systematics;
	  //TLine *line = new TLine(60,1,250,1);

	  // Iteration systematics
	  TCanvas *cIterSys = new TCanvas("cIterSys","cIterSys",1200,600);
	  cIterSys->Divide(2,1);
	  cIterSys->cd(2);
	  TH1F *hRecoIterSysPP[100];
	  TH1F *hRebinPP_tmp         = rebin(uhist[nbins_cent]->hReco, (char*)"hRebinPP_tmp");
	  TLegend *legBayesianIterPP = myLegend(0.4,0.7,0.9,0.9);
	  legBayesianIterPP->AddEntry("","PP","");
		 
	  for (int j=2;j<7;j++) {
	    hRecoIterSysPP[j] = rebin(uhist[nbins_cent]->hRecoIterSys[j],Form("hRecoIterSysPP_IterSys%d",j));
	    hRecoIterSysPP[j]->SetLineColor(colorCode[j-2]);
	    hRecoIterSysPP[j]->SetMarkerColor(colorCode[j-2]);
	    hRecoIterSysPP[j]->SetMarkerStyle(20);
	    hRecoIterSysPP[j]->Divide(hRebinPP_tmp);
	    if (j==2){
	      makeHistTitle(hRecoIterSysPP[j],(char*)"",(char*)"Jet p_{T} (GeV/c)",(char*)"Ratio (Unfolded / Nominal)");
	      hRecoIterSysPP[j]->SetTitleOffset(1.3,"Y");
	      hRecoIterSysPP[j]->SetTitleOffset(1.2,"X");
	     // if(isMC)
               hRecoIterSysPP[j]->SetAxisRange(0.5,1.5,"Y");
	    //  else hRecoIterSysPP[j]->SetAxisRange(0.,2.,"Y");
	     // hRecoIterSysPP[j]->SetAxisRange(20,500,"X");
	      hRecoIterSysPP[j]->GetXaxis()->SetRangeUser(recoJetPtCut,600);
	      hRecoIterSysPP[j]->GetXaxis()->SetNdivisions(505);
	      hRecoIterSysPP[j]->DrawCopy(); 
	    } else {
	      hRecoIterSysPP[j]->DrawCopy("same");
	    }
		 
	    checkMaximumSys(systematics.hSysIter[nbins_cent],hRecoIterSysPP[j],0,1.);
	    legBayesianIterPP->AddEntry(hRecoIterSysPP[j],Form("Iteration %d",j),"pl");     
	  }

	  legBayesianIterPP->Draw("same");
	  line->Draw("same");
	  drawEnvelope(systematics.hSysIter[nbins_cent],(char*)"hist same");


	  cIterSys->cd(1);
	  TH1F *hRecoIterSysPbPb[100];
	  TH1F *hRebinPbPb_tmp         = rebin(uhist[0]->hReco, (char*)"hRebinPbPb_tmp");
	  TLegend *legBayesianIterPbPb = myLegend(0.4,0.7,0.9,0.9);
	  legBayesianIterPbPb->AddEntry("","PPb","");
	  for (int j=2;j<7;j++) {
	    hRecoIterSysPbPb[j] = rebin(uhist[0]->hRecoIterSys[j],Form("hRecoIterSysPbPb_IterSys%d",j));
	    hRecoIterSysPbPb[j]->SetLineColor(colorCode[j-2]);
	    hRecoIterSysPbPb[j]->SetMarkerColor(colorCode[j-2]);
	    hRecoIterSysPbPb[j]->SetMarkerStyle(20);
	    hRecoIterSysPbPb[j]->Divide(hRebinPbPb_tmp);
	    if (j==2){
	      makeHistTitle(hRecoIterSysPbPb[j],(char*)"",(char*)"Jet p_{T} (GeV/c)",(char*)"Ratio (Unfolded / Nominal)");
	      hRecoIterSysPbPb[j]->SetTitleOffset(1.3,"Y");
	      hRecoIterSysPbPb[j]->SetTitleOffset(1.2,"X");
	    //  if(isMC)
                 hRecoIterSysPbPb[j]->SetAxisRange(0.5,1.5,"Y");
	    //  else hRecoIterSysPbPb[j]->SetAxisRange(0.,2.,"Y");
	 //     hRecoIterSysPbPb[j]->SetAxisRange(20,500,"X");
	      hRecoIterSysPbPb[j]->GetXaxis()->SetRangeUser(recoJetPtCut,600);
	      hRecoIterSysPbPb[j]->GetXaxis()->SetNdivisions(505);
	      hRecoIterSysPbPb[j]->DrawCopy(); 
	    } else {
	      hRecoIterSysPbPb[j]->DrawCopy("same");
	    }
		 
	    checkMaximumSys(systematics.hSysIter[0],hRecoIterSysPbPb[j],0,1.);
	    legBayesianIterPbPb->AddEntry(hRecoIterSysPbPb[j],Form("Iteration %d",j),"pl");     
	  }
	  legBayesianIterPbPb->Draw("same");
	  line->Draw("same");
	  drawEnvelope(systematics.hSysIter[0],(char*)"hist same");

	 cIterSys->Update();

	  TString data ;
	  if(isMC) data="MC";
	  else data="Data";
	 TString anaType ;
	  if(doBjets) anaType="Bjet";
	  else anaType="Inclusive";

	  if(SavePlot){
	    cMatrix->SaveAs(Form("plots/%s%s%sResponseMatrixEtaBin%.f_%.f.gif", data.Data(), anaType.Data(), algoName[algo], etamin*10,etamax*10));
	    cPbPb->SaveAs(Form("plots/%s%s%sJetSpectraEtaBin%.f_%.f.gif",  data.Data(), anaType.Data(), algoName[algo],etamin*10,etamax*10));
	    cRatio->SaveAs(Form("plots/%s%s%sJetRatioEtaBin%.f_%.f.gif",  data.Data(), anaType.Data(), algoName[algo],etamin*10,etamax*10));
	    cIterSys->SaveAs(Form("plots/%s%s%sIterationSysEtaBin%.f_%.f.gif",  data.Data(), anaType.Data(), algoName[algo],etamin*10,etamax*10));
	   }

	}





