#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <TVector3.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH2F.h>
#include <TFile.h>
#include <TCanvas.h>
#include "TPaveStats.h"
#include <TLorentzVector.h>
#include <TTree.h>
#include <TMath.h>
#include <TLegend.h>
#include <TStyle.h>
#include <THStack.h>
#include <TLatex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <stdio.h>
#include <TString.h>

#include "../binning.h"
#include "../common.h"

#include "../ptBinning.h"
#include "../alphaBinning.h"
//#include "../fitTools.h"

#include "tdrstyle_mod14.C"

void drawDataMcComparisonWithoutRatio(const string& canvasName, THStack *hsMc, TH1 *hMcG, TH1 *hMcQcd, TH1 *hMcTot, TH1 *hData, const string& Xtitle, const string& path, int inLogScale = 1, int i=0, string legendPos = "r"){
	TCanvas *cCanvas = new TCanvas(canvasName.c_str(),canvasName.c_str());
	cCanvas->cd();
	gStyle->SetOptStat(0);
	if(i==1){
		hData->Draw("pe1");
		hData->SetMarkerColor(1);
		hData->SetLineColor(1);
		hData->SetMarkerStyle(20);
		//hData->SetMarkerSize(0.8);
		hsMc->Draw("histsame");
		hsMc->GetXaxis()->SetTitle(Xtitle.c_str());
		hMcTot->SetFillColor(kBlack);
		hMcTot->SetFillStyle(3001);
		hMcTot->DrawCopy("e2same");
		hMcTot->SetFillStyle(1001);
		hMcTot->SetFillColor(47);
		hMcTot->SetLineColor(1);
		gPad-> SetLogy(inLogScale);	
		gPad->RedrawAxis();
	}
	else {
		hsMc->Draw("hist");
		hsMc->GetXaxis()->SetTitle(Xtitle.c_str());
		hData->Draw("pe1same");
		hData->SetMarkerColor(1);
		hData->SetLineColor(1);
		hData->SetMarkerStyle(20);
		//hData->SetMarkerSize(0.8);
		hMcTot->SetFillColor(kBlack);
		hMcTot->SetFillStyle(3001);
		hMcTot->DrawCopy("e2same");
		hMcTot->SetFillStyle(1001);
		hMcTot->SetFillColor(47);
		hMcTot->SetLineColor(1);
		gPad-> SetLogy(inLogScale);
		gPad->RedrawAxis();
	}
	TLegend *myLegend;
	if(legendPos == "l") {
		myLegend = new TLegend(0.13,0.69,0.33,0.88);
	}
	else {
		myLegend = new TLegend(0.74,0.70,0.88,0.88);
	}
	myLegend->SetBorderSize(0);
	myLegend->AddEntry(hMcQcd,"MC QCD","f");
	myLegend->AddEntry(hMcG,"MC #gamma + jet","f");
	myLegend->AddEntry(hData,"Data 2012","p");
	myLegend->SetFillColor(0);
	myLegend->Draw("SAME");
// 	string myPath = aPath + "_MC_woLogScale.pdf"
	cCanvas->SaveAs(path.c_str());
}


void drawDataMcComparison_Flavour(THStack *hsMc, vector<TH1F*> vFlavourHisto, TH1 *hData, const string& path, string sXTitle, string sTitle, int withNoMatched = 1, int inLogScale = 1, int withRatio = 1){
//*******************************************************************************************************
// description: function which draws the response distrition for each 2DTagging zone, 
//              both for data and MC (histogram containing the different flavour responses)
// arguments: - hsMc: THStack containning the MC responses for each flavour
//            - vFlavourHisto: vector containing each flavour response histogram
//            - hData: data response distribution
//            - path: where you want to save the plots
//            - sXTitle: title fo X axis
//            - sTitle: title of the plot
//            - withNoMatched: 1 if you want to include the no matched jets, 0 if you don't
//            - inLogScale: 1 if you want to have log scale on Y axis, 0 if you don't
//*******************************************************************************************************
	TCanvas *cCanvas = new TCanvas("cCanvas","cCanvas");
  cCanvas->cd();
	cout<<"cCanvas->cd()"<<endl;
	gStyle->SetOptStat(0);

  TH1F *hMcTot = NULL;
	if(withNoMatched == 1) {
		hMcTot = (TH1F*)vFlavourHisto[getFlavourNumber()-1]->Clone();
		cout<<"hMcTot created"<<endl;
	}
	else if(withNoMatched == 0 ) {
	hMcTot = (TH1F*)vFlavourHisto[0]->Clone();
		for(int i = 1; i<getFlavourNumber()-2; i++) {
			hMcTot->Add(vFlavourHisto[i]);
			cout<<"hMcTot created"<<endl;
		}
	}

  if(withRatio == 1) {
    // Data / MC comparison

    TPad* pad_hi = new TPad("pad_hi", "", 0., 0.33, 0.99, 0.99);
    pad_hi->Draw();
    //pad_hi->SetLogx();
    pad_hi->SetLeftMargin(0.12);
    pad_hi->SetBottomMargin(0.015);

    // Data / MC ratio
    TPad* pad_lo = new TPad("pad_lo", "", 0., 0., 0.99, 0.33);
    pad_lo->Draw();
    //pad_lo->SetLogx();
    pad_lo->SetLeftMargin(0.12);
    pad_lo->SetTopMargin(1.);
    pad_lo->SetBottomMargin(0.3);

    pad_lo->cd();
    TF1* ratioFit = new TF1("ratioFit", "[0]", hData->GetXaxis()->GetXmin(), hData->GetXaxis()->GetXmax());
    ratioFit->SetParameter(0, 0.);
    ratioFit->SetLineColor(46);
    ratioFit->SetLineWidth(2);
    TH1* h = (TH1F*)hData->Clone();
    h->Divide(hMcTot);
    h1_style_lo(h);
    h->SetStats(1);
    h->SetTitle("Data / MC");
    h->SetXTitle(sXTitle.c_str());
    h->SetYTitle("Data / MC");
    h->SetMarkerSize(1.0);
    h->SetMarkerColor(1);
    h->SetLineColor(1);
    h->SetMarkerStyle(20);
    h->SetMaximum(2);
    h->SetMinimum(0);


    //for(int i=1;i<h->GetEntries();i++){
    //h->SetBinError(i,sqrt(pow(hData->GetBinError(i)/hMc->GetBinContent(i),2)+pow(hMc->GetBinError(i)*hData->GetBinContent(i)/(pow(hMc->GetBinContent(i),2)),2)));
    //}
    h->Fit(ratioFit, "RQ");
    h->GetYaxis()->SetRangeUser(-1,3);
    double fitValue = ratioFit->GetParameter(0);
    double fitError = ratioFit->GetParError(0);

    TPaveText* fitlabel = new TPaveText(0.55, 0.77, 0.88, 0.83, "brNDC");
    fitlabel->SetTextSize(0.08);
    fitlabel->SetFillColor(0);
    fitlabel->SetTextFont(43);
    TString fitLabelText = TString::Format("Fit: %.3f #pm %.3f", fitValue, fitError);
    fitlabel->AddText(fitLabelText);
    fitlabel->Draw("same");

    gStyle->SetOptFit(0);
    h->Draw("PE1same");
    gPad->RedrawAxis();

    pad_hi->cd();
  }



  hsMc->Draw("hist"); 
  h1_style_hi(hData);
  hs_style_hi(hsMc);
  hsMc->SetMinimum(0.1);

  if(inLogScale == 0) {
    Double_t ymax = 0.;
    Double_t ymax_mc = hMcTot->GetMaximum();
    Double_t ymax_data = hData->GetMaximum();
    if(ymax_mc>=ymax_data) {
      ymax = ymax_mc;
    }
    else ymax = ymax_data;
    hsMc->SetMaximum(ymax*1.1);
    //hData->GetYaxis()->SetRangeUser(0.,ymax*1.3);
    gPad->RedrawAxis();
  }
  else if(inLogScale == 1) {
    Double_t ymax = 0.;
    Double_t ymax_mc = hMcTot->GetMaximum();
    Double_t ymax_data = hData->GetMaximum();
    if(ymax_mc>=ymax_data) {
      ymax = ymax_mc;
    }
    else ymax = ymax_data;
    hsMc->SetMaximum(ymax*100.);
    //hData->GetYaxis()->SetRangeUser(0.,ymax*1.3);
    gPad->RedrawAxis();
  }

  hsMc->GetXaxis()->SetTitle(sXTitle.c_str());
  hsMc->SetTitle(sTitle.c_str());
  //cout<<"hsMc style"<<endl;
  hData->Draw("pe1same");
  hData->SetMarkerColor(1);
  hData->SetLineColor(1);
  hData->SetMarkerStyle(20);
  //cout<<"hData style"<<endl;
  hMcTot->SetFillColor(kBlack);
  hMcTot->SetFillStyle(3001);
  hMcTot->DrawCopy("e2same");
  hMcTot->SetFillStyle(1001);
  hMcTot->SetFillColor(47);
  hMcTot->SetLineColor(1);
  //cout<<"hMcTot style"<<endl;
  //hMcTot->SetMarkerSize(0.8);
  //cCanvas->SetLogy();
  gPad-> SetLogy(inLogScale);	
  gPad->RedrawAxis();
  TLegend *aLegend = new TLegend(0.69,0.65,0.89,0.89);
  aLegend->SetBorderSize(0);
  if(withNoMatched == 1) {
    for(int i=0; i<getFlavourNumber()-1; i++) {
      aLegend->AddEntry(vFlavourHisto[i],getFlavourBinName(i).c_str(),"f");
    }
  }
  else if(withNoMatched == 0 ) {
    for(int i=0; i<getFlavourNumber()-2; i++) {
      aLegend->AddEntry(vFlavourHisto[i],getFlavourBinName(i).c_str(),"f");
    }
  }
  aLegend->AddEntry(hData,"Data 2012","p");
  aLegend->SetFillStyle(0);
  aLegend->Draw("SAME");
  cCanvas->SaveAs(path.c_str());
  hMcTot->Delete();
}


void drawRmpfPerZone(TFile* afMC, TH1* hData, vector<int> aHistoColor, string aZoneName, string aPtName, int zoneBin, int isForLumi = 1, int withNoMatched = 1, bool isMCTot = false) {
//*******************************************************************************************************
// description: function which retrieves the informations from the rootfiles to draw the response
//              distrition for each 2DTagging zone, 
//              both for data and MC (histogram containing the different flavour responses)/
//              Function to be called for each pt bin. 
//              Different routines for :
//                  - lumi/shape distributions
//									- with/without no matched jets
//                  - with/without Y log scale
// arguments: - afMC: MC rootfile containning informations we want to extract
//            - hData: data rootfile containning informations we want to extract
//            - aHistoColor: vector containing the flavour colour code
//            - aZoneName: name of the 2DTagging zone, usefull to retrieve the histo from rootfiles
//            - aPtName: name of the pt 2DTagging bin
//            - aCanvas: canvas where plots will be drawn
//            - zoneBin: number of the corresponding 2DTagging zone bin
//            - aPath: where you want to save the plot
//            - isForLumi: 1 if you want MC distributions normalized to the luminosity, 
//                         0 if you want MC distributions normalized to the number of data events (shape)
//            - withNoMatched: 1 if you want to include the no matched jets, 0 if you don't
//            - inLogScale: 1 if you want to have log scale on Y axis, 0 if you don't
//*************************************************************************************************
	vector<TH1F*> myFlavourVector;
	THStack* hsMC = new THStack();
	string aXTitle;
	aXTitle = "R_{" + aZoneName + "}^{MPF}";
	string aTitle;
	aTitle = aXTitle + " distribution";

	// first retrieve the histogramms from the rootfiles by reconstructing their name
	for(int i=0; i<getFlavourNumber(); i++){
		string aXTitleMCflavour;
		aXTitleMCflavour = "R_{" + aZoneName + "}^{MPF," + getFlavourBinName(i) + "}";
		string aTitleMCflavour;
		aTitleMCflavour = aXTitleMCflavour + "distribution";
		string histoNameMC;
		histoNameMC = "response_Zone_Flavour_Pt/Rmpf/Rmpf_" + aZoneName + "_" + getFlavourBinName(i) + "_" + aPtName;
		//cout<<histoNameMC<<endl;
		myFlavourVector.push_back((TH1F*)afMC->Get(histoNameMC.c_str()));
		h1_style(myFlavourVector[i]);
		myFlavourVector[i]->SetXTitle(aXTitleMCflavour.c_str());
		myFlavourVector[i]->SetTitle(aTitleMCflavour.c_str());
		//normalize MC distributions to the lumi
		if(isForLumi == 1) {
			myFlavourVector[i]->Scale(getLumi());
		}
	}

	for(int i=0; i<getFlavourNumber(); i++){
		myFlavourVector[i]->SetLineColor(1);
		myFlavourVector[i]->SetFillColor(aHistoColor[i]);
	}

	float NentriesMC = 0;
	int excluded = 0;
	TH1F *hMcTot = NULL;
	if(withNoMatched == 0) {
		excluded = 2;
		hMcTot = (TH1F*)myFlavourVector[0]->Clone();
		for(int i = 1; i<getFlavourNumber()-2; i++) {
			hMcTot->Add(myFlavourVector[i]);
			//cout<<"hMcTot created"<<endl;
		}
// 		for(int i=0; i<getFlavourNumber()-excluded; i++){
// 			NentriesMC = NentriesMC + myFlavourVector[i]->Integral();
// 		}
	}
	else if(withNoMatched == 1) {
		excluded = 1;
		hMcTot = (TH1F*)myFlavourVector[getFlavourNumber()-1]->Clone();
		//NentriesMC = myFlavourVector[getFlavourNumber()-1]->Integral();		
	}

	NentriesMC = hMcTot->Integral();

	//normalize MC distributions to the number of events in data
	if(isForLumi == 0) {
		float NentriesData = hData->Integral();
		for(int i=0; i<getFlavourNumber(); i++){
			myFlavourVector[i]->Scale(NentriesData/NentriesMC);	
		}
	}
			
	if(zoneBin != 4 && zoneBin != 5) {
		hsMC->Add(myFlavourVector[zoneBin])	;
		for(int i=0; i<getFlavourNumber()-excluded; i++){
			if(i != zoneBin) {					
				hsMC->Add(myFlavourVector[i])	;
			}
		}
	}
	else {
		for(int i=0; i<getFlavourNumber()-excluded; i++){
			hsMC->Add(myFlavourVector[i])	;
		}		
	}

	string aPath;
  if (isMCTot) {
    aPath = "plotsResult/RmpfPerZonePerPt/Rmpf_" + aZoneName + "_" + aPtName + "_MC_TOT_woLogScale_physics.pdf";
  } else {
    aPath = "plotsResult/RmpfPerZonePerPt/Rmpf_" + aZoneName + "_" + aPtName + "_MC_G_woLogScale_physics.pdf";
  }
	drawDataMcComparison_Flavour(hsMC, myFlavourVector, hData, aPath, aXTitle, aTitle.c_str(), withNoMatched, 0);

  if (isMCTot) {
    aPath = "plotsResult/RmpfPerZonePerPt/Rmpf_" + aZoneName + "_" + aPtName + "_MC_TOT_wLogScale_physics.pdf";
  } else {
    aPath = "plotsResult/RmpfPerZonePerPt/Rmpf_" + aZoneName + "_" + aPtName + "_MC_G_wLogScale_physics.pdf";
  }
	drawDataMcComparison_Flavour(hsMC, myFlavourVector, hData, aPath, aXTitle, aTitle.c_str(), withNoMatched, 1);

  if (isMCTot) {
    aPath = "plotsResult/RmpfPerZonePerPt/Rmpf_" + aZoneName + "_" + aPtName + "_MC_TOT_woLogScale_physics.C";
  } else {
    aPath = "plotsResult/RmpfPerZonePerPt/Rmpf_" + aZoneName + "_" + aPtName + "_MC_G_woLogScale_physics.C";
  }
  drawDataMcComparison_Flavour(hsMC, myFlavourVector, hData, aPath, aXTitle, aTitle.c_str(), withNoMatched, 0);

  if (isMCTot) {
    aPath = "plotsResult/RmpfPerZonePerPt/Rmpf_" + aZoneName + "_" + aPtName + "_MC_TOT_wLogScale_physics.C";
  } else {
    aPath = "plotsResult/RmpfPerZonePerPt/Rmpf_" + aZoneName + "_" + aPtName + "_MC_G_wLogScale_physics.C";
  }
  drawDataMcComparison_Flavour(hsMC, myFlavourVector, hData, aPath, aXTitle, aTitle.c_str(), withNoMatched, 1);

	hMcTot->Delete();
}



void drawFlavourFraction(vector<TH1F*> vHisto, vector<int> aHistoColor, string path, string zoneName = "", bool isMCTot = false) {
//*******************************************************************************************************
// description: function which draws the evolution of the flavour fractions as a function of gammapt
// arguments: - vHisto: vector containing the gammapt histo for each flavour 
//            - aHistoColor: vector containing the flavour colour code
//            - path: where you want to save the plot
//*******************************************************************************************************
	
  int nflavours = getFlavourNumber();
	THStack* hsPtGamma = new THStack("hsPtGamma","Distribution of p_{t}^{#gamma}");
	for(int i=nflavours-2; i>=0; i--) {
		hsPtGamma->Add(vHisto[i]);
		vHisto[i]->SetLineColor(1);
		vHisto[i]->SetFillColor(aHistoColor[i]);

	}
	
	//draw the histo we have in vector to check what there is in
	TCanvas *cGammaPt = new TCanvas("cGammaPt","cGammaPt");
	cGammaPt->cd();
	gStyle->SetOptStat(0);	
	hsPtGamma->Draw("hist"); 
	hsPtGamma->SetMinimum(0.00001);
	hsPtGamma->GetXaxis()->SetTitle("p_{t}^{#gamma}");
	vHisto[nflavours-1]->SetFillColor(kBlack);
	vHisto[nflavours-1]->SetFillStyle(3001);
	vHisto[nflavours-1]->DrawCopy("e2same");
	//cGammaPt->SetLogy();
	gPad-> SetLogy();	
	gPad->RedrawAxis();
	TLegend *lGammaPt = new TLegend(0.69,0.65,0.89,0.89);
	lGammaPt->SetBorderSize(0);
	for(int i=0; i<nflavours-1; i++) {
		lGammaPt->AddEntry(vHisto[i],getFlavourBinName(i).c_str(),"f");
	}
	lGammaPt->SetFillColor(0);
	lGammaPt->Draw("SAME");

  string aName;
  if (zoneName == "") {
    aName = "plotsResult/fractionEvolution/flavourGammapt_physics";
  } else {
    aName = "plotsResult/fractionEvolution/flavourGammapt_" + zoneName + "_physics";
  }

  string aPath;
  if (isMCTot) {
    aPath = aName + "_MC_TOT.pdf";
  } else {
    aPath = aName + "_MC_G.pdf";
  }
	cGammaPt->SaveAs(aPath.c_str());
  if (isMCTot) {
    aPath = aName + "_MC_TOT.C";
  } else {
    aPath = aName + "_MC_G.C";
  }
	cGammaPt->SaveAs(aPath.c_str());

//******************************************************************************************
	//compute the flavour fractions and their errors from the TH1

  vector<TH1*> vHistoFrac;
  vHistoFrac.resize(nflavours-1);
  double binning[] = {0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 390, 420, 800};
  TH1* all_rebin = vHisto[nflavours-1]->Rebin(15, "tot", binning);
	for(int i=0; i<nflavours-1; i++) {
    std::string newName = getFlavourBinName(i);
    vHistoFrac[i] = vHisto[i]->Rebin(15, newName.c_str(), binning);
    vHistoFrac[i]->Divide(all_rebin);
   	vHistoFrac[i]->SetMarkerStyle(20);
   	vHistoFrac[i]->SetMarkerColor(aHistoColor[i]);
   	vHistoFrac[i]->SetLineColor(aHistoColor[i]);
	}

  	TCanvas* cFractionPerPtBin = new TCanvas("cFractionPerPtBin","cFractionPerPtBin");
	cFractionPerPtBin->cd();
  vHistoFrac[0]->Draw("PE1");
	for(int i=1; i<nflavours-1; i++) {
    vHistoFrac[i]->Draw("PE1same");
	}

	vHistoFrac[0]->SetTitle("Flavour fractions with respect to p_{t}^{#gamma}");
	vHistoFrac[0]->GetXaxis()->SetRangeUser(30.,800.);
	vHistoFrac[0]->GetXaxis()->SetTitle("p_{t}^{#gamma} [GeV/c]");
	vHistoFrac[0]->GetYaxis()->SetTitle("flavour fraction");
  vHistoFrac[0]->GetYaxis()->SetTitleOffset(1.7);
  vHistoFrac[0]->GetXaxis()->SetTitleOffset(1.5);
	vHistoFrac[0]->GetYaxis()->SetRangeUser(0,1);
	TLegend *lFractionPerPtBin = new TLegend(0.65,0.51,0.87,0.73);
  lFractionPerPtBin->SetBorderSize(0);
	lFractionPerPtBin->SetFillColor(0);
	lFractionPerPtBin->SetFillStyle(0);
	lFractionPerPtBin->SetTextFont(42);
	lFractionPerPtBin->SetTextSizePixels(24);
	for(int i=0; i<nflavours-1; i++) {
		lFractionPerPtBin->AddEntry(vHistoFrac[i],getFlavourBinName(i).c_str(),"p");
	}
	lFractionPerPtBin->Draw("same");
  string fullPath = path + ".pdf";
	cFractionPerPtBin->SaveAs(fullPath.c_str());
  fullPath = path + ".C";
	cFractionPerPtBin->SaveAs(fullPath.c_str());

/*  int NbinPt = vHisto[nflavours-1]->GetNbinsX();*/
  //cout<<"NbinPt : "<<NbinPt<<endl;
  //float aPtBin[NbinPt];
  //float aPtBinError[NbinPt];
  //vector<float*> vFraction;
  //vector<float*> vFractionError;
  //for(int i=0; i<nflavours-1; i++) {
    //vFraction.push_back(new float[NbinPt]);
    //vFractionError.push_back(new float[NbinPt]);
  //}
  //for ( int i=1 ;  i<=NbinPt ;  i++) {
    //aPtBin[i-1] = vHisto[nflavours-1]->GetBinCenter(i);
    //aPtBinError[i-1] = 0;
    //if(vHisto[nflavours-1]->GetBinContent(i) != 0.) {
      //for(int j=0; j<nflavours-1; j++) {
        //vFraction[j][i-1] = 100.*vHisto[j]->GetBinContent(i)/vHisto[nflavours-1]->GetBinContent(i);
        ////vFractionError[j][i-1] = 0.;
        //vFractionError[j][i-1] = 100.*sqrt(pow(vHisto[j]->GetBinError(i)/vHisto[nflavours-1]->GetBinContent(i),2) + pow(vHisto[nflavours-1]->GetBinError(i)*vHisto[j]->GetBinContent(i)/(pow(vHisto[nflavours-1]->GetBinContent(i),2)),2));
        //cout<<"aPtBin["<<i-1<<"] : "<<aPtBin[i-1]<<endl;
        //cout<<"aFraction["<<j<<"]["<<i-1<<"] : "<<vFraction[j][i-1]<<endl;
        ////cout<<"vFractionError["<<j<<"]["<<i-1<<"] : "<<vFractionError[j][i-1]<<endl;
      //}
    //}
  //}
  //vector<TGraphErrors*> vGraph;
  //TMultiGraph* mgFractionPerPtBin = new TMultiGraph();
  //for(int i=0; i<nflavours-1; i++) {
    //vGraph.push_back(new TGraphErrors(NbinPt,aPtBin, vFraction[i], aPtBinError, vFractionError[i]));
    //string graphName;
    //graphName = "gFraction_" + getFlavourBinName(i);
    //vGraph[i]->SetName(graphName.c_str());
     //vGraph[i]->SetMarkerStyle(20);
     //vGraph[i]->SetMarkerColor(aHistoColor[i]);
     //vGraph[i]->SetLineColor(aHistoColor[i]);
    //mgFractionPerPtBin->Add(vGraph[i]);
  //}

  //TCanvas* cFractionPerPtBin = new TCanvas("cFractionPerPtBin","cFractionPerPtBin");
  //cFractionPerPtBin->cd();
  //mgFractionPerPtBin->Draw("APE");
  //TMultiGraph_style(mgFractionPerPtBin);
  //mgFractionPerPtBin->SetTitle("Flavour fractions with respect to p_{t}^{#gamma}");
  //mgFractionPerPtBin->GetXaxis()->SetRangeUser(30.,800.);
  //mgFractionPerPtBin->GetXaxis()->SetTitle("p_{t}^{#gamma} [GeV/c]");
  //mgFractionPerPtBin->GetYaxis()->SetTitle("flavour fraction (%)");
  //mgFractionPerPtBin->GetYaxis()->SetRangeUser(0,100);
  //TLegend *lFractionPerPtBin = new TLegend(0.65,0.51,0.87,0.73);
  //lFractionPerPtBin->SetBorderSize(0);
  //lFractionPerPtBin->SetFillColor(0);
  //lFractionPerPtBin->SetFillStyle(0);
  //lFractionPerPtBin->SetTextFont(42);
  //lFractionPerPtBin->SetTextSizePixels(24);
  //for(int i=0; i<nflavours-1; i++) {
    //lFractionPerPtBin->AddEntry(vGraph[i],getFlavourBinName(i).c_str(),"p");
  //}
  //lFractionPerPtBin->Draw("same");
  //string fullPath = path + ".pdf";
  //cFractionPerPtBin->SaveAs(fullPath.c_str());
  //fullPath = path + ".C";
  /*cFractionPerPtBin->SaveAs(fullPath.c_str());*/
}


void drawTagger(ptBinning aPtBinning, TFile* afMC, TH1* hData, vector<int> aHistoColor, string aTaggerName, string aPtName, int ptBin, int isForLumi = 1, int withNoMatched = 1, int inLogScale = 1) {
//*************************************************************************************************
// description: function which draws the flavour tagger (CSV or QGL) distributions for data and MC
// arguments: - aPtBinning: object from the class 'ptBinning' to reconstruct the name of histo tyou want to 
//              catch up from the rootfile
//            - afMC: MC rootfile containning informations we want to extract
//            - hData: tagger distribution for data
//            - aHistoColor: vector containing the flavour colour code
//            - aTaggerName: name of the tagger (CSV or QGL) you want to study
//            - aPtName: name of the 2DTagging pt bin
//            - aCanvas: canvas you want to draw the plot in
//            - ptBin: number of the corresponding 2DTagging pt bin
//            - isForLumi: 1 if you want MC distributions normalized to the luminosity, 
//                         0 if you want MC distributions normalized to the number of data events (shape)
//            - withNoMatched: 1 if you want to include the no matched jets, 0 if you don't
//            - inLogScale: 1 if you want to have log scale on Y axis, 0 if you don't
//*************************************************************************************************
	THStack* hsMC = new THStack();
	string aXTitle;
	aXTitle = aTaggerName;
	std::stringstream PtBinningName ;
	PtBinningName << (int) aPtBinning.getBinValueInf(ptBin) << " GeV < p_{t}^{#gamma} < " << (int) aPtBinning.getBinValueSup(ptBin) << " GeV";
	string aTitle;
	aTitle = aXTitle + " distribution for " + PtBinningName.str();

	vector<TH1F*> vHisto;
	int nflavours = getFlavourNumber();
	for(int i=0; i<nflavours; i++){
		string histoNameMC;
		histoNameMC = "tagger_Flavour_Pt/" + aTaggerName + "/" + aTaggerName + "_" + getFlavourBinName(i) + "_" + aPtName;
		//cout<<histoNameMC<<endl;
		vHisto.push_back((TH1F*)afMC->Get(histoNameMC.c_str()));
		h1_style(vHisto[i]);
		if(isForLumi == 1) {
			vHisto[i]->Scale(getLumi());
		}
		vHisto[i]->SetLineColor(1);
		vHisto[i]->SetFillColor(aHistoColor[i]);
	}	

	TH1F *hMcTot = NULL;
	float NentriesMC = 0;
	int excluded = 0;
	if(withNoMatched == 0) {
		excluded = 2;
		hMcTot = (TH1F*)vHisto[0]->Clone();
		for(int i = 1; i<getFlavourNumber()-excluded; i++) {
			hMcTot->Add(vHisto[i]);
		}
// 		for(int i=0; i<nflavours-excluded; i++){
// 			NentriesMC = NentriesMC + vHisto[i]->Integral();
// 		}
	}
	else if(withNoMatched == 1) {
		excluded = 1;
		hMcTot = (TH1F*)vHisto[getFlavourNumber()-1]->Clone();
// 		NentriesMC = vHisto[nflavours-1]->Integral();		
	}

	NentriesMC = hMcTot->Integral();

	if(isForLumi == 0) {
		float NentriesData = hData->Integral();
		for(int i=0; i<nflavours; i++){
			vHisto[i]->Scale(NentriesData/NentriesMC);	
		}
	}	
	for(int i=0; i<nflavours-excluded; i++){
		hsMC->Add(vHisto[i])	;
	}		
	string aPath = "plotsResult/tagger/" + aTaggerName + "_" + aPtName + "_physics.pdf";
	drawDataMcComparison_Flavour(hsMC, vHisto, hData, aPath.c_str(), aXTitle.c_str(), aTitle.c_str(), withNoMatched, inLogScale, 0);		

	aPath = "plotsResult/tagger/" + aTaggerName + "_" + aPtName + "_physics.C";
	drawDataMcComparison_Flavour(hsMC, vHisto, hData, aPath.c_str(), aXTitle.c_str(), aTitle.c_str(), withNoMatched, inLogScale, 0);	

	hMcTot->Delete();
}

void my2DTaggingPlots () 
{

 /* TString innameMC_G="../output_rootfile/output2DTagging_MC_G.root";*/
	//TString innameMC_QCD="../output_rootfile/output2DTagging_MC_QCD.root";
	//TString innameData="../output_rootfile/output2DTagging_data.root";
	/*TString innameMC_TOT="../output_rootfile/output2DTagging_MC_TOT.root";*/

	TString innameMC_G="../output_rootfile/output2DTagging_MC_G_physics.root";
	TString innameMC_QCD="../output_rootfile/output2DTagging_MC_QCD_physics.root";
	TString innameData="../output_rootfile/output2DTagging_data.root";
	TString innameMC_TOT="../output_rootfile/output2DTagging_MC_TOT_physics.root";

// 	TString innameMC_G="../output_rootfile/output2DTagging_MC_G_stageM2.root";
// 	TString innameMC_QCD="../output_rootfile/output2DTagging_MC_QCD_stageM2.root";
// 	TString innameData="../output_rootfile/output2DTagging_data_stageM2_woResiduals.root";
// 	TString innameMC_TOT="../output_rootfile/output2DTagging_MC_TOT_stageM2.root";

	//open input file
	TFile *fMC_G=TFile::Open(innameMC_G);
	TFile *fMC_QCD=TFile::Open(innameMC_QCD);
	TFile *fData=TFile::Open(innameData);
	TFile *fMC_TOT=TFile::Open(innameMC_TOT);

  setTDRStyle();

	int isFor2DTagging = 1;
	ptBinning my2DTaggingPtBinning(isFor2DTagging);
	int nptbins = my2DTaggingPtBinning.getSize();
	int nzones = getZoneNumber();
	int nflavours = getFlavourNumber(); 
	int ntaggers = getTaggerNumber();

	int isForLumi = 0 ;
	int withNoMatched = 1 ;
	//stageM2
// 	int withNoMatched = 0 ;
	int inLogScale = 1;

	vector<int> myHistoColor = HistoColor();

	vector<vector<TH1F*> > vRmpf_PtZone_data;
	vRmpf_PtZone_data.resize(nptbins);

  for(int j=0; j<nptbins; j++) {
		vRmpf_PtZone_data[j].resize(nzones);
    for(int i=0; i<nzones; i++) {			
			string ZoneName = getZoneBinName(i);
			string PtName = my2DTaggingPtBinning.getName(j);
			string histoNameData;
			histoNameData = "response_Zone_Pt/Rmpf/Rmpf_" + ZoneName + "_" + PtName;
			//cout<<histoNameData<<endl;
			vRmpf_PtZone_data[j][i] = (TH1F*)fData->Get(histoNameData.c_str());
      drawRmpfPerZone(fMC_G, vRmpf_PtZone_data[j][i], myHistoColor, ZoneName, PtName, i, isForLumi, withNoMatched);
      drawRmpfPerZone(fMC_TOT, vRmpf_PtZone_data[j][i], myHistoColor, ZoneName, PtName, i, isForLumi, withNoMatched, true);
		}
	}

	vector<TH1F*> vGammaPt_Flavour_MC_G;
	vector<TH1F*> vGammaPt_Flavour_MC_TOT;

	//first retrieve the TH1
	for(int i=0; i<nflavours; i++){
		string histoNameMC;
		histoNameMC = "gammapt_Flavour/Gammapt_" + getFlavourBinName(i);
		//cout<<histoNameMC<<endl;
		vGammaPt_Flavour_MC_G.push_back((TH1F*)fMC_G->Get(histoNameMC.c_str()));
		vGammaPt_Flavour_MC_TOT.push_back((TH1F*)fMC_TOT->Get(histoNameMC.c_str()));
		h1_style(vGammaPt_Flavour_MC_G[i]);
		h1_style(vGammaPt_Flavour_MC_TOT[i]);
	}	

  drawFlavourFraction(vGammaPt_Flavour_MC_G, myHistoColor, "plotsResult/fractionEvolution/flavourFraction_MC_G_physics");
  drawFlavourFraction(vGammaPt_Flavour_MC_TOT, myHistoColor, "plotsResult/fractionEvolution/flavourFraction_MC_TOT_physics", "", true);

	vector<vector<TH1F*> > vGammaPt_ZoneFlavour_MC_G;
	vector<vector<TH1F*> > vGammaPt_ZoneFlavour_MC_TOT;
  vGammaPt_ZoneFlavour_MC_G.resize(nzones);
  vGammaPt_ZoneFlavour_MC_TOT.resize(nzones);

	for(int j=0; j<nzones; j++) {
		vGammaPt_ZoneFlavour_MC_G[j].resize(nflavours);
		vGammaPt_ZoneFlavour_MC_TOT[j].resize(nflavours);
		for(int i=0; i<nflavours; i++) {			
			string ZoneName = getZoneBinName(j);
			string FlavourName = getFlavourBinName(i);
			string histoName;
			histoName = "variables/gammaPtPerZonePerFlavour/Gammapt_" + ZoneName + "_" + FlavourName;
			//cout<<histoName<<endl;
			vGammaPt_ZoneFlavour_MC_G[j][i] = (TH1F*)fMC_G->Get(histoName.c_str());
			vGammaPt_ZoneFlavour_MC_TOT[j][i] = (TH1F*)fMC_TOT->Get(histoName.c_str());
		}
	}

	for(int j=0; j<nzones; j++) {
			string ZoneName = getZoneBinName(j);
			string histoName = "plotsResult/fractionEvolution/flavourFraction_MC_G_physics_" + ZoneName;
      drawFlavourFraction(vGammaPt_ZoneFlavour_MC_G[j], myHistoColor, histoName, ZoneName);
			histoName = "plotsResult/fractionEvolution/flavourFraction_MC_TOT_physics_" + ZoneName;
      drawFlavourFraction(vGammaPt_ZoneFlavour_MC_TOT[j], myHistoColor, histoName, ZoneName, true);
	}


	vector<vector<TH1F*> > vTag_TaggerZone_data;
	vTag_TaggerZone_data.resize(nptbins);

	for(int j=0; j<nptbins; j++) {
		vTag_TaggerZone_data[j].resize(nzones);
		for(int i=0; i<ntaggers; i++) {			
			string TaggerName = getTaggerBinName(i);
			string PtName = my2DTaggingPtBinning.getName(j);
			string histoNameData;
			histoNameData = "tagger_Pt/" + TaggerName + "/" + TaggerName + "_" + PtName;
			//cout<<histoNameData<<endl;
			vTag_TaggerZone_data[j][i] = (TH1F*)fData->Get(histoNameData.c_str());
			if(TaggerName=="QGL") inLogScale = 0;
			else if(TaggerName=="CSV") inLogScale = 1;
      drawTagger(my2DTaggingPtBinning, fMC_G, vTag_TaggerZone_data[j][i], myHistoColor, TaggerName, PtName, j, isForLumi, withNoMatched, inLogScale);
		}
	}

//*************************************************************************************************
//
//                                      gamma pt 
//
//*************************************************************************************************

	TH1F* hGammaPt_MC_G=(TH1F*)fMC_G->Get("variables/hGammaPt");
	TH1F* hGammaPt_MC_QCD=(TH1F*)fMC_QCD->Get("variables/hGammaPt");
	TH1F* hGammaPt_data=(TH1F*)fData->Get("variables/hGammaPt");
	TH1F* hGammaPt_MC_tot=(TH1F*)fMC_TOT->Get("variables/hGammaPt");

	h1_style(hGammaPt_MC_G);	
	h1_style(hGammaPt_MC_QCD);
	h1_style(hGammaPt_data);
	h1_style(hGammaPt_MC_tot);
	
	hGammaPt_data->SetXTitle("p_{T}^{#gamma} [GeV/c]");
	hGammaPt_MC_QCD->SetXTitle("p_{T}^{#gamma} [GeV/c]");
	hGammaPt_MC_G->SetXTitle("p_{T}^{#gamma} [GeV/c]");
	hGammaPt_data->SetTitle("Distribution en p_{T}^{#gamma}");
	
	//rescale the Monte Carlo histogramm with luminosity
	hGammaPt_MC_G->Scale(getLumi());
	hGammaPt_MC_QCD->Scale(getLumi());
	hGammaPt_MC_tot->Scale(getLumi());
		
	THStack* hsMCgammapt = new THStack("hsMCgammapt","hsMCgammapt");
	hsMCgammapt->Add(hGammaPt_MC_G);
	hsMCgammapt->Add(hGammaPt_MC_QCD);	

	hsMCgammapt->SetMinimum(1);
	hGammaPt_MC_G->SetFillColor(2);
	hGammaPt_MC_G->SetLineColor(1);
	hGammaPt_MC_QCD->SetFillColor(5);
	hGammaPt_MC_QCD->SetLineColor(1);

  drawDataMcComparisonWithoutRatio("c1bis", hsMCgammapt, hGammaPt_MC_G, hGammaPt_MC_QCD, hGammaPt_MC_tot, hGammaPt_data, "p_{T}^{#gamma} [GeV/c]", "plotsResult/gammapt/cGammaPtTot2013_lumi_withoutRatio_physics.pdf");

  drawDataMcComparisonWithoutRatio("c1bis", hsMCgammapt, hGammaPt_MC_G, hGammaPt_MC_QCD, hGammaPt_MC_tot, hGammaPt_data, "p_{T}^{#gamma} [GeV/c]", "plotsResult/gammapt/cGammaPtTot2013_lumi_withoutRatio_physics.C");

//*************************************************************************************************
//
//                                      alpha
//
//*************************************************************************************************

	TH1F* hAlpha_MC_G=(TH1F*)fMC_G->Get("variables/hAlpha");
	TH1F* hAlpha_MC_QCD=(TH1F*)fMC_QCD->Get("variables/hAlpha");
	TH1F* hAlpha_data=(TH1F*)fData->Get("variables/hAlpha");
	TH1F* hAlpha_MC_tot=(TH1F*)fMC_TOT->Get("variables/hAlpha");

	h1_style(hAlpha_MC_G);	
	h1_style(hAlpha_MC_QCD);
	h1_style(hAlpha_data);
	h1_style(hAlpha_MC_tot);
	
	hAlpha_data->SetXTitle("#alpha");
	hAlpha_MC_QCD->SetXTitle("#alpha");
	hAlpha_MC_G->SetXTitle("#alpha");
	hAlpha_data->SetTitle("#alpha distribution");
	
	//rescale the Monte Carlo histogramm with luminosity
	hAlpha_MC_G->Scale(getLumi());
	hAlpha_MC_QCD->Scale(getLumi());
	hAlpha_MC_tot->Scale(getLumi());
		
	THStack* hsMCalpha = new THStack("hsMCalpha","hsMCalpha");
	hsMCalpha->Add(hAlpha_MC_QCD);	
	hsMCalpha->Add(hAlpha_MC_G);
	hsMCalpha->SetMinimum(1);
	hAlpha_MC_G->SetFillColor(2);
	hAlpha_MC_G->SetLineColor(1);
	hAlpha_MC_QCD->SetFillColor(5);
	hAlpha_MC_QCD->SetLineColor(1);

  drawDataMcComparisonWithoutRatio("c2", hsMCalpha, hAlpha_MC_G, hAlpha_MC_QCD, hAlpha_MC_tot, hAlpha_data, "#alpha", "plotsResult/alpha/cAlphaTot2013_lumi_withoutRatio_physics.pdf",0);
  drawDataMcComparisonWithoutRatio("c2", hsMCalpha, hAlpha_MC_G, hAlpha_MC_QCD, hAlpha_MC_tot, hAlpha_data, "#alpha", "plotsResult/alpha/cAlphaTot2013_lumi_withoutRatio_physics.C",0);

}




