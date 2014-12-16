#include <iomanip>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <TVector3.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH2F.h>
#include <TColor.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveStats.h>
#include <TStyle.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TRandom.h>
#include <TMultiGraph.h>
#include <THStack.h>
#include "TLatex.h"
#include "TString.h"

#pragma once

#define TITLE_FONTSIZE 26
#define LABEL_FONTSIZE 18
#define LEFT_MARGIN 0.17
#define RIGHT_MARGIN 0.03
#define TOP_MARGIN 0.05
#define BOTTOM_MARGIN 0.13

TStyle* createStyle() {
  TStyle *style = new TStyle("style", "style");

  // For the canvas:
  style->SetCanvasBorderMode(0);
  style->SetCanvasColor(kWhite);
  style->SetCanvasDefH(800); //Height of canvas
  style->SetCanvasDefW(800); //Width of canvas
  style->SetCanvasDefX(0);   //POsition on screen
  style->SetCanvasDefY(0);

  // For the Pad:
  style->SetPadBorderMode(0);
  style->SetPadColor(kWhite);
  style->SetPadGridX(false);
  style->SetPadGridY(false);
  style->SetGridColor(0);
  style->SetGridStyle(3);
  style->SetGridWidth(1);

  // For the frame:
  style->SetFrameBorderMode(0);
  style->SetFrameBorderSize(1);
  style->SetFrameFillColor(0);
  style->SetFrameFillStyle(0);
  style->SetFrameLineColor(1);
  style->SetFrameLineStyle(1);
  style->SetFrameLineWidth(1);

  // For the histo:
  style->SetHistLineStyle(1);
  style->SetHistLineWidth(2);
  style->SetEndErrorSize(2);
  style->SetMarkerStyle(20);

  //For the fit/function:
  style->SetFitFormat("5.4g");
  style->SetFuncColor(2);
  style->SetFuncStyle(1);
  style->SetFuncWidth(1);

  // For the statistics box:
  style->SetOptFile(0);
  style->SetStatColor(kWhite);
  //style->SetStatFont(43);
  //style->SetStatFontSize(0.025);
  style->SetStatTextColor(1);
  style->SetStatFormat("6.4g");
  style->SetStatBorderSize(1);
  //style->SetStatH(0.1);
  //style->SetStatW(0.15);

  //For the date:
  style->SetOptDate(0);

  // Margins:
  style->SetPadTopMargin(TOP_MARGIN);
  style->SetPadBottomMargin(BOTTOM_MARGIN);
  style->SetPadLeftMargin(LEFT_MARGIN);
  style->SetPadRightMargin(RIGHT_MARGIN);

  // For the Global title:
  style->SetOptTitle(0);
  style->SetTitleFont(63);
  style->SetTitleColor(1);
  style->SetTitleFillColor(10);
  style->SetTitleBorderSize(0);
  style->SetTitleAlign(33);
  style->SetTitleX(1);
  style->SetTitleFontSize(TITLE_FONTSIZE);

  // For the axis titles:

  style->SetTitleColor(1, "XYZ");
  style->SetTitleFont(43, "XYZ");
  style->SetTitleSize(TITLE_FONTSIZE, "XYZ");
  //style->SetTitleXOffset(3.5); //FIXME
  style->SetTitleYOffset(2.5); //FIXME
  style->SetTitleXOffset(1.5);
  //style->SetTitleYOffset(1.5);

  style->SetLabelColor(1, "XYZ");
  style->SetLabelFont(43, "XYZ");
  style->SetLabelOffset(0.01, "YZ");
  style->SetLabelOffset(0.015, "X");
  style->SetLabelSize(LABEL_FONTSIZE, "XYZ");

  style->SetAxisColor(1, "XYZ");
  style->SetStripDecimals(kTRUE);
  style->SetTickLength(0.03, "XYZ");
  style->SetNdivisions(510, "XYZ");
  style->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  style->SetPadTickY(1);

  style->SetOptLogx(0);
  style->SetOptLogy(0);
  style->SetOptLogz(0);

  style->SetHatchesSpacing(1.3);
  style->SetHatchesLineWidth(1);

  style->cd();

  return style;
}

void cms_style(bool isData = false){
  std::string status = "Simulation preliminary";
  if (isData) status = "Preliminary";
  TPaveText* pt_exp = new TPaveText(LEFT_MARGIN, 1 - 0.5 * TOP_MARGIN, 1 - RIGHT_MARGIN, 1, "brNDC");
  pt_exp->SetFillStyle(0);
  pt_exp->SetBorderSize(0);
  pt_exp->SetMargin(0);
  pt_exp->SetTextFont(62);
  pt_exp->SetTextSize(0.75 * TOP_MARGIN);
  pt_exp->SetTextAlign(13);
  TString d = TString::Format("CMS #font[52]{#scale[0.76]{%s}}", status.c_str());
  pt_exp->AddText(d);
  pt_exp->Draw();

  TString lumi_s = "19.48 fb^{-1} (8 TeV)";
  TPaveText* pt_lumi = new TPaveText(LEFT_MARGIN, 1 - 0.5 * TOP_MARGIN, 1 - RIGHT_MARGIN, 1, "brNDC");
  pt_lumi->SetFillStyle(0);
  pt_lumi->SetBorderSize(0);
  pt_lumi->SetMargin(0);
  pt_lumi->SetTextFont(42);
  pt_lumi->SetTextSize(0.6 * TOP_MARGIN);
  pt_lumi->SetTextAlign(33);
  pt_lumi->AddText(lumi_s);
  pt_lumi->Draw();

}

void h1_style(TH1 *h, int optstat=0) {
	h->SetStats(optstat);	
	h->SetLabelFont(42,"X");       // 42
	h->SetLabelFont(42,"Y");       // 42
	h->SetLabelOffset(0.005,"X");  // D=0.005
	h->SetLabelOffset(0.005,"Y");  // D=0.005
	h->SetLabelSize(0.045,"X");
	h->SetLabelSize(0.045,"Y");
	h->SetTitleOffset(1.15,"X");
	h->SetTitleOffset(1.15,"Y");
	h->SetTitleSize(0.04,"X");
	h->SetTitleSize(0.04,"Y");
	//h->SetTitle(0);
	h->SetTitleFont(42, "XYZ");
}

void drawHist(const string& canvasName, TH1 *hist,const string& path) {
	TCanvas *cCanvas = new TCanvas(canvasName.c_str(),canvasName.c_str());
	cCanvas->cd();
	hist->Draw("hist");
	cCanvas->SaveAs(path.c_str());
}

void TMultiGraph_style (TMultiGraph* h) {
	h->GetXaxis()->SetLabelOffset(0.005);
	h->GetXaxis()->SetLabelFont(42);
	h->GetXaxis()->SetLabelSize(0.055);
	h->GetXaxis()->SetTitleOffset(1.15);
	h->GetXaxis()->SetTitleSize(0.04);
	h->GetXaxis()->SetTitleFont(42);
	//h->GetYaxis()->SetRangeUser(0.,1.);
	h->GetYaxis()->SetLabelOffset(0.005);
	h->GetYaxis()->SetLabelFont(42);
	h->GetYaxis()->SetLabelSize(0.045);	
	h->GetYaxis()->SetTitleOffset(1.2);
	h->GetYaxis()->SetTitleFont(42);
	h->GetYaxis()->SetTitleSize(0.04);	
}

void TGraph_style (TGraph* h) {
	h->GetXaxis()->SetLabelOffset(0.005);
	h->GetXaxis()->SetLabelFont(42);
	h->GetXaxis()->SetLabelSize(0.055);
	h->GetXaxis()->SetTitleOffset(1.15);
	h->GetXaxis()->SetTitleSize(0.04);
	h->GetXaxis()->SetTitleFont(42);
	//h->GetYaxis()->SetRangeUser(0.,1.);
	h->GetYaxis()->SetLabelOffset(0.005);
	h->GetYaxis()->SetLabelFont(42);
	h->GetYaxis()->SetLabelSize(0.045);	
	h->GetYaxis()->SetTitleOffset(1.2);
	h->GetYaxis()->SetTitleFont(42);
	h->GetYaxis()->SetTitleSize(0.04);	
}


void fractionPerPtBin () 
{	
// 	TFile *f1=TFile::Open("output_rootfile/MC_G_out_neutrinosStudy.root");
// 	
// 	TH1F *hGammaPt_allFlavour=(TH1F*)f1->Get("neutrinos_study/Pt_neutrinos/hGammaPt_allFlavour");
// 	TH1F *hGammaPt_uds=(TH1F*)f1->Get("neutrinos_study/Pt_neutrinos/hGammaPt_uds");
// 	TH1F *hGammaPt_c=(TH1F*)f1->Get("neutrinos_study/Pt_neutrinos/hGammaPt_c");
// 	TH1F *hGammaPt_b=(TH1F*)f1->Get("neutrinos_study/Pt_neutrinos/hGammaPt_b");
// 	TH1F *hGammaPt_g=(TH1F*)f1->Get("neutrinos_study/Pt_neutrinos/hGammaPt_g");
	
// 	TFile *f1=TFile::Open("output_rootfile/MC_G_finalizer_allPtGamma_partonFlavour.root");
	TFile *f1=TFile::Open("output_rootfile/output2DTagging_MC_G_physics.root");

  //applyStyle();
  TStyle* m_style = createStyle();
  m_style->cd();

  Int_t g_color = TColor::GetColor("#542437");
  Int_t uds_color = TColor::GetColor("#D95B43");
  Int_t c_color = TColor::GetColor("#C02942");
  Int_t b_color = TColor::GetColor("#53777A");
  Int_t noMatched_color = TColor::GetColor("#69D2E7");
		
  TH1F *hGammaPt_allFlavour=(TH1F*)f1->Get("gammapt_Flavour/Gammapt_all");
	TH1F* hGammaPt_g=(TH1F*)f1->Get("gammapt_Flavour/Gammapt_g");
	TH1F* hGammaPt_c=(TH1F*)f1->Get("gammapt_Flavour/Gammapt_c");
	TH1F* hGammaPt_b=(TH1F*)f1->Get("gammapt_Flavour/Gammapt_b");
	TH1F* hGammaPt_uds=(TH1F*)f1->Get("gammapt_Flavour/Gammapt_uds");
	TH1F* hGammaPt_noMatched=(TH1F*)f1->Get("gammapt_Flavour/Gammapt_noMatched");
	
 /* h1_style(hGammaPt_allFlavour);*/
	//h1_style(hGammaPt_g);
	//h1_style(hGammaPt_uds);
	//h1_style(hGammaPt_c);
	/*h1_style(hGammaPt_b);*/
	
	THStack* hsPtGamma = new THStack("hsPtGamma","Distribution of p_{t}^{#gamma}");
	hsPtGamma->Add(hGammaPt_noMatched);
	hsPtGamma->Add(hGammaPt_b);
	hsPtGamma->Add(hGammaPt_c);
	hsPtGamma->Add(hGammaPt_g);
	hsPtGamma->Add(hGammaPt_uds);
	hGammaPt_g->SetFillColor(g_color);
	hGammaPt_g->SetLineColor(1);
	hGammaPt_c->SetFillColor(c_color);
	hGammaPt_c->SetLineColor(1);	
	hGammaPt_b->SetFillColor(b_color);
	hGammaPt_b->SetLineColor(1);
	hGammaPt_noMatched->SetFillColor(noMatched_color);
	hGammaPt_noMatched->SetLineColor(1);
	hGammaPt_uds->SetFillColor(uds_color);
	hGammaPt_uds->SetLineColor(1);	
	
	TCanvas *cGammaPt = new TCanvas("cGammaPt","cGammaPt");
	cGammaPt->cd();
	gStyle->SetOptStat(0);	
	hsPtGamma->Draw("hist"); 
	hsPtGamma->SetMinimum(0.0001);
	hsPtGamma->GetXaxis()->SetTitle("p_{t}^{#gamma}");
	hGammaPt_allFlavour->SetFillColor(kBlack);
	hGammaPt_allFlavour->SetFillStyle(3001);
	hGammaPt_allFlavour->DrawCopy("e2same");
	//cGammaPt->SetLogy();
	gPad-> SetLogy();	
	gPad->RedrawAxis();
	TLegend *lGammaPt = new TLegend(0.69,0.65,0.89,0.89);
	lGammaPt->SetBorderSize(0);
	lGammaPt->AddEntry(hGammaPt_g,"gluon","f");
	lGammaPt->AddEntry(hGammaPt_c,"c","f");
	lGammaPt->AddEntry(hGammaPt_b,"b","f");
	lGammaPt->AddEntry(hGammaPt_uds,"uds","f");
	lGammaPt->AddEntry(hGammaPt_noMatched,"no matched","f");
	lGammaPt->SetFillColor(0);
	lGammaPt->Draw("SAME");
  cms_style(false);
	cGammaPt->SaveAs("images2DTagging/fractionPerPtBin/cGammaPt_MC_G_physics.pdf");
	cGammaPt->SaveAs("images2DTagging/fractionPerPtBin/cGammaPt_MC_G_physics.C");
	
	int NbinPt = hGammaPt_allFlavour->GetNbinsX();
  std::cout << "NbinPt: " << NbinPt << std::endl;
	double *aPtBin = new double[NbinPt];
	double *aPtBinError = new double[NbinPt];
	double *aFraction_uds = new double[NbinPt];
	double *aFraction_c = new double[NbinPt];
	double *aFraction_b = new double[NbinPt];
	double *aFraction_noMatched = new double[NbinPt];
	double *aFraction_g = new double[NbinPt];
	double *aFractionError_uds = new double[NbinPt];
	double *aFractionError_c = new double[NbinPt];
	double *aFractionError_b = new double[NbinPt];
	double *aFractionError_noMatched = new double[NbinPt];
	double *aFractionError_g = new double[NbinPt];
	
	for ( int i=1 ;  i<NbinPt ;  i++) {
		aPtBin[i-1] = hGammaPt_allFlavour->GetBinCenter(i);
		aPtBinError[i-1] = 0;
		
		cout<<"hGammaPt_allFlavour->GetBinContent(i) : "<<hGammaPt_allFlavour->GetBinContent(i)<<endl;
		
		if(hGammaPt_allFlavour->GetBinContent(i) != 0.) {
			aFraction_uds[i-1] = 100.*hGammaPt_uds->GetBinContent(i)/hGammaPt_allFlavour->GetBinContent(i);
			aFraction_c[i-1] = 100.*hGammaPt_c->GetBinContent(i)/hGammaPt_allFlavour->GetBinContent(i);
			aFraction_b[i-1] = 100.*hGammaPt_b->GetBinContent(i)/hGammaPt_allFlavour->GetBinContent(i);
			aFraction_noMatched[i-1] = 100.*hGammaPt_noMatched->GetBinContent(i)/hGammaPt_allFlavour->GetBinContent(i);
			aFraction_g[i-1] = 100.*hGammaPt_g->GetBinContent(i)/hGammaPt_allFlavour->GetBinContent(i);	
			
			cout<<"aPtBin["<<i-1<<"] : "<<aPtBin[i-1]<<endl;
			cout<<"aFraction_uds["<<i-1<<"] : "<<aFraction_uds[i-1]<<endl;
			cout<<"aFraction_c["<<i-1<<"] : "<<aFraction_c[i-1]<<endl;
			cout<<"aFraction_b["<<i-1<<"] : "<<aFraction_b[i-1]<<endl;
			cout<<"aFraction_noMatched["<<i-1<<"] : "<<aFraction_noMatched[i-1]<<endl;
			cout<<"aFraction_g["<<i-1<<"] : "<<aFraction_g[i-1]<<endl;
			
	/*		aFractionError_uds[i-1] =0;
			
			aFractionError_c[i-1] =0;
			
			aFractionError_b[i-1] =0;

			aFractionError_noMatched[i-1] =0;
			
			aFractionError_g[i-1] =0;
	*/
	
			aFractionError_uds[i-1] = 100.*sqrt(pow(hGammaPt_uds->GetBinError(i)/hGammaPt_allFlavour->GetBinContent(i),2) + pow(hGammaPt_allFlavour->GetBinError(i)*hGammaPt_uds->GetBinContent(i)/(pow(hGammaPt_allFlavour->GetBinContent(i),2)),2));
			
			aFractionError_c[i-1] = 100.*sqrt(pow(hGammaPt_c->GetBinError(i)/hGammaPt_allFlavour->GetBinContent(i),2) + pow(hGammaPt_allFlavour->GetBinError(i)*hGammaPt_c->GetBinContent(i)/(pow(hGammaPt_allFlavour->GetBinContent(i),2)),2));
			
			aFractionError_b[i-1] = 100.*sqrt(pow(hGammaPt_b->GetBinError(i)/hGammaPt_allFlavour->GetBinContent(i),2) + pow(hGammaPt_allFlavour->GetBinError(i)*hGammaPt_b->GetBinContent(i)/(pow(hGammaPt_allFlavour->GetBinContent(i),2)),2));


			aFractionError_noMatched[i-1] = 100.*sqrt(pow(hGammaPt_noMatched->GetBinError(i)/hGammaPt_allFlavour->GetBinContent(i),2) + pow(hGammaPt_allFlavour->GetBinError(i)*hGammaPt_noMatched->GetBinContent(i)/(pow(hGammaPt_allFlavour->GetBinContent(i),2)),2));
			
			aFractionError_g[i-1] = 100.*sqrt(pow(hGammaPt_g->GetBinError(i)/hGammaPt_allFlavour->GetBinContent(i),2) + pow(hGammaPt_allFlavour->GetBinError(i)*hGammaPt_g->GetBinContent(i)/(pow(hGammaPt_allFlavour->GetBinContent(i),2)),2));
			
			cout<<"aFractionError_uds["<<i-1<<"] : "<<aFractionError_uds[i-1]<<endl;
			cout<<"aFractionError_c["<<i-1<<"] : "<<aFractionError_c[i-1]<<endl;
			cout<<"aFractionError_b["<<i-1<<"] : "<<aFractionError_b[i-1]<<endl;
			cout<<"aFractionError_noMatched["<<i-1<<"] : "<<aFractionError_noMatched[i-1]<<endl;
			cout<<"aFractionError_g["<<i-1<<"] : "<<aFractionError_g[i-1]<<endl;
		}
	}
	
	TGraphErrors* gFraction_uds = new TGraphErrors(NbinPt,&aPtBin[0], &aFraction_uds[0], &aPtBinError[0], &aFractionError_uds[0]);
   	gFraction_uds->SetName("gFraction_uds");
   	gFraction_uds->SetMarkerStyle(20);
   	gFraction_uds->SetMarkerColor(uds_color);
   	gFraction_uds->SetLineColor(uds_color);
   	
   	TGraphErrors* gFraction_c = new TGraphErrors(NbinPt,aPtBin, aFraction_c, aPtBinError, aFractionError_c);
   	gFraction_c->SetName("gFraction_c");
   	gFraction_c->SetMarkerStyle(20);
   	gFraction_c->SetMarkerColor(c_color);
   	gFraction_c->SetLineColor(c_color);
   	
   	TGraphErrors* gFraction_b = new TGraphErrors(NbinPt,aPtBin, aFraction_b, aPtBinError, aFractionError_b);
   	gFraction_b->SetName("gFraction_b");
   	gFraction_b->SetMarkerStyle(20);
   	gFraction_b->SetMarkerColor(b_color);
   	gFraction_b->SetLineColor(b_color);

   	TGraphErrors* gFraction_noMatched = new TGraphErrors(NbinPt, aPtBin, aFraction_noMatched, aPtBinError, aFractionError_noMatched);
   	gFraction_noMatched->SetName("gFraction_noMatched");
   	gFraction_noMatched->SetMarkerStyle(20);
   	gFraction_noMatched->SetMarkerColor(noMatched_color);
   	gFraction_noMatched->SetLineColor(noMatched_color);
   	
   	TGraphErrors* gFraction_g = new TGraphErrors(NbinPt, aPtBin, aFraction_g, aPtBinError, aFractionError_g);
   	gFraction_g->SetName("gFraction_g");
   	gFraction_g->SetMarkerStyle(20);
   	gFraction_g->SetMarkerColor(g_color);
   	gFraction_g->SetLineColor(g_color);
	
	TMultiGraph* mgFractionPerPtBin = new TMultiGraph();
	mgFractionPerPtBin->Add(gFraction_uds);
	mgFractionPerPtBin->Add(gFraction_c);
	mgFractionPerPtBin->Add(gFraction_b);
	mgFractionPerPtBin->Add(gFraction_noMatched);
	mgFractionPerPtBin->Add(gFraction_g);
	TCanvas* cFractionPerPtBin = new TCanvas("cFractionPerPtBin","cFractionPerPtBin");
	cFractionPerPtBin->cd();
	mgFractionPerPtBin->Draw("APE");
	//TMultiGraph_style(mgFractionPerPtBin);
	mgFractionPerPtBin->SetTitle("Flavour fractions with respect to p_{t}^{#gamma}");
	mgFractionPerPtBin->GetXaxis()->SetRangeUser(30.,800.);
	mgFractionPerPtBin->GetXaxis()->SetTitle("p_{t}^{#gamma} [GeV/c]");
	mgFractionPerPtBin->GetYaxis()->SetTitle("flavour fraction (%)");
	TLegend *lFractionPerPtBin = new TLegend(0.708,0.540,0.871,0.727);
  lFractionPerPtBin->SetBorderSize(0);
	lFractionPerPtBin->SetFillColor(0);
	lFractionPerPtBin->SetFillStyle(0);
	lFractionPerPtBin->SetTextFont(42);
	lFractionPerPtBin->SetTextSizePixels(24);
	lFractionPerPtBin->AddEntry(gFraction_uds,"uds","p");
	lFractionPerPtBin->AddEntry(gFraction_c,"c","p");
	lFractionPerPtBin->AddEntry(gFraction_b,"b","p");
	lFractionPerPtBin->AddEntry(gFraction_g,"g","p");
	lFractionPerPtBin->AddEntry(gFraction_noMatched,"no matched","p");
	lFractionPerPtBin->Draw("same");
  cms_style(false);
  cFractionPerPtBin->SaveAs("images2DTagging/fractionPerPtBin/cFractionPerPtBin_MC_G_physics.pdf");
	cFractionPerPtBin->SaveAs("images2DTagging/fractionPerPtBin/cFractionPerPtBin_MC_G_physics.C");

  delete[] aPtBin;
	delete[] aPtBinError;
	delete[] aFraction_uds;
	delete[] aFraction_c;
	delete[] aFraction_b;
	delete[] aFraction_noMatched;
	delete[] aFraction_g;
	delete[] aFractionError_uds;
	delete[] aFractionError_c;
	delete[] aFractionError_b;
	delete[] aFractionError_noMatched;
	delete[] aFractionError_g;

}
























