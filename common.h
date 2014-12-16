#pragma once

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <TVector3.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TF1.h>
#include <TH2F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <TMath.h>
#include <TLegend.h>
#include <TStyle.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <stdio.h>
#include <TString.h>
#include <TColor.h>
#include "TMatrixF.h"
#include "TMatrixFSym.h"
#include "TMatrixFSparse.h"
#include "TMatrixFLazy.h"
#include "TVectorF.h"
#include <TRandom.h>
#include <THStack.h>
#include "TROOT.h"
#include<math.h>

#include "binning.h"
#include "ptBinning.h"


float getLumi () {
  return 1.948*pow(10.,4.);
}

int getHistoNbinsx() {
  return 30;
}

double getHistoXlow() {
  return 0.;
}

double getHistoXup() {
  return 2.;
}

std::vector<TH1F*> buildPtVectorH1(ptBinning aPtBinning,const std::string& aName,Int_t nbinsx,Double_t xlow,Double_t xup)
{
  std::vector<TH1F*> myVector;
  for(int i=0; i<aPtBinning.getSize(); i++) {
    std::string histoName;
    histoName = aName + "_" + aPtBinning.getName(i) ;
    myVector.push_back(new TH1F(histoName.c_str(),histoName.c_str(),nbinsx, xlow, xup));
  }
  return myVector;
}

std::vector<TH1F*> buildFlavourVectorH1(const std::string& aName,Int_t nbinsx,Double_t xlow,Double_t xup)
{
  std::vector<TH1F*> myVector;
  for(int i=0; i<getFlavourNumber(); i++){
    std::string histoName;
    histoName = aName + "_" + getFlavourBinName(i);
    myVector.push_back(new TH1F(histoName.c_str(),histoName.c_str(),nbinsx, xlow, xup));
  }
  return myVector;
}

std::vector<std::vector<TH1F*> > buildFlavourPtVectorH1(ptBinning aPtBinning,const std::string& aName,Int_t nbinsx,Double_t xlow,Double_t xup)
{
  std::vector<std::vector<TH1F*> > myVector;
  for(int i=0; i<getFlavourNumber(); i++){
    std::string histoName;
    histoName = aName + "_" + getFlavourBinName(i);
    myVector.push_back(buildPtVectorH1(aPtBinning,histoName.c_str(),nbinsx, xlow, xup));
  }
  return myVector;
}

std::vector<std::vector<TH1F*> > buildZonePtVectorH1(ptBinning aPtBinning,const std::string& aName,Int_t nbinsx,Double_t xlow,Double_t xup)
{
  std::vector<std::vector<TH1F*> > myVector;
  for(int i=0; i<getZoneNumber(); i++){
    std::string histoName;
    histoName = aName + "_" + getZoneBinName(i);
    myVector.push_back(buildPtVectorH1(aPtBinning,histoName.c_str(),nbinsx, xlow, xup));
  }
  return myVector;
}

std::vector<std::vector<TH1F*> > buildZoneFlavourVectorH1(const std::string& aName,Int_t nbinsx,Double_t xlow,Double_t xup)
{
  std::vector<std::vector<TH1F*> > myVector;
  for(int i=0; i<getZoneNumber(); i++){
    std::string histoName;
    histoName = aName + "_" + getZoneBinName(i);
    myVector.push_back(buildFlavourVectorH1(histoName.c_str(),nbinsx, xlow, xup));
  }
  return myVector;
}

std::vector<std::vector<std::vector<TH1F*> > > buildZoneFlavourPtVectorH1(ptBinning aPtBinning,const std::string& aName,Int_t nbinsx,Double_t xlow,Double_t xup)
{
  std::vector<std::vector<std::vector<TH1F*> > > myVector;
  for(int i=0; i<getZoneNumber(); i++){
    std::string histoName;
    histoName = aName + "_" + getZoneBinName(i);
    myVector.push_back(buildFlavourPtVectorH1(aPtBinning, histoName.c_str(), nbinsx, xlow, xup));
  }
  return myVector;
}

std::vector<TH2F*> buildPtVectorH2(ptBinning aPtBinning,const std::string& aName,Int_t nbinsx,Double_t xlow,Double_t xup,Int_t nbinsy,Double_t ylow,Double_t yup)
{
  std::vector<TH2F*> myVector;
  for(int i=0; i<aPtBinning.getSize(); i++) {
    std::string histoName;
    histoName = aName + "_" + aPtBinning.getName(i) ;
    myVector.push_back(new TH2F(histoName.c_str(),histoName.c_str(),nbinsx, xlow, xup, nbinsy, ylow, yup));
  }
  return myVector;
}

std::vector<std::vector<TH2F*> > buildFlavourPtVectorH2(ptBinning aPtBinning,const std::string& aName,Int_t nbinsx,Double_t xlow,Double_t xup,Int_t nbinsy,Double_t ylow,Double_t yup)
{
  std::vector<std::vector<TH2F*> > myVector;
  for(int i=0; i<getFlavourNumber(); i++){
    std::string histoName;
    histoName = aName + "_" + getFlavourBinName(i);
    myVector.push_back(buildPtVectorH2(aPtBinning,histoName.c_str(),nbinsx, xlow, xup, nbinsy, ylow, yup));
  }
  return myVector;
}

std::vector<std::vector<std::vector<float> > > build_fraction_ZoneFlavourPtVector(ptBinning aPtBinning)
{
  std::vector<std::vector<std::vector<float> > > myFinalVector;
  std::vector<std::vector<float> > myVector_2;
  std::vector<float> myVector_1;

  for(int k=0; k<aPtBinning.getSize(); k++) {
    myVector_1.push_back(0.);
  }

  for(int j=0; j<getFlavourNumber(); j++) {//j=4 : noMatched ; j=5 : total
    myVector_2.push_back(myVector_1);
  }

  for(int i=0; i<getZoneNumber(); i++){//6 zones
    myFinalVector.push_back(myVector_2);
  }
  return myFinalVector;
}

std::vector<TMatrixF> buildMatrixPtVector(ptBinning aPtBinning)
{
  std::vector<TMatrixF> myVector;
  for(int i=0; i<aPtBinning.getSize(); i++) {
    int ncolumn = getFlavourNumber() - 2; //number of flavours in the matrix
    int nrow = getZoneNumber();
    TMatrixF m(nrow,ncolumn);
    myVector.push_back(m);
  }
  return myVector;
}

std::vector<TMatrixF> buildSquareMatrixPtVector(ptBinning aPtBinning)
{
  std::vector<TMatrixF> myVector;
  for(int i=0; i<aPtBinning.getSize(); i++) {
    int ncolumn = getFlavourNumber() - 2; //number of flavours in the matrix
    int nrow = getZoneNumber() - 2;//we do not use the mix zones
    TMatrixF m(nrow,ncolumn);
    myVector.push_back(m);
  }
  return myVector;
}


// //common tagger values
int getZoneBin(float x_csv, float y_qgl)
{
  if(x_csv<0.3 && y_qgl>=0.9) {
    return 0;
  }
  else if(x_csv<0.3 && y_qgl<0.1) {
    return 1;
  }
  else if(x_csv>=0.7 && x_csv<0.9) {
    return 2;
  }
  else if(x_csv>=0.9) {
    return 3;
  }
  else if(x_csv<0.3 && y_qgl>=0.1 && y_qgl<0.9) {
    return 4;
  }
  else if(x_csv>=0.3 && x_csv<0.7) {
    return 5;
  }
  else {
    return -1;
  }
}

// test purity gluon
// int getZoneBin(float x_csv, float y_qgl)
// {
// 	if(x_csv<0.25 && y_qgl>=0.9) {
// 		return 0;
// 	}
// 	else if(x_csv<0.25 && y_qgl<0.05) {
// 		return 1;
// 	}
// 	else if(x_csv>=0.7 && x_csv<0.9) {
// 		return 2;
// 	}
// 	else if(x_csv>=0.9) {
// 		return 3;
// 	}
// 	else if(x_csv<0.25 && y_qgl>=0.05 && y_qgl<0.9) {
// 		return 4;
// 	}
// 	else if(x_csv>=0.25 && x_csv<0.7) {
// 		return 5;
// 	}
// 	else {
// 		return -1;
// 	}
// }

//stage M2 previous tagger values for gammapt>100 GeV
/*int getZoneBin(float x_csv, float y_qgl)
{
  if(x_csv<0.10 && y_qgl>0.85) {
    return 0;
  }
  else if(x_csv<0.3 && y_qgl<0.05) {
    return 1;
  }
  else if(x_csv>0.7 && x_csv<0.9 && y_qgl>0.55) {
    return 2;
  }
  else if(x_csv>0.95) {
    return 3;
  }
  else if(x_csv<0.3 && y_qgl>=0.1 && y_qgl<0.9) {
    return 4;
  }
  else if(x_csv>=0.3 && x_csv<0.7) {
    return 5;
  }
  else {
    return -1;
  }
}*/

//****************************************************************************************
//                                 !!!! WARNING !!!!
//        if you change the correspondance binName-binNumber in binning.h
//            do not forget to change it in the getFlavourBin function 
//****************************************************************************************

int getFlavourBin(float afirstjetgenpdgid)
{
  if(TMath::Abs(afirstjetgenpdgid) == 1 || TMath::Abs(afirstjetgenpdgid) == 2 || TMath::Abs(afirstjetgenpdgid) == 3) {//uds quarks
    return 0;
  }	
  else if(TMath::Abs(afirstjetgenpdgid) == 21) {//gluons
    return 1;
  }
  else if(TMath::Abs(afirstjetgenpdgid) == 4) {//c
    return 2;
  }
  else if(TMath::Abs(afirstjetgenpdgid) == 5) {//b
    return 3;
  }
  else {//no matched
    return 4;
  }	
}

//*********************************************************************************************************

void countFlavour(int ZoneBin, float& N_zone1, float& N_zone2, float& N_zone3, float& N_zone4, float& N_zone5, float& N_zone6, float weight) {
  if(ZoneBin == 0) {
    N_zone1 = N_zone1 + weight;
  }
  else if(ZoneBin == 1) {
    N_zone2 = N_zone2 + weight;
  }
  else if(ZoneBin == 2) {
    N_zone3 = N_zone3 + weight;
  }
  else if(ZoneBin == 3) {
    N_zone4 = N_zone4 + weight;
  }
  else if(ZoneBin == 4) {
    N_zone5 = N_zone5 + weight;
  }
  else if(ZoneBin == 5) {
    N_zone6 = N_zone6 + weight;
  }		
}	

void printLatexTableFlavourRatio(float Nuds, float Ng, float Nc, float Nb, float NnoMatched, float Ntot,TH1 *h, const std::string& path, ptBinning aPtBinning, int zoneNumber=0, int ptNumber=0) {
  TString titleName;
  std::stringstream PtBinningName;
  PtBinningName << (int) aPtBinning.getBinValueInf(ptNumber) << " GeV < gammapt < " << (int) aPtBinning.getBinValueSup(ptNumber)<< " GeV";
  if (zoneNumber == 0) {
    std::cout<<"******************************************************************************************"<<std::endl;
    std::cout<<PtBinningName.str()<<std::endl;
    // Print latex table
    std::cout << "\\begin{tabular}{@{}|l|c|c|c|c|c|@{}}" << std::endl << "\\hline" << std::endl;
    std::cout << "\\diagbox{Zone}{Flavour} & uds & gluon & c & b & no matched" << "\\\\" <<std::endl << "\\hline" << std::endl;
  }
  std::cout << std::setprecision(3);
  std::cout << getZoneBinNiceName(zoneNumber) << " & " <<Nuds/Ntot*100.<<" \\% &	"<<Ng/Ntot*100.<<" \\% &	"<<Nc/Ntot*100.<<" \\% &	"<<Nb/Ntot*100.<<" \\% &	"<<NnoMatched/Ntot*100.<<" \\% " << "\\\\" <<std::endl;
  if (zoneNumber == getZoneNumber () - 1) {
    std::cout << "\\hline" << std::endl << "\\end{tabular}" << std::endl;;
  }

}

void computeFlavourRatio(float Nuds, float Ng, float Nc, float Nb, float NnoMatched, float Ntot,TH1 *h, const std::string& path, ptBinning aPtBinning, int zoneNumber=0, int ptNumber=0) {
  TString titleName;
  std::stringstream PtBinningName;
  PtBinningName << (int) aPtBinning.getBinValueInf(ptNumber) << " GeV < gammapt < " << (int) aPtBinning.getBinValueSup(ptNumber)<< " GeV";

  std::cout<<"******************************************************************************************"<<std::endl;
  std::cout<<PtBinningName.str()<<std::endl;
  std::cout<<"Zone Number "<<(zoneNumber+1)<<std::endl;
  titleName = TString::Format("Pourcentage de jets dans la zone%i pour %i GeV < p_{t}^{#gamma} < %i GeV",zoneNumber,(int) aPtBinning.getBinValueInf(ptNumber),(int) aPtBinning.getBinValueSup(ptNumber));

  std::cout<<"Ntot"<<"	"<<"Nuds"<<"	"<<"Ng"<<"	"<<"Nc"<<"	"<<"Nb"<<"	"<<"NnoMatched"<<std::endl;
  std::cout<<Ntot*getLumi()<<"	"<<Nuds*getLumi()<<"	"<<Ng*getLumi()<<"	"<<Nc*getLumi()<<"	"<<Nb*getLumi()<<"	"<<NnoMatched*getLumi()<<std::endl<<std::endl;
  std::cout<<"Ratio"<<"	"<<"Nuds"<<"	"<<"Ng"<<"	"<<"Nc"<<"	"<<"Nb"<<"	"<<"NnoMatched"<<std::endl;
  std::cout<<"%"<<"	"<<Nuds/Ntot*100.<<"	"<<Ng/Ntot*100.<<"	"<<Nc/Ntot*100.<<"	"<<Nb/Ntot*100.<<"	"<<NnoMatched/Ntot*100.<<std::endl;
  std::cout<<"******************************************************************************************"<<std::endl;

  h->SetBinContent(1,Nuds/Ntot*100.);
  h->SetBinContent(2,Ng/Ntot*100.);
  h->SetBinContent(3,Nc/Ntot*100.);
  h->SetBinContent(4,Nb/Ntot*100.);
  h->SetBinContent(5,NnoMatched/Ntot*100.);
  h->GetXaxis()->SetNdivisions(5,kTRUE);
  h->GetXaxis()->SetBinLabel(1,"quarks u,d,s");
  h->GetXaxis()->SetBinLabel(2,"gluon");
  h->GetXaxis()->SetBinLabel(3,"quark c");
  h->GetXaxis()->SetBinLabel(4,"quark b");
  h->GetXaxis()->SetBinLabel(5,"no matched");
  h->SetStats(0);
  //h->GetXaxis()->SetRangeUser(0.,1.);
  h->GetXaxis()->SetLabelOffset(0.005);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelSize(0.045);
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
  h->GetYaxis()->SetTitle("jet flavour/all mixture ratio (%)");
  h->SetTitle(titleName);

  TCanvas *cFlavour = new TCanvas("cFlavour","cFlavour");
  cFlavour->cd();
  cFlavour->SetLogy(1);
  h->Draw();
  cFlavour->SaveAs(path.c_str());
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

void TGraph_style2 (TGraph* h, const std::string& Xtitle, const std::string& Ytitle, int color, int style, int size=1) {
  h->GetXaxis()->SetLabelOffset(0.005);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelSize(0.055);
  h->GetXaxis()->SetTitleOffset(1.15);
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitle(Xtitle.c_str());
  //h->GetYaxis()->SetRangeUser(0.,1.);
  h->GetYaxis()->SetLabelOffset(0.005);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelSize(0.045);	
  h->GetYaxis()->SetTitleOffset(1.2);
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetTitle(Ytitle.c_str());	

  h->SetMarkerColor(color);
  h->SetMarkerStyle(style);
  h->SetMarkerSize(size);
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

void h1_style_lo(TH1 *h, int optstat=0) {
  h->SetStats(optstat);
  h->SetXTitle("Photon p_{T} [GeV/c]");
  h->SetYTitle("Data / MC");
  h->GetXaxis()->SetTitleOffset(1.2);
  h->GetYaxis()->SetTitleOffset(0.55);
  h->GetXaxis()->SetTickLength(0.06);
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  //h->GetXaxis()->SetLabelSize(0.);
  h->GetXaxis()->SetLabelSize(0.085);
  h->GetYaxis()->SetLabelSize(0.07);
  h->GetXaxis()->SetTitleSize(0.09);
  h->GetYaxis()->SetTitleSize(0.08);	
  //h->SetTitleFont(42, "XYZ");
}

void h1_style_hi(TH1 *h, int optstat=0) {
  h->GetXaxis()->SetTitleOffset(1.1);
  h->GetYaxis()->SetTitleOffset(1.2);
  h->GetYaxis()->SetTitleSize(0.045);
  //h->GetXaxis()->SetMoreLogLabels();
  //h->GetXaxis()->SetNoExponent();
  h->GetXaxis()->SetLabelSize(0.);
}

void hs_style_hi(THStack *h, int optstat=0) {
  h->GetXaxis()->SetTitleOffset(1.1);
  h->GetYaxis()->SetTitleOffset(1.2);
  h->GetYaxis()->SetTitleSize(0.045);
  //h->GetXaxis()->SetMoreLogLabels();
  //h->GetXaxis()->SetNoExponent();
  h->GetXaxis()->SetLabelSize(0.);
}

std::vector<int> MarkerStyle() {
  std::vector<int> MarkerStyle;
  MarkerStyle.resize(getFlavourNumber());
  //uds
  MarkerStyle[0] = 20;
  //g
  MarkerStyle[1] = 20;
  //c
  MarkerStyle[2] = 20;
  //b
  MarkerStyle[3] = 20;
  //no matched
  MarkerStyle[4] = 20;
  //all
  MarkerStyle[5] = 20;

  return MarkerStyle;
}

std::vector<int> HistoColor() {
  std::vector<int> HistoColor;
  HistoColor.resize(getFlavourNumber());
  //uds
  HistoColor[0] = TColor::GetColor("#D95B43");
  //g
  HistoColor[1] = TColor::GetColor("#542437");
  //c
  HistoColor[2] = TColor::GetColor("#C02942");
  //b
  HistoColor[3] = TColor::GetColor("#53777A");
  //no matched
  HistoColor[4] = TColor::GetColor("#69D2E7");
  //all
  HistoColor[5] = 6;

  return HistoColor;
}

