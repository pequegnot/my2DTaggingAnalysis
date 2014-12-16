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
#include "TMatrixF.h"
#include "TMatrixFSym.h"
#include "TMatrixFSparse.h"
#include "TMatrixFLazy.h"
#include "TVectorF.h"
#include <TRandom.h>
#include <THStack.h>
#include "TROOT.h"

#include "binning.h"
#include "common.h"

#include "ptBinning.h"
#include "alphaBinning.h"
//#include "fitTools.h"
#include "QGSyst.h"

using namespace std;

int main () 
{
//open input file
  	//TFile *f=TFile::Open("input_rootfile/PhotonJet_Photon_Run2012_residuals_PFlowAK5chs.root"); 
		TFile *f=TFile::Open("input_rootfile/PhotonJet_Photon_Run2012_residuals_alpha_0.3_PFlowAK5chs.root"); 

 	
  	

//*********************************************************************************************************
//
//                                      histograms definition
//
//*********************************************************************************************************

	ptBinning myPtBinning;
	int isFor2DTagging = 1;
	ptBinning my2DTaggingPtBinning(isFor2DTagging);
  int n2DTaggingPtBins = my2DTaggingPtBinning.getSize();
	
//usefull variables
	Double_t xlow = getHistoXlow();
	Double_t xup = getHistoXup();
	//Double_t    binrange = 0.1;
	//Int_t    nbinsx = (xup - xlow)/binrange;
	Int_t    nbinsx = getHistoNbinsx();
	
	Double_t xlow2D = 0;
	Double_t xup2D = 1;
	Int_t    nbinsx2D = 20;
	Double_t ylow2D = 0;
	Double_t yup2D = 1;
	Int_t    nbinsy2D = 20;

	int nzones = getZoneNumber();
	int nflavours = getFlavourNumber(); //uds, g, c, b, noMatched, all

//vectors for 2D tagging study
	//responses per zone per pt
	vector<vector<TH1F*> > vRmpf_ZonePt = buildZonePtVectorH1(my2DTaggingPtBinning,"Rmpf",nbinsx,xlow,xup) ;

	//QG-likelihood per pt
	vector<TH1F*> vQGL_Pt = buildPtVectorH1(my2DTaggingPtBinning,"QGL",nbinsx2D,xlow2D,xup2D) ;

	//Btag-CSV per pt
	vector<TH1F*> vCSV_Pt = buildPtVectorH1(my2DTaggingPtBinning,"CSV",nbinsx2D,xlow2D,xup2D) ;

	
	for(int j=0; j<my2DTaggingPtBinning.getSize(); j++) {
		vQGL_Pt[j]->Sumw2();
		vCSV_Pt[j]->Sumw2();
		for(int i=0; i<nzones; i++) {
			vRmpf_ZonePt[i][j]->Sumw2();
		}
	}

	TH1F* hGammaPt=new TH1F("hGammaPt","hGammaPt",40,40,800);
	hGammaPt->SetXTitle("p_{T}^{#gamma} [GeV/c]");
	hGammaPt->Sumw2();
	
	TH1F* hFirstJetPt=new TH1F("hFirstJetPt","hFirstJetPt",50,0,1400);
	hFirstJetPt->SetXTitle("1^{st} jet p_{T} [GeV/c]");
	hFirstJetPt->Sumw2();

	TH1F* hSecondJetPt=new TH1F("hSecondJetPt","hSecondJetPt",50,0,200);
	hSecondJetPt->SetXTitle("2^{nd} jet p_{T} [GeV/c]");
	hSecondJetPt->Sumw2();
	
	TH1F* hMetPt=new TH1F("hMetPt","hMetPt",50,0,620);
	hMetPt->SetXTitle("MET p_{T} [GeV/c]");
	hMetPt->Sumw2();
	
	TH1F* hParalleleMetPt=new TH1F("hParalleleMetPt","hParalleleMetPt",50,-200,200);
	hParalleleMetPt->SetXTitle("parallele MET p_{T} [GeV/c]");
	hParalleleMetPt->Sumw2();
	
	TH1F* hPerpendicularMetPt=new TH1F("hPerpendicularMetPt","hPerpendicularMetPt",50,-140,140);
	hPerpendicularMetPt->SetXTitle("perpendicular MET p_{T} [GeV/c]");
	hPerpendicularMetPt->Sumw2();
	
	TH1F* hAlpha=new TH1F("hAlpha","hAlpha",10,0.,1.);
	hAlpha->SetXTitle("p_{T}^{2^{nd} jet}/p_{T}^{#gamma}");
	hAlpha->Sumw2();
	
	TH1F* hDeltaPhi_j1gamma=new TH1F("hDeltaPhi_j1gamma","hDeltaPhi_j1gamma",40,0,4);
	hDeltaPhi_j1gamma->SetXTitle("|#Delta#phi (#gamma,1^{st}jet)|");
	hDeltaPhi_j1gamma->Sumw2();

	TH1F* hDeltaPhi_j2gamma=new TH1F("hDeltaPhi_j2gamma","hDeltaPhi_j2gamma",40,0,4);
	hDeltaPhi_j2gamma->SetXTitle("|#Delta#phi (#gamma,2^{nd}jet)|");
	hDeltaPhi_j2gamma->Sumw2();
	
	TH1F* hDeltaPhi_j1j2=new TH1F("hDeltaPhi_j1j2","hDeltaPhi_j1j2",40,0,4);
	hDeltaPhi_j1j2->SetXTitle("|#Delta#phi (1^{st} jet,2^{nd} jet)|");
	hDeltaPhi_j1j2->Sumw2();
	
	TH1F* hDeltaR_j1gamma=new TH1F("hDeltaR_j1gamma","hDeltaR_j1gamma",40,0,8);
	hDeltaR_j1gamma->SetXTitle("|#DeltaR (#gamma,1^{st} jet)|");
	hDeltaR_j1gamma->Sumw2();
	
	TH1F* hDeltaR_j2gamma=new TH1F("hDeltaR_j2gamma","hDeltaR_j2gamma",40,0,8);
	hDeltaR_j2gamma->SetXTitle("|#DeltaR (#gamma,2^{nd} jet)|");
	hDeltaR_j2gamma->Sumw2();
	
	TH1F* hDeltaR_j1j2=new TH1F("hDeltaR_j1j2","hDeltaR_j1j2",40,0,8);
	hDeltaR_j1j2->SetXTitle("|#DeltaR (1^{st} jet,2^{nd} jet)|");
	hDeltaR_j1j2->Sumw2();




//*****************************************************************************************************
//
//                                      reading the root file 
//
//*****************************************************************************************************

//retrieve the trees
	TTree* t_firstjet=(TTree*)f->Get("first_jet");
	TTree* t_firstjetraw=(TTree*)f->Get("first_jet_raw");
	TTree* t_secondjet=(TTree*)f->Get("second_jet");
	TTree* t_electron=(TTree*)f->Get("electrons");
	TTree* t_muon=(TTree*)f->Get("muons");
	TTree* t_met=(TTree*)f->Get("met");
	TTree* t_gamma=(TTree*)f->Get("photon");
	TTree* t_misc=(TTree*)f->Get("misc");
	TTree* t_rho=(TTree*)f->Get("rho");
	//TTree* t_gen_particles=(TTree*)f->Get("gen_particles");

//retrieve the variables
	//First jet
	float firstjetpt;
	t_firstjet->SetBranchAddress("pt",&firstjetpt);
	float firstjetpx;
	t_firstjet->SetBranchAddress("px",&firstjetpx);
	float firstjetpy;
	t_firstjet->SetBranchAddress("py",&firstjetpy);
	float firstjetphi;
	t_firstjet->SetBranchAddress("phi",&firstjetphi);
	float firstjeteta;
	t_firstjet->SetBranchAddress("eta",&firstjeteta);
	float firstjetbtag_csv;
	t_firstjet->SetBranchAddress("btag_csv",&firstjetbtag_csv);
	float firstjetqg_tag_mlp;
	t_firstjet->SetBranchAddress("qg_tag_mlp",&firstjetqg_tag_mlp);
	float firstjetqg_tag_likelihood;
	t_firstjet->SetBranchAddress("qg_tag_likelihood",&firstjetqg_tag_likelihood);
	
	//First jet raw
	float firstjetrawpt;
	t_firstjetraw->SetBranchAddress("pt",&firstjetrawpt);
	float firstjetrawphi;
	t_firstjetraw->SetBranchAddress("phi",&firstjetrawphi);
	float firstjetraweta;
	t_firstjetraw->SetBranchAddress("eta",&firstjetraweta);
	
	//misc
	float miscevent_weight;
	t_misc->SetBranchAddress("event_weight",&miscevent_weight);

	//Second jet
	float secondjetpt;
	t_secondjet->SetBranchAddress("pt",&secondjetpt);
	float secondjetphi;
	t_secondjet->SetBranchAddress("phi",&secondjetphi);
	float secondjeteta;
	t_secondjet->SetBranchAddress("eta",&secondjeteta);
	
	//Gamma (reco)
	float gammapt;
	t_gamma->SetBranchAddress("pt",&gammapt);
	float gammaphi;
	t_gamma->SetBranchAddress("phi",&gammaphi);
	float gammaeta;
	t_gamma->SetBranchAddress("eta",&gammaeta);
	float gammapx;
	t_gamma->SetBranchAddress("px",&gammapx);
	float gammapy;
	t_gamma->SetBranchAddress("py",&gammapy);
	float gammapz;
	t_gamma->SetBranchAddress("pz",&gammapz);
	bool gamma_has_pixel_seed;
	t_gamma->SetBranchAddress("has_pixel_seed",&gamma_has_pixel_seed);
	
	//electron
	float electronpt[100];
	t_electron->SetBranchAddress("pt",&electronpt);
	float electronphi[100];
	t_electron->SetBranchAddress("phi",&electronphi);
	float electroneta[100];
	t_electron->SetBranchAddress("eta",&electroneta);
	int electronn;
	t_electron->SetBranchAddress("n",&electronn);
	float electronisolation[100];
	t_electron->SetBranchAddress("isolation",&electronisolation);
	
	//muon
	float muonpt[100];
	t_muon->SetBranchAddress("pt",&muonpt);
	float muonphi[100];
	t_muon->SetBranchAddress("phi",&muonphi);
	float muoneta[100];
	t_muon->SetBranchAddress("eta",&muoneta);
	int muonn;
	t_muon->SetBranchAddress("n",&muonn);
	float muonisolation[100];
	t_muon->SetBranchAddress("relative_isolation",&muonisolation);
	
	//met
	float metpt;
	t_met->SetBranchAddress("pt",&metpt);
	float metpx;
	t_met->SetBranchAddress("px",&metpx);
	float metpy;
	t_met->SetBranchAddress("py",&metpy);
	float metpz;
	t_met->SetBranchAddress("pz",&metpz);
	float metphi;
	t_met->SetBranchAddress("phi",&metphi);
	
	//rho
	double rho;
	t_rho->SetBranchAddress("rho",&rho);
	

	
	//Usefull variables
	float Rmpf;
	float alpha;
	float metParal;
	float metPerpend;
	Double_t DeltaPhi_j1gamma;
	Double_t DeltaPhi_j2gamma;
	Double_t DeltaPhi_j1j2;	
	Double_t DeltaPhi_j1met;
	Double_t DeltaR_j1gamma;
	Double_t DeltaR_j2gamma;
	Double_t DeltaR_j1j2;	
	Double_t DeltaEta_j1gamma;
	Double_t DeltaEta_j2gamma;
	Double_t DeltaEta_j1j2;

	int binPt;//bin en pt
	int bin2DTaggingPt;//gros bins en pt pour l etude 2D tagging
	int binZone;// 2D tagging zone bin
	
	bool dropEvent=false;
	
	//count events in the tree
	unsigned int nEvents = (int)t_firstjet->GetEntries();

	double totalWeight;
	TFile *f_weight=TFile::Open("output_rootfile_last/additionalPrescaleWeight.root");
	TH1F* hprescaleWeight=(TH1F*)f_weight->Get("hprescaleWeight");
	
	//loop over them
	for(unsigned int ievt=0; ievt<nEvents; ievt++) {
		t_firstjet->GetEntry(ievt);
		t_firstjetraw->GetEntry(ievt);
		t_secondjet->GetEntry(ievt);
		t_electron->GetEntry(ievt);
		t_muon->GetEntry(ievt);
		t_gamma->GetEntry(ievt);
		t_met->GetEntry(ievt);
		t_misc->GetEntry(ievt);
		t_rho->GetEntry(ievt);   
	
//*********************************************************************************************************
		
		//hAlpha->Fill(alpha,miscevent_weight);

		dropEvent=false;
		
		if(gamma_has_pixel_seed == true) continue;
		  	
		binPt = myPtBinning.getPtBin(gammapt);
		if(binPt == -1) continue;

		bin2DTaggingPt = my2DTaggingPtBinning.getPtBin(gammapt);		

		//totalWeight = miscevent_weight*(hprescaleWeight->GetBinContent(hprescaleWeight->FindBin(gammapt)));
		
		binZone = getZoneBin(firstjetbtag_csv, firstjetqg_tag_likelihood);
		//if(binZone == -1) continue;		

		Rmpf = 1 + (gammapx*metpx + gammapy*metpy)/(pow(gammapt,2));		
		alpha = (secondjetpt)/(gammapt);

		
//*****************************************************************************************************
//
//                                      deltaPhi calculation
//
//*****************************************************************************************************


//DeltaPhi between the 1st jet and the gamma calculation
    	DeltaPhi_j1gamma = TMath::Abs((gammaphi) - (firstjetphi));
    	if(DeltaPhi_j1gamma>TMath::Pi()){
    	  DeltaPhi_j1gamma = 2*TMath::Pi()-DeltaPhi_j1gamma;
    	}
    	
//DeltaPhi between the 2nd jet and the gamma calculation
    	DeltaPhi_j2gamma = TMath::Abs((gammaphi) - (secondjetphi));
    	if(DeltaPhi_j2gamma>TMath::Pi()){
    	  DeltaPhi_j2gamma = 2*TMath::Pi()-DeltaPhi_j2gamma;
    	}


//DeltaPhi between the first jet and the second jet calculation
    	DeltaPhi_j1j2 = TMath::Abs((secondjetphi) - (firstjetphi));
    	if(DeltaPhi_j1j2>TMath::Pi()){
      	DeltaPhi_j1j2 = 2*TMath::Pi()-DeltaPhi_j1j2;
    	}
    	

 //DeltaPhi between the first jet and the met
    	DeltaPhi_j1met = TMath::Abs((metphi) - (firstjetphi));
    	
    	
//*****************************************************************************************************
//
//                                      deltaEta calculation 
//
//*****************************************************************************************************

//DeltaEta between the first jet and the gamma calculation
    	DeltaEta_j1gamma = TMath::Abs((gammaeta) - (firstjeteta));

//DeltaEta between the second jet and the gamma calculation
    	DeltaEta_j2gamma = TMath::Abs((gammaeta) - (secondjeteta));
    	
//DeltaEta between the first jet and the second jet calculation
    	DeltaEta_j1j2 = TMath::Abs((secondjeteta) - (firstjeteta));

    	    	
//*****************************************************************************************************
//
//                                      deltaR calculation 
//
//*****************************************************************************************************


//DeltaR between the first jet and the gamma calculation
    	DeltaR_j1gamma = sqrt ( pow(DeltaEta_j1gamma,2) + pow(DeltaPhi_j1gamma,2) );
    	
//DeltaR between the second jet and the gamma calculation
    	DeltaR_j2gamma = sqrt ( pow(DeltaEta_j2gamma,2) + pow(DeltaPhi_j2gamma,2) );

//DeltaR between the first jet and the second jet calculation
    	DeltaR_j1j2 = sqrt ( pow(DeltaEta_j1j2,2) + pow(DeltaPhi_j1j2,2) );

    	
//*****************************************************************************************************
//
//                                      met parallele calculation 
//
//*****************************************************************************************************

		metParal = metpt * TMath::Cos(DeltaPhi_j1met);

//*****************************************************************************************************
//
//                                      met perpendicular calculation 
//
//*****************************************************************************************************

		metPerpend = metpt * TMath::Sin(DeltaPhi_j1met);

//*****************************************************************************************************
//
//                                      filling histogramms 
//
//*****************************************************************************************************
	

		//if(fabs(firstjetpt)>12. && DeltaPhi_j1gamma>2.8 && fabs(firstjeteta)<1.3 && (alpha<0.3 || secondjetpt<10)) {
		if(fabs(firstjetpt)>12. && DeltaPhi_j1gamma>2.8 && fabs(firstjeteta)<1.3) {
			if(muonn==0){
				if(electronn==0) {
					if(gammapt>200.) {
						hAlpha->Fill(alpha,miscevent_weight);
					}
					if((alpha<0.3 || secondjetpt<10)) {
						hGammaPt->Fill(gammapt,miscevent_weight);
						hFirstJetPt->Fill(firstjetpt,miscevent_weight);
						hMetPt->Fill(metpt,miscevent_weight);
						hParalleleMetPt->Fill(metParal,miscevent_weight);
						hPerpendicularMetPt->Fill(metPerpend,miscevent_weight);
						hSecondJetPt->Fill(secondjetpt,miscevent_weight);
						hDeltaPhi_j1gamma->Fill(DeltaPhi_j1gamma,miscevent_weight);
						hDeltaPhi_j2gamma->Fill(DeltaPhi_j2gamma,miscevent_weight);
						hDeltaPhi_j1j2->Fill(DeltaPhi_j1j2,miscevent_weight);
						hDeltaR_j1gamma->Fill(DeltaR_j1gamma,miscevent_weight);
						hDeltaR_j2gamma->Fill(DeltaR_j2gamma,miscevent_weight);
						hDeltaR_j1j2->Fill(DeltaR_j1j2,miscevent_weight);	
		
						if((firstjetbtag_csv>=0 && firstjetbtag_csv<=1) && (firstjetqg_tag_likelihood>=0 && firstjetqg_tag_likelihood<=1)) {
						//protection contre les jets non matches
							if(binZone!=-1) {
								vRmpf_ZonePt[binZone][bin2DTaggingPt]->Fill(Rmpf,miscevent_weight);
							}
							vQGL_Pt[bin2DTaggingPt]->Fill(firstjetqg_tag_likelihood,miscevent_weight);
							vCSV_Pt[bin2DTaggingPt]->Fill(firstjetbtag_csv,miscevent_weight);
						}
					}
				}
			}
		}
	}
    
	

	

//*********************************************************************************************************
//
//                                      Output file
//
//*********************************************************************************************************

	TFile *out_mikko = new TFile("output_rootfile/output2DTagging_data_forMikko.root", "recreate");

	out_mikko->cd();	
TDirectory *response_Zone_PtDir_mikko = out_mikko->mkdir("response_perZone","response_perZone");
	TDirectory *Rmpf_Zone_PtDir_mikko = response_Zone_PtDir_mikko->mkdir("Rmpf","Rmpf");
	Rmpf_Zone_PtDir_mikko->cd();
  		for(int i=0; i<nzones; i++) {
			vRmpf_ZonePt[i][n2DTaggingPtBins-1]->Write();
		}
  out_mikko->Close();

	//create the output file and write into it
	TFile *out = new TFile("output_rootfile/output2DTagging_data.root", "recreate");
	
	out->cd();	
	TDirectory *response_Zone_PtDir = out->mkdir("response_Zone_Pt","response_Zone_Pt");
	TDirectory *Rmpf_Zone_PtDir = response_Zone_PtDir->mkdir("Rmpf","Rmpf");
	Rmpf_Zone_PtDir->cd();
	for(int j=0; j<my2DTaggingPtBinning.getSize(); j++) {
		for(int i=0; i<nzones; i++) {
			vRmpf_ZonePt[i][j]->Write();
		}
	}


	TDirectory *tagger_PtDir = out->mkdir("tagger_Pt","tagger_Pt");
	TDirectory *CSV_PtDir = tagger_PtDir->mkdir("CSV","CSV");
	CSV_PtDir->cd();
	for(int j=0; j<my2DTaggingPtBinning.getSize(); j++) {
			vCSV_Pt[j]->Write();
	}
	TDirectory *QGL_Flavour_PtDir = tagger_PtDir->mkdir("QGL","QGL");
	QGL_Flavour_PtDir->cd();
	for(int j=0; j<my2DTaggingPtBinning.getSize(); j++) {
			vQGL_Pt[j]->Write();
	}
	
	
	TDirectory *variablesDir = out->mkdir("variables","variables");
	variablesDir->cd();
	hGammaPt->Write();
	hFirstJetPt->Write();
	hSecondJetPt->Write();
	hMetPt->Write();
	hParalleleMetPt->Write();
	hPerpendicularMetPt->Write();
	hAlpha->Write();
	hDeltaPhi_j1gamma->Write();
	hDeltaPhi_j2gamma->Write();
	hDeltaPhi_j1j2->Write();
	hDeltaR_j1gamma->Write();
	hDeltaR_j2gamma->Write();
	hDeltaR_j1j2->Write();
	
	out->Close();
	
	return 0;
}









