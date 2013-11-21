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
  	TFile *f=TFile::Open("input_rootfile/PhotonJet_G_PFlowAK5chs.root");
// 		TFile *f=TFile::Open("input_rootfile/PhotonJet_QCD_PFlowAK5chs.root");
// 		TFile *f=TFile::Open("input_rootfile/PhotonJet_MC_TOT_PFlowAK5chs.root");
 	
  	

//*****************************************************************************************************
//
//                                      histograms definition
//
//*****************************************************************************************************
	QGSyst qgsyst;
	//qgsyst.ReadDatabase("SystDatabase.txt");
	qgsyst.ReadDatabaseDoubleMin("SystDatabase_doubleMin.txt");
	qgsyst.SetTagger("QGLHisto");

	ptBinning myPtBinning;
	int isFor2DTagging = 1;
	ptBinning my2DTaggingPtBinning(isFor2DTagging);
	
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
	int n2DTaggingPtBins = my2DTaggingPtBinning.getSize();

	vector<TMatrixF> v4x4MatrixPt = buildSquareMatrixPtVector(my2DTaggingPtBinning);
	vector<TMatrixF> v6x4MatrixPt = buildMatrixPtVector(my2DTaggingPtBinning);

//vectors for 2D tagging study
	//flavour fractions per pt
	vector<TH1F*> vFractionHisto_Pt = buildPtVectorH1(my2DTaggingPtBinning,"FractionHisto",nflavours-1,0,nflavours-1) ;

	//responses per zone per pt
	vector<vector<TH1F*> > vRmpf_ZonePt = buildZonePtVectorH1(my2DTaggingPtBinning,"Rmpf",nbinsx,xlow,xup) ;

	vector<vector<TH1F*> > vRtrue_ZonePt = buildZonePtVectorH1(my2DTaggingPtBinning,"Rtrue",nbinsx,xlow,xup) ;


	//QG-likelihood per flavour per pt
	vector<vector<TH1F*> > vQGL_FlavourPt = buildFlavourPtVectorH1(my2DTaggingPtBinning,"QGL",nbinsx2D,xlow2D,xup2D) ;

	//Btag-CSV per flavour per pt
	vector<vector<TH1F*> > vCSV_FlavourPt = buildFlavourPtVectorH1(my2DTaggingPtBinning,"CSV",nbinsx2D,xlow2D,xup2D) ;

	//responses per flavour per pt
	vector<vector<TH1F*> > vRmpf_FlavourPt = buildFlavourPtVectorH1(my2DTaggingPtBinning,"Rmpf",nbinsx,xlow,xup) ;

	vector<vector<TH1F*> > vRtrue_FlavourPt = buildFlavourPtVectorH1(my2DTaggingPtBinning,"Rtrue",nbinsx,xlow,xup) ;

	//responses per flavour per pt when we are in one of the 2D tagging zones
	vector<vector<TH1F*> > vRmpf_in2DTaggingZone_FlavourPt = buildFlavourPtVectorH1(my2DTaggingPtBinning,"Rmpf_in2DTaggingZone",nbinsx,xlow,xup) ;

	vector<vector<TH1F*> > vRtrue_in2DTaggingZone_FlavourPt = buildFlavourPtVectorH1(my2DTaggingPtBinning,"Rtrue_in2DTaggingZone",nbinsx,xlow,xup) ;

	
	//Gammapt per flavour 
	vector<TH1F*> vGammapt_Flavour = buildFlavourVectorH1("Gammapt",30,0,800);


	//responses per zone per flavour per pt
	vector<vector<vector<TH1F*> > > vRmpf_ZoneFlavourPt = buildZoneFlavourPtVectorH1(my2DTaggingPtBinning,"Rmpf",nbinsx,xlow,xup) ;

	vector<vector<vector<TH1F*> > > vRtrue_ZoneFlavourPt = buildZoneFlavourPtVectorH1(my2DTaggingPtBinning,"Rtrue",nbinsx,xlow,xup) ;

	
	//2D tagging plans per flavour per pt
	vector<vector<TH2F*> > v2DTaggingPlan_FlavourPt = buildFlavourPtVectorH2(my2DTaggingPtBinning,"2DTaggingPlan",nbinsx2D,xlow2D,xup2D,nbinsy2D,ylow2D,yup2D) ;

	vector<vector<TH2F*> > v2DTaggingPlan_FlavourPt_divided = buildFlavourPtVectorH2(my2DTaggingPtBinning,"2DTaggingPlan_divided",nbinsx2D,xlow2D,xup2D,nbinsy2D,ylow2D,yup2D) ;

	//flavour number per zone per flavour per pt
	vector<vector<vector<float> > > vFlavourNumber = build_fraction_ZoneFlavourPtVector(my2DTaggingPtBinning) ;

	//flavour number per zone per flavour per pt
	vector<vector<vector<float> > > vFlavourFraction = build_fraction_ZoneFlavourPtVector(my2DTaggingPtBinning) ;

	
	for(int j=0; j<n2DTaggingPtBins; j++) {
		vFractionHisto_Pt[j]->Sumw2();
		for(int i=0; i<nzones; i++) {
			vRmpf_ZonePt[i][j]->Sumw2();
			vRtrue_ZonePt[i][j]->Sumw2();
		}
		for(int k=0; k<nflavours; k++) {
			vRmpf_FlavourPt[k][j]->Sumw2();
			vRtrue_FlavourPt[k][j]->Sumw2();
			vRmpf_in2DTaggingZone_FlavourPt[k][j]->Sumw2();
			vRtrue_in2DTaggingZone_FlavourPt[k][j]->Sumw2();
			v2DTaggingPlan_FlavourPt[k][j]->Sumw2();
			vQGL_FlavourPt[k][j]->Sumw2();
			vCSV_FlavourPt[k][j]->Sumw2();
			for(int l=0; l<nzones; l++) {
				vRmpf_ZoneFlavourPt[l][k][j]->Sumw2();
				vRtrue_ZoneFlavourPt[l][k][j]->Sumw2();
			}
		}
	}
	
	for(int k=0; k<nflavours; k++) {
		vGammapt_Flavour[k]->Sumw2();
	}
	

	TH1F* hGammaPt=new TH1F("hGammaPt","hGammaPt",30,0,800);
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
	
	TH1F* hAlpha=new TH1F("hAlpha","hAlpha",nbinsx,xlow,xup);
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
	TTree* t_firstjetgen=(TTree*)f->Get("first_jet_gen");
	TTree* t_secondjetgen=(TTree*)f->Get("second_jet_gen");
	TTree* t_secondjet=(TTree*)f->Get("second_jet");
	TTree* t_electron=(TTree*)f->Get("electrons");
	TTree* t_muon=(TTree*)f->Get("muons");
	TTree* t_met=(TTree*)f->Get("met");
	TTree* t_gamma=(TTree*)f->Get("photon");
	TTree* t_gammagen=(TTree*)f->Get("photon_gen");
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
	
	//First jet generated
	gROOT->ProcessLine("#include <vector>"); 
	float firstjetgenpt;
	t_firstjetgen->SetBranchAddress("pt",&firstjetgenpt);
	float firstjetgenpx;
	t_firstjetgen->SetBranchAddress("px",&firstjetgenpx);
	float firstjetgenpy;
	t_firstjetgen->SetBranchAddress("py",&firstjetgenpy);
	float firstjetgenphi;
	t_firstjetgen->SetBranchAddress("phi",&firstjetgenphi);
	float firstjetgeneta;
	t_firstjetgen->SetBranchAddress("eta",&firstjetgeneta);
	int firstjetgenpdgid;
	t_firstjetgen->SetBranchAddress("parton_pdg_id",&firstjetgenpdgid);
// 	t_firstjetgen->SetBranchAddress("parton_flavour",&firstjetgenpdgid);
// 	TClonesArray *aneutrino_4vect = 0;
// 	t_firstjetgen->SetBranchAddress("neutrinos",&aneutrino_4vect);

	
	//Second jet generated
	float secondjetgenpt;
	t_secondjetgen->SetBranchAddress("pt",&secondjetgenpt);
	float secondjetgenphi;
	t_secondjetgen->SetBranchAddress("phi",&secondjetgenphi);
	float secondjetgeneta;
	t_secondjetgen->SetBranchAddress("eta",&secondjetgeneta);
	
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
	
	//Gamma gen
	float gammagenpt;
	t_gammagen->SetBranchAddress("pt",&gammagenpt);
	float gammagenpx;
	t_gammagen->SetBranchAddress("px",&gammagenpx);
	float gammagenpy;
	t_gammagen->SetBranchAddress("py",&gammagenpy);

	
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
	float Rtrue;
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
	Double_t DeltaPhi_j1j2_gen;
	Double_t DeltaEta_j1j2_gen;
	Double_t DeltaR_j1j2_gen;

	int binPt;//bin en pt
	int bin2DTaggingPt;//gros bins en pt pour l etude 2D tagging
	int binZone;// 2D tagging zone bin
	int binFlavour;
	
	bool dropEvent=false;
	
	//count events in the tree
	unsigned int nEvents = (int)t_firstjet->GetEntries();
	
	string type;
	
	//loop over them
	for(unsigned int ievt=0; ievt<nEvents; ievt++) {
		t_firstjet->GetEntry(ievt);
		t_firstjetraw->GetEntry(ievt);
		t_secondjet->GetEntry(ievt);
		t_electron->GetEntry(ievt);
		t_muon->GetEntry(ievt);
		t_gamma->GetEntry(ievt);
		t_gammagen->GetEntry(ievt);
		t_met->GetEntry(ievt);
		t_misc->GetEntry(ievt);
		t_firstjetgen->GetEntry(ievt);
		t_secondjetgen->GetEntry(ievt); 
		t_rho->GetEntry(ievt);   
		//t_gen_particles->GetEntry(ievt);
	
//*****************************************************************************************************
		
		dropEvent=false;
		
		if(gamma_has_pixel_seed == true) continue;
		  	
		binPt = myPtBinning.getPtBin(gammapt);
		if(binPt == -1) continue;

		bin2DTaggingPt = my2DTaggingPtBinning.getPtBin(gammapt);		

		if(TMath::Abs(firstjetgenpdgid) == 21) {
			type = "gluon";
		}
		else if(TMath::Abs(firstjetgenpdgid) == 4) {//c
			type = "quark";
		}
		else if(TMath::Abs(firstjetgenpdgid) == 5) {//b
			type = "quark";
		}
		else if(TMath::Abs(firstjetgenpdgid) == 1 || TMath::Abs(firstjetgenpdgid) == 2 || TMath::Abs(firstjetgenpdgid) == 3) {
			type = "quark";
		}
		else {
			type = "gluon";
		}
		
		//firstjetqg_tag_likelihood = qgsyst.Smear(firstjetpt, firstjeteta, rho, firstjetqg_tag_likelihood, type);

		binZone = getZoneBin(firstjetbtag_csv, firstjetqg_tag_likelihood);
		//if(binZone == -1) continue;	
		binFlavour = getFlavourBin(firstjetgenpdgid);
		//if(binFlavour == -1) continue;
		
		Rmpf = 1 + (gammapx*metpx + gammapy*metpy)/(pow(gammapt,2));		
		alpha = (secondjetpt)/(gammapt);
		Rtrue = (firstjetpt)/(firstjetgenpt);

		
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
    	
//DeltaPhi between the first jetgen and the second jetgen calculation
    	DeltaPhi_j1j2_gen = TMath::Abs((secondjetgenphi) - (firstjetgenphi));
    	if(DeltaPhi_j1j2_gen>TMath::Pi()){
    	  DeltaPhi_j1j2_gen = 2*TMath::Pi()-DeltaPhi_j1j2_gen;
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
    	
 //DeltaEta between the first jetgen and the second jetgen calculation
    	DeltaEta_j1j2_gen = TMath::Abs((secondjetgeneta) - (firstjetgeneta));
    	    	
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
    	
//DeltaR between the first jetgen and the second jetgen calculation
    	DeltaR_j1j2_gen = sqrt ( pow(DeltaEta_j1j2_gen,2) + pow(DeltaPhi_j1j2_gen,2) );
    	
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

	hAlpha->Fill(alpha,miscevent_weight);
	
		if(DeltaPhi_j1gamma>2.8 && fabs(firstjeteta)<1.3 && (alpha<0.2 || secondjetpt<10)) {
			if(muonn==0){
				if(electronn==0) {
					vGammapt_Flavour[binFlavour]->Fill(gammapt,miscevent_weight);	
					vGammapt_Flavour[nflavours-1]->Fill(gammapt,miscevent_weight);
	
	
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
						if(binFlavour != 4) {
							vQGL_FlavourPt[binFlavour][bin2DTaggingPt]->Fill(firstjetqg_tag_likelihood,miscevent_weight);
							vCSV_FlavourPt[binFlavour][bin2DTaggingPt]->Fill(firstjetbtag_csv,miscevent_weight);
							v2DTaggingPlan_FlavourPt[binFlavour][bin2DTaggingPt]->Fill(firstjetbtag_csv,firstjetqg_tag_likelihood,miscevent_weight);
							vQGL_FlavourPt[nflavours-1][bin2DTaggingPt]->Fill(firstjetqg_tag_likelihood,miscevent_weight);
							vCSV_FlavourPt[nflavours-1][bin2DTaggingPt]->Fill(firstjetbtag_csv,miscevent_weight);
							v2DTaggingPlan_FlavourPt[nflavours-1][bin2DTaggingPt]->Fill(firstjetbtag_csv,firstjetqg_tag_likelihood,miscevent_weight);//nflavours-1 : corresponds to total 
																							//i.e. Ntot = Nuds + Nc + Nb + Ng

							vRmpf_FlavourPt[binFlavour][bin2DTaggingPt]->Fill(Rmpf,miscevent_weight);
							vRtrue_FlavourPt[binFlavour][bin2DTaggingPt]->Fill(Rtrue,miscevent_weight);
						}

						if(binZone!=-1 && binFlavour != 4) {
							vRmpf_ZonePt[binZone][bin2DTaggingPt]->Fill(Rmpf,miscevent_weight);
							vRtrue_ZonePt[binZone][bin2DTaggingPt]->Fill(Rtrue,miscevent_weight);
							vRmpf_ZoneFlavourPt[binZone][binFlavour][bin2DTaggingPt]->Fill(Rmpf,miscevent_weight);
							vRtrue_ZoneFlavourPt[binZone][binFlavour][bin2DTaggingPt]->Fill(Rtrue,miscevent_weight);
							vRmpf_ZoneFlavourPt[binZone][nflavours-1][bin2DTaggingPt]->Fill(Rmpf,miscevent_weight);
							vRtrue_ZoneFlavourPt[binZone][nflavours-1][bin2DTaggingPt]->Fill(Rtrue,miscevent_weight);
							vRmpf_in2DTaggingZone_FlavourPt[binFlavour][bin2DTaggingPt]->Fill(Rmpf,miscevent_weight);
							vRtrue_in2DTaggingZone_FlavourPt[binFlavour][bin2DTaggingPt]->Fill(Rtrue,miscevent_weight);
							vFlavourNumber[binZone][binFlavour][bin2DTaggingPt] = vFlavourNumber[binZone][binFlavour][bin2DTaggingPt] + miscevent_weight;	
			
							vFlavourNumber[binZone][nflavours-1][bin2DTaggingPt] = vFlavourNumber[binZone][nflavours-1][bin2DTaggingPt] + miscevent_weight;
// 							vQGL_FlavourPt[binFlavour][bin2DTaggingPt]->Fill(firstjetqg_tag_likelihood,miscevent_weight);
// 							vCSV_FlavourPt[binFlavour][bin2DTaggingPt]->Fill(firstjetbtag_csv,miscevent_weight);
// 							v2DTaggingPlan_FlavourPt[binFlavour][bin2DTaggingPt]->Fill(firstjetbtag_csv,firstjetqg_tag_likelihood,miscevent_weight);
// 							vQGL_FlavourPt[nflavours-1][bin2DTaggingPt]->Fill(firstjetqg_tag_likelihood,miscevent_weight);
// 							vCSV_FlavourPt[nflavours-1][bin2DTaggingPt]->Fill(firstjetbtag_csv,miscevent_weight);
						}
					}	
				}
			}
		}
	}
    
	

	
//*****************************************************************************************************
//
//                                      2D Tagging Plans
//
//*****************************************************************************************************

//v2DTaggingPlan_FlavourPt[binFlavour][bin2DTaggingPt]
	std::string name2DTaggingPlan; 
	std::string name2DTaggingPlan_check; 
	std::string nameHisto;
	std::string nameHisto_check;
	TCanvas *cLikelihood_vs_csv = new TCanvas();
	cLikelihood_vs_csv->cd();
	string histoName;
	for(int j=0; j<n2DTaggingPtBins; j++) {
		for(int i=0; i<nflavours-1; i++) {//divide flavour distri by total distri
			cLikelihood_vs_csv->Clear();
			histoName = "images2DTagging/2DTaggingZones/beforeDividing/" + getFlavourBinName(i) + "_" + my2DTaggingPtBinning.getName(j) + "_stageM2.pdf";
			//name2DTaggingPlan_check = ("images2DTagging/2DTaggingZones/beforeDividing/");

			v2DTaggingPlan_FlavourPt[i][j]->Draw("colz");
			v2DTaggingPlan_FlavourPt[i][j]->SetStats(0);
			v2DTaggingPlan_FlavourPt[i][j]->GetXaxis()->SetRangeUser(0.,1.);
			v2DTaggingPlan_FlavourPt[i][j]->GetXaxis()->SetLabelOffset(0.005);
			v2DTaggingPlan_FlavourPt[i][j]->GetXaxis()->SetLabelFont(42);
			v2DTaggingPlan_FlavourPt[i][j]->GetXaxis()->SetLabelSize(0.045);
			v2DTaggingPlan_FlavourPt[i][j]->GetXaxis()->SetTitleOffset(1.15);
			v2DTaggingPlan_FlavourPt[i][j]->GetXaxis()->SetTitleSize(0.04);
			v2DTaggingPlan_FlavourPt[i][j]->GetXaxis()->SetTitleFont(42);
			v2DTaggingPlan_FlavourPt[i][j]->GetXaxis()->SetTitle("btag_csv");
			v2DTaggingPlan_FlavourPt[i][j]->GetYaxis()->SetRangeUser(0.,1.);
			v2DTaggingPlan_FlavourPt[i][j]->GetYaxis()->SetLabelOffset(0.005);
			v2DTaggingPlan_FlavourPt[i][j]->GetYaxis()->SetLabelFont(42);
			v2DTaggingPlan_FlavourPt[i][j]->GetYaxis()->SetLabelSize(0.045);	
			v2DTaggingPlan_FlavourPt[i][j]->GetYaxis()->SetTitleOffset(1.2);
			v2DTaggingPlan_FlavourPt[i][j]->GetYaxis()->SetTitleFont(42);
			v2DTaggingPlan_FlavourPt[i][j]->GetYaxis()->SetTitleSize(0.04);	
			v2DTaggingPlan_FlavourPt[i][j]->GetYaxis()->SetTitle("QG_likelihood");
			v2DTaggingPlan_FlavourPt[i][j]->GetZaxis()->SetLabelSize(0.03);
			//nameHisto_check = v2DTaggingPlan_FlavourPt[i][j]->GetTitle();
			//name2DTaggingPlan_check += nameHisto;
			//name2DTaggingPlan_check += "_stageM2.pdf";
			//cLikelihood_vs_csv->SaveAs(name2DTaggingPlan_check.c_str());
			cLikelihood_vs_csv->SaveAs(histoName.c_str());

			v2DTaggingPlan_FlavourPt_divided[i][j] = (TH2F*)v2DTaggingPlan_FlavourPt[i][j]->Clone(v2DTaggingPlan_FlavourPt[i][j]->GetTitle());

			v2DTaggingPlan_FlavourPt_divided[i][j]->Divide(v2DTaggingPlan_FlavourPt[nflavours-1][j]);

			name2DTaggingPlan = ("images2DTagging/2DTaggingZones/");

			cLikelihood_vs_csv->Clear();
	
			v2DTaggingPlan_FlavourPt_divided[i][j]->Draw("colz");
			v2DTaggingPlan_FlavourPt_divided[i][j]->SetStats(0);
			v2DTaggingPlan_FlavourPt_divided[i][j]->GetXaxis()->SetRangeUser(0.,1.);
			v2DTaggingPlan_FlavourPt_divided[i][j]->GetXaxis()->SetLabelOffset(0.005);
			v2DTaggingPlan_FlavourPt_divided[i][j]->GetXaxis()->SetLabelFont(42);
			v2DTaggingPlan_FlavourPt_divided[i][j]->GetXaxis()->SetLabelSize(0.045);
			v2DTaggingPlan_FlavourPt_divided[i][j]->GetXaxis()->SetTitleOffset(1.15);
			v2DTaggingPlan_FlavourPt_divided[i][j]->GetXaxis()->SetTitleSize(0.04);
			v2DTaggingPlan_FlavourPt_divided[i][j]->GetXaxis()->SetTitleFont(42);
			v2DTaggingPlan_FlavourPt_divided[i][j]->GetXaxis()->SetTitle("btag_csv");
			v2DTaggingPlan_FlavourPt_divided[i][j]->GetYaxis()->SetRangeUser(0.,1.);
			v2DTaggingPlan_FlavourPt_divided[i][j]->GetYaxis()->SetLabelOffset(0.005);
			v2DTaggingPlan_FlavourPt_divided[i][j]->GetYaxis()->SetLabelFont(42);
			v2DTaggingPlan_FlavourPt_divided[i][j]->GetYaxis()->SetLabelSize(0.045);	
			v2DTaggingPlan_FlavourPt_divided[i][j]->GetYaxis()->SetTitleOffset(1.2);
			v2DTaggingPlan_FlavourPt_divided[i][j]->GetYaxis()->SetTitleFont(42);
			v2DTaggingPlan_FlavourPt_divided[i][j]->GetYaxis()->SetTitleSize(0.04);	
			v2DTaggingPlan_FlavourPt_divided[i][j]->GetYaxis()->SetTitle("QG_likelihood");
			v2DTaggingPlan_FlavourPt_divided[i][j]->GetZaxis()->SetLabelSize(0.03);
			nameHisto = v2DTaggingPlan_FlavourPt[i][j]->GetTitle();
			name2DTaggingPlan += nameHisto;
			name2DTaggingPlan += "_stageM2.pdf";
			cLikelihood_vs_csv->SaveAs(name2DTaggingPlan.c_str());
			nameHisto = "";
			name2DTaggingPlan = "";
// 			nameHisto_check = "";
// 			name2DTaggingPlan_check = "";
			histoName = "";
				}
	}





//*****************************************************************************************************
//
//                                      Jet flavour fraction
//
//*****************************************************************************************************

	float Nuds,Ng,Nc,Nb,NnoMatched,Ntot;
	std::string nameFraction;
	std::string nameHistoFraction;
	
	for(int j=0;j<n2DTaggingPtBins; j++) {
		for(int i=0; i<nzones; i++) {
			nameFraction = ("images2DTagging/FlavourFractions/");
			Nuds       = vFlavourNumber[i][0][j];
			Ng         = vFlavourNumber[i][1][j];
			Nc         = vFlavourNumber[i][2][j];
			Nb         = vFlavourNumber[i][3][j];
			NnoMatched = vFlavourNumber[i][4][j];
			Ntot       = vFlavourNumber[i][5][j];
			//Ntot       = vFlavourNumber[i][0][j]+vFlavourNumber[i][1][j]+vFlavourNumber[i][2][j]+vFlavourNumber[i][3][j];
			nameHistoFraction = vFractionHisto_Pt[j]->GetName();
			nameFraction += nameHistoFraction ;
			nameFraction += "_stageM2.pdf";
			computeFlavourRatio(Nuds, Ng, Nc, Nb, NnoMatched, Ntot,vFractionHisto_Pt[j], nameFraction, my2DTaggingPtBinning, i, j);
			nameFraction = "";
			nameHistoFraction = "";
		}
	} 

	for(int k=0; k<n2DTaggingPtBins; k++) {
		for(int j=0; j<nflavours ; j++) {
			for(int i=0; i<nzones; i++) {
				vFlavourFraction[i][j][k] = vFlavourNumber[i][j][k];
			}
		}
	}

	cout<<""<<endl<<endl;

	for(int k=0; k<n2DTaggingPtBins; k++) {
		for(int j=0; j<nflavours ; j++) {
			for(int i=0; i<nzones; i++) {
				//Ntot       = vFlavourNumber[i][0][k]+vFlavourNumber[i][1][k]+vFlavourNumber[i][2][k]+vFlavourNumber[i][3][k];
				//vFlavourFraction[i][j][k] = vFlavourFraction[i][j][k]/Ntot;
				vFlavourFraction[i][j][k] = vFlavourFraction[i][j][k]/vFlavourFraction[i][nflavours-1][k];
				cout<<"vFlavourFraction["<<i<<"]["<<j<<"]["<<k<<"] : "<<vFlavourFraction[i][j][k]*100.<<" %"<<endl;
			}
		}
	}

//****************************************************************************************************
//
//                                      filling the matrix with the flavour fractions per zone
//
//****************************************************************************************************

	cout<<""<<endl<<endl;

	//filling the matrices
	for(int k=0; k<n2DTaggingPtBins; k++) {
		for(int j=0; j<nflavours-2; j++) {
			for(int i=0; i<nzones-2; i++) {
				v4x4MatrixPt[k](i,j) = vFlavourFraction[i][j][k];
				cout<<"v4x4MatrixPt["<<k<<"]("<<i<<","<<j<<") :  "<<v4x4MatrixPt[k](i,j)*100.<<" %"<<endl;
			}
		}
	}

	cout<<""<<endl<<endl;

	for(int k=0; k<n2DTaggingPtBins; k++) {
		for(int j=0; j<nflavours-2; j++) {
			for(int i=0; i<nzones; i++) {
				v6x4MatrixPt[k](i,j) = vFlavourFraction[i][j][k];
				cout<<"v6x4MatrixPt["<<k<<"]("<<i<<","<<j<<") : "<<v6x4MatrixPt[k](i,j)*100.<<" %"<<endl;
			}
		}
	}

	TFile *out_matrix = new TFile("output_rootfile/outputMatrix2DTagging_MC_G_stageM2.root", "recreate");
	out_matrix->cd();
	TDirectory *matrixDir = out_matrix->mkdir("matrix","matrix");
	matrixDir->cd();
	for(int k=0; k<n2DTaggingPtBins; k++) {
		string matrix6x4Name;
		string matrix4x4Name;
		matrix6x4Name = "matrix6x4_" + my2DTaggingPtBinning.getName(k);
		matrix4x4Name = "matrix4x4_" + my2DTaggingPtBinning.getName(k);
		v6x4MatrixPt[k].Write(matrix6x4Name.c_str());
		v4x4MatrixPt[k].Write(matrix4x4Name.c_str());
	}
	out_matrix->Close();

//*****************************************************************************************************
//
//                                      Output file
//
//*****************************************************************************************************

	//create the output file and write into it
	TFile *out = new TFile("output_rootfile/output2DTagging_MC_G_stageM2.root", "recreate");
	
	out->cd();	
	TDirectory *response_Zone_PtDir = out->mkdir("response_Zone_Pt","response_Zone_Pt");
	TDirectory *Rmpf_Zone_PtDir = response_Zone_PtDir->mkdir("Rmpf","Rmpf");
	Rmpf_Zone_PtDir->cd();
	for(int j=0; j<n2DTaggingPtBins; j++) {
		for(int i=0; i<nzones; i++) {
			vRmpf_ZonePt[i][j]->Write();
		}
	}
	TDirectory *Rtrue_Zone_PtDir = response_Zone_PtDir->mkdir("Rtrue","Rtrue");
	Rtrue_Zone_PtDir->cd();
	for(int j=0; j<n2DTaggingPtBins; j++) {
		for(int i=0; i<nzones; i++) {
			vRtrue_ZonePt[i][j]->Write();
		}
	}


	TDirectory *response_Flavour_PtDir = out->mkdir("response_Flavour_Pt","response_Flavour_Pt");
	TDirectory *Rmpf_Flavour_PtDir = response_Flavour_PtDir->mkdir("Rmpf","Rmpf");
	Rmpf_Flavour_PtDir->cd();
	for(int j=0; j<n2DTaggingPtBins; j++) {
		for(int k=0; k<nflavours; k++) {
			vRmpf_FlavourPt[k][j]->Write();
			vRmpf_in2DTaggingZone_FlavourPt[k][j]->Write();
		}
	}
	TDirectory *Rtrue_Flavour_PtDir = response_Flavour_PtDir->mkdir("Rtrue","Rtrue");
	Rtrue_Flavour_PtDir->cd();
	for(int j=0; j<n2DTaggingPtBins; j++) {
		for(int k=0; k<nflavours; k++) {
			vRtrue_FlavourPt[k][j]->Write();
			vRtrue_in2DTaggingZone_FlavourPt[k][j]->Write();
		}
	}

	TDirectory *response_Zone_Flavour_PtDir = out->mkdir("response_Zone_Flavour_Pt","response_Zone_Flavour_Pt");
	TDirectory *Rmpf_Zone_Flavour_PtDir = response_Zone_Flavour_PtDir->mkdir("Rmpf","Rmpf");
	Rmpf_Zone_Flavour_PtDir->cd();
	for(int j=0; j<n2DTaggingPtBins; j++) {
		for(int k=0; k<nflavours; k++) {
			for(int l=0; l<nzones; l++) {
				vRmpf_ZoneFlavourPt[l][k][j]->Write();
			}
		}
	}
	TDirectory *Rtrue_Zone_Flavour_PtDir = response_Zone_Flavour_PtDir->mkdir("Rtrue","Rtrue");
	Rtrue_Zone_Flavour_PtDir->cd();
	for(int j=0; j<n2DTaggingPtBins; j++) {
		for(int k=0; k<nflavours; k++) {
			for(int l=0; l<nzones; l++) {
				vRtrue_ZoneFlavourPt[l][k][j]->Write();
			}
		}
	}

	TDirectory *tagger_Flavour_PtDir = out->mkdir("tagger_Flavour_Pt","tagger_Flavour_Pt");
	TDirectory *CSV_Flavour_PtDir = tagger_Flavour_PtDir->mkdir("CSV","CSV");
	CSV_Flavour_PtDir->cd();
	for(int j=0; j<n2DTaggingPtBins; j++) {
		for(int k=0; k<nflavours; k++) {
			vCSV_FlavourPt[k][j]->Write();
		}
	}
	TDirectory *QGL_Flavour_PtDir = tagger_Flavour_PtDir->mkdir("QGL","QGL");
	QGL_Flavour_PtDir->cd();
	for(int j=0; j<n2DTaggingPtBins; j++) {
		for(int k=0; k<nflavours; k++) {
			vQGL_FlavourPt[k][j]->Write();
		}
	}
	
	TDirectory *gammapt_FlavourDir = out->mkdir("gammapt_Flavour","gammapt_Flavour");
	gammapt_FlavourDir->cd();
	for(int k=0; k<nflavours; k++) {
		vGammapt_Flavour[k]->Write();
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

	TDirectory *plan2DTaggingDir = out->mkdir("plan2DTagging","plan2DTagging");
	plan2DTaggingDir->cd();
	for(int j=0; j<n2DTaggingPtBins; j++) {
		for(int i=0; i<nflavours; i++) {//divide flavour distri by total distri
			v2DTaggingPlan_FlavourPt[i][j]->Write();
			v2DTaggingPlan_FlavourPt_divided[i][j]->Write();	
		}
	}
	
	out->Close();
	
	return 0;
}









