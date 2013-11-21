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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <stdio.h>
#include <TString.h>
#include "TMatrixF.h"
#include "TMatrixFSym.h"
#include "TMatrixFSparse.h"
#include "TMatrixFLazy.h"

#include "../../binning.h"
#include "../../common.h"

#include "../../ptBinning.h"
#include "../../alphaBinning.h"
#include "../../fitTools.h"



int getNrows() {
	//size of the square matrix
	//return getZoneNumber() - 2;//we do not use the mix zones
	return 4;//we do not use the mix zones
}

int getNcolumns() {
	//size of the square matrix
	//return getFlavourNumber() - 2; //number of flavours in the matrix
	return 4; //number of flavours in the matrix
}

void TGraph_flavourLabel(TGraph* g) {
	g->GetXaxis()->SetNdivisions(4,kTRUE);
	for(int i=0; i<getFlavourNumber()-2; i++) {
		g->GetXaxis()->SetBinLabel((g->GetXaxis()->FindBin(i+1)),getFlavourBinName(i).c_str());
	}
	g->GetXaxis()->LabelsOption("h");
	g->GetXaxis()->SetLabelSize(0.055);
}

vector<TH1F*> getResponse_Zone_Pt(TFile* f, string aPtName) {
	vector<TH1F*> vHisto;
	for(int i=0; i<getZoneNumber(); i++){
		string histoName;
		histoName = "response_Zone_Pt/Rmpf/Rmpf_" + getZoneBinName(i) + "_" + aPtName;
		//cout<<histoNameMC<<endl;
		vHisto.push_back((TH1F*)f->Get(histoName.c_str()));
	}
	return vHisto;
}

vector<TH1F*> getResponse_Flavour_Pt(TFile* f, string aPtName) {
//*******************************************************************************************************
// description: function which gets all the flavour responses for a giving 2DTagging pt bin 
//              and stores these histo in a vector<TH1F*> which is returned
// arguments: - f: MC rootfile containning informations we want to extract
//            - aPtName: name of the 2DTagging pt bin under study
//*******************************************************************************************************
	vector<TH1F*> vHisto;
	for(int i=0; i<getFlavourNumber(); i++){
		string histoName;
		histoName = "response_Flavour_Pt/Rmpf/Rmpf_in2DTaggingZone_" + getFlavourBinName(i) + "_" + aPtName;
		//cout<<histoNameMC<<endl;
		vHisto.push_back((TH1F*)f->Get(histoName.c_str()));
	}
	return vHisto;
}

TMatrixF buildSquareMatrix() {
// 	int ncolumn = getNcolumns();
// 	int nrow = getNrows();
// 	TMatrixF m(nrow,ncolumn);
	TMatrixF m(4,4);
	return m;
}



//************************************************************************************************
// Functions used in the construction of the flavour fraction fluctuation matrix
//************************************************************************************************

int getN_negFraction(const TMatrixF& myMatrix, int row) {
	int n=0;
	for(int j=0; j<4; j++) {
		if( myMatrix(row,j) < 0.) {
			n=n+1;
		}
	}
	return n;
}

int getN_nullFraction(const TMatrixF& myMatrix, int row) {
	int n=0;
	for(int j=0; j<4; j++) {
		if(myMatrix(row,j) == 0.) {
			n=n+1;
		}
	}
	return n;
}

int getN_sup1Fraction(const TMatrixF& myMatrix, int row) {
	int n=0;
	for(int j=0; j<4; j++) {
		if(myMatrix(row,j) > 1.) {
			n=n+1;
		}
	}
	return n;
}

int getIndex_negFraction(const TMatrixF& myMatrix, int row) {
	int index=0;
	for(int j=0; j<4; j++) {
		if(myMatrix(row,j) < 0.) {
			index = j;
		}
	}
	return index;
}

int getIndex_nullFraction(const TMatrixF& myMatrix, int row) {
	int index=0;
	for(int j=0; j<4; j++) {
		if(myMatrix(row,j) == 0.) {
			index = j;
		}
	}
	return index;
}

int getIndex_posFraction(const TMatrixF& myMatrix, int row, int index_max) {
	int index=0;
	for(int j=0; j<4; j++) {
		if(j != index_max) {
			if(myMatrix(row,j) > 0.) {
				index = j;
			}
		}
	}
	return index;
}

int getIndex_maxFraction(const TMatrixF& myMatrix, int row) {
	int index=0;
	float maxFraction_tmp = myMatrix(row,0);
	for(int j=0; j<4; j++) {
		if(myMatrix(row,j) > maxFraction_tmp) {
			maxFraction_tmp = myMatrix(row,j);
			index = j;
		}
	}
	return index;
}

void setMaxFracTo1(TMatrixF& myMatrix, int row) {
	int index_tmp = 0;
	for(int j=0; j<4; j++) {
		if(myMatrix(row,j)>1) {
			index_tmp = j;
			myMatrix(row,j) = 1.;
		}
	}
	for(int j=0; j<4; j++) {
		if(j != index_tmp) {
			myMatrix(row,j) = 0.;
		}
	}
}


TMatrixF getMatrix(TFile* f_matrix, string aPtName) {
//*******************************************************************************************************
// description: function which catches up the 4x4 matrix containing the flavour fractions per 2DT zone
// arguments: - f_matrix: MC rootfile containning the flavour matrix we want to retrieve
//            - aPtName: name of the 2DTagging pt bin under study
//*******************************************************************************************************
	TMatrixF m = buildSquareMatrix();
	string matrixName;
	matrixName = "matrix/matrix4x4_" + aPtName;
	m = *((TMatrixF*) f_matrix->Get(matrixName.c_str()));
	return m;
}

TMatrixF buildColumnMatrix() {
	int ncolumn = 1; //number of flavours in the matrix
	int nrow = getNrows();
	TMatrixF m(nrow,ncolumn);
	return m;
}

vector<TMatrixF> buildErrorInverseMatrix(const TMatrixF& m0) {
// build systematique error matrix from the raw flavour matrix
// i.e. construct the 8 matrices with +/- 30% error on flavour fractions
// and return the 8 corresponding inverse matrices

//********************************************************************************************

	vector<TMatrixF> vMatrix;
	vector<TMatrixF> vMatrix_inverse;
	cout<<"creating matrices..."<<endl;
	//for(int i=0; i<2*getNcolumns(); i++) {
	for(int i=0; i<8; i++) {
		vMatrix.push_back(buildSquareMatrix());
		vMatrix_inverse.push_back(buildSquareMatrix());
	}
	cout<<"creation done..."<<endl;

	for(int k=0; k<4; k++) {
		for(int i = vMatrix[k].GetRowLwb(); i <= vMatrix[k].GetRowUpb(); i++) {
			for(int j = vMatrix[k].GetColLwb(); j <= vMatrix[k].GetColUpb(); j++) {
				if(j==k) {
					vMatrix[k](i,j) = m0(i,j)*1.3;
				}
				else {
					vMatrix[k](i,j) = m0(i,j)-m0(i,k)*0.1;
				}
			}
		}
	}
	
	for(int k=4; k<8; k++) {
		for(int i = vMatrix[k].GetRowLwb(); i <= vMatrix[k].GetRowUpb(); i++) {
			for(int j = vMatrix[k].GetColLwb(); j <= vMatrix[k].GetColUpb(); j++) {
				if(j==(k-4)) {
					vMatrix[k](i,j) = m0(i,j)*0.7;
				}
				else {
					vMatrix[k](i,j) = m0(i,j)+m0(i,k-4)*0.1;
				}
			}
		}
	}

	for(int k=0; k<2*getNcolumns(); k++) {		
		cout<<"verification matrice vMatrix["<<k<<"]..."<<endl;
		for (int i = vMatrix[k].GetRowLwb(); i <= vMatrix[k].GetRowUpb(); i++) {
			for (int j = vMatrix[k].GetColLwb(); j <= vMatrix[k].GetColUpb(); j++) {
				cout<<"vMatrix["<<k<<"]("<<i<<","<<j<<") = "<<vMatrix[k](i,j)<<endl;
			}
		}
	}


	cout<<"verification matrice terminee..."<<endl;
	
	for(int k=0; k<2*getNcolumns(); k++) {	
		for(int i=0; i<4; i++) {
			cout<<"matrix "<<k<<endl;
			cout<<"ligne "<<i<<endl;

			int n_negFraction = 0;
			int n_sup1Fraction = 0;
				
			n_negFraction = getN_negFraction(vMatrix[k], i);
			n_sup1Fraction = getN_sup1Fraction(vMatrix[k], i);

			cout<<"n_negFraction: "<< n_negFraction <<endl;
			cout<<"n_sup1Fraction: "<< n_sup1Fraction <<endl;
			for(int j=0; j<4; j++) {
				cout<<"vMatrix["<<k<<"]("<<i<<","<<j<<") = "<<vMatrix[k](i,j)<<endl;
			}

			while(n_negFraction > 0 || n_sup1Fraction > 0) {
				cout<<"In the while condition..."<<endl;
				int need_to_check_neg_frac = 1;
			//first check if there is no fraction > 1.
			//if there are, then you don't need to check the condition 'fraction>0'
				if(n_sup1Fraction > 0) {
					need_to_check_neg_frac = 0;
					setMaxFracTo1(vMatrix[k], i);
				}
	
				if(need_to_check_neg_frac == 1) {
					//first find the higher fraction
					int index_maxFraction = getIndex_maxFraction(vMatrix[k], i);
			
					//if n_negFraction = 1
					if(n_negFraction == 1) {
						int index_negFraction = getIndex_negFraction(vMatrix[k], i);
						int n_nullFraction = getN_nullFraction(vMatrix[k], i);
						if(n_nullFraction == 0) {
							float value = fabs(vMatrix[k](i,index_negFraction)/2.);
							for(int j=0; j<4; j++) {
								if(j!=index_maxFraction && j!=index_negFraction) {
									vMatrix[k](i,j) = vMatrix[k](i,j) - value;
								}
								if(j==index_negFraction) {
									vMatrix[k](i,j) = 0.;
								}								
							}
						}
						else if(n_nullFraction == 1) {
							int index_nullFraction = getIndex_nullFraction(vMatrix[k], i);
							float value = fabs(vMatrix[k](i,index_negFraction));
							for(int j=0; j<4; j++) {
								if(j!=index_maxFraction && j!=index_negFraction && j!=index_nullFraction) {
									vMatrix[k](i,j) = vMatrix[k](i,j) - value;
								}
								if(j==index_negFraction) {
									vMatrix[k](i,j) = 0.;
								}		
							}
						}
						else if(n_nullFraction == 2) {
							setMaxFracTo1(vMatrix[k], i);
						}
					}
			
					
					//if n_negFraction = 2
					else if(n_negFraction == 2) {
						int n_nullFraction = getN_nullFraction(vMatrix[k], i);
						if(n_nullFraction == 0) {
							int index_negFraction1 = getIndex_negFraction(vMatrix[k], i);
							int index_negFraction2 = 0;
							for(int j=0; j<4; j++) {
								if(vMatrix[k](i,j)<0. && j!=index_negFraction1) {
									index_negFraction2 = j;
								}
							}	
		
							int index_posFraction = getIndex_posFraction(vMatrix[k], i, index_maxFraction);

							cout<<"index_negFraction1: "<<index_negFraction1<<endl;	
							cout<<"index_negFraction2: "<<index_negFraction2<<endl;	
							cout<<"index_posFraction: "<<index_posFraction<<endl;	
							cout<<"index_maxFraction: "<<index_maxFraction<<endl;	

							vMatrix[k](i,index_posFraction) = vMatrix[k](i,index_posFraction) - fabs(vMatrix[k](i,index_negFraction1)) - fabs(vMatrix[k](i,index_negFraction2));
							vMatrix[k](i,index_negFraction1) = 0.;
							vMatrix[k](i,index_negFraction2) = 0.;

// 							for(int j=0; j<4; j++) {
// 								if(j!=index_maxFraction) {
// 									if(vMatrix[k](i,j)>0.) {
// 									vMatrix[k](i,j) = vMatrix[k](i,j) - fabs(vMatrix[k](i,index_negFraction1)) - fabs(vMatrix[k](i,index_negFraction2));
// 									}
// 								}
// 							}

						}
						else if(n_nullFraction == 1) {
							setMaxFracTo1(vMatrix[k], i);
						}
					}
			
					//if n_negFraction = 3
					else if(n_negFraction == 3) {
						setMaxFracTo1(vMatrix[k], i);
					}
				}

								
				for(int j=0; j<4; j++) {
					cout<<"vMatrix["<<k<<"]("<<i<<","<<j<<") = "<<vMatrix[k](i,j)<<endl;
				}

				n_sup1Fraction = getN_sup1Fraction(vMatrix[k], i);		
				n_negFraction = getN_negFraction(vMatrix[k], i);

				cout<<"n_negFraction: "<< n_negFraction <<endl;
				cout<<"n_sup1Fraction: "<< n_sup1Fraction <<endl;
			}
		}
	}



	for(int k=0; k<2*getNcolumns(); k++) {
		vMatrix_inverse[k] = ((TMatrixF*)vMatrix[k].Clone())->Invert();	
		for (int i = vMatrix_inverse[k].GetRowLwb(); i <= vMatrix_inverse[k].GetRowUpb(); i++) {
			for (int j = vMatrix_inverse[k].GetColLwb(); j <= vMatrix_inverse[k].GetColUpb(); j++) {
				cout<<"vMatrix_inverse["<<k<<"]("<<i<<","<<j<<") = "<<vMatrix_inverse[k](i,j)<<endl;
			}
		}		
		cout<<"verification matrice vMatrix["<<k<<"] inchangee..."<<endl;
		for (int i = vMatrix[k].GetRowLwb(); i <= vMatrix[k].GetRowUpb(); i++) {
			for (int j = vMatrix[k].GetColLwb(); j <= vMatrix[k].GetColUpb(); j++) {
				cout<<"vMatrix["<<k<<"]("<<i<<","<<j<<") = "<<vMatrix[k](i,j)<<endl;
			}
		}
	}
	
	return 	vMatrix_inverse;



//********************************************************************************************
/*
	//fraction of partons in each zone
	//zone 1
	float fuds1;
	float fb1;
	float fc1;
	float fg1;
	
	//zone 2
	float fuds2;
	float fb2;
	float fc2;
	float fg2;
	
	//zone 3
	float fuds3;
	float fb3;
	float fc3;
	float fg3;
	
	//zone 4
	float fuds4;
	float fb4;
	float fc4;
	float fg4;

	fuds1 = m0(0,0);
	fc1 = m0(0,1);
	fb1 = m0(0,2);
	fg1 = m0(0,3);
	fuds2 = m0(1,0);
	fc2 = m0(1,1);
	fb2 = m0(1,2);
	fg2 = m0(1,3);
	fuds3 = m0(2,0);
	fc3 = m0(2,1);
	fb3 = m0(2,2);
	fg3 = m0(2,3);
	fuds4 = m0(3,0);
	fc4 = m0(3,1);
	fb4 = m0(3,2);
	fg4 = m0(3,3);
	
	cout<<"fuds1 : "<<fuds1<<endl;
	cout<<"fc1 : "<<fc1<<endl;
	cout<<"fb1 : "<<fb1<<endl;
	cout<<"fg1 : "<<fg1<<endl;
	cout<<"fuds2 : "<<fuds2<<endl;
	cout<<"fc2 : "<<fc2<<endl;
	cout<<"fb2 : "<<fb2<<endl;
	cout<<"fg2 : "<<fg2<<endl;
	cout<<"fuds3 : "<<fuds3<<endl;
	cout<<"fc3 : "<<fc3<<endl;
	cout<<"fb3 : "<<fb3<<endl;
	cout<<"fg3 : "<<fg3<<endl;
	cout<<"fuds4 : "<<fuds4<<endl;
	cout<<"fc4 : "<<fc4<<endl;
	cout<<"fb4 : "<<fb4<<endl;
	cout<<"fg4 : "<<fg4<<endl;

	vector<TMatrixF> vMatrix;
	vector<TMatrixF> vMatrix_inverse;
	for(int i=0; i<2*getNcolumns(); i++) {
		vMatrix.push_back(buildSquareMatrix());
		vMatrix_inverse.push_back(buildSquareMatrix());
	}

	vMatrix[0](0,0) = fuds1*1.3;
// 	vMatrix[0](0,1) = 0.;
// 	vMatrix[0](0,2) = 0.;
// 	vMatrix[0](0,3) = fg1-(fuds1*0.3)/3.+(fc1-(fuds1*0.3)/3.)+(fb1-(fuds1*0.3)/3.);
	vMatrix[0](0,1) = fc1-(fuds1*0.3)/3.;
	vMatrix[0](0,2) = fb1-(fuds1*0.3)/3.;
	vMatrix[0](0,3) = fg1-(fuds1*0.3)/3.;
// 	vMatrix[0](0,0) = 1;
// 	vMatrix[0](0,1) = 0;
// 	vMatrix[0](0,2) = 0;
// 	vMatrix[0](0,3) = 0;
	vMatrix[0](1,0) = fuds2*1.3;
	vMatrix[0](1,1) = fc2-(fuds2*0.3)/3.;
	vMatrix[0](1,2) = fb2-(fuds2*0.3)/3.;
	vMatrix[0](1,3) = fg2-(fuds2*0.3)/3.;
// 	vMatrix[0](1,1) = fc2-(fuds2*0.3)/3.-TMath::Abs((fg2-(fuds2*0.1))/2.);
// 	vMatrix[0](1,2) = fb2-(fuds2*0.3)/3.-TMath::Abs((fg2-(fuds2*0.1))/2.);
// 	vMatrix[0](1,3) = 0;
	vMatrix[0](2,0) = fuds3*1.3;
	vMatrix[0](2,1) = fc3-(fuds3*0.3)/3.;
	vMatrix[0](2,2) = fb3-(fuds3*0.3)/3.;
	vMatrix[0](2,3) = fg3-(fuds3*0.3)/3.;
	vMatrix[0](3,0) = fuds4*1.3;
// 	vMatrix[0](3,1) = fc4-(fuds4*0.3)/3.-TMath::Abs((fb4-(fuds4*0.3)/3.)/2.);
// 	vMatrix[0](3,2) = 0.;
// 	vMatrix[0](3,3) = fg4-(fuds4*0.3)/3.-TMath::Abs((fb4-(fuds4*0.3)/3.)/2.);
	vMatrix[0](3,1) = fc4-(fuds4*0.3)/3.;
	vMatrix[0](3,2) = fb4-(fuds4*0.3)/3.;
	vMatrix[0](3,3) = fg4-(fuds4*0.3)/3.;
// 	vMatrix[0](3,1) = 0;
// 	vMatrix[0](3,2) = 0;
// 	vMatrix[0](3,3) = fg4-(fuds4*0.3)/3.-TMath::Abs(fc4-(fuds4*0.3)/3.)-TMath::Abs(fb4-(fuds4*0.3)/3.);
	
	//cout<<"matrice vMatrix[0]..."<<endl;
	for (int i = vMatrix[0].GetRowLwb(); i <= vMatrix[0].GetRowUpb(); i++) {
    	for (int j = vMatrix[0].GetColLwb(); j <= vMatrix[0].GetColUpb(); j++) {
      		//cout<<"vMatrix[0]("<<i<<","<<j<<") = "<<vMatrix[0](i,j)<<endl;
      	}
	}
	
	vMatrix[1](0,0) = fuds1*0.7;
	vMatrix[1](0,1) = fc1+(fuds1*0.3)/3.;
	vMatrix[1](0,2) = fb1+(fuds1*0.3)/3.;
	vMatrix[1](0,3) = fg1+(fuds1*0.3)/3.;
	vMatrix[1](1,0) = fuds2*0.7;
	vMatrix[1](1,1) = fc2+(fuds2*0.3)/3.;
	vMatrix[1](1,2) = fb2+(fuds2*0.3)/3.;
	vMatrix[1](1,3) = fg2+(fuds2*0.3)/3.;
	vMatrix[1](2,0) = fuds3*0.7;
	vMatrix[1](2,1) = fc3+(fuds3*0.3)/3.;
	vMatrix[1](2,2) = fb3+(fuds3*0.3)/3.;
	vMatrix[1](2,3) = fg3+(fuds3*0.3)/3.;
	vMatrix[1](3,0) = fuds4*0.7;
	vMatrix[1](3,1) = fc4+(fuds4*0.3)/3.;
	vMatrix[1](3,2) = fb4+(fuds4*0.3)/3.;
	vMatrix[1](3,3) = fg4+(fuds4*0.3)/3.;
	
	//cout<<"matrice vMatrix[1]..."<<endl;
	for (int i = vMatrix[1].GetRowLwb(); i <= vMatrix[1].GetRowUpb(); i++) {
    	for (int j = vMatrix[1].GetColLwb(); j <= vMatrix[1].GetColUpb(); j++) {
      		//cout<<"vMatrix[1]("<<i<<","<<j<<") = "<<vMatrix[1](i,j)<<endl;
      	}
	}
	
	vMatrix[2](0,0) = fuds1-(fc1*0.1);
	vMatrix[2](0,1) = fc1*1.3;
	vMatrix[2](0,2) = fb1-(fc1*0.1);
	vMatrix[2](0,3) = fg1-(fc1*0.1);
// 	vMatrix[2](0,0) = fuds1-(fc1*0.1)-TMath::Abs((fb1-(fc1*0.1))/2.);
// 	vMatrix[2](0,2) = 0;
// 	vMatrix[2](0,3) = fg1-(fc1*0.1)-TMath::Abs((fb1-(fc1*0.1))/2.);
// 	vMatrix[2](1,0) = 0.;
// 	vMatrix[2](1,1) = 1.;
// 	vMatrix[2](1,2) = 0.;
// 	vMatrix[2](1,3) = 0.;
	vMatrix[2](1,0) = fuds2-(fc2*0.1);
	vMatrix[2](1,1) = fc2*1.3;
	vMatrix[2](1,2) = fb2-(fc2*0.1);
	vMatrix[2](1,3) = fg2-(fc2*0.1);
// 	vMatrix[2](1,0) = fuds2-(fc2*0.1)-TMath::Abs((fg2-(fc2*0.1))/2.);
// 	vMatrix[2](1,2) = fb2-(fc2*0.1)-TMath::Abs((fg2-(fc2*0.1))/2.);
// 	vMatrix[2](1,3) = 0;
	vMatrix[2](2,0) = fuds3-(fc3*0.1);
	vMatrix[2](2,1) = fc3*1.3;
	vMatrix[2](2,2) = fb3-(fc3*0.1);
	vMatrix[2](2,3) = fg3-(fc3*0.1);
// 	vMatrix[2](2,0) = fuds3-(fc3*0.1)-TMath::Abs((fg3-(fc3*0.1))/2.);
// 	vMatrix[2](2,2) = fb3-(fc3*0.1)-TMath::Abs((fg3-(fc3*0.1))/2.);
// 	vMatrix[2](2,3) = 0;
	vMatrix[2](3,0) = fuds4-(fc4*0.1);
	vMatrix[2](3,1) = fc4*1.3;
	vMatrix[2](3,2) = fb4-(fc4*0.1);
	vMatrix[2](3,3) = fg4-(fc4*0.1);
// 	vMatrix[2](3,0) = fuds4-(fc4*0.1)-TMath::Abs((fb4-(fc4*0.1))/2.);
// 	vMatrix[2](3,2) = 0;
// 	vMatrix[2](3,3) = fg4-(fc4*0.1)-TMath::Abs((fb4-(fc4*0.1))/2.);
	
	//cout<<"matrice vMatrix[2]..."<<endl;
	for (int i = vMatrix[2].GetRowLwb(); i <= vMatrix[2].GetRowUpb(); i++) {
    	for (int j = vMatrix[2].GetColLwb(); j <= vMatrix[2].GetColUpb(); j++) {
      		//cout<<"vMatrix[2]("<<i<<","<<j<<") = "<<vMatrix[2](i,j)<<endl;
      	}
	}
	
	vMatrix[3](0,0) = fuds1+(fc1*0.1);
	vMatrix[3](0,1) = fc1*0.7;
	vMatrix[3](0,2) = fb1+(fc1*0.1);
	vMatrix[3](0,3) = fg1+(fc1*0.1);
	vMatrix[3](1,0) = fuds2+(fc2*0.1);
	vMatrix[3](1,1) = fc2*0.7;
	vMatrix[3](1,2) = fb2+(fc2*0.1);
	vMatrix[3](1,3) = fg2+(fc2*0.1);
	vMatrix[3](2,0) = fuds3+(fc3*0.1);
	vMatrix[3](2,1) = fc3*0.7;
	vMatrix[3](2,2) = fb3+(fc3*0.1);
	vMatrix[3](2,3) = fg3+(fc3*0.1);
	vMatrix[3](3,0) = fuds4+(fc4*0.1);
	vMatrix[3](3,1) = fc4*0.7;
	vMatrix[3](3,2) = fb4+(fc4*0.1);
	vMatrix[3](3,3) = fg4+(fc4*0.1);
	
	//cout<<"matrice vMatrix[3]..."<<endl;
	for (int i = vMatrix[3].GetRowLwb(); i <= vMatrix[3].GetRowUpb(); i++) {
    	for (int j = vMatrix[3].GetColLwb(); j <= vMatrix[3].GetColUpb(); j++) {
      		//cout<<"vMatrix[3]("<<i<<","<<j<<") = "<<vMatrix[3](i,j)<<endl;
      	}
	}
	
	vMatrix[4](0,0) = fuds1-(fb1*0.1);
	vMatrix[4](0,1) = fc1-(fb1*0.1);
	vMatrix[4](0,2) = fb1*1.3;
	vMatrix[4](0,3) = fg1-(fb1*0.1);
	vMatrix[4](1,0) = fuds2-(fb2*0.1);
	vMatrix[4](1,1) = fc2-(fb2*0.1);
	vMatrix[4](1,2) = fb2*1.3;
	vMatrix[4](1,3) = fg2-(fb2*0.1);
	vMatrix[4](2,0) = fuds3-(fb3*0.1);
	vMatrix[4](2,1) = fc3-(fb3*0.1);
	vMatrix[4](2,2) = fb3*1.3;
	vMatrix[4](2,3) = fg3-(fb3*0.1);
// 	vMatrix[4](2,0) = 0.;
// 	vMatrix[4](2,1) = 0.;
// 	vMatrix[4](2,2) = 1.;
// 	vMatrix[4](2,3) = 0.;
	vMatrix[4](3,0) = fuds4-(fb4*0.1);
	vMatrix[4](3,1) = fc4-(fb4*0.1);
	vMatrix[4](3,2) = fb4*1.3;
	vMatrix[4](3,3) = fg4-(fb4*0.1);
	
	//cout<<"matrice vMatrix[4]..."<<endl;
	for (int i = vMatrix[4].GetRowLwb(); i <= vMatrix[4].GetRowUpb(); i++) {
    	for (int j = vMatrix[4].GetColLwb(); j <= vMatrix[4].GetColUpb(); j++) {
      		//cout<<"vMatrix[4]("<<i<<","<<j<<") = "<<vMatrix[4](i,j)<<endl;
      	}
	}
	
	vMatrix[5](0,0) = fuds1+(fb1*0.1);
	vMatrix[5](0,1) = fc1+(fb1*0.1);
	vMatrix[5](0,2) = fb1*0.7;
	vMatrix[5](0,3) = fg1+(fb1*0.1);
	vMatrix[5](1,0) = fuds2+(fb2*0.1);
	vMatrix[5](1,1) = fc2+(fb2*0.1);
	vMatrix[5](1,2) = fb2*0.7;
	vMatrix[5](1,3) = fg2+(fb2*0.1);
	vMatrix[5](2,0) = fuds3+(fb3*0.1);
	vMatrix[5](2,1) = fc3+(fb3*0.1);
	vMatrix[5](2,2) = fb3*0.7;
	vMatrix[5](2,3) = fg3+(fb3*0.1);
	vMatrix[5](3,0) = fuds4+(fb4*0.1);
	vMatrix[5](3,1) = fc4+(fb4*0.1);
	vMatrix[5](3,2) = fb4*0.7;
	vMatrix[5](3,3) = fg4+(fb4*0.1);
	
	//cout<<"matrice vMatrix[5]..."<<endl;
	for (int i = vMatrix[5].GetRowLwb(); i <= vMatrix[5].GetRowUpb(); i++) {
    	for (int j = vMatrix[5].GetColLwb(); j <= vMatrix[5].GetColUpb(); j++) {
      		//cout<<"vMatrix[5]("<<i<<","<<j<<") = "<<vMatrix[5](i,j)<<endl;
      	}
	}
	
	vMatrix[6](0,0) = fuds1-(fg1*0.1);
	vMatrix[6](0,1) = fc1-(fg1*0.1);
	vMatrix[6](0,2) = fb1-(fg1*0.1);
// 	vMatrix[6](0,0) = fuds1-(fg1*0.1)-TMath::Abs((fb1-(fg1*0.1))/2.);
// 	vMatrix[6](0,1) = fc1-(fg1*0.1)-TMath::Abs((fb1-(fg1*0.1))/2.);
// 	vMatrix[6](0,2) = 0;
	vMatrix[6](0,3) = fg1*1.3;
	vMatrix[6](1,0) = fuds2-(fg2*0.1);
	vMatrix[6](1,1) = fc2-(fg2*0.1);
	vMatrix[6](1,2) = fb2-(fg2*0.1);
	vMatrix[6](1,3) = fg2*1.3;
	vMatrix[6](2,0) = fuds3-(fg3*0.1);
	vMatrix[6](2,1) = fc3-(fg3*0.1);
	vMatrix[6](2,2) = fb3-(fg3*0.1);
	vMatrix[6](2,3) = fg3*1.3;
	vMatrix[6](3,0) = fuds4-(fg4*0.1);
	vMatrix[6](3,1) = fc4-(fg4*0.1);
	vMatrix[6](3,2) = fb4-(fg4*0.1);
// 	vMatrix[6](3,0) = fuds4-(fg4*0.1)-TMath::Abs(fc4-(fg4*0.1))-TMath::Abs(fb4-(fg4*0.1));
// 	vMatrix[6](3,1) = 0;
// 	vMatrix[6](3,2) = 0;
	vMatrix[6](3,3) = fg4*1.3;
// 	vMatrix[6](3,0) = 0.;
// 	vMatrix[6](3,1) = 0.;
// 	vMatrix[6](3,2) = 0.;
// 	vMatrix[6](3,3) = 1.;
	
	//cout<<"matrice vMatrix[6]..."<<endl;
	for (int i = vMatrix[6].GetRowLwb(); i <= vMatrix[6].GetRowUpb(); i++) {
    	for (int j = vMatrix[6].GetColLwb(); j <= vMatrix[6].GetColUpb(); j++) {
      		//cout<<"vMatrix[6]("<<i<<","<<j<<") = "<<vMatrix[6](i,j)<<endl;
      	}
	}
	
	vMatrix[7](0,0) = fuds1+(fg1*0.1);
	vMatrix[7](0,1) = fc1+(fg1*0.1);
	vMatrix[7](0,2) = fb1+(fg1*0.1);
	vMatrix[7](0,3) = fg1*0.7;
	vMatrix[7](1,0) = fuds2+(fg2*0.1);
	vMatrix[7](1,1) = fc2+(fg2*0.1);
	vMatrix[7](1,2) = fb2+(fg2*0.1);
	vMatrix[7](1,3) = fg2*0.7;
	vMatrix[7](2,0) = fuds3+(fg3*0.1);
	vMatrix[7](2,1) = fc3+(fg3*0.1);
	vMatrix[7](2,2) = fb3+(fg3*0.1);
	vMatrix[7](2,3) = fg3*0.7;
	vMatrix[7](3,0) = fuds4+(fg4*0.1);
	vMatrix[7](3,1) = fc4+(fg4*0.1);
	vMatrix[7](3,2) = fb4+(fg4*0.1);
	vMatrix[7](3,3) = fg4*0.7;
	
	//cout<<"matrice vMatrix[7]..."<<endl;
	for (int i = vMatrix[7].GetRowLwb(); i <= vMatrix[7].GetRowUpb(); i++) {
    	for (int j = vMatrix[7].GetColLwb(); j <= vMatrix[7].GetColUpb(); j++) {
      		//cout<<"vMatrix[7]("<<i<<","<<j<<") = "<<vMatrix[7](i,j)<<endl;
      	}
	}
	

	for(int k=0; k<2*getNcolumns(); k++) {
		vMatrix_inverse[k] = ((TMatrixF*)vMatrix[k].Clone())->Invert();	
		for (int i = vMatrix_inverse[k].GetRowLwb(); i <= vMatrix_inverse[k].GetRowUpb(); i++) {
				for (int j = vMatrix_inverse[k].GetColLwb(); j <= vMatrix_inverse[k].GetColUpb(); j++) {
						//cout<<"vMatrix_inverse["<<k<<"]("<<i<<","<<j<<") = "<<vMatrix_inverse[k](i,j)<<endl;
					}
		}		
		//cout<<"verification matrice vMatrix["<<k<<"] inchangee..."<<endl;
		for (int i = vMatrix[k].GetRowLwb(); i <= vMatrix[k].GetRowUpb(); i++) {
				for (int j = vMatrix[k].GetColLwb(); j <= vMatrix[k].GetColUpb(); j++) {
						//cout<<"vMatrix["<<k<<"]("<<i<<","<<j<<") = "<<vMatrix[k](i,j)<<endl;
					}
		}
	}
	
	return 	vMatrix_inverse;
*/
}

vector<TGraphErrors*> do2DTaggingMethod(TFile* afMC_matrix, TFile* f, string aPtName, int withSystErr = 0) {
//*******************************************************************************************************
// description: function which applies the 2D tagging method
//              i.e. reconstructs the mean and RMS from responses per zone and flavour matrix
//							vector.at(0): result TGraphErrors for the mean flavour response
//							vector.at(1): result TGraphErrors for the 'RMS' flavour response
// arguments: - afMC_matrix: MC rootfile containning the flavour fraction matrix
//            - f: rootfile containing the response distributions per zone usefull for the 2DTagging method
//						     This file can be a MC file to evaluate the performance of the 2DTagging method
//						     or a data file 
//            - aPtName: name of the 2DTagging pt bin under study
//            - withSystErr: 1 to compute systematic errors du to flavour fraction flucutation
//*******************************************************************************************************
	vector<TGraphErrors*> myVector;
	vector<TH1F*> vHisto_Zone;
	TMatrixF m0 = buildSquareMatrix();
	TMatrixF m0_inverse = buildSquareMatrix();

	//retrieve the flavour matrix and inverse it
	m0 = getMatrix(afMC_matrix, aPtName);
	m0_inverse = ((TMatrixF*)m0.Clone())->Invert();
	for (int i = m0_inverse.GetRowLwb(); i <= m0_inverse.GetRowUpb(); i++) {
    	for (int j = m0_inverse.GetColLwb(); j <= m0_inverse.GetColUpb(); j++) {
      		cout<<"m0_inverse("<<i<<","<<j<<") = "<<m0_inverse(i,j)<<endl;
      	}
	}

	//built the 8 systematic error flavour matrices
	vector<TMatrixF> vm_inverse = buildErrorInverseMatrix(m0);

	//retrieve the response distributions per zone, for a given 2DTagging pt bin
	vHisto_Zone = getResponse_Zone_Pt(f,aPtName);

	//extract mean and RMS from 2DTagging method and compute statistical errors
	TMatrixF mRzonePlusRMS(getNrows(),1);
	TMatrixF mRzoneMinusRMS(getNrows(),1);
	TMatrixF mRzonePlusMeanError(getNrows(),1);
	TMatrixF mRzoneMinusMeanError(getNrows(),1);
	TMatrixF mRzonePlusRMSError(getNrows(),1);
	TMatrixF mRzoneMinusRMSError(getNrows(),1);
	TMatrixF mRzone(getNrows(),1);

	
	for(int i=0; i<getNrows(); i++) {
		//vHisto_Zone[i]->Draw();
		cout<<"in the loop......."<<endl;
		mRzonePlusRMS(i,0) = vHisto_Zone[i]->GetMean() + vHisto_Zone[i]->GetRMS();
		cout<<"mRzonePlusRMS......."<<endl;
		mRzoneMinusRMS(i,0) = vHisto_Zone[i]->GetMean() - vHisto_Zone[i]->GetRMS();
		mRzonePlusMeanError(i,0) = vHisto_Zone[i]->GetMean() + vHisto_Zone[i]->GetMeanError();
		mRzoneMinusMeanError(i,0) = vHisto_Zone[i]->GetMean() - vHisto_Zone[i]->GetMeanError();
		mRzonePlusRMSError(i,0) = vHisto_Zone[i]->GetMean() + vHisto_Zone[i]->GetRMS() + vHisto_Zone[i]->GetRMSError();
		mRzoneMinusRMSError(i,0) = vHisto_Zone[i]->GetMean() - vHisto_Zone[i]->GetRMS() - vHisto_Zone[i]->GetRMSError();
		mRzone(i,0) = vHisto_Zone[i]->GetMean();
	}


	TMatrix mRflavour_2Dtagging(4,1);
	mRflavour_2Dtagging = m0_inverse*mRzone;
	TMatrix mRflavour_2DtaggingPlusRMS(4,1);
	mRflavour_2DtaggingPlusRMS = m0_inverse*mRzonePlusRMS;
	TMatrix mRflavour_2DtaggingMinusRMS(4,1);
	mRflavour_2DtaggingMinusRMS = m0_inverse*mRzoneMinusRMS;
	TMatrix mRflavour_2DtaggingPlusMeanError(4,1);
	mRflavour_2DtaggingPlusMeanError = m0_inverse*mRzonePlusMeanError;
	TMatrix mRflavour_2DtaggingMinusMeanError(4,1);
	mRflavour_2DtaggingMinusMeanError = m0_inverse*mRzoneMinusMeanError;
	TMatrix mRflavour_2DtaggingPlusRMSError(4,1);
	mRflavour_2DtaggingPlusRMSError = m0_inverse*mRzonePlusRMSError;
	TMatrix mRflavour_2DtaggingMinusRMSError(4,1);
	mRflavour_2DtaggingMinusRMSError = m0_inverse*mRzoneMinusRMSError;


	float aXflavour[4];
	aXflavour[0] = 1;
	aXflavour[1] = 2;
	aXflavour[2] = 3;
	aXflavour[3] = 4;
	float aXflavour_error[4];
	aXflavour_error[0] = 0;
	aXflavour_error[1] = 0;
	aXflavour_error[2] = 0;
	aXflavour_error[3] = 0;

	float aMean_2Dtagging[4];
	float aRMS_2Dtagging[4];
	float aMeanErrorStatL_2Dtagging[4];
	float aMeanErrorStatH_2Dtagging[4];
	float aRMSErrorStatL_2Dtagging[4];
	float aRMSErrorStatH_2Dtagging[4];

    for(int i=0; i<4; i++) {
    	aMean_2Dtagging[i]= mRflavour_2Dtagging(i,0);
			aRMS_2Dtagging[i] = (mRflavour_2DtaggingPlusRMS(i,0) - mRflavour_2DtaggingMinusRMS(i,0))/2.;
    	aMeanErrorStatH_2Dtagging[i]=(mRflavour_2DtaggingPlusMeanError(i,0) - mRflavour_2Dtagging(i,0));
    	aMeanErrorStatL_2Dtagging[i]=(mRflavour_2Dtagging(i,0) - mRflavour_2DtaggingMinusMeanError(i,0));
    	aRMSErrorStatH_2Dtagging[i]=(mRflavour_2DtaggingPlusRMSError(i,0)) - (mRflavour_2DtaggingPlusRMS(i,0));
    	aRMSErrorStatL_2Dtagging[i]=(mRflavour_2DtaggingMinusRMS(i,0) - mRflavour_2DtaggingMinusRMSError(i,0));   	
			cout<<"aMean_2Dtagging["<<i<<"] : "<<aMean_2Dtagging[i]<<endl;
			cout<<"aRMS_2Dtagging["<<i<<"] : "<<aRMS_2Dtagging[i]<<endl;
			cout<<"aMeanErrorStatH_2Dtagging["<<i<<"] : "<<aMeanErrorStatH_2Dtagging[i]<<endl;
			cout<<"aMeanErrorStatL_2Dtagging["<<i<<"] : "<<aMeanErrorStatL_2Dtagging[i]<<endl;
			cout<<"aRMSErrorStatH_2Dtagging["<<i<<"] : "<<aRMSErrorStatH_2Dtagging[i]<<endl;
			cout<<"aRMSErrorStatL_2Dtagging["<<i<<"] : "<<aRMSErrorStatL_2Dtagging[i]<<endl<<endl;

    }

	float aMeanErrorStat_2Dtagging[4];
	float aRMSErrorStat_2Dtagging[4];
	for(int j=0; j<4; j++) {
		aMeanErrorStat_2Dtagging[j] = (aMeanErrorStatL_2Dtagging[j]+aMeanErrorStatH_2Dtagging[j])/2.;
		cout<<"aMeanErrorStat_2Dtagging["<<j<<"] : "<<aMeanErrorStat_2Dtagging[j]<<endl;

		aRMSErrorStat_2Dtagging[j] = (aRMSErrorStatL_2Dtagging[j]+aRMSErrorStatH_2Dtagging[j])/2.;
		cout<<"aRMSErrorStat_2Dtagging["<<j<<"] : "<<aRMSErrorStat_2Dtagging[j]<<endl<<endl;
	}

	float aMeanErrorSyst_2Dtagging[4];
	float aMeanErrorTot_2Dtagging[4];
	float aRMSErrorSyst_2Dtagging[4];
	float aRMSErrorTot_2Dtagging[4];

	if(withSystErr == 1) {
		//compute systematic errors from fraction fluctuations
		vector<TMatrixF> vColumnMatrix_Mean;
		vector<TMatrixF> vColumnMatrix_RMS_plus;
		vector<TMatrixF> vColumnMatrix_RMS_minus;
		vector<float*> vMean;
		vector<TH1F*> vHistoMean;
		vector<float*> vRMS;
		vector<TH1F*> vHistoRMS;
		for(int i=0; i<2*getNcolumns(); i++) {
			vColumnMatrix_Mean.push_back(buildColumnMatrix());
			vColumnMatrix_RMS_plus.push_back(buildColumnMatrix());
			vColumnMatrix_RMS_minus.push_back(buildColumnMatrix());
	
			vColumnMatrix_Mean[i] = vm_inverse[i]*mRzone;	
			vColumnMatrix_RMS_plus[i] = vm_inverse[i]*mRzonePlusRMS;	
			vColumnMatrix_RMS_minus[i] = vm_inverse[i]*mRzoneMinusRMS;	
		}
	
		vector<TGraph*> vGraphMean;
		vector<TGraph*> vGraphRMS;
	
		for(int j=0; j<4; j++) {
			vMean.push_back(new float[2*getNcolumns()]);
			vRMS.push_back(new float[2*getNcolumns()]);
			stringstream aHistoMeanName;
			aHistoMeanName <<"h"<<j<<"Mean";
			string myMeanHistoName = aHistoMeanName.str() + aPtName;
			vHistoMean.push_back(new TH1F(myMeanHistoName.c_str(),myMeanHistoName.c_str(),200,0.7,1.2));
			stringstream aHistoRMSName;
			aHistoRMSName <<"h"<<j<<"RMS";
			string myRMSHistoName = aHistoRMSName.str() + aPtName;
			vHistoRMS.push_back(new TH1F(myRMSHistoName.c_str(),myRMSHistoName.c_str(),200,0.08,0.20));
			for(int i=0; i<2*getNcolumns(); i++) {
				vMean[j][i] = vColumnMatrix_Mean[i](j,0);
				vHistoMean[j]->Fill(vColumnMatrix_Mean[i](j,0));
				vRMS[j][i] = (vColumnMatrix_RMS_plus[i](j,0)-vColumnMatrix_RMS_minus[i](j,0))/2.;
				vHistoRMS[j]->Fill((vColumnMatrix_RMS_plus[i](j,0)-vColumnMatrix_RMS_minus[i](j,0))/2.);
				cout<<"vMean["<<j<<"]["<<i<<"] : "<<vMean[j][i]<<endl;
				cout<<"vRMS["<<j<<"]["<<i<<"] : "<<vRMS[j][i]<<endl<<endl;
			}
		}

		//compute total error (syst+stat)
		for(int j=0; j<4; j++) {
			vGraphMean.push_back(new TGraph(4,aXflavour,vMean[j]));
			vGraphMean[j]->Draw("APE");
			aMeanErrorSyst_2Dtagging[j] = vHistoMean[j]->GetRMS();
			aMeanErrorTot_2Dtagging[j] = sqrt(pow(aMeanErrorSyst_2Dtagging[j],2) + pow(aMeanErrorStat_2Dtagging[j],2));
			cout<<"aMeanErrorSyst_2Dtagging["<<j<<"] : "<<aMeanErrorSyst_2Dtagging[j]<<endl;
			cout<<"aMeanErrorTot_2Dtagging["<<j<<"] : "<<aMeanErrorTot_2Dtagging[j]<<endl;
	
			vGraphRMS.push_back(new TGraph(4,aXflavour,vRMS[j]));
			vGraphRMS[j]->Draw("APE");
			aRMSErrorSyst_2Dtagging[j] = vHistoRMS[j]->GetRMS();
			aRMSErrorTot_2Dtagging[j] = sqrt(pow(aRMSErrorSyst_2Dtagging[j],2) + pow(aRMSErrorStat_2Dtagging[j],2));
			cout<<"aRMSErrorSyst_2Dtagging["<<j<<"] : "<<aRMSErrorSyst_2Dtagging[j]<<endl;
			cout<<"vHistoRMS["<<j<<"]->GetRMS() : "<<vHistoRMS[j]->GetRMS()<<endl;
			cout<<"aRMSErrorTot_2Dtagging["<<j<<"] : "<<aRMSErrorTot_2Dtagging[j]<<endl<<endl;
		}

		vHistoMean.clear();
		vHistoRMS.clear();

	}
	
	else if(withSystErr == 0) {
		for(int j=0; j<4; j++) {
			aMeanErrorTot_2Dtagging[j] = aMeanErrorStat_2Dtagging[j];
			cout<<"aMeanErrorTot_2Dtagging["<<j<<"] : "<<aMeanErrorTot_2Dtagging[j]<<endl;
	
			aRMSErrorTot_2Dtagging[j] = aRMSErrorStat_2Dtagging[j];
			cout<<"aRMSErrorTot_2Dtagging["<<j<<"] : "<<aRMSErrorTot_2Dtagging[j]<<endl<<endl;
		}
	}

	//construct the TGraphErrors for mean and RMS results
	TGraphErrors* gRmpf_flavour_2Dtagging_Mean = new TGraphErrors(4,aXflavour, aMean_2Dtagging, aXflavour_error, aMeanErrorTot_2Dtagging);
  gRmpf_flavour_2Dtagging_Mean->SetName("gRmpf_flavour_2Dtagging_Mean");

	TGraphErrors* gRmpf_flavour_2Dtagging_RMS = new TGraphErrors(4,aXflavour, aRMS_2Dtagging, aXflavour_error, aRMSErrorTot_2Dtagging);
  gRmpf_flavour_2Dtagging_RMS->SetName("gRmpf_flavour_2Dtagging_RMS");


	myVector.push_back(gRmpf_flavour_2Dtagging_Mean);
	myVector.push_back(gRmpf_flavour_2Dtagging_RMS);
	return myVector;
}

vector<TGraphErrors*> getMCTruth(TFile* afMC, string aPtName) {
	//function which builds the MC truth TGraphErrors i.e. exctracts true mean and RMS 
	//from the response distributions per flavour obtained from MC truth
	vector<TGraphErrors*> myVector;
	vector<TH1F*> vHisto_Flavour_MC;
	vHisto_Flavour_MC = getResponse_Flavour_Pt(afMC, aPtName);

	float aRflavourRMS_MC[getNcolumns()];
	float aRflavourMeanError_MC[getNcolumns()];
	float aRflavourRMSError_MC[getNcolumns()];
	float aRflavourMean_MC[getNcolumns()];
	
	for(int i=0; i<getNcolumns(); i++) {
		aRflavourRMS_MC[i] = vHisto_Flavour_MC[i]->GetRMS();
		aRflavourMeanError_MC[i] = vHisto_Flavour_MC[i]->GetMeanError();
		aRflavourRMSError_MC[i] = vHisto_Flavour_MC[i]->GetRMSError();
		aRflavourMean_MC[i] = vHisto_Flavour_MC[i]->GetMean();
		cout<<"aRflavourRMS_MC["<<i<<"] : "<<aRflavourRMS_MC[i]<<endl;
		cout<<"aRflavourMeanError_MC["<<i<<"] : "<<aRflavourMeanError_MC[i]<<endl;
		cout<<"aRflavourRMSError_MC["<<i<<"] : "<<aRflavourRMSError_MC[i]<<endl;
		cout<<"aRflavourMean_MC["<<i<<"] : "<<aRflavourMean_MC[i]<<endl<<endl;
	}

	float aXflavour[4];
	aXflavour[0] = 1;
	aXflavour[1] = 2;
	aXflavour[2] = 3;
	aXflavour[3] = 4;
	float aXflavour_error[4];
	aXflavour_error[0] = 0;
	aXflavour_error[1] = 0;
	aXflavour_error[2] = 0;
	aXflavour_error[3] = 0;

	TGraphErrors* gRmpf_flavour_MCtruth_Mean = new TGraphErrors(4,aXflavour, aRflavourMean_MC, aXflavour_error, aRflavourMeanError_MC);
	gRmpf_flavour_MCtruth_Mean->SetName("gRmpf_flavour_MCtruth_Mean");
	gRmpf_flavour_MCtruth_Mean->SetMarkerStyle(24);
	gRmpf_flavour_MCtruth_Mean->SetMarkerColor(2);
	gRmpf_flavour_MCtruth_Mean->SetLineColor(2);

	TGraphErrors* gRmpf_flavour_MCtruth_RMS = new TGraphErrors(4,aXflavour, aRflavourRMS_MC, aXflavour_error, aRflavourRMSError_MC);
	gRmpf_flavour_MCtruth_RMS->SetName("gRmpf_flavour_MCtruth_RMS");
	gRmpf_flavour_MCtruth_RMS->SetMarkerStyle(24);
	gRmpf_flavour_MCtruth_RMS->SetMarkerColor(2);
	gRmpf_flavour_MCtruth_RMS->SetLineColor(2);

	myVector.push_back(gRmpf_flavour_MCtruth_Mean);
	myVector.push_back(gRmpf_flavour_MCtruth_RMS);

	return myVector;
}

void drawRelativGraph(TGraphErrors* gMC, TGraphErrors* gData, string aYTitle, string aPath) {
//*******************************************************************************************************
// description: function which draws the mean/RMS flavour responses comparison in relative difference
// arguments: 1) 2DTagging method performance study
//                - gMC: TGraphErrors for MC truth
//                - gData: TGraphErrors obtained from 2DTagging method for MC 
//            2) 2DTagging method study on data
//                - gMC: TGraphErrors obtained from 2DTagging method for MC 
//                - gData: TGraphErrors obtained from 2DTagging method for data
//            - aYTitle: Y axis title
//            - aPath: where you want to save the plot
//*******************************************************************************************************
	Double_t arelative[4];
	Double_t arelative_error[4];
	Double_t aMC_error[4];
	Double_t aData_error[4];
	Double_t aMCx[4];
	Double_t aDatax[4];
	Double_t aMCy[4];
	Double_t aDatay[4];
	
	for(Int_t i=0; i<4; i++) { 
// 		aMC_error[i] = gMC->GetErrorX(gMC->GetXaxis()->FindBin(i+1));
// 		aData_error[i] = gData->GetErrorX(gData->GetXaxis()->FindBin(i+1));
		aMC_error[i] = gMC->GetErrorY(i);
		aData_error[i] = gData->GetErrorY(i);
		gMC->GetPoint(i, aMCx[i], aMCy[i]);
		gData->GetPoint(i, aDatax[i], aDatay[i]);		  	

		arelative[i] = 100.*(aDatay[i]-aMCy[i])/aMCy[i];
		arelative_error[i] = 100.*sqrt(pow(aData_error[i],2)/pow(aMCy[i],2)+pow(aMC_error[i],2) * pow(aDatay[i]/(pow(aMCy[i],2)),2));
	
		cout<<"aDatay["<<i<<"] : "<<aDatay[i]<<endl;
		cout<<"aData_error["<<i<<"] : "<<aData_error[i]<<endl;
		
		cout<<"aMCy["<<i<<"] : "<<aMCy[i]<<endl;
		cout<<"aMC_error["<<i<<"] : "<<aMC_error[i]<<endl;
		
		cout<<"arelative["<<i<<"] : "<<arelative[i]<<endl;
		cout<<"arelative_error["<<i<<"] : "<<arelative_error[i]<<endl<<endl;		
	}

	Double_t aXflavour[4];
	aXflavour[0] = 1;
	aXflavour[1] = 2;
	aXflavour[2] = 3;
	aXflavour[3] = 4;
	Double_t aXflavour_error[4];
	aXflavour_error[0] = 0;
	aXflavour_error[1] = 0;
	aXflavour_error[2] = 0;
	aXflavour_error[3] = 0;

	TGraphErrors* gR_flavour_relative = new TGraphErrors(4,aXflavour, arelative, aXflavour_error, arelative_error);
	gR_flavour_relative->SetMarkerStyle(21);
	gR_flavour_relative->SetMarkerColor(2);
	gR_flavour_relative->SetLineColor(2);
	
	TGraph_style (gR_flavour_relative);


	TCanvas *crelative = new TCanvas("crelative","crelative");
	crelative->cd();
	gR_flavour_relative->Draw("APE");
	//gR_flavour_relative->GetYaxis()->SetRangeUser(0.70,1.15);
	gR_flavour_relative->GetYaxis()->SetTitle(aYTitle.c_str());
	gR_flavour_relative->GetXaxis()->SetNdivisions(4,kTRUE);
	for(int i=0; i<getFlavourNumber()-2; i++) {
		gR_flavour_relative->GetXaxis()->SetBinLabel((gR_flavour_relative->GetXaxis()->FindBin(i+1)),getFlavourBinName(i).c_str());
	}
	gR_flavour_relative->GetXaxis()->LabelsOption("h");
	gPad->RedrawAxis();
	crelative->SaveAs(aPath.c_str());
}

void drawGraph(TCanvas* aCanvas, TGraphErrors* gtruth, TGraphErrors* gMC, TGraphErrors* gdata, string YTitle, string path, int isRMS = 0, int ngraphs = 3) {
//*******************************************************************************************************
// description: function which draws the mean/RMS flavour responses comparison 
// arguments: - aCanvas: canvas where you want to draw your plots
//            - gtruth: TGraphErrors for MC truth (mean or RMS)
//            - gMC: TGraphErrors for MC 2DTagging (mean or RMS)
//            - gdata: TGraphErrors for data 2DTagging (mean or RMS)
//            - YTitle: title of Y axis (depending if it is mean or RMS)
//            - path: where you want to save the plot
//            - isRMS: you have to precize if it is RMS plot (=1) or not (=0)
//            - ngraphs: you can choose to draw only 2DTagging data and 2DTagging MC using ngraphs=2
//*******************************************************************************************************
 		gStyle->SetOptStat(0);
		TMultiGraph *mg = new TMultiGraph();
		TGraph_flavourLabel(gtruth);
		TGraph_flavourLabel(gMC);
		TGraph_flavourLabel(gdata);
   	aCanvas->cd();
		aCanvas->Clear();

		if(ngraphs == 3) {
			mg->Add(gtruth);
			mg->Add(gMC);
			mg->Add(gdata);
		}
		else if(ngraphs == 2) {
			mg->Add(gtruth);
			mg->Add(gMC);
			//mg->Add(gdata);
		}


   	mg->Draw("APE");
		if(isRMS == 0) {
    		mg->GetYaxis()->SetRangeUser(0.60,1.70);
// 				mg->GetYaxis()->SetRangeUser(0.92,1.26);
		}
		else if(isRMS == 1) {
//    		mg->GetYaxis()->SetRangeUser(0.08,0.18);
		}


   	mg->GetYaxis()->SetTitle(YTitle.c_str());
   	TMultiGraph_style(mg);
   	mg->GetXaxis()->SetNdivisions(4,kTRUE);
		for(int i=0; i<getFlavourNumber()-2; i++) {
			mg->GetXaxis()->SetBinLabel((gtruth->GetXaxis()->FindBin(i+1)),getFlavourBinName(i).c_str());
		}
 		mg->GetXaxis()->LabelsOption("h");
 		mg->GetXaxis()->SetLabelSize(0.055);

		gPad->RedrawAxis();
		TLegend *l = new TLegend(0.12,0.12,0.45,0.25);
		//l->SetBorderSize(0);
		l->SetFillColor(0);
		l->SetFillStyle(0);
		l->SetTextFont(42);
		l->SetTextSizePixels(24);
		if(ngraphs == 3) {
			l->AddEntry(gtruth,"MC truth","p");
			l->AddEntry(gMC,"2D tagging MC","p");
			l->AddEntry(gdata,"2D tagging data","p");
		}
		else if(ngraphs == 2) {
			l->AddEntry(gtruth,"MC truth","p");
			l->AddEntry(gMC,"2D tagging MC","p");
			//l->AddEntry(gdata,"2D tagging data","p");
		}

		l->Draw("same");
		aCanvas->SaveAs(path.c_str());
}


void my2DTaggingMethod_4x4 () 
{

	TString innameMC_G_matrix="../../output_rootfile/outputMatrix2DTagging_MC_TOT.root";
	//TString innameMC_QCD_matrix="../../output_rootfile/outputMatrix2DTagging_MC_QCD.root";
	TString innameData="../../output_rootfile/output2DTagging_data.root";
// 	TString innameMC_TOT_matrix="../../output_rootfile/outputMatrix2DTagging_MC_TOT.root";
	TString innameMC_G="../../output_rootfile/output2DTagging_MC_TOT.root";
	//TString innameMC_QCD="../../output_rootfile/output2DTagging_MC_QCD.root";
// 	TString innameMC_TOT="../../output_rootfile/output2DTagging_MC_TOT.root";

	//open input file
	TFile *fMC_G_matrix=TFile::Open(innameMC_G_matrix);
	//TFile *fMC_QCD_matrix=TFile::Open(innameMC_QCD_matrix);
	TFile *fData=TFile::Open(innameData);
	//TFile *fMC_TOT_matrix=TFile::Open(innameMC_TOT_matrix);
	TFile *fMC_G=TFile::Open(innameMC_G);
	//TFile *fMC_QCD=TFile::Open(innameMC_QCD);
	//TFile *fMC_TOT=TFile::Open(innameMC_TOT);

	int isFor2DTagging = 1;
	ptBinning my2DTaggingPtBinning(isFor2DTagging);
	int nptbins = my2DTaggingPtBinning.getSize();
	int nzones = getZoneNumber();
	int nflavours = getFlavourNumber(); 
	int ntaggers = getTaggerNumber();

	int isForLumi = 0 ;
	int withNoMatched = 1 ;
	int inLogScale = 1;
	int withSystErrors = 0;

	vector<int> myHistoColor = HistoColor();

	TCanvas* c2DTagging = new TCanvas();
	TGraphErrors* gRmpf_flavour_MCtruth_Mean = NULL;
	TGraphErrors* gRmpf_flavour_MCtruth_RMS = NULL;
	TGraphErrors* gRmpf_flavour_2Dtagging_Mean_data = NULL;
	TGraphErrors* gRmpf_flavour_2Dtagging_RMS_data = NULL;
	TGraphErrors* gRmpf_flavour_2Dtagging_Mean_MC = NULL;
	TGraphErrors* gRmpf_flavour_2Dtagging_RMS_MC = NULL;

	vector<TGraphErrors*> vData;
	vector<TGraphErrors*> vMC;
	vector<TGraphErrors*> vMCTruth;
	
	string path;	

	for(int j=0; j<nptbins; j++) {
		string PtName = my2DTaggingPtBinning.getName(j);
		cout<<"*********************************************************"<<endl;
		cout<<PtName<<endl;
		cout<<"*********************************************************"<<endl<<endl;
		cout<<"*********************************************************"<<endl;
		cout<<"do2DTaggingMethod for MC:"<<endl;
		cout<<"*********************************************************"<<endl<<endl;
		vMC = do2DTaggingMethod(fMC_G_matrix, fMC_G, PtName,withSystErrors);
		cout<<"*********************************************************"<<endl;
		cout<<"do2DTaggingMethod for data:"<<endl;
		cout<<"*********************************************************"<<endl<<endl;
		vData = do2DTaggingMethod(fMC_G_matrix, fData, PtName, withSystErrors);
		cout<<"*********************************************************"<<endl;
		cout<<"getMCTruth:"<<endl;
		cout<<"*********************************************************"<<endl<<endl;
		vMCTruth = getMCTruth(fMC_G, PtName);

		gRmpf_flavour_MCtruth_Mean =  vMCTruth[0];
		gRmpf_flavour_MCtruth_RMS = vMCTruth[1];
		gRmpf_flavour_2Dtagging_Mean_data =  vData[0];
		gRmpf_flavour_2Dtagging_RMS_data = vData[1];
		gRmpf_flavour_2Dtagging_Mean_MC = vMC[0]; 
		gRmpf_flavour_2Dtagging_RMS_MC = vMC[1];

		gRmpf_flavour_2Dtagging_Mean_MC->SetMarkerStyle(20);
		gRmpf_flavour_2Dtagging_Mean_MC->SetMarkerColor(4);
		gRmpf_flavour_2Dtagging_Mean_MC->SetLineColor(4);
	
		gRmpf_flavour_2Dtagging_RMS_MC->SetMarkerStyle(20);
		gRmpf_flavour_2Dtagging_RMS_MC->SetMarkerColor(4);
		gRmpf_flavour_2Dtagging_RMS_MC->SetLineColor(4);
	
		gRmpf_flavour_2Dtagging_Mean_data->SetMarkerStyle(20);
		gRmpf_flavour_2Dtagging_Mean_data->SetMarkerColor(1);
		gRmpf_flavour_2Dtagging_Mean_data->SetLineColor(1);
	
		gRmpf_flavour_2Dtagging_RMS_data->SetMarkerStyle(20);
		gRmpf_flavour_2Dtagging_RMS_data->SetMarkerColor(1);
		gRmpf_flavour_2Dtagging_RMS_data->SetLineColor(1);

		path = "plotsResult/Mean_" + PtName + "_MCtot.pdf";
		drawGraph(c2DTagging, gRmpf_flavour_MCtruth_Mean, gRmpf_flavour_2Dtagging_Mean_MC,gRmpf_flavour_2Dtagging_Mean_data, "Mean +/- Stat. Error", path.c_str(),0);

		path = "plotsResult/RMS_" + PtName + "_MCtot.pdf";
		drawGraph(c2DTagging, gRmpf_flavour_MCtruth_RMS, gRmpf_flavour_2Dtagging_RMS_MC,gRmpf_flavour_2Dtagging_RMS_data, "RMS +/- Stat. Error", path.c_str(),1);

		path = "plotsResult/MeanRelative_" + PtName + "_MC_MCTruth_MCtot.pdf";
		drawRelativGraph(gRmpf_flavour_MCtruth_Mean, gRmpf_flavour_2Dtagging_Mean_MC, "Mean relative difference (%)", path.c_str());
		
		path = "plotsResult/RMSRelative_" + PtName + "_MC_MCTruth_MCtot.pdf";
		drawRelativGraph(gRmpf_flavour_MCtruth_RMS, gRmpf_flavour_2Dtagging_RMS_MC, "RMS relative difference (%)", path.c_str());
		
		path = "plotsResult/MeanRelative_" + PtName + "_MC_data_MCtot.pdf";
		drawRelativGraph(gRmpf_flavour_2Dtagging_Mean_MC, gRmpf_flavour_2Dtagging_Mean_data, "Mean relative difference (%)", path.c_str());
		
		path = "plotsResult/RMSRelative_" + PtName + "_MC_data_MCtot.pdf";
		drawRelativGraph(gRmpf_flavour_2Dtagging_RMS_MC, gRmpf_flavour_2Dtagging_RMS_data, "RMS relative difference (%)", path.c_str());

	}	




}


