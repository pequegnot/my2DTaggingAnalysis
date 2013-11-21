int getN_negFraction(TMatrix* m, int row) {
	int n=0;
	for(int j=0; j<4; j++) {
		if(m(row,j) < 0.) {
			n=n+1;
		}
	}
	return n;
}

int getN_nullFraction(TMatrix* m, int row) {
	int n=0;
	for(int j=0; j<4; j++) {
		if(m(row,j) == 0.) {
			n=n+1;
		}
	}
	return n;
}

int getIndex_negFraction(TMatrix* m, int row) {
	int index=0;
	for(int j=0; j<4; j++) {
		if(m(row,j) < 0.) {
			index = j;
		}
	}
	return index;
}

int getIndex_nullFraction(TMatrix* m, int row) {
	int index=0;
	for(int j=0; j<4; j++) {
		if(m(row,j) == 0.) {
			index = j;
		}
	}
	return index;
}

int getIndex_maxFraction(TMatrix* m, int row) {
	int index=0;
	float maxFraction_tmp = m(row,0);
	for(int j=0; j<4; j++) {
		if(m(row,j) > maxFraction_tmp) {
			maxFraction_tmp = m(row,j);
			index = j;
		}
	}
	return j;
}

void setMaxFracTo1(TMatrix* m, int row, int& need_to_check_neg_frac) {
		int index_tmp = 0;
		for(int j=0; j<4; j++) {
		if(m(row,j)>1) {
			index_tmp = j;
			m(row,j) = 1.;
		}
	}
	for(int j=0; j<4; j++) {
		if(j != index_tmp) {
			m(row,j) = 0.;
		}
	}
}

for(int k=0; k<4; k++) {
	for(int i=0; i<4; i++) {
		if(j==k) {
			vMatrix[k](i,j) = m0(i,j)*1.3;
		}
		else {
			vMatrix[k](i,j) = m0(i,j)-m0(i,k)*0.1;
		}
	}
}

for(int k=4; k<8; k++) {
	for(int i=0; i<4; i++) {
		if(j==k) {
			vMatrix[k](i,j) = m0(i,j)*0.7;
		}
		else {
			vMatrix[k](i,j) = m0(i,j)+m0(i,k)*0.1;
		}
	}
}


for(int i=0; i<4< i++) {
	int need_to_check_neg_frac = 1;
	int n_negFraction = 0;
	int isSup1 = 0;

	//first check if there is no fraction > 1.
	//if there are, then you don't need to check the condition 'fraction>0'
	for(int j=0; j<4; j++) {
		if(vMatrix[0](i,j)>1) {
			need_to_check_neg_frac = 0;
			isSup1 = 1;
		}
	}
	if(isSup1 == 1) {
		setMaxFracTo1(vMatrix[0], i);
	}
	
// 	int index_tmp = 0;
// 	for(int j=0; j<4; j++) {
// 		if(vMatrix[0](i,j)>1) {
// 			need_to_check_neg_frac = 0;
// 			index_tmp = j;
// 			vMatrix[0](i,j) = 1.;
// 		}
// 	}
// 	for(int j=0; j<4; j++) {
// 		if(j != index_tmp) {
// 			vMatrix[0](i,j) = 0.;
// 		}
// 	}

	n_negFraction = getN_negFraction(vMatrix[0], i);
	while(n_negFraction != 0) {
		//first find the higher fraction
		int index_maxFraction = getIndex_maxFraction(vMatrix[0], i);

		//if n_negFraction = 1
		if(n_negFraction == 1) {
			int index_negFraction = getIndex_negFraction(vMatrix[0], i);
			int n_nullFraction = getN_nullFraction(vMatrix[0], i);
			if(n_nullFraction == 0) {
				for(int j=0; j<4; j++) {
					if(j!=index_maxFraction && j!=index_negFraction) {
						vMatrix[0](i,j) = vMatrix[0](i,j) - fabs(vMatrix[0](i,index_negFraction)/2.);
					}
				}
			}
			else if(n_nullFraction == 1) {
				int index_nullFraction = getIndex_nullFraction(vMatrix[0], i);
				for(int j=0; j<4; j++) {
					if(j!=index_maxFraction && j!=index_negFraction && j!=index_nullFraction) {
						vMatrix[0](i,j) = vMatrix[0](i,j) - fabs(vMatrix[0](i,index_negFraction));
					}
				}
			}
			else if(n_nullFraction == 2) {
				setMaxFracTo1(vMatrix[0], i);
			}
		}

		//if n_negFraction = 2
		else if(n_negFraction == 2) {
			int n_posFraction = 0;
			if(n_nullFraction == 0) {
				int index_negFraction1 = getIndex_negFraction(vMatrix[0], i);
				int index_negFraction2 = 0;
				for(int j=0; j<4; j++) {
					if(Matrix[0](i,j)<0. && j!=index_negFraction1) {
						index_negFraction2 = j;
					}
				}				
				
				for(int j=0; j<4; j++) {
					if(j!=index_maxFraction) {
						if(Matrix[0](i,j)>0.) {
						vMatrix[0](i,j) = vMatrix[0](i,j) - fabs(vMatrix[0](i,index_negFraction1)) - fabs(vMatrix[0](i,index_negFraction2));
						}
					}
				}
			}
			else if(n_nullFraction == 1) {
				setMaxFracTo1(vMatrix[0], i);
			}
		}

		//if n_negFraction = 3
		else if(n_negFraction == 3) {
			setMaxFracTo1(vMatrix[0], i, need_to_check_neg_frac);
		}

		n_negFraction = getN_negFraction(vMatrix[0], i);
	}



}



































