#include "ptBinning.h"

ptBinning::ptBinning()
{
	fillPtBins();
}

ptBinning::ptBinning(int isFor2DTagging)
{	
	if(isFor2DTagging == 1) {
		fill2DTaggingPtBins();
	}
	else {
		fillPtBins();
	}
}

ptBinning::~ptBinning()
{

}

int ptBinning::getPtBin(float pt)
{
	for(int i=0; i<mPtBins.size(); i++)
	{
		std::pair<float, float> bin = mPtBins[i];
		if (pt >= bin.first && pt < bin.second)
		{
			return i;
		}	
	}
	
	return -1;
}

size_t ptBinning::getSize()
{
	return mPtBins.size();
}

std::pair<float, float> ptBinning::getBinValue(int bin)  
{
	return mPtBins[bin];
}

float ptBinning::getBinValueInf(int bin)  
{
	return mPtBins[bin].first;
}

float ptBinning::getBinValueSup(int bin)  
{
	return mPtBins[bin].second;
}

std::string ptBinning::getName(int bin) {
	std::stringstream PtBinningName ;
	PtBinningName << "ptPhot_" << (int) getBinValueInf(bin) << "_" << (int) getBinValueSup(bin);
	return PtBinningName.str();
}


void ptBinning::fillPtBins() 
{
	mPtBins.push_back(std::make_pair(40., 50.));
	mPtBins.push_back(std::make_pair(50., 60.));
	mPtBins.push_back(std::make_pair(60., 75.));
	mPtBins.push_back(std::make_pair(75., 100.));
	mPtBins.push_back(std::make_pair(100., 125.));
	mPtBins.push_back(std::make_pair(125., 155.));
	mPtBins.push_back(std::make_pair(155., 180.));
	mPtBins.push_back(std::make_pair(180., 210.));
	mPtBins.push_back(std::make_pair(210., 250.));
	mPtBins.push_back(std::make_pair(250., 300.));
	mPtBins.push_back(std::make_pair(300., 350.));
	mPtBins.push_back(std::make_pair(350., 400.));
	mPtBins.push_back(std::make_pair(400., 500.));
	mPtBins.push_back(std::make_pair(500., 600.));
	mPtBins.push_back(std::make_pair(600., 800.));
}

void ptBinning::fill2DTaggingPtBins() 
{
  mPtBins.push_back(std::make_pair(40., 100.));
  mPtBins.push_back(std::make_pair(100., 200.));
  mPtBins.push_back(std::make_pair(200., 800.));

	//stage M2
/* 	mPtBins.push_back(std::make_pair(40., 100.));*/
 	/*mPtBins.push_back(std::make_pair(100., 800.));*/
}




