#include "alphaBinning.h"

alphaBinning::alphaBinning()
{
	fillAlphaBins();
}

alphaBinning::~alphaBinning()
{

}

int alphaBinning::getAlphaBin(float alpha)
{
	for(int i=0; i<mAlphaBins.size(); i++)
	{
		std::pair<float, float> bin = mAlphaBins[i];
		if (alpha >= bin.first && alpha < bin.second)
		{
			return i;
		}	
	}
	
	return -1;
}

size_t alphaBinning::getSize()
{
	return mAlphaBins.size();
}

std::pair<float, float> alphaBinning::getBinValue(int bin)  
{
	return mAlphaBins[bin];
}

float alphaBinning::getBinValueInf(int bin)  
{
	return mAlphaBins[bin].first;
}

float alphaBinning::getBinValueSup(int bin)  
{
	return mAlphaBins[bin].second;
}

void alphaBinning::fillAlphaBins() 
{
	mAlphaBins.push_back(std::make_pair(0., 0.05));
	mAlphaBins.push_back(std::make_pair(0.05, 0.10));
	mAlphaBins.push_back(std::make_pair(0.10, 0.15));
	mAlphaBins.push_back(std::make_pair(0.15, 0.20));
	mAlphaBins.push_back(std::make_pair(0.20, 0.25));
	mAlphaBins.push_back(std::make_pair(0.25, 0.30));
}




