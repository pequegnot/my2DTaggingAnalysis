#ifndef DEF_ALPHABINNING
#define DEF_ALPHABINNING

#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <utility>

class alphaBinning
{
    public:

    alphaBinning();
    ~alphaBinning();
    void fillAlphaBins();
    int getAlphaBin(float alpha);
    size_t getSize();
    std::pair<float, float> getBinValue(int bin);
    float getBinValueInf(int bin);
    float getBinValueSup(int bin);


    private:
    std::vector<std::pair<float,float> > mAlphaBins;


};

#endif


