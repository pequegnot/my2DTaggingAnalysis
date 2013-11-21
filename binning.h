#pragma once

#include<string>
#include<math.h>

std::string getFlavourBinName(int flavourBin) {
//****************************************************************************************
//                   !!!! WARNING !!!!
//        if you change the correspondance binName-binNumber
//        do not forget to change it in the getFlavourBin function of common.h
//****************************************************************************************
		if(flavourBin==0) {
			return "uds";
		}
		else if(flavourBin==1) {
			return "g";
		}
		else if(flavourBin==2) {
			return "c";
		}
		else if(flavourBin==3) {
			return "b";
		}
		else if(flavourBin==4) {
			return "noMatched";
		}
		else if(flavourBin==5) {
			return "all";
		}
		return ""; 
}

std::string getZoneBinName(int zoneBin) {
		if(zoneBin==0) {
			return "zone1Q";
		}
		else if(zoneBin==1) {
			return "zone2G";
		}
		else if(zoneBin==2) {
			return "zone3C";
		}
		else if(zoneBin==3) {
			return "zone4B";
		}
		else if(zoneBin==4) {
			return "zone5QG";
		}
		else if(zoneBin==5) {
			return "zone6LC";
		}
		return ""; 
}

std::string getTaggerBinName(int taggerBin) {
		if(taggerBin==0) {
			return "CSV";
		}
		else if(taggerBin==1) {
			return "QGL";
		}
		return ""; 
}

int getFlavourNumber () {
	return 6;
}

int getZoneNumber () {
	return 6;
}

int getTaggerNumber () {
	return 2;
}



