#include <iostream>
#include "file_io.hpp"
#include "mask.hpp"

void generateLineSpace(vector<float>& msk, int sizeX, int sizeY, LineSpace& lineSpace) {
  int period = lineSpace.lineWidth + lineSpace.spaceWidth;
  if(lineSpace.isHorizontal){
	for (int y = 0; y < sizeY; ++y) {
	  for (int x = 0; x < sizeX; ++x) {
		if ((y % period) < lineSpace.lineWidth) {
		  msk[y * sizeX + x] = 1.0; // line
		} else {
		  msk[y * sizeX + x] = 0.0; // space
		}
	  }
	}
  
  }else{
	for (int y = 0; y < sizeY; ++y) {
	  for (int x = 0; x < sizeX; ++x) {
		if ((x % period) < lineSpace.lineWidth) {
		  msk[y * sizeX + x] = 1.0; // line
		} else {
		  msk[y * sizeX + x] = 0.0; // space
		}
	  }
	}
  }
}

void createMask(vector<float>& msk, int Lx, int Ly, int maskSizeX, int maskSizeY,
				string maskType, LineSpace& lineSpace,
				string maskInFile, stringstream& ss_mskParam, float dose) {
  vector<float> inMsk(maskSizeX * maskSizeY);
  
  ss_mskParam <<"Mask Period: " << Lx << " x " << Ly << "\nMask Size: " << maskSizeX << " x " << maskSizeY << "\nMask Type: " << maskType << "\n";
 
  if (maskType == "LineSpace") {
    ss_mskParam << "Line Width: " << lineSpace.lineWidth << "\nSpace Width: " << lineSpace.spaceWidth << "\n";
    generateLineSpace(inMsk, maskSizeX, maskSizeY, lineSpace);
  }
  else if (maskType == "Import") {
    // Read mask data
    cout << "Read Mask: " << maskInFile <<endl;
    string extension = maskInFile.substr(maskInFile.find_last_of(".") + 1);
    if(extension == "bin"){
	  readFloatArrayFromBinary(maskInFile, inMsk, maskSizeX * maskSizeY);
    } else {
	  readFloatArrayFromTxt(maskInFile, inMsk, maskSizeX * maskSizeY);
	
    }
  }
  int difYh = (Ly - maskSizeY) / 2;
  int difXh = (Lx - maskSizeX) / 2;
  for (int y = 0; y < maskSizeY; y++) {
    for (int x = 0; x < maskSizeX; x++) {
	  msk[(y + difYh) * Lx + x + difXh] = inMsk[y * maskSizeX + x] * dose;
    }
  }
}
