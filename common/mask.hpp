#ifndef MASK_HPP
#define MASK_HPP

#include <string>
#include <vector>
#include <sstream>

using namespace std;

struct LineSpace {
  bool isHorizontal;
  int lineWidth;
  int spaceWidth;
};

void generateLineSpace(vector<float>& msk, int sizeX, int sizeY, LineSpace& lineSpace) ;

void createMask(vector<float>& msk, int Lx, int Ly, int maskSizeX, int maskSizeY,
				string maskType, LineSpace& lineSpace,
				string maskInFile, stringstream& ss_mskParam, float dose) ;

#endif
