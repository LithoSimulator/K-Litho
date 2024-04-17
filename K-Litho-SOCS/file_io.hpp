#pragma once
#include <vector>
#include <complex>
#include <string>
using namespace std;
typedef std::complex<double> Complex;

/* #include "stb_image_write.h" */

/* #define STB_IMAGE_WRITE_IMPLEMENTATION */
/* #include "stb_image_write.h" */

using namespace std;
typedef std::complex<double> Complex;

float writeFloatArrayToPNGWithGrayscale(int sizeX, int sizeY, vector<float>& data, float min, float max, string filename) ;

bool writeFloatArrayToPNGWith7Colors(int sizeX, int sizeY, vector<float>& data, string filename, float min, float max) ;

bool writeFloatArrayToPNGWithContinuousColor(int sizeX, int sizeY, vector<float>& data, string filename, float min, float max) ;


bool writeFloatArrayToBinary(const string& filename, const vector<float>& data, int size) ;

bool writeComplexArrayToBinary(string realFile, string imagFile, vector<Complex>& data, int size) ;

bool writeFloatArrayToTxt(string filename, vector<float>& data, int sizeX, int sizeY) ;

bool readFloatArrayFromText(const string& filename, vector<float>& data, int size) ;

bool readFloatArrayFromBinary(const string& filename, vector<float>& data, int size) ;

void writeKernelInfo(string filename, int krnSizeX, int krnSizeY, int nk, vector<float>& scales) ;

void readKernelInfo(string filename, int& krnSizeX, int& krnSizeY, int& nk, vector<float>& scales) ;
