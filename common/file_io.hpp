#ifndef FILE_IO_HPP
#define FILE_IO_HPP

#include <vector>
#include <complex>
#include <string>
using namespace std;
typedef complex<double> Complex;

void writeFloatArrayToPNGWithGrayscale(int sizeX, int sizeY, vector<float>& data, float min, float max, string filename) ;

void writeFloatArrayToPNGWith7Colors(int sizeX, int sizeY, vector<float>& data, string filename, float min, float max) ;

void writeFloatArrayToPNGWithContinuousColor(int sizeX, int sizeY, vector<float>& data, string filename, float min, float max) ;

void writeFloatArrayToBinary( string filename,  vector<float>& data, int size) ;

void writeComplexArrayToBinary(string realFile, string imagFile, vector<Complex>& data, int size) ;

void writeFloatArrayToTxt(string filename, vector<float>& data, int sizeX, int sizeY) ;

void readFloatArrayFromTxt( string filename, vector<float>& data, int size) ;

void readFloatArrayFromBinary( string filename, vector<float>& data, int size) ;

void writeKernelInfo(string filename, int krnSizeX, int krnSizeY, int nk, vector<float>& scales) ;

void readKernelInfo(string filename, int& krnSizeX, int& krnSizeY, int& nk, vector<float>& scales) ;

#endif
