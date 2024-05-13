#include <iostream>
#include <complex> 
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <chrono>
#include "fftw3.h"
#include "file_io.hpp"
#include "mask.hpp"

using namespace std;
using namespace std::chrono;
typedef std::complex<double> Complex;

float getMax(vector<float>& data){
  float val, max=-1e+30;
  for (int i=0; i< (int) data.size(); i++) {
	val = data[i];
	if (max < val){
	  max = val;
	}
  }
  return max;
}
float getMin(vector<float>& data){
  float val, min=1e+30;
  for (int i=0; i< (int) data.size(); i++) {
	val = data[i];
	if (min > val) {
	  min = val;
	}
  }
  return min;
}

bool extractCentralRegion(vector<float>& in, int inSizeX, int inSizeY, vector<float>& out, int outSizeX, int outSizeY) {
  if (outSizeX > inSizeX || outSizeY > inSizeY) {
	cerr << "Error: outSizeX and outSizeY must be less than or equal to inSizeX and inSizeY respectively." << endl;
	return false;
  }
  // Compute the starting indices for the center of the in array
  int startIdxX = inSizeX / 2 - outSizeX / 2;
  int startIdxY = inSizeY / 2 - outSizeY / 2;

  // Copy the central outSizeX * outSizeY data from in to out
  for (int y = 0; y < outSizeY; ++y) {
	for (int x = 0; x < outSizeX; ++x) {
	  out[y * outSizeX + x] = in[(startIdxY + y) * inSizeX + (startIdxX + x)];
	}
  }
  return true;
}

template <typename T1, typename T2>
void myShift(vector<T1>& in, vector<T2>& out, int sizeX, int sizeY, bool shiftTypeX, bool shiftTypeY) {
  int sx, sy;
  int xh = shiftTypeX ? sizeX / 2 : (sizeX + 1) / 2;
  int yh = shiftTypeY ? sizeY / 2 : (sizeY + 1) / 2;

  for (int y = 0; y < sizeY; y++) {
	for (int x = 0; x < sizeX; x++) {
	  sy = (y + yh) % sizeY;
	  sx = (x + xh) % sizeX;
	  out[sy * sizeX + sx] = in[y * sizeX + x];
	}
  }
}

//Fast Fourier transform from real data to complex data
void FT_r2c(vector<float>& in, vector<Complex>& out, int sizeX, int sizeY,
			fftw_plan &fwd, vector<double>& datar, vector<Complex>& datac){
  myShift(in, datar, sizeX, sizeY, false, false);
  fftw_execute(fwd);
  int sizeXh = sizeX / 2;
  int wt = (sizeX / 2) + 1;
  int AY = (sizeY + 1) / 2;
  int BY = sizeY - AY;
  vector<Complex> R(wt * sizeY, 0);
  vector<Complex> L(wt * sizeY, 0);
  for (int y = 0, y2 = AY; y < BY; y++, y2++) {
	for (int x = 0; x < wt; x++) {
	  R[y * wt + x] = datac[y2 * wt + x];
	}
  }
  for (int y = BY, y2 = 0; y < sizeY; y++, y2++) {
	for (int x = 0; x < wt; x++) {
	  R[y * wt + x] = datac[y2 * wt + x];
	}
  }
  if (sizeY % 2 == 0) {
	for (int x = 0, x2 = wt - 1; x < wt; x++, x2--) {
	  L[x] = conj(R[x2]);
	}
	for (int y = 1, y2 = sizeY - 1; y < sizeY; y++, y2--) {
	  for (int x = 0, x2 = wt - 1; x < wt; x++, x2--) {
        L[y * wt + x] = conj(R[y2 * wt + x2]);
	  }
	}
  }
  else {
	for (int y = 0, y2 = sizeY - 1; y < sizeY; y++, y2--) {
	  for (int x = 0, x2 = wt - 1; x < wt; x++, x2--) {
		L[y * wt + x] = conj(R[y2 * wt + x2]);
	  }
	}
  }
  for (int y = 0; y < sizeY; y++) {
	for (int x = 0; x < wt; x++) {
	  out[y * sizeX + x] = L[y * wt + x];
	}
  }
  for (int y = 0; y < sizeY; y++) {
	for (int x = sizeXh, x2 = 0; x < sizeX; x++, x2++) {
	  out[y * sizeX + x] = R[y * wt + x2];
	}
  }
  for (int i = 0; i < sizeX * sizeY; i++) {
	out[i] /= (sizeX * sizeY);
  }
}

//Fast Fourier transform from complex data to real data
void FT_c2r(vector<Complex>& in, vector<float>& out, int sizeX, int sizeY,
			fftw_plan &bwd, vector<Complex>& datac, vector<double>& datar){
  int wt = (sizeX / 2) + 1;
  int sx, sy;
  int xh = (sizeX + 1) / 2;
  int yh = (sizeY + 1) / 2;
  for (int y = 0; y < sizeY; y++) {
	for (int x = 0; x < wt; x++) {
	  datac[y * wt + x] = 0;
	}
  }
  for (int y = 0; y < sizeY; y++) {
	for (int x = xh-1; x < sizeX; x++) {
	  sx = (x + xh) % sizeX;
	  sy = (y + yh) % sizeY;
	  datac[sy * wt + sx] = in[y * sizeX + x];
	}
  }
  if (sizeX % 2 == 0) {
	vector<Complex> R(sizeY);
	R[0] = conj(in[0]);
	for (int y = 1, y2 = sizeY - 1; y < sizeY; y++, y2--) {
	  R[y] = conj(in[y2 * sizeX]);
	}

	for (int y = 0; y < sizeY; y++) {
	  datac[y * wt + wt - 1] = R[y];
	}
  }
  
  fftw_execute(bwd);
  myShift(datar, out, sizeX, sizeY, false, false);
}

//Fourier interpolation
template <typename T1, typename T2>
void FI(vector<T1>& in, int in_sizeX, int in_sizeY,
		vector<T2>& out, int out_sizeX, int out_sizeY,
		fftw_plan &R2C_plan, vector<double>& R2CDataReal, vector<Complex>& R2CDataComplex,
		fftw_plan &C2R_plan, vector<Complex>& C2RDataComplex, vector<double>& C2RDataReal){
  int difX= out_sizeX-in_sizeX;
  int difY= out_sizeY-in_sizeY;
  if(difX<0 || difY<0){
    cout<<"ERROR in FI()"<<endl;
    exit(1);
  }
  int in_sizeXY= in_sizeX*in_sizeY;
  fill(R2CDataReal.begin(), R2CDataReal.end(), 0);
  fill(R2CDataComplex.begin(), R2CDataComplex.end(), 0);
  myShift(in, R2CDataReal, in_sizeX,in_sizeY,false,false);
  fftw_execute(R2C_plan);
  for(int i=0; i<in_sizeY*(in_sizeX/2+1); i++){
    R2CDataComplex[i]/=(in_sizeXY);
  }
  fill(C2RDataComplex.begin(), C2RDataComplex.end(), 0);
  for(int y=0; y<in_sizeY/2+1; y++){
    for(int x=0; x<in_sizeX/2+1; x++){
	  C2RDataComplex[(y)*(out_sizeX/(2)+(1))+(x)]=R2CDataComplex[y*(in_sizeX/2+1)+x];
    }
  }
  
  for(int y=in_sizeY/2+1; y< in_sizeY; y++){
    for(int x=0; x<in_sizeX/2+1; x++){
	  C2RDataComplex[(y+difY)*(out_sizeX/(2)+(1))+(x)]=R2CDataComplex[y*(in_sizeX/2+1)+x];
    }
  }
  fill(C2RDataReal.begin(), C2RDataReal.end(), 0);
  fftw_execute(C2R_plan);
  myShift(C2RDataReal, out, out_sizeX,out_sizeY,false,false);
  
}

void calcSOCS(vector<float>& image, vector<Complex>& msk, int Lx, int Ly,
			  vector<vector<Complex>>& krns, vector<float>& scales, int nk, int Nx, int Ny,
			  fftw_plan& LxyC2R_plan, vector<Complex>& complexDataLxy, vector<double>& realDataLxy) {
  int sizeX = 4 * Nx + 1;
  int sizeY = 4 * Ny + 1;
  int difX = sizeX - (2 * Nx + 1);
  int difY = sizeY - (2 * Ny + 1);
  int Lxh = Lx / 2;
  int Lyh = Ly / 2;

  vector<Complex> bwdDataIn(sizeX * sizeY, Complex(0, 0));
  vector<Complex> bwdDataOut(sizeX * sizeY, Complex(0, 0));
  fftw_plan bwd_plan = fftw_plan_dft_2d(sizeY, sizeX, reinterpret_cast<fftw_complex*>(bwdDataIn.data()), reinterpret_cast<fftw_complex*>(bwdDataOut.data()), FFTW_BACKWARD, FFTW_ESTIMATE);

  vector<double> realDataNxy(sizeX * sizeY);
  vector<Complex> complexDataNxy((sizeX / 2 + 1) * sizeY);
  fftw_plan NxyR2C = fftw_plan_dft_r2c_2d(sizeY, sizeX, realDataNxy.data(), reinterpret_cast<fftw_complex*>(complexDataNxy.data()), FFTW_ESTIMATE);
  
  vector<double> tmpImg(sizeX * sizeY, 0);
  vector<double> tmpImgp(sizeX * sizeY, 0);

  // Compute the optical image by SOCS
  for (int k = 0; k < nk; ++k) {
    for (int y = -Ny; y <= Ny; ++y) {
      for (int x = -Nx; x <= Nx; ++x) {
        bwdDataIn[(difY + y + Ny) * sizeX + difX + x + Nx] = krns[k][(y + Ny) * (2 * Nx + 1) + x + Nx] * msk[(y + Lyh) * Lx + x + Lxh];
      }
    }
	
    // Inverse Fourier transform of the product of kernel and mask
    fftw_execute(bwd_plan);

    for (int i = 0; i < sizeX * sizeY; ++i) {
      tmpImg[i] += scales[k] * (bwdDataOut[i].real() * bwdDataOut[i].real() + bwdDataOut[i].imag() * bwdDataOut[i].imag());
    }
  }

  myShift(tmpImg, tmpImgp, sizeX, sizeY, true, true);

  // Apply Fourier interpolation to enlarge the optical image
  FI(tmpImgp, sizeX, sizeY, image, Lx, Ly, NxyR2C, realDataNxy, complexDataNxy, LxyC2R_plan, complexDataLxy, realDataLxy);

  // Destroy the FFTW plans
  fftw_destroy_plan(bwd_plan);
  fftw_destroy_plan(NxyR2C);
}

int main(int argc, char* argv[]) {
  cout << "************************************" << endl;
  cout << "* Lithography Simulation initiated..." << endl;
  cout << "************************************" << endl;

  auto start = system_clock::now();

  // Check the number of command-line arguments
  if (argc != 11) {
	cerr << "Error: Invalid number of arguments: " << argc << ". Expected 11 arguments." << endl;
    return 1;
  }
  // Period of mask
  int Lx = stoi(argv[1]);
  int Ly = stoi(argv[2]);
  // Size of mask
  int maskSizeX = stoi(argv[3]);
  int maskSizeY = stoi(argv[4]);
  string maskType = argv[5];
  LineSpace lineSpace;
  // Parameters for the line and space pattern mask
  lineSpace.lineWidth = stoi(argv[6]);
  lineSpace.spaceWidth = stoi(argv[7]);
  lineSpace.isHorizontal = stoi(argv[8]);
  // Input file for the mask
  string maskInFile = argv[9];
  // Directory containing the kernel files
  string KrnDir = argv[10];
  // Other parameters
  int Lxy = Lx * Ly;
  float dose = 1.0;
  int Nx, Ny, nk;
  float TARGET_INTENSITY = 0.225;

  // Validate input parameters
  if (Lx <= 0 || Ly <= 0 || maskSizeX <= 0 || maskSizeY <= 0 || Lx < maskSizeX || Ly < maskSizeY){ 
    cerr << "Error: Invalid input parameters. Please check the values." << endl;
    return 1;
  }

  vector<Complex> complexDataLxy(Ly * (Lx / 2 + 1));
  vector<double> realDataLxy(Lxy);
  
  fftw_plan LxyC2R_plan = fftw_plan_dft_c2r_2d(Ly, Lx, reinterpret_cast<fftw_complex*>(complexDataLxy.data()), realDataLxy.data(), FFTW_ESTIMATE);
  fftw_plan LxyR2C_plan = fftw_plan_dft_r2c_2d(Ly, Lx, realDataLxy.data(), reinterpret_cast<fftw_complex*>(complexDataLxy.data()), FFTW_ESTIMATE);

  vector<Complex> mskf(Lxy);
  vector<float> msk(Lxy, 0);
  vector<float> img(Lxy);
  float min, max;
  int krnSizeX=0, krnSizeY=0;
 
  vector<float> scales;
  // Read eigenvalues
  readKernelInfo(KrnDir + "kernel_info.txt",krnSizeX, krnSizeY, nk, scales);
  if (krnSizeX % 2 == 0 || krnSizeY % 2 == 0) {
	cerr << "Error: Kernel size must be odd." << endl;
	return 1;
  }
  Nx = (krnSizeX - 1)/2;
  Ny = (krnSizeY - 1)/2;
  
  // Read kernels
  vector<vector<Complex>> krns(nk, vector<Complex>(krnSizeX * krnSizeY));
  vector<float> krnReal, krnImag;
  for (int i = 0; i < nk; i++) {
	string krnRealInFile = KrnDir + "krn_" + to_string(i) + "_r.bin";
	string krnImagInFile = KrnDir + "krn_" + to_string(i) + "_i.bin";
	cout << "Read Kernel " << i <<" : " << krnRealInFile <<endl;
	readFloatArrayFromBinary(krnRealInFile, krnReal, krnSizeX * krnSizeY);
	cout << "Read Kernel " << i <<" : " << krnImagInFile <<endl;
	readFloatArrayFromBinary(krnImagInFile, krnImag, krnSizeX * krnSizeY);
	for (int j = 0; j <krnSizeX * krnSizeY ; j++) {
	  krns[i][j]=Complex(krnReal[j],krnImag[j]);
	}
  }
  
  // Create mask
  stringstream ss_mskParam;
  createMask(msk, Lx, Ly, maskSizeX, maskSizeY, maskType, lineSpace, maskInFile, ss_mskParam, dose);

  // Output mask to PNG file
  writeFloatArrayToPNGWithGrayscale(Lx, Ly, msk, getMin(msk), getMax(msk), "./out/mask.png");

  
  // Output simulation parameters
  cout << endl << "Simulation Parameters:" << endl;
  cout << "-------------------" << endl;
  cout << ss_mskParam.str();
  cout << "Nx: " << Nx << ",  Ny: " << Ny << endl;
  cout << "Number of Kernels: " << nk << endl << endl;
  
  // Calculate FFT of mask
  cout << endl << "------------------------------------" << endl;
  cout << "Calculating FFT Mask..." << endl;
  FT_r2c(msk, mskf, Lx, Ly, LxyR2C_plan, realDataLxy, complexDataLxy);

  // Calculate optical image
  cout << endl << "------------------------------------" << endl;
  cout << "Calc Image..." << endl;
  calcSOCS(img, mskf, Lx, Ly, krns, scales, nk, Nx, Ny, LxyC2R_plan, complexDataLxy, realDataLxy);
 
  max = getMax(img);
  min = getMin(img);
  cout << "Max. Intensity: " << max << endl;
  
  // Output optical image to binary file
  writeFloatArrayToBinary("./out/image.bin", img, Lx * Ly);
  // Output optical image to PNG file
  writeFloatArrayToPNGWithContinuousColor(Lx, Ly, img, "./out/image.png", min, max);
  writeFloatArrayToPNGWith7Colors(Lx, Ly, img, "./out/image_s.png", 0, TARGET_INTENSITY);
  
  auto end = system_clock::now();
  auto dur = end - start;
  auto microsec = duration_cast<chrono::microseconds>(dur).count();
  cout << endl << "------------------------------------" << endl;
  cout << "Simulation completed." << endl;
  cout << "Total Execution Time: " << microsec * 1.0e-06 << " sec." << endl;
  cout << "************************************" << endl;

  // Destroy FFTW plans
  fftw_destroy_plan(LxyC2R_plan);
  fftw_destroy_plan(LxyR2C_plan);

  return 0;
}
