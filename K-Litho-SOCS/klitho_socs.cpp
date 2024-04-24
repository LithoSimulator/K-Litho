#include <iostream>
#include <complex> 
#include <string>
#include <fstream>
#include <chrono>
#include "fftw3.h"
#include "file_io.hpp"
#include <vector>

using namespace std;
using namespace std::chrono;
typedef std::complex<double> Complex;

float getMax(vector<float>& data, int size){
  float val, max=-1e+30;
  
  for(int i=0; i<size; i++){
	val = data[i];
	if(max<val){max= val;}
  }
  return max;
}

float getMin(vector<float>& data, int size){
  float val, min=1e+30;
  
  for(int i=0; i<size; i++){
	val = data[i];
	if(min>val){min= val;}
  }
  return min;
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
void FT_r2c(vector<float>& in, vector<Complex>& out, int sizeX, int sizeY, fftw_plan &fwd, vector<double>& datar, vector<Complex>& datac){
  myShift(in, datar,sizeX,sizeY,false,false);
 
  fftw_execute(fwd);
  int sizeXh=sizeX/2;
  int wt=(sizeX/(2)+(1));
  int AY=(sizeY+1)/2;
  int BY=sizeY-AY;

  vector<Complex> R(wt*sizeY, 0);
  vector<Complex> L(wt*sizeY, 0);
  
  for(int y=0, y2=AY; y<BY; y++, y2++){
	for(int x=0; x<wt; x++){
	  R[y*wt+x]=datac[y2*wt+x];
	}
  }
  for(int y=BY, y2=0; y<sizeY; y++, y2++){
	for(int x=0; x<wt; x++){
	  R[y*wt+x]=datac[y2*wt+x];
	}
  }
  if(sizeY%2==0){
	for(int x=0, x2=wt-1; x<wt; x++, x2--){
	  L[x]=conj(R[x2]);
	}
    
	for(int y=1, y2=sizeY-1; y<sizeY; y++, y2--){
	  for(int x=0, x2=wt-1; x<wt; x++, x2--){
		L[y*wt+x]=conj(R[y2*wt+x2]);
	  }
	}
  }else{
	for(int y=0, y2=sizeY-1; y<sizeY; y++, y2--){
	  for(int x=0, x2=wt-1; x<wt; x++, x2--){
		L[y*wt+x]=conj(R[y2*wt+x2]);
	  }
	}
  }
  for(int y=0; y<sizeY; y++){
	for(int x=0; x<wt; x++){
	  out[y*sizeX+x]=L[y*wt+x];
	}
  }
  for(int y=0; y<sizeY; y++){
	for(int x=sizeXh, x2=0; x<sizeX; x++, x2++){
	  out[y*sizeX+x]=R[y*wt+x2];
	}
  }
  
  for(int i=0; i<sizeX*sizeY; i++){
	out[i] /=(sizeX*sizeY);
  }
}

//Fast Fourier transform from complex data to real data
void FT_c2r(vector<Complex>& in, vector<float>& out, int sizeX, int sizeY, fftw_plan &bwd, vector<Complex>& datac, vector<double>& datar){
  int wt=(sizeX/(2)+(1));
  int sx, sy;
  int xh = (sizeX+1)/2;
  int yh = (sizeY+1)/2;
  for (int y=0; y<sizeY; y++) {
    for(int x=0; x<wt; x++) {
      datac[y*wt+x] = 0;
    }
  }
  for (int y=0; y<sizeY; y++) {
    for(int x=xh; x<sizeX; x++) {
      sx = (x+xh) % sizeX;
      sy = (y+yh) % sizeY;
      datac[sy*wt+sx] = in[y*sizeX+x];
    }
  }
  if(sizeX%2==0){
	vector<Complex> R(sizeY);
	R[0]=conj(in[0]);
	for(int y=1, y2=sizeY-1; y<sizeY; y++, y2--){
	  R[y]=conj(in[y2*sizeX]);
	}
	for(int y=0; y<sizeY; y++){
	  datac[y*wt+wt-1]=R[y];
	}
  }

  fftw_execute(bwd);
 
  myShift(datar,out,sizeX,sizeY,false,false);
}

//Fourier interpolation
template <typename T1, typename T2>
void FI(vector<T1>& in, int in_sizeX, int in_sizeY, vector<T2>& out, int out_sizeX, int out_sizeY, fftw_plan &R2C_plan, vector<double>& R2CDataReal, vector<Complex>& R2CDataComplex, fftw_plan &C2R_plan, vector<Complex>& C2RDataComplex, vector<double>& C2RDataReal){
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

bool extractCentralRegion(vector<float>& in, int inSizeX, int inSizeY, vector<float>& out, int outSizeX, int outSizeY) {
  if (outSizeX > inSizeX || outSizeY > inSizeY) {
	std::cerr << "Error: outSizeX and outSizeX must be less than or equal to inSizeX and inSizeY respectively." << std::endl;
	return false;
  }
  // Compute the starting indices for the center of the in array
  int startIdxX = inSizeX / 2 - outSizeX / 2;
  int startIdxY = inSizeY / 2 - outSizeX / 2;

  // Copy the central outSizeX * outSizeY data from in to out
  for (int y = 0; y < outSizeY; ++y) {
	for (int x = 0; x < outSizeX; ++x) {
	  out[y * outSizeX + x] = in[(startIdxY + y) * inSizeX + (startIdxX + x)];
	}
  }
  return true;
}

void calcSOCS(vector<float>& image, int Nx, int Ny, int nk, const vector<vector<Complex>>& krns, const vector<float>& scales, const vector<Complex>& msk, int Lx, int Ly, fftw_plan& LxyC2R_plan, vector<Complex>& complexDataLxy, vector<double>& realDataLxy) {
  const int sizeX = 4 * Nx + 1;
  const int sizeY = 4 * Ny + 1;
  const int difX = sizeX - (2 * Nx + 1);
  const int difY = sizeY - (2 * Ny + 1);
  const int Lxh = Lx / 2;
  const int Lyh = Ly / 2;

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
  cout << "*****************************************************************" << endl;
  cout << "* Lithography Simulation initiated..." << endl;
  cout << "*****************************************************************" << endl;
  auto start = system_clock::now();

  const float TARGET_INTENSITY = 0.225;

  // Check the number of command-line arguments
  if (argc != 7) {
    cerr << "Error: Incorrect number of command-line arguments." << endl;
    return 1;
  }
  
  const int Lx = stoi(argv[1]);
  const int Ly = stoi(argv[2]);
  const int Lxy = Lx * Ly;
  const int maskSizeX = stoi(argv[3]);
  const int maskSizeY = stoi(argv[4]);
  const string maskInFile = argv[5];
  const string KrnDir = argv[6];
  const float dose = 1.0;
  int Nx, Ny, nk;
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
  vector<float> inMsk;
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
  vector<vector<Complex>> krns(nk, vector<Complex>((2 * Nx + 1) * (2 * Ny + 1)));
  vector<float> tmpDataReal, tmpDataImag;
  for (int i = 0; i < nk; i++) {
	string krnRealInFile = KrnDir + "krn_" + to_string(i) + "_r.bin";
	string krnImagInFile = KrnDir + "krn_" + to_string(i) + "_i.bin";
	cout << "Read Kernel " << i <<" : " << krnRealInFile <<endl;
	if (!readFloatArrayFromBinary(krnRealInFile, tmpDataReal, (2 * Nx + 1) * (2 * Ny + 1))) {
	  cerr << "Error: Failed to read kernel " << i << " from file: " << krnRealInFile <<endl; 
	  return 1;
	}
	cout << "Read Kernel " << i <<" : " << krnImagInFile <<endl;
	if (!readFloatArrayFromBinary(krnImagInFile, tmpDataImag, (2 * Nx + 1) * (2 * Ny + 1))) {
	  cerr << "Error: Failed to read kernel " << i << " from file: " << krnImagInFile <<endl; 
	  return 1;
	}
	for (int j = 0; j <(2 * Nx + 1) * (2 * Ny + 1) ; j++) {
	  krns[i][j]=Complex(tmpDataReal[j],tmpDataImag[j]);
	}
  }
  
  // Output simulation parameters
  cout << endl << "Simulation Parameters:" << endl;
  cout << "-------------------" << endl;
  cout << "Mask Period: " << Lx << " x " << Ly << endl;
  cout << "Mask Size: " << maskSizeX << " x " << maskSizeY << endl;
  cout << "Nx: " << Nx << ",  Ny: " << Ny << endl;
  cout << "Number of Kernels: " << nk << endl << endl;
  
  // Read mask data
  cout << "Read Mask: " << maskInFile <<endl;
  string extension = maskInFile.substr(maskInFile.find_last_of(".") + 1);
  if (extension == "bin") {
    if (!readFloatArrayFromBinary(maskInFile, inMsk, maskSizeX * maskSizeY)) {
      cerr << "Error: Failed to read mask file: " << maskInFile << endl;
      return 1;
    } 
  } else if (extension == "txt") {
    if (!readFloatArrayFromText(maskInFile, inMsk, maskSizeX * maskSizeY)) {
      cerr << "Error: Failed to read mask file: " << maskInFile << endl;
      return 1;
    } 
  } else {
    cerr << "Error: Unsupported file extension: " << extension << endl;
    return 1;
  }

  max = getMax(inMsk, maskSizeX * maskSizeY);
  min = getMin(inMsk, maskSizeX * maskSizeY);
  
  if (min < 0 || max > 1) {
    cerr << "Error: Mask values must be in the range [0, 1]." << endl;
    return 1;
  }
  writeFloatArrayToPNGWithGrayscale(maskSizeX, maskSizeY, inMsk, min, max, "./out/mask.png");

  const int difYh = (Ly - maskSizeY) / 2;
  const int difXh = (Lx - maskSizeX) / 2;
  for (int y = 0; y < maskSizeY; ++y) {
    for (int x = 0; x < maskSizeX; ++x) {
      if (inMsk[y * maskSizeX + x] != 0) {
        msk[(y + difYh) * Lx + x + difXh] = inMsk[y * maskSizeX + x] * dose;
      }
    }
  }
  
  // Calculate FFT of mask
  cout << endl << "-----------------------------------------------------------------" << endl;
  cout << "Calculating FFT Mask..." << endl;
  FT_r2c(msk, mskf, Lx, Ly, LxyR2C_plan, realDataLxy, complexDataLxy);
  
  cout << endl << "-----------------------------------------------------------------" << endl;
  cout << "Calc Image..." << endl;
  calcSOCS(img, Nx, Ny, nk, krns, scales, mskf, Lx, Ly, LxyC2R_plan, complexDataLxy, realDataLxy);
 
  max = getMax(img, Lx * Ly);
  min = getMin(img, Lx * Ly);
  cout << "Max. Intensity: " << max << endl;
  cout << "Min. Intensity: " << min << endl;
  
  if (!writeFloatArrayToBinary("./out/image.bin", img, Lx * Ly)) {
    cerr << "Error: Failed to write image data to file: ./out/image.bin" << endl;
    return 1;
  }
  writeFloatArrayToPNGWithContinuousColor(Lx, Ly, img, "./out/image.png", min, max);
  writeFloatArrayToPNGWith7Colors(Lx, Ly, img, "./out/image_s.png", 0, TARGET_INTENSITY);

  //output image only mask region
  // vector<float> imgMaskSize(maskSizeX * maskSizeY);
  // extractCentralRegion(img, Lx, Ly, imgMaskSize, maskSizeX, maskSizeY);
  // if (!writeFloatArrayToBinary("./out/image.bin", imgMaskSize, maskSizeX*maskSizeY)) {
  //   cerr << "Error: Failed to write image data to file: ./out/image.bin" << endl;
  //   return 1;
  // }
  // max = getMax(imgMaskSize, maskSizeX * maskSizeY);
  // min = getMin(imgMaskSize, maskSizeX * maskSizeY);
  // writeFloatArrayToPNGWithContinuousColor(maskSizeX, maskSizeY, imgMaskSize, "./out/image.png", min, max);
  // writeFloatArrayToPNGWith7Colors(maskSizeX, maskSizeY, imgMaskSize, "./out/image_s.png", 0, TARGET_INTENSITY);
  
  
  auto end = system_clock::now();
  auto dur = end - start;
  auto microsec = duration_cast<chrono::microseconds>(dur).count();
  cout << endl << "-----------------------------------------------------------------" << endl;
  cout << "Simulation completed." << endl;
  cout << "Total Execution Time: " << microsec * 1.0e-06 << " sec." << endl;
  cout << "*****************************************************************" << endl;

  // Destroy FFTW plans
  fftw_destroy_plan(LxyC2R_plan);
  fftw_destroy_plan(LxyR2C_plan);

  return 0;
}
