#include <iostream>
#include <complex> 
#include <string>
#include <chrono>
#include <vector>
#include <sstream>
#include "fftw3.h"
#include "file_io.hpp"
#include "source.hpp"
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

float getOuterSigma(vector<float>& src, int srcSize){
  // Find the distance of the coordinate farthest from the center coordinate
  int center = (srcSize - 1) / 2;
  float max_distance = 0.0;
  for (int y = 0; y < srcSize; y++) {
	for (int x = 0; x < srcSize; x++) {
	  if (src[y*srcSize+x] != 0.0) {
		float distance = sqrt((y - center)*(y - center) + (x - center)*(x - center));
		if (distance > max_distance) {
		  max_distance = distance;
		}
	  }
	}
  }
  return max_distance/center;
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

void complexToFloat(vector<Complex>& dataComplex, vector<float>& dataReal, vector<float>& dataImag){
  int size = dataComplex.size();
  dataReal.resize(size);
  dataImag.resize(size);
  for (int i=0; i<size; i++) {
	dataReal[i]=dataComplex[i].real();
	dataImag[i]=dataComplex[i].imag();
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


//Calculate the TCC matrix from the source
void calcTCC(vector<Complex>& tcc, int tccSize, vector<float>& src, int srcSize, float outerSigma,
			 float lambda, float NA, float defocus, int Lx, int Ly, int Nx, int Ny) {
  int sh = ((srcSize - 1) / 2.0);
  int oSgm = ceil(sh * outerSigma);
  float dz = defocus / (NA * NA / (float)lambda), k = 2 * M_PI / lambda;
  int nxMinTmp = -Nx, nxMaxTmp = Nx, nyMinTmp = -Ny, nyMaxTmp = Ny;
  int startmx, srcID, ID1, ID2;
  Complex tmpValComplex;
  float tmpValFloat;
  float sx, sy, fx, fy, rho2;
  float Lx_norm = Lx * NA / lambda, Ly_norm = Ly * NA / lambda;
  Complex srcVal;
  vector<vector<Complex>> pupil(tccSize, vector<Complex>(srcSize * srcSize, Complex(0, 0)));

  //calculate pupil function
  for (int q = -oSgm; q <= oSgm; q++) {
	for (int p = -oSgm; p <= oSgm; p++) {
	  srcID = (q + sh) * srcSize + p + sh;
	  if (src[srcID] != 0 && (p * p + q * q) <= sh * sh) {
		sx = (float)p / sh;
		sy = (float)q / sh;
		nxMinTmp = (-1 - sx) * Lx_norm;
		nxMaxTmp = (1 - sx) * Lx_norm;
		nyMinTmp = (-1 - sy) * Ly_norm;
		nyMaxTmp = (1 - sy) * Ly_norm;

		for (int ny = nyMinTmp; ny <= nyMaxTmp; ny++) {
		  fy = (float)ny / Ly_norm + sy;
		  if (fy * fy <= 1) {
			for (int nx = nxMinTmp; nx <= nxMaxTmp; nx++) {
			  fx = (float)nx / Lx_norm + sx;
			  rho2 = fx * fx + fy * fy;
			  if (rho2 <= 1.0) {
				tmpValFloat = dz * k * (sqrt(1.0 - (rho2)*NA * NA));
				pupil[(ny + Ny) * (Nx * 2 + 1) + nx + Nx][srcID] = Complex(cos(tmpValFloat), sin(tmpValFloat));
			  }
			}
		  }
		}
	  }
	}
  }

  //calculate TCC matrix (Calculate only upper triangular component)
  for (int q = -oSgm; q <= oSgm; q++) {
	for (int p = -oSgm; p <= oSgm; p++) {
	  srcID = (q + sh) * srcSize + p + sh;
	  if (src[srcID] != 0  && (p * p + q * q) <= sh * sh) {
		srcVal = Complex(src[srcID], 0);
		sx = (float)p / sh;
		sy = (float)q / sh;
		nxMinTmp = (-1 - sx) * Lx_norm;
		nxMaxTmp = (1 - sx) * Lx_norm;
		nyMinTmp = (-1 - sy) * Ly_norm;
		nyMaxTmp = (1 - sy) * Ly_norm;

		for (int ny = nyMinTmp; ny <= nyMaxTmp; ny++) {
		  for (int nx = nxMinTmp; nx <= nxMaxTmp; nx++) {
			ID1 = (ny + Ny) * (Nx * 2 + 1) + (nx + Nx);
			if (pupil[ID1][srcID] != Complex(0, 0)) {
			  tmpValComplex= pupil[ID1][srcID] * srcVal;
			  for (int my = ny; my <= nyMaxTmp; my++) {
				startmx = (ny == my) ? nx : nxMinTmp;
				for (int mx = startmx; mx <= nxMaxTmp; mx++) {
				  ID2 = (my + Ny) * (Nx * 2 + 1) + (mx + Nx);
				  tcc[ID1 * tccSize + ID2] += tmpValComplex * conj(pupil[ID2][srcID]);
				}
			  }
			}
		  }
		}
	  }
	}
  }
  //calculate lower triangular component from upper triangular component
  for (int i = 0; i < tccSize; i++) {
	for (int j = i + 1; j < tccSize; j++) {
	  tcc[j * tccSize + i] = conj(tcc[i * tccSize + j]);
	}
  }
}

//Calculate the optical image in the frequency domain using TCC and mask
void calcImage(vector<Complex>& imgf, vector<Complex>& msk, int Lx, int Ly, vector<Complex>& tcc, int tccSize, int Nx, int Ny){
  int tccSizeh = (tccSize - 1) / 2;
  int Lxh = Lx / 2;
  int Lyh = Ly / 2;
  Complex val = 0;
  for (int ny2 = -2 * Ny; ny2 <= 2 * Ny; ny2++) {
	for (int nx2 =  0 * Nx; nx2 <= 2 * Nx; nx2++) {
	  val = 0;
	  for (int ny1 = -1 * Ny; ny1 <= Ny; ny1++) {
		for (int nx1 = -1 * Nx; nx1 <= Nx; nx1++) {
		  if (abs(nx2 + nx1) <= Nx && abs(ny2 + ny1) <= Ny) {
			val += msk[(ny2 + ny1 + Lyh) * Lx + (nx2 + nx1 + Lxh)]
			  * conj(msk[(ny1 + Lyh) * Lx + (nx1 + Lxh)])
			  * tcc[((ny1) * (2 * Nx + 1) + nx1 + tccSizeh) * tccSize + ((ny2 + ny1) * (2 * Nx + 1) + (nx2 + nx1) + tccSizeh)];
		  }
		}
	  }
	  imgf[(ny2 + Lyh) * Lx + (nx2 + Lxh)] = val;
	}
  }
}

extern "C" {
  // LAPACK's eigenvalue problem driver
  void zheevr_(const char& jobz, const char& range, const char& uplo, int* n, Complex* a,
			   int* lda, double* vl, double* vu, int* il, int* iu, double* abstol,
			   int* m, double* w, Complex* z, int* ldz, int* isuppz, Complex* work, int* lwork,
			   double* rwork, int* lrwork, int* iwork, int* liwork, int* info);
}

//Derive the SOCS kernels from TCC
int calcKernels(vector<vector<Complex>>& krns, vector<float>& scales, int nk, vector<Complex>& tcc, int tccSize, int Nx, int Ny) {
  int n = tccSize, lda = tccSize, ldz = tccSize, il = tccSize - (nk - 1), iu = tccSize, m = iu - il + 1, info;
  double abstol, vl, vu;
  Complex* a = nullptr;
  Complex* z = nullptr;
  double* w = nullptr;
  int* isuppz = nullptr;
  try {
	isuppz = new int[2 * m];
	a = new Complex[tccSize * tccSize];
	z = new Complex[tccSize * m];
	w = new double[m];
	for (int i = 0; i < tccSize; i++) {
	  for (int j = 0; j < tccSize; j++) {
        a[i * tccSize + j] = tcc[j * tccSize + i];
	  }
	}

	// Negative abstol means using the default value
	abstol = -1.0;

	// Query and allocate the optimal workspace
	int lwork = -1, lrwork = -1, liwork = -1;
	Complex wkopt;
	double rwkopt;
	int iwkopt;
	zheevr_('V', 'I', 'U', &n, a, &lda, &vl, &vu, &il, &iu,
			&abstol, &m, w, z, &ldz, isuppz, &wkopt, &lwork, &rwkopt, &lrwork, &iwkopt, &liwork, &info);

	lwork = static_cast<int>(wkopt.real());
	lrwork = static_cast<int>(rwkopt);
	liwork = iwkopt;

	Complex* work = static_cast<Complex*>(malloc(lwork * sizeof(Complex)));
	double* rwork = static_cast<double*>(malloc(lrwork * sizeof(double)));
	int* iwork = static_cast<int*>(malloc(liwork * sizeof(int)));

	// Solve eigenproblem
	cout << "Solving eigenproblem..." << endl;
	zheevr_('V', 'I', 'U', &n, a, &lda, &vl, &vu, &il, &iu,
			&abstol, &m, w, z, &ldz, isuppz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);

	cout << "Eigenvalues: " << endl;
	for (int j = m - 1; j >= 0; j--) {
	  cout << w[j] << " ";
	  scales[m - 1 - j] = w[j];
	}
	cout << endl;

	for (int i = 0; i < m; i++) {
	  for (int y = 2 * Ny; y >= 0; y--) {
        for (int x = 2 * Nx; x >= 0; x--) {
		  krns[i][((2 * Ny + 1) - 1 - y) * (2 * Nx + 1) + ((2 * Nx + 1) - 1 - x)] = z[y + x * (2 * Ny + 1) + ((m - 1 - i) * n)];
        }
	  }
	}

	for (int i = 0; i < m; i++) {
	  for (int j = 0; j < n; j++) {
        krns[i][j] = z[(m - 1 - i) * n + (n - 1 - j)];
	  }
	}

	free(work);
	free(rwork);
	free(iwork);
  }
  catch (const bad_alloc&) {
	cerr << "Memory allocation failed." << endl;
	return -1;
  }
  delete[] a;
  delete[] z;
  delete[] w;
  delete[] isuppz;
  return info;
}

int main(int argc, char* argv[]) {
  cout << "************************************" << endl;
  cout << "* Lithography Simulation initiated..." << endl;
  cout << "************************************" << endl;

  auto start = system_clock::now();
  // Input parameters
  if (argc != 23) {
    cerr << "Error: Invalid number of arguments: " << argc << ". Expected 23 arguments." << endl;
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
  // Grid Size of source
  int srcSize = stoi(argv[10]);
  string srcType = argv[11]; 
  Annular annular;
  Dipole dipole;
  CrossQuadrupole crossQuadrupole;
  Point point;
  // Parameters for the annular source
  annular.innerRadius = stof(argv[12]);
  annular.outerRadius = stof(argv[13]);
  // Parameters for the Dipole/Closs-Quadrupole source
  dipole.radius = crossQuadrupole.radius = stof(argv[14]);
  dipole.offset = crossQuadrupole.offset = stof(argv[15]);
  dipole.onXAxis = static_cast<bool>(stoi(argv[16]));
  // Parameters for the point source
  point.x = stof(argv[17]);
  point.y = stof(argv[18]);
  // Input file for the source
  string sourceInFile= argv[19];
  // Numerical aperture
  double NA = stof(argv[20]);
  // Normalized defocus
  float defocus = stof(argv[21]);
  // Number of output kernels
  int nk = stoi(argv[22]);

  // Validate input parameters
  if (Lx <= 0 || Ly <= 0 || maskSizeX <= 0 || maskSizeY <= 0 || Lx < maskSizeX || Ly < maskSizeY || srcSize <= 0) {
    cerr << "Error: Invalid input parameters. Please check the values." << endl;
    return 1;
  }

  // Other parameters
  int Lxy = Lx * Ly;
  float wavelength = 193;
  float lambda = wavelength;
  float dose = 1.0;
  float TARGET_INTENSITY = 0.225;
  float outerSigma;
  float min, max;
  
  if (NA < 0 || NA > 1) {
    cerr << "Error: NA must be in the range [0, 1]." << endl;
    return 1;
  }

  vector<float> src(srcSize * srcSize);
  vector<float> msk(Lxy, 0);
  vector<Complex> mskf(Lxy);
  vector<Complex> imgf(Lxy);
  vector<float> img(Lxy);
  vector<Complex> tcc;
  
  vector<double> realDataLxy(Lxy, 0);
  vector<Complex> complexDataLxy(Ly * (Lx / 2 + 1), 0);
  // Create FFTW plans
  fftw_plan LxyR2C_plan = fftw_plan_dft_r2c_2d(Ly, Lx, realDataLxy.data(), reinterpret_cast<fftw_complex*>(complexDataLxy.data()), FFTW_ESTIMATE);
  fftw_plan LxyC2R_plan = fftw_plan_dft_c2r_2d(Ly, Lx, reinterpret_cast<fftw_complex*>(complexDataLxy.data()), realDataLxy.data(), FFTW_ESTIMATE);

  // Create source
  stringstream ss_srcParam;
  createSource(src, srcSize, srcType, annular, dipole, crossQuadrupole, point, sourceInFile, ss_srcParam);
  
  // Output source to PNG file
  writeFloatArrayToPNGWithGrayscale(srcSize, srcSize, src, getMin(src), getMax(src), "./out/source.png");
  
  // Create mask
  stringstream ss_mskParam;
  createMask(msk, Lx, Ly, maskSizeX, maskSizeY, maskType, lineSpace, maskInFile, ss_mskParam, dose);
 
  // Output mask to PNG file
  writeFloatArrayToPNGWithGrayscale(Lx, Ly, msk, getMin(msk), getMax(msk), "./out/mask.png");
  
  // Calculate Number of elements in TCC matrix
  outerSigma = getOuterSigma(src, srcSize);
  if (outerSigma > 1) {
    cout << "Warning: Outer Sigma is greater than 1. Setting outerSigma = 1." << endl;
	outerSigma = 1;
  }
  int Nx = floor(static_cast<float>(NA) * Lx * (1 + outerSigma) / lambda);
  int Ny = floor(static_cast<float>(NA) * Ly * (1 + outerSigma) / lambda);
  int tccSize = (2 * Nx + 1) * (2 * Ny + 1);
  tcc.resize(tccSize * tccSize, 0);
  if (tccSize < nk) {
    cout << "Warning: Number of kernels (nk) is greater than TCC size (tccSize). Setting nk = tccSize." << endl;
    nk = tccSize;
  }
  
  // Output simulation parameters
  cout << endl << "Simulation Parameters:" << endl;
  cout << "-------------------" << endl;
  cout << ss_mskParam.str();
  cout << ss_srcParam.str();
  cout << "Wavelength: " << wavelength << " nm" << endl;
  cout << "NA: " << NA << endl;
  cout << "Defocus: " << defocus << endl;
  cout << "Dose: " << dose << endl;
  cout << "Nx: " << Nx << ", Ny: " << Ny << endl;
  cout << "TCC Size: " << tccSize << " x " << tccSize << endl;
  cout << "Number of Kernels: " << nk << endl << endl;
  

  // Calculate FFT of mask
  cout << endl << "------------------------------------" << endl;
  cout << "Calculating FFT Mask..." << endl;
  FT_r2c(msk, mskf, Lx, Ly, LxyR2C_plan, realDataLxy, complexDataLxy);

  // Calculate TCC matrix
  cout << endl << "------------------------------------" << endl;
  cout << "Starting calcTCC..." << endl;
  calcTCC(tcc, tccSize, src, srcSize, outerSigma, lambda, NA, defocus, Lx, Ly, Nx, Ny);

  // Output Real/imaginary part of TCC matrix to PNG file
  vector<float> tccReal, tccImag;
  complexToFloat(tcc, tccReal, tccImag);

  max = getMax(tccReal);
  min = getMin(tccReal);
  cout << "Max. TCC Real Value: " << max << endl;
  writeFloatArrayToPNGWithContinuousColor(tccSize, tccSize, tccReal, "./out/tcc_r.png", min, max);

  max = getMax(tccImag);
  min = getMin(tccImag);
  cout << "Max. TCC Imag Value: " << max << endl;
  writeFloatArrayToPNGWithContinuousColor(tccSize, tccSize, tccImag, "./out/tcc_i.png", min, max);

  // Calculate optical image in Fourier domain
  cout << endl << "------------------------------------" << endl;
  cout << "Starting calcImage..." << endl;
  calcImage(imgf, mskf, Lx, Ly, tcc, tccSize, Nx, Ny);
  
  // Inverse Fourier transform of optical image
  FT_c2r(imgf, img, Lx, Ly, LxyC2R_plan, complexDataLxy, realDataLxy);
  
  max = getMax(img);
  min = getMin(img);
  cout << "Max. Intensity: " << max << endl;

  // Output optical image to binary file
  writeFloatArrayToBinary("./out/image.bin", img, Lx*Ly);
  // Output optical image to PNG file
  writeFloatArrayToPNGWithContinuousColor(Lx, Ly, img, "./out/image.png", min, max);
  writeFloatArrayToPNGWith7Colors(Lx, Ly, img, "./out/image_s.png", 0, TARGET_INTENSITY);

  if (nk > 0) {
    cout << endl << "------------------------------------" << endl;
    cout << "Starting calcKernels..." << endl;

    vector<vector<Complex>> krns(nk, vector<Complex>(tccSize));
    vector<float> scales(nk);
	// Calculate kernels
    int info = calcKernels(krns, scales, nk, tcc, tccSize, Nx, Ny);
    if (info != 0) {
	  cerr << "Kernels derivation failed with info = " << info << endl;
	  return 1;
    }

	// Output Real/imaginary part of each kernel to binary/PNG file
    for (int i = 0; i < nk; i++) {
	  vector<float> krnReal, krnImag;
	  complexToFloat(krns[i], krnReal, krnImag);
	  writeFloatArrayToBinary("./out/kernels/krn_" + to_string(i) + "_r.bin", krnReal, tccSize);
	  writeFloatArrayToBinary("./out/kernels/krn_" + to_string(i) + "_i.bin", krnImag, tccSize);
	  writeFloatArrayToPNGWithContinuousColor(2 * Nx + 1, 2 * Ny + 1, krnReal, "./out/kernels/png/krn_" + to_string(i) + "_r.png", getMin(krnReal), getMax(krnReal));
	  writeFloatArrayToPNGWithContinuousColor(2 * Nx + 1, 2 * Ny + 1, krnImag, "./out/kernels/png/krn_" + to_string(i) + "_i.png", getMin(krnImag), getMax(krnImag));
	}
	// Output kernel parameters and eigenvalues
	writeKernelInfo("./out/kernels/kernel_info.txt", 2 * Nx + 1, 2 * Ny + 1, nk, scales);
  }

  // Destroy FFTW plans
  fftw_destroy_plan(LxyR2C_plan);
  fftw_destroy_plan(LxyC2R_plan);

  // Output execution time
  auto end = system_clock::now();
  auto dur = end - start;
  auto microsec = duration_cast<chrono::microseconds>(dur).count();
  cout << endl << "------------------------------------" << endl;
  cout << "Simulation completed." << endl;
  cout << "Total Execution Time: " << microsec * 1.0e-06 << " sec." << endl;
  cout << "************************************" << endl;

  return 0;
}
