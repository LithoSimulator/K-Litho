#include <iostream>
#include <complex> 
#include <string>
#include <fstream>
#include <chrono>
#include "fftw3.h"
#include "file_io.hpp"
#include <vector>
#include <cmath>

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
void myShift(T1 *in, T2 *out, int sizeX, int sizeY, bool shiftTypeX, bool shiftTypeY) {
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

//Create a parametric source
void createAnnular(double innerRadius, double outerRadius, int matrixSize, vector<float>& matrix) {
  int innerRadiusSize = innerRadius * ((matrixSize - 1) / 2);
  int outerRadiusSize = outerRadius * ((matrixSize - 1) / 2);
  int centerX = ((matrixSize - 1) / 2);
  int centerY = ((matrixSize - 1) / 2);

  for (int y = 0; y < matrixSize; ++y) {
	for (int x = 0; x < matrixSize; ++x) {
	  // Calculate the distance from the current point to the center point
	  double distance = std::sqrt(std::pow(x - centerX, 2) + std::pow(y - centerY, 2));

	  // Check if the distance is within the range of inner and outer radii
	  if (distance >= innerRadiusSize && distance <= outerRadiusSize) {
		matrix[y * matrixSize + x] = 1; // Inside
	  } else {
		matrix[y * matrixSize + x] = 0; // Outside
	  }
	}
  }
}
void createCrossQuadrupole(double radius, double offset, int matrixSize, vector<float>& matrix) {
  int offsetSize = offset * ((matrixSize - 1) / 2);
  int radiusSize = radius * ((matrixSize - 1) / 2);
  int centerX = ((matrixSize - 1) / 2);
  int centerY = ((matrixSize - 1) / 2);

  for (int y = 0; y < matrixSize; ++y) {
	for (int x = 0; x < matrixSize; ++x) {
	  // Calculate the center coordinates of the four circles
	  int centers[4][2] = {
		{centerX, centerY + offsetSize},       // Top
		{centerX, centerY - offsetSize},       // Bottom
		{centerX + offsetSize, centerY},       // Right
		{centerX - offsetSize, centerY}        // Left
	  };

	  // Calculate the distance from each circle's center and check if it's within the radius
	  bool insideAnyCircle = false;
	  for (int i = 0; i < 4; ++i) {
		double dx = x - centers[i][0];
		double dy = y - centers[i][1];
		double distance = std::sqrt(dx * dx + dy * dy);
		if (distance <= radiusSize) {
		  insideAnyCircle = true;
		  break;  // If inside at least one circle, no need to check the others
		}
	  }

	  // Store 1.0 for inside, 0.0 for outside
	  matrix[y * matrixSize + x] = insideAnyCircle ? 1.0 : 0.0;
	}
  }
}
void createDipole(double radius, double offset, int matrixSize, vector<float>& matrix, bool onXAxis) {
  int offsetSize = offset * ((matrixSize - 1) / 2);
  int radiusSize = radius * ((matrixSize - 1) / 2);
  int centerX = matrixSize / 2;
  int centerY = matrixSize / 2;

  for (int y = 0; y < matrixSize; ++y) {
	for (int x = 0; x < matrixSize; ++x) {
	  // 2つの円の中心座標を計算
	  int centers[2][2] = {
		{centerX + (onXAxis ? offsetSize : 0), centerY + (onXAxis ? 0 : offsetSize)},  // 1つ目の円の中心
		{centerX - (onXAxis ? offsetSize : 0), centerY - (onXAxis ? 0 : offsetSize)}   // 2つ目の円の中心
	  };

	  // 各円の中心からの距離を計算し、半径内にあるかどうか判定
	  bool insideAnyCircle = false;
	  for (int i = 0; i < 2; ++i) {
		double dx = x - centers[i][0];
		double dy = y - centers[i][1];
		double distance = std::sqrt(dx * dx + dy * dy);
		if (distance <= radiusSize) {
		  insideAnyCircle = true;
		  break;  // 少なくとも1つの円の内側にある場合は他の円はチェック不要
		}
	  }

	  // 内側は1.0、外側は0.0で格納
	  matrix[y * matrixSize + x] = insideAnyCircle ? 1.0 : 0.0;
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

void normalizeMatrix(vector<float>& matrix, int matrixSize) {
  float sum = 0.0;

  // First, calculate the sum of all elements
  for (int i = 0; i < matrixSize; ++i) {
	sum += matrix[i];
  }
  
  // If the sum is not zero, normalize each element by dividing it by the sum
  if (sum != 0.0) {
	for (int i = 0; i < matrixSize; ++i) {
	  matrix[i] /= sum;
	}
  }
}

//Calculate the TCC matrix from the source
void calcTCC(int srcSize, vector<float>& src, int tccSize, vector<Complex>& tcc, float outerSigma, float NA, float defocus,  bool srcSimX, bool srcSimY, int Nx, int Ny, int Lx, int Ly, float lambda){
  int sh = ((srcSize-1)/2.0);
  int oSgm= ceil(sh*outerSigma);
  Complex val;
  float kx1_2=0, ky1_2=0, kx2_2=0, ky2_2=0;
  float sx,sy;

  float dz=defocus/(NA*NA/(float)lambda), k= 2*M_PI/lambda;

  vector<vector< float > > sxkx_vals(srcSize, vector<float>(2*Nx+1)), syky_vals(srcSize, vector<float>(2*Ny+1));

  int startX, starY;
  vector<float> tmpSrc(srcSize*srcSize,0);
  vector<Complex> tmpTcc(tccSize*tccSize,0);

  startX=-oSgm; starY=-oSgm;
  if(srcSimX){startX=0;}
  if(srcSimY){starY=0;}
  
  for(int q=starY; q<=oSgm; q++){
	for(int p=startX; p<=oSgm; p++){
	  tmpSrc[(q+sh)*srcSize+sh+p]=src[(q+sh)*srcSize+sh+p];
	}
  }
  if(srcSimX){
	for(int q=starY; q<=oSgm; q++){
	  tmpSrc[(q+sh)*srcSize+sh]/=2.0;
	}
  }
  if(srcSimY){
	for(int p=startX; p<=oSgm; p++){
	  tmpSrc[(sh)*srcSize+sh+p]/=2.0;
	}
  }
  for(int p=startX; p<=oSgm; p++){
	sx=(float)p/sh;
	for(int i=-Nx; i<=Nx; i++){
	  sxkx_vals[p+sh][i+Nx]= pow( i*lambda/(NA*Lx)+sx,2);
	}
  }
  for(int q=starY; q<=oSgm; q++){
	sy=(float)q/sh;
	for(int i=-Ny; i<=Ny; i++){
	  syky_vals[q+sh][i+Ny]= pow(i*lambda/(NA*Ly)+sy,2);
	}
  }
  for(int q=starY; q<=oSgm; q++){
	sy=(float)q/sh;
	for(int p=startX; p<=oSgm; p++){
		
	  int sc=(q+sh)*srcSize+(p+sh);
	  float sr=tmpSrc[sc];
	  if(sr==0)continue;
	  sx= (float)p/sh;

	  for(int ny=-Ny; ny<=Ny; ny++){
		ky2_2=syky_vals[q+sh][ny+Ny];
		for(int nx=-Nx; nx<=Nx; nx++){
		  kx2_2=sxkx_vals[p+sh][nx+Nx];
		  if((kx2_2 + ky2_2)>1.0) continue;
		  int idy=((ny+Ny)*(Nx*2+1) + nx+Nx)*tccSize;
		  float SqrtK2 = sqrt(1-(kx2_2 + ky2_2)*NA*NA);
		  for(int my=-Ny; my<=Ny; my++){
			ky1_2=syky_vals[q+sh][my+Ny];
			for(int mx=-Nx; mx<=Nx ; mx++){
			  if(((ny+Ny)*(Nx*2+1) + nx+Nx) > ((my+Ny)*(Nx*2+1) + mx+Nx)){continue;}
			  kx1_2=sxkx_vals[p+sh][mx+Nx];
			  if((kx1_2 + ky1_2)<=1.0){
				//cout << "\r" << ++i << "/" << tccSize*tccSize << std::flush;
				float th=dz*k*(sqrt(1-(kx1_2 + ky1_2)*NA*NA)-SqrtK2);
				val=Complex(cos(th),-sin(th)); val*=sr;
				tmpTcc[idy+((my+Ny)*(Nx*2+1) + mx+Nx)]+= val;  
			  }
			}
		  }	    
		}
      
	  }
	}
  }

  for(int ny=-Ny, ny2=Ny; ny<=Ny; ny++, ny2--){
	for(int nx=-Nx, nx2=Nx; nx<=Nx; nx++, nx2--){
	  for(int my=-Ny, my2=Ny; my<=Ny; my++, my2--){
		for(int mx=-Nx, mx2=Nx; mx<=Nx; mx++, mx2--){
		  if(((ny+Ny)*(Nx*2+1) + nx+Nx) > ((my+Ny)*(Nx*2+1) + mx+Nx)){continue;}
			
		  int idx1=((ny+Ny)*(Nx*2+1) + nx+Nx)*tccSize+((my+Ny)*(Nx*2+1) + mx+Nx); int idx1T=((my+Ny)*(Nx*2+1) + mx+Nx)*tccSize +((ny+Ny)*(Nx*2+1) + nx+Nx);
		  Complex val=tmpTcc[idx1];
		  Complex valc=conj(tmpTcc[idx1]);

		  tcc[idx1] += val;
		  if(((ny+Ny)*(Nx*2+1) + nx+Nx) != ((my+Ny)*(Nx*2+1) + mx+Nx)){
			tcc[idx1T] += valc;
		  }
		  if(srcSimX){
			int idx2=((ny+Ny)*(Nx*2+1) + nx2+Nx)*tccSize+((my+Ny)*(Nx*2+1) + mx2+Nx); int idx2T=((my+Ny)*(Nx*2+1) + mx2+Nx)*tccSize + ((ny+Ny)*(Nx*2+1) + nx2+Nx);
			tcc[idx2] += val;
			if(((ny+Ny)*(Nx*2+1) + nx+Nx) != ((my+Ny)*(Nx*2+1) + mx+Nx)){
			  tcc[idx2T] += valc;
			}
		  }
		  if(srcSimY){
			int idx3=((ny2+Ny)*(Nx*2+1) + nx+Nx)*tccSize+((my2+Ny)*(Nx*2+1) + mx+Nx); int idx3T=((my2+Ny)*(Nx*2+1) + mx+Nx)*tccSize + ((ny2+Ny)*(Nx*2+1) + nx+Nx);
			tcc[idx3] += val;
			if(((ny+Ny)*(Nx*2+1) + nx+Nx) != ((my+Ny)*(Nx*2+1) + mx+Nx)){
			  tcc[idx3T] += valc;
			}
		  }
		  if(srcSimX && srcSimY){
			int idx4=((ny2+Ny)*(Nx*2+1) + nx2+Nx)*tccSize+((my2+Ny)*(Nx*2+1) + mx2+Nx); int idx4T=((my2+Ny)*(Nx*2+1) + mx2+Nx)*tccSize+((ny2+Ny)*(Nx*2+1) + nx2+Nx);
			tcc[idx4] += val;
			if(((ny+Ny)*(Nx*2+1) + nx+Nx) != ((my+Ny)*(Nx*2+1) + mx+Nx)){
			  tcc[idx4T] += valc;
			}
		  }
			
		}
		
	  }
	}
  }
}

//Calculate the optical image in the frequency domain using TCC and mask
void calcImage(vector<Complex>& imgf, vector<Complex>& msk, int tccSize, vector<Complex>& tcc, int Nx, int Ny, int Lx, int Ly){
  int ht = (tccSize-1)/2;
  int Lxh = Lx/2;
  int Lyh = Ly/2;
  Complex val=0;
  for(int ny2=-2*Ny; ny2<=2*Ny; ny2++){
	for(int nx2=0; nx2<=2*Nx; nx2++){	  
	  val=0;
	  for(int ny1=-1*Ny; ny1<=Ny; ny1++){
		for(int nx1=-1*Nx; nx1<=Nx; nx1++){
		  if(abs(nx2+nx1)>Nx || abs(ny2+ny1)>Ny)continue;
		  else{
			val+= msk[(ny2+ny1+Lyh)*Lx+(nx2+nx1+Lxh)]*conj(msk[(ny1+Lyh)*Lx+(nx1+Lxh)])*tcc[((ny1)*(2*Nx+1)+nx1+ht)*tccSize+((ny2+ny1)*(2*Nx+1)+(nx2+nx1)+ht)];
		  }
		}
	  }
	  imgf[(ny2+Lyh)*Lx+(nx2+Lxh)]=val;
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
int calcKernels(vector<Complex>& tcc, int tccSize, int nk, vector<vector<Complex>>& krns, vector<float>& scales, string fname, int Nx, int Ny) {
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
		  krns[i][((2 * Ny + 1) - 1 - y) * (2 * Nx + 1) + ((2 * Nx + 1) - 1 - x)] = z[y + (x) * (2 * Ny + 1) + ((m - 1 - i) * n)];
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

//Fast Fourier transform from real data to complex data
void FT_r2c(vector<float>& in, vector<Complex>& out, int sizeX, int sizeY, fftw_plan &fwd, vector<double>& datar, vector<Complex>& datac){
  myShift(in.data(), datar.data(),sizeX,sizeY,false,false);
 
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
 
  myShift(datar.data(),out.data(),sizeX,sizeY,false,false);
}

bool isSymmetricY(vector<float>& matrix, int width, int height) {
  for (int y = 0; y < height / 2; ++y) {
	for (int x = 0; x < width; ++x) {
	  if (matrix[y * width + x] != matrix[(height - 1 - y) * width + x]) {
		return false;
	  }
	}
  }
  return true;
}

bool isSymmetricX(vector<float>& matrix, int width, int height) {
  for (int y = 0; y < height; ++y) {
	for (int x = 0; x < width / 2; ++x) {
	  if (matrix[y * width + x] != matrix[y * width + (width - 1 - x)]) {
		return false;
	  }
	}
  }
  return true;
}

int main(int argc, char* argv[]) {
  cout << "*****************************************************************" << endl;
  cout << "* Lithography Simulation initiated..." << endl;
  cout << "*****************************************************************" << endl;

  auto start = system_clock::now();
  // Input parameters
  if (argc != 19) {
    cerr << "Error: Invalid number of arguments. Expected 18 arguments." << endl;
    return 1;
  }
  
  const int Lx = stoi(argv[1]);
  const int Ly = stoi(argv[2]);
  const int maskSizeX = stoi(argv[3]);
  const int maskSizeY = stoi(argv[4]);
  const string maskInFile = argv[5];
  const int srcSize = stoi(argv[6]);
  const string srctype = argv[7];
  const float innerRad = stof(argv[8]);
  const float outerRad = stof(argv[9]);
  const float Rad = stof(argv[10]);
  const float Offset = stof(argv[11]);
  const bool onXAxis = static_cast<bool>(stoi(argv[12]));
  const float ptX = stof(argv[13]);
  const float ptY = stof(argv[14]);
  const string sourceInFile= argv[15];
  const double NA = stof(argv[16]);
  const float defocus = stof(argv[17]);
  const int nk = stoi(argv[18]);

  // Validate input parameters
  if (Lx <= 0 || Ly <= 0 || maskSizeX <= 0 || maskSizeY <= 0 || Lx < maskSizeX || Ly < maskSizeY || srcSize <= 0 || NA <= 0) {
    cerr << "Error: Invalid input parameters. Please check the values." << endl;
    return 1;
  }

  // Other parameters
  const int Lxy = Lx * Ly;
  const float wavelength = 193;
  const float lambda = wavelength;
  const float dose = 1.0;
  const float TARGET_INTENSITY = 0.225;
  
  // Source parameters
  float outerSigma = 0, innerSigma = 0;
  int ptXId = 0, ptYId = 0;

  if (srctype == "Annular") {
    innerSigma = innerRad;
    outerSigma = outerRad;
  } else if (srctype == "CrossQuadrupole") {
    outerSigma = Offset + Rad;
    innerSigma = Offset - Rad;
  } else if (srctype == "Dipole") {
    outerSigma = Offset + Rad;
    innerSigma = Offset - Rad;
  } else if (srctype == "Point") {
    outerSigma = sqrt(ptX * ptX + ptY * ptY);
    float sh = srcSize / 2.0;
    ptXId = round(ptX * sh + sh);
    ptYId = round(-ptY * sh + sh);
    if (ptX >= 0) ptXId -= 1;
    if (-ptY >= 0) ptYId -= 1;
    innerSigma = 0;
  } else if (srctype == "Import") {
	outerSigma = 1;
	innerSigma = 0;
  }else {
    cerr << "Error: Invalid source type. Supported types: Annular, CrossQuadrupole, Dipole, Point, Import." << endl;
    return 1;
  }

  float min, max;
  const int Nx = floor(static_cast<float>(NA) * Lx * (1 + outerSigma) / lambda);
  const int Ny = floor(static_cast<float>(NA) * Ly * (1 + outerSigma) / lambda);
  const int tccSize = (2 * Nx + 1) * (2 * Ny + 1);

  // Output simulation parameters
  cout << endl << "Simulation Parameters:" << endl;
  cout << "-------------------" << endl;
  cout << "Source Type: " << srctype << endl;
  cout << "Source Size: " << srcSize << "x" << srcSize << endl;

  if (srctype == "Annular") {
    cout << "Inner Radius: " << innerRad << endl;
    cout << "Outer Radius: " << outerRad << endl;
  } else if (srctype == "CrossQuadrupole") {
    cout << "Radius: " << Rad << endl;
    cout << "Offset: " << Offset << endl;
  } else if (srctype == "Dipole") {
    cout << "Radius: " << Rad << endl;
    cout << "Offset: " << Offset << endl;
    cout << "onXAxis: " << boolalpha << onXAxis << endl;
  } else if (srctype == "Point") {
    cout << "Point Coordinate: (" << ptX << ", " << ptY << ")" << endl;
  }
  
  cout << "Wavelength: " << wavelength << " nm" << endl;
  cout << "NA: " << NA << endl;
  cout << "Defocus: " << defocus << endl;
  cout << "Dose: " << dose << endl;
  cout << "Mask Period: " << Lx << " x " << Ly << endl;
  cout << "Mask Size: " << maskSizeX << " x " << maskSizeY << endl;
  cout << "Nx: " << Nx << ", Ny: " << Ny << endl;
  cout << "TCC Size: " << tccSize << " x " << tccSize << endl;
  cout << "Number of Kernels: " << nk << endl << endl;

  // Validate parameters
  if (tccSize < nk) {
    cerr << "Error: Number of kernels cannot be greater than TCC size." << endl;
    return 1;
  }

  if (NA < 0 || NA > 1) {
    cerr << "Error: NA must be in the range [0, 1]." << endl;
    return 1;
  }
  
  // Allocate memory for data structures
  vector<double> realDataLxy(Lxy, 0);
  vector<Complex> complexDataLxy(Ly * (Lx / 2 + 1), 0);
  vector<Complex> mskf(Lxy);
  vector<float> inMsk(maskSizeX * maskSizeY);
  vector<float> msk(Lxy, 0);
  vector<float> src(srcSize * srcSize);
  vector<Complex> imgf(Lxy);
  vector<float> img(Lxy);
  vector<Complex> tcc(tccSize * tccSize, 0);
  
  // Create FFTW plans
  fftw_plan LxyR2C_plan = fftw_plan_dft_r2c_2d(Ly, Lx, realDataLxy.data(), reinterpret_cast<fftw_complex*>(complexDataLxy.data()), FFTW_ESTIMATE);
  fftw_plan LxyC2R_plan = fftw_plan_dft_c2r_2d(Ly, Lx, reinterpret_cast<fftw_complex*>(complexDataLxy.data()), realDataLxy.data(), FFTW_ESTIMATE);

  // Create source
  if (srctype == "Annular") {
    createAnnular(innerRad, outerRad, srcSize, src);
  } else if (srctype == "CrossQuadrupole") {
    createCrossQuadrupole(Rad, Offset, srcSize, src);
  } else if (srctype == "Dipole") {
    createDipole(Rad, Offset, srcSize, src, onXAxis);
  } else if (srctype == "Point") {
    src[ptYId * srcSize + ptXId] = 1.0;
  } else if (srctype == "Import") {
	cout << "Read Source: " << sourceInFile <<endl;
	string extension = sourceInFile.substr(sourceInFile.find_last_of(".") + 1);
	if(extension == "bin"){
	  if (!readFloatArrayFromBinary(sourceInFile, src, srcSize * srcSize)) {
		cerr << "Error: Failed to read source file: " << sourceInFile << endl;
		return 1;
	  } 
	}else if(extension == "txt"){
	  if (!readFloatArrayFromText(sourceInFile, src, srcSize * srcSize)) {
		cerr << "Error: Failed to read source file: " << sourceInFile << endl;
		return 1;
	  } 
	}else {
	  cerr << "Error: Unsupported file extension " << extension << endl;
	  return 1;
	}
    outerSigma = getOuterSigma(src, srcSize);
  }
  // cout << "Outer Sigma: " << outerSigma <<endl;

  // Normalizing source elements to sum to 1.
  normalizeMatrix(src, srcSize * srcSize);

  max = getMax(src, srcSize * srcSize);
  min = getMin(src, srcSize * srcSize);

  writeFloatArrayToPNGWithGrayscale(srcSize, srcSize, src, min, max, "./out/source.png");
  
  if (min < 0) {
    cerr << "Error: Source elements must be greater than or equal to zero." << std::endl;
    return 1;
  }
  if (outerSigma < 0 || outerSigma > 1) {
    cerr << "Error: Outer sigma must be in the range [0, 1]." << endl;
    return 1;
  }
  if (innerSigma < 0 || innerSigma > 1) {
    cerr << "Error: Inner sigma must be in the range [0, 1]." << endl;
    return 1;
  }
  if (innerSigma > outerSigma) {
    cerr << "Error: outer Sigma must be greater than inner Sigma." << endl;
    return 1;
  }
  const bool srcSymX = isSymmetricX(src, srcSize, srcSize);
  const bool srcSymY = isSymmetricY(src, srcSize, srcSize);
  // cout << "Source SymmetricX: " << boolalpha << srcSymX <<endl;
  // cout << "Source SymmetricY: " << boolalpha << srcSymY <<endl;

  // Read mask data
  cout << "Read Mask: " << maskInFile <<endl;
  string extension = maskInFile.substr(maskInFile.find_last_of(".") + 1);
  if(extension == "bin"){
	if (!readFloatArrayFromBinary(maskInFile, inMsk, maskSizeX * maskSizeY)) {
	  cerr << "Error: Failed to read mask file: " << maskInFile << endl;
	  return 1;
	} 
  }else if(extension == "txt"){
	if (!readFloatArrayFromText(maskInFile, inMsk, maskSizeX * maskSizeY)) {
	  cerr << "Error: Failed to read mask file: " << maskInFile << endl;
	  return 1;
	} 
  }else {
	cerr << "Error: Unsupported file extension " << extension << endl;
	return 1;
  }
  const int difYh = (Ly - maskSizeY) / 2;
  const int difXh = (Lx - maskSizeX) / 2;
  for (int y = 0; y < maskSizeY; y++) {
    for (int x = 0; x < maskSizeX; x++) {
      if (inMsk[y * maskSizeX + x] != 0) {
        msk[(y + difYh) * Lx + x + difXh] = inMsk[y * maskSizeX + x] * dose;
      }
    }
  }
  max = getMax(inMsk, maskSizeX * maskSizeY);
  min = getMin(inMsk, maskSizeX * maskSizeY);
  if (min < 0 || max > 1) {
    std::cerr << "Error: Mask values must be in the range [0, 1]." << std::endl;
    return 1;
  }
  writeFloatArrayToPNGWithGrayscale(maskSizeX, maskSizeY, inMsk, min, max, "./out/mask.png");
  
  // Calculate FFT of mask
  cout << endl << "-----------------------------------------------------------------" << endl;
  cout << "Calculating FFT Mask..." << endl;
  FT_r2c(msk, mskf, Lx, Ly, LxyR2C_plan, realDataLxy, complexDataLxy);

  // Calculate TCC matrix
  cout << endl << "-----------------------------------------------------------------" << endl;
  cout << "Starting calcTCC..." << endl;
  for(int i=0; i<1; i++){
	calcTCC(srcSize, src, tccSize, tcc, outerSigma, NA, defocus, srcSymX, srcSymY, Nx, Ny, Lx, Ly, lambda);
  }
  vector<float> tccReal(tccSize * tccSize), tccImag(tccSize * tccSize);
  for (int i = 0; i < tccSize * tccSize; i++) {
    tccReal[i] = tcc[i].real();
    tccImag[i] = tcc[i].imag();
  }
  max = getMax(tccReal, tccSize * tccSize);
  min = getMin(tccReal, tccSize * tccSize);
  cout << "Max. TCC Real Value: " << max << endl;
  cout << "Min. TCC Real Value: " << min << endl;
  writeFloatArrayToPNGWithContinuousColor(tccSize, tccSize, tccReal, "./out/tcc_r.png", min, max);
  
  max = getMax(tccImag, tccSize * tccSize);
  min = getMin(tccImag, tccSize * tccSize);
  cout << "Man. TCC Imag Value: " << max << endl;
  cout << "Min. TCC Imag Value: " << min << endl;
  writeFloatArrayToPNGWithContinuousColor(tccSize, tccSize, tccImag, "./out/tcc_i.png", min, max);
  
  // Calculate optical image in Fourier domain
  cout << endl << "-----------------------------------------------------------------" << endl;
  cout << "Starting calcImage..." << endl;
  calcImage(imgf, mskf, tccSize, tcc, Nx, Ny, Lx, Ly);

  // Inverse Fourier transform of optical image
  FT_c2r(imgf, img, Lx, Ly, LxyC2R_plan, complexDataLxy, realDataLxy);

  max = getMax(img, Lx * Ly);
  min = getMin(img, Lx * Ly);
  cout << "Max. Intensity: " << max << endl;
  cout << "Min. Intensity: " << min << endl;
  
  if (!writeFloatArrayToBinary("./out/image.bin", img, Lx*Ly)) {
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
  
  if (nk > 0) {
	cout << endl << "-----------------------------------------------------------------" << endl;
	cout << "Starting calcKernels..." << endl;

	vector<vector<Complex>> krns(nk, vector<Complex>(tccSize));

	vector<float> scales(nk);

	int info = calcKernels(tcc, tccSize, nk, krns, scales, "./out/kernels/", Nx, Ny);
	
    if (info != 0) {
	  cerr << "Calculation failed with info = " << info << endl;
	  return 1;
    }
	
	for (int i = 0; i < nk; i++) {
	  vector<float> dataReal(tccSize), dataImag(tccSize);
	  for (int j = 0; j < tccSize; j++) {
		dataReal[j] = krns[i][j].real();
		dataImag[j] = krns[i][j].imag();
	  }

	  writeFloatArrayToPNGWithContinuousColor(2 * Nx + 1, 2 * Ny + 1, dataReal, "./out/kernels/png/krn_" + to_string(i) + "_r.png", getMin(dataReal, tccSize), getMax(dataReal, tccSize));
	  writeFloatArrayToPNGWithContinuousColor(2 * Nx + 1, 2 * Ny + 1, dataImag, "./out/kernels/png/krn_" + to_string(i) + "_i.png", getMin(dataImag, tccSize), getMax(dataImag, tccSize));
	  // writeFloatArrayToTxt("./out/kernels/krn_" + to_string(i) + "_r.txt", dataReal, (2 * Nx + 1), (2 * Ny + 1)) ;
	  // writeFloatArrayToTxt("./out/kernels/krn_" + to_string(i) + "_i.txt", dataImag, (2 * Nx + 1), (2 * Ny + 1)) ;
	  if (!writeFloatArrayToBinary("./out/kernels/krn_" + to_string(i) + "_r.bin", dataReal, tccSize) || !writeFloatArrayToBinary("./out/kernels/krn_" + to_string(i) + "_i.bin", dataImag, tccSize)) {
		cerr << "Error: Failed to write kernel " << i << " to file." << endl; 
		return 1;
	  }
	}
	
	writeKernelInfo("./out/kernels/kernel_info.txt", 2 * Nx + 1, 2 * Ny + 1, nk, scales);
	
	
  }
  
  // Destroy FFTW plans
  fftw_destroy_plan(LxyR2C_plan);
  fftw_destroy_plan(LxyC2R_plan);

  // Output execution time
  auto end = system_clock::now();
  auto dur = end - start;
  auto microsec = duration_cast<chrono::microseconds>(dur).count();
  cout << endl << "-----------------------------------------------------------------" << endl;
  cout << "Simulation completed." << endl;
  cout << "Total Execution Time: " << microsec * 1.0e-06 << " sec." << endl;
  cout << "*****************************************************************" << endl;

  return 0;
}
