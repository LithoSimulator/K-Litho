#ifndef SOURCE_HPP
#define SOURCE_HPP

#include <string>
#include <vector>
#include <sstream>
using namespace std;

struct Annular {
  float innerRadius;
  float outerRadius;
};
struct Dipole {
  float radius;
  float offset;
  bool onXAxis;
};
struct CrossQuadrupole {
  float radius;
  float offset;
};
struct Point {
  float x;
  float y;
};

void createAnnular(vector<float>& matrix, int matrixSize, Annular& annular) ;
void createDipole(vector<float>& matrix, int matrixSize, Dipole& dipole) ;
void createCrossQuadrupole(vector<float>& matrix, int matrixSize, CrossQuadrupole& crossQuadrupole) ;
void normalizeMatrix(vector<float>& matrix, int matrixSize);
void createSource(vector<float>& src, int srcSize, string srcType,
				  Annular& annular, Dipole& dipole, CrossQuadrupole& crossQuadrupole, Point& point,
				  string sourceInFile, stringstream& ss_srcParam);


#endif
