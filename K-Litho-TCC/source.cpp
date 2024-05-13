#include <iostream>
#include "file_io.hpp"
#include "source.hpp"

void createAnnular(vector<float>& matrix, int matrixSize, Annular& annular) {
  int innerRadiusSize = annular.innerRadius * ((matrixSize - 1) / 2);
  int outerRadiusSize = annular.outerRadius * ((matrixSize - 1) / 2);
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
void createCrossQuadrupole(vector<float>& matrix, int matrixSize, CrossQuadrupole& crossQuadrupole) {
  int offsetSize = crossQuadrupole.offset * ((matrixSize - 1) / 2);
  int radiusSize = crossQuadrupole.radius * ((matrixSize - 1) / 2);
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
void createDipole(vector<float>& matrix, int matrixSize, Dipole& dipole) {
  int offsetSize = dipole.offset * ((matrixSize - 1) / 2);
  int radiusSize = dipole.radius * ((matrixSize - 1) / 2);
  int centerX = matrixSize / 2;
  int centerY = matrixSize / 2;

  for (int y = 0; y < matrixSize; ++y) {
	for (int x = 0; x < matrixSize; ++x) {
	  // 2つの円の中心座標を計算
	  int centers[2][2] = {
		{centerX + (dipole.onXAxis ? offsetSize : 0), centerY + (dipole.onXAxis ? 0 : offsetSize)},  // 1つ目の円の中心
		{centerX - (dipole.onXAxis ? offsetSize : 0), centerY - (dipole.onXAxis ? 0 : offsetSize)}   // 2つ目の円の中心
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

void createSource(vector<float>& src, int srcSize, string srcType,
				  Annular& annular, Dipole& dipole, CrossQuadrupole& crossQuadrupole, Point& point,
				  string sourceInFile, stringstream& ss_srcParam){
  ss_srcParam << "Source Size: " << srcSize << " x " << srcSize << "\nSource Type: " << srcType << "\n";
  // Create source
  if (srcType == "Annular") {	
    ss_srcParam << "Inner Radius: " << annular.innerRadius << "\nOuter Radius: " << annular.outerRadius << "\n";
    createAnnular(src, srcSize, annular);
  }
  else if (srcType == "Dipole") {
    ss_srcParam << "Radius: " << dipole.radius << "\nOffset: " << dipole.offset << "\nOnXAxis: " << dipole.onXAxis << "\n";
    createDipole(src, srcSize, dipole);
  }
  else if (srcType == "CrossQuadrupole") {
    ss_srcParam << "Radius: " << crossQuadrupole.radius << "\nOffset: " << crossQuadrupole.offset << "\n";
    createCrossQuadrupole(src, srcSize, crossQuadrupole);
  }
  else if (srcType == "Point") {
    ss_srcParam << "Point Position: (" << point.x << ", " << point.y << ")\n";
    int ptXId = round(point.x * srcSize / 2.0 + srcSize / 2.0);
    int ptYId = round(-1 * point.y * srcSize / 2.0 + srcSize / 2.0);
    if (point.x >= 0) ptXId -= 1;
    if (-1 * point.y >= 0) ptYId -= 1;
    src[ptYId * srcSize + ptXId] = 1.0;
  }
  else if (srcType == "Import") {
    cout << "Read Source: " << sourceInFile <<endl;
    string extension = sourceInFile.substr(sourceInFile.find_last_of(".") + 1);
    if(extension == "bin"){
	  readFloatArrayFromBinary(sourceInFile, src, srcSize * srcSize);
    }
    else {
	  readFloatArrayFromTxt(sourceInFile, src, srcSize * srcSize);
    }
  }
  // Normalizing source elements to sum to 1.
  normalizeMatrix(src, srcSize * srcSize);
}
