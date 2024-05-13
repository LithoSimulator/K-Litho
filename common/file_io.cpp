#include <iostream>
#include <fstream>
#include "file_io.hpp"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

void writeFloatArrayToPNGWithGrayscale(int sizeX, int sizeY, vector<float>& data, float min, float max, string filename){
  // cout << "write " << fixed<< filename << ": black= "<<min<<", white= "<<max<<endl;
  vector<unsigned char> grayscaleData(sizeX * sizeY);
  
  // Convert float values to grayscale (0-255)
  int size = sizeX * sizeY;
  for (int i = 0; i < size; ++i) {
	float value = std::max(0.0f, std::min(1.0f, (data[i] - min)/(max - min)));
	grayscaleData[i] = static_cast<unsigned char>(value * 255.0f);
  }
    
  // Write grayscale data to PNG file
  if (!stbi_write_png(filename.c_str(), sizeX, sizeY, 1, grayscaleData.data(), sizeX)) {
	cerr << "Error: Failed to write PNG: " << filename << endl;
	exit(1); // Error
  }
}

void writeFloatArrayToPNGWith7Colors(int sizeX, int sizeY, vector<float>& data, string filename, float min, float max) {
  // cout << "write " << filename << ": blue= "<<min<<", red= "<<max<<endl;
  int size = sizeX * sizeY;
  vector<unsigned char> imgData(size * 3); // 3 for RGB
  if (min == max) {
	for(int i = 0; i < size; ++i){
	  imgData[i * 3]=0; // R
	  imgData[i * 3 + 1]=255; // G
	  imgData[i * 3 + 2]=0; // B
	}
  }
  else {
	for (int i = 0; i < size; ++i) {
	  double value = std::max(0.0f, std::min(1.0f, (data[i] - min)/(max - min)));
	  unsigned char R = 0, G = 0, B = 0; // Initialize all colors to 0

	  if (value < 1.0 / 6.0) {
		// blue
		B = 255;
		G = 0;
	  } else if (value < 2.0 / 6.0) {
		// aqua
		G = 127;
		B = 255;
	  } else if (value < 3.0 / 6.0) {
		// bright green
		G = 255;
		B = 127;
	  } else if (value < 4.0 / 6.0) {
		// pure green
		R = 0;
		G = 255;
	  } else if (value < 5.0 / 6.0) {
		// yellow
		R = 127;
		G = 255;
	  } else if (value < 6.0 / 6.0) {
		// orange
		R = 255;
		G = 127;
	  }else{
		// red
		R=255;
	  }

	  imgData[i * 3] = R;
	  imgData[i * 3 + 1] = G;
	  imgData[i * 3 + 2] = B;
	}
  }
  // Write to PNG
  if (!stbi_write_png(filename.c_str(), sizeX, sizeY, 3, imgData.data(), sizeX * 3)) {
	cerr << "Error: Failed to write PNG: " << filename << endl;
	exit(1);
  }
}

void writeFloatArrayToPNGWithContinuousColor(int sizeX, int sizeY, vector<float>& data, string filename, float min, float max) {
  // cout << "write " << fixed<< filename << ": blue= "<<min<<", red= "<<max<<endl;
  int size = sizeX * sizeY;
  vector<unsigned char> imgData(size * 3); // 3 for RGB
  if(min == max){
	for(int i = 0; i < size; ++i){
	  imgData[i * 3]=0; // R
	  imgData[i * 3 + 1]=255; // G
	  imgData[i * 3 + 2]=0; // B
	}
  }
  else{
	for (int i = 0; i < size; ++i) {
	  double value = std::max(0.0f, std::min(1.0f, (data[i] - min)/(max - min)));
	  // Ensure the value is in [0, 1]
	  unsigned char R = 0, G = 0, B = 0; // Initialize all colors to 0

	  if (value < 1.0 / 6.0) {
		// From blue to aqua (Increase G)
		B = 255;
		G = static_cast<unsigned char>(value * 6 * 127); // 0 to 127
	  } else if (value < 2.0 / 6.0) {
		// From aqua to bright green (Increase G, decrease B)
		G = 127 + static_cast<unsigned char>((value - 1.0 / 6.0) * 6 * 128); // 127 to 255
		B = static_cast<unsigned char>(255 - (value - 1.0 / 6.0) * 6 * 127); // 255 to 128
	  } else if (value < 3.0 / 6.0) {
		// From bright green to pure green (Decrease B)
		G = 255;
		B = static_cast<unsigned char>(128 - (value - 2.0 / 6.0) * 6 * 128); // 128 to 0
	  } else if (value < 4.0 / 6.0) {
		// From green to yellow (Increase R)
		R = static_cast<unsigned char>((value - 3.0 / 6.0) * 6 * 127); // 0 to 127
		G = 255;
	  } else if (value < 5.0 / 6.0) {
		// From yellow to orange (Increase R, decrease G)
		R = 127 + static_cast<unsigned char>((value - 4.0 / 6.0) * 6 * 128); // 127 to 255
		G = static_cast<unsigned char>(255 - (value - 4.0 / 6.0) * 6 * 127); // 255 to 128
	  } else {
		// From orange to red (Decrease G)
		R = 255;
		G = static_cast<unsigned char>(128 - (value - 5.0 / 6.0) * 6 * 128); // 128 to 0
	  }

	  imgData[i * 3] = R;
	  imgData[i * 3 + 1] = G;
	  imgData[i * 3 + 2] = B;
	}
  }
  // Write to PNG
  if (!stbi_write_png(filename.c_str(), sizeX, sizeY, 3, imgData.data(), sizeX * 3)) {
	cerr << "Error: Failed to write PNG: " << filename << endl;
	exit(1);
  }
}

void writeFloatArrayToBinary( string filename,  vector<float>& data, int size) {
  ofstream file(filename, ios::binary);
  if (!file.is_open()) {
    cerr << "Error: Failed to open file " << filename << endl;
    exit(1);
  }
  file.write(reinterpret_cast< char*>(data.data()), size * sizeof(float));
  if (!file) {
    cerr << "Error: Failed to write data to file " << filename << endl;
    file.close();
    exit(1);
  }
  file.close();
}

void writeComplexArrayToBinary(string realFile, string imagFile, vector<Complex>& data, int size) {
  ofstream fileReal(realFile, ios::binary);
  ofstream fileImag(imagFile, ios::binary);

  if (!fileReal.is_open() || !fileImag.is_open()) {
    cerr << "Error: Failed to open files " << realFile << " or " << imagFile << endl;
    exit(1);
  }

  for (int i = 0; i < size; ++i) {
    float realValue = data[i].real();
    float imagValue = data[i].imag();

    fileReal.write(reinterpret_cast< char*>(&realValue), sizeof(float));
    fileImag.write(reinterpret_cast< char*>(&imagValue), sizeof(float));

    if (!fileReal || !fileImag) {
      cerr << "Error: Failed to write data to files " << realFile << " or " << imagFile << endl;
      fileReal.close();
      fileImag.close();
      exit(1);
    }
  }

  fileReal.close();
  fileImag.close();
}

void writeFloatArrayToTxt(string filename, vector<float>& data, int sizeX, int sizeY) {
  ofstream file(filename);

  if (!file.is_open()) {
    cerr << "Error: Failed to open file " << filename << endl;
    exit(1);
  }

  for (int y = 0; y < sizeY; ++y) {
	for (int x = 0; x < sizeX; ++x) {
	  file << data[y*sizeX+x] << " ";
	}
	file << endl;
  }

  if (!file) {
    cerr << "Error: Failed to write data to file " << filename << endl;
    file.close();
	exit(1);
  }

  file.close();
}

void readFloatArrayFromTxt( string filename, vector<float>& data, int size) {
  ifstream file(filename);
  if (!file.is_open()) {
	cerr << "Error: Failed to open file " << filename << endl;
	exit(1);
  }

  // Read the data
  data.clear();
  float value;
  while (file >> value) {
	data.push_back(value);
  }

  if (!file.eof()) {
	cerr << "Error: Failed to read data from file " << filename << endl;
	file.close();
	exit(1);
  }
  file.close();
  
  // Check if the data size matches the expected size
  if ((int)data.size() != size) {
	cerr << "Warning: Data size (" << data.size() << ") does not match the expected size (" << size << ")" << endl;
	 data.resize(size);
  }
}

void readFloatArrayFromBinary( string filename, vector<float>& data, int size) {
  ifstream file(filename, ios::binary);
  if (!file.is_open()) {
    cerr << "Error: Failed to open file " << filename << endl;
	exit(1);
  }
  // Get the file size
  file.seekg(0, ios::end);
  //  streamsize fileSize = file.tellg();
  int fileSize = (int)  file.tellg();
  file.seekg(0, ios::beg);

  // Read the data
  data.resize(fileSize);
  file.read(reinterpret_cast<char*>(data.data()), fileSize);
  if (!file) {
    cerr << "Error: Failed to read data from file " << filename << endl;
    file.close();
	exit(1);
  }
  file.close();
  
   // Check if the data size matches the file size
  if (fileSize != (int)(size * sizeof(float))) {
	cerr << "Warning: File size (" << fileSize << " bytes) does not match the expected data size (" << size * sizeof(float) << " bytes)" << endl;
	 data.resize(size);
  }
  
}

void writeKernelInfo(string filename, int krnSizeX, int krnSizeY, int nk, vector<float>& scales) {
  ofstream outfile(filename);
  if (!outfile) {
	cerr << "Error opening file for writing: " << filename << endl;
	exit(1);
  }

  outfile << "# Kernel Information\n";
  outfile << "## Kernel Size and Count\n";
  outfile << "- Kernel Size: " << krnSizeX << "x" << krnSizeY << "\n";
  outfile << "- Number of Kernels: " << nk << "\n";
  outfile << "## Eigenvalues of Each Kernel\n";

  for (int i = 0; i < nk; ++i) {
	outfile << "Eigenvalue " << i << ": " << scales[i] << "\n";
  }

  outfile.close();
}

void readKernelInfo(string filename, int& krnSizeX, int& krnSizeY, int& nk, vector<float>& scales) {
  ifstream infile(filename);
  if (!infile) {
	cerr << "Error opening file for reading: " << filename << endl;
	exit(1);
  }

  string line;
  while (getline(infile, line)) {
	if (line.find("Kernel Size:") != string::npos) {
	  sscanf(line.c_str(), "- Kernel Size: %dx%d", &krnSizeX, &krnSizeY);
	} else if (line.find("Number of Kernels:") != string::npos) {
	  sscanf(line.c_str(), "- Number of Kernels: %d", &nk);
	} else if (line.find("Eigenvalue") != string::npos) {
	  int index;
	  double value;
	  sscanf(line.c_str(), "Eigenvalue %d: %lf", &index, &value);
	  if (index >= 0 && index < nk) {
		scales.resize(nk);
		scales[index] = value;
	  }
	}
  }

  infile.close();
}
