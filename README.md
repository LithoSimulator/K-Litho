# K-Litho
K-Litho is a toolset for simulating optical lithography.

K-Litho includes two tools: K-Litho-TCC and K-Litho-SOCS.

- K-Litho-TCC: Calculates the optical image using the Transmission Cross Coefficient (TCC) and derives the Sum of Coherent Systems (SOCS) kernels. K-Litho-TCC allows you to set optical conditions such as mask transmittance, source shape and intensity data, NA, and defocus.

- K-Litho-SOCS: Calculates the optical image using the Sum of Coherent Systems (SOCS) kernels derived by K-Litho-TCC.

## Prerequisites

Before implementing K-Litho, make sure that [FFTW3](http://www.fftw.org/) and [LAPACK](http://www.netlib.org/lapack/) are installed.

If these packages are not installed, install them using the appropriate commands for your Linux distribution:

Below are the installation methods for Ubuntu and CentOS.

For Ubuntu/Debian:
```
sudo apt-get install liblapack-dev libfftw3-dev
```

For CentOS/RHEL:
```
sudo yum install lapack-devel fftw-devel
```

## K-Litho-TCC

Calculates the optical image using TCC and derives SOCS kernels.

You can set optical conditions such as mask, source, NA, and defocus.

### Compilation

Move to the K-Litho-TCC directory and compile using the make command.
```
cd K-Litho-TCC
make
```

After compilation, an executable file named "klitho_tcc" will be created.

### Execution

You can specify various simulation parameters when running the program.

You can also use text or binary files containing mask transmittance data and source intensity data as input files.

These input files must contain data in Float format.

The "mask" directory contains transmittance data of 1024x1024 for 10 test benches (T1 to T10) provided by ICCAD2013 as examples of input mask data.

Without the need for file input, you can utilize masks with Line and Space patterns, as well as parametric sources such as Annular or Dipole shapes by adjusting their parameters.

A sample script for running K-Litho-TCC is provided in "test_tcc.sh". When running this script for the first time, grant execution permission with the following command:
```
chmod +x test_tcc.sh
```

Then, run it with the following command:
```
./test_tcc.sh
```

### Output Files

This program outputs multiple files. by default, the following files are output:

- `./out/mask.png`: PNG file of the input mask
- `./out/src.png`: PNG file of the source used
- `./out/tcc_r.png`, `./out/tcc_i.png`: PNG files of the real/imaginary parts of the calculated TCC matrix
- `./out/image.bin`, `./out/image.png`: binary/PNG files of the calculated optical image
- `./out/kernels/kernel_info.txt`: Text file containing information on the number, size, and eigenvalues of SOCS kernels
- `./out/kernels/krn_*_r.bin`, `./out/kernels/krn_*_i.bin`: binary files of the real/imaginary parts of the eigenfunctions of SOCS kernels
- `./out/kernels/png/krn_*_r.png`, `./out/kernels/png/krn_*_i.png`: PNG files of the real/imaginary parts of the eigenfunctions of SOCS kernels

## K-Litho-SOCS

Calculates the optical image using SOCS. Inputs the mask and kernels and calculates the optical image at high speed.

You can use the kernel files created by K-Litho-TCC as input kernels.

### Compilation

Move to the K-Litho-SOCS directory and compile using the make command.
```
cd K-Litho-SOCS
make
```

After compilation, an executable file named "klitho_socs" will be created.

### Execution

You can specify simulation parameters when running the program. A sample script for running K-Litho-SOCS is provided in "test_socs.sh".

A sample script for running K-Litho-SOCS is provided in "test_socs.sh". When running this script for the first time, grant execution permission with the following command:
```
chmod +x test_socs.sh
```

Then, run it with the following command:
```
./test_socs.sh
```

### Output Files

This program outputs multiple files. by default, the following files are output:

- `./out/mask.png`: PNG file of the input mask
- `./out/image.bin`, `./out/image.png`: binary/PNG files of the calculated optical image


## Publications

- Masaki KURAMOCHI, Yukihide KOHIRA, Hiroyoshi TANABE, Tetsuaki MATSUNAWA,  Chikaaki KODAMA,
"Development of a Lithography Simulation Tool Set in Various Optical Conditions for Source Mask Optimization",
IEEE Access, vol. 12, pp. 58490-58501, 2024. (doi: 10.1109/ACCESS.2024.3390936)  (<a href="https://ieeexplore.ieee.org/document/10504808">IEEE Xplore</a>, <a href="https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=10504808">paper</a>)
