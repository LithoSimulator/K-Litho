Lx=2048 # Period of mask (x direction)
Ly=2048 # Period of mask (y direction)
maskSizeX=1024 # Size of mask (x direction)
maskSizeY=1024 # Size of mask (x direction)

maskPath='../mask/T1.bin' # filename of the mask to read

krnDir='../K-Litho-TCC/out/kernels/'
# execution
./klitho_socs ${Lx} ${Ly} ${maskSizeX} ${maskSizeY} ${maskPath} ${krnDir}

