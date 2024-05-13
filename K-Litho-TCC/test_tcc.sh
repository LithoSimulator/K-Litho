#!/bin/bash

# Mask parameter
Lx=2048  # Period of mask (x direction)
Ly=2048 # Period of mask (y direction)
maskSizeX=1024 # Size of mask (x direction)
maskSizeY=1024 # Size of mask (x direction)

maskTypeList=(LineSpace Import)
maskType=${maskTypeList[1]} # Select mask type (0: LineSpace, 1: Import file)

## for Line/Space pattern mask
LineWidth=128
SpaceWidth=128
isHorizontal=0 # 0: vertical direction, 1: horizontal direction

## for import mask file
MaskFile='../mask/T1.bin'  # Filename of the mask to read

# Source parameter
srcSize=101 # Grid size of source
srcTypeList=(Annular Dipole CrossQuadrupole Point Import)
srcType=${srcTypeList[0]} # Select source type (0: Annular, 1: Dipole, 2: CrossQuadrupole, 3: Point, 4: Import)

## for annular source
InnerRadius=0.6 # [0:1]
OuterRadius=0.9 # [0:1]

## for quadrupole, dipole source
Radius=0.1 # [0:1]
Offset=0.5 # [0:1]

## for dipole source
OnXAxis=1 # 0: dipole on Y-axis, 1: dipole on X-axis

## for point source
ptX=0.6 # [-1:1]
ptY=0.2 # [-1:1]

## for import source file
SourceFile="./source/src_test1_size101.bin"  # Filename of the source to read

# Other parameter
NA=0.8 # [0:1]
Defocus=0.2 # Defocus normalized by dividing by wavelength/(NA^2)
NumKernels=3 # Number of output kernels

# Execution
./klitho_tcc ${Lx} ${Ly} ${maskSizeX} ${maskSizeY} ${maskType} ${LineWidth} ${SpaceWidth} ${isHorizontal} ${MaskFile} ${srcSize} ${srcType} ${InnerRadius} ${OuterRadius} ${Radius} ${Offset} ${OnXAxis} ${ptX} ${ptY} ${SourceFile} ${NA} ${Defocus} ${NumKernels}
