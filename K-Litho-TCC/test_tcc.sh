#!/bin/bash

Lx=2048 # Period of mask (x direction)
Ly=2048 # Period of mask (y direction)
maskSizeX=1024 # Size of mask (x direction)
maskSizeY=1024 # Size of mask (x direction)

maskPath='../mask/T1.bin' # filename of the mask to read

# source parameter
srcSize=201 # grid size of source
srcTypeList=(Annular CrossQuadrupole Dipole Point Import)
srcType=${srcTypeList[0]}

## for annular source
InnerRadius=0.6
OuterRadius=0.9

## for quadrupole, dipole source
Radius=0.1
Offset=0.5

## for dipole source
OnXAxis=1

## for point source
ptX=0.6
ptY=0.2

InputSourcePath="./source/src_test1_size101.bin"

NA=0.83
Defocus=0.5 # Defocus normalized by dividing by wavelength/(NA^2)
NumKernels=3

# execution
./klitho_tcc ${Lx} ${Ly} ${maskSizeX} ${maskSizeY} ${maskPath} ${srcSize} ${srcType} ${InnerRadius} ${OuterRadius} ${Radius} ${Offset} ${OnXAxis} ${ptX} ${ptY} ${InputSourcePath} ${NA} ${Defocus} ${NumKernels}
