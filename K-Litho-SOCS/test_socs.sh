#!/bin/bash

# Mask parameter
Lx=2048 # Period of mask (x direction)
Ly=2048 # Period of mask (y direction)
maskSizeX=1024 # Size of mask (x direction)
maskSizeY=1024 # Size of mask (x direction)

maskTypeList=(LineSpace Import)
maskType=${maskTypeList[1]}  # Select mask type (0: LineSpace, 1: Import file)

## for Line/Space pattern mask
LineWidth=128
SpaceWidth=128
isHorizontal=0 # 0: vertical direction, 1: horizontal direction

## for import mask file
MaskFile='../mask/T1.bin' # filename of the mask to read

# Kernel
krnDir='../K-Litho-TCC/out/kernels/' # Directory containing kernel files

# Execution
./klitho_socs ${Lx} ${Ly} ${maskSizeX} ${maskSizeY}  ${maskType} ${LineWidth} ${SpaceWidth} ${isHorizontal} ${MaskFile} ${krnDir}

