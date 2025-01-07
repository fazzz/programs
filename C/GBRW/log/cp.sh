#!/bin/sh

cp ~/work/programs/MGaussian_reweight_oneD/src/main_oneD_Gaussian_MP.c ~/work/programs/GBRW/src/main.c
cp ~/work/programs/MGaussian_reweight_oneD/src/Simpson_integ_oneD_Gaussianbase.c ~/work/programs/GBRW/src/Simpson.c
cp ~/work/programs/MGaussian_reweight_oneD/src/EF.c ~/work/programs/GBRW/src/EF.c
cp ~/work/programs/MGaussian_reweight_oneD/src/Simpson_integ_oneD_Gaussianbase.h ~/work/programs/GBRW/src/Simpson.h
cp ~/work/programs/MGaussian_reweight_oneD/src/EF.h ~/work/programs/GBRW/src/EF.h
cp ~/work/programs/MGaussian_reweight_oneD/src/Makefile2 ~/work/programs/GBRW/src/Makefile

cp ~/work/programs/MGaussian_reweight_oneD/test/input/MGaussian_reweight_36_2.in ~/work/programs/GBRW/test/input
cp ~/work/programs/MGaussian_reweight_oneD/test/input/umb_36_2.in ~/work/programs/GBRW/test/input
