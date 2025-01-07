#!/bin/sh

files=( ~/work/programs/efunc/EF.h ~/work/programs/efunc/efunc.c \
~/work/programs/LA/LA.h ~/work/programs/LA/LA.c \
~/work/programs/readParmtop/PTL.h ~/work/programs/readParmtop/readParmtopL.c \
~/work/programs/CFF/FFL.h ~/work/programs/CFF/calcFFL.c \
~/work/programs/TOPO/TOPO.h ~/work/programs/TOPO/TOPO.c \
~/work/programs/RAND/RAND.h \
~/work/programs/RAND/BOXMULL.h ~/work/programs/RAND/BOXMULL.c \
~/work/programs/ABA/quaternion.h ~/work/programs/ABA/quaternion.c \
~/work/programs/netcdf/netcdf_mineL.h ~/work/programs/netcdf/netcdf_mineL.c )

grep -nH -e include ${files[*]}
