#!/bin/sh

dir=~/work/programs
srcdir=${dir}/MD_pep/src

files=( ${dir}/MD/MD_pep_NH_MP1996_AAFF_Amber.c \
    ${dir}/MD/MD_NHC_MP1996.c ${dir}/MD/MD_NHC_MP1996.h \
    ${dir}/MD/MD.c            ${dir}/MD/MD.h \
    ${dir}/UMBP/UMBP.c ${dir}/UMBP/UMBP.h \
    ${dir}/efunc/efunc.c                  ${dir}/efunc/EF.h \
    ${dir}/LA/LA.c                        ${dir}/LA/LA.h \
    ${dir}/readParmtop/readParmtopL.c     ${dir}/readParmtop/PTL.h \
    ${dir}/CFF/calcFFL.c                  ${dir}/CFF/FFL.h \
    ${dir}/TOPO/TOPO.c                    ${dir}/TOPO/TOPO.h \
    ${dir}/RAND/BOXMULL.c                 ${dir}/RAND/BOXMULL.h \
    ${dir}/RAND/mt19937ar.c               ${dir}/RAND/BOXMULL.h \
    ${dir}/netcdf/netcdf_mineL.c          ${dir}/netcdf/netcdf_mineL.h \
    ${dir}/ABA_makefile/MB.c              ${dir}/ABA_makefile/MB.h \
    ${dir}/MYMATH/mymath.c                ${dir}/MYMATH/mymath.h \
    ${dir}/GOLM/GOLMAA_PROTEINS2008.c     ${dir}/GOLM/GOLMAA_PROTEINS2008.h \
    ${dir}/GOLM/GOLMAA_PROTEINS2008_set.c ${dir}/GOLM/GOLMAA_PROTEINS2008_set.h )

liblaries=( ~/software/netcdf-4.1.1 )

cp ${files[*]} ${srcdir}

cp -r ${liblaries[*]} ${srcdir}

cd ${srcdir}

mv ${srcdir}/efunc.c ${srcdir}/EF.c
mv ${srcdir}/readParmtopL.c ${srcdir}/PTL.c
mv ${srcdir}/calcFFL.c ${srcdir}/FFL.c
mv ${srcdir}/calcFF.c ${srcdir}/FF.c
