#!/bin/sh

dir=~/work/programs/
srcdir=${dir}/ABAMD_2014_05_01/src/

files=( ${dir}/ABA/mainABAMD_NH_new_2014-01-09_force_checker_2014-01-29.c \
    ${dir}/ABA/ABA.c                     ${dir}/ABA/ABA.h \
    ${dir}/ABA/ABAb.c                    ${dir}/ABA/ABAb.h \
    ${dir}/ABA/ABA_prepass.c             ${dir}/ABA/ABA_prepass.h \
    ${dir}/ABA/ABA_mainpass.c            ${dir}/ABA/ABA_mainpass.h \
    ${dir}/ABA/ABA_backpass.c            ${dir}/ABA/ABA_backpass.h \
    ${dir}/ABA/ABA_Nose-Hoover.c         ${dir}/ABA/ABA_Nose-Hoover.h \
    ${dir}/ABA/ABA_Nose-Hoover_chain.c   ${dir}/ABA/ABA_Nose-Hoover_chain.h \
    ${dir}/ABA/ABA_Nose-Hoover_new.c     ${dir}/ABA/ABA_Nose-Hoover_new.h \
    ${dir}/ABA/ABA_Nose-Hoover_new_mvV.c ${dir}/ABA/ABA_Nose-Hoover_new_mvV.h \
    ${dir}/ABA/ABA_pick_data.c           ${dir}/ABA/ABA_pick_data.h \
    ${dir}/ABA/ABA_set_frc.c             ${dir}/ABA/ABA_set_frc.h \
    ${dir}/ABA/ABA_set_imat.c            ${dir}/ABA/ABA_set_imat.h \
    ${dir}/ABA/ABA_set_imatb.c           ${dir}/ABA/ABA_set_imatb.h \
    ${dir}/ABA/ABA_set_lref.c            ${dir}/ABA/ABA_set_lref.h \
    ${dir}/ABA/ABA_set_rst.c             ${dir}/ABA/ABA_set_rst.h \
    ${dir}/ABA/ABA_set_tmat.c            ${dir}/ABA/ABA_set_tmat.h \
    ${dir}/ABA/ABA_set_trans.c           ${dir}/ABA/ABA_set_trans.h \
    ${dir}/ABA/ABA_integ.c               ${dir}/ABA/ABA_integ.h \
    ${dir}/ABA/ABA_update.c              ${dir}/ABA/ABA_update.h \
    ${dir}/ABA/ABA_calcattfrc.c          ${dir}/ABA/ABA_calcattfrc.h \
    ${dir}/ABA/ABA_calcKineE.c           ${dir}/ABA/ABA_calcKineE.h \
    ${dir}/ABA/ABA_Inverse.c             ${dir}/ABA/ABA_Inverse.h \
    ${dir}/ABA/ABA_Inverse_mainpass.c    ${dir}/ABA/ABA_Inverse_mainpass.h \
    ${dir}/ABA/ABA_Inverse_backpass.c    ${dir}/ABA/ABA_Inverse_backpass.h \
    ${dir}/ABA/ABA_gtree.c               ${dir}/ABA/ABA_gtree.h \
    ${dir}/ABA/ABA_hosoku.c              ${dir}/ABA/ABA_hosoku.h \
    ${dir}/ABA/quaternion.c              ${dir}/ABA/quaternion.h \
    ${dir}/efunc/efunc.c                 ${dir}/efunc/EF.h \
    ${dir}/LA/LA.c                       ${dir}/LA/LA.h \
    ${dir}/readParmtop/readParmtopL.c    ${dir}/readParmtop/PTL.h \
    ${dir}/CFF/calcFFL.c                 ${dir}/CFF/FFL.h \
    ${dir}/TOPO/TOPO.c                   ${dir}//TOPO/TOPO.h \
    ${dir}/RAND/BOXMULL.c                ${dir}/RAND/BOXMULL.h \
    ${dir}/RAND/mt19937ar.c              ${dir}/RAND/BOXMULL.h \
    ${dir}/netcdf/netcdf_mineL.c         ${dir}/netcdf/netcdf_mineL.h \
#    ${dir}/MolBas/MB.c                   ${dir}/MolBas/MB.h \
    ${dir}/ABA_makefile/MB.c             ${dir}//ABA_makefile/MB.h \
    ${dir}/MYMATH/mymath.c               ${dir}/MYMATH/mymath.h \
    ${dir}/SBFF/SBFF.c                   ${dir}/SBFF/SBFF.h \
    ${dir}/CFF/calcFF.c                  ${dir}/CFF/FF.h \
    ${dir}/readParmtop/readParmtop.c     ${dir}/readParmtop/PT.h )

liblaries=( ~/software/netcdf-4.1.1 ~/software/CLAPACK-3.2.1 ~/software/gc6.7 )

echo ${files[*]}

cp ${files[*]} ${srcdir}

echo ${liblaries[*]}

cp -r ${liblaries[*]} ${srcdir}

cd ${srcdir}

mv ${srcdir}/mainABAMD_NH_new_2014-01-09_force_checker_2014-01-29.c ${srcdir}/mainABAMD.c
mv ${srcdir}/efunc.c ${srcdir}/EF.c
mv ${srcdir}/readParmtopL.c ${srcdir}/PTL.c
mv ${srcdir}/readParmtop.c ${srcdir}/PT.c
mv ${srcdir}/calcFFL.c ${srcdir}/FFL.c
mv ${srcdir}/calcFF.c ${srcdir}/FF.c
