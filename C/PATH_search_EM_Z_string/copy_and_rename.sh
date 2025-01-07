#!/bin/sh

dir=~/work/programs/
srcdir=${dir}/PATH_search_EM_Z_string/src/

files=( ${dir}/PATH_SE_MUSTERMD/main_path_search_MGaussian_calc_KLdiv_detailed_FEL_non_square_range.c \
    ${dir}/EMalg/EMalg.c ${dir}/EMalg/EMalg.h \
    ${dir}/EMalg/EMalg_non_square_range.c ${dir}/EMalg/EMalg_non_square_range.h \
    ${dir}/Gaussian/Gaussian.c ${dir}/Gaussian/Gaussian.h \
    ${dir}/K_means/K_means.c ${dir}/K_means/K_means.h \
    ${dir}/K_means/K_means_non_square_range.c ${dir}/K_means/K_means_non_square_range.h \
    ${dir}/STRING/zp_string.c ${dir}/STRING/STRING.h \
    ${dir}/STRING/zp_string_CV.c ${dir}/STRING/STRING_CV.h \
    ${dir}/CSI/cubic_spline_interpolation.c ${dir}/CSI/CSI.h \
    ${dir}/efunc/efunc.c ${dir}/efunc/EF.h \
    ${dir}/CFF/calcFF.c                 ${dir}/CFF/FF.h \
    ${dir}/TOPO/TOPO.c                  ${dir}/TOPO/TOPO.h \
    ${dir}/ABA_makefile/MB.c            ${dir}/ABA_makefile/MB.h \
    ${dir}/MYMATH/mymath.c              ${dir}/MYMATH/mymath.h \
    ${dir}/LA/LA.c                      ${dir}/LA/LA.h \
    ${dir}/readParmtop/readParmtop.c    ${dir}/readParmtop/PT.h \
    ${dir}/IO/inoutdata.c ${dir}/IO/IO.h )

echo ${files[*]}

cp ${files[*]} ${srcdir}

mv ${srcdir}/main_path_search_MGaussian_calc_KLdiv_detailed_FEL_non_square_range.c ${srcdir}/main.c
mv ${srcdir}/cubic_spline_interpolation.c ${srcdir}/CSI.c
mv ${srcdir}/zp_string.c ${srcdir}/STRING.c
mv ${srcdir}/zp_string_CV.c ${srcdir}/STRING_CV.c
mv ${srcdir}/efunc.c ${srcdir}/EF.c
mv ${srcdir}/calcFF.c ${srcdir}/FF.c
mv ${srcdir}/readParmtop.c ${srcdir}/PT.c
mv ${srcdir}/inoutdata.c ${srcdir}/IO.c
