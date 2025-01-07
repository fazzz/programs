#bin/sh

FEL_MS=~/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=1.0/mZ=100.00/TZ=1000/SB_KZMAX=1000_NR8_woeljd0.001_2013-08-28/nEX=1000/fq=10/pmf/pmf_pymbar_TAA=300_TCG=300_TZ=1000_KZAAo=1000_KZCGo=0_AA_2013-05-26_20_2

FEL_REMD=~/calspa/refcalc/REMD/AD/s_REVAC_2012-07-23_ff99SB/f300t1200/nEX=100000/freq=1ps/pmf/pmf_py3

FEL_TAMD=~/calspa/refcalc/REUS/AD/s_REUSVAC_2013-08-30_ff99SB/UmbAt_Nbin_4x4_K_0.5/nEX_500/freq_10ps/pmf/pmf_UmbMD_vac_nbx=40_nby=40.txt_2

FEL_REUS=~/calspa/refcalc/REUS/AD/s_REUSVAC_2013-08-30_ff99SB/UmbAt_Nbin_4x4_K_0.5/nEX_500/freq_10ps/pmf/pmf_UmbMD_vac_nbx=40_nby=40.txt_2

FEL_MS_GBSA=~/calspa/InV2InW/AD/e_GBSA_2012-08-06/e_CG-FG_NH_2012-07-27//tau=1.0/mZ=100.00/TZ=1000/SB_KZMAX=1000_NR8_woeljd0.001_2013-08-28/nEX=1000/fq=10/pmf/pmf_TAA=300_TCG_300_TZ_1000_0.3_0_1000_CG_InW

FEL_MS2=~/papers/CG-FG_TACCM_REMD/FEL_MS_2013-09-04
FEL_REMD2=~/papers/CG-FG_TACCM_REMD/FEL_REMD_2013-09-04
FEL_TAMD2=~/papers/CG-FG_TACCM_REMD/FEL_TAMD_2013-09-04
FEL_REUS2=~/papers/CG-FG_TACCM_REMD/FEL_REUS_2013-09-04
FEL_MS_GBSA2=~/papers/CG-FG_TACCM_REMD/FEL_MS_GBSA_2013-09-04

mim_FEL_MS=~/papers/CG-FG_TACCM_REMD/min_FEL_MS_2013-09-04
mim_FEL_REMD=~/papers/CG-FG_TACCM_REMD/min_FEL_REMD_2013-09-04
mim_FEL_TAMD=~/papers/CG-FG_TACCM_REMD/min_FEL_TAMD_2013-09-04
mim_FEL_REUS=~/papers/CG-FG_TACCM_REMD/min_FEL_REUS_2013-09-04
mim_FEL_MS_GBSA=~/papers/CG-FG_TACCM_REMD/min_FEL_MS_GBSA_2013-09-04

mim_FEL_MS2=~/papers/CG-FG_TACCM_REMD/min_FEL_MS_s_2013-09-04
mim_FEL_REMD2=~/papers/CG-FG_TACCM_REMD/min_FEL_REMD_s_2013-09-04
mim_FEL_TAMD2=~/papers/CG-FG_TACCM_REMD/min_FEL_TAMD_s_2013-09-04
mim_FEL_REUS2=~/papers/CG-FG_TACCM_REMD/min_FEL_REUS_s_2013-09-04
mim_FEL_MS_GBSA2=~/papers/CG-FG_TACCM_REMD/min_FEL_MS_GBSA_s_2013-09-04

cp ${FEL_MS} ${FEL_MS2}
cp ${FEL_REMD} ${FEL_REMD2}
cp ${FEL_TAMD} ${FEL_TAMD2}
cp ${FEL_REUS} ${FEL_REUS2}
cp ${FEL_MS_GBSA} ${FEL_MS_GBSA2}

sed -i -e "s/\*\*\*\*\*\*/100/g" ${FEL_MS2}
sed -i -e "s/\*\*\*\*\*\*/100/g" ${FEL_REMD2}
sed -i -e "s/\*\*\*\*\*\*/100/g" ${FEL_TAMD2}
sed -i -e "s/\*\*\*\*\*\*/100/g" ${FEL_REMD2}
sed -i -e "s/\*\*\*\*\*\*/100/g" ${FEL_MS_GBSA2}

~/mybin/find_mimg 20 20 ${FEL_MS2} > ${mim_FEL_MS}	
~/mybin/find_mimg 30 30 ${FEL_REMD2} > ${mim_FEL_REMD}
~/mybin/find_mimg 40 40   ${F40_TAMD2} > ${mim_FEL_TAMD}
~/mybin/find_mimg 40 40 ${FEL_REUS2} > ${mim_FEL_REUS}
~/mybin/find_mimg 21 21 ${FEL_MS_GBSA2} > ${mim_FEL_MS_GBSA}

sort -n -k4  ${mim_FEL_MS} > ${mim_FEL_MS2}
sort -n -k4  ${mim_FEL_REMD} > ${mim_FEL_REMD2}
sort -n -k4  ${mim_FEL_TAMD} > ${mim_FEL_TAMD2}
sort -n -k4  ${mim_FEL_MS_GBSA} > ${mim_FEL_MS_GBSA2}

