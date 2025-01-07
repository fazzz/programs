#bin/sh

FEL_MS=~/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=1.0/mZ=100.00/TZ=700/SB_KZMAX=500_NR4_woeljd0.001/nEX=1000/fq=10/pmf/pmf_pymbar_TAA=300_TCG=300_TZ=700_KZAAo=500_KZCGo=0_AA_2013-05-26_20_2

FEL_REMD=~/calspa/refcalc/REMD/AD/s_REVAC_2012-07-23_ff99SB/f300t600_rep4/nEX=100000/freq=1ps/pmf/pmf_py3

FEL_TAMD=~/calspa/TACCM/AD/e_TACCM_NH_2012-07-31_99SB/tau=1.0/TB=700/KZ=1000/mZ=100.00/pmf/pmf_T=300.Zhist_0.2

FEL_REUS=~/calspa/refcalc/REUS/AD/s_REUSVAC_2013-08-12_ff99SB/UmbAt_Nbin_4x4_K_1/nEX_250/freq_10ps/pmf/pmf_UmbMD_vac_nbx=20_nby=20.txt

FEL_MS2=~/papers/CG-FG_TACCM_REMD/FEL_MS_2013-08-22
FEL_REMD2=~/papers/CG-FG_TACCM_REMD/FEL_REMD_2013-08-22
FEL_TAMD2=~/papers/CG-FG_TACCM_REMD/FEL_TAMD_2013-08-22
FEL_REUS2=~/papers/CG-FG_TACCM_REMD/FEL_REUS_2013-08-22

mim_FEL_MS=~/papers/CG-FG_TACCM_REMD/min_FEL_MS_2013-08-22
mim_FEL_REMD=~/papers/CG-FG_TACCM_REMD/min_FEL_REMD_2013-08-22
mim_FEL_TAMD=~/papers/CG-FG_TACCM_REMD/min_FEL_TAMD_2013-08-22
mim_FEL_REUS=~/papers/CG-FG_TACCM_REMD/min_FEL_REUS_2013-08-22

mim_FEL_MS2=~/papers/CG-FG_TACCM_REMD/min_FEL_MS_s_2013-08-22
mim_FEL_REMD2=~/papers/CG-FG_TACCM_REMD/min_FEL_REMD_s_2013-08-22
mim_FEL_TAMD2=~/papers/CG-FG_TACCM_REMD/min_FEL_TAMD_s_2013-08-22
mim_FEL_REUS2=~/papers/CG-FG_TACCM_REMD/min_FEL_REUS_s_2013-08-22

cp ${FEL_MS} ${FEL_MS2}
cp ${FEL_REMD} ${FEL_REMD2}
cp ${FEL_TAMD} ${FEL_TAMD2}
cp ${FEL_REUS} ${FEL_REUS2}

sed -i -e "s/\*\*\*\*\*\*/100/g" ${FEL_MS2}
sed -i -e "s/\*\*\*\*\*\*/100/g" ${FEL_REMD2}
sed -i -e "s/\*\*\*\*\*\*/100/g" ${FEL_TAMD2}
sed -i -e "s/\*\*\*\*\*\*/100/g" ${FEL_REMD2}

~/mybin/find_mimg 20 20 ${FEL_MS2} > ${mim_FEL_MS}	
~/mybin/find_mimg 20 20 ${FEL_REMD2} > ${mim_FEL_REMD}
~/mybin/find_mimg 32 32   ${FEL_TAMD2} > ${mim_FEL_TAMD}
~/mybin/find_mimg 20 20 ${FEL_REUS2} > ${mim_FEL_REUS}

sort -n -k4  ${mim_FEL_MS} > ${mim_FEL_MS2}
sort -n -k4  ${mim_FEL_REMD} > ${mim_FEL_REMD2}
sort -n -k4  ${mim_FEL_TAMD} > ${mim_FEL_TAMD2}

