#bin/sh

FEL_REMD=~/calspa/refcalc/REMD/AD/s_REVAC_2012-07-23_ff99SB/f300t1000/nEX=100000/freq=1ps/pmf/pmf_ADv_0.01
FEL_REUS=/home/yamamori/calspa/refcalc/REUS/AD/s_REUSVAC_2013-08-12_ff99SB/UmbAt_Nbin_4x4_K_1/nEX_500/freq_10ps/pmf/pmf_UmbMD_vac_nbx=40_nby=40.txt

FEL_REMD2=~/papers/CG-FG_TACCM_REMD/FEL_REMD_2
FEL_REUS2=~/papers/CG-FG_TACCM_REMD/FEL_REUS

mim_FEL_REMD2=~/papers/CG-FG_TACCM_REMD/min_FEL_REMD_s_2
mim_FEL_REUS=~/papers/CG-FG_TACCM_REMD/min_FEL_REUS

mim_FEL_REMD2=~/papers/CG-FG_TACCM_REMD/min_FEL_REMD_s_2
mim_FEL_REUS2=~/papers/CG-FG_TACCM_REMD/min_FEL_REUS_s

cp ${FEL_REMD} ${FEL_REMD2}
cp ${FEL_REUS} ${FEL_REUS2}

sed -i -e "s/\*\*\*\*\*\*/100/g" ${FEL_REMD2}
sed -i -e "s/inf/100/g" ${FEL_REUS2}

~/mybin/find_mimg 40 40 ${FEL_REMD2} > ${mim_FEL_REMD}
~/mybin/find_mimg 40 40 ${FEL_REUS2} > ${mim_FEL_REUS}

sort -n -k4  ${mim_FEL_REMD} > ${mim_FEL_REMD2}
sort -n -k4  ${mim_FEL_REUS} > ${mim_FEL_REUS2}

