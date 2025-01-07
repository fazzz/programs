#bin/sh

FEL_MS_GBSA=~/calspa/InV2InW/AD/e_GBSA_2012-08-06/e_CG-FG_NH_2012-07-27/tau=1.0/mZ=100.00/TZ=700/SB_KZMAX=500_NR4_woeljd0.001/nEX=1000/fq=10/pmf/pmf_TAA=300_TCG_300_TZ_700_0.3_500_0_AA_InW
FEL_REMD_GBSA=~/calspa/refcalc/REMD/AD/s_REGB_2012-07-23_ff99SB/f300t500_rep4/nEX=100000/freq=1ps/pmf/pmf_ADg

FEL_MS_GBSA2=~/papers/CG-FG_TACCM_REMD/FEL_MS_GBSA_2013-08-24
FEL_REMD_GBSA2=~/papers/CG-FG_TACCM_REMD/FEL_REMD_GBSA_2013-08-24

min_FEL_MS_GBSA=~/papers/CG-FG_TACCM_REMD/min_FEL_MS_GBSA_2013-08-24
min_FEL_REMD_GBSA=~/papers/CG-FG_TACCM_REMD/min_FEL_REMD_GBSA_2013-08-24

min_FEL_MS_GBSA2=~/papers/CG-FG_TACCM_REMD/min_FEL_MS_GBSA_s_2013-08-24
min_FEL_REMD_GBSA2=~/papers/CG-FG_TACCM_REMD/min_FEL_REMD_GBSA_s_2013-08-24

cp ${FEL_MS_GBSA} ${FEL_MS_GBSA2}
cp ${FEL_REMD_GBSA} ${FEL_REMD_GBSA2}

sed -i -e "s/\*\*\*\*\*\*/100/g" ${FEL_MS_GBSA2}
sed -i -e "s/\*\*\*\*\*\*/100/g" ${FEL_REMD_GBSA2}

~/mybin/find_mimg 21 21  ${FEL_MS_GBSA2} > ${min_FEL_MS_GBSA}
~/mybin/find_mimg 21 21  ${FEL_REMD_GBSA2} > ${min_FEL_REMD_GBSA}

sort -n -k4  ${mim_FEL_MS_GBSA} > ${min_FEL_MS_GBSA2}
sort -n -k4  ${mim_FEL_REMD_GBSA} > ${min_FEL_REMD_GBSA2}



