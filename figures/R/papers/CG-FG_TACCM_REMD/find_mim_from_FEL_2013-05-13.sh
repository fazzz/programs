#bin/sh

FEL_MS=~/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=1.0/mZ=100.00/TZ=750/SB_KZMAX=5000_NR8-2_woeljd0.001/nEX=1000/fq=10/pmf/pmf_TAA=300_TCG_300_TZ_750_0.01_0_5000_CG_4_2012-08-21

FEL_REMD=~/calspa/refcalc/REMD/AD/s_REVAC_2012-07-23_ff99SB/f300t400/nEX=100000/freq=1ps/pmf/pmf_ADv_0.01

FEL_US=~/calspa/refcalc/UmbSam/AD/s_UmbSam_vac_2012-11-12_ff99SB/Umb_Nbin=12x12_K=10/pmf/pmf_UmbMD_vac_10ns_nbx=628_nby=628.txt

FEL_TAMD=~/calspa/TACCM/AD/e_TACCM_NH_2012-07-31_99SB/tau=1.0/TB=750/KZ=1000/mZ=50000.00/pmf/pmf_T=300.Zhist_0.1

FEL_MS_GBSA=~/calspa/InV2InW/AD/e_GBSA_2012-08-06/e_CG-FG_NH_2012-07-27/tau=1.0/mZ=50000.00/TZ=750/T1_wrefd0.1/nEX=1000/fq=10/pmf/pmf_TAA=300_TCG_300_TZ_750_0.3_1000_0_AA_InW

FEL_REMD_GBSA=~/calspa/refcalc/REMD/AD/s_REGB_2012-07-23_ff99/f300t400/nEX=100000/freq=1ps/pmf/pmf_ADg_radi

FEL_MS_ME=~/calspa/TACCM_CGAAREMD/MetEnk/e_CG-FG_CMD_NH_2012-06-04/1FG2CG/tau=1.0/mZ=50000.00/ep=0.2/cutoff=4.7/TZ=1400/2CG21_KAA=1000/freq=1/pmf/pmf_TAA=300_TCG1_370_TCG2_370_TZ_1400_0.01_1000_0_0_AA_bo10000ps_2

FEL_MS2=~/papers/CG-FG_TACCM_REMD/FEL_MS
FEL_REMD2=~/papers/CG-FG_TACCM_REMD/FEL_REMD
FEL_US2=~/papers/CG-FG_TACCM_REMD/FEL_US
FEL_TAMD2=~/papers/CG-FG_TACCM_REMD/FEL_TAMD
FEL_MS_GBSA2=~/papers/CG-FG_TACCM_REMD/FEL_MS_GBSA
FEL_REMD_GBSA2=~/papers/CG-FG_TACCM_REMD/FEL_REMD_GBSA
FEL_MS_ME2=~/papers/CG-FG_TACCM_REMD/FEL_MS_ME

mim_FEL_MS=~/papers/CG-FG_TACCM_REMD/min_FEL_MS
mim_FEL_REMD=~/papers/CG-FG_TACCM_REMD/min_FEL_REMD
mim_FEL_US=~/papers/CG-FG_TACCM_REMD/min_FEL_US
mim_FEL_TAMD=~/papers/CG-FG_TACCM_REMD/min_FEL_TAMD
mim_FEL_MS_GBSA=~/papers/CG-FG_TACCM_REMD/min_FEL_MS_GBSA
min_FEL_REMD_GBSA=~/papers/CG-FG_TACCM_REMD/min_FEL_REMD_GBSA
min_FEL_MS_ME=~/papers/CG-FG_TACCM_REMD/min_FEL_MS_ME

mim_FEL_MS2=~/papers/CG-FG_TACCM_REMD/min_FEL_MS_s
mim_FEL_REMD2=~/papers/CG-FG_TACCM_REMD/min_FEL_REMD_s
mim_FEL_US2=~/papers/CG-FG_TACCM_REMD/min_FEL_US_s
mim_FEL_TAMD2=~/papers/CG-FG_TACCM_REMD/min_FEL_TAMD_s
mim_FEL_MS_GBSA2=~/papers/CG-FG_TACCM_REMD/min_FEL_MS_GBSA_s
min_FEL_REMD_GBSA2=~/papers/CG-FG_TACCM_REMD/min_FEL_REMD_GBSA_s
min_FEL_MS_ME2=~/papers/CG-FG_TACCM_REMD/min_FEL_MS_ME_s

cp ${FEL_MS} ${FEL_MS2}
cp ${FEL_REMD} ${FEL_REMD2}
cp ${FEL_US} ${FEL_US2}
cp ${FEL_TAMD} ${FEL_TAMD2}
cp ${FEL_MS_GBSA} ${FEL_MS_GBSA2}
cp ${FEL_REMD_GBSA} ${FEL_REMD_GBSA2}
cp ${FEL_MS_ME} ${FEL_MS_ME2}

sed -i -e "s/\*\*\*\*\*\*/100/g" ${FEL_MS2}
sed -i -e "s/\*\*\*\*\*\*/100/g" ${FEL_REMD2}
sed -i -e "s/\*\*\*\*\*\*/100/g" ${FEL_US2}
sed -i -e "s/\*\*\*\*\*\*/100/g" ${FEL_TAMD2}
sed -i -e "s/\*\*\*\*\*\*/100/g" ${FEL_MS_GBSA2}
sed -i -e "s/\*\*\*\*\*\*/100/g" ${FEL_REMD_GBSA2}
sed -i -e "s/\*\*\*\*\*\*/100/g" ${FEL_MS_ME2}

~/mybin/find_mimg 629 629 ${FEL_MS2} > ${mim_FEL_MS}	
~/mybin/find_mimg 629 629 ${FEL_REMD2} > ${mim_FEL_REMD}
~/mybin/find_mimg 628 628 ${FEL_US2} > ${mim_FEL_US}	
~/mybin/find_mimg 63 63   ${FEL_TAMD2} > ${mim_FEL_TAMD}
~/mybin/find_mimg 21 21   ${FEL_MS_GBSA2} > ${mim_FEL_MS_GBSA}
~/mybin/find_mimg 60 60   ${FEL_REMD_GBSA2} > ${min_FEL_REMD_GBSA}
~/mybin/find_mimg 31 31   ${FEL_MS_ME2} > ${min_FEL_MS_ME}	  

sort -n -k4  ${mim_FEL_MS} > ${mim_FEL_MS2}
sort -n -k4  ${mim_FEL_REMD} > ${mim_FEL_REMD2}
sort -n -k4  ${mim_FEL_US} > ${mim_FEL_US2}
sort -n -k4  ${mim_FEL_TAMD} > ${mim_FEL_TAMD2}
sort -n -k4  ${mim_FEL_MS_GBSA} > ${mim_FEL_MS_GBSA2}
sort -n -k4  ${min_FEL_REMD_GBSA} > ${min_FEL_REMD_GBSA2}
sort -n -k4  ${min_FEL_MS_ME} > ${min_FEL_MS_ME2}
