#!/bin/sh
#$ -S /bin/sh 
#$ -cwd 
#$ -V
#$ -N moSA
#$ -q all.q@ajisai01
#$ -q all.q@ajisai02
#$ -e /home/yamamori/work/programs/SA/log/moMCSA.e
#$ -o /home/yamamori/work/programs/SA/log/moMCSA.o
 
cd /home/yamamori/work/programs/SA/command
LD_LIBRARY_PATH=/opt/intel/Compiler/11.0/083/mkl/lib/em64t:/opt/intel/Compiler/11.0/083/lib/intel64:/opt/intel/Compiler/11.0/083/lib/intel64:/opt/intel/impi/3.2.1.009/lib64:/opt/intel/Compiler/11.0/083/mkl/lib/em64t:/opt/intel/Compiler/11.0/083/lib/intel64:/opt/intel/Compiler/11.0/083/lib/intel64:/home/yamamori/mylib
export LD_LIBRARY_PATH
/home/yamamori/mybin/MC_get_ene_term ../ADv_out_C7ax.nc ../ene_C7ax
/home/yamamori/mybin/CD_ncg -K -p /home/yamamori/calspa/input/ADv.top ../ADv_out_C7ax.nc ../ADv_C7ax.dtrj 
