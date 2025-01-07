#!/bin/sh
#$ -S /bin/sh 
#$ -cwd 
#$ -V
#$ -N SA
#$ -q all.q
#$ -l hostname=ajisai02
#$ -e /home/yamamori/work/programs/SA/log/MCSA.e
#$ -o /home/yamamori/work/programs/SA/log/MCSA.o
 
cd /home/yamamori/work/programs/SA/command
LD_LIBRARY_PATH=/opt/intel/Compiler/11.0/083/mkl/lib/em64t:/opt/intel/Compiler/11.0/083/lib/intel64:/opt/intel/Compiler/11.0/083/lib/intel64:/opt/intel/impi/3.2.1.009/lib64:/opt/intel/Compiler/11.0/083/mkl/lib/em64t:/opt/intel/Compiler/11.0/083/lib/intel64:/opt/intel/Compiler/11.0/083/lib/intel64:/home/yamamori/mylib
export LD_LIBRARY_PATH
/home/yamamori/mybin/MC_simulated_annealing -a -n 10000000 -i 10000 -j 10000 -k 10000 -d 0.01 -X 1000 -F 10 -O 100 /home/yamamori/calspa/input/ADv.crd /home/yamamori/calspa/input/ADv.top ../ADv_out_C7ax.nc  /home/yamamori/work/programs/MC/AD_C7ax.drst
