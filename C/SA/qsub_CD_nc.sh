#!/bin/sh
#$ -S /bin/sh 
#$ -cwd 
#$ -V
#$ -N CD
#$ -q all.q
#$ -l hostname=ajisai02
#$ -e /home/yamamori/work/programs/SA/CD.e
#$ -o /home/yamamori/work/programs/SA/CD.o
 
cd /home/yamamori/work/programs/SA
LD_LIBRARY_PATH=/opt/intel/Compiler/11.0/083/mkl/lib/em64t:/opt/intel/Compiler/11.0/083/lib/intel64:/opt/intel/Compiler/11.0/083/lib/intel64:/opt/intel/impi/3.2.1.009/lib64:/opt/intel/Compiler/11.0/083/mkl/lib/em64t:/opt/intel/Compiler/11.0/083/lib/intel64:/opt/intel/Compiler/11.0/083/lib/intel64:/home/yamamori/mylib
export LD_LIBRARY_PATH
/home/yamamori/mybin/CD_nc -K -i ADv_trj.nc -t 100 -o ADv.dtrj -p /home/yamamori/calspa/input/ADv.top
