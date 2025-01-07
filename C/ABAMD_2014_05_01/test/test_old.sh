#!/bin/csh
#$ -S /bin/sh 
#$ -cwd 
#$ -V
#$ -N test
#$ -q all.q@ajisai01
#$ -e /home/yamamori/work/programs/ABAMD_2014_05_01/test/output_old/test.e
#$ -o /home/yamamori/work/programs/ABAMD_2014_05_01/test/output_old/test.o
 
cd /home/yamamori/work/programs/ABAMD_2014_05_01/test

LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/yamamori/work/programs/ABAMD_2014_05_01/lib
export LD_LIBRARY_PATH
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/yamamori/lib
export LD_LIBRARY_PATH

/home/yamamori/mybin/ABAMD_NH_new_12-02-23 \
    --nve --dt 0.001 --termon --temp 300 --tau 1.0 \
    --nums 100000 --int 10 --intout 10 --intnc 10 \
    --rst output_old/ALADv.rst --rstv output_old/ALADv.rve \
    input/ALADv_min.rst input/ALADv.clt input/ALADv.top \
    output_old/ALADv.out output_old/ALADv.out2_old output_old/ALADv.trj
