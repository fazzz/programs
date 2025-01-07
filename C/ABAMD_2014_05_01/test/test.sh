#!/bin/csh
#$ -S /bin/sh 
#$ -cwd 
#$ -V
#$ -N test
#$ -q all.q@ajisai01
#$ -e /home/yamamori/work/programs/ABAMD_2014_05_01/test/output/test.e
#$ -o /home/yamamori/work/programs/ABAMD_2014_05_01/test/output/test.o
 
cd /home/yamamori/work/programs/ABAMD_2014_05_01/test

LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/yamamori/work/programs/ABAMD_2014_05_01/lib
export LD_LIBRARY_PATH
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/yamamori/lib
export LD_LIBRARY_PATH

/home/yamamori/work/programs/ABAMD_2014_05_01/src/ABAMD \
    --nve --dt 0.001 --termon --temp 300 --tau 1.0 \
    --nums 100000 --int 10 --intout 10 --intnc 10 \
    --rst output/ALADv.rst --rstv output/ALADv.rve \
    input/ALADv_min.rst input/ALADv.clt input/ALADv.top \
    output/ALADv.out output/ALADv.out2 output/ALADv.trj
