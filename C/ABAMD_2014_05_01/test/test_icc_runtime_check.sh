#!/bin/sh

LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/yamamori/work/programs/ABAMD_2014_05_01/lib
export LD_LIBRARY_PATH
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/yamamori/lib
export LD_LIBRARY_PATH

/home/yamamori/work/programs/ABAMD_2014_05_01/src/ABAMD \
    --nve --dt 0.001 --termon --temp 300 --tau 1.0 \
    --nums 100000 --int 10 --intout 10 --intnc 10 \
    --rst output_icc_runtime_check/ALADv.rst --rstv output_icc_runtime_check/ALADv.rve \
    input/ALADv_min.rst input/ALADv.clt input/ALADv.top \
    output_icc_runtime_check/ALADv.out output_icc_runtime_check/ALADv.out2 output_icc_runtime_check/ALADv.trj > output_icc_runtime_check/log.txt
