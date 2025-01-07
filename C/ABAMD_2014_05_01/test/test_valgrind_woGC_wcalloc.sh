#!/bin/sh

LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/yamamori/work/programs/ABAMD_2014_05_01/lib
export LD_LIBRARY_PATH
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/yamamori/lib
export LD_LIBRARY_PATH
#10485760
valgrind -v --track-origins=yes /home/yamamori/work/programs/ABAMD_2014_05_01/src_woGC_wcalloc/ABAMD_woGC_wcalloc \
    --nve --dt 0.001 --termon --temp 300 --tau 1.0 \
    --nums 100 --int 10 --intout 10 --intnc 10 \
    --rst output_valgrind_woGC_wcalloc/ALADv.rst --rstv output_valgrind_woGC_wcalloc/ALADv.rve \
    input/ALADv_min.rst input/ALADv.clt input/ALADv.top \
    output_valgrind_woGC_wcalloc/ALADv.out output_valgrind_woGC_wcalloc/ALADv.out2 output_valgrind_woGC_wcalloc/ALADv.trj > output_valgrind_woGC_wcalloc/log.txt
