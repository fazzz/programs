#!/bin/sh

LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/yamamori/work/programs/ABAMD_2014_05_01/lib
export LD_LIBRARY_PATH
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/yamamori/lib
export LD_LIBRARY_PATH
#10485760
valgrind -v --track-origins=yes /home/yamamori/work/programs/ABAMD_2014_05_01/src_woGC/ABAMD_woGC \
    --nve --dt 0.001 --termon --temp 300 --tau 1.0 \
    --nums 100 --int 10 --intout 10 --intnc 10 \
    --rst output_valgrind_woGC/ALADv.rst --rstv output_valgrind_woGC/ALADv.rve \
    input/ALADv_min.rst input/ALADv.clt input/ALADv.top \
    output_valgrind_woGC/ALADv.out output_valgrind_woGC/ALADv.out2 output_valgrind_woGC/ALADv.trj > output_valgrind_woGC/log.txt

#    --nums 1 --int 1 --intout 1 --intnc 1 \
#    --nums 100000 --int 10 --intout 10 --intnc 10 \

#==24015== ERROR SUMMARY: 866 errors from 108 contexts (suppressed: 4 from 4)
#==24067== ERROR SUMMARY: 848 errors from 108 contexts (suppressed: 4 from 4)
#==24405== ERROR SUMMARY: 806 errors from 107 contexts (suppressed: 4 from 4)
#==24686== ERROR SUMMARY: 805 errors from 106 contexts (suppressed: 4 from 4)
#==25174== ERROR SUMMARY: 763 errors from 105 contexts (suppressed: 4 from 4)
#==25328== ERROR SUMMARY: 745 errors from 104 contexts (suppressed: 4 from 4)
#==27190== ERROR SUMMARY: 181 errors from 75 contexts (suppressed: 4 from 4)
#==5236== ERROR SUMMARY: 127 errors from 21 contexts (suppressed: 4 from 4)
#==20423== ERROR SUMMARY: 139 errors from 25 contexts (suppressed: 4 from 4)
#==21088== ERROR SUMMARY: 132 errors from 24 contexts (suppressed: 4 from 4)
#==21998== ERROR SUMMARY: 120 errors from 20 contexts (suppressed: 4 from 4)
#==21911== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 4 from 4)
