#!/bin/sh

LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/yamamori/work/programs/ABAMD_2014_05_01/lib
export LD_LIBRARY_PATH
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/yamamori/lib
export LD_LIBRARY_PATH

RES=( dummy ALA ARG ASN ASP CYS GLN GLU GLY HIE ILE LEU LYS MET SER TRP THR TYR PHE PRO VAL ASH CYM CYX GLH HIE HIP LYN )
#RES=( dummy HIE )
#RES=( dummy ARG )
#RES=( dummy ALA )
NRES=$(expr ${#RES[*]} - 1)
nvet=( dummy --nve "" )
nnvet=$(expr ${#nvet[*]} - 1)

for i in `seq 1 ${NRES}`; do
    for j in `seq 1 ${nnvet}`; do

	prog=/home/yamamori/work/programs/ABAMD_2014_05_01/src_woGC_wcalloc/ABAMD_woGC_wcalloc

	crd=~/calspa/input/dipep/${RES[$i]}Dv.crd
	top=~/calspa/input/dipep/${RES[$i]}Dv.top
	clt=~/calspa/input/dipep/${RES[$i]}Dv.clt

	rst=output_valgrind_woGC_wcalloc/${RES[$i]}Dv.rst
	rve=output_valgrind_woGC_wcalloc/${RES[$i]}Dv.rve
	out=output_valgrind_woGC_wcalloc/${RES[$i]}Dv.out
	out2=output_valgrind_woGC_wcalloc/${RES[$i]}Dv.out2
	trj=output_valgrind_woGC_wcalloc/${RES[$i]}Dv.trj
	
	log=output_valgrind_woGC_wcalloc/${RES[$i]}Dv.log

	valgrind -v --track-origins=yes ${prog} \
	    ${nvet[$j]} --dt 0.001 --termon --temp 300 --tau 1.0 \
	    --nums 100 --int 10 --intout 10 --intnc 10 \
	    --rst ${rst} --rstv ${rve} \
	    ${crd} ${clt} ${top} \
	    ${out} ${out2} ${trj} > ${log}
	done
done

#	--nums 1 --int 1 --intout 1 --intnc 1 \
#	--nums 100 --int 10 --intout 10 --intnc 10 \
#	--nums 100000 --int 1000 --intout 1000 --intnc 1000 \

#run --nve --dt 0.001 --termon --temp 300 --tau 1.0 --nums 100 --int 10 --intout 10 --intnc 10 --rst output_valgrind_woGC_wcalloc/ARGDv.rst --rstv output_valgrind_woGC_wcalloc/ARGDv.rve /home/yamamori/calspa/input/dipep/ARGDv.crd /home/yamamori/calspa/input/dipep/ARGDv.clt /home/yamamori/calspa/input/dipep/ARGDv.top output_valgrind_woGC_wcalloc/ARGDv.out output_valgrind_woGC_wcalloc/ARGDv.out2 output_valgrind_woGC_wcalloc/ARGDv.trj

#ARG ASN MET TRP TYR PHE PRO HIP

