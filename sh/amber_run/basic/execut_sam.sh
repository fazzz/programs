#!/bin/sh
#./execute_sam.sh 4 4 2 10 20

if [ -z "$5" ]
then
   echo :Usage ./$( basename $0 ) ndt ntu ic initial max_run
   exit
fi

ndt=$1
ntu=$2
ic=$3

i=1
while [ $i -le $ndt ] 
do
j=1
while [ $j -le $ntu ] 
do
k=1
while [ $k -le $ic ] 
do
cd ../dt${i}tu${j}st${k}

set dummy 0.002 0.004 0.006 0.008
shift ${i}
dt=$1
echo ${dt}
set dummy 500000 250000 1670000 125000
#set dummy 50 25 16 12
shift ${j}
step=$1
set dummy 500 250 167 125
out=$1
#out=1
set dummy 500.0 1000.0 1500.0 2000.0
shift ${k} 
tu=$1
tu=$(expr ${tu}*${dt})
echo $tu
cat << eof1 > md.in
1 0 300.0 0.0 $dt $step 1 20 $out $out 0 0 0.0 $tu 1.0
eof1

cat << eof2 > qsub_equ${i}_${j}_${k}.sh
#!/bin/sh
#$ -N TAMDsam${i}_${j}_${k}
#$ -e /home/yamamori/calspa/AAN1002_4/dt${i}tu${j}st${k}/TAMD.e
#$ -o /home/yamamori/calspa/AAN1002_4/dt${i}tu${j}st${k}/TAMD.o
 
cd /home/yamamori/calspa/AAN1002_4/dt${i}tu${j}st${k}
./mass_MD.exe -c crd_rst.in
eof2

remain=$(qstat | grep -G yamamori | wc -l) 
let remain="expr ${max_run} - ${remain}"
while [ -z $remain ] ; do
sleep 10
remain=$(qstat | grep -G yamamori | wc -l) 
let remain="expr ${max_run} - ${remain}"
done

chmod +x qsub_equ${i}_${j}_${k}.sh
rm *.out *.pdb
#qsub qsub_equ${i}_${j}_${k}.sh
#./mass_MD.exe

let k="$k+1"

done
let j="$j+1"
done
let i="$i+1"
done

cd ../comand
