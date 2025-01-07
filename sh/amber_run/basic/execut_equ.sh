#!/bin/sh
#./execute_sim.sh 0.001 10000 300.0 100,000 100

if [ -z "$5" ]
then
   echo :Usage ./$( basename $0 ) dt tu temp step out
   exit
fi

dt=$1
tu=$2
temp=$3
step=$4
out=$5

cat << eof1 > md.in
1 1 $temp 100.0 $dt $step 1 20 $out $out 0 0 0.0 $tu 1.0
eof1

dir=$(pwd)
cat md.in

cat << eof2 > qsub_sim.sh
#!/bin/csh -f
#$ -N tamd
#$ -q all.q
#$ -V
#$ -e ${dir}/TAMD.e
#$ -o ${dir}/TAMD.o
 
cd ${dir}
pwd
./mass_MD.exe
eof2

chmod +x qsub_sim.sh
rm *.o *.e *.out *.pdb velo.in crd_rst.in
qsub qsub_sim.sh

