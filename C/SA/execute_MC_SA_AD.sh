#!/bin/sh

opt=(dummy betaini betafin numdt numstep interval  )
nopt=${#opt[*]}
if [ $# -le `expr ${nopt} - 2` ]; then
    echo "USAGE ./execute_min.sh" ${opt[*]:1:${nopt}}
    echo $*
    exit
fi

num=1
while [ $num -le `expr ${nopt} - 1` ]; do
    eval ${opt[$num]}=$1
    shift 1
    num=`expr $num + 1`
done

state=( dummy C7eq C7ax )

for i in `seq 1 2 `; do
    dir=/home/yamamori/work/programs/SA
    cd ${dir}
    mkdir log
    mkdir command

    cat << eof2 > command/qsub_MCSA_${state[$i]}.sh
#!/bin/sh
#$ -S /bin/sh 
#$ -cwd 
#$ -V
#$ -N SA
#$ -q all.q
#$ -l hostname=ajisai02
#$ -e ${dir}/log/MCSA.e
#$ -o ${dir}/log/MCSA.o
 
cd ${dir}/command
LD_LIBRARY_PATH=/opt/intel/Compiler/11.0/083/mkl/lib/em64t:/opt/intel/Compiler/11.0/083/lib/intel64:/opt/intel/Compiler/11.0/083/lib/intel64:/opt/intel/impi/3.2.1.009/lib64:/opt/intel/Compiler/11.0/083/mkl/lib/em64t:/opt/intel/Compiler/11.0/083/lib/intel64:/opt/intel/Compiler/11.0/083/lib/intel64:/home/yamamori/mylib
export LD_LIBRARY_PATH
/home/yamamori/mybin/MC_simulated_annealing -a -n ${numstep} -i ${interval} -j ${interval} -k ${interval} -d 0.01 -X ${betaini} -F ${betafin} -O ${numdt} $INPUT/ADv.crd $INPUT/ADv.top ../ADv_out_${state[$i]}.nc  /home/yamamori/work/programs/MC/AD_${state[$i]}.drst
eof2

    cat command/qsub_MCSA_${state[$i]}.sh
    chmod +x command/qsub_MCSA_${state[$i]}.sh
    qsub command/qsub_MCSA_${state[$i]}.sh
done
