#!/bin/sh


opt=(dummy name  )				        
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

flag=`qstat -u yamamori | grep -G ${name} |wc -l` 	  
echo ${flag}
while [ ${flag} -gt 0  ]; do	       	  	  

    for i in `seq 1 2 `; do
	dir=/home/yamamori/work/programs/SA
	cd ${dir}
	mkdir log
	mkdir command

	cat << eof2 > command/qsub_moMCSA_${state[$i]}.sh
#!/bin/sh
#$ -S /bin/sh 
#$ -cwd 
#$ -V
#$ -N moSA
#$ -q all.q@ajisai01
#$ -q all.q@ajisai02
#$ -e ${dir}/log/moMCSA.e
#$ -o ${dir}/log/moMCSA.o
 
cd ${dir}/command
LD_LIBRARY_PATH=/opt/intel/Compiler/11.0/083/mkl/lib/em64t:/opt/intel/Compiler/11.0/083/lib/intel64:/opt/intel/Compiler/11.0/083/lib/intel64:/opt/intel/impi/3.2.1.009/lib64:/opt/intel/Compiler/11.0/083/mkl/lib/em64t:/opt/intel/Compiler/11.0/083/lib/intel64:/opt/intel/Compiler/11.0/083/lib/intel64:/home/yamamori/mylib
export LD_LIBRARY_PATH
/home/yamamori/mybin/MC_get_ene_term ../ADv_out_${state[$i]}.nc ../ene_${state[$i]}
/home/yamamori/mybin/CD_ncg -K -p $INPUT/ADv.top ../ADv_out_${state[$i]}.nc ../ADv_${state[$i]}.dtrj 
/home/yamamori/mybin/MC_get_trj -l 1 /home/yamamori/work/programs/SA/ADv_out_${state[$i]}.nc ADv_${state[$i]}.ini
eof2

	cat command/qsub_moMCSA_${state[$i]}.sh
	chmod +x command/qsub_moMCSA_${state[$i]}.sh
	qsub command/qsub_moMCSA_${state[$i]}.sh

	flag=`qstat -u yamamori | grep -G moSA |wc -l` 	  
	while [ ${flag} -gt 0  ]; do	       	  	  
	    sleep 2			        	  
	    flag=`qstat -u yamamori | grep -G moSA |wc -l` 	  
	done		            
    done

    flag=`qstat -u yamamori | grep -G ${name} |wc -l` 	  
done       				            
