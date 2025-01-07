#!/bin/sh

if [ -z $3 ]; then
    echo "USAGE: ./'$0' o/g/b  programname sourcefilename1 sourcefilename2 ......"
    exit
fi

source ${PROG}/compile/lib.sh
#source ${WORK}/programs_trunk/compile/lib.sh

opt=$1
programname=$2
comsh=${PROG}/compile/compile.sh

if [ ${opt} == "g" -o ${opt} == "b" ]; then 
    echo "#!/bin/sh" > ${HOME}/mybin/compileg_${programname}.sh
    echo "chmod +x" ${comsh} >> ${HOME}/mybin/compileg_${programname}.sh
    echo "cd " $(pwd) >> ${HOME}/mybin/compileg_${programname}.sh
    echo "sh " ${comsh}  $@ >> ${HOME}/mybin/compileg_${programname}.sh
fi
if [ ${opt} == "o" -o ${opt} == "b" ]; then 
    echo "#!/bin/sh" > ${HOME}/mybin/compile_${programname}.sh
    echo "chmod +x" ${comsh} >> ${HOME}/mybin/compile_${programname}.sh
    echo "cd " $(pwd) >> ${HOME}/mybin/compile_${programname}.sh
    echo "sh " ${comsh}  $@ >> ${HOME}/mybin/compile_${programname}.sh
fi

nfile=`expr $# - 1 `
for num in `seq 0 ${nfile}`; do
    sourcefilename[$num]=$3
    shift 1
done

nfile=`expr ${nfile} - 2 `
if [ ${opt} == "g" -o ${opt} == "b" ]; then 
    gcc -o ${programname}g -pg -g ${sourcefilename[@]} -L${HOME}/lib -I${HOME}/include  -L${HOME}/mylib -I${HOME}/myinclude -Xlinker -rpath -Xlinker ${HOME}/mylib -Xlinker -rpath -Xlinker ${HOME}/lib -L${HOME}/lib/glib-2.0 -I${HOME}/include/glib-2.0 -I${HOME}/lib/glib-2.0/include -I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include  -L/lib64  ${libg[@]} ;
    mv ${programname}g ${HOME}/mybin
    cat ${HOME}/mybin/compileg_${programname}.sh
    echo 

    if [ ! -f ${HOME}/mybin/gdbini_${programname} ]; then
	cat << eof > ${HOME}/mybin/gdbcdc_${programname}
cd $(pwd)
pwd
break main
run
eof
	cat << eof > ${HOME}/mybin/gdbsave_${programname}
set history filename ${HOME}/mybin/gdbhist_${programname}
set history save on
quit
eof
    fi
fi

if [ ${opt} == "o" -o ${opt} == "b" ]; then 
    gcc -o ${programname} -O4 ${sourcefilename[@]} -L${HOME}/lib -I${HOME}/include -L${HOME}/mylib -I${HOME}/myinclude -Xlinker -rpath -Xlinker -I${HOME}/mylib -L${HOME}/lib/glib-2.0 -I${HOME}/include/glib-2.0 -I${HOME}/lib/glib-2.0/include -I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include  -L/lib64 ${libo[@]} ;
    mv ${programname} ${HOME}/mybin
    cat ${HOME}/mybin/compile_${programname}.sh
    echo

    if [ ! -f ${HOME}/mybin/execute_${programname}.sh  ]; then
	cat << eof > ${HOME}/mybin/execute_${programname}.sh
#!/bin/sh

opt=(dummy )
nopt='${#opt[*]}'
if [ $# -le '`expr ${nopt} - 2`'  ]; then
    echo "USAGE " '$0'  '${opt[*]:1:${nopt}}'
    echo $*
    exit
fi

num=1
while [ $num -le `expr ${nopt} - 1` ]; do
    eval ${opt[$num]}=$1
    shift 1
    num=`expr $num + 1`
done

dir=$(pwd)
cat << eof2 > qsub_'${programname}'.sh
#!/bin/sh
#PBS -N '${programname}'
#PBS -q serial
#PBS -l nodes=1
#PBS -l ncpus=1
#PBS -e ${dir}/'${programname}'.e
#PBS -o ${dir}/'${programname}'.o
 
cd ${dir}
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/yamamori/mylib
export LD_LIBRARY_PATH
/home/yamamori/mybin/'${programname}' 
eof2

cat qsub_'${programname}'.sh
chmod +x qsub_'${programname}'.sh
qsub qsub_'${programname}'.sh
eof
    fi
    
fi

mkdir -p ~/test/${programname}/log
cat <<EOF >> ~/test/${programname}/log/log.txt
#$(date +"%Y-%m-%d")
EOF

cat <<EOF >> ~/test/${programname}/${programname}_debug_base.sh
#$(date +"%Y-%m-%d")

EOF

