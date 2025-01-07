#!/bin/sh

opt=( dummy match ca tmpfiles )
nopt=${#opt[*]}
if [ $# -le `expr ${nopt} - 2` ]; then
    echo "USAGE: $0" ${opt[*]:1:${nopt}}
    echo $*
    exit
fi

num=1
while [ $num -le `expr ${nopt} - 1` ]; do
    eval ${opt[$num]}=$1
    shift 1
    num=`expr $num + 1`
done

cat <<EOF > ${tmpfiles}/conv.tmp.awk
BEGIN{

EOF

cat ${ca} | while read line; do
    name_atm=$( echo ${line} | awk '{print $1}' )
    num_atm=$( echo ${line} | awk '{print $2}' )
    num_res=$( echo ${line} | awk '{print $3}' )
    
    echo "MA[\"${num_res}\"]=\"${num_atm}\"; " >> ${tmpfiles}/conv.tmp.awk
done

cat <<EOF >> ${tmpfiles}/conv.tmp.awk
}

EOF

cat <<EOF >> ${tmpfiles}/conv.tmp.awk
{ printf("%4d CA %4d\n", \$1, MA[\$1]) }
EOF

awk -f ${tmpfiles}/conv.tmp.awk ${match}












