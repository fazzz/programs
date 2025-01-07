#!/bin/sh

opt=( dummy prot_name pdbid )
nopt=${#opt[*]}
if [ $# -le `expr ${nopt} - 2` ]; then
    echo "USAGE $0" ${opt[*]:1:${nopt}}
    echo $*
    exit
fi

num=1
while [ $num -le `expr ${nopt} - 1` ]; do
    eval ${opt[$num]}=$1
    shift 1
    num=`expr $num + 1`
done

echo "prot_name="${prot_name}

echo "pdbid="${pdbid}

mkdir ~/work/data/pdb/pdb_${pdbid}

cd ~/work/data/pdb/pdb_${pdbid}

if [ -f ${pdbid}.pdb.gz ]; then
    rm ${pdbid}.pdb.gz
fi

wget http://rcsb.org/pdb/files/${pdbid}.pdb.gz

gunzip ${pdbid}.pdb.gz

gawk '$1 ~ /ATOM/ { print $0} ' ${pdbid}.pdb > ${prot_name}.pdb

echo ${pdbid} ${prot_name} > ~/work/data/pdb/pdbid_proteinname.txt
