#!/bin/sh

opt=( dummy pdbref pdbtgt pdbfit )
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

if [ ! -d tempfiles ]; then
    mkdir tempfiles
fi

pdbref=$( basename ${pdbref} .pdb)
pdbtgt=$( basename ${pdbtgt} .pdb)

awk '$1~/ATOM/{print $0}' ${pdbref}.pdb > tempfiles/${pdbref}_ATOM.pdb
awk '$1~/ATOM/{print $0}' ${pdbtgt}.pdb > tempfiles/${pdbtgt}_ATOM.pdb

TMalign/TMalign tempfiles/${pdbref}_ATOM.pdb tempfiles/${pdbtgt}_ATOM.pdb > tempfiles/TMalign.txt

tail -4 tempfiles/TMalign.txt | head -1 | fold -w 1 > tempfiles/temp01_TMalign
tail -3 tempfiles/TMalign.txt | head -1 | fold -w 1 > tempfiles/temp02_TMalign
tail -2 tempfiles/TMalign.txt | head -1 | fold -w 1 > tempfiles/temp03_TMalign

paste -d "," tempfiles/temp02_TMalign tempfiles/temp01_TMalign tempfiles/temp03_TMalign > tempfiles/TMalign_inv.txt

awk -F"," -f bestfit/script/parse_for_fit_01.awk tempfiles/TMalign_inv.txt | awk -f bestfit/script/parse_for_fit_03.awk > tempfiles/MatchResList_pdbref.txt

awk -F"," -f bestfit/script/parse_for_fit_02.awk tempfiles/TMalign_inv.txt | awk -f bestfit/script/parse_for_fit_03.awk > tempfiles/MatchResList_pdbtgt.txt

awk -f bestfit/script/parse_for_fit_04.awk tempfiles/${pdbref}_ATOM.pdb > tempfiles/CAList_pdbref.txt

awk -f bestfit/script/parse_for_fit_04.awk tempfiles/${pdbtgt}_ATOM.pdb > tempfiles/CAList_pdbtgt.txt

bestfit/script/CAlist_Matchlist.sh tempfiles/MatchResList_pdbref.txt tempfiles/CAList_pdbref.txt tempfiles | awk -f bestfit/script/parse_for_fit_05.awk > tempfiles/pdbref.ndx
bestfit/script/CAlist_Matchlist.sh tempfiles/MatchResList_pdbtgt.txt tempfiles/CAList_pdbtgt.txt tempfiles | awk -f bestfit/script/parse_for_fit_05.awk > tempfiles/pdbtgt.ndx

bestfit/bin/bestfit tempfiles/${pdbref}_ATOM.pdb tempfiles/${pdbtgt}_ATOM.pdb tempfiles/pdbref.ndx tempfiles/pdbtgt.ndx > ${pdbfit}








