#!/bin/sh

files=( dummy FF PT )

Nfiles=$(expr ${#files[*]} - 1)

for i in `seq 1 ${Nfiles}` ; do
    ls ~/work/programs/*/${files[$i]}.h
    same_name_files=$( ls ~/work/programs/*/${files[$i]}.h )
    Nsame_name_files=$(expr ${#same_name_files[*]} - 1)

    for j in `seq 1 ${Nsame_name_files}` ; do
	diff -c ${same_name_files[0]} ${same_name_files[$j]}
    done
done
