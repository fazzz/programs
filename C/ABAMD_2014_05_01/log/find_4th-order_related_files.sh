#!/bin/sh

files=( ~/work/programs/CFF/FF.h ~/work/programs/CFF/calcFF.c \
~/work/programs/readParmtop/PT.h ~/work/programs/readParmtop/readParmtop.c )

grep -nH -e include ${files[*]}
