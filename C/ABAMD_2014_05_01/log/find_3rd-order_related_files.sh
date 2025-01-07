#!/bin/sh

files=( ~/work/programs/MolBas/MB.h ~/work/programs/MolBas/MB.c
~/work/programs/MYMATH/mymath.h ~/work/programs/MYMATH/mymath.c
~/work/programs/SBFF/SBFF.h ~/work/programs/SBFF/SBFF.c )

grep -nH -e include ${files[*]}
