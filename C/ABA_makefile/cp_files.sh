#!/bin/sh

cp $PROG/ABA/*.[ch] $PROG/LA/*.[ch]  $PROG/readParmtop/*.[ch] $PROG/netcdf/*.[ch] $PROG/TOPO/*.[ch] $PROG/efunc/*.[ch] $PROG/MYMATH/*.[ch] $PROG/CFF/*.[ch] $PROG/GOLM/*.[ch] $PROG/CD/*.[ch]  $PROG/ABA_makefile 

mv readParmtopL.c PTL.c
mv calcFFL.c FFL.c
mv efunc.c EF.c



