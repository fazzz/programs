#!/bin/csh
#$ -S /bin/sh 
#$ -cwd 
#$ -V
#$ -N leapACE
#$ -q all.q
#$ -l hostname=ajisai01
#$ -e /home/yamamori/work/programs/yuMD_ECEPP/src/pre_MD/src/leap.e
#$ -o /home/yamamori/work/programs/yuMD_ECEPP/src/pre_MD/src/leap.o
 
cd /home/yamamori/work/programs/yuMD_ECEPP/src/pre_MD/src
export AMBERHOME=/home/appl/amber10
/home/appl/amber10/exe/tleap -s -f //home/appl/amber10/dat/leap/cmd/leaprc.ff99SB -f /home/yamamori/work/programs/yuMD_ECEPP/src/pre_MD/src/A10v.cmd
