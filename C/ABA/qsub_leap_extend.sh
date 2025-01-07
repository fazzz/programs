#!/bin/csh
#$ -S /bin/sh 
#$ -cwd 
#$ -V
#$ -N leap
#$ -q all.q
#$ -l hostname=ajisai01
#$ -l hostname=ajisai03
#$ -e /home/yamamori/work/programs/ABA/leap.e
#$ -o /home/yamamori/work/programs/ABA/leap.o
 
cd /home/yamamori/work/programs/ABA
export AMBERHOME=/home/appl/amber10
/home/appl/amber10/exe/tleap -s -f //home/appl/amber10/dat/leap/cmd/leaprc.ff99SB -f /home/yamamori/work/programs/ABA/CHIvv_extend.cmd
