#!/bin/csh -f
#$ -N tamd
#$ -q all.q
#$ -V
#$ -e /home/yamamori/calspa/AAN1003/TAMD.e
#$ -o /home/yamamori/calspa/AAN1003/TAMD.o
 
cd /home/yamamori/calspa/AAN1003
pwd
./mass_MD.exe
