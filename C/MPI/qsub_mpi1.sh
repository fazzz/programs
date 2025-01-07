#!/bin/csh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -q all.q
#$ -pe mpi 2
#$ -e job.eo
#$ -o job.eo
#$ -N job
#$ -l hostname=ajisai02

cd ~/work/programs/MPI/
mpirun -np 2 ~/work/programs/MPI/mpi1



