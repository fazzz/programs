#!/bin/csh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -q all.q
#$ -pe mpi 2
#$ -e mpi2.e
#$ -o mpi2.o
#$ -N job
#$ -l hostname=ajisai02

cd ~/work/programs/MPI/
mpirun -np 2 ~/work/programs/MPI/mpi2



