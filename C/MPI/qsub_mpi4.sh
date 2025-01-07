#!/bin/csh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -q all.q
#$ -pe mpi 4
#$ -e mpi4.e
#$ -o mpi4.o
#$ -N job
#$ -l hostname=ajisai02

cd ~/work/programs/MPI/
mpirun -np 4 ~/work/programs/MPI/mpi4



