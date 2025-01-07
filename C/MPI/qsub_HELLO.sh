#!/bin/csh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -q all.q@ajisai01
#$ -q all.q@ajisai02
#$ -pe mpi 2
#$ -e ~/work/programs/MPI/HELLO.e
#$ -o ~/work/programs/MPI/HELLO.o
#$ -N job

cd ~/work/programs/MPI/
mpirun -np 2 ~/work/programs/MPI/hello
