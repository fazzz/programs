#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mpi.h"

int main(int argc, char *argv[])
{
int num_procs, myrank;
 double a[128], b[128];

 MPI_Init(&argc, &argv);
 MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
 MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
 a[myrank] = num_procs - myrank;
 b[myrank] = myrank;

 printf("%e + %e = %e\n", a[myrank], b[myrank], a[myrank] + b[myrank]);

 MPI_Finalize();

 return EXIT_SUCCESS;
}
