#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mpi.h"

int main(int argc, char *argv[])
{
int namelen, num_procs, myrank;
 char processor_name[MPI_MAX_PROCESSOR_NAME];

 MPI_Init(&argc, &argv);

 MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
 MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
 MPI_Get_processor_name(processor_name, &namelen);

 printf("Hellow, MPI! at Process %d of %d on %s\n", myrank, num_procs, processor_name);

 MPI_Finalize();

 return EXIT_SUCCESS;
}
