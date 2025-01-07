#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mpi.h"

int main(int argc, char *argv[]) {
  int num_procs, myrank;
  int tag=0;
  MPI_Status status;
  double a=0.0, b=0.0;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  a=myrank;

  if(myrank%2==0) {
    if (num_procs>myrank) 
      MPI_Send((void *)&a, 1, MPI_DOUBLE, myrank+1, tag, MPI_COMM_WORLD);
  }
  else if(myrank%2==1) {
    MPI_Recv((void *)&b, 1, MPI_DOUBLE, myrank-1, tag, MPI_COMM_WORLD, &status);
  }

 printf("Process:%d %e + %e = %e\n", myrank,a, b, a + b);

 MPI_Finalize();

 return EXIT_SUCCESS;
}
