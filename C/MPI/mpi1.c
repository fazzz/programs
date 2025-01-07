#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mpi.h"

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);

 printf("Hellow, MPI!\n");

 MPI_Finalize();

 return EXIT_SUCCESS;
 }
