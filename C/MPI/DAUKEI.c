
#include <stdio.h>
#include "mpi.h"
#include "daikei.h"

float function (float);
int assign (int, int, float, float, int, float *, float *, int *);

#define X_MIN 0.0		
#define X_MAX 100.0		
#define DIV_NUM 120		

int
main (int argc, char **argv)
{
    int my_rank;		
    int p_num;			
    float x_min = X_MIN;	
    float x_max = X_MAX;	
    int div_num = DIV_NUM;	
    float result;		
    float h;
    float local_x_min;		
    float local_x_max;		
    int local_n;		
    float local_result;		
    int source;
    int dest = 0;
    int tag = 0;
    MPI_Status status;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size (MPI_COMM_WORLD, &p_num);

    assign (p_num, my_rank, x_min, x_max, div_num, &local_x_min, &local_x_max,
	    &local_n);
    local_result =
	daikei_integral (function, local_x_min, local_x_max, local_n);


    if (my_rank == 0)		
      {
	  fprintf (stdout, "P%d :: x = [%f,%f] : n = %d : s = %f\n", 0,
		   local_x_min, local_x_max, local_n, local_result);
	  result = local_result;
	  for (source = 1; source < p_num; source++)
	    {
		MPI_Recv (&local_result, 1, MPI_FLOAT, source, tag,
			  MPI_COMM_WORLD, &status);

		assign (p_num, source, x_min, x_max, div_num, &local_x_min,
			&local_x_max, &local_n);
		fprintf (stdout, "P%d :: x = [%f,%f] : n = %d : s = %f\n",
			 source, local_x_min, local_x_max, local_n,
			 local_result);

		result += local_result;
	    }
      }
    else			
      {
	  MPI_Send (&local_result, 1, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
      }

    if (my_rank == 0)
      {
	  fprintf (stdout, "TOTAL :: x = [%f,%f] : n = %d : S = %f\n", X_MIN,
		   X_MAX, DIV_NUM, result);
      }

    MPI_Finalize ();
}

float
function (float x)
{
    float result;

    result = x;

    return result;
}

int
assign (int p_num, int rank, float x_min, float x_max, int div_num,
	float *local_x_min, float *local_x_max, int *local_n)
{
    float h;
    *local_n = div_num / p_num;
    h = (x_max - x_min) / div_num;
    *local_x_min = x_min + rank * (*local_n) * h;
    *local_x_max = (*local_x_min) + (*local_n) * h;
    return 0;
}
