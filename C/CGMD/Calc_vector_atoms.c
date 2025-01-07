#include <stdio.h>
#include <math.h>
#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "MD.h"
#include "force.h"
#include "physics.h"

// 原子間のベクトルの計算を行う関数
void Calc_vector_atoms(void)
{
	int i,j,alpha;
	FILE *out;
	double d_i_j;

	for(i=0;i<prot.num_atom;++i){
	  for(j=/*0*/i+1;j<prot.num_atom;++j){
	    d_i_j = 0.0;
	    for(alpha=0;alpha<3;++alpha){
		/**********************/
                /* if (i != j)	      */
		/*   {		      */
                /**********************/
	      //					q[i][j][alpha]/*A*/ = prot.coord[i][alpha]/*A*/
	      //					                    - prot.coord[j][alpha]/*A*/;
	      q[i][j][alpha]/*A*/ = prot.coord[j][alpha]/*A*/- prot.coord[i][alpha]/*A*/;
	      q[j][i][alpha] -= q[i][j][alpha];
	      d_i_j/*A^2*/ += q[i][j][alpha]/*A*/*q[i][j][alpha]/*A*/;
	    }
	    /************/
	    /* }	    */
	    /************/
	    /**********************/
            /* if (i != j)	  */
	    /*   {		  */
            /**********************/
		d_i_j/*A*/ = sqrt(d_i_j/*A^2*/);
		len_q[i][j]/*A*/ = d_i_j/*A*/;
		len_q[i][j] = d_i_j;
	    /***********************************************/
            /*   }					   */
	    /* else					   */
	    /*   {					   */
	    /* 	len_q[i][j]/\*A*\/ = 0.0/\*A*\/;	   */
	    /*   }					   */
            /***********************************************/
	  }
	}

//	out = fopen("len_q.out","w");
//	for(i=0;i<prot.num_atom;++i)
//	{
//		fprintf(out, "%d \n", i);
//		for(j=0;j<prot.num_atom;++j)
//		{
//			fprintf(out, "%d ", j);
//			fprintf(out, "%lf \n", len_q[i][j]);
//		}
//		fprintf(out, "%d \n", i);
//	}
//	fclose(out);
}

int which_calc_non_bonding_pote(int nNumClut, int i_c, int j_a) {
  int kkk;
  
  int flag = 1;

  for(kkk=0;;++kkk){
    if (clust[nNumClut].pairs.not_interacting[i_c][kkk] == 0){
      break;
    }
    else if (clust[nNumClut].pairs.not_interacting[i_c][kkk] == j_a+1){
      flag = 0;
      break;
    }
    else{
      ;
    }
  }

  if (flag == 0){
    return 0;
  }
  else{
    return 1;
  }
}

int which_calc_1_4_non_bonding_pote(int nNumClut, int i_c, int j_a) {
  int kkk;

  int flag = 1;
  
  for(kkk=0;;++kkk){
    if (clust[nNumClut].o_f_pairs.o_f_not_interacting[i_c][kkk] == 0) {
      break;
    }
    else if (clust[nNumClut].o_f_pairs.o_f_not_interacting[i_c][kkk]== j_a+1){
      flag = 0;
      break;
    }
    else {
      ;
    }
  }

  if (flag == 0) {
    return 1;
  }
  else {
    return 0;
  }
}

double power(double num,int n)
{
	int i;

	double pop=1.0;

	for(i = 0; i < n; ++i)
	{
		pop = pop*num;
	}

	return pop;
}
