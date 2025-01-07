#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "MD.h"
#include "force.h"
#include "BD.h"
#include "EF.h"

//#include "GOLM.h"

// 現在の構造でのタンパク質のポテンシャルエネルギー、
// タンパク質に及ぼす力の計算を行う関数
void calc_force2_CG(int pflag,
		    double *ele, double *ALJ, double *BLJ,
		    double *p_e, double *p_1_4_e, 
		    double *p_LJ,double *p_1_4_LJ,
		    double *f_e, double *f_1_4_e, 
		    double *f_LJ,double *f_1_4_LJ,
		    int numnb, int *indexnb,
		    int num14, int *index14,double *eigen,double *eig_14,double *eig_dihed) {
  int i,ii,i_c,alpha,alpha2,alpha3,j;

  int num;
  int nNumClut;
  
  int nNumOrigc;
  int nNumAbsoc;
  
  int origin_of_this_branch = 0;
  
  double **q_c/*[MAXNAC][3]*/;
  
  double qq[3], qq2[3];
  
  double **f_clust/*[MAXDOF][3]*/;
  double **N_clust/*[MAXDOF][3]*/;
  double ***f_clust3/*[MAXDOF][20][3]*/;
  double ***f_clust4/*[MAXDOF][20][3]*/;
  double f_i_c[3];
  
  double RotatnNumtonNumMiOn[3][3];
  
  double *cord;
  
  FILE *out;
  
  int num_a_prot;
  int tnum_atom_clust,tDOF=prot.DOF;
  
  q_c=(double **)gcemalloc(sizeof(double *)*prot.num_atom);
  for (i=0;i<prot.num_atom;++i)
    q_c[i]=(double *)gcemalloc(sizeof(double)*3);
  f_clust=(double **)gcemalloc(sizeof(double *)*prot.DOF);
  N_clust=(double **)gcemalloc(sizeof(double *)*prot.DOF);
  for (i=0;i<prot.DOF;++i) {
    f_clust[i]=(double *)gcemalloc(sizeof(double)*3);
    N_clust[i]=(double *)gcemalloc(sizeof(double)*3);
  }
  f_clust3=(double **)gcemalloc(sizeof(double *)*prot.DOF);
  f_clust4=(double **)gcemalloc(sizeof(double *)*prot.DOF);
  for (i=0;i<prot.DOF;++i) {
    f_clust3[i]=(double **)gcemalloc(sizeof(double *)*20);
    f_clust4[i]=(double **)gcemalloc(sizeof(double *)*20);
    for (j=0;j<20;++j) {
      f_clust3[i][j]=(double *)gcemalloc(sizeof(double)*20);
      f_clust4[i][j]=(double *)gcemalloc(sizeof(double)*20);
    }
  }
  
  
  // spatial_force の初期化を行う
  for(nNumClut = 0; nNumClut < prot.DOF; ++nNumClut)
    {
      for(alpha=0;alpha<3;++alpha)
	{
	  clust[nNumClut].f_c.sp_f_clust[0].f_clust[alpha] = 0.0;
	  clust[nNumClut].f_c.sp_f_clust[0].N_clust[alpha] = 0.0;
	  f_clust[nNumClut][alpha] = 0.0;
	  N_clust[nNumClut][alpha] = 0.0;
	}
	}
  
  for(nNumClut = 0;nNumClut < prot.DOF; ++nNumClut)
    {
      for(i=0;i < clust[nNumClut].num_atom_clust;++i)
	{
	  for(alpha=0;alpha < 3; ++alpha)
	    {
	      clust[nNumClut].f_c.f_elesta[i][alpha] = 0.0;
	      clust[nNumClut].f_c.f_L_J[i][alpha] = 0.0;
	      clust[nNumClut].f_c.f_1_4_elesta[i][alpha] = 0.0;
	      clust[nNumClut].f_c.f_1_4_L_J[i][alpha] = 0.0;
	      f_clust3[nNumClut][i][alpha] = 0.0;
	      f_clust4[nNumClut][i][alpha] = 0.0;
	    }
	}
    }
  
  if(nbstopflag!=ON) {
    ///////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////
    // VDW 相互作用の計算

    if ((cord=(double *)malloc(sizeof(double)*prot.num_atom*3))==NULL) {
      printf("error of allocation cord\n");
    }
    for (i=0;i<prot.num_atom;++i) {
      for (alpha=0;alpha<3;++alpha) {
	cord[i*3+alpha]=prot.coord[i][alpha];
      }
    }
    GOLMff_calcff(cord,prot.num_atom,&ene,1,1,1,1,1);
    free(cord);
      
  }
  
  for(nNumClut=0;nNumClut<prot.DOF;++nNumClut){
    for(i_c=0;i_c<clust[nNumClut].num_atom_clust;++i_c){
      for(alpha=0;alpha < 3;++alpha){
	f_clust[nNumClut][alpha]
	  +=
	  clust[nNumClut].f_c.f_L_J[i_c][alpha]
	  + clust[nNumClut].f_c.f_elesta[i_c][alpha]
	  + clust[nNumClut].f_c.f_1_4_L_J[i_c][alpha]
	  + clust[nNumClut].f_c.f_1_4_elesta[i_c][alpha]
	  ;
	
	f_clust3[nNumClut][i_c][alpha]
	  =
	  clust[nNumClut].f_c.f_L_J[i_c][alpha]
	  + clust[nNumClut].f_c.f_elesta[i_c][alpha]
	  + clust[nNumClut].f_c.f_1_4_L_J[i_c][alpha]
	  + clust[nNumClut].f_c.f_1_4_elesta[i_c][alpha];
      }
    }
  }

  for(nNumClut=0;nNumClut<prot.DOF;++nNumClut)
    {
      for (alpha=0;alpha<3;++alpha)
	{
	  for (alpha2=0;alpha2<3;++alpha2)
	    {
	      // Rot_0_to_N * Rot_n-1_to_0
	      clust[nNumClut].f_c.sp_f_clust[0].f_clust[alpha]
		+= 	clust[nNumClut].trans_A_to_CN[0][alpha][alpha2]
		*f_clust[nNumClut][alpha2];
	    }
	}
    }

  for(nNumClut=0;nNumClut<prot.DOF;++nNumClut)
    {
      for(i_c=0;i_c<clust[nNumClut].num_atom_clust;++i_c)
	{
	  for (alpha=0;alpha<3;++alpha)
	    {
	      for (alpha2=0;alpha2<3;++alpha2)
		{
		  // Rot_0_to_N * Rot_n-1_to_0
		  f_clust4[nNumClut][i_c][alpha]
		    += 	clust[nNumClut].trans_A_to_CN[0][alpha][alpha2]
		    *f_clust3[nNumClut][i_c][alpha2];
		}
	    }
	}
    }

  for(nNumClut=0;nNumClut<prot.DOF;++nNumClut)
    {
      num = 0/*0410*//*clust[nNumClut].origin_xoord_a*/;
      for(i_c=0;i_c<clust[nNumClut].num_atom_clust;++i_c)
	{
	  clust[nNumClut].f_c.sp_f_clust[0].N_clust[0]
	    += -clust[nNumClut].xoord_clust/*[0]*/[i_c+num][2]*f_clust4[nNumClut][i_c][1]
	    +clust[nNumClut].xoord_clust/*[0]*/[i_c+num][1]*f_clust4[nNumClut][i_c][2];
	  
	  clust[nNumClut].f_c.sp_f_clust[0].N_clust[1] 
	    += clust[nNumClut].xoord_clust/*[0]*/[i_c+num][2]*f_clust4[nNumClut][i_c][0]
	    -clust[nNumClut].xoord_clust/*[0]*/[i_c+num][0]*f_clust4[nNumClut][i_c][2];
	  
	  clust[nNumClut].f_c.sp_f_clust[0].N_clust[2]
	    += -clust[nNumClut].xoord_clust/*[0]*/[i_c+num][1]*f_clust4[nNumClut][i_c][0]
	    +clust[nNumClut].xoord_clust/*[0]*/[i_c+num][0]*f_clust4[nNumClut][i_c][1];
	  
	}
    }
}
