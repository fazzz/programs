#include <stdio.h>
#include <math.h>

#include "const.h"
#include "dPCA.h"

double pick_dihed_one_clust(double coord[4][3]);

void pick_dihed(/*double *dihed_traj ,*/double *traj,int numatom,int numdihed, int atom_dihed_pair[MAXNUMDIHED][4], int time)
{
  int i,j,k;
  double coord[4][3];

  dihed_traj=malloc(sizeof(double)*time*numdihed);
 
  for (i=0;i<time;++i) {
    for (j=0;j<numdihed;++j) {
      for (k=0;k<3;++k) {
	coord[0][k]=traj[i*numatom*3+atom_dihed_pair[j][0]*3+k];
	coord[1][k]=traj[i*numatom*3+atom_dihed_pair[j][1]*3+k];
	coord[2][k]=traj[i*numatom*3+atom_dihed_pair[j][2]*3+k];
	coord[3][k]=traj[i*numatom*3+atom_dihed_pair[j][3]*3+k];
      }
      dihed_traj[i*numdihed+j]/*rad*/ = pick_dihed_one_clust(coord);
    }
  }
}

// 剛体中の全ての 2 面角の計算を行う関数
double pick_dihed_one_clust(double coord[4][3])
{
  int alpha;
  int i;
  int num_ATOM_N_Clut;
  int num_ATOM_N_Clut_1;
  int nNumDihedClut;
  int nNumDihedClut_1;
  int nNumDihedClut_2;
  int nNumAtomOrigin;
  int nNumAtomTarget;
  
  int num_dihed_clust;
  
  double CN_A[3];
  double HN_A[3];
  double CN_1_A[3];
  double HN_1_A[3];
  
  double cc1[3];
  double cc2[3];
  double cc1cc2[3];

  double d_CN_1_A_CN_A;
  double d_HN_A_CN_A;
  double d_HN_1_A_CN_1_A ;

  double cs=0.0;
  double sn=0.0;

  double det;

  double vect[3];

  double sintheta=0.0;
  double sinphi=0.0;

  double theta;


  // CN_A と CN_A_1 の座標の取得
  for(alpha=0;alpha<3;++alpha) {
    CN_A[alpha]/*A*/ = coord[1][alpha]/*A*/;
    CN_1_A[alpha]/*A*/ = coord[2][alpha]/*A*/;
  }

  // HN_A の座標の取得
  for(alpha=0;alpha<3;++alpha){
		HN_A[alpha]/*A*/ = coord[0][alpha]/*A*/;
	}

	d_CN_1_A_CN_A    =0.0;
	d_HN_A_CN_A		 =0.0;
	d_HN_1_A_CN_1_A  =0.0;

	cs=0.0;
	sn=0.0;

 	sintheta = 0.0;
	sinphi = 0.0;

	// HN_A_1 の座標の取得
	for(alpha = 0 ; alpha < 3 ; ++alpha)
	{
		HN_1_A[alpha]/*A*/
		 = coord[3][alpha]/*A*/;
	}

	for(alpha = 0 ; alpha < 3 ; ++alpha)
	{
		d_CN_1_A_CN_A/*A^2*/   +=  (CN_1_A[alpha]-CN_A[alpha])/*A*/
		                          *(CN_1_A[alpha]-CN_A[alpha])/*A*/;
		d_HN_A_CN_A/*A^2*/     +=  (HN_A[alpha]-CN_A[alpha])/*A*/
		                          *(HN_A[alpha]-CN_A[alpha])/*A*/;
		d_HN_1_A_CN_1_A/*A^2*/ +=  (HN_1_A[alpha]-CN_1_A[alpha])/*A*/
		                          *(HN_1_A[alpha]-CN_1_A[alpha])/*A*/;
	}

	d_CN_1_A_CN_A/*A*/ = sqrt(d_CN_1_A_CN_A)/*A*/;
	d_HN_A_CN_A/*A*/ = sqrt(d_HN_A_CN_A)/*A*/;
	d_HN_1_A_CN_1_A/*A*/ = sqrt(d_HN_1_A_CN_1_A)/*A*/;

	for(alpha=0;alpha<3;++alpha)
	{
		sintheta/*rad*/ += (CN_1_A[alpha]-CN_A[alpha])
		                  *(HN_A[alpha]-CN_A[alpha]);
		sinphi/*rad*/   += (CN_A[alpha]-CN_1_A[alpha])
		                  *(HN_1_A[alpha]-CN_1_A[alpha]);
	}

	sintheta/*costheta*//*rad*/ = sintheta/*rad*/
	                              /(d_CN_1_A_CN_A*d_HN_A_CN_A);
	sintheta/*rad*/ = 1.0 - sintheta/*rad*/*sintheta/*rad*/;
	if (sintheta != 0)
		sintheta/*rad*/ = sqrt(sintheta)/*rad*/;
	else 
		sintheta/*rad*/ = 0.0;

	sinphi/*rad*/ = sinphi/(d_CN_1_A_CN_A*d_HN_1_A_CN_1_A);
	sinphi/*rad*/ = 1.0 - sinphi*sinphi;
	if (sinphi != 0)
		sinphi/*rad*/ = sqrt(sinphi)/*rad*/;
	else
		sinphi/*rad*/ = 0.0;


	cc1[0]=(  (HN_A[1]-CN_A[1])
		 *(CN_1_A[2]-CN_A[2])
		   -  
		  (HN_A[2]-CN_A[2])
		 *(CN_1_A[1]-CN_A[1])
	        )
		   /(d_CN_1_A_CN_A*d_HN_A_CN_A*sintheta);
		
	cc1[1]=(  
		   (HN_A[2]-CN_A[2])
		  *(CN_1_A[0]-CN_A[0])
		   -  
		   (HN_A[0]-CN_A[0])
		  *(CN_1_A[2]-CN_A[2])
		)
		   /(d_CN_1_A_CN_A*d_HN_A_CN_A*sintheta);
			
	cc1[2]=(  (HN_A[0]-CN_A[0])*(CN_1_A[1]-CN_A[1])
		   -(HN_A[1]-CN_A[1])*(CN_1_A[0]-CN_A[0])
		 )
		   /(d_CN_1_A_CN_A*d_HN_A_CN_A*sintheta);
			
	cc2[0]=(  (HN_1_A[1]-CN_1_A[1])
		    *(CN_A[2]-CN_1_A[2])
		  -  
		  (HN_1_A[2]-CN_1_A[2])
	            *(CN_A[1]-CN_1_A[1])
		)
	    /(d_CN_1_A_CN_A*d_HN_1_A_CN_1_A*sinphi);
			
	cc2[1]=(   (HN_1_A[2]-CN_1_A[2])
	      *(CN_A[0]-CN_1_A[0])
		   -    
		   (HN_1_A[0]-CN_1_A[0])
		      *(CN_A[2]-CN_1_A[2]) 
		   )
		    /(d_CN_1_A_CN_A*d_HN_1_A_CN_1_A*sinphi);
			
	cc2[2]=(    (HN_1_A[0]-CN_1_A[0])
		       *(CN_A[1]-CN_1_A[1])
		    -   
		    (HN_1_A[1]-CN_1_A[1])
		      * (CN_A[0]-CN_1_A[0]))
	    /(d_CN_1_A_CN_A*d_HN_1_A_CN_1_A*sinphi);

	for(alpha=0;alpha<3;++alpha)
	{
		cs += cc1[alpha]*cc2[alpha];
	}

	if (cs < -1.0 )
	{
		cs = -1.0;
	}
	else if (cs > 1.0 )
	{
		cs = 1.0;
	}

	cc1cc2[0] = cc1[1]*cc2[2]-cc1[2]*cc2[1];
	cc1cc2[1] = cc1[2]*cc2[0]-cc1[0]*cc2[2];
	cc1cc2[2] = cc1[0]*cc2[1]-cc1[1]*cc2[0];

	det = cc1cc2[0]*(CN_1_A[0]-CN_A[0])
		 +cc1cc2[1]*(CN_1_A[1]-CN_A[1])
		 +cc1cc2[2]*(CN_1_A[2]-CN_A[2]);

	if (det > 0.0)
	{
	  theta/*rad*/ = pi/*rad*/ + acos(cs)/*rad*/;
	}
	else
	{
	  theta/*rad*/ = pi/*rad*/ - acos(cs)/*rad*/;
	}

	return theta/*rad*/;

}

