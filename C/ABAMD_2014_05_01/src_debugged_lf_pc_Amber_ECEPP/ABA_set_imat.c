
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#include "ABA.h" // 2014-06-18
#include "ABAb.h"  // 2014-06-18
//#include "EF.h"
#include "EF.h" // 2014-06-17

void ABAs_inertia_matrix(CLT *clt,int nNumClut_all,int num_atom_all,double *crd,double *mass) {
  int i;
  int nNumClut;
  int natomtotal=0;
  int joinflag=0;

  for(nNumClut=0; nNumClut<nNumClut_all; ++nNumClut) {
    //    clt[nNumClut].mass=(double *)gcemalloc(sizeof(double)*clt[nNumClut].num_atom_clust); // 2014-07-22
    //    clt[nNumClut].mass=(double *)emalloc(sizeof(double)*clt[nNumClut].num_atom_clust); // 2014-07-22 // 2014-09-05
    clt[nNumClut].mass=(double *)calloc(clt[nNumClut].num_atom_clust,sizeof(double)); // 2014-07-22 // 2014-09-05
    for(i=0; i<clt[nNumClut].num_atom_clust; ++i) {
      clt[nNumClut].mass[i]=mass[natomtotal+i];
    }
    natomtotal+=clt[nNumClut].num_atom_clust;
    ABAs_Inertia_clust(clt[nNumClut].I,
		       clt[nNumClut].num_atom_clust,
		       clt[nNumClut].xoord,clt[nNumClut].mass);
    ABAs_InertiaMatrix(clt[nNumClut].IM,
		       clt[nNumClut].I,
		       clt[nNumClut].qCOM,&(clt[nNumClut].sum_mass),
		       clt[nNumClut].num_atom_clust,
		       clt[nNumClut].xoord,
		       clt[nNumClut].mass);
    //    clt[nNumClut].join = joinflag;                      // 2011-11
    //    if (clt[nNumClut].num_branch > 1) joinflag +=1;     // 2011-11
    //    if (clt[nNumClut].terminal == 0) joinflag -=1;      // 2011-11
  }

  joinflag=0;
  for(nNumClut=0; nNumClut<nNumClut_all; ++nNumClut) {
    joinflag=ABA_setJoin(clt,nNumClut,joinflag);
  }

}

void ABAs_Inertia_clust(double Inertia_clust[3][3],
			int num_atom_clust,double *crd_local,double *mass) {
  int i,j,k;

  Inertia_clust[0][0]=0.0;
  Inertia_clust[1][1]=0.0;
  Inertia_clust[2][2]=0.0;
  Inertia_clust[0][1]=0.0;
  Inertia_clust[0][2]=0.0;
  Inertia_clust[1][2]=0.0;
	
  for(i=0;i<num_atom_clust;++i) {
    Inertia_clust[0][0] += 
      mass[i]*(crd_local[i*3+1]*crd_local[i*3+1]+crd_local[i*3+2]*crd_local[i*3+2]);
    
    Inertia_clust[1][1] += 
      mass[i]*(crd_local[i*3+0]*crd_local[i*3+0]+crd_local[i*3+2]*crd_local[i*3+2]);
    
    Inertia_clust[2][2] += 
      mass[i]*(crd_local[i*3+0]*crd_local[i*3+0]+crd_local[i*3+1]*crd_local[i*3+1]);
    
    Inertia_clust[0][1] -=  mass[i]*crd_local[i*3+0]*crd_local[i*3+1];
    
    Inertia_clust[0][2] -=  mass[i]*crd_local[i*3+0]*crd_local[i*3+2];
    
    Inertia_clust[1][2] -=  mass[i]*crd_local[i*3+1]*crd_local[i*3+2];
  }
  
  Inertia_clust[1][0]=Inertia_clust[0][1];
  Inertia_clust[2][0]=Inertia_clust[0][2];
  Inertia_clust[2][1]=Inertia_clust[1][2];
}

void ABAs_InertiaMatrix(double InertiaMatrix[6][6],
			double Inertia_clust[3][3],
			double qcom[3], double *summass,
			int num_atom_clust,
			double *crd_local,double *mass) {
  int i;
  int nNumAtom;
  int alpha, alpha2;
  double mq[3];

  *summass=0.0;
  for(nNumAtom=0;nNumAtom<num_atom_clust;++nNumAtom) *summass+=mass[nNumAtom];
  
  for(alpha=0;alpha<3;++alpha) 
    for(alpha2=0;alpha2<3;++alpha2) 
      InertiaMatrix[alpha][alpha2] = Inertia_clust[alpha][alpha2];

  InertiaMatrix[3][3]=*summass;
  InertiaMatrix[3][4]=0.0;
  InertiaMatrix[3][5]=0.0;
  InertiaMatrix[4][3]=0.0;
  InertiaMatrix[4][4]=InertiaMatrix[3][3];
  InertiaMatrix[4][5]=0.0;
  InertiaMatrix[5][3]=0.0;
  InertiaMatrix[5][4]=0.0;
  InertiaMatrix[5][5]=InertiaMatrix[3][3];

  for (alpha=0;alpha<3;++alpha)  mq[alpha] = 0.0;  
  for(nNumAtom=0;nNumAtom<num_atom_clust;++nNumAtom)
    mq[0] +=  mass[nNumAtom]*crd_local[nNumAtom*3+0];
  
  for(nNumAtom=0;nNumAtom<num_atom_clust;++nNumAtom)
    mq[1] +=  mass[nNumAtom]*crd_local[nNumAtom*3+1];

  for(nNumAtom=0;nNumAtom<num_atom_clust;++nNumAtom)
    mq[2] +=  mass[nNumAtom]*crd_local[nNumAtom*3+2];

  InertiaMatrix[0][3]=0.0;
  InertiaMatrix[0][4]=-mq[2];
  InertiaMatrix[0][5]=mq[1];
  InertiaMatrix[1][3]=mq[2];
  InertiaMatrix[1][4]=0.0;
  InertiaMatrix[1][5]=-mq[0];
  InertiaMatrix[2][3]=-mq[1];
  InertiaMatrix[2][4]=mq[0];
  InertiaMatrix[2][5]=0.0;

  InertiaMatrix[3][0]=0.0;
  InertiaMatrix[3][1]=mq[2];
  InertiaMatrix[3][2]=-mq[1];
  InertiaMatrix[4][0]=-mq[2];
  InertiaMatrix[4][1]=0.0;
  InertiaMatrix[4][2]=mq[0];
  InertiaMatrix[5][0]=mq[1];
  InertiaMatrix[5][1]=-mq[0];
  InertiaMatrix[5][2]=0.0;

  for (i=0;i<3;++i) qcom[i]=mq[i]/(*summass);
  
}


int ABA_setJoin(CLT *clt,int nNumClut, int joinflag) {

  clt[nNumClut].join = joinflag;
  
  if (clt[nNumClut].num_branch > 1) {
    joinflag +=1;
  }
  if (clt[nNumClut].terminal == 0) {
    joinflag -=1;
  }

  return joinflag;
}
