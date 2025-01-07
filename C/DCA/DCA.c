
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "DCA.h"
#include "EF.h"
#include "LA.h"
#include "f2c.h"
#include "clapack.h"

#define ON 1
#define OFF 0

double solverDCA(double *qacc,double *qvel,CLT *clt,double *Q,double *frc,AST *assembletree,
		 int numclt,int numatom,int numtree){
  int i,j,k,l;
  int leafflag;
  int na,nb,nc;
  int naref,nbref;
  double *temp;
  double *PA1,*PA2,*PA12,*PA21,*bA1,*bA2;
  DCA *dca;

  double *test,invtest[6][6];

  dca=(DCA *)gcemalloc(sizeof(DCA)*numtree);
  for (i=0;i<numtree;++i) {
    dca[i].P1=(double *)gcemalloc(sizeof(double)*6*6);
    dca[i].P2=(double *)gcemalloc(sizeof(double)*6*6);
    dca[i].P12=(double *)gcemalloc(sizeof(double)*6*6);
    dca[i].P21=(double *)gcemalloc(sizeof(double)*6*6);
    dca[i].f1=(double *)gcemalloc(sizeof(double)*6);
    dca[i].f2=(double *)gcemalloc(sizeof(double)*6);
    dca[i].b1=(double *)gcemalloc(sizeof(double)*6);
    dca[i].b2=(double *)gcemalloc(sizeof(double)*6);
    dca[i].W=(double *)gcemalloc(sizeof(double)*6*6);
    dca[i].V=(double *)gcemalloc(sizeof(double)*6*6);
    dca[i].gamma=(double *)gcemalloc(sizeof(double)*6);
    dca[i].beta=(double *)gcemalloc(sizeof(double)*6);
  }
  temp=(double *)gcemalloc(sizeof(double)*6*6);
  PA1=(double *)gcemalloc(sizeof(double)*6*6);
  PA2=(double *)gcemalloc(sizeof(double)*6*6);
  PA12=(double *)gcemalloc(sizeof(double)*6*6);
  PA21=(double *)gcemalloc(sizeof(double)*6*6);
  bA1=(double *)gcemalloc(sizeof(double)*6);
  bA2=(double *)gcemalloc(sizeof(double)*6*6);

  //  DCAp_prepass(clt,qvel,numclt,numatom);

  for (i=numtree-1;i>=0;--i) {
    leafflag=assembletree[i].leafflag;
    if (leafflag==ON) {
      nc=assembletree[i].right-1;
      for (j=0;j<6;++j)	for (k=0;k<6;++k) temp[j*6+k]=clt[nc].IM[j][k];
      invm2(temp,dca[i].P1,6);

      for (j=0;j<6;++j) {
	for (k=0;k<6;++k) {
	  dca[i].P2[j*6+k]=dca[i].P1[j*6+k];
	  dca[i].P12[j*6+k]=dca[i].P1[j*6+k];
	  dca[i].P21[j*6+k]=dca[i].P1[j*6+k];
	}
      }
      mvmult(dca[i].P1,clt[i].Spfrc,dca[i].b1,6);
    }
    else {
      na=assembletree[i].right-1;
      nb=assembletree[i].left-1;
      nc=assembletree[i].num-1;
      naref=assembletree[i].refright-1;
      nbref=assembletree[i].refleft-1;
      DCAs_trans_mainpass(naref,nbref,clt,
			  PA1,PA2,PA21,bA1,bA2,
			  dca[na].P1,dca[na].P2,dca[na].P21,dca[na].b1,dca[na].b2);
      DCAm_mainpass(dca[nc].P1,  dca[nc].P2, 
		    dca[nc].P12, dca[nc].P21,
		    dca[nc].b1,  dca[nc].b2,
		    PA1,         PA2,
		    PA12,        PA21,
		    bA1,         bA2,
 		    dca[nb].P1,  dca[nb].P2, 
		    dca[nb].P12, dca[nb].P21, 
		    dca[nb].b1,  dca[nb].b2,
		    clt[nc].Coacc,Q[nc]);
    }
  }

  for (i=0;assembletree[i]!=-1;++i) {
    na=assembletree[i].right-1;
    nb=assembletree[i].left-1
    nc=assembletree[i].num-1;
    naref=assembletree[i].refright-1;
    nbref=assembletree[i].refleft-1;
    DCAs_trans_mainpass(naref,nbref,clt,
			PA1,PA2,PA21,bA1,bA2,
			dca[na].P1,dca[na].P2,dca[na].P21,dca[na].b1,dca[na].b2);
    qacc[nc]=DCAb_backpass(dca[na].f1,dca[nb].f2,dca[na].f2,dca[nb].f1,
  			   dca[na].P21,dca[nb].P12,dca[nc].V,
  			   clt[nc].Q,dca[nc].beta,dca[nc].gamma,dca[nc].W);
  }
}
