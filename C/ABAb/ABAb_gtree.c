
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "EF.h"
#include "RAND.h"
#include "ABAb.h"

void SampleRandom(double *qrand,double *crdbase,CLTb *clt,int numclut,int numatom){
  /*****************************************************/
  /* int i,j;					       */
  /* double u,pi;				       */
  /* double *qrot;				       */
  /* 						       */
  /* pi=acos(-1.0);				       */
  /* 						       */
  /* qrot=(double *)gcemalloc(sizeof(double)*numclut); */
  /* 						       */
  /* ABAb_set_ini(qrot,clt,crdbase,numclut);	       */
  /* 						       */
  /* for (i=0;i<numclut;++i) {			       */
  /*   u=genrand_real2();			       */
  /*   qrand[i]=u*2.0*pi;			       */
  /*   qrot[i]-=qrand[i];			       */
  /* }						       */
  /* 						       */
  /* 						       */
  /* ABAb_update(clt,crdrand,qrot,numclut,numatom);     */
  /*****************************************************/

}

void Createnewcrd(double *crdnew,double *crdbase,double *qrand,CLTb *clt,int numclut,int numatom,double f){
  int i,j;
  double *qrot,coff;
  double maxq;
  double pi;

  pi=acos(-1.0);

  maxq=qrand[0];
  for (i=1;i<numclut;++i) {
    if (maxq < qrand[i]) {
      maxq=qrand[i];
    }
  }
  coff=f/maxq;
  qrot=(double *)gcemalloc(sizeof(double)*numclut);
  for (i=0;i<numclut;++i) qrot[i]=coff*qrand[i];

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j)
      crdnew[i*3+j]=crdbase[i*3+j];
  ABAb_update(clt,crdnew,qrot,numclut,numatom);
}

//void Findnearconfig(double **tree,double *crdrand,double *crdnode,int treesize){
//  int i;
//
//
//}

