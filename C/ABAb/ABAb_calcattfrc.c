
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ABAb.h"
#include "EF.h"
#include "TOPO.h"
#include "PTL.h"
//#include "PT.h"

double ABAb_calcAttTorque_H(double *Q, CLTb *clt,double *crd,double *q,double *qgoal,double kc,int numclut,FILE *outputfile2,int intervalflag){
  int i,j,k,t;
  double atom[4][3];
  int *num;
  int **atomdihedpairs;
  int *pairs;
  double d;
  double e=0.0;
  double /**q,*/pi;
  double UNIT=418.4070;

  pi=acos(-1.0);

  ABAb_set_ini(q,clt,crd,numclut);

  for (i=1;i<numclut;++i){
    d=q[i]-qgoal[i];

    if (fabs(d)>pi) {
      if (d>0.0) {
	d=-(2.0*pi-fabs(d));
      }
      else {
	d=2.0*pi-fabs(d);
      }
    }
    Q[i]=kc*d*UNIT;
    e+=0.5*kc*d*d;
    
    if (intervalflag==0) {
      fprintf(outputfile2,"%lf ",d);
    }
    //    e+=kc*(1.0+cos(q[i]-qgoal[i]));
    //    Q[i]=-kc*(sin(q[i]-qgoal[i]))*UNIT;
  }

  //  for (i=0;i<numclut;++i) fprintf(outputfile2,"%lf ",q[i]);
  if (intervalflag==0) {
    fprintf(outputfile2,"\n");
  }

  return e;
}

double ABAb_calcAttTorque_Hwd(double *Q,CLTb *clt,double *crd,double *q,double *qvel,double *qgoal,double kc,double kd,int numclut,FILE *outputfile2){
  int i,j,k,t;
  double atom[4][3];
  int *num;
  int **atomdihedpairs;
  int *pairs;
  double d;
  double e=0.0;
  double /**q,*/pi;
  double UNIT=418.4070;

  pi=acos(-1.0);

  ABAb_set_ini(q,clt,crd,numclut);

  for (i=1;i<numclut;++i){
    d=q[i]-qgoal[i];

    if (fabs(d)>pi) {
      if (d>0.0) {
	d=-(2.0*pi-fabs(d));
      }
      else {
	d=2.0*pi-fabs(d);
      }
    }
    Q[i]=kc*d*UNIT-kd*qvel[i]*UNIT;
    //    e+=0.5*kc*d*d;
    
    fprintf(outputfile2,"%lf ",d);
  }

  fprintf(outputfile2,"\n");

  return e;
}

double ABAb_calcAttTorque_C(double *Q,CLTb *clt,double *crd,double *q,double *qgoal,double kc,int numclut,FILE *outputfile2,double GoalReachedThershold){
  int i,j,k,t;
  double atom[4][3];
  int *num;
  int **atomdihedpairs;
  int *pairs;
  double d;
  double e=0.0;
  double /**q,*/pi;
  double UNIT=418.4070;

  pi=acos(-1.0);

  ABAb_set_ini(q,clt,crd,numclut);

  for (i=1;i<numclut;++i){
    d=q[i]-qgoal[i];

    if (fabs(d)>pi) {
      if (d>0.0) {
	d=-(2.0*pi-fabs(d));
      }
      else {
	d=2.0*pi-fabs(d);
      }
    }
    if (d<0.0 && fabs(d)>GoalReachedThershold)
      Q[i]=-1.0*kc*UNIT;
    else if (d>=0.0 && fabs(d)>GoalReachedThershold)
      Q[i]=kc*UNIT;
    else 
      Q[i]=0.0;
    //    e+=0.5*kc*d*d;
    
    fprintf(outputfile2,"%lf ",d);
  }
  fprintf(outputfile2,"\n");

  return e;
}


void ABAb_set_ini(double *q,CLTb *clt,double *crd,int numclut){
  int i,j,k,t,p,l,ll;
  double atom[4][3];
  int *num;
  int **atomdihedpairs;
  int *pairs;
  double pi;

  pi=acos(-1.0);

  pairs=(int *)gcemalloc(sizeof(int)*numclut*4);

  atomdihedpairs=(int **)gcemalloc(sizeof(int *)*5);

  atomdihedpairs[0]=(int *)gcemalloc(sizeof(int)*4); // PHI,PSI
  atomdihedpairs[1]=(int *)gcemalloc(sizeof(int)*4); // OMEGA
  atomdihedpairs[2]=(int *)gcemalloc(sizeof(int)*4); // KI
  atomdihedpairs[3]=(int *)gcemalloc(sizeof(int)*8); // ACE
  atomdihedpairs[4]=(int *)gcemalloc(sizeof(int)*8); // NME
  num=(int *)gcemalloc(sizeof(int)*5);
  readdihedpairsL(atomdihedpairs,num);
  //  readdihedpairs(atomdihedpairs,num);

  for (t=0;t<5;++t) {
    for (i=0;i<num[t];++i) {
      for (j=0;j<numclut;++j) {
	p = clt[j].nNumClutOfParent-1;
	l = 1;
	ll = 2;
	if (atomdihedpairs[t][i*4+1] > atomdihedpairs[t][i*4+2]) {
	  l = 2;
	  ll = 1;
	}
	//	if (clt[j].terminal_atom_a[0]-1==atomdihedpairs[t][i*4+1]) {
	if (atomdihedpairs[t][i*4+l]==clt[p].terminal_atom_a[0]-1 && 
	    atomdihedpairs[t][i*4+ll]==clt[j].origin_atom_a-1){
	  for (k=0;k<4;++k) {
	    pairs[(j)*4+k]=atomdihedpairs[t][i*4+k];
	  }
	  break;
	}
      }
    }
  }

  for (i=1;i<numclut;++i) {
    for (j=0;j<4;++j) for (k=0;k<3;++k) atom[j][k]=crd[(pairs[i*4+j])*3+k];
    q[i]=dih(atom[0],atom[1],atom[2],atom[3]);
    //    if (q[i]>pi) q[i]-=2.0*pi;
    //    else if (q[i]<-1.0*pi) q[i]+=2.0*pi;
  }
}
