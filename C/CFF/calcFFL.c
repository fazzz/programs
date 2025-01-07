
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "FFL.h"
#include "MB.h"
#include "PTL.h"
#include "LA.h"
#include "TOPO.h"
#include "mymath.h"
#include "EF.h"

#define YES 1
#define NO  0

double calcANGKE_force(double atomi[3],double atomj[3],double atomk[3],double kang,double ang_eq,double *f);

int ffL_calcFFNB(double *ele, double *ALJ, double *BLJ,
		double *p_e,double *p_LJ,
		double *f_e,double *f_LJ,
		int numnb, int *indexnb,
		int num_atom,
		double *cord,
		int flagp, int flagf) {
  int i,j;
  int num_a_prot,NUM_A_PROT;
  double potedummy;
  double f[3];
  double len,len2,len6;
  double vec[3],fdummy[3];
  double p12,p6;

  if (flagf != 0 && flagf != 1 && flagf !=2 ) {
    printf("error !\n");exit(1);
  }
  if (flagp != 0 && flagp != 1 && flagp !=2 ) {
    printf("error !\n");exit(1);
  }

  if (flagp==1 ) {
    for(i=0;i<num_atom;++i) {
      p_e[i]=0.0;p_LJ[i]=0.0;
      if (flagf==1) {
	for(j=0;j<3;++j) {
	  f_LJ[i*3+j]=0.0;
	  f_e[i*3+j]=0.0;
	}
      }
    }
  }

  if (flagp==2 ) {
    for(i=0;i<numnb;++i) {
      p_e[i]=0.0;
      p_LJ[i]=0.0; 
      if (flagf==2) {
	for(j=0;j<3;++j) {
	  f_LJ[i*3+j]=0.0;
	  f_e[i*3+j]=0.0;
	}
      }
    }
  }

  for(i=0;i<numnb;++i){
    num_a_prot=indexnb[i*2];
    NUM_A_PROT=indexnb[i*2+1];
    len2 = 0.0;
    for(j=0;j<3;++j){
      vec[j] = cord[NUM_A_PROT*3+j]-cord[num_a_prot*3+j];
      len2 += vec[j]*vec[j];
    }
    len = sqrt(len2);
    potedummy=ele[num_a_prot]*ele[NUM_A_PROT]/(len);
    if (flagp==1) {
      p_e[num_a_prot] += potedummy;
      p_e[NUM_A_PROT] += potedummy;
    }
    else if (flagp==2) p_e[i] = potedummy;

    if (flagf==1) {
      for(j=0; j<3; ++j) {
	fdummy[j] = -potedummy*UNIT*vec[j]/len2;
	f_e[num_a_prot*3+j] += fdummy[j];
	f_e[NUM_A_PROT*3+j] -= fdummy[j];
      }
    }
    else if (flagf==2) {
      for(j=0; j<3; ++j) {
	fdummy[j] = -potedummy*UNIT*vec[j]/len2;
	f_e[i*3+j] += fdummy[j];
      }
    }
    len6=len2;
    for (j=0;j<2;++j)  len6 = len6*len2;
    p12 = ALJ[num_a_prot*num_atom+NUM_A_PROT]/(len6*len6);
    p6  = BLJ[num_a_prot*num_atom+NUM_A_PROT]/len6;
    if (flagp==1) {
      p_LJ[num_a_prot] += p12-p6;
      p_LJ[NUM_A_PROT] += p12-p6;
    }
    else if (flagp==2) p_LJ[i] = p12-p6;

    if (flagf==1) {
      for(j=0;j<3;++j) {
	fdummy[j]=-6.0*(2.0*p12-p6)/(len2)*vec[j]*UNIT;
	f_LJ[num_a_prot*3+j] +=fdummy[j];
	f_LJ[NUM_A_PROT*3+j] -= fdummy[j];
      }
    }
    else if (flagf==2) {
      for(j=0;j<3;++j) {
	fdummy[j]=-6.0*(2.0*p12-p6)/(len2)*vec[j]*UNIT;
	f_LJ[i*3+j] +=fdummy[j];
      }
    }
  }

  return 0;
}

int ffL_calcFFNB_wtune(double *ele, double *ALJ, double *BLJ,
		       double *p_e,double *p_LJ,
		       double *f_e,double *f_LJ,
		       int numnb, int *indexnb,
		       int num_atom,
		       double *cord,
		       int flagp, int flagf,
		       int *atom_tune_pairs, double *tune_val, int numtune) {
  int i,j;
  int num_a_prot,NUM_A_PROT;
  double potedummy;
  double f[3];
  double len,len2,len6;
  double vec[3],fdummy[3];
  double p12,p6;

  int ON=1;
  int OFF=0;

  int tuneflag=OFF,tuneindex;

  if (flagf != 0 && flagf != 1 && flagf !=2 ) {
    printf("error !\n");exit(1);
  }
  if (flagp != 0 && flagp != 1 && flagp !=2 ) {
    printf("error !\n");exit(1);
  }

  if (flagp==1 ) {
    for(i=0;i<num_atom;++i) {
      p_e[i]=0.0;p_LJ[i]=0.0;
      if (flagf==1) {
	for(j=0;j<3;++j) {
	  f_LJ[i*3+j]=0.0;
	  f_e[i*3+j]=0.0;
	}
      }
    }
  }

  if (flagp==2 ) {
    for(i=0;i<numnb;++i) {
      p_e[i]=0.0;
      p_LJ[i]=0.0; 
      if (flagf==2) {
	for(j=0;j<3;++j) {
	  f_LJ[i*3+j]=0.0;
	  f_e[i*3+j]=0.0;
	}
      }
    }
  }

  for(i=0;i<numnb;++i){
    num_a_prot=indexnb[i*2];
    NUM_A_PROT=indexnb[i*2+1];
    tuneflag=OFF;
    for (j=0;j<numtune;++j) {
      if (num_a_prot==atom_tune_pairs[j*2] && NUM_A_PROT==atom_tune_pairs[j*2+1]) {
	tuneindex=j;
	tuneflag=ON;
      }
    }

    len2 = 0.0;
    for(j=0;j<3;++j){
      vec[j] = cord[NUM_A_PROT*3+j]-cord[num_a_prot*3+j];
      len2 += vec[j]*vec[j];
    }
    len = sqrt(len2);
    if (tuneflag==ON)
      potedummy=tune_val[tuneindex]*ele[num_a_prot]*ele[NUM_A_PROT]/(len);
    else
      potedummy=ele[num_a_prot]*ele[NUM_A_PROT]/(len);
    if (flagp==1) {
      p_e[num_a_prot] += potedummy;
      p_e[NUM_A_PROT] += potedummy;
    }
    else if (flagp==2) p_e[i] = potedummy;

    if (flagf==1) {
      for(j=0; j<3; ++j) {
	fdummy[j] = -potedummy*UNIT*vec[j]/len2;
	f_e[num_a_prot*3+j] += fdummy[j];
	f_e[NUM_A_PROT*3+j] -= fdummy[j];
      }
    }
    else if (flagf==2) {
      for(j=0; j<3; ++j) {
	fdummy[j] = -potedummy*UNIT*vec[j]/len2;
	f_e[i*3+j] += fdummy[j];
      }
    }
    len6=len2;
    for (j=0;j<2;++j)  len6 = len6*len2;
    if (tuneflag==ON) {
      p12 = tune_val[tuneindex]*ALJ[num_a_prot*num_atom+NUM_A_PROT]/(len6*len6);
      p6  = tune_val[tuneindex]*BLJ[num_a_prot*num_atom+NUM_A_PROT]/len6;
    }
    else {
      p12 = ALJ[num_a_prot*num_atom+NUM_A_PROT]/(len6*len6);
      p6  = BLJ[num_a_prot*num_atom+NUM_A_PROT]/len6;
    }
    if (flagp==1) {
      p_LJ[num_a_prot] += p12-p6;
      p_LJ[NUM_A_PROT] += p12-p6;
    }
    else if (flagp==2) p_LJ[i] = p12-p6;

    if (flagf==1) {
      for(j=0;j<3;++j) {
	fdummy[j]=-6.0*(2.0*p12-p6)/(len2)*vec[j]*UNIT;
	f_LJ[num_a_prot*3+j] +=fdummy[j];
	f_LJ[NUM_A_PROT*3+j] -= fdummy[j];
      }
    }
    else if (flagf==2) {
      for(j=0;j<3;++j) {
	fdummy[j]=-6.0*(2.0*p12-p6)/(len2)*vec[j]*UNIT;
	f_LJ[i*3+j] +=fdummy[j];
      }
    }
  }

  return 0;
}

int ffL_calcFFNB_wtuneb(double *ele, double *ALJ, double *BLJ,
			double *p_e,double *p_LJ,
			double *f_e,double *f_LJ,
			int numnb, int *indexnb,
			int num_atom,
			double *cord,
			int flagp, int flagf,
			int *atom_tune_pairs_es, int *atom_tune_pairs_LJ,
			double *tune_val_es, double *tune_val_LJ,
			int numtune_es, int numtune_LJ) {
  int i,j;
  int num_a_prot,NUM_A_PROT;
  double potedummy;
  double f[3];
  double len,len2,len6;
  double vec[3],fdummy[3];
  double p12,p6;

  int ON=1;
  int OFF=0;

  int tuneflag_es=OFF,tuneindex_es;
  int tuneflag_LJ=OFF,tuneindex_LJ;

  if (flagf != 0 && flagf != 1 && flagf !=2 ) {
    printf("error !\n");exit(1);
  }
  if (flagp != 0 && flagp != 1 && flagp !=2 ) {
    printf("error !\n");exit(1);
  }

  if (flagp==1 ) {
    for(i=0;i<num_atom;++i) {
      p_e[i]=0.0;p_LJ[i]=0.0;
      if (flagf==1) {
	for(j=0;j<3;++j) {
	  f_LJ[i*3+j]=0.0;
	  f_e[i*3+j]=0.0;
	}
      }
    }
  }

  if (flagp==2 ) {
    for(i=0;i<numnb;++i) {
      p_e[i]=0.0;
      p_LJ[i]=0.0;
      if (flagf==2) {
	for(j=0;j<3;++j) {
	  f_LJ[i*3+j]=0.0;
	  f_e[i*3+j]=0.0;
	}
      }
    }
  }

  for(i=0;i<numnb;++i){
    num_a_prot=indexnb[i*2];
    NUM_A_PROT=indexnb[i*2+1];
    tuneflag_es=OFF;
    for (j=0;j<numtune_es;++j) {
      if (num_a_prot==atom_tune_pairs_es[j*2] && NUM_A_PROT==atom_tune_pairs_es[j*2+1]) {
	tuneindex_es=j;
	tuneflag_es=ON;
      }
    }
    tuneflag_LJ=OFF;
    for (j=0;j<numtune_LJ;++j) {
      if (num_a_prot==atom_tune_pairs_LJ[j*2] && NUM_A_PROT==atom_tune_pairs_LJ[j*2+1]) {
	tuneindex_LJ=j;
	tuneflag_LJ=ON;
      }
    }

    len2 = 0.0;
    for(j=0;j<3;++j){
      vec[j] = cord[NUM_A_PROT*3+j]-cord[num_a_prot*3+j];
      len2 += vec[j]*vec[j];
    }
    len = sqrt(len2);
    if (tuneflag_es==ON)
      potedummy=tune_val_es[tuneindex_es]*ele[num_a_prot]*ele[NUM_A_PROT]/(len);
    else
      potedummy=ele[num_a_prot]*ele[NUM_A_PROT]/(len);
    if (flagp==1) {
      p_e[num_a_prot] += potedummy;
      p_e[NUM_A_PROT] += potedummy;
    }
    else if (flagp==2) p_e[i] = potedummy;

    if (flagf==1) {
      for(j=0; j<3; ++j) {
	fdummy[j] = -potedummy*UNIT*vec[j]/len2;
	f_e[num_a_prot*3+j] += fdummy[j];
	f_e[NUM_A_PROT*3+j] -= fdummy[j];
      }
    }
    else if (flagf==2) {
      for(j=0; j<3; ++j) {
	fdummy[j] = -potedummy*UNIT*vec[j]/len2;
	f_e[i*3+j] += fdummy[j];
      }
    }
    len6=len2;
    for (j=0;j<2;++j)  len6 = len6*len2;
    if (tuneflag_LJ==ON) {
      p12 = tune_val_LJ[tuneindex_LJ]*ALJ[num_a_prot*num_atom+NUM_A_PROT]/(len6*len6);
      p6  = tune_val_LJ[tuneindex_LJ]*BLJ[num_a_prot*num_atom+NUM_A_PROT]/len6;
    }
    else {
      p12 = ALJ[num_a_prot*num_atom+NUM_A_PROT]/(len6*len6);
      p6  = BLJ[num_a_prot*num_atom+NUM_A_PROT]/len6;
    }
    if (flagp==1) {
      p_LJ[num_a_prot] += p12-p6;
      p_LJ[NUM_A_PROT] += p12-p6;
    }
    else if (flagp==2) p_LJ[i] = p12-p6;

    if (flagf==1) {
      for(j=0;j<3;++j) {
	fdummy[j]=-6.0*(2.0*p12-p6)/(len2)*vec[j]*UNIT;
	f_LJ[num_a_prot*3+j] +=fdummy[j];
	f_LJ[NUM_A_PROT*3+j] -= fdummy[j];
      }
    }
    else if (flagf==2) {
      for(j=0;j<3;++j) {
	fdummy[j]=-6.0*(2.0*p12-p6)/(len2)*vec[j]*UNIT;
	f_LJ[i*3+j] +=fdummy[j];
      }
    }
  }

  return 0;
}

int ffL_calcFFNB_14(double *ele, double *ALJ, double *BLJ,
		    double *p_e,double *p_LJ,
		    double *f_e,double *f_LJ,
		    int numnb, int *indexnb,
		    int num_atom,
		    double *cord,
		    int flagp, int flagf) {
  int i,j;
  int num_a_prot,NUM_A_PROT;
  double potedummy;
  double f[3];
  double len,len2,len6;
  double vec[3],fdummy[3];
  double p12,p6;

  if (flagf != 0 && flagf != 1 && flagf !=2 ) {
    printf("error !\n");exit(1);
  }
  if (flagp != 0 && flagp != 1 && flagp !=2 ) {
    printf("error !\n");exit(1);
  }

  if (flagp==1 ) {
    for(i=0;i<num_atom;++i) {
      p_e[i]=0.0;p_LJ[i]=0.0;
      if (flagf==1) {
	for(j=0;j<3;++j) {
	  f_LJ[i*3+j]=0.0;
	  f_e[i*3+j]=0.0;
	}
      }
    }
  }

  if (flagp==2 ) {
    for(i=0;i<numnb;++i) {
      p_e[i]=0.0;
      p_LJ[i]=0.0; 
      if (flagf==2) {
	for(j=0;j<3;++j) {
	  f_LJ[i*3+j]=0.0;
	  f_e[i*3+j]=0.0;
	}
      }
    }
  }

  for(i=0;i<numnb;++i){
    num_a_prot=indexnb[i*2];
    NUM_A_PROT=indexnb[i*2+1];
    len2 = 0.0;
    for(j=0;j<3;++j){
      vec[j] = cord[NUM_A_PROT*3+j]-cord[num_a_prot*3+j];
      len2 += vec[j]*vec[j];
    }
    len = sqrt(len2);
    potedummy=/*1.0/1.2**/ele[num_a_prot]*ele[NUM_A_PROT]/(len);
    if (flagp==1) {
      p_e[num_a_prot] += 1.0/1.2*potedummy;
      p_e[NUM_A_PROT] += 1.0/1.2*potedummy;
    }
    else if (flagp==2) p_e[i] = 1.0/1.2*potedummy;

    if (flagf==1) {
      for(j=0; j<3; ++j) {
	fdummy[j] = -1.0/1.2*potedummy*/*UNIT**/4.184070*100.0*vec[j]/len2;
	f_e[num_a_prot*3+j] += fdummy[j];
	f_e[NUM_A_PROT*3+j] -= fdummy[j];
      }
    }
    else if (flagf==2) {
      for(j=0; j<3; ++j) {
	fdummy[j] = -1.0/1.2*potedummy*/*UNIT**/4.184070*100.0*vec[j]/len2;
	f_e[i*3+j] += fdummy[j];
      }
    }
    len6=len2;
    for (j=0;j<2;++j)  len6 = len6*len2;
    p12 = ALJ[num_a_prot*num_atom+NUM_A_PROT]/(len6*len6);
    p6  = BLJ[num_a_prot*num_atom+NUM_A_PROT]/len6;
    if (flagp==1) {
      p_LJ[num_a_prot] += 0.5*(p12-p6);
      p_LJ[NUM_A_PROT] += 0.5*(p12-p6);
    }
    else if (flagp==2) p_LJ[i] = 0.5*(p12-p6);

    if (flagf==1) {
      for(j=0;j<3;++j) {
	fdummy[j]=-3.0*(2.0*p12-p6)/(len2)*vec[j]*UNIT;
	f_LJ[num_a_prot*3+j] +=fdummy[j];
	f_LJ[NUM_A_PROT*3+j] -= fdummy[j];
      }
    }
    else if (flagf==2) {
      for(j=0;j<3;++j) {
	fdummy[j]=-3.0*(2.0*p12-p6)/(len2)*vec[j]*UNIT;
	f_LJ[i*3+j] +=fdummy[j];
      }
    }
  }

  return 0;
}

int ffL_calcFFNB_14_6(double *ele, double *ALJ, double *BLJ, // 0911
		     double *p_e,double *p_LJ,
		     double *f_e,double *f_LJ,
		     int numnb, int *indexnb,
		     int num_atom,
		     double *cord,
		     int flagp, int flagf,
		     double scale,
		     int numrepul) {
  int i,j;
  int num_a_prot,NUM_A_PROT;
  double potedummy;
  double f[3];
  double len,len2,len6,len14;
  double vec[3],fdummy[3];
  double p14,p6;

  if (flagf != 0 && flagf != 1 && flagf !=2 ) {
    printf("error !\n");exit(1);
  }
  if (flagp != 0 && flagp != 1 && flagp !=2 ) {
    printf("error !\n");exit(1);
  }

  if (flagp==1 ) {
    for(i=0;i<num_atom;++i) {
      p_e[i]=0.0;p_LJ[i]=0.0;
      if (flagf==1) {
	for(j=0;j<3;++j) {
	  f_LJ[i*3+j]=0.0;
	  f_e[i*3+j]=0.0;
	}
      }
    }
  }

  if (flagp==2 ) {
    for(i=0;i<numnb;++i) {
      p_e[i]=0.0;
      p_LJ[i]=0.0; 
      if (flagf==2) {
	for(j=0;j<3;++j) {
	  f_LJ[i*3+j]=0.0;
	  f_e[i*3+j]=0.0;
	}
      }
    }
  }

  for(i=0;i<numnb;++i){
    num_a_prot=indexnb[i*2];
    NUM_A_PROT=indexnb[i*2+1];
    len2 = 0.0;
    for(j=0;j<3;++j){
      vec[j] = cord[NUM_A_PROT*3+j]-cord[num_a_prot*3+j];
      len2 += vec[j]*vec[j];
    }
    len = sqrt(len2);
    potedummy=ele[num_a_prot]*ele[NUM_A_PROT]/(len);
    if (flagp==1) {
      p_e[num_a_prot] += potedummy;
      p_e[NUM_A_PROT] += potedummy;
    }
    else if (flagp==2) p_e[i] = potedummy;

    if (flagf==1) {
      for(j=0; j<3; ++j) {
	fdummy[j] = -potedummy*UNIT*vec[j]/len2;
	f_e[num_a_prot*3+j] += fdummy[j];
	f_e[NUM_A_PROT*3+j] -= fdummy[j];
      }
    }
    else if (flagf==2) {
      for(j=0; j<3; ++j) {
	fdummy[j] = -potedummy*UNIT*vec[j]/len2;
	f_e[i*3+j] += fdummy[j];
      }
    }
    len6=len2;
    for (j=0;j<2;++j)  len6 = len6*len2;
    len14 = len;
    for (j=0;j<numrepul-1;++j)  len14 = len14*len;
    p14 = scale*ALJ[num_a_prot*num_atom+NUM_A_PROT]/(len14);
    p6  = BLJ[num_a_prot*num_atom+NUM_A_PROT]/len6;
    if (flagp==1) {
      p_LJ[num_a_prot] += p14-p6;
      p_LJ[NUM_A_PROT] += p14-p6;
    }
    else if (flagp==2) p_LJ[i] = p14-p6;

    if (flagf==1) {
      for(j=0;j<3;++j) {
	fdummy[j]=-6.0*(2.0*p14-p6)/(len2)*vec[j]*UNIT;
	f_LJ[num_a_prot*3+j] +=fdummy[j];
	f_LJ[NUM_A_PROT*3+j] -= fdummy[j];
      }
    }
    else if (flagf==2) {
      for(j=0;j<3;++j) {
	fdummy[j]=-6.0*(2.0*p14-p6)/(len2)*vec[j]*UNIT;
	f_LJ[i*3+j] +=fdummy[j];
      }
    }
  }

  return 0;
}


int ffL_calcFFNB_wpep(double *ele, double *ALJ, double *BLJ,
		     double *p_e,double *p_LJ,
		     double *f_e,double *f_LJ,
		     int numnb, int *indexnb,
		     int num_atom,
		     double *cord,
		     double *u_1_5, double fact,
		     int flagp, int flagf) {
  int i,j;
  int num_a_prot,NUM_A_PROT;
  double potedummy;
  double f[3];
  double len,len2,len6;
  double vec[3],fdummy[3];
  double p12,p6;

  if (flagf != 0 && flagf != 1 && flagf !=2 ) {
    printf("error !\n");exit(1);
  }
  if (flagp != 0 && flagp != 1 && flagp !=2 ) {
    printf("error !\n");exit(1);
  }

  if (flagp==1 ) {
    for(i=0;i<num_atom;++i) {
      p_e[i]=0.0;p_LJ[i]=0.0;
      if (flagf==1) {
	for(j=0;j<3; ++j) {
	  f_LJ[i*3+j]=0.0;f_e[i*3+j]=0.0;
	}
      }
    }
  }

  if (flagp==2 ) {
    for(i=0;i<numnb;++i) {
      p_e[i]=0.0;p_LJ[i]=0.0; 
      if (flagf==2) {
	for(j=0;j<3; ++j) {
	  f_LJ[i*3+j]=0.0;f_e[i*3+j]=0.0;
	}
      }
    }
  }

  for(i=0;i<numnb;++i){
    num_a_prot=indexnb[i*2];
    NUM_A_PROT=indexnb[i*2+1];
    len2 = 0.0;
    for(j=0;j<3;++j){
      vec[j] = cord[NUM_A_PROT*3+j]-cord[num_a_prot*3+j];
      len2 += vec[j]*vec[j];
    }
    len = sqrt(len2);
    potedummy=ele[num_a_prot]*ele[NUM_A_PROT]/(len);
    potedummy=potedummy*u_1_5[i*2]*fact;
    if (flagp==1) {
      p_e[num_a_prot] += potedummy;
      p_e[NUM_A_PROT] += potedummy;
    }
    else if (flagp==2)
      p_e[i] = potedummy;
    if (flagf==1) {
      for(j=0; j<3; ++j) {
	fdummy[j] = -potedummy*UNIT*vec[j]/len2;
	f_e[num_a_prot*3+j] += fdummy[j];
	f_e[NUM_A_PROT*3+j] -= fdummy[j];
      }
    }
    else if (flagf==2) {
      for(j=0; j<3; ++j) {
	fdummy[j] = -potedummy*UNIT*vec[j]/len2;
	f_e[i*3+j] += fdummy[j];
      }
    }
    len6=len2;
    for (j=0;j<2;++j) {
      len6 = len6*len2;
    }
    p12 = ALJ[num_a_prot*num_atom+NUM_A_PROT]/(len6*len6);
    p6  = BLJ[num_a_prot*num_atom+NUM_A_PROT]/len6;
    p12 = p12*u_1_5[i*2+1]*fact;
    p6  = p6*u_1_5[i*2+1]*fact;
    if (flagp==1) {
      p_LJ[num_a_prot] += p12-p6;
      p_LJ[NUM_A_PROT] += p12-p6;
    }
    else if (flagp==2)
      p_LJ[i] = p12-p6;

    if (flagf==1) {
      for(j=0;j<3;++j) {
	fdummy[j]=-6.0*(2.0*p12-p6)/(len2)*vec[j]*UNIT;
	f_LJ[num_a_prot*3+j] +=fdummy[j];
	f_LJ[NUM_A_PROT*3+j] -= fdummy[j];
      }
    }
    else if (flagf==2) {
      for(j=0;j<3;++j) {
	fdummy[j]=-6.0*(2.0*p12-p6)/(len2)*vec[j]*UNIT;
	f_LJ[i*3+j] +=fdummy[j];
      }
    }
  }

  return 0;
}


void ffL_set_NB_PARM(double *ele, double *ALJ, double *BLJ, int numatom){
  int i,j,index;
  
  for (i=0;i<numatom;++i) {
    ele[i]= AP.CHRG[i];
    for (j=0;j<numatom;++j) {
      index=AP.ICO[(AP.IAC[i]-1)*AP.NTYPES+(AP.IAC[j]-1)]-1;
      ALJ[i*numatom+j] = AP.CN1[index];
      BLJ[i*numatom+j] = AP.CN2[index];
    }
  }
}

int ffL_set_NB_index(int *indexnb,int numnb, int numatom) {
  int i,j,k,l,m,n;
  int flag;
  FILE *log;

  l=0;n=0;
  for(i=0;i<numatom;++i){
    for (j=i+1;j<numatom;++j) {
      flag=1;
      for (k=0;k<AP.NUMEX[i];++k) {
	if (AP.NATEX[(n+k)]-1==j) {
	  flag=0;
	  break;
	}
      }
      if (flag==1) {
	indexnb[l*2]=i;
	indexnb[l*2+1]=j;
	++l;
      }
    }
    n+=AP.NUMEX[i];
  }

  if((log=fopen("log.txt","w"))==NULL)
    printf("error\n");
  for(i=0;i<numnb;++i){
    fprintf(log,"%d - %d \n",indexnb[i*2]+1,indexnb[i*2+1]+1);
  }
  fprintf(log,"\n  ");
  fclose(log);

  return l;
}

int ffL_set_numnb(void) {
  int i,numnb;

  numnb=0;
  for(i=0;i<AP.NATOM;++i)
    numnb+=AP.NATOM-i-1-AP.NUMEX[i];

  return numnb+1;
}

int ffL_calcDIHE(double *p_d,
		double *n_d,
		double *cord,
		int flagp, int flagf, int flaginp) {
  int i,j,k,l,flag;
  int dtype;
  
  double atom[4][3];
  double dihedang;
  
  for (i=0;i<AP.NPHIH;++i){
    p_d[i] = 0.0;
    if (flagf==1)
      n_d[i] = 0.0;
  }
  
  for (i=0;i<AP.MPHIA;++i){
    p_d[i+AP.NPHIH] = 0.0;
    if (flagf==1)
      n_d[i+AP.NPHIH] = 0.0;
  }
  
  for (i=0;i<AP.NPHIH;++i) {
    dtype = AP.PH[i][4]-1;
    flag=1;
    /************************************/
    /* for (j=0;j<inpnum;++j) {	        */
    /*   if (i == inpindex[j]) {        */
    /* 	flag=0;break;		        */
    /*   }			        */
    /* }			        */
    /************************************/
    if (flag==1 || flaginp==0) {
      for (j=0;j<4;++j) {
  	for (k=0;k<3;++k) {
  	  atom[j][k]=cord[abs(AP.PH[i][j])+k];
  	}
      }
  
      dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
      p_d[i] = AP.PK[dtype]*(1.0+cos(AP.PN[dtype]*dihedang-AP.PHASE[dtype]));
      if (flagf==1)
  	n_d[i] = -AP.PK[dtype]*4.18407*100.0*(sin(AP.PN[dtype]*dihedang-AP.PHASE[dtype])*AP.PN[dtype]);
    }
  }
  
  for (i=0;i<AP.MPHIA;++i) {
    dtype = AP.PA[i][4]-1;
    flag=1;
    /************************************/
    /* for (j=0;j<inpnum;++j) {	        */
    /*   if (i == inpindex[j]) {        */
    /* 	flag=0;break;		        */
    /*   }			        */
    /* }			        */
    /************************************/
    if (flag==1 || flaginp==0) {
      for (j=0;j<4;++j) {
  	for (k=0;k<3;++k) {
  	  atom[j][k]=cord[abs(AP.PA[i][j])+k];
  	}
      }
  
      dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
      p_d[i+AP.NPHIH] = AP.PK[dtype]*(1.0+cos(AP.PN[dtype]*dihedang-AP.PHASE[dtype]));
      if (flagf==1)
  	n_d[i+AP.NPHIH] = -AP.PK[dtype]*4.18407*100.0*(sin(AP.PN[dtype]*dihedang-AP.PHASE[dtype])*AP.PN[dtype]);
    }
  }
}

int ffL_calcDIHE_force_Cartesian(double *f_d,double *cord) {
  int i,j,k,l,flag;
  int dtype;

  double fa,fb[3],fc[3];
  double *n1,*n2,ln1,ln2;
  double vij[3],vkj[3],vki[3],vjl[3],vji[3],vik[3],vkl[3],vlk[3];
  double op1[3],op2[3],op3[3],op4[3],op5[3],op6[3];
  
  double atom[4][3];
  double cosdih,sindih;
  double dihedang;

  n1=(double *)gcemalloc(sizeof(double)*3);
  n2=(double *)gcemalloc(sizeof(double)*3);

  //  for (i=0;i<(AP.NPHIH+AP.MPHIA)*3;++i) f_d[i] = 0.0;
  for (i=0;i<AP.NATOM*3;++i) f_d[i] = 0.0;
  
  for (i=0;i<AP.NPHIH;++i) {
    dtype = AP.PH[i][4]-1;
    for (j=0;j<4;++j) for (k=0;k<3;++k) atom[j][k]=cord[abs(AP.PH[i][j])+k];

    for (j=0;j<3;++j) {
      vij[j] = atom[1][j]-atom[0][j];
      vkj[j] = atom[1][j]-atom[2][j];
      vki[j] = atom[0][j]-atom[2][j];
      vjl[j] = atom[3][j]-atom[1][j];
      vji[j] = atom[0][j]-atom[1][j];
      vik[j] = atom[2][j]-atom[0][j];
      vkl[j] = atom[3][j]-atom[2][j];
      vlk[j] = atom[2][j]-atom[3][j];
    }

    outprod(vij,vkj,n1);
    outprod(vkj,vkl,n2);
    ln1=sqrt(inprod(n1,n1,3));
    ln2=sqrt(inprod(n2,n2,3));
  
    csdih(atom[0],atom[1],atom[2],atom[3],&cosdih,&sindih);

    if (AP.PN[dtype]==1) fa=-AP.PK[dtype]*AP.PN[dtype]*cos(AP.PHASE[dtype]);
    else if (AP.PN[dtype]==2) fa=-AP.PK[dtype]*AP.PN[dtype]*2.0*cosdih*cos(AP.PHASE[dtype]);
    else if (AP.PN[dtype]==3) fa=-AP.PK[dtype]*AP.PN[dtype]*(-4.0*sindih*sindih+3.0)*cos(AP.PHASE[dtype]);
    else if (AP.PN[dtype]==4) fa=-AP.PK[dtype]*AP.PN[dtype]*4.0*(cosdih*(2.0*cosdih*cosdih-1.0))*cos(AP.PHASE[dtype]);
    else {
      printf("error:phase must be 1~4\n");
      exit(1);
    }

    for (j=0;j<3;++j) fb[j]=(n2[j]/ln2-cosdih*n1[j]/ln1)/ln1;
    for (j=0;j<3;++j) fc[j]=(n1[j]/ln1-cosdih*n2[j]/ln2)/ln2;

    outprod(fb,vkj,op1);
    //    outprod(fc,vki,op2); // 1111
    //    outprod(fb,vik,op3); // 1111
    outprod(fb,vki,op2);       // 1111
    outprod(fc,vlk,op3);       // 1111
    outprod(fb,vij,op4);
    outprod(fc,vjl,op5);
    outprod(fc,vkj,op6);

    for (j=0;j<3;++j) {
      f_d[abs(AP.PH[i][0])+j] += fa*op1[j]*UNIT;
      f_d[abs(AP.PH[i][1])+j] += fa*(-op2[j]+op3[j])*UNIT;
      f_d[abs(AP.PH[i][2])+j] += fa*(-op4[j]+op5[j])*UNIT;
      f_d[abs(AP.PH[i][3])+j] += fa*op6[j]*UNIT;
    }
  }
  
  for (i=0;i<AP.MPHIA;++i) {
    dtype = AP.PA[i][4]-1;
    for (j=0;j<4;++j) for (k=0;k<3;++k)	atom[j][k]=cord[abs(AP.PA[i][j])+k];

    for (j=0;j<3;++j) {
      vij[j] = atom[1][j]-atom[0][j];
      vkj[j] = atom[1][j]-atom[2][j];
      vki[j] = atom[0][j]-atom[2][j];
      vjl[j] = atom[3][j]-atom[1][j];
      vji[j] = atom[0][j]-atom[1][j];
      vik[j] = atom[2][j]-atom[0][j];
      vkl[j] = atom[3][j]-atom[2][j];
      vlk[j] = atom[2][j]-atom[3][j];
    }

    outprod(vij,vkj,n1);
    outprod(vkj,vkl,n2);
    ln1=sqrt(inprod(n1,n1,3));
    ln2=sqrt(inprod(n2,n2,3));
  
    csdih(atom[0],atom[1],atom[2],atom[3],&cosdih,&sindih);

    if (AP.PN[dtype]==1) fa=-AP.PK[dtype]*AP.PN[dtype]*cos(AP.PHASE[dtype]);
    else if (AP.PN[dtype]==2) fa=-AP.PK[dtype]*AP.PN[dtype]*2.0*cosdih*cos(AP.PHASE[dtype]);
    else if (AP.PN[dtype]==3) fa=-AP.PK[dtype]*AP.PN[dtype]*(-4.0*sindih*sindih+3.0)*cos(AP.PHASE[dtype]);
    else if (AP.PN[dtype]==4) fa=-AP.PK[dtype]*AP.PN[dtype]*4.0*(cosdih*(2.0*cosdih*cosdih-1.0))*cos(AP.PHASE[dtype]);
    else {
      printf("error:periodicity must be 1~4\n");
      exit(1);
    }

    for (j=0;j<3;++j) fb[j]=(n2[j]/ln2-cosdih*n1[j]/ln1)/ln1;
    for (j=0;j<3;++j) fc[j]=(n1[j]/ln1-cosdih*n2[j]/ln2)/ln2;

    outprod(fb,vkj,op1);
    //    outprod(fc,vki,op2); // 1111
    //    outprod(fb,vik,op3); // 1111
    outprod(fb,vki,op2);       // 1111
    outprod(fc,vlk,op3);       // 1111
    outprod(fb,vij,op4);
    outprod(fc,vjl,op5);
    outprod(fc,vkj,op6);

    for (j=0;j<3;++j) {
      f_d[abs(AP.PA[i][0])+j] += fa*op1[j]*UNIT;
      f_d[abs(AP.PA[i][1])+j] += fa*(-op2[j]+op3[j])*UNIT;
      f_d[abs(AP.PA[i][2])+j] += fa*(-op4[j]+op5[j])*UNIT;
      f_d[abs(AP.PA[i][3])+j] += fa*op6[j]*UNIT;
    }
    // for debug
    //    FILE *db;
    //    db=efopen("db_z_fa","a");
    //    fprintf(db,"%e\n",fa);
    //    fclose(db);
    //    FILE *db2;
    //    db2=efopen("db_z_sd_cd","a");
    //    fprintf(db2,"%e %e\n",sindih,cosdih);
    //    fclose(db2);
    // for debug
  }
}

int ffL_calc_spe_type_DIHE(double *p,
			  double *n,
			  double angle,
			  int type) {

  *p = AP.PK[type]*(1.0+cos(AP.PN[type]*angle-AP.PHASE[type]));
  *n = AP.PK[type]*4.18407*100.0*(sin(AP.PN[type]*angle-AP.PHASE[type])*AP.PN[type]);
  
}


int ffL_calcDIHE_wpep(double *p_d,
		     double *n_d,
		     double *cord,
		     double *u_d,
		     int flagp, int flagf, int flaginp) {
  int i,j,k,l,flag;
  int dtype;
  
  double atom[4][3];
  double dihedang;
  
  for (i=0;i<AP.NPHIH;++i){
    p_d[i] = 0.0;
    if (flagf==1)
      n_d[i] = 0.0;
  }
  
  for (i=0;i<AP.MPHIA;++i){
    p_d[i+AP.NPHIH] = 0.0;
    if (flagf==1)
      n_d[i+AP.NPHIH] = 0.0;
  }
  
  for (i=0;i<AP.NPHIH;++i) {
    dtype = AP.PH[i][4]-1;
    flag=1;
    /************************************/
    /* for (j=0;j<inpnum;++j) {	        */
    /*   if (i == inpindex[j]) {        */
    /* 	flag=0;break;		        */
    /*   }			        */
    /* }			        */
    /************************************/
    if (flag==1 || flaginp==0) {
      for (j=0;j<4;++j) {
  	for (k=0;k<3;++k) {
  	  atom[j][k]=cord[abs(AP.PH[i][j])+k];
  	}
      }
  
     dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
      p_d[i] = AP.PK[dtype]*(1.0+cos(AP.PN[dtype]*dihedang-AP.PHASE[dtype]))*u_d[i];
      if (flagf==1)
  	n_d[i] = -AP.PK[dtype]*4.18407*100.0*(sin(AP.PN[dtype]*dihedang-AP.PHASE[dtype])*AP.PN[dtype])*u_d[i];
    }
  }
  
  for (i=0;i<AP.MPHIA;++i) {
    dtype = AP.PA[i][4]-1;
    flag=1;
    /************************************/
    /* for (j=0;j<inpnum;++j) {	        */
    /*   if (i == inpindex[j]) {        */
    /* 	flag=0;break;		        */
    /*   }			        */
    /* }			        */
    /************************************/
    if (flag==1 || flaginp==0) {
      for (j=0;j<4;++j) {
  	for (k=0;k<3;++k) {
  	  atom[j][k]=cord[abs(AP.PA[i][j])+k];
  	}
      }
  
      dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
      p_d[i] = AP.PK[dtype]*(1.0+cos(AP.PN[dtype]*dihedang-AP.PHASE[dtype]))*u_d[i+AP.NPHIH];
      if (flagf==1)
  	n_d[i] = -AP.PK[dtype]*4.18407*100.0*(sin(AP.PN[dtype]*dihedang-AP.PHASE[dtype])*AP.PN[dtype])*u_d[i+AP.NPHIH];
    }
  }
}

int ffL_calcANGLE(double *p_a,double *cord){
  int i,j,k,l;
  int type;
  
  double atom[3][3];
  double ang;
  
  for (i=0;i<AP.NTHETH+AP.MTHETA;++i)
    p_a[i] = 0.0;
  
  for (i=0;i<AP.NTHETH;++i) {
    type = AP.TH[i][3]-1;
     for (j=0;j<3;++j)
      for (k=0;k<3;++k)
	atom[j][k]=cord[abs(AP.TH[i][j])+k];
  
    ang = pick_angle(atom[0],atom[1],atom[2],0,0.0);
    p_a[i] = AP.TK[type]*(ang-AP.TEQ[type])*(ang-AP.TEQ[type]);
  }
  
  for (i=0;i<AP.MTHETA;++i) {
    type = AP.TA[i][3]-1;
     for (j=0;j<3;++j)
      for (k=0;k<3;++k)
	atom[j][k]=cord[abs(AP.TA[i][j])+k];
  
    ang = pick_angle(atom[0],atom[1],atom[2],0,0.0);
    p_a[i+AP.NTHETH] = AP.TK[type]*(ang-AP.TEQ[type])*(ang-AP.TEQ[type]);
  }

  return 0;
}

int ffL_calcANGLE_force_Cartesian(double *f_a,double *cord){
  int i,j,k,l;
  int numatom;
  int type;
  double kang,ang_eq;

  double atom[3][3];
  double *f_temp;

  numatom=AP.NATOM;
  for (i=0;i<numatom*3;++i) f_a[i] = 0.0;
  f_temp=(double *)gcemalloc(sizeof(double)*3*3);

  for (i=0;i<AP.NTHETH;++i) {
    for (j=0;j<3;++j)  for (k=0;k<3;++k) atom[j][k]=cord[abs(AP.TH[i][j])+k];
    type = AP.TH[i][3]-1;
    kang = AP.TK[type];
    ang_eq = AP.TEQ[type];
    calcANGKE_force(atom[0],atom[1],atom[2],kang,ang_eq,f_temp);
    for (j=0;j<3;++j) {
      f_a[abs(AP.TH[i][0])+j] += f_temp[j];
      f_a[abs(AP.TH[i][1])+j] += f_temp[3+j];
      f_a[abs(AP.TH[i][2])+j] += f_temp[6+j];
    }
  }

  for (i=0;i<AP.MTHETA;++i) {
    for (j=0;j<3;++j)  for (k=0;k<3;++k) atom[j][k]=cord[abs(AP.TA[i][j])+k];
    type = AP.TA[i][3]-1;
    kang = AP.TK[type];
    ang_eq = AP.TEQ[type];
    calcANGKE_force(atom[0],atom[1],atom[2],kang,ang_eq,f_temp);
    for (j=0;j<3;++j) {
      f_a[abs(AP.TA[i][0])+j] += f_temp[j];
      f_a[abs(AP.TA[i][1])+j] += f_temp[3+j];
      f_a[abs(AP.TA[i][2])+j] += f_temp[6+j];
    }
  }

  return 0;
}

double calcANGKE_force(double atomi[3],double atomj[3],double atomk[3],double kang,double ang_eq,double *f) {
  int i,j,k;
  double lenij,lenkj;
  double vij[3],vkj[3];
  double cosijk,angijk;
  double f1,f2;

  lenij = len(atomi,atomj);
  lenkj = len(atomk,atomj);
  for (j=0;j<3;++j) {
    vij[j]=atomj[j]-atomi[j];
    vkj[j]=atomj[j]-atomk[j];
  }
  cosijk=inprod(vij,vkj,3);
  cosijk=cosijk/lenij/lenkj;
  angijk = acos(cosijk);
  //  angijk=ang(atomi,atomj,atomk);
  //  cosijk=cos(angijk);

  for (j=0;j<3;++j) {
    /*******************************************************************************************/
    /* f1 = -kang*(angijk-ang_eq)/(lenij*sin(angijk))*(vkj[j]/lenkj-cosijk*vij[j]/lenij)*UNIT; */
    /* f2 = -kang*(angijk-ang_eq)/(lenkj*sin(angijk))*(vij[j]/lenij-cosijk*vkj[j]/lenkj)*UNIT; */
    /*******************************************************************************************/
    f1 = -2.0*kang*(angijk-ang_eq)/(lenij*sin(angijk))*(vkj[j]/lenkj-cosijk*vij[j]/lenij)*UNIT;
    f2 = -2.0*kang*(angijk-ang_eq)/(lenkj*sin(angijk))*(vij[j]/lenij-cosijk*vkj[j]/lenkj)*UNIT;

    f[j] = f1;
    f[6+j] = f2;
    f[3+j] = -f1-f2;
  }
}

int ffL_calcBOND(double *p_b,double *cord){
  int i,j,k;
  int type;
  double len;
  double atom[2][3];

  for (i=0;i<AP.NBONH+AP.MBONA;++i)
    p_b[i] = 0.0;
  
  for (i=0;i<AP.NBONH;++i) {
    type = AP.BH[i][2]-1;
    for (j=0;j<2;++j)
      for (k=0;k<3;++k)
	atom[j][k]=cord[abs(AP.BH[i][j])+k];
  
    len = pick_bond_leng(atom[0],atom[1]);
    p_b[i] = AP.RK[type]*(len-AP.REQ[type])*(len-AP.REQ[type]);
  }

  for (i=0;i<AP.MBONA;++i) {
    type = AP.BA[i][2]-1;
    for (j=0;j<2;++j)
      for (k=0;k<3;++k)
	atom[j][k]=cord[abs(AP.BA[i][j])+k];  

    len = pick_bond_leng(atom[0],atom[1]);
    p_b[i+AP.NBONH] = AP.RK[type]*(len-AP.REQ[type])*(len-AP.REQ[type]);
  }

  return 0;
}

int ffL_calcBOND_force_Cartesian(double *f_b,double *cord){
  int i,j,k;
  int numatom;
  int type;
  double f;
  double lenij;
  double atom[2][3];

  numatom=AP.NATOM;
  for (i=0;i<numatom*3;++i) f_b[i] = 0.0;
  
  for (i=0;i<AP.NBONH;++i) {
    type = AP.BH[i][2]-1;
    for (j=0;j<2;++j) for (k=0;k<3;++k) atom[j][k]=cord[abs(AP.BH[i][j])+k];
  
    lenij = len(atom[0],atom[1]);
    for (j=0;j<3;++j) {
      f = -2.0*AP.RK[type]*(lenij-AP.REQ[type])*(atom[1][j]-atom[0][j])/lenij*UNIT;
      f_b[abs(AP.BH[i][0])+j] += f;
      f_b[abs(AP.BH[i][1])+j] += -f;
    }
  }

  for (i=0;i<AP.MBONA;++i) {
    type = AP.BA[i][2]-1;
    for (j=0;j<2;++j) for (k=0;k<3;++k) atom[j][k]=cord[abs(AP.BA[i][j])+k];
  
    lenij = len(atom[0],atom[1]);
    for (j=0;j<3;++j) {
      f = -2.0*AP.RK[type]*(lenij-AP.REQ[type])*(atom[1][j]-atom[0][j])/lenij*UNIT;
      f_b[abs(AP.BA[i][0])+j] += f;
      f_b[abs(AP.BA[i][1])+j] += -f;
    }
  }

  return 0;
}

int ffL_set_calcff(int numnb, int num14,FILE *inputfile, struct potential *ene){
  int i,numatom,f;
  char dummy;

  (*ene).parm.numnb=numnb;
  (*ene).parm.num14=num14;

  numatom=AP.NATOM;
  (*ene).parm.indexnb=(int *)gcemalloc(sizeof(int)*numnb*2);
  (*ene).parm.ele=(double *)gcemalloc(sizeof(double)*numatom);
  (*ene).parm.ALJ=(double *)gcemalloc(sizeof(double)*numatom*numatom);
  (*ene).parm.BLJ=(double *)gcemalloc(sizeof(double)*numatom*numatom);
  (*ene).parm.index14=(int *)gcemalloc(sizeof(int)*num14*2);
  (*ene).p_e=(double *)gcemalloc(sizeof(double)*numnb);
  (*ene).p_LJ=(double *)gcemalloc(sizeof(double)*numnb);
  (*ene).p_e_14=(double *)gcemalloc(sizeof(double)*numnb);
  (*ene).p_LJ_14=(double *)gcemalloc(sizeof(double)*numnb);
  (*ene).p_d=(double *)gcemalloc(sizeof(double)*(AP.NPHIH+AP.MPHIA));
  (*ene).p_a=(double *)gcemalloc(sizeof(double)*(AP.NTHETH+AP.MTHETA));
  (*ene).p_b=(double *)gcemalloc(sizeof(double)*(AP.NBONH+AP.MBONA));

  ffL_set_NB_PARM((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,numatom);

  for (i=0;i<numnb;++i) {
    fscanf(inputfile,"%d",&f);(*ene).parm.indexnb[i*2]=f-1;
    fscanf(inputfile,"%s",&dummy);
    fscanf(inputfile,"%d",&f);(*ene).parm.indexnb[i*2+1]=f-1;
  }
  for (i=0;i<num14;++i) {
    fscanf(inputfile,"%d",&f);(*ene).parm.index14[i*2]=f-1;
    fscanf(inputfile,"%s",&dummy);
    fscanf(inputfile,"%d",&f);(*ene).parm.index14[i*2+1]=f-1;
  }
}

int ffL_set_calcffsp(struct potential *ene){
  int i,numatom,fd;
  int numnb,num14;
  char dummy;

  ffL_set_non_bonding_index_1(&numnb,&num14);
  (*ene).parm.numnb=numnb;
  (*ene).parm.num14=num14;
  (*ene).parm.indexnb=(int *)gcemalloc(sizeof(int)*numnb*2);
  (*ene).parm.index14=(int *)gcemalloc(sizeof(int)*num14*2);
  ffL_set_non_bonding_index_2((*ene).parm.indexnb,(*ene).parm.index14);

  numatom=AP.NATOM;
  (*ene).parm.ele=(double *)gcemalloc(sizeof(double)*numatom);
  (*ene).parm.ALJ=(double *)gcemalloc(sizeof(double)*numatom*numatom);
  (*ene).parm.BLJ=(double *)gcemalloc(sizeof(double)*numatom*numatom);
  (*ene).p_e=(double *)gcemalloc(sizeof(double)*numnb);
  (*ene).p_LJ=(double *)gcemalloc(sizeof(double)*numnb);
  (*ene).p_e_14=(double *)gcemalloc(sizeof(double)*numnb);
  (*ene).p_LJ_14=(double *)gcemalloc(sizeof(double)*numnb);
  (*ene).p_d=(double *)gcemalloc(sizeof(double)*(AP.NPHIH+AP.MPHIA));
  (*ene).p_a=(double *)gcemalloc(sizeof(double)*(AP.NTHETH+AP.MTHETA));
  (*ene).p_b=(double *)gcemalloc(sizeof(double)*(AP.NBONH+AP.MBONA));

  ffL_set_NB_PARM((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,numatom);

}



int ffL_set_calcffandforce(struct potential *ene, struct force *f){
  int i,numatom,fd;
  int numnb,num14;
  char dummy;

  ffL_set_non_bonding_index_1(&numnb,&num14);
  (*ene).parm.numnb=numnb;
  (*ene).parm.num14=num14;
  (*ene).parm.indexnb=(int *)gcemalloc(sizeof(int)*numnb*2);
  (*ene).parm.index14=(int *)gcemalloc(sizeof(int)*num14*2);
  ffL_set_non_bonding_index_2((*ene).parm.indexnb,(*ene).parm.index14);

  numatom=AP.NATOM;
  (*ene).parm.ele=(double *)gcemalloc(sizeof(double)*numatom);
  (*ene).parm.ALJ=(double *)gcemalloc(sizeof(double)*numatom*numatom);
  (*ene).parm.BLJ=(double *)gcemalloc(sizeof(double)*numatom*numatom);
  (*ene).p_e=(double *)gcemalloc(sizeof(double)*numnb);
  (*ene).p_LJ=(double *)gcemalloc(sizeof(double)*numnb);
  (*ene).p_e_14=(double *)gcemalloc(sizeof(double)*numnb);
  (*ene).p_LJ_14=(double *)gcemalloc(sizeof(double)*numnb);
  (*ene).p_d=(double *)gcemalloc(sizeof(double)*(AP.NPHIH+AP.MPHIA));
  (*ene).p_a=(double *)gcemalloc(sizeof(double)*(AP.NTHETH+AP.MTHETA));
  (*ene).p_b=(double *)gcemalloc(sizeof(double)*(AP.NBONH+AP.MBONA));

  ffL_set_NB_PARM((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,numatom);

  (*f).f_t=(double *)gcemalloc(sizeof(double)*numatom*3);
  (*f).f_e=(double *)gcemalloc(sizeof(double)*numatom*3);
  (*f).f_LJ=(double *)gcemalloc(sizeof(double)*numatom*3);
  (*f).f_e_14=(double *)gcemalloc(sizeof(double)*numatom*3);
  (*f).f_LJ_14=(double *)gcemalloc(sizeof(double)*numatom*3);
  (*f).f_d=(double *)gcemalloc(sizeof(double)*numatom*3);
  (*f).f_a=(double *)gcemalloc(sizeof(double)*numatom*3);
  (*f).f_b=(double *)gcemalloc(sizeof(double)*numatom*3);

}

double ffL_calcff(double *crd, int numatom,struct potential *ene) {
  int i;
  int numnb,num14;
  double *f_e,*f_LJ,*n_d;
  double *f_e_14,*f_LJ_14;

  numnb=(*ene).parm.numnb;
  num14=(*ene).parm.num14;

  ffL_calcDIHE((*ene).p_d,n_d,crd,1,0,0);
  ffL_calcANGLE((*ene).p_a,crd);
  ffL_calcBOND((*ene).p_b,crd);


  ffL_calcFFNB((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e,(*ene).p_LJ,f_e,f_LJ,numnb,(*ene).parm.indexnb,numatom,crd,2,0);
  ffL_calcFFNB((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e_14,(*ene).p_LJ_14,f_e_14,f_LJ_14,num14,(*ene).parm.index14,numatom,crd,2,0);
  
  (*ene).p_t=0.0;
  (*ene).p_e_t=0.0;
  (*ene).p_LJ_t=0.0;
  (*ene).p_e_14_t=0.0;
  (*ene).p_LJ_14_t=0.0;
  (*ene).p_d_t=0.0;
  (*ene).p_a_t=0.0;
  (*ene).p_b_t=0.0;
    
  for (i=0;i<numnb;++i) {
    (*ene).p_t+=(*ene).p_e[i]+(*ene).p_LJ[i];
    (*ene).p_e_t+=(*ene).p_e[i];
    (*ene).p_LJ_t+=(*ene).p_LJ[i];
  }
  for (i=0;i<num14;++i) {
    (*ene).p_t+=1.0/1.2*(*ene).p_e_14[i]+0.5*(*ene).p_LJ_14[i];
    (*ene).p_e_14_t+=1.0/1.2*(*ene).p_e_14[i];
    (*ene).p_LJ_14_t+=0.5*(*ene).p_LJ_14[i];
  }
  for (i=0;i<AP.NPHIH+AP.MPHIA;++i) {
    (*ene).p_t+=(*ene).p_d[i];
    (*ene).p_d_t+=(*ene).p_d[i];
  }
  for (i=0;i<AP.NTHETH+AP.MTHETA;++i) {
    (*ene).p_t+=(*ene).p_a[i];
    (*ene).p_a_t+=(*ene).p_a[i];
  }
  for (i=0;i<AP.NBONH+AP.MBONA;++i) {
    (*ene).p_t+=(*ene).p_b[i];
    (*ene).p_b_t+=(*ene).p_b[i];
  }
  
  return (*ene).p_t;
}

double ffL_calcffandforce(double *crd, int numatom,struct potential *ene,struct force *f) {
  int i;
  int numnb,num14;
  double *n_d;

  numnb=(*ene).parm.numnb;
  num14=(*ene).parm.num14;

  ffL_calcFFNB((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e,(*ene).p_LJ,(*f).f_e,(*f).f_LJ,numnb,(*ene).parm.indexnb,numatom,crd,1,1);
  ffL_calcFFNB/*_14*/((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e_14,(*ene).p_LJ_14,(*f).f_e_14,(*f).f_LJ_14,num14,(*ene).parm.index14,numatom,crd,1,1);

  ffL_calcDIHE((*ene).p_d,n_d,crd,1,0,0);
  ffL_calcDIHE_force_Cartesian((*f).f_d,crd); // 1111

  ffL_calcANGLE((*ene).p_a,crd);
  ffL_calcANGLE_force_Cartesian((*f).f_a,crd);

  ffL_calcBOND((*ene).p_b,crd);
  ffL_calcBOND_force_Cartesian((*f).f_b,crd);
  
  (*ene).p_t=0.0;
  (*ene).p_e_t=0.0;
  (*ene).p_LJ_t=0.0;
  (*ene).p_e_14_t=0.0;
  (*ene).p_LJ_14_t=0.0;
  (*ene).p_d_t=0.0;
  (*ene).p_a_t=0.0;
  (*ene).p_b_t=0.0;
  for (i=0;i<numatom*3;++i) (*f).f_t[i]=0.0;

  for (i=0;i<numnb;++i) {
    (*ene).p_t+=(*ene).p_e[i]+(*ene).p_LJ[i];
    (*ene).p_e_t+=(*ene).p_e[i];
    (*ene).p_LJ_t+=(*ene).p_LJ[i];
  }
  for (i=0;i<num14;++i) {
    //    (*ene).p_t+=1.0/1.2*(*ene).p_e_14[i]+0.5*(*ene).p_LJ_14[i]; // 1011
    //    (*ene).p_e_14_t+=/*1.0/1.2**/(*ene).p_e_14[i];                  // 1011
    //    (*ene).p_LJ_14_t+=/*0.5**/(*ene).p_LJ_14[i];                    // 1011
    (*ene).p_e_14_t+=(*ene).p_e_14[i];             
    (*ene).p_LJ_14_t+=(*ene).p_LJ_14[i];               
  }
  (*ene).p_e_14_t=1.0/1.2*(*ene).p_e_14_t;
  (*ene).p_LJ_14_t=0.5*(*ene).p_LJ_14_t;
  (*ene).p_t=(*ene).p_e_14_t+(*ene).p_LJ_14_t;

  for (i=0;i<AP.NPHIH+AP.MPHIA;++i) {
    (*ene).p_t+=(*ene).p_d[i];
    (*ene).p_d_t+=(*ene).p_d[i];
  }
  for (i=0;i<AP.NTHETH+AP.MTHETA;++i) {
    (*ene).p_t+=(*ene).p_a[i];
    (*ene).p_a_t+=(*ene).p_a[i];             // 0911
  }                                          // 0911
  for (i=0;i<AP.NBONH+AP.MBONA;++i) {        // 0911
    (*ene).p_t+=(*ene).p_b[i];               // 0911
    (*ene).p_b_t+=(*ene).p_b[i];
  }
  for (i=0;i<numatom*3;++i) {
    (*f).f_e_14[i]=1.0/1.2*(*f).f_e_14[i];
    (*f).f_LJ_14[i]=0.5*(*f).f_LJ_14[i];
  }

  for (i=0;i<numatom*3;++i)
    (*f).f_t[i]+=(*f).f_e[i]+(*f).f_LJ[i]+(*f).f_e_14[i]+(*f).f_LJ_14[i]+(*f).f_d[i]+(*f).f_a[i]+(*f).f_b[i];
  
  return (*ene).p_t;
}

double ffL_calcffandforce_w14tune(double *crd, int numatom,struct potential *ene,struct force *f, int *atom_tune_pairs, double *tune_val, int numtune) {
  int i;
  int numnb,num14;
  double *n_d;

  numnb=(*ene).parm.numnb;
  num14=(*ene).parm.num14;

  ffL_calcFFNB((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e,(*ene).p_LJ,(*f).f_e,(*f).f_LJ,numnb,(*ene).parm.indexnb,numatom,crd,1,1);
  ffL_calcFFNB_wtune/*_14*/((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e_14,(*ene).p_LJ_14,(*f).f_e_14,(*f).f_LJ_14,num14,(*ene).parm.index14,numatom,crd,1,1,atom_tune_pairs,tune_val,numtune);

  ffL_calcDIHE((*ene).p_d,n_d,crd,1,0,0);
  ffL_calcDIHE_force_Cartesian((*f).f_d,crd); // 1111

  ffL_calcANGLE((*ene).p_a,crd);
  ffL_calcANGLE_force_Cartesian((*f).f_a,crd);

  ffL_calcBOND((*ene).p_b,crd);
  ffL_calcBOND_force_Cartesian((*f).f_b,crd);
  
  (*ene).p_t=0.0;
  (*ene).p_e_t=0.0;
  (*ene).p_LJ_t=0.0;
  (*ene).p_e_14_t=0.0;
  (*ene).p_LJ_14_t=0.0;
  (*ene).p_d_t=0.0;
  (*ene).p_a_t=0.0;
  (*ene).p_b_t=0.0;
  for (i=0;i<numatom*3;++i) (*f).f_t[i]=0.0;

  for (i=0;i<numnb;++i) {
    (*ene).p_t+=(*ene).p_e[i]+(*ene).p_LJ[i];
    (*ene).p_e_t+=(*ene).p_e[i];
    (*ene).p_LJ_t+=(*ene).p_LJ[i];
  }
  for (i=0;i<num14;++i) {
    //    (*ene).p_t+=1.0/1.2*(*ene).p_e_14[i]+0.5*(*ene).p_LJ_14[i]; // 1011
    //    (*ene).p_e_14_t+=/*1.0/1.2**/(*ene).p_e_14[i];                  // 1011
    //    (*ene).p_LJ_14_t+=/*0.5**/(*ene).p_LJ_14[i];                    // 1011
    (*ene).p_e_14_t+=(*ene).p_e_14[i];             
    (*ene).p_LJ_14_t+=(*ene).p_LJ_14[i];               
  }
  (*ene).p_e_14_t=1.0/1.2*(*ene).p_e_14_t;
  (*ene).p_LJ_14_t=0.5*(*ene).p_LJ_14_t;
  (*ene).p_t=(*ene).p_e_14_t+(*ene).p_LJ_14_t;

  for (i=0;i<AP.NPHIH+AP.MPHIA;++i) {
    (*ene).p_t+=(*ene).p_d[i];
    (*ene).p_d_t+=(*ene).p_d[i];
  }
  for (i=0;i<AP.NTHETH+AP.MTHETA;++i) {
    (*ene).p_t+=(*ene).p_a[i];
    (*ene).p_a_t+=(*ene).p_a[i];             // 0911
  }                                          // 0911
  for (i=0;i<AP.NBONH+AP.MBONA;++i) {        // 0911
    (*ene).p_t+=(*ene).p_b[i];               // 0911
    (*ene).p_b_t+=(*ene).p_b[i];
  }
  for (i=0;i<numatom*3;++i) {
    (*f).f_e_14[i]=1.0/1.2*(*f).f_e_14[i];
    (*f).f_LJ_14[i]=0.5*(*f).f_LJ_14[i];
  }

  for (i=0;i<numatom*3;++i)
    (*f).f_t[i]+=(*f).f_e[i]+(*f).f_LJ[i]+(*f).f_e_14[i]+(*f).f_LJ_14[i]+(*f).f_d[i]+(*f).f_a[i]+(*f).f_b[i];
  
  return (*ene).p_t;
}

double ffL_calcffandforce_w14tuneb(double *crd, int numatom,struct potential *ene,struct force *f, int *atom_tune_pairs_es, int *atom_tune_pairs_LJ, double *tune_val_es,double *tune_val_LJ, int numtune_es, int numtune_LJ) {
  int i;
  int numnb,num14;
  double *n_d;

  numnb=(*ene).parm.numnb;
  num14=(*ene).parm.num14;

  ffL_calcFFNB((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e,(*ene).p_LJ,(*f).f_e,(*f).f_LJ,numnb,(*ene).parm.indexnb,numatom,crd,1,1);
  ffL_calcFFNB_wtuneb((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e_14,(*ene).p_LJ_14,(*f).f_e_14,(*f).f_LJ_14,num14,(*ene).parm.index14,numatom,crd,1,1,atom_tune_pairs_es,atom_tune_pairs_LJ,tune_val_es,tune_val_LJ,numtune_es,numtune_LJ);

  ffL_calcDIHE((*ene).p_d,n_d,crd,1,0,0);
  ffL_calcDIHE_force_Cartesian((*f).f_d,crd); // 1111

  ffL_calcANGLE((*ene).p_a,crd);
  ffL_calcANGLE_force_Cartesian((*f).f_a,crd);

  ffL_calcBOND((*ene).p_b,crd);
  ffL_calcBOND_force_Cartesian((*f).f_b,crd);

  (*ene).p_t=0.0;
  (*ene).p_e_t=0.0;
  (*ene).p_LJ_t=0.0;
  (*ene).p_e_14_t=0.0;
  (*ene).p_LJ_14_t=0.0;
  (*ene).p_d_t=0.0;
  (*ene).p_a_t=0.0;
  (*ene).p_b_t=0.0;
  for (i=0;i<numatom*3;++i) (*f).f_t[i]=0.0;

  for (i=0;i<numnb;++i) {
    (*ene).p_t+=(*ene).p_e[i]+(*ene).p_LJ[i];
    (*ene).p_e_t+=(*ene).p_e[i];
    (*ene).p_LJ_t+=(*ene).p_LJ[i];
  }
  for (i=0;i<num14;++i) {
    //    (*ene).p_t+=1.0/1.2*(*ene).p_e_14[i]+0.5*(*ene).p_LJ_14[i]; // 1011
    //    (*ene).p_e_14_t+=/*1.0/1.2**/(*ene).p_e_14[i];                  // 1011
    //    (*ene).p_LJ_14_t+=/*0.5**/(*ene).p_LJ_14[i];                    // 1011
    (*ene).p_e_14_t+=(*ene).p_e_14[i];
    (*ene).p_LJ_14_t+=(*ene).p_LJ_14[i];
  }
  (*ene).p_e_14_t=1.0/1.2*(*ene).p_e_14_t;
  (*ene).p_LJ_14_t=0.5*(*ene).p_LJ_14_t;
  (*ene).p_t=(*ene).p_e_14_t+(*ene).p_LJ_14_t;

  for (i=0;i<AP.NPHIH+AP.MPHIA;++i) {
    (*ene).p_t+=(*ene).p_d[i];
    (*ene).p_d_t+=(*ene).p_d[i];
  }
  for (i=0;i<AP.NTHETH+AP.MTHETA;++i) {
    (*ene).p_t+=(*ene).p_a[i];
    (*ene).p_a_t+=(*ene).p_a[i];             // 0911
  }                                          // 0911
  for (i=0;i<AP.NBONH+AP.MBONA;++i) {        // 0911
    (*ene).p_t+=(*ene).p_b[i];               // 0911
    (*ene).p_b_t+=(*ene).p_b[i];
  }
  for (i=0;i<numatom*3;++i) {
    (*f).f_e_14[i]=1.0/1.2*(*f).f_e_14[i];
    (*f).f_LJ_14[i]=0.5*(*f).f_LJ_14[i];
  }

  for (i=0;i<numatom*3;++i)
    (*f).f_t[i]+=(*f).f_e[i]+(*f).f_LJ[i]+(*f).f_e_14[i]+(*f).f_LJ_14[i]+(*f).f_d[i]+(*f).f_a[i]+(*f).f_b[i];

  return (*ene).p_t;
}

double ffL_calcffandforce_14D_woH(double *crd, int numatom,struct potential *ene,struct force *f) {
  int i;
  int num14;
  double *n_d;

  num14=(*ene).parm.num14;

  ffL_calcFFNB_woH((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e_14,(*ene).p_LJ_14,(*f).f_e_14,(*f).f_LJ_14,num14,(*ene).parm.index14,numatom,crd,1,1);

  ffL_calcDIHE_woH((*ene).p_d,n_d,crd,1,0,0);

  (*ene).p_t=0.0;
  (*ene).p_e_14_t=0.0;
  (*ene).p_LJ_14_t=0.0;
  (*ene).p_d_t=0.0;
  for (i=0;i<numatom*3;++i) (*f).f_t[i]=0.0;

  for (i=AP.NPHIH;i<AP.NPHIH+AP.MPHIA;++i) {
    (*ene).p_t+=(*ene).p_d[i];
    (*ene).p_d_t+=(*ene).p_d[i];
  }

  for (i=0;i<num14;++i) {
    (*ene).p_e_14_t+=(*ene).p_e_14[i];             
    (*ene).p_LJ_14_t+=(*ene).p_LJ_14[i];               
  }
  (*ene).p_e_14_t=1.0/1.2*(*ene).p_e_14_t;
  (*ene).p_LJ_14_t=0.5*(*ene).p_LJ_14_t;
  (*ene).p_t=(*ene).p_e_14_t+(*ene).p_LJ_14_t;

  for (i=0;i<numatom*3;++i) {
    (*f).f_e_14[i]=1.0/1.2*(*f).f_e_14[i];
    (*f).f_LJ_14[i]=0.5*(*f).f_LJ_14[i];
  }

  for (i=0;i<numatom*3;++i) (*f).f_t[i]+=(*f).f_e_14[i]+(*f).f_LJ_14[i];
  
  return (*ene).p_t;
}

double ffL_calcffandforce_14DAB_woH(double *crd, int numatom,struct potential *ene,struct force *f) {
  int i;
  int num14;
  double *n_d;

  num14=(*ene).parm.num14;

  ffL_calcFFNB_woH((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e_14,(*ene).p_LJ_14,(*f).f_e_14,(*f).f_LJ_14,num14,(*ene).parm.index14,numatom,crd,1,1);

  ffL_calcDIHE_woH((*ene).p_d,n_d,crd,1,0,0);
  ffL_calcDIHE_force_Cartesian_woH((*f).f_d,crd);

  /*******************************************************/
  /* ffL_calcDIHE((*ene).p_d,n_d,crd,1,0,0);		 */
  /* ffL_calcDIHE_force_Cartesian((*f).f_d,crd); // 1111 */
  /*******************************************************/


  ffL_calcANGLE_woH((*ene).p_a,crd);
  ffL_calcANGLE_force_Cartesian_woH((*f).f_a,crd);

  ffL_calcBOND_woH((*ene).p_b,crd);
  ffL_calcBOND_force_Cartesian_woH((*f).f_b,crd);

  (*ene).p_t=0.0;
  (*ene).p_e_14_t=0.0;
  (*ene).p_LJ_14_t=0.0;
  (*ene).p_d_t=0.0;
  (*ene).p_a_t=0.0;
  (*ene).p_b_t=0.0;
  for (i=0;i<numatom*3;++i) (*f).f_t[i]=0.0;

  for (i=0;i<num14;++i) {
    (*ene).p_e_14_t+=(*ene).p_e_14[i];             
    (*ene).p_LJ_14_t+=(*ene).p_LJ_14[i];               
  }
  (*ene).p_e_14_t=1.0/1.2*(*ene).p_e_14_t;
  (*ene).p_LJ_14_t=0.5*(*ene).p_LJ_14_t;
  (*ene).p_t=(*ene).p_e_14_t+(*ene).p_LJ_14_t;

  for (i=0;i<numatom*3;++i) {
    (*f).f_e_14[i]=1.0/1.2*(*f).f_e_14[i];
    (*f).f_LJ_14[i]=0.5*(*f).f_LJ_14[i];
  }

  for (i=AP.NPHIH;i<AP.NPHIH+AP.MPHIA;++i) {
    (*ene).p_t+=(*ene).p_d[i];
    (*ene).p_d_t+=(*ene).p_d[i];
  }

  for (i=AP.NTHETH;i<AP.NTHETH+AP.MTHETA;++i) {
    (*ene).p_t+=(*ene).p_a[i];
    (*ene).p_a_t+=(*ene).p_a[i];
  }                                     

  for (i=AP.NBONH;i<AP.NBONH+AP.MBONA;++i) {    
    (*ene).p_t+=(*ene).p_b[i];
    (*ene).p_b_t+=(*ene).p_b[i];
  }

  for (i=0;i<numatom*3;++i) (*f).f_t[i]+=(*f).f_e_14[i]+(*f).f_LJ_14[i]+(*f).f_d[i]+(*f).f_a[i]+(*f).f_b[i];
  
  return (*ene).p_t;
}

double ffL_calcffandforce_14vdWDAB_woH(double *crd, int numatom,struct potential *ene,struct force *f) {
  int i;
  int num14;
  double *n_d;

  num14=(*ene).parm.num14;

  ffL_calcFFNBvdW_woH((*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_LJ_14,(*f).f_LJ_14,num14,(*ene).parm.index14,numatom,crd,1,1);

  ffL_calcDIHE_woH((*ene).p_d,n_d,crd,1,0,0);
  ffL_calcDIHE_force_Cartesian_woH((*f).f_d,crd);

  ffL_calcANGLE_woH((*ene).p_a,crd);
  ffL_calcANGLE_force_Cartesian_woH((*f).f_a,crd);

  ffL_calcBOND_woH((*ene).p_b,crd);
  ffL_calcBOND_force_Cartesian_woH((*f).f_b,crd);

  (*ene).p_t=0.0;
  (*ene).p_LJ_14_t=0.0;
  (*ene).p_d_t=0.0;
  (*ene).p_a_t=0.0;
  (*ene).p_b_t=0.0;
  for (i=0;i<numatom*3;++i) (*f).f_t[i]=0.0;

  for (i=0;i<num14;++i) {
    (*ene).p_LJ_14_t+=(*ene).p_LJ_14[i];               
  }
  (*ene).p_LJ_14_t=0.5*(*ene).p_LJ_14_t;
  (*ene).p_t=(*ene).p_LJ_14_t;

  for (i=0;i<numatom*3;++i) {
    (*f).f_LJ_14[i]=0.5*(*f).f_LJ_14[i];
  }

  for (i=AP.NPHIH;i<AP.NPHIH+AP.MPHIA;++i) {
    (*ene).p_t+=(*ene).p_d[i];
    (*ene).p_d_t+=(*ene).p_d[i];
  }

  for (i=AP.NTHETH;i<AP.NTHETH+AP.MTHETA;++i) {
    (*ene).p_t+=(*ene).p_a[i];
    (*ene).p_a_t+=(*ene).p_a[i];
  }                                     

  for (i=AP.NBONH;i<AP.NBONH+AP.MBONA;++i) {    
    (*ene).p_t+=(*ene).p_b[i];
    (*ene).p_b_t+=(*ene).p_b[i];
  }

  for (i=0;i<numatom*3;++i) (*f).f_t[i]+=(*f).f_e_14[i]+(*f).f_LJ_14[i]+(*f).f_d[i]+(*f).f_a[i]+(*f).f_b[i];
  
  return (*ene).p_t;
}



double ffL_calcffandforce_14_6(double *crd, int numatom,struct potential *ene,struct force *f,double scale,int numrepul) { // 0911
  int i;
  int numnb,num14;
  double *n_d;

  numnb=(*ene).parm.numnb;
  num14=(*ene).parm.num14;

  ffL_calcFFNB_14_6((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e,(*ene).p_LJ,(*f).f_e,(*f).f_LJ,numnb,(*ene).parm.indexnb,numatom,crd,1,1,scale,numrepul);
  ffL_calcFFNB_14_6((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e_14,(*ene).p_LJ_14,(*f).f_e_14,(*f).f_LJ_14,num14,(*ene).parm.index14,numatom,crd,1,1,scale,numrepul);

  ffL_calcDIHE((*ene).p_d,n_d,crd,1,0,0);
  ffL_calcDIHE_force_Cartesian((*f).f_d,crd);

  ffL_calcANGLE((*ene).p_a,crd);
  ffL_calcANGLE_force_Cartesian((*f).f_a,crd);

  ffL_calcBOND((*ene).p_b,crd);
  ffL_calcBOND_force_Cartesian((*f).f_b,crd);
  
  (*ene).p_t=0.0;
  (*ene).p_e_t=0.0;
  (*ene).p_LJ_t=0.0;
  (*ene).p_e_14_t=0.0;
  (*ene).p_LJ_14_t=0.0;
  (*ene).p_d_t=0.0;
  (*ene).p_a_t=0.0;
  (*ene).p_b_t=0.0;
  for (i=0;i<numatom*3;++i) (*f).f_t[i]=0.0;

  for (i=0;i<numnb;++i) {
    (*ene).p_t+=(*ene).p_e[i]+(*ene).p_LJ[i];
    (*ene).p_e_t+=(*ene).p_e[i];
    (*ene).p_LJ_t+=(*ene).p_LJ[i];
  }
  for (i=0;i<num14;++i) {
    (*ene).p_t+=1.0/1.2*(*ene).p_e_14[i]+0.5*(*ene).p_LJ_14[i];
    (*ene).p_e_14_t+=(*ene).p_e_14[i];
    (*ene).p_LJ_14_t+=(*ene).p_LJ_14[i];
  }
  for (i=0;i<AP.NPHIH+AP.MPHIA;++i) {
    (*ene).p_t+=(*ene).p_d[i];
    (*ene).p_d_t+=(*ene).p_d[i];
  }
  for (i=0;i<AP.NTHETH+AP.MTHETA;++i) {
    (*ene).p_t+=(*ene).p_a[i];
    (*ene).p_a_t+=(*ene).p_a[i];             
  }                                          
  for (i=0;i<AP.NBONH+AP.MBONA;++i) {        
    (*ene).p_t+=(*ene).p_b[i];               
    (*ene).p_b_t+=(*ene).p_b[i];
  }
  for (i=0;i<numatom*3;++i) {
    (*f).f_e_14[i]=1.0/1.2*(*f).f_e_14[i];
    (*f).f_LJ_14[i]=0.5*(*f).f_LJ_14[i];
  }
  for (i=0;i<numatom*3;++i)
    (*f).f_t[i]+=(*f).f_e[i]+(*f).f_LJ[i]+(*f).f_e_14[i]+(*f).f_LJ_14[i]+(*f).f_d[i]+(*f).f_a[i]+(*f).f_b[i];
  
  return (*ene).p_t;
}


void ffL_set_non_bonding_index_1(int *numindexnb, int *numindex14) {
  int i,j,k,l;
  int numatom;

  int gnumnb=0;
  int gnum14=0;

  numatom=AP.NATOM;

  gnumnb=0;
  gnum14=0;

  for(i=0;i<numatom;++i) {
    for(j=i+1;j<numatom;++j) {
      if (which_calc_nb(i,j) == YES)
	++gnumnb;
      else if(which_calc_14_nb(i,j)==YES)
      	++gnum14;
    }
  }

  (*numindexnb)=gnumnb;
  (*numindex14)=gnum14;

}

void ffL_set_non_bonding_index_2(int *gindexnb,int *gindex14) {
  int i,j,k,l;
  int numatom;

  int gnumnb=0;
  int gnum14=0;

  numatom=AP.NATOM;

  k=0;
  l=0;
  for(i=0;i<numatom;++i)
    for(j=i+1;j<numatom;++j)
      if (which_calc_nb(i,j) == YES) {
  	gindexnb[k*2]=i;
  	gindexnb[k*2+1]=j;
  	++k;
      }
      else if(which_calc_14_nb(i,j)==YES){
  	gindex14[l*2]=i;
  	gindex14[l*2+1]=j;
  	++l;
      }

}

int which_calc_nb(int atomi,int atomj) {
  int k;
  int sum=0;

  for(k=0;k<atomi;++k) sum+=AP.NUMEX[k];

  for(k=sum;k<sum+AP.NUMEX[atomi];++k){
    if (AP.NATEX[k]-1 == atomj)
      return NO;
  }

  return YES;
}

int which_calc_14_nb(int atomi,int atomj) {
  int k,l;
  int flag;

  if ((flag=which_calc_nb(atomi,atomj))==YES)
    return NO;
  else {
    for(l=0;l<AP.NBONH;++l){
      if ((AP.BH[l][0])/3 == atomi ){
  	if ((AP.BH[l][1])/3 == atomj)
  	  return NO;
      }
      else if ((AP.BH[l][1])/3 == atomi){
  	if ((AP.BH[l][0])/3 == atomj)
  	  return NO;
      }
      else if ((AP.BH[l][0])/3 == atomj){
  	if ((AP.BH[l][1])/3 == atomi)
  	  return NO;
      }
      else if ((AP.BH[l][1])/3 == atomj){
  	if ((AP.BH[l][0])/3 == atomi)
  	  return NO;
      }
    }
  
    for(l=0;l<AP.MBONA;++l){
      if ((AP.BA[l][0])/3 == atomi ){
  	if ((AP.BA[l][1])/3 == atomj)
  	  return NO;
      }
      else if ((AP.BA[l][1])/3 == atomi){
  	if ((AP.BA[l][0])/3 == atomj)
  	  return NO;
      }
      else if ((AP.BA[l][0])/3 == atomj){
  	if ((AP.BA[l][1])/3 == atomi)
  	  return NO;
      }
      else if ((AP.BA[l][1])/3 == atomj){
  	if ((AP.BA[l][0])/3 == atomi)
  	  return NO;
      }
    }
  
    for(l=0;l<AP.NTHETH;++l){
      if ((AP.TH[l][0])/3 == atomi ){
  	if ((AP.TH[l][2])/3 == atomj)
  	  return NO;
      }
      else if ((AP.TH[l][2])/3 == atomi){
  	if ((AP.TH[l][0])/3 == atomj)
  	  return NO;
      }
      else if ((AP.TH[l][0])/3 == atomj){
  	if ((AP.TH[l][2])/3 == atomi)
  	  return NO;
      }
      else if ((AP.TH[l][2])/3 == atomj){
  	if ((AP.TH[l][0])/3 == atomi)
  	  return NO;
      }
    }
  
    for(l=0;l<AP.MTHETA;++l){
      if ((AP.TA[l][0])/3 == atomi ) {
  	if ((AP.TA[l][2])/3 == atomj)
  	  return NO;
      }
      //      else if ((AP.TH[l][2])/3 == atomi){
      else if ((AP.TA[l][2])/3 == atomi){ //1011
  	if ((AP.TA[l][0])/3 == atomj)
  	  return NO;
      }
      else if ((AP.TA[l][0])/3 == atomj){
  	if ((AP.TA[l][2])/3 == atomi)
  	  return NO;
      }
      else if ((AP.TA[l][2])/3 == atomj){
  	if ((AP.TA[l][0])/3 == atomi)
  	  return NO;
      }
    }
  }
  
  return YES;
}


///////////////////////////////////////////////////////////////////

double ffL_calcffandforce_woH(double *crd, int numatom,struct potential *ene,struct force *f) {
  int i;
  int numnb,num14;
  double *n_d;

  numnb=(*ene).parm.numnb;
  num14=(*ene).parm.num14;

  ffL_calcFFNB_woH((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e,(*ene).p_LJ,(*f).f_e,(*f).f_LJ,numnb,(*ene).parm.indexnb,numatom,crd,1,1);
  ffL_calcFFNB_woH((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e_14,(*ene).p_LJ_14,(*f).f_e_14,(*f).f_LJ_14,num14,(*ene).parm.index14,numatom,crd,1,1);

  ffL_calcDIHE_woH((*ene).p_d,n_d,crd,1,0,0);
  ffL_calcDIHE_force_Cartesian_woH((*f).f_d,crd);

  ffL_calcANGLE_woH((*ene).p_a,crd);
  ffL_calcANGLE_force_Cartesian_woH((*f).f_a,crd);

  ffL_calcBOND_woH/*_wFC100*/((*ene).p_b,crd);
  ffL_calcBOND_force_Cartesian_woH/*_wFC100*/((*f).f_b,crd);
  
  (*ene).p_t=0.0;
  (*ene).p_e_t=0.0;
  (*ene).p_LJ_t=0.0;
  (*ene).p_e_14_t=0.0;
  (*ene).p_LJ_14_t=0.0;
  (*ene).p_d_t=0.0;
  (*ene).p_a_t=0.0;
  (*ene).p_b_t=0.0;
  for (i=0;i<numatom*3;++i) (*f).f_t[i]=0.0;

  for (i=0;i<numnb;++i) {
    (*ene).p_t+=(*ene).p_e[i]+(*ene).p_LJ[i];
    (*ene).p_e_t+=(*ene).p_e[i];
    (*ene).p_LJ_t+=(*ene).p_LJ[i];
  }
  for (i=0;i<num14;++i) {
    (*ene).p_t+=1.0/1.2*(*ene).p_e_14[i]+0.5*(*ene).p_LJ_14[i];
    (*ene).p_e_14_t+=(*ene).p_e_14[i];
    (*ene).p_LJ_14_t+=(*ene).p_LJ_14[i];
  }
  for (i=0;i<AP.NPHIH+AP.MPHIA;++i) {
    (*ene).p_t+=(*ene).p_d[i];
    (*ene).p_d_t+=(*ene).p_d[i];
  }
  for (i=0;i<AP.NTHETH+AP.MTHETA;++i) {
    (*ene).p_t+=(*ene).p_a[i];
    (*ene).p_a_t+=(*ene).p_a[i];
  }                                      // 0911
  for (i=0;i<AP.NBONH+AP.MBONA;++i) {    // 0911
    (*ene).p_t+=(*ene).p_b[i];
    (*ene).p_b_t+=(*ene).p_b[i];
  }
  for (i=0;i<numatom*3;++i) {
    (*f).f_e_14[i]=1.0/1.2*(*f).f_e_14[i];
    (*f).f_LJ_14[i]=0.5*(*f).f_LJ_14[i];
  }
  for (i=0;i<numatom*3;++i)
    (*f).f_t[i]+=(*f).f_e[i]+(*f).f_LJ[i]+(*f).f_e_14[i]+(*f).f_LJ_14[i]+(*f).f_d[i]+(*f).f_a[i]+(*f).f_b[i];
  
  return (*ene).p_t;
}

int ffL_calcDIHE_woH(double *p_d,
		double *n_d,
		double *cord,
		int flagp, int flagf, int flaginp) {
  int i,j,k,l,flag;
  int dtype;
  
  double atom[4][3];
  double dihedang;
  
  for (i=0;i<AP.NPHIH;++i){
    p_d[i] = 0.0;
    if (flagf==1)
      n_d[i] = 0.0;
  }
  
  for (i=0;i<AP.MPHIA;++i){
    p_d[i+AP.NPHIH] = 0.0;
    if (flagf==1)
      n_d[i+AP.NPHIH] = 0.0;
  }
  
  for (i=0;i<AP.MPHIA;++i) {
    dtype = AP.PA[i][4]-1;
    flag=1;
    if (flag==1 || flaginp==0) {
      for (j=0;j<4;++j) {
  	for (k=0;k<3;++k) {
  	  atom[j][k]=cord[abs(AP.PA[i][j])+k];
  	}
      }
  
      dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
      p_d[i+AP.NPHIH] = AP.PK[dtype]*(1.0+cos(AP.PN[dtype]*dihedang-AP.PHASE[dtype]));
      if (flagf==1)
  	n_d[i+AP.NPHIH] = -AP.PK[dtype]*4.18407*100.0*(sin(AP.PN[dtype]*dihedang-AP.PHASE[dtype])*AP.PN[dtype]);
    }
  }
}

int ffL_calcDIHE_force_Cartesian_woH(double *f_d,double *cord) {
  int i,j,k,l,flag;
  int dtype;

  double fa,fb[3],fc[3];
  double *n1,*n2,ln1,ln2;
  double vij[3],vkj[3],vki[3],vjl[3],vji[3],vik[3],vkl[3],vlk[3];
  double op1[3],op2[3],op3[3],op4[3],op5[3],op6[3];
  
  double atom[4][3];
  double cosdih,sindih;
  double dihedang;

  n1=(double *)gcemalloc(sizeof(double)*3);
  n2=(double *)gcemalloc(sizeof(double)*3);

  //  for (i=0;i<(AP.NPHIH+AP.MPHIA)*3;++i) f_d[i] = 0.0;
  for (i=0;i<AP.NATOM*3;++i) f_d[i] = 0.0;
  
  for (i=0;i<AP.MPHIA;++i) {
    dtype = AP.PA[i][4]-1;
    for (j=0;j<4;++j) for (k=0;k<3;++k)	atom[j][k]=cord[abs(AP.PA[i][j])+k];

    for (j=0;j<3;++j) {
      vij[j] = atom[1][j]-atom[0][j];
      vkj[j] = atom[1][j]-atom[2][j];
      vki[j] = atom[0][j]-atom[2][j];
      vjl[j] = atom[3][j]-atom[1][j];
      vji[j] = atom[0][j]-atom[1][j];
      vik[j] = atom[2][j]-atom[0][j];
      vkl[j] = atom[3][j]-atom[2][j];
      vlk[j] = atom[2][j]-atom[3][j];
    }

    outprod(vij,vkj,n1);
    outprod(vkj,vkl,n2);
    ln1=sqrt(inprod(n1,n1,3));
    ln2=sqrt(inprod(n2,n2,3));
  
    csdih(atom[0],atom[1],atom[2],atom[3],&cosdih,&sindih);

    if (AP.PN[dtype]==1) fa=-AP.PK[dtype]*AP.PN[dtype]*cos(AP.PHASE[dtype]);
    else if (AP.PN[dtype]==2) fa=-AP.PK[dtype]*AP.PN[dtype]*2.0*cosdih*cos(AP.PHASE[dtype]);
    else if (AP.PN[dtype]==3) fa=-AP.PK[dtype]*AP.PN[dtype]*(-4.0*sindih*sindih+3.0)*cos(AP.PHASE[dtype]);
    else if (AP.PN[dtype]==4) fa=-AP.PK[dtype]*AP.PN[dtype]*4.0*(cosdih*(2.0*cosdih*cosdih-1.0))*cos(AP.PHASE[dtype]);
    else {
      printf("error:periodicity must be 1~4\n");
      exit(1);
    }

    for (j=0;j<3;++j) fb[j]=(n2[j]/ln2-cosdih*n1[j]/ln1)/ln1;
    for (j=0;j<3;++j) fc[j]=(n1[j]/ln1-cosdih*n2[j]/ln2)/ln2;

    outprod(fb,vkj,op1);
    //    outprod(fc,vki,op2);  // 1111
    //    outprod(fb,vik,op3);  // 1111
    outprod(fb,vki,op2);        // 1111
    outprod(fc,vlk,op3);        // 1111
    outprod(fb,vij,op4);
    outprod(fc,vjl,op5);
    outprod(fc,vkj,op6);

    for (j=0;j<3;++j) {
      f_d[abs(AP.PA[i][0])+j] += fa*op1[j]*UNIT;
      f_d[abs(AP.PA[i][1])+j] += fa*(-op2[j]+op3[j])*UNIT;
      f_d[abs(AP.PA[i][2])+j] += fa*(-op4[j]+op5[j])*UNIT;
      f_d[abs(AP.PA[i][3])+j] += fa*op6[j]*UNIT;
    }
  }
}

int ffL_calcffandforce_speD(double *f_d,double *p_d,double *cord,int numspedihed,
			    int *atom1,int *atom2,int *atom3,int *atom4) {
  int i,j,k,l,flag;
  int dtype;

  double fa,fb[3],fc[3];
  double *n1,*n2,ln1,ln2;
  double vij[3],vkj[3],vki[3],vjl[3],vji[3],vik[3],vkl[3],vlk[3];
  double op1[3],op2[3],op3[3],op4[3],op5[3],op6[3];
  
  double atom[4][3];
  double cosdih,sindih;
  double dihedang;

  n1=(double *)gcemalloc(sizeof(double)*3);
  n2=(double *)gcemalloc(sizeof(double)*3);

  //  for (i=0;i<(AP.NPHIH+AP.MPHIA)*3;++i) f_d[i] = 0.0;
  for (i=0;i<4*3;++i) f_d[i] = 0.0;
  
  dtype = AP.PA[numspedihed][4]-1;
  for (j=0;j<4;++j) for (k=0;k<3;++k)	atom[j][k]=cord[abs(AP.PA[numspedihed][j])+k];
  
  for (j=0;j<3;++j) {
    vij[j] = atom[1][j]-atom[0][j];
    vkj[j] = atom[1][j]-atom[2][j];
    vki[j] = atom[0][j]-atom[2][j];
    vjl[j] = atom[3][j]-atom[1][j];
    vji[j] = atom[0][j]-atom[1][j];
    vik[j] = atom[2][j]-atom[0][j];
    vkl[j] = atom[3][j]-atom[2][j];
    vlk[j] = atom[2][j]-atom[3][j];
  }
  
  outprod(vij,vkj,n1);
  outprod(vkj,vkl,n2);
  ln1=sqrt(inprod(n1,n1,3));
  ln2=sqrt(inprod(n2,n2,3));
  
  csdih(atom[0],atom[1],atom[2],atom[3],&cosdih,&sindih);

  dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
  *p_d = AP.PK[dtype]*(1.0+cos(AP.PN[dtype]*dihedang-AP.PHASE[dtype]));

  if (AP.PN[dtype]==1) fa=-AP.PK[dtype]*AP.PN[dtype]*cos(AP.PHASE[dtype]);
  else if (AP.PN[dtype]==2) fa=-AP.PK[dtype]*AP.PN[dtype]*2.0*cosdih*cos(AP.PHASE[dtype]);
  else if (AP.PN[dtype]==3) fa=-AP.PK[dtype]*AP.PN[dtype]*(-4.0*sindih*sindih+3.0)*cos(AP.PHASE[dtype]);
  else if (AP.PN[dtype]==4) fa=-AP.PK[dtype]*AP.PN[dtype]*4.0*(cosdih*(2.0*cosdih*cosdih-1.0))*cos(AP.PHASE[dtype]);
  else {
    printf("error:periodicity must be 1~4\n");
    exit(1);
  }

  for (j=0;j<3;++j) fb[j]=(n2[j]/ln2-cosdih*n1[j]/ln1)/ln1;
  for (j=0;j<3;++j) fc[j]=(n1[j]/ln1-cosdih*n2[j]/ln2)/ln2;
    
  outprod(fb,vkj,op1);
  //  outprod(fc,vki,op2);
  //  outprod(fb,vik,op3);
  outprod(fb,vki,op2);
  outprod(fc,vlk,op3);
  outprod(fb,vij,op4);
  outprod(fc,vjl,op5);
  outprod(fc,vkj,op6);
  
  for (j=0;j<3;++j) {
    f_d[0*3+j] += fa*op1[j]*UNIT;
    f_d[1*3+j] += fa*(-op2[j]+op3[j])*UNIT;
    f_d[2*3+j] += fa*(-op4[j]+op5[j])*UNIT;
    f_d[3*3+j] += fa*op6[j]*UNIT;
  }
  *atom1=AP.PA[numspedihed][0];
  *atom2=AP.PA[numspedihed][1];
  *atom3=AP.PA[numspedihed][2];
  *atom4=AP.PA[numspedihed][3];

}

int ffL_calcANGLE_woH(double *p_a,double *cord){
  int i,j,k,l;
  int type;
  
  double atom[3][3];
  double ang;
  
  for (i=0;i<AP.NTHETH+AP.MTHETA;++i)
    p_a[i] = 0.0;
  
  for (i=0;i<AP.MTHETA;++i) {
    type = AP.TA[i][3]-1;
     for (j=0;j<3;++j)
      for (k=0;k<3;++k)
	atom[j][k]=cord[abs(AP.TA[i][j])+k];
  
    ang = pick_angle(atom[0],atom[1],atom[2],0,0.0);
    p_a[i+AP.NTHETH] = AP.TK[type]*(ang-AP.TEQ[type])*(ang-AP.TEQ[type]);
  }

  return 0;
}

int ffL_calcANGLE_force_Cartesian_woH(double *f_a,double *cord){
  int i,j,k,l;
  int numatom;
  int type;
  double kang,ang_eq;

  double atom[3][3];
  double *f_temp;

  numatom=AP.NATOM;
  for (i=0;i<numatom*3;++i) f_a[i] = 0.0;
  f_temp=(double *)gcemalloc(sizeof(double)*3*3);

  for (i=0;i<AP.MTHETA;++i) {
    for (j=0;j<3;++j)  for (k=0;k<3;++k) atom[j][k]=cord[abs(AP.TA[i][j])+k];
    type = AP.TA[i][3]-1;
    kang = AP.TK[type];
    ang_eq = AP.TEQ[type];
    calcANGKE_force(atom[0],atom[1],atom[2],kang,ang_eq,f_temp);
    for (j=0;j<3;++j) {
      f_a[abs(AP.TA[i][0])+j] += f_temp[j];
      f_a[abs(AP.TA[i][1])+j] += f_temp[3+j];
      f_a[abs(AP.TA[i][2])+j] += f_temp[6+j];
    }
  }

  return 0;
}

int ffL_calcBOND_woH(double *p_b,double *cord){
  int i,j,k;
  int type;
  double len;
  double atom[2][3];

  for (i=0;i<AP.NBONH+AP.MBONA;++i)
    p_b[i] = 0.0;
  
  for (i=0;i<AP.MBONA;++i) {
    type = AP.BA[i][2]-1;
    for (j=0;j<2;++j)
      for (k=0;k<3;++k)
	atom[j][k]=cord[abs(AP.BA[i][j])+k];  

    len = pick_bond_leng(atom[0],atom[1]);
    p_b[i+AP.NBONH] = AP.RK[type]*(len-AP.REQ[type])*(len-AP.REQ[type]);
  }

  return 0;
}

int ffL_calcBOND_woH_wFC100(double *p_b,double *cord){
  ;
}


int ffL_calcBOND_force_Cartesian_woH(double *f_b,double *cord){
  int i,j,k;
  int numatom;
  int type;
  double f;
  double lenij;
  double atom[2][3];

  numatom=AP.NATOM;
  for (i=0;i<numatom*3;++i) f_b[i] = 0.0;
  
  for (i=0;i<AP.MBONA;++i) {
    type = AP.BA[i][2]-1;
    for (j=0;j<2;++j) for (k=0;k<3;++k) atom[j][k]=cord[abs(AP.BA[i][j])+k];
  
    lenij = len(atom[0],atom[1]);
    for (j=0;j<3;++j) {
      f = -2.0*AP.RK[type]*(lenij-AP.REQ[type])*(atom[1][j]-atom[0][j])/lenij*UNIT;
      f_b[abs(AP.BA[i][0])+j] += f;
      f_b[abs(AP.BA[i][1])+j] += -f;
    }
  }

  return 0;
}

int ffL_calcBOND_force_Cartesian_woH_wFC100(double *f_b,double *cord){
  ;
}


int ffL_calcFFNB_woH(double *ele, double *ALJ, double *BLJ,
		    double *p_e,double *p_LJ,
		    double *f_e,double *f_LJ,
		    int numnb, int *indexnb,
		    int num_atom,
		    double *cord,
		    int flagp, int flagf) {
  int i,j;
  int num_a_prot,NUM_A_PROT;
  double potedummy;
  double f[3];
  double len,len2,len6;
  double vec[3],fdummy[3];
  double p12,p6;

  if (flagf != 0 && flagf != 1 && flagf !=2 ) {
    printf("error !\n");exit(1);
  }
  if (flagp != 0 && flagp != 1 && flagp !=2 ) {
    printf("error !\n");exit(1);
  }

  if (flagp==1 ) {
    for(i=0;i<num_atom;++i) {
      p_e[i]=0.0;p_LJ[i]=0.0;
      if (flagf==1) {
	for(j=0;j<3;++j) {
	  f_LJ[i*3+j]=0.0;
	  f_e[i*3+j]=0.0;
	}
      }
    }
  }

  if (flagp==2 ) {
    for(i=0;i<numnb;++i) {
      p_e[i]=0.0;
      p_LJ[i]=0.0; 
      if (flagf==2) {
	for(j=0;j<3;++j) {
	  f_LJ[i*3+j]=0.0;
	  f_e[i*3+j]=0.0;
	}
      }
    }
  }

  for(i=0;i<numnb;++i){
    num_a_prot=indexnb[i*2];
    NUM_A_PROT=indexnb[i*2+1];
    if (strncmp(AP.IGRAPH[num_a_prot],"H",1)!=0 && strncmp(AP.IGRAPH[NUM_A_PROT],"H",1)!=0 ) {
      len2 = 0.0;
      for(j=0;j<3;++j){
	vec[j] = cord[NUM_A_PROT*3+j]-cord[num_a_prot*3+j];
	len2 += vec[j]*vec[j];
      }
      len = sqrt(len2);
      potedummy=ele[num_a_prot]*ele[NUM_A_PROT]/(len);
      if (flagp==1) {
	p_e[num_a_prot] += potedummy;
	p_e[NUM_A_PROT] += potedummy;
      }
      else if (flagp==2) p_e[i] = potedummy;

      if (flagf==1) {
	for(j=0; j<3; ++j) {
	  fdummy[j] = -potedummy*UNIT*vec[j]/len2;
	  f_e[num_a_prot*3+j] += fdummy[j];
	  f_e[NUM_A_PROT*3+j] -= fdummy[j];
	}
      }
      else if (flagf==2) {
	for(j=0; j<3; ++j) {
	  fdummy[j] = -potedummy*UNIT*vec[j]/len2;
	  f_e[i*3+j] += fdummy[j];
	}
      }
      len6=len2;
      for (j=0;j<2;++j)  len6 = len6*len2;
      p12 = ALJ[num_a_prot*num_atom+NUM_A_PROT]/(len6*len6);
      p6  = BLJ[num_a_prot*num_atom+NUM_A_PROT]/len6;
      if (flagp==1) {
	p_LJ[num_a_prot] += p12-p6;
	p_LJ[NUM_A_PROT] += p12-p6;
      }
      else if (flagp==2) p_LJ[i] = p12-p6;

      if (flagf==1) {
	for(j=0;j<3;++j) {
	  fdummy[j]=-6.0*(2.0*p12-p6)/(len2)*vec[j]*UNIT;
	  f_LJ[num_a_prot*3+j] +=fdummy[j];
	  f_LJ[NUM_A_PROT*3+j] -= fdummy[j];
	}
      }
      else if (flagf==2) {
	for(j=0;j<3;++j) {
	  fdummy[j]=-6.0*(2.0*p12-p6)/(len2)*vec[j]*UNIT;
	  f_LJ[i*3+j] +=fdummy[j];
	}
      }
    }
  }

  return 0;
}

int ffL_calcFFNBvdW_woH( double *ALJ, double *BLJ,
			 double *p_LJ,
			 double *f_LJ,
			 int numnb, int *indexnb,
			 int num_atom,
			 double *cord,
			 int flagp, int flagf) {
  int i,j;
  int num_a_prot,NUM_A_PROT;
  double potedummy;
  double f[3];
  double len,len2,len6;
  double vec[3],fdummy[3];
  double p12,p6;

  if (flagf != 0 && flagf != 1 && flagf !=2 ) {
    printf("error !\n");exit(1);
  }
  if (flagp != 0 && flagp != 1 && flagp !=2 ) {
    printf("error !\n");exit(1);
  }

  if (flagp==1 ) {
    for(i=0;i<num_atom;++i) {
      p_LJ[i]=0.0;
      if (flagf==1) {
	for(j=0;j<3;++j) {
	  f_LJ[i*3+j]=0.0;
	}
      }
    }
  }

  if (flagp==2 ) {
    for(i=0;i<numnb;++i) {
      p_LJ[i]=0.0; 
      if (flagf==2) {
	for(j=0;j<3;++j) {
	  f_LJ[i*3+j]=0.0;
	}
      }
    }
  }

  for(i=0;i<numnb;++i){
    num_a_prot=indexnb[i*2];
    NUM_A_PROT=indexnb[i*2+1];
    if (strncmp(AP.IGRAPH[num_a_prot],"H",1)!=0 && strncmp(AP.IGRAPH[NUM_A_PROT],"H",1)!=0 ) {
      len2 = 0.0;
      for(j=0;j<3;++j){
	vec[j] = cord[NUM_A_PROT*3+j]-cord[num_a_prot*3+j];
	len2 += vec[j]*vec[j];
      }
      len = sqrt(len2);
      len6=len2;
      for (j=0;j<2;++j)  len6 = len6*len2;
      p12 = ALJ[num_a_prot*num_atom+NUM_A_PROT]/(len6*len6);
      p6  = BLJ[num_a_prot*num_atom+NUM_A_PROT]/len6;
      if (flagp==1) {
	p_LJ[num_a_prot] += p12-p6;
	p_LJ[NUM_A_PROT] += p12-p6;
      }
      else if (flagp==2) p_LJ[i] = p12-p6;

      if (flagf==1) {
	for(j=0;j<3;++j) {
	  fdummy[j]=-6.0*(2.0*p12-p6)/(len2)*vec[j]*UNIT;
	  f_LJ[num_a_prot*3+j] +=fdummy[j];
	  f_LJ[NUM_A_PROT*3+j] -= fdummy[j];
	}
      }
      else if (flagf==2) {
	for(j=0;j<3;++j) {
	  fdummy[j]=-6.0*(2.0*p12-p6)/(len2)*vec[j]*UNIT;
	  f_LJ[i*3+j] +=fdummy[j];
	}
      }
    }
  }

  return 0;
}

int ffL_calcTorque(double *Q,double *cord,int numclut,int *nNumClutOfParent,int *terminal_atom_a,int *origin_atom_a) {
  int i,j,k,l,n;
  int dtype,flag;  
  double atom[4][3];
  double dihedang;
  int *inpindex,inpnumH,inpnumA,*indexclut;
  int ON=1,OFF=0;

  indexclut=(int *)gcemalloc(sizeof(int)*(AP.NPHIH+AP.MPHIA));
  //  inpindex=(int *)gcemalloc(sizeof(int)*1);
  inpindex=ffL_make_inpindex(&inpnumH,&inpnumA,indexclut,numclut,nNumClutOfParent,terminal_atom_a,origin_atom_a);

  for (i=0;i<AP.NPHIH;++i) {
    flag=ON;
    for (j=0;j<inpnumH;++j) {
      if (i == inpindex[j]) {
	flag=OFF;
	break;
      }
    }

    if (flag==ON) {
      dtype = AP.PH[i][4]-1;
      n = indexclut[i]-1;
      for (j=0;j<4;++j) for (k=0;k<3;++k) atom[j][k]=cord[abs(AP.PH[i][j])+k];
      dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
      //      Q[n] += -AP.PK[dtype]*4.18407*100.0*(sin(AP.PN[dtype]*dihedang-AP.PHASE[dtype])*AP.PN[dtype]);
      Q[n] += -AP.PK[dtype]*(sin(AP.PN[dtype]*dihedang-AP.PHASE[dtype])*AP.PN[dtype]);
    }
  }

  for (i=0;i<numclut;++i) {
    Q[i] = 4.18407*100.0*Q[i];
  } 

  for (i=0;i<AP.MPHIA;++i) {
    flag=ON;
    for (j=0;j<inpnumA;++j) {
      if (i == inpindex[j+inpnumH]) {
	flag=OFF;
	break;
      }
    }

    if (flag==ON) {
      dtype = AP.PA[i][4]-1;
      n = indexclut[i+AP.NPHIH]-1;
      for (j=0;j<4;++j) for (k=0;k<3;++k) atom[j][k]=cord[abs(AP.PA[i][j])+k];
  
      dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
      Q[n] += -AP.PK[dtype]*4.18407*100.0*(sin(AP.PN[dtype]*dihedang-AP.PHASE[dtype])*AP.PN[dtype]);
    }
  }
}

double ffL_calcTorque_wtune(double *Q,double *p_d,double *cord,int numclut,int *nNumClutOfParent,int *terminal_atom_a,int *origin_atom_a,int *atom_tune_pairs, double *tune_val_by_period, int numtune) {
  int i,j,k,l,n;
  int dtype,flag;  
  double atom[4][3];
  double dihedang;
  int *inpindex,inpnumH,inpnumA,*indexclut;
  int ON=1,OFF=0;

  int tuneflag,tuneindex;
  double p_d_t=0.0;

  indexclut=(int *)gcemalloc(sizeof(int)*(AP.NPHIH+AP.MPHIA));
  //  inpindex=(int *)gcemalloc(sizeof(int)*1);
  inpindex=ffL_make_inpindex(&inpnumH,&inpnumA,indexclut,numclut,nNumClutOfParent,terminal_atom_a,origin_atom_a);

  for (i=0;i<AP.NPHIH;++i) {
    flag=ON;
    for (j=0;j<inpnumH;++j) {
      if (i == inpindex[j]) {
	flag=OFF;
	break;
      }
    }

    if (flag==ON) {
      dtype = AP.PH[i][4]-1;
      n = indexclut[i]-1;
      for (j=0;j<4;++j) for (k=0;k<3;++k) atom[j][k]=cord[abs(AP.PH[i][j])+k];

      tuneflag=OFF;
      for (j=0;j<numtune;++j) {
	if (atom_tune_pairs[j*4+0]==abs(AP.PH[i][0])/3 && atom_tune_pairs[j*4+1]==abs(AP.PH[i][1])/3 
	    && atom_tune_pairs[j*4+2]==abs(AP.PH[i][2])/3 && atom_tune_pairs[j*4+3]==abs(AP.PH[i][3])/3) {
	  tuneindex=j;
	  tuneflag=ON;
	  break;
	}
      }

      dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
      //      Q[n] += -AP.PK[dtype]*4.18407*100.0*(sin(AP.PN[dtype]*dihedang-AP.PHASE[dtype])*AP.PN[dtype]);
      if (tuneflag==OFF) {
	p_d[i] = AP.PK[dtype]*(1.0+cos(AP.PN[dtype]*dihedang-AP.PHASE[dtype]));
	p_d_t+=p_d[i];
	Q[n] += -AP.PK[dtype]*(sin(AP.PN[dtype]*dihedang-AP.PHASE[dtype])*AP.PN[dtype]);
      }
      else {
	p_d[i] = 1.0/(tune_val_by_period[tuneindex]*tune_val_by_period[tuneindex])
	  *AP.PK[dtype]*(1.0+cos(AP.PN[dtype]*dihedang-AP.PHASE[dtype]));
	p_d_t+=p_d[i];
	Q[n] += -1.0/(tune_val_by_period[tuneindex]*tune_val_by_period[tuneindex])
	  *AP.PK[dtype]*(sin(AP.PN[dtype]*dihedang-AP.PHASE[dtype])*AP.PN[dtype]);
      }
    }
  }

  for (i=0;i<numclut;++i) {
    Q[i] = 4.18407*100.0*Q[i];
  } 

  for (i=0;i<AP.MPHIA;++i) {
    flag=ON;
    for (j=0;j<inpnumA;++j) {
      if (i == inpindex[j+inpnumH]) {
	flag=OFF;
	break;
      }
    }

    if (flag==ON) {
      dtype = AP.PA[i][4]-1;
      n = indexclut[i+AP.NPHIH]-1;
      for (j=0;j<4;++j) for (k=0;k<3;++k) atom[j][k]=cord[abs(AP.PA[i][j])+k];

      tuneflag=OFF;
      for (j=0;j<numtune;++j) {
	if (atom_tune_pairs[j*4+0]==abs(AP.PA[i][0])/3 && atom_tune_pairs[j*4+1]==abs(AP.PA[i][1])/3 
	    && atom_tune_pairs[j*4+2]==abs(AP.PA[i][2])/3 && atom_tune_pairs[j*4+3]==abs(AP.PA[i][3])/3) {
	  tuneindex=j;
	  tuneflag=ON;
	  break;
	}
      }
  
      dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
      if (tuneflag==OFF) {
	p_d[i+AP.NPHIH] = AP.PK[dtype]*(1.0+cos(AP.PN[dtype]*dihedang-AP.PHASE[dtype]));
	p_d_t+=p_d[i];
	Q[n] += -AP.PK[dtype]*4.18407*100.0*(sin(AP.PN[dtype]*dihedang-AP.PHASE[dtype])*AP.PN[dtype]);
      }
      else { 
	p_d[i+AP.NPHIH] = 1.0/(tune_val_by_period[tuneindex]*tune_val_by_period[tuneindex])
	  *AP.PK[dtype]*(1.0+cos(AP.PN[dtype]*dihedang-AP.PHASE[dtype]));
	p_d_t+=p_d[i];
	Q[n] += -1.0/(tune_val_by_period[tuneindex]*tune_val_by_period[tuneindex])
	  *AP.PK[dtype]*4.18407*100.0*(sin(AP.PN[dtype]*dihedang-AP.PHASE[dtype])*AP.PN[dtype]);
      }
    }
  }

  return p_d_t;
}

double ffL_calcTorque_wtuneb(double *Q,double *p_d,double *cord,int numclut,int *nNumClutOfParent,int *terminal_atom_a,int *origin_atom_a,int *atom_tune_pairs, double *tune_val_V_n, double *tune_val_n_phase, int numtune, double pi) {
  int i,j,k,l,n;
  int dtype,flag;  
  double atom[4][3];
  double dihedang;
  int *inpindex,inpnumH,inpnumA,*indexclut;
  int ON=1,OFF=0;

  int tuneflag,tuneindex;
  double p_d_t=0.0;

  double N,PHASE;

  indexclut=(int *)gcemalloc(sizeof(int)*(AP.NPHIH+AP.MPHIA));
  //  inpindex=(int *)gcemalloc(sizeof(int)*1);
  inpindex=ffL_make_inpindex(&inpnumH,&inpnumA,indexclut,numclut,nNumClutOfParent,terminal_atom_a,origin_atom_a);

  for (i=0;i<AP.NPHIH;++i) {
    flag=ON;
    for (j=0;j<inpnumH;++j) {
      if (i == inpindex[j]) {
	flag=OFF;
	break;
      }
    }

    if (flag==ON) {
      dtype = AP.PH[i][4]-1;
      n = indexclut[i]-1;
      for (j=0;j<4;++j) for (k=0;k<3;++k) atom[j][k]=cord[abs(AP.PH[i][j])+k];

      tuneflag=OFF;
      for (j=0;j<numtune;++j) {
	if (atom_tune_pairs[j*4+0]==abs(AP.PH[i][0])/3 && atom_tune_pairs[j*4+1]==abs(AP.PH[i][1])/3 
	    && atom_tune_pairs[j*4+2]==abs(AP.PH[i][2])/3 && atom_tune_pairs[j*4+3]==abs(AP.PH[i][3])/3) {
	  tuneindex=j;
	  tuneflag=ON;
	  break;
	}
      }

      dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
      //      Q[n] += -AP.PK[dtype]*4.18407*100.0*(sin(AP.PN[dtype]*dihedang-AP.PHASE[dtype])*AP.PN[dtype]);
      if (tuneflag==OFF) {
	p_d[i] = AP.PK[dtype]*(1.0+cos(AP.PN[dtype]*dihedang-AP.PHASE[dtype]));
	p_d_t+=p_d[i];
	Q[n] += -AP.PK[dtype]*(sin(AP.PN[dtype]*dihedang-AP.PHASE[dtype])*AP.PN[dtype]);
      }
      else {
	N=tune_val_n_phase[tuneindex]*AP.PN[dtype];
	PHASE=2.0*pi/N;
 	if (PHASE >= 2.0*pi) PHASE-=2.0*pi;
	p_d[i] = tune_val_V_n[tuneindex]*AP.PK[dtype]*(1.0+cos(N*dihedang-PHASE));
	p_d_t+=p_d[i];
	Q[n] += -tune_val_V_n[tuneindex]*AP.PK[dtype]*(sin(N*dihedang-PHASE)*N);
      }
    }
  }

  for (i=0;i<numclut;++i) {
    Q[i] = 4.18407*100.0*Q[i];
  } 

  for (i=0;i<AP.MPHIA;++i) {
    flag=ON;
    for (j=0;j<inpnumA;++j) {
      if (i == inpindex[j+inpnumH]) {
	flag=OFF;
	break;
      }
    }

    if (flag==ON) {
      dtype = AP.PA[i][4]-1;
      n = indexclut[i+AP.NPHIH]-1;
      for (j=0;j<4;++j) for (k=0;k<3;++k) atom[j][k]=cord[abs(AP.PA[i][j])+k];

      tuneflag=OFF;
      for (j=0;j<numtune;++j) {
	if (atom_tune_pairs[j*4+0]==abs(AP.PA[i][0])/3 && atom_tune_pairs[j*4+1]==abs(AP.PA[i][1])/3 
	    && atom_tune_pairs[j*4+2]==abs(AP.PA[i][2])/3 && atom_tune_pairs[j*4+3]==abs(AP.PA[i][3])/3) {
	  tuneindex=j;
	  tuneflag=ON;
	  break;
	}
      }
  
      dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
      if (tuneflag==OFF) {
	p_d[i+AP.NPHIH] = AP.PK[dtype]*(1.0+cos(AP.PN[dtype]*dihedang-AP.PHASE[dtype]));
	p_d_t+=p_d[i];
	Q[n] += -AP.PK[dtype]*4.18407*100.0*(sin(AP.PN[dtype]*dihedang-AP.PHASE[dtype])*AP.PN[dtype]);
      }
      else { 
	N=tune_val_n_phase[tuneindex]*AP.PN[dtype];
	PHASE=2.0*pi/N;
 	if (PHASE >= 2.0*pi) PHASE-=2.0*pi;
	p_d[i] = tune_val_V_n[tuneindex]*AP.PK[dtype]*(1.0+cos(N*dihedang-PHASE));
	p_d_t+=p_d[i];
	Q[n] += -tune_val_V_n[tuneindex]*AP.PK[dtype]*(sin(N*dihedang-PHASE)*N);
      }
    }
  }

  return p_d_t;
}



int ffL_calcTorque_woH(double *Q,double *cord,int numclut,int *nNumClutOfParent,int *terminal_atom_a,int *origin_atom_a) {
  int i,j,k,l,n;
  int dtype,flag;  
  double atom[4][3];
  double dihedang;
  int *inpindex,inpnumH,inpnumA,*indexclut;
  int ON=1,OFF=0;

  indexclut=(int *)gcemalloc(sizeof(int)*(AP.NPHIH+AP.MPHIA));
  //  inpindex=(int *)gcemalloc(sizeof(int)*1);
  inpindex=ffL_make_inpindex(&inpnumH,&inpnumA,indexclut,numclut,nNumClutOfParent,terminal_atom_a,origin_atom_a);

  for (i=0;i<AP.MPHIA;++i) {
    flag=ON;
    for (j=0;j<inpnumA;++j) {
      if (i == inpindex[j+inpnumH]) {
	flag=OFF;
	break;
      }
    }

    if (flag==ON) {
      dtype = AP.PA[i][4]-1;
      n = indexclut[i+AP.NPHIH]-1;
      for (j=0;j<4;++j) for (k=0;k<3;++k) atom[j][k]=cord[abs(AP.PA[i][j])+k];
  
      dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
      Q[n] += -AP.PK[dtype]*4.18407*100.0*(sin(AP.PN[dtype]*dihedang-AP.PHASE[dtype])*AP.PN[dtype]);
    }
  }
}



int *ffL_make_inpindex(int *inpnumH,int *inpnumA,int *indexclut,int numclut,int *nNumClutOfParent,int *terminal_atom_a,int *origin_atom_a) {
  int i,j,k,l,ll,p;
  int flag;
  int *atom_dihed_pair;
  int ON=1,OFF=0;
  int *inpindexH,*inpindexA;
  int *inpindex;

  (*inpnumH)=0;
  (*inpnumA)=0;
  
  atom_dihed_pair=(int *)gcemalloc(sizeof(int)*(AP.NPHIH+AP.MPHIA)*6);
  inpindexH=(int *)gcemalloc(sizeof(int)*1);
  inpindexA=(int *)gcemalloc(sizeof(int)*1);
  
  for (k=0;k<AP.NPHIH;++k) {
    flag=OFF;
    for (i=0;i<AP.NBONH;++i) {
      if (((AP.BH[i][0] == abs(AP.PH[k][1]) && AP.BH[i][1] == abs(AP.PH[k][0]))
  	   || (AP.BH[i][0] == abs(AP.PH[k][0]) && AP.BH[i][1] == abs(AP.PH[k][1])))) {
  	flag = ON;
  	break;
      }
    }
    if (flag==OFF) {
      for (i=0;i<AP.NBONA;++i) {
  	if (((AP.BA[i][0] == abs(AP.PH[k][1]) && AP.BA[i][1] == abs(AP.PH[k][0]))
  	     || (AP.BA[i][0] == abs(AP.PH[k][0]) && AP.BA[i][1] == abs(AP.PH[k][1])))) {
  	  flag = ON;
  	  break;
  	}
      }
    }
    if (flag==ON) {
      flag=OFF;
      for (i=0;i<AP.NBONH;++i) {
  	if (((AP.BH[i][0] == abs(AP.PH[k][2]) && AP.BH[i][1] == abs(AP.PH[k][3]))
  	     || (AP.BH[i][0] == abs(AP.PH[k][3]) && AP.BH[i][1] == abs(AP.PH[k][2]))) ) {
  	  flag = ON;
  	  break;
  	}
      }
      if (flag==OFF ) {
  	for (i=0;i<AP.NBONA;++i) {
  	  if (((AP.BA[i][0] == abs(AP.PH[k][2]) && AP.BA[i][1] == abs(AP.PH[k][3]))
  	       || (AP.BA[i][0] == abs(AP.PH[k][3]) && AP.BA[i][1] == abs(AP.PH[k][2]))) ) {
  	    flag = ON;
  	    break;
  	  }
  	}
      }
    }
    if (flag==OFF) {
      inpindexH[(*inpnumH)]=k;
      ++(*inpnumH);
      inpindexH=(int *)gcerealloc(inpindexH,sizeof(int)*(*inpnumH));
    }
  
    for (i=0;i<4;++i) {
      atom_dihed_pair[k/*][*/*6+i] = abs(AP.PH[k][i])/3+1;
    }
    atom_dihed_pair[k/*][*/*6+4] = AP.PH[k][i];
    l = 1;
    ll = 2;
    if (atom_dihed_pair[k/*][*/*6+1] > atom_dihed_pair[k/*][*/*6+2]) {
      l = 2;
      ll = 1;
    }
  
    for (i=0;i<numclut;++i){
      p = nNumClutOfParent[i]-1;
      if (atom_dihed_pair[k/*][*/*6+l]==terminal_atom_a[p] && atom_dihed_pair[k/*][*/*6+ll]==origin_atom_a[i]) {
  	atom_dihed_pair[k/*][*/*6+5] = i+1;
      }
    }
  }
  
  for (k=0;k<AP.MPHIA;++k) {
    // check for improper dihed
    flag=OFF;
    for (i=0;i<AP.NBONH;++i) {
      if (   ((AP.BH[i][0] == abs(AP.PA[k][1]) && AP.BH[i][1] == abs(AP.PA[k][0])) || (AP.BH[i][0] == abs(AP.PA[k][0]) && AP.BH[i][1] == abs(AP.PA[k][1])))) {
  	flag = ON;
  	break;
      }
    }
    if (flag==OFF) {
      for (i=0;i<AP.NBONA;++i) {
  	if (   ((AP.BA[i][0] == abs(AP.PA[k][1]) && AP.BA[i][1] == abs(AP.PA[k][0])) || (AP.BA[i][0] == abs(AP.PA[k][0]) && AP.BA[i][1] == abs(AP.PA[k][1])))) {
  	  flag = ON;
  	  break;
  	}
      }
    }
    if (flag==ON) {
      flag=OFF;
      for (i=0;i<AP.NBONH;++i) {
  	if (((AP.BH[i][0] == abs(AP.PA[k][2]) && AP.BH[i][1] == abs(AP.PA[k][3])) || (AP.BH[i][0] == abs(AP.PA[k][3]) && AP.BH[i][1] == abs(AP.PA[k][2]))) ) {
  	  flag = ON;
  	  break;
  	}
      }
      if (flag==OFF ) {
  	for (i=0;i<AP.NBONA;++i) {
  	  if (((AP.BA[i][0] == abs(AP.PA[k][2]) && AP.BA[i][1] == abs(AP.PA[k][3])) || (AP.BA[i][0] == abs(AP.PA[k][3]) && AP.BA[i][1] == abs(AP.PA[k][2]))) ) {
  	    flag = ON;
  	    break;
  	  }
  	}
      }
    }
    if (flag==OFF) {
      inpindexA[(*inpnumA)]=k;
      ++(*inpnumA);
      inpindexA=(int *)gcerealloc(inpindexA,sizeof(int)*(*inpnumA));
    }
  
    for (i=0;i<4;++i) {
      atom_dihed_pair[(k+AP.NPHIH)/*][*/*6+i] = abs(AP.PA[k][i])/3+1;
    }
    atom_dihed_pair[(k+AP.NPHIH)/*][*/*6+4] = AP.PA[k][4];
    l = 1;
    ll = 2;
    if (atom_dihed_pair[(k+AP.NPHIH)/*][*/*6+1] > atom_dihed_pair[(k+AP.NPHIH)/*][*/*6+2]) {
      l = 2;
      ll = 1;
    }
    for (i=0;i<numclut;++i){
      p = nNumClutOfParent[i]-1;
      if (atom_dihed_pair[(k+AP.NPHIH)/*][*/*6+l]==terminal_atom_a[p] && atom_dihed_pair[(k+AP.NPHIH)/*][*/*6+ll]==origin_atom_a[i]){
  
  	atom_dihed_pair[(k+AP.NPHIH)/*][*/*6+5] = i+1;
      }
    }
  }
  
  for (i=0;i<AP.NPHIH+AP.MPHIA;++i) {
    indexclut[i]=atom_dihed_pair[i*6+5];
  }

  //  inpindex=(int *)gcerealloc(inpindex,sizeof(int)*((*inpnumH)+(*inpnumA)));
  inpindex=(int *)gcemalloc(sizeof(int)*((*inpnumH)+(*inpnumA)));
  for (i=0;i<(*inpnumH);++i)
    inpindex[i]=inpindexH[i];
  for (i=0;i<(*inpnumA);++i)
    inpindex[i+(*inpnumH)]=inpindexA[i];

  return inpindex;
}

int *ffL_make_inpindex_A(int *inpnumA) {
  int i,j,k,l,ll,p;
  int flag;
  int ON=1,OFF=0;
  int *inpindexA;

  (*inpnumA)=0;
  
  inpindexA=(int *)gcemalloc(sizeof(int)*AP.MPHIA);

  for (k=0;k<AP.MPHIA;++k) {
    // check for improper dihed
    flag=OFF;
    for (i=0;i<AP.NBONH;++i) {
      if ( ((AP.BH[i][0] == abs(AP.PA[k][1]) && AP.BH[i][1] == abs(AP.PA[k][0])) || (AP.BH[i][0] == abs(AP.PA[k][0]) && AP.BH[i][1] == abs(AP.PA[k][1])))) {
  	flag = ON;
  	break;
      }
    }
    if (flag==OFF) {
      for (i=0;i<AP.NBONA;++i) {
  	if (((AP.BA[i][0] == abs(AP.PA[k][1]) && AP.BA[i][1] == abs(AP.PA[k][0])) || (AP.BA[i][0] == abs(AP.PA[k][0]) && AP.BA[i][1] == abs(AP.PA[k][1])))) {
  	  flag = ON;
  	  break;
  	}
      }
    }
    if (flag==ON) {
      flag=OFF;
      for (i=0;i<AP.NBONH;++i) {
  	if (((AP.BH[i][0] == abs(AP.PA[k][2]) && AP.BH[i][1] == abs(AP.PA[k][3])) || (AP.BH[i][0] == abs(AP.PA[k][3]) && AP.BH[i][1] == abs(AP.PA[k][2]))) ) {
  	  flag = ON;
  	  break;
  	}
      }
      if (flag==OFF ) {
  	for (i=0;i<AP.NBONA;++i) {
  	  if (((AP.BA[i][0] == abs(AP.PA[k][2]) && AP.BA[i][1] == abs(AP.PA[k][3])) || (AP.BA[i][0] == abs(AP.PA[k][3]) && AP.BA[i][1] == abs(AP.PA[k][2]))) ) {
  	    flag = ON;
  	    break;
  	  }
  	}
      }
    }
    if (flag==OFF) {
      inpindexA[k]=-1;
      ++(*inpnumA);
    }
  }

  return inpindexA;
}


void ffL_out_formated(FILE *outputfile,struct potential e,double KE,double KEv,double PEv,double T,int i,double dt) {

  fprintf(outputfile,"/***********************************************/\n");
  fprintf(outputfile,"steps            = %d  \n",i);
  fprintf(outputfile,"total time       = %10.3lf ps  \n" ,dt*(double)i);
  fprintf(outputfile,"T_kelvin         = %e K  \n",T);
  fprintf(outputfile,"toal_energy      = %e kcal/mol  \n",e.p_t+KE);
  fprintf(outputfile,"toal_vertial_energy      = %e kcal/mol  \n",e.p_t+KE+KEv+PEv);
  fprintf(outputfile,"kinetic_energy   = %e kcal/mol  \n",KE);
  fprintf(outputfile,"kinetic_energy_vertial   = %e kcal/mol  \n",KEv);
  fprintf(outputfile,"potential_energy_real = %e kcal/mol  \n",e.p_t);
  fprintf(outputfile,"potential_energy_vertial   = %e kcal/mol  \n",PEv);
  fprintf(outputfile,"dihedral_energy  = %e kcal/mol  \n",e.p_d_t);
  fprintf(outputfile,"elect_energy     = %e kcal/mol  \n",e.p_e_t);
  fprintf(outputfile,"VDW_energy       = %e kcal/mol  \n",e.p_LJ_t);
  fprintf(outputfile,"1_4_elect_energy = %e kcal/mol  \n",e.p_e_14_t);
  fprintf(outputfile,"1_4_VDW_energy   = %e kcal/mol  \n",e.p_LJ_14_t);

}

///////////////////////////////////////////////////////////////////
