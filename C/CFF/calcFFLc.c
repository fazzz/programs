
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "FFLc.h"
#include "FFL.h"
#include "MB.h"

#include "PTLb.h"

#include "LA.h"
#include "TOPO.h"
#include "mymath.h"
#include "EF.h"

#define YES 1
#define NO  0

double calcANGKE_forcec(double atomi[3],double atomj[3],double atomk[3],double kang,double ang_eq,double *f);

int ffLc_calcFFNB(double *ele, double *ALJ, double *BLJ,
		double *p_e,double *p_LJ,
		double *f_e,double *f_LJ,
		int numnb, int *indexnb,
		int num_atom,
		double *cord,
		int flagp, int flagf,struct AmberParmL ap) {
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

int ffLc_calcFFNB_wtune(double *ele, double *ALJ, double *BLJ,
		       double *p_e,double *p_LJ,
		       double *f_e,double *f_LJ,
		       int numnb, int *indexnb,
		       int num_atom,
		       double *cord,
		       int flagp, int flagf,
		       int *atom_tune_pairs, double *tune_val, int numtune,struct AmberParmL ap) {
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

int ffLc_calcFFNB_wtuneb(double *ele, double *ALJ, double *BLJ,
			double *p_e,double *p_LJ,
			double *f_e,double *f_LJ,
			int numnb, int *indexnb,
			int num_atom,
			double *cord,
			int flagp, int flagf,
			int *atom_tune_pairs_es, int *atom_tune_pairs_LJ,
			double *tune_val_es, double *tune_val_LJ,
			int numtune_es, int numtune_LJ,struct AmberParmL ap) {
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

int ffLc_calcFFNB_14(double *ele, double *ALJ, double *BLJ,
		    double *p_e,double *p_LJ,
		    double *f_e,double *f_LJ,
		    int numnb, int *indexnb,
		    int num_atom,
		    double *cord,
		    int flagp, int flagf,struct AmberParmL ap) {
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

int ffLc_calcFFNB_14_6(double *ele, double *ALJ, double *BLJ, // 0911
		     double *p_e,double *p_LJ,
		     double *f_e,double *f_LJ,
		     int numnb, int *indexnb,
		     int num_atom,
		     double *cord,
		     int flagp, int flagf,
		     double scale,
		     int numrepul,struct AmberParmL ap) {
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


int ffLc_calcFFNB_wpep(double *ele, double *ALJ, double *BLJ,
		     double *p_e,double *p_LJ,
		     double *f_e,double *f_LJ,
		     int numnb, int *indexnb,
		     int num_atom,
		     double *cord,
		     double *u_1_5, double fact,
		     int flagp, int flagf,struct AmberParmL ap) {
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


void ffLc_set_NB_PARM(double *ele, double *ALJ, double *BLJ, int numatom,struct AmberParmL ap){
  int i,j,index;
  
  for (i=0;i<numatom;++i) {
    ele[i]= ap.CHRG[i];
    for (j=0;j<numatom;++j) {
      index=ap.ICO[(ap.IAC[i]-1)*ap.NTYPES+(ap.IAC[j]-1)]-1;
      ALJ[i*numatom+j] = ap.CN1[index];
      BLJ[i*numatom+j] = ap.CN2[index];
    }
  }
}

int ffLc_set_NB_index(int *indexnb,int numnb, int numatom,struct AmberParmL ap) {
  int i,j,k,l,m,n;
  int flag;
  FILE *log;

  l=0;n=0;
  for(i=0;i<numatom;++i){
    for (j=i+1;j<numatom;++j) {
      flag=1;
      for (k=0;k<ap.NUMEX[i];++k) {
	if (ap.NATEX[(n+k)]-1==j) {
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
    n+=ap.NUMEX[i];
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

int ffLc_set_numnb(struct AmberParmL ap) {
  int i,numnb;

  numnb=0;
  for(i=0;i<ap.NATOM;++i)
    numnb+=ap.NATOM-i-1-ap.NUMEX[i];

  return numnb+1;
}

int ffLc_calcDIHE(double *p_d,
		double *n_d,
		double *cord,
		int flagp, int flagf, int flaginp,struct AmberParmL ap) {
  int i,j,k,l,flag;
  int dtype;
  
  double atom[4][3];
  double dihedang;
  
  for (i=0;i<ap.NPHIH;++i){
    p_d[i] = 0.0;
    if (flagf==1)
      n_d[i] = 0.0;
  }
  
  for (i=0;i<ap.MPHIA;++i){
    p_d[i+ap.NPHIH] = 0.0;
    if (flagf==1)
      n_d[i+ap.NPHIH] = 0.0;
  }
  
  for (i=0;i<ap.NPHIH;++i) {
    dtype = ap.PH[i][4]-1;
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
  	  atom[j][k]=cord[abs(ap.PH[i][j])+k];
  	}
      }
  
      dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
      p_d[i] = ap.PK[dtype]*(1.0+cos(ap.PN[dtype]*dihedang-ap.PHASE[dtype]));
      if (flagf==1)
  	n_d[i] = -ap.PK[dtype]*4.18407*100.0*(sin(ap.PN[dtype]*dihedang-ap.PHASE[dtype])*ap.PN[dtype]);
    }
  }
  
  for (i=0;i<ap.MPHIA;++i) {
    dtype = ap.PA[i][4]-1;
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
  	  atom[j][k]=cord[abs(ap.PA[i][j])+k];
  	}
      }
  
      dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
      p_d[i+ap.NPHIH] = ap.PK[dtype]*(1.0+cos(ap.PN[dtype]*dihedang-ap.PHASE[dtype]));
      if (flagf==1)
  	n_d[i+ap.NPHIH] = -ap.PK[dtype]*4.18407*100.0*(sin(ap.PN[dtype]*dihedang-ap.PHASE[dtype])*ap.PN[dtype]);
    }
  }
}

int ffLc_calcDIHE_force_Cartesian(double *f_d,double *cord,struct AmberParmL ap) {
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

  //  for (i=0;i<(ap.NPHIH+ap.MPHIA)*3;++i) f_d[i] = 0.0;
  for (i=0;i<ap.NATOM*3;++i) f_d[i] = 0.0;
  
  for (i=0;i<ap.NPHIH;++i) {
    dtype = ap.PH[i][4]-1;
    for (j=0;j<4;++j) for (k=0;k<3;++k) atom[j][k]=cord[abs(ap.PH[i][j])+k];

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

    if (ap.PN[dtype]==1) fa=-ap.PK[dtype]*ap.PN[dtype]*cos(ap.PHASE[dtype]);
    else if (ap.PN[dtype]==2) fa=-ap.PK[dtype]*ap.PN[dtype]*2.0*cosdih*cos(ap.PHASE[dtype]);
    else if (ap.PN[dtype]==3) fa=-ap.PK[dtype]*ap.PN[dtype]*(-4.0*sindih*sindih+3.0)*cos(ap.PHASE[dtype]);
    else if (ap.PN[dtype]==4) fa=-ap.PK[dtype]*ap.PN[dtype]*4.0*(cosdih*(2.0*cosdih*cosdih-1.0))*cos(ap.PHASE[dtype]);
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
      f_d[abs(ap.PH[i][0])+j] += fa*op1[j]*UNIT;
      f_d[abs(ap.PH[i][1])+j] += fa*(-op2[j]+op3[j])*UNIT;
      f_d[abs(ap.PH[i][2])+j] += fa*(-op4[j]+op5[j])*UNIT;
      f_d[abs(ap.PH[i][3])+j] += fa*op6[j]*UNIT;
    }
  }
  
  for (i=0;i<ap.MPHIA;++i) {
    dtype = ap.PA[i][4]-1;
    for (j=0;j<4;++j) for (k=0;k<3;++k)	atom[j][k]=cord[abs(ap.PA[i][j])+k];

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

    if (ap.PN[dtype]==1) fa=-ap.PK[dtype]*ap.PN[dtype]*cos(ap.PHASE[dtype]);
    else if (ap.PN[dtype]==2) fa=-ap.PK[dtype]*ap.PN[dtype]*2.0*cosdih*cos(ap.PHASE[dtype]);
    else if (ap.PN[dtype]==3) fa=-ap.PK[dtype]*ap.PN[dtype]*(-4.0*sindih*sindih+3.0)*cos(ap.PHASE[dtype]);
    else if (ap.PN[dtype]==4) fa=-ap.PK[dtype]*ap.PN[dtype]*4.0*(cosdih*(2.0*cosdih*cosdih-1.0))*cos(ap.PHASE[dtype]);
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
      f_d[abs(ap.PA[i][0])+j] += fa*op1[j]*UNIT;
      f_d[abs(ap.PA[i][1])+j] += fa*(-op2[j]+op3[j])*UNIT;
      f_d[abs(ap.PA[i][2])+j] += fa*(-op4[j]+op5[j])*UNIT;
      f_d[abs(ap.PA[i][3])+j] += fa*op6[j]*UNIT;
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

int ffLc_calc_spe_type_DIHE(double *p,
			  double *n,
			  double angle,
			  int type,struct AmberParmL ap) {

  *p = ap.PK[type]*(1.0+cos(ap.PN[type]*angle-ap.PHASE[type]));
  *n = ap.PK[type]*4.18407*100.0*(sin(ap.PN[type]*angle-ap.PHASE[type])*ap.PN[type]);
  
}


int ffLc_calcDIHE_wpep(double *p_d,
		     double *n_d,
		     double *cord,
		     double *u_d,
		     int flagp, int flagf, int flaginp,struct AmberParmL ap) {
  int i,j,k,l,flag;
  int dtype;
  
  double atom[4][3];
  double dihedang;
  
  for (i=0;i<ap.NPHIH;++i){
    p_d[i] = 0.0;
    if (flagf==1)
      n_d[i] = 0.0;
  }
  
  for (i=0;i<ap.MPHIA;++i){
    p_d[i+ap.NPHIH] = 0.0;
    if (flagf==1)
      n_d[i+ap.NPHIH] = 0.0;
  }
  
  for (i=0;i<ap.NPHIH;++i) {
    dtype = ap.PH[i][4]-1;
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
  	  atom[j][k]=cord[abs(ap.PH[i][j])+k];
  	}
      }
  
     dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
      p_d[i] = ap.PK[dtype]*(1.0+cos(ap.PN[dtype]*dihedang-ap.PHASE[dtype]))*u_d[i];
      if (flagf==1)
  	n_d[i] = -ap.PK[dtype]*4.18407*100.0*(sin(ap.PN[dtype]*dihedang-ap.PHASE[dtype])*ap.PN[dtype])*u_d[i];
    }
  }
  
  for (i=0;i<ap.MPHIA;++i) {
    dtype = ap.PA[i][4]-1;
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
  	  atom[j][k]=cord[abs(ap.PA[i][j])+k];
  	}
      }
  
      dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
      p_d[i] = ap.PK[dtype]*(1.0+cos(ap.PN[dtype]*dihedang-ap.PHASE[dtype]))*u_d[i+ap.NPHIH];
      if (flagf==1)
  	n_d[i] = -ap.PK[dtype]*4.18407*100.0*(sin(ap.PN[dtype]*dihedang-ap.PHASE[dtype])*ap.PN[dtype])*u_d[i+ap.NPHIH];
    }
  }
}

int ffLc_calcANGLE(double *p_a,double *cord,struct AmberParmL ap){
  int i,j,k,l;
  int type;
  
  double atom[3][3];
  double ang;
  
  for (i=0;i<ap.NTHETH+ap.MTHETA;++i)
    p_a[i] = 0.0;
  
  for (i=0;i<ap.NTHETH;++i) {
    type = ap.TH[i][3]-1;
     for (j=0;j<3;++j)
      for (k=0;k<3;++k)
	atom[j][k]=cord[abs(ap.TH[i][j])+k];
  
    ang = pick_angle(atom[0],atom[1],atom[2],0,0.0);
    p_a[i] = ap.TK[type]*(ang-ap.TEQ[type])*(ang-ap.TEQ[type]);
  }
  
  for (i=0;i<ap.MTHETA;++i) {
    type = ap.TA[i][3]-1;
     for (j=0;j<3;++j)
      for (k=0;k<3;++k)
	atom[j][k]=cord[abs(ap.TA[i][j])+k];
  
    ang = pick_angle(atom[0],atom[1],atom[2],0,0.0);
    p_a[i+ap.NTHETH] = ap.TK[type]*(ang-ap.TEQ[type])*(ang-ap.TEQ[type]);
  }

  return 0;
}

int ffLc_calcANGLE_force_Cartesian(double *f_a,double *cord,struct AmberParmL ap){
  int i,j,k,l;
  int numatom;
  int type;
  double kang,ang_eq;

  double atom[3][3];
  double *f_temp;

  numatom=ap.NATOM;
  for (i=0;i<numatom*3;++i) f_a[i] = 0.0;
  f_temp=(double *)gcemalloc(sizeof(double)*3*3);

  for (i=0;i<ap.NTHETH;++i) {
    for (j=0;j<3;++j)  for (k=0;k<3;++k) atom[j][k]=cord[abs(ap.TH[i][j])+k];
    type = ap.TH[i][3]-1;
    kang = ap.TK[type];
    ang_eq = ap.TEQ[type];
    calcANGKE_force(atom[0],atom[1],atom[2],kang,ang_eq,f_temp);
    for (j=0;j<3;++j) {
      f_a[abs(ap.TH[i][0])+j] += f_temp[j];
      f_a[abs(ap.TH[i][1])+j] += f_temp[3+j];
      f_a[abs(ap.TH[i][2])+j] += f_temp[6+j];
    }
  }

  for (i=0;i<ap.MTHETA;++i) {
    for (j=0;j<3;++j)  for (k=0;k<3;++k) atom[j][k]=cord[abs(ap.TA[i][j])+k];
    type = ap.TA[i][3]-1;
    kang = ap.TK[type];
    ang_eq = ap.TEQ[type];
    calcANGKE_force(atom[0],atom[1],atom[2],kang,ang_eq,f_temp);
    for (j=0;j<3;++j) {
      f_a[abs(ap.TA[i][0])+j] += f_temp[j];
      f_a[abs(ap.TA[i][1])+j] += f_temp[3+j];
      f_a[abs(ap.TA[i][2])+j] += f_temp[6+j];
    }
  }

  return 0;
}

double calcANGKE_forcec(double atomi[3],double atomj[3],double atomk[3],double kang,double ang_eq,double *f) {
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

int ffLc_calcBOND(double *p_b,double *cord,struct AmberParmL ap){
  int i,j,k;
  int type;
  double len;
  double atom[2][3];

  for (i=0;i<ap.NBONH+ap.MBONA;++i)
    p_b[i] = 0.0;
  
  for (i=0;i<ap.NBONH;++i) {
    type = ap.BH[i][2]-1;
    for (j=0;j<2;++j)
      for (k=0;k<3;++k)
	atom[j][k]=cord[abs(ap.BH[i][j])+k];
  
    len = pick_bond_leng(atom[0],atom[1]);
    p_b[i] = ap.RK[type]*(len-ap.REQ[type])*(len-ap.REQ[type]);
  }

  for (i=0;i<ap.MBONA;++i) {
    type = ap.BA[i][2]-1;
    for (j=0;j<2;++j)
      for (k=0;k<3;++k)
	atom[j][k]=cord[abs(ap.BA[i][j])+k];  

    len = pick_bond_leng(atom[0],atom[1]);
    p_b[i+ap.NBONH] = ap.RK[type]*(len-ap.REQ[type])*(len-ap.REQ[type]);
  }

  return 0;
}

int ffLc_calcBOND_force_Cartesian(double *f_b,double *cord,struct AmberParmL ap){
  int i,j,k;
  int numatom;
  int type;
  double f;
  double lenij;
  double atom[2][3];

  numatom=ap.NATOM;
  for (i=0;i<numatom*3;++i) f_b[i] = 0.0;
  
  for (i=0;i<ap.NBONH;++i) {
    type = ap.BH[i][2]-1;
    for (j=0;j<2;++j) for (k=0;k<3;++k) atom[j][k]=cord[abs(ap.BH[i][j])+k];
  
    lenij = len(atom[0],atom[1]);
    for (j=0;j<3;++j) {
      f = -2.0*ap.RK[type]*(lenij-ap.REQ[type])*(atom[1][j]-atom[0][j])/lenij*UNIT;
      f_b[abs(ap.BH[i][0])+j] += f;
      f_b[abs(ap.BH[i][1])+j] += -f;
    }
  }

  for (i=0;i<ap.MBONA;++i) {
    type = ap.BA[i][2]-1;
    for (j=0;j<2;++j) for (k=0;k<3;++k) atom[j][k]=cord[abs(ap.BA[i][j])+k];
  
    lenij = len(atom[0],atom[1]);
    for (j=0;j<3;++j) {
      f = -2.0*ap.RK[type]*(lenij-ap.REQ[type])*(atom[1][j]-atom[0][j])/lenij*UNIT;
      f_b[abs(ap.BA[i][0])+j] += f;
      f_b[abs(ap.BA[i][1])+j] += -f;
    }
  }

  return 0;
}

int ffLc_set_calcff(int numnb, int num14,FILE *inputfile, struct potential *ene,struct AmberParmL ap){
  int i,numatom,f;
  char dummy;

  (*ene).parm.numnb=numnb;
  (*ene).parm.num14=num14;

  numatom=ap.NATOM;
  (*ene).parm.indexnb=(int *)gcemalloc(sizeof(int)*numnb*2);
  (*ene).parm.ele=(double *)gcemalloc(sizeof(double)*numatom);
  (*ene).parm.ALJ=(double *)gcemalloc(sizeof(double)*numatom*numatom);
  (*ene).parm.BLJ=(double *)gcemalloc(sizeof(double)*numatom*numatom);
  (*ene).parm.index14=(int *)gcemalloc(sizeof(int)*num14*2);
  (*ene).p_e=(double *)gcemalloc(sizeof(double)*numnb);
  (*ene).p_LJ=(double *)gcemalloc(sizeof(double)*numnb);
  (*ene).p_e_14=(double *)gcemalloc(sizeof(double)*numnb);
  (*ene).p_LJ_14=(double *)gcemalloc(sizeof(double)*numnb);
  (*ene).p_d=(double *)gcemalloc(sizeof(double)*(ap.NPHIH+ap.MPHIA));
  (*ene).p_a=(double *)gcemalloc(sizeof(double)*(ap.NTHETH+ap.MTHETA));
  (*ene).p_b=(double *)gcemalloc(sizeof(double)*(ap.NBONH+ap.MBONA));

  ffLc_set_NB_PARM((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,numatom,ap);

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

int ffLc_set_calcffsp(struct potential *ene,struct AmberParmL ap){
  int i,numatom,fd;
  int numnb,num14;
  char dummy;

  ffLc_set_non_bonding_index_1(&numnb,&num14,ap);
  (*ene).parm.numnb=numnb;
  (*ene).parm.num14=num14;
  (*ene).parm.indexnb=(int *)gcemalloc(sizeof(int)*numnb*2);
  (*ene).parm.index14=(int *)gcemalloc(sizeof(int)*num14*2);
  ffLc_set_non_bonding_index_2((*ene).parm.indexnb,(*ene).parm.index14,ap);

  numatom=ap.NATOM;
  (*ene).parm.ele=(double *)gcemalloc(sizeof(double)*numatom);
  (*ene).parm.ALJ=(double *)gcemalloc(sizeof(double)*numatom*numatom);
  (*ene).parm.BLJ=(double *)gcemalloc(sizeof(double)*numatom*numatom);
  (*ene).p_e=(double *)gcemalloc(sizeof(double)*numnb);
  (*ene).p_LJ=(double *)gcemalloc(sizeof(double)*numnb);
  (*ene).p_e_14=(double *)gcemalloc(sizeof(double)*numnb);
  (*ene).p_LJ_14=(double *)gcemalloc(sizeof(double)*numnb);
  (*ene).p_d=(double *)gcemalloc(sizeof(double)*(ap.NPHIH+ap.MPHIA));
  (*ene).p_a=(double *)gcemalloc(sizeof(double)*(ap.NTHETH+ap.MTHETA));
  (*ene).p_b=(double *)gcemalloc(sizeof(double)*(ap.NBONH+ap.MBONA));

  ffLc_set_NB_PARM((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,numatom,ap);

}



int ffLc_set_calcffandforce(struct potential *ene, struct force *f,struct AmberParmL ap){
  int i,numatom,fd;
  int numnb,num14;
  char dummy;

  ffLc_set_non_bonding_index_1(&numnb,&num14,ap);
  (*ene).parm.numnb=numnb;
  (*ene).parm.num14=num14;
  (*ene).parm.indexnb=(int *)gcemalloc(sizeof(int)*numnb*2);
  (*ene).parm.index14=(int *)gcemalloc(sizeof(int)*num14*2);
  ffLc_set_non_bonding_index_2((*ene).parm.indexnb,(*ene).parm.index14,ap);

  numatom=ap.NATOM;
  (*ene).parm.ele=(double *)gcemalloc(sizeof(double)*numatom);
  (*ene).parm.ALJ=(double *)gcemalloc(sizeof(double)*numatom*numatom);
  (*ene).parm.BLJ=(double *)gcemalloc(sizeof(double)*numatom*numatom);
  (*ene).p_e=(double *)gcemalloc(sizeof(double)*numnb);
  (*ene).p_LJ=(double *)gcemalloc(sizeof(double)*numnb);
  (*ene).p_e_14=(double *)gcemalloc(sizeof(double)*numnb);
  (*ene).p_LJ_14=(double *)gcemalloc(sizeof(double)*numnb);
  (*ene).p_d=(double *)gcemalloc(sizeof(double)*(ap.NPHIH+ap.MPHIA));
  (*ene).p_a=(double *)gcemalloc(sizeof(double)*(ap.NTHETH+ap.MTHETA));
  (*ene).p_b=(double *)gcemalloc(sizeof(double)*(ap.NBONH+ap.MBONA));

  ffLc_set_NB_PARM((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,numatom,ap);

  (*f).f_t=(double *)gcemalloc(sizeof(double)*numatom*3);
  (*f).f_e=(double *)gcemalloc(sizeof(double)*numatom*3);
  (*f).f_LJ=(double *)gcemalloc(sizeof(double)*numatom*3);
  (*f).f_e_14=(double *)gcemalloc(sizeof(double)*numatom*3);
  (*f).f_LJ_14=(double *)gcemalloc(sizeof(double)*numatom*3);
  (*f).f_d=(double *)gcemalloc(sizeof(double)*numatom*3);
  (*f).f_a=(double *)gcemalloc(sizeof(double)*numatom*3);
  (*f).f_b=(double *)gcemalloc(sizeof(double)*numatom*3);

}

double ffLc_calcff(double *crd, int numatom,struct potential *ene,struct AmberParmL ap) {
  int i;
  int numnb,num14;
  double *f_e,*f_LJ,*n_d;
  double *f_e_14,*f_LJ_14;

  numnb=(*ene).parm.numnb;
  num14=(*ene).parm.num14;

  ffLc_calcDIHE((*ene).p_d,n_d,crd,1,0,0,ap);
  ffLc_calcANGLE((*ene).p_a,crd,ap);
  ffLc_calcBOND((*ene).p_b,crd,ap);


  ffLc_calcFFNB((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e,(*ene).p_LJ,f_e,f_LJ,numnb,(*ene).parm.indexnb,numatom,crd,2,0,ap);
  ffLc_calcFFNB((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e_14,(*ene).p_LJ_14,f_e_14,f_LJ_14,num14,(*ene).parm.index14,numatom,crd,2,0,ap);
  
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
  for (i=0;i<ap.NPHIH+ap.MPHIA;++i) {
    (*ene).p_t+=(*ene).p_d[i];
    (*ene).p_d_t+=(*ene).p_d[i];
  }
  for (i=0;i<ap.NTHETH+ap.MTHETA;++i) {
    (*ene).p_t+=(*ene).p_a[i];
    (*ene).p_a_t+=(*ene).p_a[i];
  }
  for (i=0;i<ap.NBONH+ap.MBONA;++i) {
    (*ene).p_t+=(*ene).p_b[i];
    (*ene).p_b_t+=(*ene).p_b[i];
  }
  
  return (*ene).p_t;
}

double ffLc_calcffandforce(double *crd, int numatom,struct potential *ene,struct force *f,struct AmberParmL ap) {
  int i;
  int numnb,num14;
  double *n_d;

  numnb=(*ene).parm.numnb;
  num14=(*ene).parm.num14;

  ffLc_calcFFNB((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e,(*ene).p_LJ,(*f).f_e,(*f).f_LJ,numnb,(*ene).parm.indexnb,numatom,crd,1,1,ap);
  ffLc_calcFFNB/*_14*/((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e_14,(*ene).p_LJ_14,(*f).f_e_14,(*f).f_LJ_14,num14,(*ene).parm.index14,numatom,crd,1,1,ap);

  ffLc_calcDIHE((*ene).p_d,n_d,crd,1,0,0,ap);
  ffLc_calcDIHE_force_Cartesian((*f).f_d,crd,ap); // 1111

  ffLc_calcANGLE((*ene).p_a,crd,ap);
  ffLc_calcANGLE_force_Cartesian((*f).f_a,crd,ap);

  ffLc_calcBOND((*ene).p_b,crd,ap);
  ffLc_calcBOND_force_Cartesian((*f).f_b,crd,ap);
  
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

  for (i=0;i<ap.NPHIH+ap.MPHIA;++i) {
    (*ene).p_t+=(*ene).p_d[i];
    (*ene).p_d_t+=(*ene).p_d[i];
  }
  for (i=0;i<ap.NTHETH+ap.MTHETA;++i) {
    (*ene).p_t+=(*ene).p_a[i];
    (*ene).p_a_t+=(*ene).p_a[i];             // 0911
  }                                          // 0911
  for (i=0;i<ap.NBONH+ap.MBONA;++i) {        // 0911
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

double ffLc_calcffandforce_w14tune(double *crd, int numatom,struct potential *ene,struct force *f, int *atom_tune_pairs, double *tune_val, int numtune,struct AmberParmL ap) {
  int i;
  int numnb,num14;
  double *n_d;

  numnb=(*ene).parm.numnb;
  num14=(*ene).parm.num14;

  ffLc_calcFFNB((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e,(*ene).p_LJ,(*f).f_e,(*f).f_LJ,numnb,(*ene).parm.indexnb,numatom,crd,1,1,ap);
  ffLc_calcFFNB_wtune/*_14*/((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e_14,(*ene).p_LJ_14,(*f).f_e_14,(*f).f_LJ_14,num14,(*ene).parm.index14,numatom,crd,1,1,atom_tune_pairs,tune_val,numtune,ap);

  ffLc_calcDIHE((*ene).p_d,n_d,crd,1,0,0,ap);
  ffLc_calcDIHE_force_Cartesian((*f).f_d,crd,ap); // 1111

  ffLc_calcANGLE((*ene).p_a,crd,ap);
  ffLc_calcANGLE_force_Cartesian((*f).f_a,crd,ap);

  ffLc_calcBOND((*ene).p_b,crd,ap);
  ffLc_calcBOND_force_Cartesian((*f).f_b,crd,ap);
  
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

  for (i=0;i<ap.NPHIH+ap.MPHIA;++i) {
    (*ene).p_t+=(*ene).p_d[i];
    (*ene).p_d_t+=(*ene).p_d[i];
  }
  for (i=0;i<ap.NTHETH+ap.MTHETA;++i) {
    (*ene).p_t+=(*ene).p_a[i];
    (*ene).p_a_t+=(*ene).p_a[i];             // 0911
  }                                          // 0911
  for (i=0;i<ap.NBONH+ap.MBONA;++i) {        // 0911
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

double ffLc_calcffandforce_w14tuneb(double *crd, int numatom,struct potential *ene,struct force *f, int *atom_tune_pairs_es, int *atom_tune_pairs_LJ, double *tune_val_es,double *tune_val_LJ, int numtune_es, int numtune_LJ,struct AmberParmL ap) {
  int i;
  int numnb,num14;
  double *n_d;

  numnb=(*ene).parm.numnb;
  num14=(*ene).parm.num14;

  ffLc_calcFFNB((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e,(*ene).p_LJ,(*f).f_e,(*f).f_LJ,numnb,(*ene).parm.indexnb,numatom,crd,1,1,ap);
  ffLc_calcFFNB_wtuneb((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e_14,(*ene).p_LJ_14,(*f).f_e_14,(*f).f_LJ_14,num14,(*ene).parm.index14,numatom,crd,1,1,atom_tune_pairs_es,atom_tune_pairs_LJ,tune_val_es,tune_val_LJ,numtune_es,numtune_LJ,ap);

  ffLc_calcDIHE((*ene).p_d,n_d,crd,1,0,0,ap);
  ffLc_calcDIHE_force_Cartesian((*f).f_d,crd,ap); // 1111

  ffLc_calcANGLE((*ene).p_a,crd,ap);
  ffLc_calcANGLE_force_Cartesian((*f).f_a,crd,ap);

  ffLc_calcBOND((*ene).p_b,crd,ap);
  ffLc_calcBOND_force_Cartesian((*f).f_b,crd,ap);

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

  for (i=0;i<ap.NPHIH+ap.MPHIA;++i) {
    (*ene).p_t+=(*ene).p_d[i];
    (*ene).p_d_t+=(*ene).p_d[i];
  }
  for (i=0;i<ap.NTHETH+ap.MTHETA;++i) {
    (*ene).p_t+=(*ene).p_a[i];
    (*ene).p_a_t+=(*ene).p_a[i];             // 0911
  }                                          // 0911
  for (i=0;i<ap.NBONH+ap.MBONA;++i) {        // 0911
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

double ffLc_calcffandforce_14D_woH(double *crd, int numatom,struct potential *ene,struct force *f,struct AmberParmL ap) {
  int i;
  int num14;
  double *n_d;

  num14=(*ene).parm.num14;

  ffLc_calcFFNB_woH((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e_14,(*ene).p_LJ_14,(*f).f_e_14,(*f).f_LJ_14,num14,(*ene).parm.index14,numatom,crd,1,1,ap);

  ffLc_calcDIHE_woH((*ene).p_d,n_d,crd,1,0,0,ap);

  (*ene).p_t=0.0;
  (*ene).p_e_14_t=0.0;
  (*ene).p_LJ_14_t=0.0;
  (*ene).p_d_t=0.0;
  for (i=0;i<numatom*3;++i) (*f).f_t[i]=0.0;

  for (i=ap.NPHIH;i<ap.NPHIH+ap.MPHIA;++i) {
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

double ffLc_calcffandforce_14DAB_woH(double *crd, int numatom,struct potential *ene,struct force *f,struct AmberParmL ap) {
  int i;
  int num14;
  double *n_d;

  num14=(*ene).parm.num14;

  ffLc_calcFFNB_woH((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e_14,(*ene).p_LJ_14,(*f).f_e_14,(*f).f_LJ_14,num14,(*ene).parm.index14,numatom,crd,1,1,ap);

  ffLc_calcDIHE_woH((*ene).p_d,n_d,crd,1,0,0,ap);
  ffLc_calcDIHE_force_Cartesian_woH((*f).f_d,crd,ap);

  /*******************************************************/
  /* ffLc_calcDIHE((*ene).p_d,n_d,crd,1,0,0);		 */
  /* ffLc_calcDIHE_force_Cartesian((*f).f_d,crd); // 1111 */
  /*******************************************************/


  ffLc_calcANGLE_woH((*ene).p_a,crd,ap);
  ffLc_calcANGLE_force_Cartesian_woH((*f).f_a,crd,ap);

  ffLc_calcBOND_woH((*ene).p_b,crd,ap);
  ffLc_calcBOND_force_Cartesian_woH((*f).f_b,crd,ap);

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

  for (i=ap.NPHIH;i<ap.NPHIH+ap.MPHIA;++i) {
    (*ene).p_t+=(*ene).p_d[i];
    (*ene).p_d_t+=(*ene).p_d[i];
  }

  for (i=ap.NTHETH;i<ap.NTHETH+ap.MTHETA;++i) {
    (*ene).p_t+=(*ene).p_a[i];
    (*ene).p_a_t+=(*ene).p_a[i];
  }                                     

  for (i=ap.NBONH;i<ap.NBONH+ap.MBONA;++i) {    
    (*ene).p_t+=(*ene).p_b[i];
    (*ene).p_b_t+=(*ene).p_b[i];
  }

  for (i=0;i<numatom*3;++i) (*f).f_t[i]+=(*f).f_e_14[i]+(*f).f_LJ_14[i]+(*f).f_d[i]+(*f).f_a[i]+(*f).f_b[i];
  
  return (*ene).p_t;
}

double ffLc_calcffandforce_14vdWDAB_woH(double *crd, int numatom,struct potential *ene,struct force *f,struct AmberParmL ap) {
  int i;
  int num14;
  double *n_d;

  num14=(*ene).parm.num14;

  ffLc_calcFFNBvdW_woH((*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_LJ_14,(*f).f_LJ_14,num14,(*ene).parm.index14,numatom,crd,1,1,ap);

  ffLc_calcDIHE_woH((*ene).p_d,n_d,crd,1,0,0,ap);
  ffLc_calcDIHE_force_Cartesian_woH((*f).f_d,crd,ap);

  ffLc_calcANGLE_woH((*ene).p_a,crd,ap);
  ffLc_calcANGLE_force_Cartesian_woH((*f).f_a,crd,ap);

  ffLc_calcBOND_woH((*ene).p_b,crd,ap);
  ffLc_calcBOND_force_Cartesian_woH((*f).f_b,crd,ap);

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

  for (i=ap.NPHIH;i<ap.NPHIH+ap.MPHIA;++i) {
    (*ene).p_t+=(*ene).p_d[i];
    (*ene).p_d_t+=(*ene).p_d[i];
  }

  for (i=ap.NTHETH;i<ap.NTHETH+ap.MTHETA;++i) {
    (*ene).p_t+=(*ene).p_a[i];
    (*ene).p_a_t+=(*ene).p_a[i];
  }                                     

  for (i=ap.NBONH;i<ap.NBONH+ap.MBONA;++i) {    
    (*ene).p_t+=(*ene).p_b[i];
    (*ene).p_b_t+=(*ene).p_b[i];
  }

  for (i=0;i<numatom*3;++i) (*f).f_t[i]+=(*f).f_e_14[i]+(*f).f_LJ_14[i]+(*f).f_d[i]+(*f).f_a[i]+(*f).f_b[i];
  
  return (*ene).p_t;
}



double ffLc_calcffandforce_14_6(double *crd, int numatom,struct potential *ene,struct force *f,double scale,int numrepul,struct AmberParmL ap) { // 0911
  int i;
  int numnb,num14;
  double *n_d;

  numnb=(*ene).parm.numnb;
  num14=(*ene).parm.num14;

  ffLc_calcFFNB_14_6((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e,(*ene).p_LJ,(*f).f_e,(*f).f_LJ,numnb,(*ene).parm.indexnb,numatom,crd,1,1,scale,numrepul,ap);
  ffLc_calcFFNB_14_6((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e_14,(*ene).p_LJ_14,(*f).f_e_14,(*f).f_LJ_14,num14,(*ene).parm.index14,numatom,crd,1,1,scale,numrepul,ap);

  ffLc_calcDIHE((*ene).p_d,n_d,crd,1,0,0,ap);
  ffLc_calcDIHE_force_Cartesian((*f).f_d,crd,ap);

  ffLc_calcANGLE((*ene).p_a,crd,ap);
  ffLc_calcANGLE_force_Cartesian((*f).f_a,crd,ap);

  ffLc_calcBOND((*ene).p_b,crd,ap);
  ffLc_calcBOND_force_Cartesian((*f).f_b,crd,ap);
  
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
  for (i=0;i<ap.NPHIH+ap.MPHIA;++i) {
    (*ene).p_t+=(*ene).p_d[i];
    (*ene).p_d_t+=(*ene).p_d[i];
  }
  for (i=0;i<ap.NTHETH+ap.MTHETA;++i) {
    (*ene).p_t+=(*ene).p_a[i];
    (*ene).p_a_t+=(*ene).p_a[i];             
  }                                          
  for (i=0;i<ap.NBONH+ap.MBONA;++i) {        
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


void ffLc_set_non_bonding_index_1(int *numindexnb, int *numindex14,struct AmberParmL ap) {
  int i,j,k,l;
  int numatom;

  int gnumnb=0;
  int gnum14=0;

  numatom=ap.NATOM;

  gnumnb=0;
  gnum14=0;

  for(i=0;i<numatom;++i) {
    for(j=i+1;j<numatom;++j) {
      if (ffLc_which_calc_nb(i,j,ap) == YES)
	++gnumnb;
      else if(ffLc_which_calc_14_nb(i,j,ap)==YES)
      	++gnum14;
    }
  }

  (*numindexnb)=gnumnb;
  (*numindex14)=gnum14;

}

void ffLc_set_non_bonding_index_2(int *gindexnb,int *gindex14,struct AmberParmL ap) {
  int i,j,k,l;
  int numatom;

  int gnumnb=0;
  int gnum14=0;

  numatom=ap.NATOM;

  k=0;
  l=0;
  for(i=0;i<numatom;++i)
    for(j=i+1;j<numatom;++j)
      if (ffLc_which_calc_nb(i,j,ap) == YES) {
  	gindexnb[k*2]=i;
  	gindexnb[k*2+1]=j;
  	++k;
      }
      else if(ffLc_which_calc_14_nb(i,j,ap)==YES){
  	gindex14[l*2]=i;
  	gindex14[l*2+1]=j;
  	++l;
      }

}

int ffLc_which_calc_nb(int atomi,int atomj,struct AmberParmL ap) {
  int k;
  int sum=0;

  for(k=0;k<atomi;++k) sum+=ap.NUMEX[k];

  for(k=sum;k<sum+ap.NUMEX[atomi];++k){
    if (ap.NATEX[k]-1 == atomj)
      return NO;
  }

  return YES;
}

int ffLc_which_calc_14_nb(int atomi,int atomj,struct AmberParmL ap) {
  int k,l;
  int flag;

  if ((flag=ffLc_which_calc_nb(atomi,atomj,ap))==YES)
    return NO;
  else {
    for(l=0;l<ap.NBONH;++l){
      if ((ap.BH[l][0])/3 == atomi ){
  	if ((ap.BH[l][1])/3 == atomj)
  	  return NO;
      }
      else if ((ap.BH[l][1])/3 == atomi){
  	if ((ap.BH[l][0])/3 == atomj)
  	  return NO;
      }
      else if ((ap.BH[l][0])/3 == atomj){
  	if ((ap.BH[l][1])/3 == atomi)
  	  return NO;
      }
      else if ((ap.BH[l][1])/3 == atomj){
  	if ((ap.BH[l][0])/3 == atomi)
  	  return NO;
      }
    }
  
    for(l=0;l<ap.MBONA;++l){
      if ((ap.BA[l][0])/3 == atomi ){
  	if ((ap.BA[l][1])/3 == atomj)
  	  return NO;
      }
      else if ((ap.BA[l][1])/3 == atomi){
  	if ((ap.BA[l][0])/3 == atomj)
  	  return NO;
      }
      else if ((ap.BA[l][0])/3 == atomj){
  	if ((ap.BA[l][1])/3 == atomi)
  	  return NO;
      }
      else if ((ap.BA[l][1])/3 == atomj){
  	if ((ap.BA[l][0])/3 == atomi)
  	  return NO;
      }
    }
  
    for(l=0;l<ap.NTHETH;++l){
      if ((ap.TH[l][0])/3 == atomi ){
  	if ((ap.TH[l][2])/3 == atomj)
  	  return NO;
      }
      else if ((ap.TH[l][2])/3 == atomi){
  	if ((ap.TH[l][0])/3 == atomj)
  	  return NO;
      }
      else if ((ap.TH[l][0])/3 == atomj){
  	if ((ap.TH[l][2])/3 == atomi)
  	  return NO;
      }
      else if ((ap.TH[l][2])/3 == atomj){
  	if ((ap.TH[l][0])/3 == atomi)
  	  return NO;
      }
    }
  
    for(l=0;l<ap.MTHETA;++l){
      if ((ap.TA[l][0])/3 == atomi ) {
  	if ((ap.TA[l][2])/3 == atomj)
  	  return NO;
      }
      //      else if ((ap.TH[l][2])/3 == atomi){
      else if ((ap.TA[l][2])/3 == atomi){ //1011
  	if ((ap.TA[l][0])/3 == atomj)
  	  return NO;
      }
      else if ((ap.TA[l][0])/3 == atomj){
  	if ((ap.TA[l][2])/3 == atomi)
  	  return NO;
      }
      else if ((ap.TA[l][2])/3 == atomj){
  	if ((ap.TA[l][0])/3 == atomi)
  	  return NO;
      }
    }
  }
  
  return YES;
}


///////////////////////////////////////////////////////////////////

double ffLc_calcffandforce_woH(double *crd, int numatom,struct potential *ene,struct force *f,struct AmberParmL ap) {
  int i;
  int numnb,num14;
  double *n_d;

  numnb=(*ene).parm.numnb;
  num14=(*ene).parm.num14;

  ffLc_calcFFNB_woH((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e,(*ene).p_LJ,(*f).f_e,(*f).f_LJ,numnb,(*ene).parm.indexnb,numatom,crd,1,1,ap);
  ffLc_calcFFNB_woH((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e_14,(*ene).p_LJ_14,(*f).f_e_14,(*f).f_LJ_14,num14,(*ene).parm.index14,numatom,crd,1,1,ap);

  ffLc_calcDIHE_woH((*ene).p_d,n_d,crd,1,0,0,ap);
  ffLc_calcDIHE_force_Cartesian_woH((*f).f_d,crd,ap);

  ffLc_calcANGLE_woH((*ene).p_a,crd,ap);
  ffLc_calcANGLE_force_Cartesian_woH((*f).f_a,crd,ap);

  ffLc_calcBOND_woH/*_wFC100*/((*ene).p_b,crd,ap);
  ffLc_calcBOND_force_Cartesian_woH/*_wFC100*/((*f).f_b,crd,ap);
  
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
  for (i=0;i<ap.NPHIH+ap.MPHIA;++i) {
    (*ene).p_t+=(*ene).p_d[i];
    (*ene).p_d_t+=(*ene).p_d[i];
  }
  for (i=0;i<ap.NTHETH+ap.MTHETA;++i) {
    (*ene).p_t+=(*ene).p_a[i];
    (*ene).p_a_t+=(*ene).p_a[i];
  }                                      // 0911
  for (i=0;i<ap.NBONH+ap.MBONA;++i) {    // 0911
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

int ffLc_calcDIHE_woH(double *p_d,
		double *n_d,
		double *cord,
		int flagp, int flagf, int flaginp,struct AmberParmL ap) {
  int i,j,k,l,flag;
  int dtype;
  
  double atom[4][3];
  double dihedang;
  
  for (i=0;i<ap.NPHIH;++i){
    p_d[i] = 0.0;
    if (flagf==1)
      n_d[i] = 0.0;
  }
  
  for (i=0;i<ap.MPHIA;++i){
    p_d[i+ap.NPHIH] = 0.0;
    if (flagf==1)
      n_d[i+ap.NPHIH] = 0.0;
  }
  
  for (i=0;i<ap.MPHIA;++i) {
    dtype = ap.PA[i][4]-1;
    flag=1;
    if (flag==1 || flaginp==0) {
      for (j=0;j<4;++j) {
  	for (k=0;k<3;++k) {
  	  atom[j][k]=cord[abs(ap.PA[i][j])+k];
  	}
      }
  
      dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
      p_d[i+ap.NPHIH] = ap.PK[dtype]*(1.0+cos(ap.PN[dtype]*dihedang-ap.PHASE[dtype]));
      if (flagf==1)
  	n_d[i+ap.NPHIH] = -ap.PK[dtype]*4.18407*100.0*(sin(ap.PN[dtype]*dihedang-ap.PHASE[dtype])*ap.PN[dtype]);
    }
  }
}

int ffLc_calcDIHE_force_Cartesian_woH(double *f_d,double *cord,struct AmberParmL ap) {
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

  //  for (i=0;i<(ap.NPHIH+ap.MPHIA)*3;++i) f_d[i] = 0.0;
  for (i=0;i<ap.NATOM*3;++i) f_d[i] = 0.0;
  
  for (i=0;i<ap.MPHIA;++i) {
    dtype = ap.PA[i][4]-1;
    for (j=0;j<4;++j) for (k=0;k<3;++k)	atom[j][k]=cord[abs(ap.PA[i][j])+k];

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

    if (ap.PN[dtype]==1) fa=-ap.PK[dtype]*ap.PN[dtype]*cos(ap.PHASE[dtype]);
    else if (ap.PN[dtype]==2) fa=-ap.PK[dtype]*ap.PN[dtype]*2.0*cosdih*cos(ap.PHASE[dtype]);
    else if (ap.PN[dtype]==3) fa=-ap.PK[dtype]*ap.PN[dtype]*(-4.0*sindih*sindih+3.0)*cos(ap.PHASE[dtype]);
    else if (ap.PN[dtype]==4) fa=-ap.PK[dtype]*ap.PN[dtype]*4.0*(cosdih*(2.0*cosdih*cosdih-1.0))*cos(ap.PHASE[dtype]);
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
      f_d[abs(ap.PA[i][0])+j] += fa*op1[j]*UNIT;
      f_d[abs(ap.PA[i][1])+j] += fa*(-op2[j]+op3[j])*UNIT;
      f_d[abs(ap.PA[i][2])+j] += fa*(-op4[j]+op5[j])*UNIT;
      f_d[abs(ap.PA[i][3])+j] += fa*op6[j]*UNIT;
    }
  }
}

int ffLc_calcffandforce_speD(double *f_d,double *p_d,double *cord,int numspedihed,
			    int *atom1,int *atom2,int *atom3,int *atom4,struct AmberParmL ap) {
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

  //  for (i=0;i<(ap.NPHIH+ap.MPHIA)*3;++i) f_d[i] = 0.0;
  for (i=0;i<4*3;++i) f_d[i] = 0.0;
  
  dtype = ap.PA[numspedihed][4]-1;
  for (j=0;j<4;++j) for (k=0;k<3;++k)	atom[j][k]=cord[abs(ap.PA[numspedihed][j])+k];
  
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
  *p_d = ap.PK[dtype]*(1.0+cos(ap.PN[dtype]*dihedang-ap.PHASE[dtype]));

  if (ap.PN[dtype]==1) fa=-ap.PK[dtype]*ap.PN[dtype]*cos(ap.PHASE[dtype]);
  else if (ap.PN[dtype]==2) fa=-ap.PK[dtype]*ap.PN[dtype]*2.0*cosdih*cos(ap.PHASE[dtype]);
  else if (ap.PN[dtype]==3) fa=-ap.PK[dtype]*ap.PN[dtype]*(-4.0*sindih*sindih+3.0)*cos(ap.PHASE[dtype]);
  else if (ap.PN[dtype]==4) fa=-ap.PK[dtype]*ap.PN[dtype]*4.0*(cosdih*(2.0*cosdih*cosdih-1.0))*cos(ap.PHASE[dtype]);
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
  *atom1=ap.PA[numspedihed][0];
  *atom2=ap.PA[numspedihed][1];
  *atom3=ap.PA[numspedihed][2];
  *atom4=ap.PA[numspedihed][3];

}

int ffLc_calcANGLE_woH(double *p_a,double *cord,struct AmberParmL ap){
  int i,j,k,l;
  int type;
  
  double atom[3][3];
  double ang;
  
  for (i=0;i<ap.NTHETH+ap.MTHETA;++i)
    p_a[i] = 0.0;
  
  for (i=0;i<ap.MTHETA;++i) {
    type = ap.TA[i][3]-1;
     for (j=0;j<3;++j)
      for (k=0;k<3;++k)
	atom[j][k]=cord[abs(ap.TA[i][j])+k];
  
    ang = pick_angle(atom[0],atom[1],atom[2],0,0.0);
    p_a[i+ap.NTHETH] = ap.TK[type]*(ang-ap.TEQ[type])*(ang-ap.TEQ[type]);
  }

  return 0;
}

int ffLc_calcANGLE_force_Cartesian_woH(double *f_a,double *cord,struct AmberParmL ap){
  int i,j,k,l;
  int numatom;
  int type;
  double kang,ang_eq;

  double atom[3][3];
  double *f_temp;

  numatom=ap.NATOM;
  for (i=0;i<numatom*3;++i) f_a[i] = 0.0;
  f_temp=(double *)gcemalloc(sizeof(double)*3*3);

  for (i=0;i<ap.MTHETA;++i) {
    for (j=0;j<3;++j)  for (k=0;k<3;++k) atom[j][k]=cord[abs(ap.TA[i][j])+k];
    type = ap.TA[i][3]-1;
    kang = ap.TK[type];
    ang_eq = ap.TEQ[type];
    calcANGKE_force(atom[0],atom[1],atom[2],kang,ang_eq,f_temp);
    for (j=0;j<3;++j) {
      f_a[abs(ap.TA[i][0])+j] += f_temp[j];
      f_a[abs(ap.TA[i][1])+j] += f_temp[3+j];
      f_a[abs(ap.TA[i][2])+j] += f_temp[6+j];
    }
  }

  return 0;
}

int ffLc_calcBOND_woH(double *p_b,double *cord,struct AmberParmL ap){
  int i,j,k;
  int type;
  double len;
  double atom[2][3];

  for (i=0;i<ap.NBONH+ap.MBONA;++i)
    p_b[i] = 0.0;
  
  for (i=0;i<ap.MBONA;++i) {
    type = ap.BA[i][2]-1;
    for (j=0;j<2;++j)
      for (k=0;k<3;++k)
	atom[j][k]=cord[abs(ap.BA[i][j])+k];  

    len = pick_bond_leng(atom[0],atom[1]);
    p_b[i+ap.NBONH] = ap.RK[type]*(len-ap.REQ[type])*(len-ap.REQ[type]);
  }

  return 0;
}

int ffLc_calcBOND_woH_wFC100(double *p_b,double *cord,struct AmberParmL ap){
  ;
}


int ffLc_calcBOND_force_Cartesian_woH(double *f_b,double *cord,struct AmberParmL ap){
  int i,j,k;
  int numatom;
  int type;
  double f;
  double lenij;
  double atom[2][3];

  numatom=ap.NATOM;
  for (i=0;i<numatom*3;++i) f_b[i] = 0.0;
  
  for (i=0;i<ap.MBONA;++i) {
    type = ap.BA[i][2]-1;
    for (j=0;j<2;++j) for (k=0;k<3;++k) atom[j][k]=cord[abs(ap.BA[i][j])+k];
  
    lenij = len(atom[0],atom[1]);
    for (j=0;j<3;++j) {
      f = -2.0*ap.RK[type]*(lenij-ap.REQ[type])*(atom[1][j]-atom[0][j])/lenij*UNIT;
      f_b[abs(ap.BA[i][0])+j] += f;
      f_b[abs(ap.BA[i][1])+j] += -f;
    }
  }

  return 0;
}

int ffLc_calcBOND_force_Cartesian_woH_wFC100(double *f_b,double *cord,struct AmberParmL ap){
  ;
}


int ffLc_calcFFNB_woH(double *ele, double *ALJ, double *BLJ,
		    double *p_e,double *p_LJ,
		    double *f_e,double *f_LJ,
		    int numnb, int *indexnb,
		    int num_atom,
		    double *cord,
		    int flagp, int flagf,struct AmberParmL ap) {
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
    if (strncmp(ap.IGRAPH[num_a_prot],"H",1)!=0 && strncmp(ap.IGRAPH[NUM_A_PROT],"H",1)!=0 ) {
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

int ffLc_calcFFNBvdW_woH( double *ALJ, double *BLJ,
			 double *p_LJ,
			 double *f_LJ,
			 int numnb, int *indexnb,
			 int num_atom,
			 double *cord,
			 int flagp, int flagf,struct AmberParmL ap) {
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
    if (strncmp(ap.IGRAPH[num_a_prot],"H",1)!=0 && strncmp(ap.IGRAPH[NUM_A_PROT],"H",1)!=0 ) {
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

int ffLc_calcTorque(double *Q,double *cord,int numclut,int *nNumClutOfParent,int *terminal_atom_a,int *origin_atom_a,struct AmberParmL ap) {
  int i,j,k,l,n;
  int dtype,flag;  
  double atom[4][3];
  double dihedang;
  int *inpindex,inpnumH,inpnumA,*indexclut;
  int ON=1,OFF=0;

  indexclut=(int *)gcemalloc(sizeof(int)*(ap.NPHIH+ap.MPHIA));
  //  inpindex=(int *)gcemalloc(sizeof(int)*1);
  inpindex=ffLc_make_inpindex(&inpnumH,&inpnumA,indexclut,numclut,nNumClutOfParent,terminal_atom_a,origin_atom_a,ap);

  for (i=0;i<ap.NPHIH;++i) {
    flag=ON;
    for (j=0;j<inpnumH;++j) {
      if (i == inpindex[j]) {
	flag=OFF;
	break;
      }
    }

    if (flag==ON) {
      dtype = ap.PH[i][4]-1;
      n = indexclut[i]-1;
      for (j=0;j<4;++j) for (k=0;k<3;++k) atom[j][k]=cord[abs(ap.PH[i][j])+k];
      dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
      //      Q[n] += -ap.PK[dtype]*4.18407*100.0*(sin(ap.PN[dtype]*dihedang-ap.PHASE[dtype])*ap.PN[dtype]);
      Q[n] += -ap.PK[dtype]*(sin(ap.PN[dtype]*dihedang-ap.PHASE[dtype])*ap.PN[dtype]);
    }
  }

  for (i=0;i<numclut;++i) {
    Q[i] = 4.18407*100.0*Q[i];
  } 

  for (i=0;i<ap.MPHIA;++i) {
    flag=ON;
    for (j=0;j<inpnumA;++j) {
      if (i == inpindex[j+inpnumH]) {
	flag=OFF;
	break;
      }
    }

    if (flag==ON) {
      dtype = ap.PA[i][4]-1;
      n = indexclut[i+ap.NPHIH]-1;
      for (j=0;j<4;++j) for (k=0;k<3;++k) atom[j][k]=cord[abs(ap.PA[i][j])+k];
  
      dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
      Q[n] += -ap.PK[dtype]*4.18407*100.0*(sin(ap.PN[dtype]*dihedang-ap.PHASE[dtype])*ap.PN[dtype]);
    }
  }
}

double ffLc_calcTorque_wtune(double *Q,double *p_d,double *cord,int numclut,int *nNumClutOfParent,int *terminal_atom_a,int *origin_atom_a,int *atom_tune_pairs, double *tune_val_by_period, int numtune,struct AmberParmL ap) {
  int i,j,k,l,n;
  int dtype,flag;  
  double atom[4][3];
  double dihedang;
  int *inpindex,inpnumH,inpnumA,*indexclut;
  int ON=1,OFF=0;

  int tuneflag,tuneindex;
  double p_d_t=0.0;

  indexclut=(int *)gcemalloc(sizeof(int)*(ap.NPHIH+ap.MPHIA));
  //  inpindex=(int *)gcemalloc(sizeof(int)*1);
  inpindex=ffLc_make_inpindex(&inpnumH,&inpnumA,indexclut,numclut,nNumClutOfParent,terminal_atom_a,origin_atom_a,ap);

  for (i=0;i<ap.NPHIH;++i) {
    flag=ON;
    for (j=0;j<inpnumH;++j) {
      if (i == inpindex[j]) {
	flag=OFF;
	break;
      }
    }

    if (flag==ON) {
      dtype = ap.PH[i][4]-1;
      n = indexclut[i]-1;
      for (j=0;j<4;++j) for (k=0;k<3;++k) atom[j][k]=cord[abs(ap.PH[i][j])+k];

      tuneflag=OFF;
      for (j=0;j<numtune;++j) {
	if (atom_tune_pairs[j*4+0]==abs(ap.PH[i][0])/3 && atom_tune_pairs[j*4+1]==abs(ap.PH[i][1])/3 
	    && atom_tune_pairs[j*4+2]==abs(ap.PH[i][2])/3 && atom_tune_pairs[j*4+3]==abs(ap.PH[i][3])/3) {
	  tuneindex=j;
	  tuneflag=ON;
	  break;
	}
      }

      dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
      //      Q[n] += -ap.PK[dtype]*4.18407*100.0*(sin(ap.PN[dtype]*dihedang-ap.PHASE[dtype])*ap.PN[dtype]);
      if (tuneflag==OFF) {
	p_d[i] = ap.PK[dtype]*(1.0+cos(ap.PN[dtype]*dihedang-ap.PHASE[dtype]));
	p_d_t+=p_d[i];
	Q[n] += -ap.PK[dtype]*(sin(ap.PN[dtype]*dihedang-ap.PHASE[dtype])*ap.PN[dtype]);
      }
      else {
	p_d[i] = 1.0/(tune_val_by_period[tuneindex]*tune_val_by_period[tuneindex])
	  *ap.PK[dtype]*(1.0+cos(ap.PN[dtype]*dihedang-ap.PHASE[dtype]));
	p_d_t+=p_d[i];
	Q[n] += -1.0/(tune_val_by_period[tuneindex]*tune_val_by_period[tuneindex])
	  *ap.PK[dtype]*(sin(ap.PN[dtype]*dihedang-ap.PHASE[dtype])*ap.PN[dtype]);
      }
    }
  }

  for (i=0;i<numclut;++i) {
    Q[i] = 4.18407*100.0*Q[i];
  } 

  for (i=0;i<ap.MPHIA;++i) {
    flag=ON;
    for (j=0;j<inpnumA;++j) {
      if (i == inpindex[j+inpnumH]) {
	flag=OFF;
	break;
      }
    }

    if (flag==ON) {
      dtype = ap.PA[i][4]-1;
      n = indexclut[i+ap.NPHIH]-1;
      for (j=0;j<4;++j) for (k=0;k<3;++k) atom[j][k]=cord[abs(ap.PA[i][j])+k];

      tuneflag=OFF;
      for (j=0;j<numtune;++j) {
	if (atom_tune_pairs[j*4+0]==abs(ap.PA[i][0])/3 && atom_tune_pairs[j*4+1]==abs(ap.PA[i][1])/3 
	    && atom_tune_pairs[j*4+2]==abs(ap.PA[i][2])/3 && atom_tune_pairs[j*4+3]==abs(ap.PA[i][3])/3) {
	  tuneindex=j;
	  tuneflag=ON;
	  break;
	}
      }
  
      dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
      if (tuneflag==OFF) {
	p_d[i+ap.NPHIH] = ap.PK[dtype]*(1.0+cos(ap.PN[dtype]*dihedang-ap.PHASE[dtype]));
	p_d_t+=p_d[i];
	Q[n] += -ap.PK[dtype]*4.18407*100.0*(sin(ap.PN[dtype]*dihedang-ap.PHASE[dtype])*ap.PN[dtype]);
      }
      else { 
	p_d[i+ap.NPHIH] = 1.0/(tune_val_by_period[tuneindex]*tune_val_by_period[tuneindex])
	  *ap.PK[dtype]*(1.0+cos(ap.PN[dtype]*dihedang-ap.PHASE[dtype]));
	p_d_t+=p_d[i];
	Q[n] += -1.0/(tune_val_by_period[tuneindex]*tune_val_by_period[tuneindex])
	  *ap.PK[dtype]*4.18407*100.0*(sin(ap.PN[dtype]*dihedang-ap.PHASE[dtype])*ap.PN[dtype]);
      }
    }
  }

  return p_d_t;
}

double ffLc_calcTorque_wtuneb(double *Q,double *p_d,double *cord,int numclut,int *nNumClutOfParent,int *terminal_atom_a,int *origin_atom_a,int *atom_tune_pairs, double *tune_val_V_n, double *tune_val_n_phase, int numtune, double pi,struct AmberParmL ap) {
  int i,j,k,l,n;
  int dtype,flag;  
  double atom[4][3];
  double dihedang;
  int *inpindex,inpnumH,inpnumA,*indexclut;
  int ON=1,OFF=0;

  int tuneflag,tuneindex;
  double p_d_t=0.0;

  double N,PHASE;

  indexclut=(int *)gcemalloc(sizeof(int)*(ap.NPHIH+ap.MPHIA));
  //  inpindex=(int *)gcemalloc(sizeof(int)*1);
  inpindex=ffLc_make_inpindex(&inpnumH,&inpnumA,indexclut,numclut,nNumClutOfParent,terminal_atom_a,origin_atom_a,ap);

  for (i=0;i<ap.NPHIH;++i) {
    flag=ON;
    for (j=0;j<inpnumH;++j) {
      if (i == inpindex[j]) {
	flag=OFF;
	break;
      }
    }

    if (flag==ON) {
      dtype = ap.PH[i][4]-1;
      n = indexclut[i]-1;
      for (j=0;j<4;++j) for (k=0;k<3;++k) atom[j][k]=cord[abs(ap.PH[i][j])+k];

      tuneflag=OFF;
      for (j=0;j<numtune;++j) {
	if (atom_tune_pairs[j*4+0]==abs(ap.PH[i][0])/3 && atom_tune_pairs[j*4+1]==abs(ap.PH[i][1])/3 
	    && atom_tune_pairs[j*4+2]==abs(ap.PH[i][2])/3 && atom_tune_pairs[j*4+3]==abs(ap.PH[i][3])/3) {
	  tuneindex=j;
	  tuneflag=ON;
	  break;
	}
      }

      dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
      //      Q[n] += -ap.PK[dtype]*4.18407*100.0*(sin(ap.PN[dtype]*dihedang-ap.PHASE[dtype])*ap.PN[dtype]);
      if (tuneflag==OFF) {
	p_d[i] = ap.PK[dtype]*(1.0+cos(ap.PN[dtype]*dihedang-ap.PHASE[dtype]));
	p_d_t+=p_d[i];
	Q[n] += -ap.PK[dtype]*(sin(ap.PN[dtype]*dihedang-ap.PHASE[dtype])*ap.PN[dtype]);
      }
      else {
	N=tune_val_n_phase[tuneindex]*ap.PN[dtype];
	PHASE=2.0*pi/N;
 	if (PHASE >= 2.0*pi) PHASE-=2.0*pi;
	p_d[i] = tune_val_V_n[tuneindex]*ap.PK[dtype]*(1.0+cos(N*dihedang-PHASE));
	p_d_t+=p_d[i];
	Q[n] += -tune_val_V_n[tuneindex]*ap.PK[dtype]*(sin(N*dihedang-PHASE)*N);
      }
    }
  }

  for (i=0;i<numclut;++i) {
    Q[i] = 4.18407*100.0*Q[i];
  } 

  for (i=0;i<ap.MPHIA;++i) {
    flag=ON;
    for (j=0;j<inpnumA;++j) {
      if (i == inpindex[j+inpnumH]) {
	flag=OFF;
	break;
      }
    }

    if (flag==ON) {
      dtype = ap.PA[i][4]-1;
      n = indexclut[i+ap.NPHIH]-1;
      for (j=0;j<4;++j) for (k=0;k<3;++k) atom[j][k]=cord[abs(ap.PA[i][j])+k];

      tuneflag=OFF;
      for (j=0;j<numtune;++j) {
	if (atom_tune_pairs[j*4+0]==abs(ap.PA[i][0])/3 && atom_tune_pairs[j*4+1]==abs(ap.PA[i][1])/3 
	    && atom_tune_pairs[j*4+2]==abs(ap.PA[i][2])/3 && atom_tune_pairs[j*4+3]==abs(ap.PA[i][3])/3) {
	  tuneindex=j;
	  tuneflag=ON;
	  break;
	}
      }
  
      dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
      if (tuneflag==OFF) {
	p_d[i+ap.NPHIH] = ap.PK[dtype]*(1.0+cos(ap.PN[dtype]*dihedang-ap.PHASE[dtype]));
	p_d_t+=p_d[i];
	Q[n] += -ap.PK[dtype]*4.18407*100.0*(sin(ap.PN[dtype]*dihedang-ap.PHASE[dtype])*ap.PN[dtype]);
      }
      else { 
	N=tune_val_n_phase[tuneindex]*ap.PN[dtype];
	PHASE=2.0*pi/N;
 	if (PHASE >= 2.0*pi) PHASE-=2.0*pi;
	p_d[i] = tune_val_V_n[tuneindex]*ap.PK[dtype]*(1.0+cos(N*dihedang-PHASE));
	p_d_t+=p_d[i];
	Q[n] += -tune_val_V_n[tuneindex]*ap.PK[dtype]*(sin(N*dihedang-PHASE)*N);
      }
    }
  }

  return p_d_t;
}



int ffLc_calcTorque_woH(double *Q,double *cord,int numclut,int *nNumClutOfParent,int *terminal_atom_a,int *origin_atom_a,struct AmberParmL ap) {
  int i,j,k,l,n;
  int dtype,flag;  
  double atom[4][3];
  double dihedang;
  int *inpindex,inpnumH,inpnumA,*indexclut;
  int ON=1,OFF=0;

  indexclut=(int *)gcemalloc(sizeof(int)*(ap.NPHIH+ap.MPHIA));
  //  inpindex=(int *)gcemalloc(sizeof(int)*1);
  inpindex=ffLc_make_inpindex(&inpnumH,&inpnumA,indexclut,numclut,nNumClutOfParent,terminal_atom_a,origin_atom_a,ap);

  for (i=0;i<ap.MPHIA;++i) {
    flag=ON;
    for (j=0;j<inpnumA;++j) {
      if (i == inpindex[j+inpnumH]) {
	flag=OFF;
	break;
      }
    }

    if (flag==ON) {
      dtype = ap.PA[i][4]-1;
      n = indexclut[i+ap.NPHIH]-1;
      for (j=0;j<4;++j) for (k=0;k<3;++k) atom[j][k]=cord[abs(ap.PA[i][j])+k];
  
      dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
      Q[n] += -ap.PK[dtype]*4.18407*100.0*(sin(ap.PN[dtype]*dihedang-ap.PHASE[dtype])*ap.PN[dtype]);
    }
  }
}



int *ffLc_make_inpindex(int *inpnumH,int *inpnumA,int *indexclut,int numclut,int *nNumClutOfParent,int *terminal_atom_a,int *origin_atom_a,struct AmberParmL ap) {
  int i,j,k,l,ll,p;
  int flag;
  int *atom_dihed_pair;
  int ON=1,OFF=0;
  int *inpindexH,*inpindexA;
  int *inpindex;

  (*inpnumH)=0;
  (*inpnumA)=0;
  
  atom_dihed_pair=(int *)gcemalloc(sizeof(int)*(ap.NPHIH+ap.MPHIA)*6);
  inpindexH=(int *)gcemalloc(sizeof(int)*1);
  inpindexA=(int *)gcemalloc(sizeof(int)*1);
  
  for (k=0;k<ap.NPHIH;++k) {
    flag=OFF;
    for (i=0;i<ap.NBONH;++i) {
      if (((ap.BH[i][0] == abs(ap.PH[k][1]) && ap.BH[i][1] == abs(ap.PH[k][0]))
  	   || (ap.BH[i][0] == abs(ap.PH[k][0]) && ap.BH[i][1] == abs(ap.PH[k][1])))) {
  	flag = ON;
  	break;
      }
    }
    if (flag==OFF) {
      for (i=0;i<ap.NBONA;++i) {
  	if (((ap.BA[i][0] == abs(ap.PH[k][1]) && ap.BA[i][1] == abs(ap.PH[k][0]))
  	     || (ap.BA[i][0] == abs(ap.PH[k][0]) && ap.BA[i][1] == abs(ap.PH[k][1])))) {
  	  flag = ON;
  	  break;
  	}
      }
    }
    if (flag==ON) {
      flag=OFF;
      for (i=0;i<ap.NBONH;++i) {
  	if (((ap.BH[i][0] == abs(ap.PH[k][2]) && ap.BH[i][1] == abs(ap.PH[k][3]))
  	     || (ap.BH[i][0] == abs(ap.PH[k][3]) && ap.BH[i][1] == abs(ap.PH[k][2]))) ) {
  	  flag = ON;
  	  break;
  	}
      }
      if (flag==OFF ) {
  	for (i=0;i<ap.NBONA;++i) {
  	  if (((ap.BA[i][0] == abs(ap.PH[k][2]) && ap.BA[i][1] == abs(ap.PH[k][3]))
  	       || (ap.BA[i][0] == abs(ap.PH[k][3]) && ap.BA[i][1] == abs(ap.PH[k][2]))) ) {
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
      atom_dihed_pair[k/*][*/*6+i] = abs(ap.PH[k][i])/3+1;
    }
    atom_dihed_pair[k/*][*/*6+4] = ap.PH[k][i];
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
  
  for (k=0;k<ap.MPHIA;++k) {
    // check for improper dihed
    flag=OFF;
    for (i=0;i<ap.NBONH;++i) {
      if (   ((ap.BH[i][0] == abs(ap.PA[k][1]) && ap.BH[i][1] == abs(ap.PA[k][0])) || (ap.BH[i][0] == abs(ap.PA[k][0]) && ap.BH[i][1] == abs(ap.PA[k][1])))) {
  	flag = ON;
  	break;
      }
    }
    if (flag==OFF) {
      for (i=0;i<ap.NBONA;++i) {
  	if (   ((ap.BA[i][0] == abs(ap.PA[k][1]) && ap.BA[i][1] == abs(ap.PA[k][0])) || (ap.BA[i][0] == abs(ap.PA[k][0]) && ap.BA[i][1] == abs(ap.PA[k][1])))) {
  	  flag = ON;
  	  break;
  	}
      }
    }
    if (flag==ON) {
      flag=OFF;
      for (i=0;i<ap.NBONH;++i) {
  	if (((ap.BH[i][0] == abs(ap.PA[k][2]) && ap.BH[i][1] == abs(ap.PA[k][3])) || (ap.BH[i][0] == abs(ap.PA[k][3]) && ap.BH[i][1] == abs(ap.PA[k][2]))) ) {
  	  flag = ON;
  	  break;
  	}
      }
      if (flag==OFF ) {
  	for (i=0;i<ap.NBONA;++i) {
  	  if (((ap.BA[i][0] == abs(ap.PA[k][2]) && ap.BA[i][1] == abs(ap.PA[k][3])) || (ap.BA[i][0] == abs(ap.PA[k][3]) && ap.BA[i][1] == abs(ap.PA[k][2]))) ) {
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
      atom_dihed_pair[(k+ap.NPHIH)/*][*/*6+i] = abs(ap.PA[k][i])/3+1;
    }
    atom_dihed_pair[(k+ap.NPHIH)/*][*/*6+4] = ap.PA[k][4];
    l = 1;
    ll = 2;
    if (atom_dihed_pair[(k+ap.NPHIH)/*][*/*6+1] > atom_dihed_pair[(k+ap.NPHIH)/*][*/*6+2]) {
      l = 2;
      ll = 1;
    }
    for (i=0;i<numclut;++i){
      p = nNumClutOfParent[i]-1;
      if (atom_dihed_pair[(k+ap.NPHIH)/*][*/*6+l]==terminal_atom_a[p] && atom_dihed_pair[(k+ap.NPHIH)/*][*/*6+ll]==origin_atom_a[i]){
  
  	atom_dihed_pair[(k+ap.NPHIH)/*][*/*6+5] = i+1;
      }
    }
  }
  
  for (i=0;i<ap.NPHIH+ap.MPHIA;++i) {
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

int *ffLc_make_inpindex_A(int *inpnumA,struct AmberParmL ap) {
  int i,j,k,l,ll,p;
  int flag;
  int ON=1,OFF=0;
  int *inpindexA;

  (*inpnumA)=0;
  
  inpindexA=(int *)gcemalloc(sizeof(int)*ap.MPHIA);

  for (k=0;k<ap.MPHIA;++k) {
    // check for improper dihed
    flag=OFF;
    for (i=0;i<ap.NBONH;++i) {
      if ( ((ap.BH[i][0] == abs(ap.PA[k][1]) && ap.BH[i][1] == abs(ap.PA[k][0])) || (ap.BH[i][0] == abs(ap.PA[k][0]) && ap.BH[i][1] == abs(ap.PA[k][1])))) {
  	flag = ON;
  	break;
      }
    }
    if (flag==OFF) {
      for (i=0;i<ap.NBONA;++i) {
  	if (((ap.BA[i][0] == abs(ap.PA[k][1]) && ap.BA[i][1] == abs(ap.PA[k][0])) || (ap.BA[i][0] == abs(ap.PA[k][0]) && ap.BA[i][1] == abs(ap.PA[k][1])))) {
  	  flag = ON;
  	  break;
  	}
      }
    }
    if (flag==ON) {
      flag=OFF;
      for (i=0;i<ap.NBONH;++i) {
  	if (((ap.BH[i][0] == abs(ap.PA[k][2]) && ap.BH[i][1] == abs(ap.PA[k][3])) || (ap.BH[i][0] == abs(ap.PA[k][3]) && ap.BH[i][1] == abs(ap.PA[k][2]))) ) {
  	  flag = ON;
  	  break;
  	}
      }
      if (flag==OFF ) {
  	for (i=0;i<ap.NBONA;++i) {
  	  if (((ap.BA[i][0] == abs(ap.PA[k][2]) && ap.BA[i][1] == abs(ap.PA[k][3])) || (ap.BA[i][0] == abs(ap.PA[k][3]) && ap.BA[i][1] == abs(ap.PA[k][2]))) ) {
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


void ffLc_out_formated(FILE *outputfile,struct potential e,double KE,double KEv,double PEv,double T,int i,double dt,struct AmberParmL ap) {

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
