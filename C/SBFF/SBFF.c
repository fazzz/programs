
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SBFF.h"
#include "FF.h"
#include "MB.h"
#include "PT.h"
#include "LA.h"
#include "TOPO.h"
#include "mymath.h"
#include "EF.h"

#define YES 1
#define NO  0

#define ON 1
#define OFF 0

double SBAAff_calcff(double *crd, int numatom,struct potential_SBAA *ene) {
  int i;

  SBAAff_calcDIHE((*ene).p_d,crd,(*ene).DEQ,(*ene).EBB,(*ene).ESC,(*ene).EI);
  SBAAff_calcANGLE((*ene).p_a,crd,(*ene).EA,(*ene).AEQ);
  SBAAff_calcBOND((*ene).p_b,crd,(*ene).EB,(*ene).BEQ);
  
  SBAAff_calcFFCNB((*ene).ALJ,(*ene).BLJ,(*ene).p_cnb,(*ene).numcnb,(*ene).indexcnb,numatom,crd);
  SBAAff_calcFFNNB((*ene).ALJnnb,(*ene).p_nnb,(*ene).numnnb,(*ene).indexnnb,numatom,crd);
  (*ene).p_t=0.0;
  (*ene).p_cnb_t=0.0;
  (*ene).p_nnb_t=0.0;
  (*ene).p_d_t=0.0;
  (*ene).p_a_t=0.0;
  (*ene).p_b_t=0.0;
    
  for (i=0;i<numatom;++i) {
    (*ene).p_t+=(*ene).p_cnb[i];
    (*ene).p_cnb_t+=(*ene).p_cnb[i];
  }
  for (i=0;i<numatom;++i) {
    (*ene).p_t+=(*ene).p_nnb[i];
    (*ene).p_nnb_t+=(*ene).p_nnb[i];
  }
  for (i=0;i<AP.MPHIA;++i) {
    (*ene).p_t+=(*ene).p_d[i];
    (*ene).p_d_t+=(*ene).p_d[i];
  }
  for (i=0;i<AP.MTHETA;++i) {
    (*ene).p_t+=(*ene).p_a[i];
    (*ene).p_a_t+=(*ene).p_a[i];
  }
  for (i=0;i<AP.MBONA;++i) {
    (*ene).p_t+=(*ene).p_b[i];
    (*ene).p_b_t+=(*ene).p_b[i];
  }
  
  return (*ene).p_t;
}

int SBAAff_set_calcff(struct potential_SBAA *ene, double *cord,int numatom){
  int i;
  double *indexcnb_cradii;
    
  (*ene).BEQ=(double *)gcemalloc(sizeof(double)*AP.MBONA);
  (*ene).AEQ=(double *)gcemalloc(sizeof(double)*AP.MTHETA);
  (*ene).DEQ=(double *)gcemalloc(sizeof(double)*AP.MPHIA);

  (*ene).p_cnb=(double *)gcemalloc(sizeof(double)*numatom);
  (*ene).p_nnb=(double *)gcemalloc(sizeof(double)*numatom);
  (*ene).p_d=(double *)gcemalloc(sizeof(double)*AP.MPHIA);
  (*ene).p_a=(double *)gcemalloc(sizeof(double)*AP.MTHETA);
  (*ene).p_b=(double *)gcemalloc(sizeof(double)*AP.MBONA);

  indexcnb_cradii=SBAAff_make_nc_list(&((*ene).numcnb),cord,numatom,(*ene).criteria);
  (*ene).indexnnb=SBAAff_make_nn_list(&((*ene).numnnb),cord,numatom,(*ene).criteria);

  (*ene).cradii=(double *)gcemalloc(sizeof(double)*(*ene).numcnb*2);
  (*ene).indexcnb=(int *)gcemalloc(sizeof(int)*(*ene).numcnb*2);

  for (i=0;i<(*ene).numcnb;++i) {
    (*ene).cradii[i]=indexcnb_cradii[i*3+2];
    (*ene).indexcnb[i*2]=(int)(indexcnb_cradii[i*3]);
    (*ene).indexcnb[i*2+1]=(int)(indexcnb_cradii[i*3+1]);
  }

  SBAAff_get_eq_val(cord,(*ene).BEQ,(*ene).AEQ,(*ene).DEQ);

  (*ene).ALJ=(double *)gcemalloc(sizeof(double)*(*ene).numcnb);
  (*ene).BLJ=(double *)gcemalloc(sizeof(double)*(*ene).numcnb);

  SBAAff_set_CNB_PARM((*ene).cradii,(*ene).numcnb,(*ene).ecnb,(*ene).ALJ,(*ene).BLJ,numatom);
  SBAAff_set_NNB_PARM((*ene).cradiinnb,(*ene).numnnb,(*ene).ennb,&((*ene).ALJnnb),numatom);

}

int SBAAff_set_protein_calcff(struct potential_SBAA *ene, double *cord,int numatom){
  int i;
  double *indexcnb_cradii;
    
  (*ene).BEQ=(double *)gcemalloc(sizeof(double)*AP.MBONA);
  (*ene).AEQ=(double *)gcemalloc(sizeof(double)*AP.MTHETA);
  (*ene).DEQ=(double *)gcemalloc(sizeof(double)*AP.MPHIA);

  (*ene).p_cnb=(double *)gcemalloc(sizeof(double)*numatom);
  (*ene).p_nnb=(double *)gcemalloc(sizeof(double)*numatom);
  (*ene).p_d=(double *)gcemalloc(sizeof(double)*AP.MPHIA);
  (*ene).p_a=(double *)gcemalloc(sizeof(double)*AP.MTHETA);
  (*ene).p_b=(double *)gcemalloc(sizeof(double)*AP.MBONA);

  indexcnb_cradii=SBAAff_make_nc_protein_list(&((*ene).numcnb),cord,numatom,(*ene).criteria);
  (*ene).indexnnb=SBAAff_make_nn_protein_list(&((*ene).numnnb),cord,numatom,(*ene).criteria);

  (*ene).cradii=(double *)gcemalloc(sizeof(double)*(*ene).numcnb*2);
  (*ene).indexcnb=(int *)gcemalloc(sizeof(int)*(*ene).numcnb*2);

  for (i=0;i<(*ene).numcnb;++i) {
    (*ene).cradii[i]=indexcnb_cradii[i*3+2];
    (*ene).indexcnb[i*2]=(int)(indexcnb_cradii[i*3]);
    (*ene).indexcnb[i*2+1]=(int)(indexcnb_cradii[i*3+1]);
  }

  SBAAff_get_eq_val(cord,(*ene).BEQ,(*ene).AEQ,(*ene).DEQ);

  (*ene).ALJ=(double *)gcemalloc(sizeof(double)*(*ene).numcnb);
  (*ene).BLJ=(double *)gcemalloc(sizeof(double)*(*ene).numcnb);

  SBAAff_set_CNB_PARM((*ene).cradii,(*ene).numcnb,(*ene).ecnb,(*ene).ALJ,(*ene).BLJ,numatom);
  SBAAff_set_NNB_PARM((*ene).cradiinnb,(*ene).numnnb,(*ene).ennb,&((*ene).ALJnnb),numatom);

}


int SBAAff_calcFFCNB(double *ALJ, double *BLJ,double *p_cnb,int numcnb, int *indexcnb,int num_atom,double *cord) {
  int i,j;
  int num_a_prot,NUM_A_PROT;
  double potedummy;
  double vec[3];
  double len,len2,len6;
  double p12,p6;

  for(i=0;i<num_atom;++i) p_cnb[i]=0.0;

  for(i=0;i<numcnb;++i){
    num_a_prot=indexcnb[i*2];
    NUM_A_PROT=indexcnb[i*2+1];
    len2 = 0.0;
    for(j=0;j<3;++j){
      vec[j] = cord[NUM_A_PROT*3+j]-cord[num_a_prot*3+j];
      len2 += vec[j]*vec[j];
    }
    len = sqrt(len2);
    len6=len2;
    for (j=0;j<2;++j)  len6 = len6*len2;
    p12 = ALJ[i]/(len6*len6);
    p6  = BLJ[i]/len6;
    p_cnb[num_a_prot] += 0.5*(p12-2.0*p6);
    p_cnb[NUM_A_PROT] += 0.5*(p12-2.0*p6);
  }

  return 0;
}

int SBAAff_calcFFNNB(double ALJ,double *p_nnb,int numnnb, int *indexnnb,int num_atom,double *cord) {
  int i,j;
  int num_a_prot,NUM_A_PROT;
  double potedummy;
  double vec[3];
  double len,len2,len6;
  double p12;

  for(i=0;i<num_atom;++i) p_nnb[i]=0.0;

  for(i=0;i<numnnb;++i){
    num_a_prot=indexnnb[i*2];
    NUM_A_PROT=indexnnb[i*2+1];
    len2 = 0.0;
    for(j=0;j<3;++j){
      vec[j] = cord[NUM_A_PROT*3+j]-cord[num_a_prot*3+j];
      len2 += vec[j]*vec[j];
    }
    len = sqrt(len2);
    len6=len2;
    for (j=0;j<2;++j)  len6 = len6*len2;
    p12 = ALJ/(len6*len6);
    p_nnb[num_a_prot] += 0.5*p12;
    p_nnb[NUM_A_PROT] += 0.5*p12;
  }

  return 0;
}


int SBAAff_calcDIHE(double *p_d,double *cord,double *DEQ,double EBB,double ESC, double EI) {
  int i,j,k;
  double atom[4][3];
  double dihedang;
  double dang,pi;

  pi=acos(-1);
  
  for (i=0;i<AP.MPHIA;++i)  p_d[i] = 0.0;

  for (i=0;i<AP.MPHIA;++i) {
    for (j=0;j<4;++j) for (k=0;k<3;++k)	atom[j][k]=cord[abs(AP.PA[i][j])+k];
  
    dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
    if (dihedang>pi) dihedang-=2.0*pi;
    else if (dihedang<-1.0*pi) dihedang+=2.0*pi;

    if ((dang=dihedang-DEQ[i])>pi) dang-=2.0*pi;
    else if ((dang==dihedang-DEQ[i])<-1.0*pi) dang+=2.0*pi;

    if (AP.PK[AP.PA[i][3]-1]<1.0) {
      p_d[i] = EBB*((1.0-cos(dang))+0.5*(1.0-cos(3.0*(dang))));
    }
    else {
      p_d[i] = EI*dang*dang;
    }
  }
}

int SBAAff_calcANGLE(double *p_a,double *cord, double EA, double *AEQ){
  int i,j,k,l;
  double atom[3][3];
  double ang,dang;
  double pi;
  
  pi=acos(-1.0);
  
  for (i=0;i<AP.MTHETA;++i) p_a[i] = 0.0;
  
  for (i=0;i<AP.MTHETA;++i) {
     for (j=0;j<3;++j)
      for (k=0;k<3;++k)
	atom[j][k]=cord[abs(AP.TA[i][j])+k];
  
     ang = pick_angle(atom[0],atom[1],atom[2],0,0.0);
    if (ang>pi) ang-=2.0*pi;
    else if (ang<-1.0*pi) ang+=2.0*pi;

    if ((dang=ang-AEQ[i])>pi) dang-=2.0*pi;
    else if ((dang=ang-AEQ[i])<-1.0*pi) dang+=2.0*pi;

     p_a[i] = EA*dang*dang;
  }

  return 0;
}

int SBAAff_calcBOND(double *p_b,double *cord, double EB, double *BEQ){
  int i,j,k;
  double len;
  double atom[2][3];

  for (i=0;i<AP.MBONA;++i) {
    for (j=0;j<2;++j)
      for (k=0;k<3;++k)
	atom[j][k]=cord[abs(AP.BA[i][j])+k];  

    len = pick_bond_leng(atom[0],atom[1]);
    p_b[i] = EB*(len-BEQ[i])*(len-BEQ[i]);
  }

  return 0;
}

double *SBAAff_make_nc_list( int *numnc, double *cord, int numatom, double criteria) {
  int i,j,k;
  double *indexcnb_cradii;
  double len,atom[2][3];

  *numnc=0;
  indexcnb_cradii=(double *)gcemalloc(sizeof(double)*3);

  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      for (k=0;k<3;++k) {
	atom[0][k]=cord[i*3+k];
	atom[1][k]=cord[j*3+k];
      }

      len = pick_bond_leng(atom[0],atom[1]);
      if (which_calc_nb(i,j)==YES) {
	if (len < criteria) {
	  indexcnb_cradii=(double *)gcerealloc(indexcnb_cradii,sizeof(double)*((*numnc)+1)*3);
	  indexcnb_cradii[(*numnc)*3]=i;
	  indexcnb_cradii[(*numnc)*3+1]=j;
	  indexcnb_cradii[(*numnc)*3+2]=len;
	  ++(*numnc);
	}	
      }
    }
  }
  return indexcnb_cradii;
}

double *SBAAff_make_nc_protein_list( int *numnc, 
				     double *cord, 
				     int numatom, 
				     double criteria) {
  int i,j,k;
  double *indexcnb_cradii;
  double len,atom[2][3];

  *numnc=0;
  indexcnb_cradii=(double *)gcemalloc(sizeof(double)*3);

  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      for (k=0;k<3;++k) {
	atom[0][k]=cord[i*3+k];
	atom[1][k]=cord[j*3+k];
      }

      len = pick_bond_leng(atom[0],atom[1]);
      if (within_3_neibor_res(i,j)==NO) {
	if (len < criteria) {
	  indexcnb_cradii=(double *)gcerealloc(indexcnb_cradii,sizeof(double)*((*numnc)+1)*3);
	  indexcnb_cradii[(*numnc)*3]=i;
	  indexcnb_cradii[(*numnc)*3+1]=j;
	  indexcnb_cradii[(*numnc)*3+2]=len;
	  ++(*numnc);
	}	
      }
    }
  }
  return indexcnb_cradii;
}


int *SBAAff_make_nn_list(int *numnn, double *cord, int numatom, double criteria) {
  int i,j,k;
  int *indexnnb;
  double len,atom[2][3];

  indexnnb=(int *)gcemalloc(sizeof(int)*2);

  *numnn=0;

  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      for (k=0;k<3;++k) {
	atom[0][k]=cord[i*3+k];
	atom[1][k]=cord[j*3+k];
      }

      len = pick_bond_leng(atom[0],atom[1]);
      if (which_calc_nb(i,j)==YES) {
	if (len >= criteria) {
	  indexnnb=(int *)gcerealloc(indexnnb,sizeof(int)*((*numnn)+1)*3);
	  indexnnb[(*numnn)*2]=i;
	  indexnnb[(*numnn)*2+1]=j;
	  ++(*numnn);
	}
      }
    }
  }
  return indexnnb;
}

int *SBAAff_make_nn_protein_list(int *numnn, double *cord, int numatom, double criteria) {
  int i,j,k;
  int flag;
  int *indexnnb;
  double len,atom[2][3];

  indexnnb=(int *)gcemalloc(sizeof(int)*2);

  *numnn=0;

  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      for (k=0;k<3;++k) {
	atom[0][k]=cord[i*3+k];
	atom[1][k]=cord[j*3+k];
      }

      len = pick_bond_leng(atom[0],atom[1]);
      if (which_calc_nb(i,j)==YES) {
	flag=ON;
	if (within_3_neibor_res(i,j)==NO) {
	  if (len < criteria) {
	    flag==OFF;
	  }
	  if (flag==ON) {
	    indexnnb=(int *)gcerealloc(indexnnb,sizeof(int)*((*numnn)+1)*3);
	    indexnnb[(*numnn)*2]=i;
	    indexnnb[(*numnn)*2+1]=j;
	    ++(*numnn);
	  }
	}
      }
    }
  }
  return indexnnb;
}


void SBAAff_set_CNB_PARM(double *cradii, int numcnb,double ecnb, double *ALJ, double *BLJ, int numatom){
  int i,j,k;
  double len6,len12;

  for (i=0;i<numcnb;++i) {
    len6=cradii[i]*cradii[i]*cradii[i]*cradii[i]*cradii[i]*cradii[i];
    len12=len6*len6;
    ALJ[i] = ecnb*len12;
    BLJ[i] = ecnb*len6;
  }
}
 
void SBAAff_set_NNB_PARM(double cradiinnb,int numnnb,double ennb,double *ALJ,int numatom){
  int i,j,k;
  double len12;
  
  for (i=0;i<numnnb;++i) {
    len12=cradiinnb*cradiinnb*cradiinnb
      *cradiinnb*cradiinnb*cradiinnb
      *cradiinnb*cradiinnb*cradiinnb
      *cradiinnb*cradiinnb*cradiinnb;
    
    *ALJ = ennb*len12;
  }
}


int SBAAff_get_eq_val(double *cord, double *BEQ, double *AEQ, double *DEQ) {
  int i,j,k;
  int numatom;
  double len,atom[4][3];
  double pi;

  pi=acos(-1.0);

  for (i=0;i<AP.MBONA;++i) {
    for (j=0;j<2;++j)
      for (k=0;k<3;++k)
	atom[j][k]=cord[abs(AP.BA[i][j])+k];  

    BEQ[i] = pick_bond_leng(atom[0],atom[1]);
  }

  for (i=0;i<AP.MTHETA;++i) {
     for (j=0;j<3;++j)
      for (k=0;k<3;++k)
	atom[j][k]=cord[abs(AP.TA[i][j])+k];
  
     AEQ[i] = pick_angle(atom[0],atom[1],atom[2],0,0.0);
     if (AEQ[i]>pi) AEQ[i]-=2.0*pi;
     else if (AEQ[i]<-1.0*pi)  AEQ[i]+=2.0*pi;
  }

  for (i=0;i<AP.MPHIA;++i) {
    for (j=0;j<4;++j) for (k=0;k<3;++k)	atom[j][k]=cord[abs(AP.PA[i][j])+k];
  
    DEQ[i] = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
     if (DEQ[i]>pi) DEQ[i]-=2.0*pi;
     else if (DEQ[i]<-1.0*pi) DEQ[i]+=2.0*pi;
  }

}

void SBAAff_set_parameters_default(struct potential_SBAA *ene) {
  (*ene).cradiinnb=2.5;
  (*ene).criteria=3.0;
  (*ene).EB=100; // Tue May 31 18:24:56 2011
  (*ene).EA=20;  // Tue May 31 18:24:56 2011
  //  (*ene).EB=500;
  //  (*ene).EA=100;
  (*ene).EI=10;

  (*ene).EBB=1.0;
  (*ene).ESC=1.0;

  (*ene).ennb=0.01;
  (*ene).ecnb=1.0;

}

void SBAAff_set_parameters_protein_default(struct potential_SBAA *ene) {
  (*ene).cradiinnb=2.5;
  (*ene).criteria=6.5;
  (*ene).EB=100;
  (*ene).EA=20;
  (*ene).EI=10;

  (*ene).EBB=0.5;
  (*ene).ESC=0.5;

  (*ene).ennb=0.18;
  (*ene).ecnb=0.18;

}


int SBAAff_check_parameters(struct potential_SBAA *ene, int numatom) {
  int i;
  FILE *log;

  log=efopen("log_indexcnb.txt","w");
  for(i=0;i<(*ene).numcnb;++i) fprintf(log,"%d - %d %e\n",(*ene).indexcnb[i*2]+1,(*ene).indexcnb[i*2+1]+1,(*ene).cradii[i]);
  fprintf(log,"\n  ");
  fclose(log);

  log=efopen("log_indexnnb.txt","w");
  for(i=0;i<(*ene).numnnb;++i) fprintf(log,"%d - %d\n",(*ene).indexnnb[i*2]+1,(*ene).indexnnb[i*2+1]+1);
  fprintf(log,"\n  ");
  fclose(log);

  log=efopen("log_BEQ.txt","w");
  for(i=0;i<AP.MBONA;++i) 
    fprintf(log,"%d - %d %e\n",(abs(AP.BA[i][0]))/3,(abs(AP.BA[i][1]))/3,(*ene).BEQ[i]);
  fprintf(log,"\n  ");
  fclose(log);

  log=efopen("log_AEQ.txt","w");
  for(i=0;i<AP.MTHETA;++i) 
    fprintf(log,"%d - %d - %d %e\n",(abs(AP.TA[i][0]))/3,(abs(AP.TA[i][1]))/3,(abs(AP.TA[i][2]))/3,(*ene).AEQ[i]);
  fprintf(log,"\n  ");
  fclose(log);

  log=efopen("log_DEQ.txt","w");
  for(i=0;i<AP.MPHIA;++i) 
    fprintf(log,"%d - %d - %d -%d %e\n",(abs(AP.PA[i][0]))/3,(abs(AP.PA[i][1]))/3,(abs(AP.PA[i][2]))/3,(abs(AP.TA[i][2]))/3,(*ene).DEQ[i]);
  fprintf(log,"\n  ");
  fclose(log);

  return 0;
}

int SBAAMBff_set_calcff(struct potential_SBAA *ene1,struct potential_SBAA *ene2,double *deltaV, double *cord1,double *cord2,int numatom){

  SBAAff_set_calcff(ene1,cord1,numatom);
  SBAAff_set_calcff(ene2,cord2,numatom);

  SBAAff_calcff(cord1,numatom,ene1);
  SBAAff_calcff(cord2,numatom,ene2);

  *deltaV=ene1->p_t-ene2->p_t;
}

int SBAAMBff_set_protein_calcff(struct potential_SBAA *ene1,struct potential_SBAA *ene2,double *deltaV, double *cord1,double *cord2,int numatom){

  SBAAff_set_protein_calcff(ene1,cord1,numatom);
  SBAAff_set_protein_calcff(ene2,cord2,numatom);

  SBAAff_calcff(cord1,numatom,ene1);
  SBAAff_calcff(cord2,numatom,ene2);

  *deltaV=ene1->p_t-ene2->p_t;
}


double SBAAMBff_calcff(double *crd,int numatom,struct potential_SBAA *ene1,struct potential_SBAA *ene2, double delta, double deltaV){
  double p_t;
  double p_1,p_2;
  double p_12d;
  double p_sub1,p_sub2;

  SBAAff_calcff(crd,numatom,ene1);
  SBAAff_calcff(crd,numatom,ene2);
  
  p_1=ene1->p_t;
  p_2=ene2->p_t;

  p_sub1=0.5*(p_1+p_2+deltaV);
  p_sub2=0.5*(p_1-p_2-deltaV);

  p_t=p_sub1-sqrt(p_sub2*p_sub2+delta*delta);

  return p_t;
}

int within_3_neibor_res(int num1, int num2) {
  int numres1,numres2;

  numres1=resnum(num1);
  numres2=resnum(num2);

  if (numres2 < numres1+3)
    return YES;
  else 
    return NO;
}

int resnum(int numatom) {
  int i;

  for (i=0;i<AP.NRES;++i) {
    if (AP.IPRES[i]-1 > numatom)
      return i;
  }
  return -1;
}
