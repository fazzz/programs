
#include <stdio.h>
#include <math.h>

#include "PT.h"
#include "GOLMAA_set.h"
#include "TOPO.h"
#include "NC_check.h"
#include "MB.h"
#include "EF.h"

int GOLMAAff_set_calcff(struct potential_GOLMAA *ene, double *refcrd,int numatom,int **nb_matrix,double R_C_D,double constant){
  int i,j;
  int ii,jj,temp;
  int numca;
  int delnum;
  double atom1[3],atom2[3],atom3[3],atom4[3];
  double criteria;
  double length,len6,len12;
  double *indexcnb_cradii;
  int numres,**ncmap,*indexncb,numnc;
  int **bp_f,*numb;

  //  double constant=1.0;
  double pi;

  pi=acos(-1.0);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  //  (*ene).FC_dihed=1.0/AP.MPHIA;
  (*ene).FC_dihed=constant/AP.MPHIA*(1.0/(R_C_D+1.0));
  (*ene).DEQ=(double *)gcemalloc(sizeof(double)*AP.MPHIA);
  (*ene).p_d=(double *)gcemalloc(sizeof(double)*AP.MPHIA);
  // (*ene).f_d=(double **)gcemalloc(sizeof(double *)*numatom);
  //  for (i=0;i<numatom;++i) (*ene).f_d[i]=(double *)gcemalloc(sizeof(double)*3);

  for (i=0;i<AP.MPHIA;++i) {
    for (j=0;j<3;++j) {
      atom1[j]=refcrd[abs(AP.PA[i][0])+j];
      atom2[j]=refcrd[abs(AP.PA[i][1])+j];
      atom3[j]=refcrd[abs(AP.PA[i][2])+j];
      atom4[j]=refcrd[abs(AP.PA[i][3])+j];
    }

    (*ene).DEQ[i] = pick_dihed(atom1,atom2,atom3,atom4,0,0.0);

    if ((*ene).DEQ[i]>pi) (*ene).DEQ[i]-=2.0*pi;
    else if ((*ene).DEQ[i]<-1.0*pi) (*ene).DEQ[i]+=2.0*pi;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  numres=AP.NRES;
  (*ene).ncmap=(int **)gcemalloc(sizeof(int *)*numatom);
  for (i=0;i<numatom;++i) (*ene).ncmap[i]=(int *)gcemalloc(sizeof(int)*numatom);
  (*ene).index_natatt=make_native_contact_list_aa(&((*ene).num_natatt),refcrd,numatom,criteria_NC,(*ene).ncmap,ON);

  (*ene).ALJ_natatt=(double *)gcemalloc(sizeof(double)*(*ene).num_natatt);
  (*ene).BLJ_natatt=(double *)gcemalloc(sizeof(double)*(*ene).num_natatt);
  (*ene).cradii_natatt=(double *)gcemalloc(sizeof(double)*(*ene).num_natatt);
  //  (*ene).ep_natatt=(*ene).FC_dihed*R_C_D;
  (*ene).ep_natatt=constant/(*ene).num_natatt*(R_C_D/(R_C_D+1.0));


  for (i=0;i<(*ene).num_natatt;++i) {
    for (j=0;j<3;++j) {
      atom1[j]=refcrd[((*ene).index_natatt[i*2])*3+j];
      atom2[j]=refcrd[((*ene).index_natatt[i*2+1])*3+j];
    }
    length=len(atom1,atom2);
    (*ene).cradii_natatt[i]=len(atom1,atom2);
    len6=(*ene).cradii_natatt[i];
    for (j=0;j<5;++j)  len6 = len6*(*ene).cradii_natatt[i];
    len12=(*ene).cradii_natatt[i];
    for (j=0;j<11;++j)  len12 = len12*(*ene).cradii_natatt[i];
    (*ene).ALJ_natatt[i]=len12;
    (*ene).BLJ_natatt[i]=len6;
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  bp_f=(int **)gcemalloc(sizeof(int *)*numatom);
  for (i=0;i<numatom;++i) bp_f[i]=(int *)gcemalloc(sizeof(int)*1);
  numb=(int *)gcemalloc(sizeof(int)*numatom);
  for (i=0;i<numatom;++i) {
    numb[i]=0;
  }
  for (i=0;i<AP.MBONA;++i) {
    ii=AP.BA[i][0]/3;
    jj=AP.BA[i][1]/3;
    if ( ii > jj ) {
      temp=ii;
      ii=jj;
      jj=temp;
    }
    bp_f[ii]=(int *)gcerealloc(bp_f[ii],sizeof(int)*(numb[ii]+1));
    bp_f[ii][numb[ii]]=jj;
    numb[ii]+=1;
  }
  for (i=0;i<numatom;++i) nb_matrix[i]=(int *)gcemalloc(sizeof(int)*numatom);
  make_nb_matrix(bp_f,numb,4,nb_matrix,numatom);

  (*ene).ep_repul=ep_repul_Sanbonmatsu;
  (*ene).cradii_repul=cradii_repul_Sanbonmatsu;
  len12=(*ene).cradii_repul;
  for (i=0;i<11;++i) len12=len12*(*ene).cradii_repul;
  (*ene).ALJ_repul=(*ene).ep_repul*len12;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  (*ene).p_natatt=(double *)gcemalloc(sizeof(double)*numatom);
  (*ene).p_repul=(double *)gcemalloc(sizeof(double)*numatom);

  (*ene).f_natatt=(double **)gcemalloc(sizeof(double *)*numatom);
  (*ene).f_repul=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) {
    (*ene).f_natatt[i]=(double *)gcemalloc(sizeof(double)*3);
    (*ene).f_repul[i]=(double *)gcemalloc(sizeof(double)*3);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////


  (*ene).f_t=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) {
    (*ene).f_t[i]=(double *)gcemalloc(sizeof(double)*3);
  }
}



