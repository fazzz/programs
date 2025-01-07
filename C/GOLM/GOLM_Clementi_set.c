
#include <stdio.h>
#include <math.h>

#include "GOLM_Clementi_set.h"

#include "PTL.h"
#include "TOPO.h"
#include "NC.h"

#include "MB.h"
#include "EF.h"

int GOLM_Clementi_ff_set_calcff(struct potential_GOLM_Clementi *ene, double *refcrd,double *refcrdAA, int numatom,int numatomAA){
  int i,j,k,nc;
  int ii,jj,temp;
  int numca;

  double atom1[3],atom2[3],atom3[3],atom4[3];
  double criteria=criteria_Clementi;
  double length,len10,len12;
  double *indexcnb_cradii;
  int **ncmap,*indexncb,numnc;
  int **bp_f,*numb,**nb_matrix;

  double pi;

  pi=acos(-1.0);

  /////////////////////////////////////////////////////////////////////////////////////////////////////

  (*ene).bon_equ=(double *)gcemalloc(sizeof(double)*(numatom-1));
  (*ene).Kb=100.0*epsilon_Clementi;
  for (i=0;i<numatom-1;++i) {
    for (j=0;j<3;++j) {
      atom1[j]=refcrd[i*3+j];
      atom2[j]=refcrd[(i+1)*3+j];
    }
    (*ene).bon_equ[i] = len(atom1,atom2);
  }

  (*ene).ang_equ=(double *)gcemalloc(sizeof(double)*(numatom-2));
  (*ene).Ka=20.0*epsilon_Clementi;
  for (i=0;i<numatom-2;++i) {
    for (j=0;j<3;++j) {
      atom1[j]=refcrd[i*3+j];
      atom2[j]=refcrd[(i+1)*3+j];
      atom3[j]=refcrd[(i+2)*3+j];
    }
    (*ene).ang_equ[i] = /*pick_*/ang/*le*/(atom1,atom2,atom3/*,0,0.0*/);
  }

  (*ene).dih_equ=(double *)gcemalloc(sizeof(double)*(numatom-3));
  (*ene).Kd1=1.0*epsilon_Clementi;
  (*ene).Kd2=0.5*epsilon_Clementi;
  for (i=0;i<numatom-3;++i) {
    for (j=0;j<3;++j) {
      atom1[j]=refcrd[i*3+j];
      atom2[j]=refcrd[(i+1)*3+j];
      atom3[j]=refcrd[(i+2)*3+j];
      atom4[j]=refcrd[(i+3)*3+j];
    }
    (*ene).dih_equ[i] = pick_dihed(atom1,atom2,atom3,atom4,0,0.0);
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////

  ncmap=GOLM_Clementi_make_native_contact(refcrdAA,criteria,&((*ene).num_natatt),numatomAA,numatom);

  nc=0;
  (*ene).index_natatt=(int *)gcemalloc(sizeof(int)*(*ene).num_natatt*2);
  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      if (ncmap[i][j]==0) {
	(*ene).index_natatt[nc*2+0]=i;
	(*ene).index_natatt[nc*2+1]=j;
	++nc;
      }
    }
  }

  (*ene).cradii_natatt=(double *)gcemalloc(sizeof(double)*(*ene).num_natatt);
  (*ene).ALJ_natatt=(double *)gcemalloc(sizeof(double)*(*ene).num_natatt);
  (*ene).BLJ_natatt=(double *)gcemalloc(sizeof(double)*(*ene).num_natatt);
  for (i=0;i<(*ene).num_natatt;++i) {
    for (j=0;j<3;++j) {
      atom1[j]=refcrd[((*ene).index_natatt[i*2+0])*3+j];
      atom2[j]=refcrd[((*ene).index_natatt[i*2+1])*3+j];
    }
    length=len(atom1,atom2);
    (*ene).cradii_natatt[i]=len(atom1,atom2);
    len10=(*ene).cradii_natatt[i];
    for (j=0;j<9;++j)  len10 = len10*(*ene).cradii_natatt[i];
    len12=(*ene).cradii_natatt[i];
    for (j=0;j<11;++j)  len12 = len12*(*ene).cradii_natatt[i];
    (*ene).ALJ_natatt[i]=len12;
    (*ene).BLJ_natatt[i]=len10;
  }
  (*ene).ep_natatt=epsilon_Clementi;

  (*ene).index_repul=(int *)gcemalloc(sizeof(int)*2);
  k=0;
  for (i=0;i<numatom;++i) {
    for (j=i+4;j<numatom;++j) {
      if ( ncmap[i][j]!=0 ) {
	++k;
	(*ene).index_repul=(int *)gcerealloc((*ene).index_repul,sizeof(int)*2*k);
	(*ene).index_repul[(k-1)*2+0]=i;
	(*ene).index_repul[(k-1)*2+1]=j;
      }
    }
  }
  (*ene).num_repul=k;

  (*ene).ep_repul=epsilon_Clementi;
  (*ene).cradii_repul=cradii_repul_Clementi;
  len12=(*ene).cradii_repul;
  for (i=0;i<11;++i) len12=len12*(*ene).cradii_repul;
  (*ene).ALJ_repul=(*ene).ep_repul*len12;

  /////////////////////////////////////////////////////////////////////////////////////////////////////

  (*ene).p_natatt=(double *)gcemalloc(sizeof(double)*(*ene).num_natatt);
  (*ene).p_repul=(double *)gcemalloc(sizeof(double)*(*ene).num_repul);
  (*ene).p_d=(double *)gcemalloc(sizeof(double)*(numatom-3));
  (*ene).p_a=(double *)gcemalloc(sizeof(double)*(numatom-2));
  (*ene).p_b=(double *)gcemalloc(sizeof(double)*(numatom-1));

  (*ene).f_natatt=(double **)gcemalloc(sizeof(double *)*numatom);
  (*ene).f_repul=(double **)gcemalloc(sizeof(double *)*numatom);
  (*ene).f_d=(double **)gcemalloc(sizeof(double *)*numatom);
  (*ene).f_a=(double **)gcemalloc(sizeof(double *)*numatom);
  (*ene).f_b=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) {
    (*ene).f_natatt[i]=(double *)gcemalloc(sizeof(double)*3);
    (*ene).f_repul[i]=(double *)gcemalloc(sizeof(double)*3);
    (*ene).f_d[i]=(double *)gcemalloc(sizeof(double)*3);
    (*ene).f_a[i]=(double *)gcemalloc(sizeof(double)*3);
    (*ene).f_b[i]=(double *)gcemalloc(sizeof(double)*3);
  }

  (*ene).f_t=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) (*ene).f_t[i]=(double *)gcemalloc(sizeof(double)*3);

  /////////////////////////////////////////////////////////////////////////////////////////////////////
}

int GOLM_Clementi_ff_set_calcff2(struct potential_GOLM_Clementi *ene, double *refcrd,double *refcrdAA, int numatom,int numatomAA, double ep){
  int i,j,k,nc;
  int ii,jj,temp;
  int numca;

  double atom1[3],atom2[3],atom3[3],atom4[3];
  double criteria=criteria_Clementi;
  double length,len10,len12;
  double *indexcnb_cradii;
  int **ncmap,*indexncb,numnc;
  int **bp_f,*numb,**nb_matrix;

  double pi;

  pi=acos(-1.0);

  /////////////////////////////////////////////////////////////////////////////////////////////////////

  (*ene).bon_equ=(double *)gcemalloc(sizeof(double)*(numatom-1));
  (*ene).Kb=100.0*epsilon_Clementi;
  for (i=0;i<numatom-1;++i) {
    for (j=0;j<3;++j) {
      atom1[j]=refcrd[i*3+j];
      atom2[j]=refcrd[(i+1)*3+j];
    }
    (*ene).bon_equ[i] = len(atom1,atom2);
  }

  (*ene).ang_equ=(double *)gcemalloc(sizeof(double)*(numatom-2));
  (*ene).Ka=20.0*epsilon_Clementi;
  for (i=0;i<numatom-2;++i) {
    for (j=0;j<3;++j) {
      atom1[j]=refcrd[i*3+j];
      atom2[j]=refcrd[(i+1)*3+j];
      atom3[j]=refcrd[(i+2)*3+j];
    }
    (*ene).ang_equ[i] = /*pick_*/ang/*le*/(atom1,atom2,atom3/*,0,0.0*/);
  }

  (*ene).dih_equ=(double *)gcemalloc(sizeof(double)*(numatom-3));
  (*ene).Kd1=1.0*epsilon_Clementi;
  (*ene).Kd2=0.5*epsilon_Clementi;
  for (i=0;i<numatom-3;++i) {
    for (j=0;j<3;++j) {
      atom1[j]=refcrd[i*3+j];
      atom2[j]=refcrd[(i+1)*3+j];
      atom3[j]=refcrd[(i+2)*3+j];
      atom4[j]=refcrd[(i+3)*3+j];
    }
    (*ene).dih_equ[i] = pick_dihed(atom1,atom2,atom3,atom4,0,0.0);
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////

  ncmap=GOLM_Clementi_make_native_contact(refcrdAA,criteria,&((*ene).num_natatt),numatomAA,numatom);

  nc=0;
  (*ene).index_natatt=(int *)gcemalloc(sizeof(int)*(*ene).num_natatt*2);
  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      if (ncmap[i][j]==0) {
	(*ene).index_natatt[nc*2+0]=i;
	(*ene).index_natatt[nc*2+1]=j;
	++nc;
      }
    }
  }

  (*ene).cradii_natatt=(double *)gcemalloc(sizeof(double)*(*ene).num_natatt);
  (*ene).ALJ_natatt=(double *)gcemalloc(sizeof(double)*(*ene).num_natatt);
  (*ene).BLJ_natatt=(double *)gcemalloc(sizeof(double)*(*ene).num_natatt);
  for (i=0;i<(*ene).num_natatt;++i) {
    for (j=0;j<3;++j) {
      atom1[j]=refcrd[((*ene).index_natatt[i*2+0])*3+j];
      atom2[j]=refcrd[((*ene).index_natatt[i*2+1])*3+j];
    }
    length=len(atom1,atom2);
    (*ene).cradii_natatt[i]=len(atom1,atom2);
    len10=(*ene).cradii_natatt[i];
    for (j=0;j<9;++j)  len10 = len10*(*ene).cradii_natatt[i];
    len12=(*ene).cradii_natatt[i];
    for (j=0;j<11;++j)  len12 = len12*(*ene).cradii_natatt[i];
    (*ene).ALJ_natatt[i]=len12;
    (*ene).BLJ_natatt[i]=len10;
  }
  (*ene).ep_natatt=ep;

  (*ene).index_repul=(int *)gcemalloc(sizeof(int)*2);
  k=0;
  for (i=0;i<numatom;++i) {
    for (j=i+4;j<numatom;++j) {
      if ( ncmap[i][j]!=0 ) {
	++k;
	(*ene).index_repul=(int *)gcerealloc((*ene).index_repul,sizeof(int)*2*k);
	(*ene).index_repul[(k-1)*2+0]=i;
	(*ene).index_repul[(k-1)*2+1]=j;
      }
    }
  }
  (*ene).num_repul=k;

  (*ene).ep_repul=epsilon_Clementi;
  (*ene).cradii_repul=cradii_repul_Clementi;
  len12=(*ene).cradii_repul;
  for (i=0;i<11;++i) len12=len12*(*ene).cradii_repul;
  (*ene).ALJ_repul=(*ene).ep_repul*len12;

  /////////////////////////////////////////////////////////////////////////////////////////////////////

  (*ene).p_natatt=(double *)gcemalloc(sizeof(double)*(*ene).num_natatt);
  (*ene).p_repul=(double *)gcemalloc(sizeof(double)*(*ene).num_repul);
  (*ene).p_d=(double *)gcemalloc(sizeof(double)*(numatom-3));
  (*ene).p_a=(double *)gcemalloc(sizeof(double)*(numatom-2));
  (*ene).p_b=(double *)gcemalloc(sizeof(double)*(numatom-1));

  (*ene).f_natatt=(double **)gcemalloc(sizeof(double *)*numatom);
  (*ene).f_repul=(double **)gcemalloc(sizeof(double *)*numatom);
  (*ene).f_d=(double **)gcemalloc(sizeof(double *)*numatom);
  (*ene).f_a=(double **)gcemalloc(sizeof(double *)*numatom);
  (*ene).f_b=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) {
    (*ene).f_natatt[i]=(double *)gcemalloc(sizeof(double)*3);
    (*ene).f_repul[i]=(double *)gcemalloc(sizeof(double)*3);
    (*ene).f_d[i]=(double *)gcemalloc(sizeof(double)*3);
    (*ene).f_a[i]=(double *)gcemalloc(sizeof(double)*3);
    (*ene).f_b[i]=(double *)gcemalloc(sizeof(double)*3);
  }

  (*ene).f_t=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) (*ene).f_t[i]=(double *)gcemalloc(sizeof(double)*3);

  /////////////////////////////////////////////////////////////////////////////////////////////////////
}

int **GOLM_Clementi_make_native_contact(double *refcrdAA,double criteria,int *numnc,int numatom, int numCAatom) {
  int i,j,k,l,d,na,ii,jj;
  int numres,numstat,numatomstop;
  int *AAtoCAnum;
  int **ncmap;

  double atom1[3],atom2[3];
  double length,pi;

  pi=acos(-1.0);

  AAtoCAnum=(int **)gcemalloc(sizeof(int)*numatom);

  numres=AP.NRES;
  if (strncmp(AP.LABERES[0],"ACE",3)!=0) {
    numstat=0;
    for (i=0;i<numatom;++i) AAtoCAnum[i]=PTL_resnum(i);
  }
  else {
    numstat=AP.IPRES[1]-1;
    for (i=0;i<numatom;++i) AAtoCAnum[i]=PTL_resnum(i)-1;
  }
  if ( strncmp(AP.LABERES[numres-1],"NME",3)!=0) {
    numatomstop=numatom;
  }
  else {
    numatomstop=AP.IPRES[numres-1]-1;
  }

  ncmap=(int **)gcemalloc(sizeof(int *)*numCAatom);
  for (i=0;i<numCAatom;++i) ncmap[i]=(int *)gcemalloc(sizeof(int)*numCAatom);

  for (i=0;i<numCAatom;++i) 
    for (j=0;j<numCAatom;++j) 
      ncmap[i][j]=-100;

  *numnc=0;
  for (i=numstat;i<numatomstop;++i) {
    for (j=i+1;j<numatomstop;++j) {
      ii=AAtoCAnum[i]-1;
      jj=AAtoCAnum[j]-1;
      if ( ii < jj-3  && strncmp(AP.IGRAPH[i],"H",1)!=0 && strncmp(AP.IGRAPH[j],"H",1)!=0 ) {
	for (k=0;k<3;++k) {
	  atom1[k]=refcrdAA[i*3+k];
	  atom2[k]=refcrdAA[j*3+k];
	}
	length=len(atom1,atom2);
	if (length <= criteria && (ncmap[ii][jj]!=0)) {
	  ++(*numnc);
	  ncmap[ii][jj]=0;
	}
      }
    }
  }

  return ncmap;
}
