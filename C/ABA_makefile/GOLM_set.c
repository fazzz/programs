
#include <stdio.h>

#include "PTL.h"
#include "GOLM.h"
#include "TOPO.h"
#include "GOLM_set.h"
#include "NC.h"
#include "MB.h"
#include "EF.h"

int GOLMff_set_calcff(struct potential_GOLM *ene, double *refcrd,int numatom){
  int i,j;
  int numca;
  int delnum;
  double atom1[3],atom2[3],atom3[3],atom4[3];
  double criteria;
  double length,len10,len12;
  double *indexcnb_cradii;
  int numres,**ncmap,*index_natatt_res,numnc;

  double pi;

  numca=AP.NRES;

  pi=acos(-1.0);

  (*ene).FC_bond=FC_bond_Clementi;
  (*ene).FC_angle=FC_angl_Clementi;
  (*ene).FC_dihed1=FC_dihe1_Clementi;
  (*ene).FC_dihed2=FC_dihe2_Clementi;

  //  cradii_natatt=cradii_natatt_Clemnti;
  //  cradii_repul=cradii_repul_Clemnti;

  (*ene).index_angl=(int *)gcemalloc(sizeof(int)*2);
  (*ene).index_dihe=(int *)gcemalloc(sizeof(int)*2);

  (*ene).index_bond=set_BOND_param(&((*ene).num_bond),refcrd);
  (*ene).index_angl=set_ANGL_param(&((*ene).num_angl),refcrd);
  (*ene).index_dihe=set_DIHE_param(&((*ene).num_dihe),refcrd);

  (*ene).BEQ=(double *)gcemalloc(sizeof(double)*(*ene).num_bond);
  (*ene).AEQ=(double *)gcemalloc(sizeof(double)*(*ene).num_angl);
  (*ene).DEQ=(double *)gcemalloc(sizeof(double)*(*ene).num_dihe);

  for (i=0;i<(*ene).num_bond;++i) {
    for (j=0;j<3;++j) {
      atom1[j]=refcrd[(*ene).index_bond[i*2]*3+j];
      atom2[j]=refcrd[(*ene).index_bond[i*2+1]*3+j];
    }
    (*ene).BEQ[i]=len(atom1,atom2);
  }
  for (i=0;i<(*ene).num_angl;++i) {
    for (j=0;j<3;++j) {
      atom1[j]=refcrd[(*ene).index_angl[i*3]*3+j];
      atom2[j]=refcrd[(*ene).index_angl[i*3+1]*3+j];
      atom3[j]=refcrd[(*ene).index_angl[i*3+2]*3+j];
    }
    (*ene).AEQ[i] = ang(atom1,atom2,atom3);
    if ((*ene).AEQ[i]>pi) (*ene).AEQ[i]-=2.0*pi;
    else if ((*ene).AEQ[i]<-1.0*pi) (*ene).AEQ[i]+=2.0*pi;
  }
  for (i=0;i<(*ene).num_dihe;++i) {
    for (j=0;j<3;++j) {
      atom1[j]=refcrd[(*ene).index_dihe[i*4]*3+j];
      atom2[j]=refcrd[(*ene).index_dihe[i*4+1]*3+j];
      atom3[j]=refcrd[(*ene).index_dihe[i*4+2]*3+j];
      atom4[j]=refcrd[(*ene).index_dihe[i*4+3]*3+j];
    }
  
    (*ene).DEQ[i] = pick_dihed(atom1,atom2,atom3,atom4,0,0.0);
    if ((*ene).DEQ[i]>pi) (*ene).DEQ[i]-=2.0*pi;
    else if ((*ene).DEQ[i]<-1.0*pi) (*ene).DEQ[i]+=2.0*pi;
  }

  (*ene).p_b=(double *)gcemalloc(sizeof(double)*(*ene).num_bond);
  (*ene).p_a=(double *)gcemalloc(sizeof(double)*(*ene).num_angl);
  (*ene).p_d=(double *)gcemalloc(sizeof(double)*(*ene).num_dihe);

  (*ene).f_b=(double **)gcemalloc(sizeof(double *)*numatom);
  (*ene).f_a=(double **)gcemalloc(sizeof(double *)*numatom);
  (*ene).f_d=(double **)gcemalloc(sizeof(double *)*numatom);

  for (i=0;i<numatom;++i) {
    (*ene).f_b[i]=(double *)gcemalloc(sizeof(double)*3);
    (*ene).f_a[i]=(double *)gcemalloc(sizeof(double)*3);
    (*ene).f_d[i]=(double *)gcemalloc(sizeof(double)*3);
  }

  numres=AP.NRES;
  ncmap=(int **)gcemalloc(sizeof(int *)*numres);
  for (i=0;i<numres;++i) ncmap[i]=(int *)gcemalloc(sizeof(int)*numres);
  index_natatt_res=make_native_contact_list(&((*ene).num_natatt),refcrd,numatom,numres,criteria_NC,ncmap);

  (*ene).index_natatt=(int *)gcemalloc(sizeof(int)*(*ene).num_natatt*2);
  (*ene).cradii_natatt=(double *)gcemalloc(sizeof(double)*(*ene).num_natatt);
  j=0;
  for (i=0;i<(*ene).num_natatt;++i) {
    if (PTL_res_ca(index_natatt_res[i*2]-1)!=-1 && PTL_res_ca(index_natatt_res[i*2]-1!=-1)) {
      ++delnum;
    }
    else {
      (*ene).index_natatt[j*2]=PTL_res_ca(index_natatt_res[i*2]-1);
      (*ene).index_natatt[j*2+1]=PTL_res_ca(index_natatt_res[i*2+1]-1);
      ++j;
    }
  }
  (*ene).num_natatt-=delnum;

  (*ene).ALJ_natatt=(double *)gcemalloc(sizeof(double)*(*ene).num_natatt);
  (*ene).BLJ_natatt=(double *)gcemalloc(sizeof(double)*(*ene).num_natatt);
  (*ene).ep_natatt=ep_natatt_Clementi;
  for (i=0;i<(*ene).num_natatt;++i) {
    for (j=0;j<3;++j) {
      atom1[j]=refcrd[(*ene).index_natatt[i*2]*3+j];
      atom2[j]=refcrd[(*ene).index_natatt[i*2+1]*3+j];
    }
    length=len(atom1,atom2);
    (*ene).cradii_natatt[i]=len(atom1,atom2);
    len10=(*ene).cradii_natatt[i];
    for (j=0;j<9;++j)  len10 = len10*(*ene).cradii_natatt[i];
    len12=(*ene).cradii_natatt[i];
    for (j=0;j<11;++j)  len12 = len12*(*ene).cradii_natatt[i];
    (*ene).ALJ_natatt[i]=/*(*ene).ep_natatt**/len12;
    (*ene).BLJ_natatt[i]=/*(*ene).ep_natatt**/len10;
  }

  (*ene).ep_repul=ep_repul_Clementi;
  (*ene).cradii_repul=cradii_repul_Clementi;
  len12=(*ene).cradii_repul;
  for (i=0;i<11;++i) len12=len12*(*ene).cradii_repul;
  (*ene).ALJ_repul=(*ene).ep_repul*len12;
  (*ene).cradii_repul=cradii_repul_Clementi;

  (*ene).p_natatt=(double *)gcemalloc(sizeof(double)*numatom);
  (*ene).p_repul=(double *)gcemalloc(sizeof(double)*numatom);

  (*ene).f_natatt=(double **)gcemalloc(sizeof(double *)*numatom);
  (*ene).f_repul=(double **)gcemalloc(sizeof(double *)*numatom);

  for (i=0;i<numatom;++i) {
    (*ene).f_natatt[i]=(double *)gcemalloc(sizeof(double)*3);
    (*ene).f_repul[i]=(double *)gcemalloc(sizeof(double)*3);
  }

  (*ene).f_t=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) {
    (*ene).f_t[i]=(double *)gcemalloc(sizeof(double)*3);
  }

}

int *set_BOND_param(int *num_bond,double *refcrd) {
  int i;
  int flag1[2],flag2[2];
  int num_bond1,num_bond2;
  int numatom;

  int *index_bond;

  index_bond=(int *)gcemalloc(sizeof(int)*2);

  numatom=AP.NATOM;

  flag1[0]=-1;flag1[1]=-1;
  flag2[0]=-1;flag2[1]=-1;

  *num_bond=-1;
  for (i=0;i<numatom;++i) {
    if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
      if (flag1[0]==-1 && flag1[1]==-1 && flag2[0]==-1 && flag2[1]==-1 ) {
	flag1[0]=i;
      }
      else if (flag1[0]!=1 && flag1[1]==-1 && flag2[0]==-1 && flag2[1]==-1 ) {
	flag1[1]=i;
	flag2[0]=i;
	++(*num_bond);
	num_bond1=(*num_bond);
	index_bond=(int *)gcerealloc(index_bond,sizeof(int)*(*num_bond)*2);
      }
      else if (flag1[0]==-1 && flag1[1]==-1 && flag2[0]!=-1 && flag2[1]==-1 ) {
	flag1[0]=i;
	flag2[1]=i;
	++(*num_bond);
	num_bond2=(*num_bond);
	index_bond=(int *)gcerealloc(index_bond,sizeof(int)*(*num_bond)*2);
      }

      if (flag1[1]!=-1) {
	index_bond[num_bond1*2]=flag1[0];
	index_bond[num_bond1*2+1]=flag1[1];
	flag1[0]=-1;flag1[1]=-1;
      }

      if (flag2[1]!=-1) {
	index_bond[num_bond2*2]=flag2[0];
	index_bond[num_bond2*2+1]=flag2[1];
	flag2[0]=-1;flag2[1]=-1;
      }

    }
  }
  return index_bond;
}

int *set_ANGL_param(int *num_angl,double *refcrd) {
  int i,j,k;
  int iniflag=OFF;
  int flag[3][3];
  int num_anglj[3];
  int numatom;

  int *index_angl;

  index_angl=(int *)gcemalloc(sizeof(int)*3);

  numatom=AP.NATOM;

  for (j=0;j<3;++j) 
    for (k=0;k<3;++k) 
      flag[j][k]=-1;

  *num_angl=0;
  for (i=0;i<numatom;++i) {
    if (strncmp(AP.IGRAPH[i],"CA",2)==0) {

      for (j=0;j<3;++j) {
	for (k=0;k<3;++k) {
	  if (flag[j][k]==-1) {
	    iniflag=ON;
	  }
	  else {
	    iniflag=OFF;
	    break;
	  }
	}
	if (iniflag==OFF) {
	  break;
	}
      }

      for (j=0;j<3;++j) {
	for (k=1;k>=0;--k) {
	  if (flag[j][k]!=-1) {
	    flag[j][k+1]=i;
	    if (k+1==2) {
	      ++(*num_angl);
	      num_anglj[j]=(*num_angl);
	      index_angl=(int *)gcerealloc(index_angl,sizeof(int)*(*num_angl)*3);
	    }
	    break;
	  }
	}
	if (j>0 && flag[j-1][0]!=-1 && flag[j][0]==-1 && flag[j-1][0]!=i) {
	  flag[j][0]=i;
	}
      }

      if (iniflag==ON) flag[0][0]=i;
      iniflag=OFF;

      for (j=0;j<3;++j) {
	if (flag[j][2]!=-1) {
	  index_angl[(num_anglj[j]-1)*3]=flag[j][0];
	  index_angl[(num_anglj[j]-1)*3+1]=flag[j][1];
	  index_angl[(num_anglj[j]-1)*3+2]=flag[j][2];
	  flag[j][0]=-1;flag[j][1]=-1;flag[j][2]=-1;
	}
      }
    }
  }

  return index_angl;
}

int *set_DIHE_param(int *num_dihe,double *refcrd) {
  int i,j,k;
  int iniflag=OFF;
  int flag[4][4];
  int num_dihej[4];
  int numatom;

  int *index_dihe;

  index_dihe=(int *)gcemalloc(sizeof(int)*4);

  numatom=AP.NATOM;

  for (j=0;j<4;++j) 
    for (k=0;k<4;++k) 
      flag[j][k]=-1;

  *num_dihe=0;
  for (i=0;i<numatom;++i) {
    if (strncmp(AP.IGRAPH[i],"CA",2)==0) {

      for (j=0;j<4;++j) {
	for (k=0;k<4;++k) {
	  if (flag[j][k]==-1) {
	    iniflag=ON;
	  }
	  else {
	    iniflag=OFF;
	    break;
	  }
	}
	if (iniflag==OFF) {
	  break;
	}
      }

      for (j=0;j<4;++j) {
	for (k=2;k>=0;--k) {
	  if (flag[j][k]!=-1) {
	    flag[j][k+1]=i;
	    if (k+1==3) {
	      ++(*num_dihe);
	      num_dihej[j]=(*num_dihe);
	      index_dihe=(int *)gcerealloc(index_dihe,sizeof(int)*(*num_dihe)*4);
	    }
	    break;
	  }
	}
	if (j>0 && flag[j-1][0]!=-1 && flag[j][0]==-1 && flag[j-1][0]!=i) {
	  flag[j][0]=i;
	}
      }

      if (iniflag==ON) flag[0][0]=i;
      iniflag=OFF;

      for (j=0;j<4;++j) {
	if (flag[j][3]!=-1) {
	  index_dihe[(num_dihej[j]-1)*4]=flag[j][0];
	  index_dihe[(num_dihej[j]-1)*4+1]=flag[j][1];
	  index_dihe[(num_dihej[j]-1)*4+2]=flag[j][2];
	  index_dihe[(num_dihej[j]-1)*4+3]=flag[j][3];
	  flag[j][0]=-1;flag[j][1]=-1;flag[j][2]=-1;flag[j][3]=-1;
	}
      }
    }
  }

  return index_dihe;
}


