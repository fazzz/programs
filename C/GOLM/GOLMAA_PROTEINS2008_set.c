
#include <stdio.h>
#include <math.h>

#include "GOLMAA_PROTEINS2008_set.h"

#include "PTL.h"
#include "TOPO.h"
#include "NC.h"
#include "FFL.h"

#include "MB.h"
#include "EF.h"

int GOLMAA_PROTEINS2008_ff_set_calcff(struct potential_GOLMAA_PROTEINS2008 *ene, double *refcrd, int numatom,int numres,int *non_bonding_index, int numnonbond, double ep, int nibnum,double criteria){
  int i,j,k,nc;
  int ii,jj,kk,ll,temp;
  int numca;
  int inpnumA;

  double atom1[3],atom2[3],atom3[3],atom4[3];
  double length,len6,len12;
  double *indexcnb_cradii;
  int **ncmap,*indexncb,numnc;
  int **bp_f,*numb,**nb_matrix;

  double pi;

  pi=acos(-1.0);

  ncmap=GOLMAA_PROTEINS2008_ff_set_make_native_contact(refcrd,criteria,&((*ene).numNC),numatom,numres,nibnum);

  nc=0;
  (*ene).NC_index=(int *)gcemalloc(sizeof(int)*(*ene).numNC*2);
  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      if (ncmap[i][j]==0) {
	(*ene).NC_index[nc*2+0]=i;
	(*ene).NC_index[nc*2+1]=j;
	++nc;
      }
    }
  }

  (*ene).cradii_natatt=(double *)gcemalloc(sizeof(double)*(*ene).numNC);
  (*ene).ALJ_natatt=(double *)gcemalloc(sizeof(double)*(*ene).numNC);
  (*ene).BLJ_natatt=(double *)gcemalloc(sizeof(double)*(*ene).numNC);
  for (i=0;i<(*ene).numNC;++i) {
    for (j=0;j<3;++j) {
      atom1[j]=refcrd[((*ene).NC_index[i*2])*3+j];
      atom2[j]=refcrd[((*ene).NC_index[i*2+1])*3+j];
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
  (*ene).ep_natatt=ep;

  ///////////////////////////////////////////////////////////////////////////////////
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
  nb_matrix=(int **)gcemalloc(sizeof(int *)*numatom);
  for (i=0;i<numatom;++i) nb_matrix[i]=(int *)gcemalloc(sizeof(int)*numatom);
  make_nb_matrix(bp_f,numb,3,nb_matrix,numatom);

  (*ene).NotNC_index=(int *)gcemalloc(sizeof(int)*2);
  k=0;
  for (i=0;i<numnonbond;++i) {
    ii=non_bonding_index[i*2+0];
    jj=non_bonding_index[i*2+1];
    if ( ncmap[ii][jj]!=0 && strncmp(AP.IGRAPH[ii],"H",1)!=0 && strncmp(AP.IGRAPH[jj],"H",1)!=0 ) {
      ++k;
      (*ene).NotNC_index=(int *)gcerealloc((*ene).NotNC_index,sizeof(int)*2*k);
      (*ene).NotNC_index[(k-1)*2+0]=ii;
      (*ene).NotNC_index[(k-1)*2+1]=jj;
    }
  }
  (*ene).numNotNC=k;

  (*ene).ep_repul=ep_repul_PROTEINS2008;
  (*ene).cradii_repul=cradii_repul_PROTEINS2008;
  len12=(*ene).cradii_repul;
  for (i=0;i<11;++i) len12=len12*(*ene).cradii_repul;
  (*ene).ALJ_repul=(*ene).ep_repul*len12;

  (*ene).p_natatt=(double *)gcemalloc(sizeof(double)*(*ene).numNC);
  (*ene).p_repul =(double *)gcemalloc(sizeof(double)*(*ene).numNotNC);

  (*ene).f_natatt=(double **)gcemalloc(sizeof(double *)*numatom);
  (*ene).f_repul =(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) {
    (*ene).f_natatt[i]=(double *)gcemalloc(sizeof(double)*3);
    (*ene).f_repul[i] =(double *)gcemalloc(sizeof(double)*3);
  }

  (*ene).f_t=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) (*ene).f_t[i]=(double *)gcemalloc(sizeof(double)*3);

  ////////////////////////////////////////////////////////////////////////////////////////////////

  (*ene).num_bond=AP.MBONA;
  (*ene).pairs_bond=(int **)gcemalloc(sizeof(int *)*(*ene).num_bond);
  for (i=0;i<(*ene).num_bond;++i) (*ene).pairs_bond[i]=(int *)gcemalloc(sizeof(int)*2);
  (*ene).bon_equ = (double *)gcemalloc(sizeof(double)*(*ene).num_bond);
  (*ene).p_b = (double *)gcemalloc(sizeof(double)*(*ene).num_bond);
  (*ene).f_b = (double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) (*ene).f_b[i]=(double *)gcemalloc(sizeof(double)*3);
  for (i=0;i<(*ene).num_bond;++i) {
    (*ene).pairs_bond[i][0]=abs(AP.BA[i][0])/3;
    (*ene).pairs_bond[i][1]=abs(AP.BA[i][1])/3;
  }
  (*ene).Kb=ep_bond_PROTEINS2008*epsilon_PROTEINS2008;
  for (i=0;i<(*ene).num_bond;++i) {
    ii=(*ene).pairs_bond[i][0];
    jj=(*ene).pairs_bond[i][1];
    for (j=0;j<3;++j) {
      atom1[j]=refcrd[ii*3+j];
      atom2[j]=refcrd[jj*3+j];
    }
    (*ene).bon_equ[i] = len(atom1,atom2);
  }

  (*ene).num_angl=AP.MTHETA;
  (*ene).pairs_angl=(int **)gcemalloc(sizeof(int *)*(*ene).num_angl);
  for (i=0;i<(*ene).num_angl;++i) (*ene).pairs_angl[i]=(int *)gcemalloc(sizeof(int)*3);
  (*ene).ang_equ = (double *)gcemalloc(sizeof(double)*(*ene).num_angl);
  (*ene).p_a = (double *)gcemalloc(sizeof(double)*(*ene).num_angl);
  (*ene).f_a = (double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) (*ene).f_a[i]=(double *)gcemalloc(sizeof(double)*3);
  for (i=0;i<(*ene).num_angl;++i) {
    (*ene).pairs_angl[i][0]=abs(AP.TA[i][0])/3;
    (*ene).pairs_angl[i][1]=abs(AP.TA[i][1])/3;
    (*ene).pairs_angl[i][2]=abs(AP.TA[i][2])/3;
  }
  (*ene).Ka=ep_angl_PROTEINS2008*epsilon_PROTEINS2008;
  for (i=0;i<(*ene).num_angl;++i) {
    ii=(*ene).pairs_angl[i][0];
    jj=(*ene).pairs_angl[i][1];
    kk=(*ene).pairs_angl[i][2];
    for (j=0;j<3;++j) {
      atom1[j]=refcrd[ii*3+j];
      atom2[j]=refcrd[jj*3+j];
      atom3[j]=refcrd[kk*3+j];
    }
    (*ene).ang_equ[i] = ang(atom1,atom2,atom3);
  }

  (*ene).num_dihe=AP.MPHIA;
  (*ene).pairs_dihe=(int **)gcemalloc(sizeof(int *)*(*ene).num_dihe);
  for (i=0;i<(*ene).num_dihe;++i) (*ene).pairs_dihe[i]=(int *)gcemalloc(sizeof(int)*4);
  (*ene).dih_equ = (double *)gcemalloc(sizeof(double)*(*ene).num_dihe);
  (*ene).p_d = (double *)gcemalloc(sizeof(double)*(*ene).num_dihe);
  (*ene).f_d = (double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) (*ene).f_d[i]=(double *)gcemalloc(sizeof(double)*3);
  for (i=0;i<(*ene).num_dihe;++i) {
    (*ene).pairs_dihe[i][0]=abs(AP.PA[i][0])/3;
    (*ene).pairs_dihe[i][1]=abs(AP.PA[i][1])/3;
    (*ene).pairs_dihe[i][2]=abs(AP.PA[i][2])/3;
    (*ene).pairs_dihe[i][3]=abs(AP.PA[i][3])/3;
  }
  (*ene).Kd1=ep_dih1_PROTEINS2008*epsilon_PROTEINS2008;
  (*ene).Kd2=ep_dih2_PROTEINS2008*epsilon_PROTEINS2008;
  (*ene).Ki=ep_impd_PROTEINS2008*epsilon_PROTEINS2008;

  (*ene).impindex=ffL_make_inpindex_A(&inpnumA);

  for (i=0;i<(*ene).num_dihe;++i) {
    ii=(*ene).pairs_dihe[i][0];
    jj=(*ene).pairs_dihe[i][1];
    kk=(*ene).pairs_dihe[i][2];
    ll=(*ene).pairs_dihe[i][3];
    for (j=0;j<3;++j) {
      atom1[j]=refcrd[ii*3+j];
      atom2[j]=refcrd[jj*3+j];
      atom3[j]=refcrd[kk*3+j];
      atom4[j]=refcrd[ll*3+j];
    }
    (*ene).dih_equ[i] = pick_dihed(atom1,atom2,atom3,atom4,0,0.0);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////

}

int GOLMAA_PROTEINS2008_ff_set_calcff_b(struct potential_GOLMAA_PROTEINS2008 *ene, double *refcrd, int numatom,int numres,int *non_bonding_index, int numnonbond, double ep, int nibnum,double criteria){
  int i,j,k,nc;
  int ii,jj,kk,ll,temp;
  int numca;
  int inpnumA;

  double atom1[3],atom2[3],atom3[3],atom4[3];
  double length,len6,len12;
  double *indexcnb_cradii;
  int **ncmap,*indexncb,numnc;
  int **bp_f,*numb,**nb_matrix;

  double pi;

  pi=acos(-1.0);

  ncmap=GOLMAA_PROTEINS2008_ff_set_make_native_contact(refcrd,criteria,&((*ene).numNC),numatom,numres,nibnum);

  nc=0;
  (*ene).NC_index=(int *)gcemalloc(sizeof(int)*(*ene).numNC*2);
  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      if (ncmap[i][j]==0) {
	(*ene).NC_index[nc*2+0]=i;
	(*ene).NC_index[nc*2+1]=j;
	++nc;
      }
    }
  }

  (*ene).cradii_natatt=(double *)gcemalloc(sizeof(double)*(*ene).numNC);
  (*ene).ALJ_natatt=(double *)gcemalloc(sizeof(double)*(*ene).numNC);
  (*ene).BLJ_natatt=(double *)gcemalloc(sizeof(double)*(*ene).numNC);
  for (i=0;i<(*ene).numNC;++i) {
    for (j=0;j<3;++j) {
      atom1[j]=refcrd[((*ene).NC_index[i*2])*3+j];
      atom2[j]=refcrd[((*ene).NC_index[i*2+1])*3+j];
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
  (*ene).ep_natatt=ep;

  ///////////////////////////////////////////////////////////////////////////////////
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
  nb_matrix=(int **)gcemalloc(sizeof(int *)*numatom);
  for (i=0;i<numatom;++i) nb_matrix[i]=(int *)gcemalloc(sizeof(int)*numatom);
  make_nb_matrix(bp_f,numb,3,nb_matrix,numatom);

  (*ene).NotNC_index=(int *)gcemalloc(sizeof(int)*2);
  k=0;
  for (i=0;i<numnonbond;++i) {
    ii=non_bonding_index[i*2+0];
    jj=non_bonding_index[i*2+1];
    if ( ncmap[ii][jj]!=0 && strncmp(AP.IGRAPH[ii],"H",1)!=0 && strncmp(AP.IGRAPH[jj],"H",1)!=0 ) {
      ++k;
      (*ene).NotNC_index=(int *)gcerealloc((*ene).NotNC_index,sizeof(int)*2*k);
      (*ene).NotNC_index[(k-1)*2+0]=ii;
      (*ene).NotNC_index[(k-1)*2+1]=jj;
    }
  }
  (*ene).numNotNC=k;

  (*ene).ep_repul=ep_repul_PROTEINS2008;
  (*ene).cradii_repul=cradii_repul_PROTEINS2008;
  len12=(*ene).cradii_repul;
  for (i=0;i<11;++i) len12=len12*(*ene).cradii_repul;
  (*ene).ALJ_repul=(*ene).ep_repul*len12;

  //  (*ene).p_natatt=(double *)gcemalloc(sizeof(double)*(*ene).numNC);
  //  (*ene).p_repul =(double *)gcemalloc(sizeof(double)*(*ene).numNotNC);

  (*ene).f_natatt=(double **)gcemalloc(sizeof(double *)*numatom);
  (*ene).f_repul =(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) {
    (*ene).f_natatt[i]=(double *)gcemalloc(sizeof(double)*3);
    (*ene).f_repul[i] =(double *)gcemalloc(sizeof(double)*3);
  }

  (*ene).f_t=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) (*ene).f_t[i]=(double *)gcemalloc(sizeof(double)*3);

  ////////////////////////////////////////////////////////////////////////////////////////////////

  (*ene).num_bond=AP.MBONA;
  (*ene).pairs_bond=(int **)gcemalloc(sizeof(int *)*(*ene).num_bond);
  for (i=0;i<(*ene).num_bond;++i) (*ene).pairs_bond[i]=(int *)gcemalloc(sizeof(int)*2);
  (*ene).bon_equ = (double *)gcemalloc(sizeof(double)*(*ene).num_bond);
  (*ene).p_b = (double *)gcemalloc(sizeof(double)*(*ene).num_bond);
  (*ene).f_b = (double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) (*ene).f_b[i]=(double *)gcemalloc(sizeof(double)*3);
  for (i=0;i<(*ene).num_bond;++i) {
    (*ene).pairs_bond[i][0]=abs(AP.BA[i][0])/3;
    (*ene).pairs_bond[i][1]=abs(AP.BA[i][1])/3;
  }
  (*ene).Kb=ep_bond_PROTEINS2008*epsilon_PROTEINS2008;
  for (i=0;i<(*ene).num_bond;++i) {
    ii=(*ene).pairs_bond[i][0];
    jj=(*ene).pairs_bond[i][1];
    for (j=0;j<3;++j) {
      atom1[j]=refcrd[ii*3+j];
      atom2[j]=refcrd[jj*3+j];
    }
    (*ene).bon_equ[i] = len(atom1,atom2);
  }

  (*ene).num_angl=AP.MTHETA;
  (*ene).pairs_angl=(int **)gcemalloc(sizeof(int *)*(*ene).num_angl);
  for (i=0;i<(*ene).num_angl;++i) (*ene).pairs_angl[i]=(int *)gcemalloc(sizeof(int)*3);
  (*ene).ang_equ = (double *)gcemalloc(sizeof(double)*(*ene).num_angl);
  (*ene).p_a = (double *)gcemalloc(sizeof(double)*(*ene).num_angl);
  (*ene).f_a = (double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) (*ene).f_a[i]=(double *)gcemalloc(sizeof(double)*3);
  for (i=0;i<(*ene).num_angl;++i) {
    (*ene).pairs_angl[i][0]=abs(AP.TA[i][0])/3;
    (*ene).pairs_angl[i][1]=abs(AP.TA[i][1])/3;
    (*ene).pairs_angl[i][2]=abs(AP.TA[i][2])/3;
  }
  (*ene).Ka=ep_angl_PROTEINS2008*epsilon_PROTEINS2008;
  for (i=0;i<(*ene).num_angl;++i) {
    ii=(*ene).pairs_angl[i][0];
    jj=(*ene).pairs_angl[i][1];
    kk=(*ene).pairs_angl[i][2];
    for (j=0;j<3;++j) {
      atom1[j]=refcrd[ii*3+j];
      atom2[j]=refcrd[jj*3+j];
      atom3[j]=refcrd[kk*3+j];
    }
    (*ene).ang_equ[i] = ang(atom1,atom2,atom3);
  }

  (*ene).num_dihe=AP.MPHIA;
  (*ene).pairs_dihe=(int **)gcemalloc(sizeof(int *)*(*ene).num_dihe);
  for (i=0;i<(*ene).num_dihe;++i) (*ene).pairs_dihe[i]=(int *)gcemalloc(sizeof(int)*4);
  (*ene).dih_equ = (double *)gcemalloc(sizeof(double)*(*ene).num_dihe);
  (*ene).p_d = (double *)gcemalloc(sizeof(double)*(*ene).num_dihe);
  (*ene).f_d = (double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) (*ene).f_d[i]=(double *)gcemalloc(sizeof(double)*3);
  for (i=0;i<(*ene).num_dihe;++i) {
    (*ene).pairs_dihe[i][0]=abs(AP.PA[i][0])/3;
    (*ene).pairs_dihe[i][1]=abs(AP.PA[i][1])/3;
    (*ene).pairs_dihe[i][2]=abs(AP.PA[i][2])/3;
    (*ene).pairs_dihe[i][3]=abs(AP.PA[i][3])/3;
  }
  (*ene).Kd1=ep_dih1_PROTEINS2008*epsilon_PROTEINS2008;
  (*ene).Kd2=ep_dih2_PROTEINS2008*epsilon_PROTEINS2008;
  (*ene).Ki=ep_impd_PROTEINS2008*epsilon_PROTEINS2008;

  (*ene).impindex=ffL_make_inpindex_A(&inpnumA);

  for (i=0;i<(*ene).num_dihe;++i) {
    ii=(*ene).pairs_dihe[i][0];
    jj=(*ene).pairs_dihe[i][1];
    kk=(*ene).pairs_dihe[i][2];
    ll=(*ene).pairs_dihe[i][3];
    for (j=0;j<3;++j) {
      atom1[j]=refcrd[ii*3+j];
      atom2[j]=refcrd[jj*3+j];
      atom3[j]=refcrd[kk*3+j];
      atom4[j]=refcrd[ll*3+j];
    }
    (*ene).dih_equ[i] = pick_dihed(atom1,atom2,atom3,atom4,0,0.0);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////

}


int **GOLMAA_PROTEINS2008_ff_set_make_native_contact(double *refcrd,double criteria,int *numNC,int numatom,int numres,int nibnum) {
  int i,j,k,l,d,na,ii,jj,nc;
  int resi,resj,resk;
  double pi;
  
  double dij,dik,djk,ajik,aijk,ajki;
  int dummy;
  
  int numncaa,numncres;
  int  *index_natatt,**ncmap,**ncmapexc,**ncmapres;
  double atom1[3],atom2[3],atomi[3],atomj[3],atomk[3];
  
  double vec[3];
  double length;

  pi=acos(-1.0);

  ncmap=(int **)gcemalloc(sizeof(int *)*numatom);
  for (i=0;i<numatom;++i) ncmap[i]=(int *)gcemalloc(sizeof(int)*numatom);
  ncmapres=(int **)gcemalloc(sizeof(int *)*numres);
  for (i=0;i<numres;++i) ncmapres[i]=(int *)gcemalloc(sizeof(int)*numres);
  
  index_natatt=make_native_contact_list_aa_wnibnum(&numncaa,&numncres,refcrd,numatom,numres,criteria,ncmap,ncmapres,EXC,nibnum);

  *numNC=numncaa;

  return ncmap;
}

