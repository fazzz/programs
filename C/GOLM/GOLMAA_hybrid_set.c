
#include <stdio.h>
#include <math.h>

//#include "GOLMAA_hybrid.h"
#include "GOLMAA_hybrid_set.h"

#include "PTL.h"
#include "TOPO.h"
#include "NC.h"

#include "MB.h"
#include "EF.h"

int GOLMAA_hybrid_ff_set_calcff(struct potential_GOLMAA_hybrid *ene, double *refcrd, int numatom,int numres,int *non_bonding_index, int numnonbond){
  int i,j,k,nc;
  int ii,jj,temp;
  int numca;

  double atom1[3],atom2[3],atom3[3],atom4[3];
  double criteria=criteria_hybrid;
  double length,len6,len12;
  double *indexcnb_cradii;
  int **ncmap,*indexncb,numnc;
  int **bp_f,*numb,**nb_matrix;

  double pi;

  pi=acos(-1.0);

  ncmap=GOLMAA_hybrid_ff_set_make_native_contact(refcrd,criteria,&((*ene).numNC),numatom,numres);

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
  (*ene).ep_natatt=ep_natatt_hybrid/4.184;

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

  //  (*ene).NotNC_index=(int *)gcemalloc(sizeof(int)*2);
  //  k=0;
  //  for (i=0;i<numatom;++i) {
  //    for (j=i+1;j<numatom;++j) {
  //      if (nb_matrix[i][j]==-1 && ncmap[i][j]!=0 && 
  //	  strncmp(AP.IGRAPH[i],"H",1)!=0 && strncmp(AP.IGRAPH[j],"H",1)!=0) {
  //	++k;
  //	(*ene).NotNC_index=(int *)gcerealloc((*ene).NotNC_index,sizeof(int)*2*k);
  //	(*ene).NotNC_index[(k-1)*2+0]=i;
  //	(*ene).NotNC_index[(k-1)*2+1]=j;
  //      }
  //    }
  //  }
  //  (*ene).numNotNC=k;
  ///////////////////////////////////////////////////////////////////////////////////

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

  (*ene).ep_repul=ep_repul_hybrid/4.184;
  (*ene).cradii_repul=cradii_repul_hybrid;
  len12=(*ene).cradii_repul;
  for (i=0;i<11;++i) len12=len12*(*ene).cradii_repul;
  (*ene).ALJ_repul=(*ene).ep_repul*len12;

  (*ene).p_natatt=(double *)gcemalloc(sizeof(double)*numatom);
  (*ene).p_repul=(double *)gcemalloc(sizeof(double)*numatom);

  (*ene).f_natatt=(double **)gcemalloc(sizeof(double *)*numatom);
  (*ene).f_repul=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) {
    (*ene).f_natatt[i]=(double *)gcemalloc(sizeof(double)*3);
    (*ene).f_repul[i]=(double *)gcemalloc(sizeof(double)*3);
  }

  (*ene).f_t=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) (*ene).f_t[i]=(double *)gcemalloc(sizeof(double)*3);
}

int GOLMAA_hybrid_ff_set_calcff_wtune(struct potential_GOLMAA_hybrid *ene, double *refcrd, int numatom,int numres,int *non_bonding_index, int numnonbond, double ep){
  int i,j,k,nc;
  int ii,jj,temp;
  int numca;

  double atom1[3],atom2[3],atom3[3],atom4[3];
  double criteria=criteria_hybrid;
  double length,len6,len12;
  double *indexcnb_cradii;
  int **ncmap,*indexncb,numnc;
  int **bp_f,*numb,**nb_matrix;

  double pi;

  pi=acos(-1.0);

  ncmap=GOLMAA_hybrid_ff_set_make_native_contact(refcrd,criteria,&((*ene).numNC),numatom,numres);

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
  (*ene).ep_natatt=ep/4.184;

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

  //  (*ene).NotNC_index=(int *)gcemalloc(sizeof(int)*2);
  //  k=0;
  //  for (i=0;i<numatom;++i) {
  //    for (j=i+1;j<numatom;++j) {
  //      if (nb_matrix[i][j]==-1 && ncmap[i][j]!=0 && 
  //	  strncmp(AP.IGRAPH[i],"H",1)!=0 && strncmp(AP.IGRAPH[j],"H",1)!=0) {
  //	++k;
  //	(*ene).NotNC_index=(int *)gcerealloc((*ene).NotNC_index,sizeof(int)*2*k);
  //	(*ene).NotNC_index[(k-1)*2+0]=i;
  //	(*ene).NotNC_index[(k-1)*2+1]=j;
  //      }
  //    }
  //  }
  //  (*ene).numNotNC=k;
  ///////////////////////////////////////////////////////////////////////////////////

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

  (*ene).ep_repul=ep_repul_hybrid/4.184;
  (*ene).cradii_repul=cradii_repul_hybrid;
  len12=(*ene).cradii_repul;
  for (i=0;i<11;++i) len12=len12*(*ene).cradii_repul;
  (*ene).ALJ_repul=(*ene).ep_repul*len12;

  (*ene).p_natatt=(double *)gcemalloc(sizeof(double)*numatom);
  (*ene).p_repul=(double *)gcemalloc(sizeof(double)*numatom);

  (*ene).f_natatt=(double **)gcemalloc(sizeof(double *)*numatom);
  (*ene).f_repul=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) {
    (*ene).f_natatt[i]=(double *)gcemalloc(sizeof(double)*3);
    (*ene).f_repul[i]=(double *)gcemalloc(sizeof(double)*3);
  }

  (*ene).f_t=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) (*ene).f_t[i]=(double *)gcemalloc(sizeof(double)*3);
}

int **GOLMAA_hybrid_ff_set_make_native_contact(double *refcrd,double criteria,int *numNC,int numatom,int numres) {
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

  index_natatt=make_native_contact_list_aa_3_nadjacent(&numncaa,&numncres,refcrd,numatom,numres,criteria,ncmap,ncmapres,EXC);

  ncmapexc=(int **)gcemalloc(sizeof(int *)*numatom);
  for (i=0;i<numatom;++i) ncmapexc[i]=(int *)gcemalloc(sizeof(int)*numatom);
  for (i=0;i<numatom;++i) for (j=0;j<numatom;++j) ncmapexc[i][j]=0;

  dummy=0;
  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      if (ncmap[i][j]==0) {
	for (k=0;k<numatom;++k) {
	  if (k!=i && k!=j && strncmp(AP.IGRAPH[k],"H",1)!=0) {
	    for (l=0;l<3;++l) {
	      atomi[l]=refcrd[i*3+l];
	      atomj[l]=refcrd[j*3+l];
	      atomk[l]=refcrd[k*3+l];
	    }
	    dij=len(atomi,atomj);
	    dik=len(atomi,atomk);
	    djk=len(atomj,atomk);
	    ajik=ang(atomj,atomi,atomk);
	    aijk=ang(atomi,atomj,atomk);
	    ajki=ang(atomj,atomk,atomi);
	    if (((dik < dij)&& (ajik < 35.0*pi/180.0)) || ((djk < dij) && (aijk < 35.0*pi/180.0)) ) {
	      ncmapexc[i][j]=-1;
	      ++dummy;
	      break;
	    }
	  }
	}
      }
    }
  }

  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      if (ncmapexc[i][j]==-1 ) {
	ncmap[i][j]=-1;
      }
    }
  }

  *numNC=numncaa-dummy;

  return ncmap;
}

int GOLMAA_hybrid_ff_set_calcff_2(struct potential_GOLMAA_hybrid *ene, double *refcrd, int numatom,int numres,int *non_bonding_index, int numnonbond){
  int i,j,k,nc;
  int ii,jj,temp;
  int numca;

  double atom1[3],atom2[3],atom3[3],atom4[3];
  double criteria=criteria_hybrid;
  double length,len6,len12;
  double *indexcnb_cradii;
  int **ncmap,*indexncb,numnc;
  int **bp_f,*numb,**nb_matrix;

  double pi;

  pi=acos(-1.0);

  ncmap=GOLMAA_hybrid_ff_set_make_native_contact_2(refcrd,criteria,&((*ene).numNC),numatom,numres);

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
  (*ene).ep_natatt=ep_natatt_hybrid/4.184;

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

  //  (*ene).NotNC_index=(int *)gcemalloc(sizeof(int)*2);
  //  k=0;
  //  for (i=0;i<numatom;++i) {
  //    for (j=i+1;j<numatom;++j) {
  //      if (nb_matrix[i][j]==-1 && ncmap[i][j]!=0 && 
  //	  strncmp(AP.IGRAPH[i],"H",1)!=0 && strncmp(AP.IGRAPH[j],"H",1)!=0) {
  //	++k;
  //	(*ene).NotNC_index=(int *)gcerealloc((*ene).NotNC_index,sizeof(int)*2*k);
  //	(*ene).NotNC_index[(k-1)*2+0]=i;
  //	(*ene).NotNC_index[(k-1)*2+1]=j;
  //      }
  //    }
  //  }
  //  (*ene).numNotNC=k;
  ///////////////////////////////////////////////////////////////////////////////////

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

  (*ene).ep_repul=ep_repul_hybrid/4.184;
  (*ene).cradii_repul=cradii_repul_hybrid;
  len12=(*ene).cradii_repul;
  for (i=0;i<11;++i) len12=len12*(*ene).cradii_repul;
  (*ene).ALJ_repul=(*ene).ep_repul*len12;

  (*ene).p_natatt=(double *)gcemalloc(sizeof(double)*numatom);
  (*ene).p_repul=(double *)gcemalloc(sizeof(double)*numatom);

  (*ene).f_natatt=(double **)gcemalloc(sizeof(double *)*numatom);
  (*ene).f_repul=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) {
    (*ene).f_natatt[i]=(double *)gcemalloc(sizeof(double)*3);
    (*ene).f_repul[i]=(double *)gcemalloc(sizeof(double)*3);
  }

  (*ene).f_t=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) (*ene).f_t[i]=(double *)gcemalloc(sizeof(double)*3);
}

int GOLMAA_hybrid_ff_set_calcff_2_wtune(struct potential_GOLMAA_hybrid *ene, double *refcrd, int numatom,int numres,int *non_bonding_index, int numnonbond, double ep){
  int i,j,k,nc;
  int ii,jj,temp;
  int numca;

  double atom1[3],atom2[3],atom3[3],atom4[3];
  double criteria=criteria_hybrid;
  double length,len6,len12;
  double *indexcnb_cradii;
  int **ncmap,*indexncb,numnc;
  int **bp_f,*numb,**nb_matrix;

  double pi;

  pi=acos(-1.0);

  ncmap=GOLMAA_hybrid_ff_set_make_native_contact_2(refcrd,criteria,&((*ene).numNC),numatom,numres);

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
  (*ene).ep_natatt=ep/4.184;

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

  (*ene).ep_repul=ep_repul_hybrid/4.184;
  (*ene).cradii_repul=cradii_repul_hybrid;
  len12=(*ene).cradii_repul;
  for (i=0;i<11;++i) len12=len12*(*ene).cradii_repul;
  (*ene).ALJ_repul=(*ene).ep_repul*len12;

  (*ene).p_natatt=(double *)gcemalloc(sizeof(double)*numatom);
  (*ene).p_repul=(double *)gcemalloc(sizeof(double)*numatom);

  (*ene).f_natatt=(double **)gcemalloc(sizeof(double *)*numatom);
  (*ene).f_repul=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) {
    (*ene).f_natatt[i]=(double *)gcemalloc(sizeof(double)*3);
    (*ene).f_repul[i]=(double *)gcemalloc(sizeof(double)*3);
  }

  (*ene).f_t=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) (*ene).f_t[i]=(double *)gcemalloc(sizeof(double)*3);
}

int **GOLMAA_hybrid_ff_set_make_native_contact_2(double *refcrd,double criteria,int *numNC,int numatom,int numres) {
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

  index_natatt=make_native_contact_list_aa_3_nadjacent_2(&numncaa,&numncres,refcrd,numatom,numres,criteria,ncmap,ncmapres,EXC);

  ncmapexc=(int **)gcemalloc(sizeof(int *)*numatom);
  for (i=0;i<numatom;++i) ncmapexc[i]=(int *)gcemalloc(sizeof(int)*numatom);
  for (i=0;i<numatom;++i) for (j=0;j<numatom;++j) ncmapexc[i][j]=0;

  dummy=0;
  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      if (ncmap[i][j]==0) {
	for (k=0;k<numatom;++k) {
	  if (k!=i && k!=j && strncmp(AP.IGRAPH[k],"H",1)!=0) {
	    for (l=0;l<3;++l) {
	      atomi[l]=refcrd[i*3+l];
	      atomj[l]=refcrd[j*3+l];
	      atomk[l]=refcrd[k*3+l];
	    }
	    dij=len(atomi,atomj);
	    dik=len(atomi,atomk);
	    djk=len(atomj,atomk);
	    ajik=ang(atomj,atomi,atomk);
	    aijk=ang(atomi,atomj,atomk);
	    ajki=ang(atomj,atomk,atomi);
	    if (((dik < dij)&& (ajik < 35.0*pi/180.0)) || ((djk < dij) && (aijk < 35.0*pi/180.0)) ) {
	      ncmapexc[i][j]=-1;
	      ++dummy;
	      break;
	    }
	  }
	}
      }
    }
  }

  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      if (ncmapexc[i][j]==-1 ) {
	ncmap[i][j]=-1;
      }
    }
  }

  *numNC=numncaa-dummy;

  return ncmap;
}

int GOLMAA_hybrid_ff_set_calcff_3(struct potential_GOLMAA_hybrid *ene, double *refcrd, int numatom,int numres,int *non_bonding_index, int numnonbond){
  int i,j,k,nc;
  int ii,jj,temp;
  int numca;

  double atom1[3],atom2[3],atom3[3],atom4[3];
  double criteria=criteria_hybrid;
  double length,len6,len12;
  double *indexcnb_cradii;
  int **ncmap,*indexncb,numnc;
  int **bp_f,*numb,**nb_matrix;

  double pi;

  pi=acos(-1.0);

  ncmap=GOLMAA_hybrid_ff_set_make_native_contact_3(refcrd,criteria,&((*ene).numNC),numatom,numres);

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
  (*ene).ep_natatt=ep_natatt_hybrid/4.184;

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

  //  (*ene).NotNC_index=(int *)gcemalloc(sizeof(int)*2);
  //  k=0;
  //  for (i=0;i<numatom;++i) {
  //    for (j=i+1;j<numatom;++j) {
  //      if (nb_matrix[i][j]==-1 && ncmap[i][j]!=0 && 
  //	  strncmp(AP.IGRAPH[i],"H",1)!=0 && strncmp(AP.IGRAPH[j],"H",1)!=0) {
  //	++k;
  //	(*ene).NotNC_index=(int *)gcerealloc((*ene).NotNC_index,sizeof(int)*2*k);
  //	(*ene).NotNC_index[(k-1)*2+0]=i;
  //	(*ene).NotNC_index[(k-1)*2+1]=j;
  //      }
  //    }
  //  }
  //  (*ene).numNotNC=k;
  ///////////////////////////////////////////////////////////////////////////////////

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

  (*ene).ep_repul=ep_repul_hybrid/4.184;
  (*ene).cradii_repul=cradii_repul_hybrid;
  len12=(*ene).cradii_repul;
  for (i=0;i<11;++i) len12=len12*(*ene).cradii_repul;
  (*ene).ALJ_repul=(*ene).ep_repul*len12;

  (*ene).p_natatt=(double *)gcemalloc(sizeof(double)*numatom);
  (*ene).p_repul=(double *)gcemalloc(sizeof(double)*numatom);

  (*ene).f_natatt=(double **)gcemalloc(sizeof(double *)*numatom);
  (*ene).f_repul=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) {
    (*ene).f_natatt[i]=(double *)gcemalloc(sizeof(double)*3);
    (*ene).f_repul[i]=(double *)gcemalloc(sizeof(double)*3);
  }

  (*ene).f_t=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) (*ene).f_t[i]=(double *)gcemalloc(sizeof(double)*3);
}

int GOLMAA_hybrid_ff_set_calcff_3_wtune(struct potential_GOLMAA_hybrid *ene, double *refcrd, int numatom,int numres,int *non_bonding_index, int numnonbond, double ep){
  int i,j,k,nc;
  int ii,jj,temp;
  int numca;

  double atom1[3],atom2[3],atom3[3],atom4[3];
  double criteria=criteria_hybrid;
  double length,len6,len12;
  double *indexcnb_cradii;
  int **ncmap,*indexncb,numnc;
  int **bp_f,*numb,**nb_matrix;

  double pi;

  pi=acos(-1.0);

  ncmap=GOLMAA_hybrid_ff_set_make_native_contact_3(refcrd,criteria,&((*ene).numNC),numatom,numres);

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
  (*ene).ep_natatt=ep/4.184;

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

  //  (*ene).NotNC_index=(int *)gcemalloc(sizeof(int)*2);
  //  k=0;
  //  for (i=0;i<numatom;++i) {
  //    for (j=i+1;j<numatom;++j) {
  //      if (nb_matrix[i][j]==-1 && ncmap[i][j]!=0 && 
  //	  strncmp(AP.IGRAPH[i],"H",1)!=0 && strncmp(AP.IGRAPH[j],"H",1)!=0) {
  //	++k;
  //	(*ene).NotNC_index=(int *)gcerealloc((*ene).NotNC_index,sizeof(int)*2*k);
  //	(*ene).NotNC_index[(k-1)*2+0]=i;
  //	(*ene).NotNC_index[(k-1)*2+1]=j;
  //      }
  //    }
  //  }
  //  (*ene).numNotNC=k;
  ///////////////////////////////////////////////////////////////////////////////////

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

  (*ene).ep_repul=ep_repul_hybrid/4.184;
  (*ene).cradii_repul=cradii_repul_hybrid;
  len12=(*ene).cradii_repul;
  for (i=0;i<11;++i) len12=len12*(*ene).cradii_repul;
  (*ene).ALJ_repul=(*ene).ep_repul*len12;

  (*ene).p_natatt=(double *)gcemalloc(sizeof(double)*numatom);
  (*ene).p_repul=(double *)gcemalloc(sizeof(double)*numatom);

  (*ene).f_natatt=(double **)gcemalloc(sizeof(double *)*numatom);
  (*ene).f_repul=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) {
    (*ene).f_natatt[i]=(double *)gcemalloc(sizeof(double)*3);
    (*ene).f_repul[i]=(double *)gcemalloc(sizeof(double)*3);
  }

  (*ene).f_t=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) (*ene).f_t[i]=(double *)gcemalloc(sizeof(double)*3);
}

int **GOLMAA_hybrid_ff_set_make_native_contact_3(double *refcrd,double criteria,int *numNC,int numatom,int numres) {
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

  index_natatt=make_native_contact_list_aa_3_nadjacent(&numncaa,&numncres,refcrd,numatom,numres,criteria,ncmap,ncmapres,EXC);

  ncmapexc=(int **)gcemalloc(sizeof(int *)*numatom);
  for (i=0;i<numatom;++i) ncmapexc[i]=(int *)gcemalloc(sizeof(int)*numatom);
  for (i=0;i<numatom;++i) for (j=0;j<numatom;++j) ncmapexc[i][j]=0;

  dummy=0;
  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      if ( ncmap[i][j]==0 ) {
	for (k=0;k<numatom;++k) {
	  resi=PTL_resnum(i);
	  resj=PTL_resnum(j);
	  resk=PTL_resnum(k);
	  if (/*resi!=resk*/abs(resi-resk)>1) {
	    if (k!=i && k!=j && strncmp(AP.IGRAPH[k],"H",1)!=0) {
	      for (l=0;l<3;++l) {
		atomi[l]=refcrd[i*3+l];
		atomj[l]=refcrd[j*3+l];
		atomk[l]=refcrd[k*3+l];
	      }
	      dij=len(atomi,atomj);
	      dik=len(atomi,atomk);
	      djk=len(atomj,atomk);
	      ajik=ang(atomj,atomi,atomk);
	      aijk=ang(atomi,atomj,atomk);
	      ajki=ang(atomj,atomk,atomi);
	      if (((dik < dij) && (ajik < 35.0*pi/180.0))) {
		ncmapexc[i][j]=-1;
		++dummy;
		break;
	      }
	    }
	  }
	  if (/*resj!=resk*/abs(resj-resk)>1) {
	    if (k!=i && k!=j && strncmp(AP.IGRAPH[k],"H",1)!=0) {
	      for (l=0;l<3;++l) {
		atomi[l]=refcrd[i*3+l];
		atomj[l]=refcrd[j*3+l];
		atomk[l]=refcrd[k*3+l];
	      }
	      dij=len(atomi,atomj);
	      dik=len(atomi,atomk);
	      djk=len(atomj,atomk);
	      ajik=ang(atomj,atomi,atomk);
	      aijk=ang(atomi,atomj,atomk);
	      ajki=ang(atomj,atomk,atomi);
	      if (((djk < dij)&& (ajik < 35.0*pi/180.0)) ) {
		ncmapexc[i][j]=-1;
		++dummy;
		break;
	      }
	    }
	  }
	}
      }
    }
  }    

  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      if (ncmapexc[i][j]==-1 ) {
	ncmap[i][j]=-1;
      }
    }
  }

  *numNC=numncaa-dummy;

  return ncmap;
}



int GOLMAA_hybrid_ff_set_calcff_4_wtune(struct potential_GOLMAA_hybrid *ene, double *refcrd, int numatom,int numres,int *non_bonding_index, int numnonbond, double ep){
  int i,j,k,nc;
  int ii,jj,temp;
  int numca;

  double atom1[3],atom2[3],atom3[3],atom4[3];
  double criteria=criteria_hybrid;
  double length,len6,len12;
  double *indexcnb_cradii;
  int **ncmap,*indexncb,numnc;
  int **bp_f,*numb,**nb_matrix;

  double pi;

  pi=acos(-1.0);

  ncmap=GOLMAA_hybrid_ff_set_make_native_contact_4(refcrd,criteria,&((*ene).numNC),numatom,numres);

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
  (*ene).ep_natatt=ep/4.184;

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

  (*ene).ep_repul=ep_repul_hybrid/4.184;
  (*ene).cradii_repul=cradii_repul_hybrid;
  len12=(*ene).cradii_repul;
  for (i=0;i<11;++i) len12=len12*(*ene).cradii_repul;
  (*ene).ALJ_repul=(*ene).ep_repul*len12;

  (*ene).p_natatt=(double *)gcemalloc(sizeof(double)*numatom);
  (*ene).p_repul=(double *)gcemalloc(sizeof(double)*numatom);

  (*ene).f_natatt=(double **)gcemalloc(sizeof(double *)*numatom);
  (*ene).f_repul=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) {
    (*ene).f_natatt[i]=(double *)gcemalloc(sizeof(double)*3);
    (*ene).f_repul[i]=(double *)gcemalloc(sizeof(double)*3);
  }

  (*ene).f_t=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) (*ene).f_t[i]=(double *)gcemalloc(sizeof(double)*3);
}


int **GOLMAA_hybrid_ff_set_make_native_contact_4(double *refcrd,double criteria,int *numNC,int numatom,int numres) {
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

  index_natatt=make_native_contact_list_aa_3_nadjacent_2(&numncaa,&numncres,refcrd,numatom,numres,criteria,ncmap,ncmapres,EXC);

  ncmapexc=(int **)gcemalloc(sizeof(int *)*numatom);
  for (i=0;i<numatom;++i) ncmapexc[i]=(int *)gcemalloc(sizeof(int)*numatom);
  for (i=0;i<numatom;++i) for (j=0;j<numatom;++j) ncmapexc[i][j]=0;

  dummy=0;
  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      if ( ncmap[i][j]==0 ) {
	for (k=0;k<numatom;++k) {
	  resi=PTL_resnum(i);
	  resj=PTL_resnum(j);
	  resk=PTL_resnum(k);
	  if (/*resi!=resk*/abs(resi-resk)>1) {
	    if (k!=i && k!=j && strncmp(AP.IGRAPH[k],"H",1)!=0) {
	      for (l=0;l<3;++l) {
		atomi[l]=refcrd[i*3+l];
		atomj[l]=refcrd[j*3+l];
		atomk[l]=refcrd[k*3+l];
	      }
	      dij=len(atomi,atomj);
	      dik=len(atomi,atomk);
	      djk=len(atomj,atomk);
	      ajik=ang(atomj,atomi,atomk);
	      aijk=ang(atomi,atomj,atomk);
	      ajki=ang(atomj,atomk,atomi);
	      if (((dik < dij) && (ajik < 35.0*pi/180.0))) {
		ncmapexc[i][j]=-1;
		++dummy;
		break;
	      }
	    }
	  }
	  if (/*resj!=resk*/abs(resj-resk)>1) {
	    if (k!=i && k!=j && strncmp(AP.IGRAPH[k],"H",1)!=0) {
	      for (l=0;l<3;++l) {
		atomi[l]=refcrd[i*3+l];
		atomj[l]=refcrd[j*3+l];
		atomk[l]=refcrd[k*3+l];
	      }
	      dij=len(atomi,atomj);
	      dik=len(atomi,atomk);
	      djk=len(atomj,atomk);
	      ajik=ang(atomj,atomi,atomk);
	      aijk=ang(atomi,atomj,atomk);
	      ajki=ang(atomj,atomk,atomi);
	      if (((djk < dij)&& (ajik < 35.0*pi/180.0)) ) {
		ncmapexc[i][j]=-1;
		++dummy;
		break;
	      }
	    }
	  }
	}
      }
    }
  }    

  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      if (ncmapexc[i][j]==-1 ) {
	ncmap[i][j]=-1;
      }
    }
  }

  *numNC=numncaa-dummy;

  return ncmap;
}

int **GOLMAA_hybrid_ff_set_make_native_contact_5(double *refcrd,double criteria,int *numNC,int numatom,int numres,int nibnum) {
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

int GOLMAA_hybrid_ff_set_calcff_6_wtune(struct potential_GOLMAA_hybrid *ene, double *refcrd, int numatom,int numres,int *non_bonding_index, int numnonbond, double ep, int nibnum,double criteria){
  int i,j,k,nc;
  int ii,jj,temp;
  int numca;

  double atom1[3],atom2[3],atom3[3],atom4[3];
  //  double criteria=criteria_hybrid;
  double length,len6,len12;
  double *indexcnb_cradii;
  int **ncmap,*indexncb,numnc;
  int **bp_f,*numb,**nb_matrix;

  double pi;

  pi=acos(-1.0);

  ncmap=GOLMAA_hybrid_ff_set_make_native_contact_6(refcrd,criteria,&((*ene).numNC),numatom,numres,nibnum);

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
  (*ene).ep_natatt=ep/4.184;

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

  (*ene).ep_repul=ep_repul_hybrid/4.184;
  (*ene).cradii_repul=cradii_repul_hybrid;
  len12=(*ene).cradii_repul;
  for (i=0;i<11;++i) len12=len12*(*ene).cradii_repul;
  (*ene).ALJ_repul=(*ene).ep_repul*len12;

  (*ene).p_natatt=(double *)gcemalloc(sizeof(double)*numatom);
  (*ene).p_repul=(double *)gcemalloc(sizeof(double)*numatom);

  (*ene).f_natatt=(double **)gcemalloc(sizeof(double *)*numatom);
  (*ene).f_repul=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) {
    (*ene).f_natatt[i]=(double *)gcemalloc(sizeof(double)*3);
    (*ene).f_repul[i]=(double *)gcemalloc(sizeof(double)*3);
  }

  (*ene).f_t=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) (*ene).f_t[i]=(double *)gcemalloc(sizeof(double)*3);
}


int **GOLMAA_hybrid_ff_set_make_native_contact_6(double *refcrd,double criteria,int *numNC,int numatom,int numres,int nibnum) {
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

  ncmapexc=(int **)gcemalloc(sizeof(int *)*numatom);
  for (i=0;i<numatom;++i) ncmapexc[i]=(int *)gcemalloc(sizeof(int)*numatom);
  for (i=0;i<numatom;++i) for (j=0;j<numatom;++j) ncmapexc[i][j]=0;
  
  dummy=0;
  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      for (k=j+1;k<numatom;++k) {
	if ( strncmp(AP.IGRAPH[i],"H",1)!=0 && strncmp(AP.IGRAPH[j],"H",1)!=0 && strncmp(AP.IGRAPH[k],"H",1)!=0  ) {
	  resi=PTL_resnum(i);
	  resj=PTL_resnum(j);
	  resk=PTL_resnum(k);
	  if ( ncmap[i][j]==0 &&  ncmap[i][k]==0 ) {
	    for (l=0;l<3;++l) {
	      atomi[l]=refcrd[i*3+l];
	      atomj[l]=refcrd[j*3+l];
	      atomk[l]=refcrd[k*3+l];
	    }
	    dij=len(atomi,atomj);
	    dik=len(atomi,atomk);
	    ajik=ang(atomj,atomi,atomk);
	    if (ajik < 35.0*pi/180.0) {
	      if (dik < dij) {
		ncmap/*exc*/[i][j]=-1;
		++dummy;
	      }
	      else {
		ncmap/*exc*/[i][k]=-1;
		++dummy;
	      }
	    }
	    if ( ncmap[i][j]==0 &&  ncmap[j][k]==0 ) {
	      for (l=0;l<3;++l) {
		atomi[l]=refcrd[i*3+l];
		atomj[l]=refcrd[j*3+l];
		atomk[l]=refcrd[k*3+l];
	      }
	      dij=len(atomi,atomj);
	      djk=len(atomj,atomk);
	      aijk=ang(atomi,atomj,atomk);
	      if (aijk < 35.0*pi/180.0) {
		if (dij < djk) {
		  ncmap/*exc*/[j][k]=-1;
		  ++dummy;
		}
		else {
		  ncmap/*exc*/[i][j]=-1;
		  ++dummy;
		}
	      }
	    }
	    if ( ncmap[i][k]==0 &&  ncmap[j][k]==0 ) {
	      for (l=0;l<3;++l) {
		atomi[l]=refcrd[i*3+l];
		atomj[l]=refcrd[j*3+l];
		atomk[l]=refcrd[k*3+l];
	      }
	      dik=len(atomi,atomk);
	      djk=len(atomj,atomk);
	      ajki=ang(atomj,atomk,atomi);
	      if (ajki < 35.0*pi/180.0) {
		if (dik < djk) {
		  ncmap/*exc*/[j][k]=-1;
		  ++dummy;
		}
		else {
		  ncmap/*exc*/[i][k]=-1;
		  ++dummy;
		}
	      }
	    }
	  }
	}
      }
    }
  }

  /*****************************************/
  /* for (i=0;i<numatom;++i) {		   */
  /*   for (j=i+1;j<numatom;++j) {	   */
  /*     if (ncmapexc[i][j]==-1 ) {	   */
  /* 	ncmap[i][j]=-1;			   */
  /*     }				   */
  /*   }				   */
  /* }					   */
  /*****************************************/

  *numNC=numncaa-dummy;

  return ncmap;
}
