
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "RBD.h"
#include "LA.h"
#include "EF.h"

double RBD_setcom(double *crd,double *mas,int numatom,double *com){
  int i,j;
  double summass=0.0;;

  for (i=0;i<3;++i)
    com[i]=0.0;

  for (i=0;i<numatom;++i) {
    summass+=mas[i];
    for (j=0;j<3;++j) {
      com[j]+=mas[i]*crd[i*3+j];
    }
  }
  for (i=0;i<3;++i)
    com[i]=com[i]/summass;

  return summass;
}

double RBD_removecom(double *crd,double *mas,int numatom){
  int i,j;
  double *com;

  com=(double *)ecalloc(sizeof(double),3);
  RBD_setcom(crd,mas,numatom,com);
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j)
      crd[i*3+j]-=com[j];

}

int RBD_setmassweightedcrd(double *crd,double *mas,int numatom,double *mwcrd){
  int i,j;

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j)
      mwcrd[i*3+j]=sqrt(mas[i])*crd[i*3+j];
}

int RBD_setInertia(double *crd,double *crd_org, double *mass,int numatom,double *Inertia){
  int i,j,k;
  double *crd_temp,*mwcrd,*mwcrd_temp,*mwcrd_tempx,*mwcrd_tempxT,*inertia_temp,com[3];

  crd_temp     = emalloc(sizeof(double)*numatom*3);
  mwcrd        = emalloc(sizeof(double)*numatom*3);
  mwcrd_temp   = emalloc(sizeof(double)*3);
  mwcrd_tempx  = emalloc(sizeof(double)*3*3);
  mwcrd_tempxT = emalloc(sizeof(double)*3*3);
  inertia_temp = emalloc(sizeof(double)*3*3);

  RBD_setcom(crd,mass,numatom,com);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j)
      crd_temp[i*3+j]=crd[i*3+j]-crd_org[j];

  RBD_setmassweightedcrd(crd_temp,mass,numatom,mwcrd);
  msetzero(Inertia,3);
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j)
      mwcrd_temp[j]=mwcrd[i*3+j];
    v_product(mwcrd_temp,mwcrd_tempx);
    mtrans(mwcrd_tempx,mwcrd_tempxT,3);
    mmult(mwcrd_tempxT,mwcrd_tempx,inertia_temp,3);
    for (j=0;j<3;++j)
      for (k=0;k<3;++k)
	Inertia[j*3+k]+=inertia_temp[j*3+k];
  }

  free(crd_temp);
  free(mwcrd);
  free(mwcrd_temp);
  free(mwcrd_tempx);
  free(mwcrd_tempxT);
  free(inertia_temp);

  return 0;
}

void RBD_trans_body_coordinate(int natom1,int natom2,int natom3,int numatom,double *coord,double *Rotmat, double *coord_body){
  int i,j;
  
  double atom1[3],atom2[3],atom3[3];
  double *coord_temp;
  double ii[3],jj[3],kk[3];  
  double d_3_1=0.0,d_2_1=0.0;
  double Sn3_1_2=0.0,Cs3_1_2=0.0;

  coord_temp=(double *)emalloc(sizeof(double)*numatom*3);
  for(i=0;i<3;++i){
    atom1[i]=coord[natom1*3+i];
    atom2[i]=coord[natom2*3+i];
    if (natom3==-1)
      atom3[i]=0.0;
    else
      atom3[i]=coord[natom3*3+i];
  }

  for(i=0;i<3;++i){
    kk[i]    = atom3[i]-atom1[i];
    d_3_1   += (atom3[i]-atom1[i])*(atom3[i]-atom1[i]);
    d_2_1   += (atom2[i]-atom1[i])*(atom2[i]-atom1[i]);
    Cs3_1_2 += (atom3[i]-atom1[i])*(atom2[i]-atom1[i]);
  }

  d_3_1 = sqrt(d_3_1);
  d_2_1 = sqrt(d_2_1);
  Cs3_1_2 = Cs3_1_2/(d_3_1*d_2_1);
  Sn3_1_2 = 1.0-Cs3_1_2*Cs3_1_2;
  Sn3_1_2 = sqrt(Sn3_1_2);
  
  for(i=0;i<3;++i)
    kk[i]=kk[i]/d_3_1;

  jj[0]=((atom3[1]-atom1[1])*(atom2[2]-atom1[2])-(atom3[2]-atom1[2])*(atom2[1]-atom1[1]))/(d_3_1*d_2_1*Sn3_1_2);
  jj[1]=((atom3[2]-atom1[2])*(atom2[0]-atom1[0])-(atom3[0]-atom1[0])*(atom2[2]-atom1[2]))/(d_3_1*d_2_1*Sn3_1_2);
  jj[2]=((atom3[0]-atom1[0])*(atom2[1]-atom1[1])-(atom3[1]-atom1[1])*(atom2[0]-atom1[0]))/(d_3_1*d_2_1*Sn3_1_2);

  ii[0]=jj[1]*kk[2]-jj[2]*kk[1];
  ii[1]=jj[2]*kk[0]-jj[0]*kk[2];
  ii[2]=jj[0]*kk[1]-jj[1]*kk[0];

  for(i=0;i<numatom;++i)
    for(j=0;j<3;++j)
      coord_temp[i*3+j] = coord[(natom1+i)*3+j]-atom1[j];
  
  for(i=0;i<3;++i){
    Rotmat[0*3+i]=ii[i];
    Rotmat[1*3+i]=jj[i];
    Rotmat[2*3+i]=kk[i];
  }

  for(i=0;i<numatom;++i){
    coord_body[i*3+0]=Rotmat[0*3+0]*coord_temp[i*3+0]+Rotmat[0*3+1]*coord_temp[i*3+1]+Rotmat[0*3+2]*coord_temp[i*3+2];
    coord_body[i*3+1]=Rotmat[1*3+0]*coord_temp[i*3+0]+Rotmat[1*3+1]*coord_temp[i*3+1]+Rotmat[1*3+2]*coord_temp[i*3+2];
    coord_body[i*3+2]=Rotmat[2*3+0]*coord_temp[i*3+0]+Rotmat[2*3+1]*coord_temp[i*3+1]+Rotmat[2*3+2]*coord_temp[i*3+2];
  }
}

void RBD_set_trans_mat(int natom1,int natom2,double *coord,double *Rotmat1, double *Rotmat2,double *Transmat1_2){
  int i,j,k;
  double *Rotmat1_2;
  double vect[3],vect2[3],*vectprod;

  Rotmat1_2=(double *)ecalloc(sizeof(double),3*3);
  vectprod=(double *)ecalloc(sizeof(double),3*3);

  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      Rotmat1_2[i*3+j] = 0.0;
  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      for (k=0;k<3;++k)
	Rotmat1_2[i*3+j]+=Rotmat1[i*3+k]*Rotmat2[j*3+k];

  for(i=0;i<3;++i)
    vect[i]=coord[natom1*3+i]-coord[natom2*3+i];

  for(i=0;i<3;++i)
    vect2[i]=0.0;
  for(i=0;i<3;++i)
    for(j=0;j<3;++j)
      vect2[i]+=Rotmat1[i*3+j]*vect[j];

  for(i=0;i<3;++i)
    for(j=0;j<3;++j)
      Transmat1_2[i*6+j]=Rotmat1_2[i*3+j];

  for(i=3;i<6;++i)
    for(j=3;j<6;++j)
      Transmat1_2[i*6+j]
	=Rotmat1_2[i-3*3+j-3];

  for(i=3;i<6;++i)
    for(j=0;j<3;++j)
      Transmat1_2[i*6+j]=0.0;

  vectprod[0*3+0]= 0.0;
  vectprod[0*3+1]=-vect2[2];
  vectprod[0*3+2]= vect2[1];
  vectprod[1*3+0]= vect2[2];
  vectprod[1*3+1]= 0.0;
  vectprod[1*3+2]=-vect2[0];
  vectprod[2*3+0]=-vect2[1];
  vectprod[2*3+1]= vect2[0];
  vectprod[2*3+2]= 0.0;

  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      Transmat1_2[i*6+(j+3)] = 0.0;

  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      for (k=0;k<3;++k)
	Transmat1_2[i*6+(j+3)] += vectprod[i*3+k]*Rotmat1_2[k*3+j];

  free(Rotmat1_2);
  free(vectprod);

}

int RMD_readRB_numbody(FILE *inputfile,RigidBodyColl RBC){
  int numbody;
  fscanf(inputfile,"%d",&numbody);
  return numbody;
}

int RBD_raedRBdata(FILE *inputfile,RigidBodyColl RBC){
  int i,j;
  
  for(i=0;i<RBC.numbody;++i)
    fscanf(inputfile,"%d",&RBC.RB[i].num_origin);
  for(i=0;i<RBC.numbody;++i)
    fscanf(inputfile,"%d",&RBC.RB[i].flag_terminal);
  for(i=0;i<RBC.numbody;++i)
    fscanf(inputfile,"%d",&RBC.RB[i].numatom_body);
  for(i=0;i<RBC.numbody;++i)
    fscanf(inputfile,"%d",&RBC.RB[i].num_branch);
  for(i=0;i<RBC.numbody;++i)
    fscanf(inputfile,"%d",&RBC.RB[i].hingmat);
  for(i=0;i<RBC.numbody;++i)
    for(j=0;j<RBC.RB[i].num_branch;++j)
      fscanf(inputfile,"%d",&RBC.RB[i].num_terminal[j]);
  for (i=0;i<RBC.numbody;++i)
    fscanf(inputfile, "%d", &RBC.RB[i].num_parent_body);
  for (i=0;i<RBC.numbody;++i)
    for(j=0;j<RBC.RB[i].num_branch;++j)
      fscanf(inputfile, "%d", &RBC.RB[i].num_child_body[j]);  
  for (i=0;i<RBC.numbody;++i)
    fscanf(inputfile, "%d", &RBC.indexofABAcyc[i]);
}

/************************************************************************/
/* int RBD_diag_Inertia(double *Inertia,double *eigenval,int numatom) { */
/*   int i,j;							        */
/*   char jobz='V';						        */
/*   char uplo='U';						        */
/*   double *w,*work,*Inertia_diag;				        */
/*   long int n,lda,lwork,info;					        */
/*   double var;						        */
/* 								        */
/*   n     = 3;							        */
/*   lda   = 3;							        */
/*   lwork = 3*3;						        */
/*   w     = (double *)ecalloc(sizeof(double),3*3);		        */
/*   work  = (double *)ecalloc(sizeof(double),3*3);		        */
/*   Inertia_diag = (double *)ecalloc(sizeof(double),3*3);	        */
/*   mtrans(cov,covm,3*3);					        */
/*   dsyev_(&jobz,&uplo,&n,covm,&lda,w,work,&lwork,&info);	        */
/*   if (info!=0){						        */
/*     printf("error; eigen vector calculation!!\n info=%d\n",info);    */
/*     exit(1);							        */
/*   }								        */
/*   for (i=0;i<3*3;++i)					        */
/*     eigenval[3*3-1-i]=w[i];					        */
/*   for (i=0;i<3*3;++i)					        */
/*     for (j=0;j<3*3;++j)					        */
/*       cov[i*3*3+j]=covm[i*3*3+j];				        */
/* 								        */
/*   return 1;							        */
/* }								        */
/************************************************************************/
