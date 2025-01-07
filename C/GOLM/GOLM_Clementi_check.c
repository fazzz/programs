#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "GOLM_Clementi_set.h"
#include "GOLM_Clementi_check.h"

#include "TOPO.h"
#include "LA.h"
#include "mymath.h"

#include "PTL.h"

#define UNIT 4.184070*100.0

#define ON 1
#define OFF 0

double GOLM_Clementi_calcff_check(char *inputfilename,char *reffilename,char *parmfilename,
				  int numspatom,double dx,
				  double f_natatt[3],double f_repul[3],double f_d[3],double f_a[3],double f_b[3],
				  double f_d1[4][3],double f_d2[4][3], int nums,
				  double f_natatt1[2][3],double f_natatt2[2][3], int numa) {
  int i,j,k,d;
  double *crd,*crddx,*crddy,*crddz,*refcrd,*refcrdAA;
  int numatom,numCAatom;
  struct potential_GOLM_Clementi e,edx,edy,edz;

  double x[3];

  double **f1,**f2,p_d,p_dx,p_dy,p_dz;
  double f_d3[4][3];

  double p_natatt,p_natatt_dx,p_natatt_dy,p_natatt_dz;
  double f_natatt3[2][3];
  int atomi,atomj;

  char *line;
  size_t len=0;

  FILE *inputfile,*reffile,*outputfile,*parmfile;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  j=0;
  for (i=0;i<numatom;++i) {
    if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
      ++j;
    }
  }
  numCAatom=j;

  crd=(double *)gcemalloc(sizeof(double)*numCAatom*3);
  refcrd=(double *)gcemalloc(sizeof(double)*numCAatom*3);
  refcrdAA=(double *)gcemalloc(sizeof(double)*numatom*3);

  crddx=(double *)gcemalloc(sizeof(double)*numCAatom*3);
  crddy=(double *)gcemalloc(sizeof(double)*numCAatom*3);
  crddz=(double *)gcemalloc(sizeof(double)*numCAatom*3);

  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&d);
  j=0;
  for (i=0;i<numatom;++i) {
    for (k=0;k<3;++k) fscanf(inputfile,"%lf",&x[k]);
    if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
      for (k=0;k<3;++k) crd[j*3+k]=x[k];
      ++j;
    }
  }
  fclose(inputfile);

  reffile=efopen(reffilename,"r");
  getline(&line,&len,reffile);
  fscanf(reffile,"%d",&d);
  j=0;
  for (i=0;i<numatom;++i) {
    for (k=0;k<3;++k) fscanf(reffile,"%lf",&refcrdAA[i*3+k]);
    if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
      for (k=0;k<3;++k) refcrd[j*3+k]=refcrdAA[i*3+k];
      ++j;
    }
  }
  fclose(reffile);

  for (i=0;i<numCAatom;++i) {
    for (j=0;j<3;++j) {
      crddx[i*3+j]=crd[i*3+j];
      crddy[i*3+j]=crd[i*3+j];
      crddz[i*3+j]=crd[i*3+j];
    }
  }

  crddx[numspatom*3]+=dx;
  crddy[numspatom*3+1]+=dx;
  crddz[numspatom*3+2]+=dx;

  GOLM_Clementi_ff_set_calcff(&e,refcrd,refcrdAA,numCAatom,numatom);
  GOLM_Clementi_ff_set_calcff(&edx,refcrd,refcrdAA,numCAatom,numatom);
  GOLM_Clementi_ff_set_calcff(&edy,refcrd,refcrdAA,numCAatom,numatom);
  GOLM_Clementi_ff_set_calcff(&edz,refcrd,refcrdAA,numCAatom,numatom);

  GOLM_Clementi_ff_calcff(crd,numCAatom,&e);
  GOLM_Clementi_ff_calcff(crddx,numCAatom,&edx);
  GOLM_Clementi_ff_calcff(crddy,numCAatom,&edy);
  GOLM_Clementi_ff_calcff(crddz,numCAatom,&edz);

  f_natatt[0]=-(edx.p_natatt_t-e.p_natatt_t)/dx*4.184070*100.0;
  f_natatt[1]=-(edy.p_natatt_t-e.p_natatt_t)/dx*4.184070*100.0;
  f_natatt[2]=-(edz.p_natatt_t-e.p_natatt_t)/dx*4.184070*100.0;

  f_repul[0]=-(edx.p_repul_t-e.p_repul_t)/dx*4.184070*100.0;
  f_repul[1]=-(edy.p_repul_t-e.p_repul_t)/dx*4.184070*100.0;
  f_repul[2]=-(edz.p_repul_t-e.p_repul_t)/dx*4.184070*100.0;

  f_d[0]=-(edx.p_d_t-e.p_d_t)/dx*4.184070*100.0;
  f_d[1]=-(edy.p_d_t-e.p_d_t)/dx*4.184070*100.0;
  f_d[2]=-(edz.p_d_t-e.p_d_t)/dx*4.184070*100.0;

  f_a[0]=-(edx.p_a_t-e.p_a_t)/dx*4.184070*100.0;
  f_a[1]=-(edy.p_a_t-e.p_a_t)/dx*4.184070*100.0;
  f_a[2]=-(edz.p_a_t-e.p_a_t)/dx*4.184070*100.0;

  f_b[0]=-(edx.p_b_t-e.p_b_t)/dx*4.184070*100.0;
  f_b[1]=-(edy.p_b_t-e.p_b_t)/dx*4.184070*100.0;
  f_b[2]=-(edz.p_b_t-e.p_b_t)/dx*4.184070*100.0;

  p_d=GOLM_Clementi_ff_calcDIHE_s(crd,numCAatom,f_d1,e.Kd1,e.Kd2,e.dih_equ,nums);

  for (i=0;i<numCAatom;++i) {
    for (j=0;j<3;++j) {
      crddx[i*3+j]=crd[i*3+j];
      crddy[i*3+j]=crd[i*3+j];
      crddz[i*3+j]=crd[i*3+j];
    }
  }

  crddx[nums*3]+=dx;
  crddy[nums*3+1]+=dx;
  crddz[nums*3+2]+=dx;

  p_dx=GOLM_Clementi_ff_calcDIHE_s(crddx,numCAatom,f_d3,e.Kd1,e.Kd2,e.dih_equ,nums);
  p_dy=GOLM_Clementi_ff_calcDIHE_s(crddy,numCAatom,f_d3,e.Kd1,e.Kd2,e.dih_equ,nums);
  p_dz=GOLM_Clementi_ff_calcDIHE_s(crddz,numCAatom,f_d3,e.Kd1,e.Kd2,e.dih_equ,nums);

  f_d2[0][0]=-(p_dx-p_d)/dx*4.184070*100.0;
  f_d2[0][1]=-(p_dy-p_d)/dx*4.184070*100.0;
  f_d2[0][2]=-(p_dz-p_d)/dx*4.184070*100.0;

  for (i=0;i<numCAatom;++i) {
    for (j=0;j<3;++j) {
      crddx[i*3+j]=crd[i*3+j];
      crddy[i*3+j]=crd[i*3+j];
      crddz[i*3+j]=crd[i*3+j];
    }
  }

  crddx[(nums+1)*3]+=dx;
  crddy[(nums+1)*3+1]+=dx;
  crddz[(nums+1)*3+2]+=dx;

  p_dx=GOLM_Clementi_ff_calcDIHE_s(crddx,numCAatom,f_d3,e.Kd1,e.Kd2,e.dih_equ,nums);
  p_dy=GOLM_Clementi_ff_calcDIHE_s(crddy,numCAatom,f_d3,e.Kd1,e.Kd2,e.dih_equ,nums);
  p_dz=GOLM_Clementi_ff_calcDIHE_s(crddz,numCAatom,f_d3,e.Kd1,e.Kd2,e.dih_equ,nums);

  f_d2[1][0]=-(p_dx-p_d)/dx*4.184070*100.0;
  f_d2[1][1]=-(p_dy-p_d)/dx*4.184070*100.0;
  f_d2[1][2]=-(p_dz-p_d)/dx*4.184070*100.0;

  for (i=0;i<numCAatom;++i) {
    for (j=0;j<3;++j) {
      crddx[i*3+j]=crd[i*3+j];
      crddy[i*3+j]=crd[i*3+j];
      crddz[i*3+j]=crd[i*3+j];
    }
  }

  crddx[(nums+2)*3]+=dx;
  crddy[(nums+2)*3+1]+=dx;
  crddz[(nums+2)*3+2]+=dx;

  p_dx=GOLM_Clementi_ff_calcDIHE_s(crddx,numCAatom,f_d3,e.Kd1,e.Kd2,e.dih_equ,nums);
  p_dy=GOLM_Clementi_ff_calcDIHE_s(crddy,numCAatom,f_d3,e.Kd1,e.Kd2,e.dih_equ,nums);
  p_dz=GOLM_Clementi_ff_calcDIHE_s(crddz,numCAatom,f_d3,e.Kd1,e.Kd2,e.dih_equ,nums);

  f_d2[2][0]=-(p_dx-p_d)/dx*4.184070*100.0;
  f_d2[2][1]=-(p_dy-p_d)/dx*4.184070*100.0;
  f_d2[2][2]=-(p_dz-p_d)/dx*4.184070*100.0;

  for (i=0;i<numCAatom;++i) {
    for (j=0;j<3;++j) {
      crddx[i*3+j]=crd[i*3+j];
      crddy[i*3+j]=crd[i*3+j];
      crddz[i*3+j]=crd[i*3+j];
    }
  }

  crddx[(nums+3)*3]+=dx;
  crddy[(nums+3)*3+1]+=dx;
  crddz[(nums+3)*3+2]+=dx;

  p_dx=GOLM_Clementi_ff_calcDIHE_s(crddx,numCAatom,f_d3,e.Kd1,e.Kd2,e.dih_equ,nums);
  p_dy=GOLM_Clementi_ff_calcDIHE_s(crddy,numCAatom,f_d3,e.Kd1,e.Kd2,e.dih_equ,nums);
  p_dz=GOLM_Clementi_ff_calcDIHE_s(crddz,numCAatom,f_d3,e.Kd1,e.Kd2,e.dih_equ,nums);

  f_d2[3][0]=-(p_dx-p_d)/dx*4.184070*100.0;
  f_d2[3][1]=-(p_dy-p_d)/dx*4.184070*100.0;
  f_d2[3][2]=-(p_dz-p_d)/dx*4.184070*100.0;

  p_natatt=GOLM_Clementi_ff_calcff_natatt_s(crd,numatom,e.index_natatt,e.num_natatt,e.ALJ_natatt,e.BLJ_natatt,e.ep_natatt,p_natatt,f_natatt1,numa);

  atomi=e.index_natatt[numa*2+0];
  atomj=e.index_natatt[numa*2+1];

  for (i=0;i<numCAatom;++i) {
    for (j=0;j<3;++j) {
      crddx[i*3+j]=crd[i*3+j];
      crddy[i*3+j]=crd[i*3+j];
      crddz[i*3+j]=crd[i*3+j];
    }
  }

  crddx[(atomi)*3]+=dx;
  crddy[(atomi)*3+1]+=dx;
  crddz[(atomi)*3+2]+=dx;

  p_natatt_dx=GOLM_Clementi_ff_calcff_natatt_s(crddx,numatom,(edx).index_natatt,(edx).num_natatt,(edx).ALJ_natatt,(edx).BLJ_natatt,(edx).ep_natatt,p_natatt,f_natatt3,numa);
  p_natatt_dy=GOLM_Clementi_ff_calcff_natatt_s(crddy,numatom,(edy).index_natatt,(edy).num_natatt,(edy).ALJ_natatt,(edy).BLJ_natatt,(edy).ep_natatt,p_natatt,f_natatt3,numa);
  p_natatt_dz=GOLM_Clementi_ff_calcff_natatt_s(crddz,numatom,(edz).index_natatt,(edz).num_natatt,(edz).ALJ_natatt,(edz).BLJ_natatt,(edz).ep_natatt,p_natatt,f_natatt3,numa);

  f_natatt2[0][0]=-(p_natatt_dx-p_natatt)/dx*4.184070*100.0;
  f_natatt2[0][1]=-(p_natatt_dy-p_natatt)/dx*4.184070*100.0;
  f_natatt2[0][2]=-(p_natatt_dz-p_natatt)/dx*4.184070*100.0;

  for (i=0;i<numCAatom;++i) {
    for (j=0;j<3;++j) {
      crddx[i*3+j]=crd[i*3+j];
      crddy[i*3+j]=crd[i*3+j];
      crddz[i*3+j]=crd[i*3+j];
    }
  }

  crddx[(atomj)*3]+=dx;
  crddy[(atomj)*3+1]+=dx;
  crddz[(atomj)*3+2]+=dx;

  p_natatt_dx=GOLM_Clementi_ff_calcff_natatt_s(crddx,numatom,(edx).index_natatt,(edx).num_natatt,(edx).ALJ_natatt,(edx).BLJ_natatt,(edx).ep_natatt,p_natatt,f_natatt3,numa);
  p_natatt_dy=GOLM_Clementi_ff_calcff_natatt_s(crddy,numatom,(edy).index_natatt,(edy).num_natatt,(edy).ALJ_natatt,(edy).BLJ_natatt,(edy).ep_natatt,p_natatt,f_natatt3,numa);
  p_natatt_dz=GOLM_Clementi_ff_calcff_natatt_s(crddz,numatom,(edz).index_natatt,(edz).num_natatt,(edz).ALJ_natatt,(edz).BLJ_natatt,(edz).ep_natatt,p_natatt,f_natatt3,numa);

  f_natatt2[1][0]=-(p_natatt_dx-p_natatt)/dx*4.184070*100.0;
  f_natatt2[1][1]=-(p_natatt_dy-p_natatt)/dx*4.184070*100.0;
  f_natatt2[1][2]=-(p_natatt_dz-p_natatt)/dx*4.184070*100.0;

}

double GOLM_Clementi_ff_calcDIHE_s(double *crd,int numatom,double f_d[4][3],double Kd1,double Kd2,double *dih_equ, int nums){
  int i,j,k,l;

  double fi[3],fj[3],fk[3],fl[3];
  double m[3],n[3],m_n[3],n_n[3],lm,ln;
  double vij[3],vkj[3],vkl[3];
  double lkj;
  double vijvkj,vklvkj;

  double dvdpsi;
  
  double atom[4][3];
  double dihed;
  double p_d_t=0.0;

  double pi;

  pi=acos(-1.0);

  for (i=0;i<4;++i) for (j=0;j<3;++j) f_d[i][j] = 0.0;
  
  i=nums;
  for (j=0;j<3;++j) {
    atom[0][j]=crd[i*3+j];
    atom[1][j]=crd[(i+1)*3+j];
    atom[2][j]=crd[(i+2)*3+j];
    atom[3][j]=crd[(i+3)*3+j];
  }

  for (j=0;j<3;++j) {
    vij[j] = atom[1][j]-atom[0][j];
    vkj[j] = atom[1][j]-atom[2][j];
    vkl[j] = atom[3][j]-atom[2][j];
  }
  lkj=sqrt(inprod(vkj,vkj,3));
  
  outprod(vij,vkj,m);
  outprod(vkj,vkl,n);
  lm=sqrt(inprod(m,m,3));
  ln=sqrt(inprod(n,n,3));
  for (j=0;j<3;++j) {
    m_n[j]=m[j]/lm;
    n_n[j]=n[j]/ln;
  }

  dihed=acos(inprod(m_n,n_n,3));
  if (inprod(vij,n,3)>0) dihed=-dihed;
  if (dihed<0.0) dihed=2.0*pi+dihed;

  vijvkj=inprod(vij,vkj,3);
  vklvkj=inprod(vkl,vkj,3);

  dvdpsi=Kd1*sin(dihed-dih_equ[i])+3.0*Kd2*sin(3.0*(dihed-dih_equ[i]));
  
  p_d_t = Kd1*(1.0-cos(dihed-dih_equ[i]))+Kd2*(1.0-cos(3.0*(dihed-dih_equ[i])));

  for (j=0;j<3;++j) {
    fi[j] = -dvdpsi*lkj*m[j]/(lm*lm);
    fl[j] =  dvdpsi*lkj*n[j]/(ln*ln);
    fj[j] = (-fi[j]+(vijvkj/(lkj*lkj))*fi[j]-(vklvkj/(lkj*lkj))*fl[j]);
    fk[j] = (-fl[j]-(vijvkj/(lkj*lkj))*fi[j]+(vklvkj/(lkj*lkj))*fl[j]);

    f_d[0][j] = fi[j]*4.184070*100.0;
    f_d[1][j] = fj[j]*4.184070*100.0;
    f_d[2][j] = fk[j]*4.184070*100.0;
    f_d[3][j] = fl[j]*4.184070*100.0;
  }

  return p_d_t;

}

double GOLM_Clementi_ff_calcff_natatt_s(double *crd, int numatom,int *index_numatt,int numatt,
					double *ALJ_natatt,double *BLJ_natatt,double ep_natatt,
					double p_natatt,double f_natatt[2][3],int nums) {
  int i,j,k;
  int atomi,atomj;

  double vec[3];
  double len,len2,len10,len12;
  double p12,p10,f;
  double p_natatt_t=0.0;


  for (i=0;i<2;++i) for (j=0;j<3;++j) f_natatt[i][j] = 0.0;

  i=nums;
  atomi=index_numatt[i*2+0];
  atomj=index_numatt[i*2+1];

  len2 = 0.0;
  for(j=0;j<3;++j){
    vec[j] = crd[atomi*3+j]-crd[atomj*3+j];
    len2 += vec[j]*vec[j];
  }
  len = sqrt(len2);
  len10=len2;
  len12=len2;
  for (j=0;j<4;++j)  len10 = len10*len2;
  for (j=0;j<5;++j)  len12 = len12*len2;
  p12 = ALJ_natatt[i]/len12;
  p10 = BLJ_natatt[i]/len10;
  p_natatt = ep_natatt*(5.0*p12-6.0*p10);
  for (j=0;j<3;++j) {
    f = ep_natatt*(12.0*5.0*p12-10.0*6.0*p10)/(len2)*vec[j]*4.184070*100.0;
    f_natatt[0][j] =  f;
    f_natatt[1][j] = -f;
  }

  return p_natatt;
}
