/***************************/
/* One Rigid Body Dynamics */
/***************************/
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "RBD.h"
#include "PT.h"
#include "EF.h"
#include "IO.h"
#include "QUA.h"
#include "LA.h"
#include "INTG.h"

void Lseof(double *vel,double *omg,double *Inertia,double Mass,double *acc,double *domg);
void Lcvec(double *crd,int numatom,double com[3], double *vec);
void LsetcoriolisCOM(double *omg, double *Inertia, double *corif);
void LsetmatCOM(double *Inertia,double Mass,double *mat);
double Lcenergy(double *vel ,double *omg ,double Mass,double *Inertia);
int Lwriteene(FILE *outfil,double ene);
int Lwritecrd(FILE *outfil,double *crd,int numatom);
void LreadAmberCoord(FILE *input, double *crd, int numatom);
int outputconf(FILE *outputfile,int numatom,double *data,int flag);

int main(int argc, char *argv[]) {
  int i,j,k,numstep,numatom,n;

  double dt;
  double *crd,*vec;
  double *Inertia,Mass,*mas,lenq;
  double ene;
  double *com;
  double *vel,*omg,*acc,*domg;
  double *rotmat;
  double omgp[3][5],omgc[3][5],qp[4][5],qc[4][5];
  double q[4],r[4],dq[4],roted[4],omg_temp[3];

  char *infilname1,*infilname2,*infilname3,*outfilname1,*outfilname2;
  FILE *infil,*infil2,*infil3,*outfil1,*outfil2;

  q[0]=1.0;q[1]=0.0;q[2]=0.0;q[3]=0.0;
  Inertia = (double *)emalloc(sizeof(double)*3*3);
  com     = (double *)emalloc(sizeof(double)*3);
  vel     = (double *)emalloc(sizeof(double)*3);
  omg     = (double *)emalloc(sizeof(double)*3);
  acc     = (double *)emalloc(sizeof(double)*3);
  domg    = (double *)emalloc(sizeof(double)*3);
  rotmat  = (double *)emalloc(sizeof(double)*3*3);

  if (argc < 6) {
    printf("USAGE: ./RBD1.exe inputfile(crd) inputfile(parmtop) inputfile(mdin) outputfile(traj) outputfile(ene)\n");
    exit(1);
  }
  infilname1  = *++argv;
  infilname2  = *++argv;
  infilname3  = *++argv;
  outfilname1 = *++argv;
  outfilname2 = *++argv;

  infil2=efopen(infilname2,"r");
  readParmtop(infil2);
  numatom=AP.NATOM;
  mas=(double *)emalloc(sizeof(double)*numatom);
  crd=(double *)emalloc(sizeof(double)*numatom*3);
  vec=(double *)emalloc(sizeof(double)*numatom*3);
  for (i=0;i<numatom;++i)
    mas[i]=AP.AMASS[i];
  infil =efopen(infilname1 ,"r");
  LreadAmberCoord(infil,crd,numatom);
  infil3 =efopen(infilname3,"r");
  fscanf(infil3,"%lf",&dt);
  fscanf(infil3,"%d",&numstep);
  for(i=0;i<3;++i)
    fscanf(infil3,"%lf",&omg[i]);
  fclose(infil);
  fclose(infil2);
  fclose(infil3);

  Mass=RBD_setcom(crd,mas,numatom,com);
  Lcvec(crd,numatom,com,vec);
  RBD_setInertia(crd,com,mas,numatom,Inertia);

  outfil1=efopen(outfilname1,"w");
  outfil2=efopen(outfilname2,"w");

  for (i=0;i<3;++i) {
    for (j=0;j<5;++j) {
      omgp[i][j]=0.0;omgc[i][j]=0.0;
    }
    omgc[i][1]=omg[i]*dt;
  }
  for (j=0;j<3;++j)
    omg_temp[j]=omg[j];
  qua_trans_omgtodqua(omg_temp,q,dq);
  for (i=0;i<4;++i) {
    for (j=0;j<5;++j) {
      qp[i][j]=0.0;qc[i][j]=0.0;
    }
    qc[i][1]=dq[i]*dt;
    qc[i][0]=q[i];
  }

  intg_set();
  for (i=0;i<numstep;++i){
    for (j=0;j<3;++j)
      intg_pc4(1,omgp[j],omgc[j],dt,0.0);
    for (j=0;j<4;++j)
      intg_pc4(1,qp[j],qc[j],dt,0.0);

    for (j=0;j<3;++j) {
      omg[j]=omgp[j][0];
      vel[j]=0.0;
    }
    Lseof(vel,omg,Inertia,Mass,acc,domg);

    for (j=0;j<3;++j)
      intg_pc4(0,omgp[j],omgc[j],dt,domg[j]);
    for (j=0;j<3;++j){
      omg_temp[j]=omgc[j][0]/dt;
      omg[j]=omgc[j][0];
    }
    qua_trans_omgtodqua(omg_temp,q,dq);

    for (j=0;j<4;++j)
      intg_pc4(0,qp[j],qc[j],dt,dq[j]);
    lenq=inprod(q,q,4);
    for (j=0;j<4;++j)
      qc[j][0]=qc[j][0]/lenq;
    for (j=0;j<4;++j)
      q[j]=qc[j][0];

    ene=Lcenergy(vel,omg,Mass,Inertia);

    for (j=0;j<numatom;++j) {
      for (k=0;k<3;++k)
	r[k+1]=vec[j*3+k];
      r[0]=0.0;
      qua_rot(r,q,roted);
      for (k=0;k<3;++k)
	crd[j*3+k]=com[k]+roted[k+1];
    }

    outputconf(outfil1,numatom,crd,'x');
    Lwriteene(outfil2,ene);
  }
  fclose(outfil1);
  fclose(outfil2);

  free(crd);
  free(mas);
  free(vec);

  free(Inertia);
  free(com);
  free(vel);    
  free(omg);    
  free(acc);    
  free(domg);   
  free(rotmat); 

  return 0;
}

void Lseof(double *vel,double *omg,double *Inertia,double Mass,double *acc,double *domg){
  int i;
  double *mat,*invmat;
  double *corif,*a;

  mat    = (double *)ecalloc(6*6,sizeof(double));
  invmat = (double *)emalloc(sizeof(double)*6*6);  
  a      = (double *)ecalloc(6,sizeof(double));
  corif  = (double *)ecalloc(6,sizeof(double));

  LsetmatCOM(Inertia,Mass,mat);
  LsetcoriolisCOM(omg,Inertia,corif);
  invm(mat, invmat, 6);
  for (i=0;i<6;++i)
    corif[i]=-corif[i];
  mvmult(invmat,corif,a,6);
  for (i=0;i<3;++i) {
    domg[i]=a[i];
    acc[i]=a[i+3];
  }

  free(mat);
  free(invmat);
  free(a);
  free(corif);
}

double Lcenergy(double *vel ,double *omg ,double Mass ,double *Inertia) {
  double *omgx,*omgxvec,*rotomgxvec,*Ixomg;
  double ene;

  Ixomg       = (double *)emalloc(sizeof(double)*3*3);

  mvmult(Inertia,omg,Ixomg,3);
  ene=0.5*Mass*inprod(vel,vel,3)+0.5*inprod(omg,Ixomg,3);

  free(Ixomg);

  return ene;
}


void Lcvec(double *crd,int numatom,double com[3], double *vec) {
  int i,j;

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j)
      vec[i*3+j]=crd[i*3+j]-com[j];
}

void Lsetmat(double *Inertia, double *com, double Mass,double *rotmat, double *mat) {
  int i,j;
  double mattemp[6][6];
  double *Vectil,*rotVectil,*rotVectilT,*VectilT;

  Vectil     = (double *)emalloc(sizeof(double)*3*3);
  VectilT    = (double *)emalloc(sizeof(double)*3*3);
  rotVectil  = (double *)emalloc(sizeof(double)*3*3);
  rotVectilT = (double *)emalloc(sizeof(double)*3*3);

  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      mat[i*6+j]=Inertia[i*3+j];
  mat[3*6+3]=Mass;mat[4*6+4]=Mass;mat[5*6+5]=Mass;
  v_product(com,Vectil);
  mmult(rotmat,Vectil,rotVectil,3);
  mtrans(rotVectil,rotVectilT,3);
  mtrans(Vectil,VectilT,3);
  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      mat[(i+3)*6+j]=Mass*rotVectilT[i*3+j];
  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      mat[i*6+(j+3)]=Mass*rotVectil[i*3+j];

  for (i=0;i<6;++i)
    for (j=0;j<6;++j)
      mattemp[i][j]=mat[i*6+j];

  free(Vectil);
  free(rotVectil);
  free(rotVectilT);
}

void LsetmatCOM(double *Inertia,double Mass,double *mat) {
  int i,j;

  for (i = 0; i < 36; ++i)
    mat[i]=0.0;   	

  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      mat[i*6+j]=Inertia[i*3+j];
  mat[3*6+3]=Mass;mat[4*6+4]=Mass;mat[5*6+5]=Mass;

}

void LsetcoriolisCOM(double *omg, double *Inertia, double *corif) {
  int i;
  double *omgtil,*omgInertia,*omgInertiaomg;

  omgtil             = (double *)emalloc(sizeof(double)*3*3);
  omgInertia         = (double *)emalloc(sizeof(double)*3*3);
  omgInertiaomg      = (double *)emalloc(sizeof(double)*3*3);

  v_product(omg,omgtil);
  mmult(omgtil,Inertia,omgInertia,3);
  mvmult(omgInertia,omg,omgInertiaomg,3);

  for (i=0;i<3;++i)
    corif[i]=omgInertiaomg[i];
  for (i=3;i<6;++i)
    corif[i]=0.0;

  free(omgInertiaomg);
  free(omgInertia);
  free(omgtil);
}

int Lwriteene(FILE *out,double ene) {
   fprintf(out,"%12.8lf \n",ene);
   return 1;
 }
 
void LreadAmberCoord(FILE *input, double *crd, int numatom) {
  int i,j;
  int d;
  double x;
  char *line;
  size_t len=0;
  
  getline(&line,&len,input);
  
  fscanf(input,"%d",&d);
  for(i=0;i<numatom;++i) {
    for(j=0;j<3;++j) {
      if (fscanf(input,"%lf",&x) == EOF) {
	printf("error: initial coordinate file eror !!! \n");
	exit(1);
      }
      else {
	crd[i*3+j]=x;
      }
    }
  }
}

int outputconf(FILE *outputfile,int numatom,double *data,int flag){
  int i,j,d;
  double f;
  
  for (i=0;i<numatom;++i) {
    if (flag=='c')
      fprintf(outputfile,"%d ",i);
    for (j=0;j<3;++j)
      fprintf(outputfile,"%10.6lf",data[i*3+j]);
    fprintf(outputfile,"\n");
  }
  fprintf(outputfile,"\n");

  return 0;
}
