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
#include "INTG.h"
#include "QUA.h"
#include "LA.h"

void Lseof(double *vel,double *omg,double *Inertia,double Mass,double *Vec,double *rotmat,double *acc,double *domg,int flag );
void Lcvec(double *crd,int numatom,double *vec);
void Lsetcoriolis(double *omg, double *vel,double *Inertia, double *Vec, double Mass, double *rotmat, double *corif);
void Lsetmat(double *Inertia, double *Vec, double Mass,double *rotmat, double *mat);
void LsetcoriolisCOM(double *omg, double *Inertia, double *Vec, double Mass, double *rotmat, double *corif);
void LsetmatCOM(double *Inertia, double *Vec, double Mass,double *rotmat, double *mat);
double Lcenergy(double *vel ,double *omg ,double *Inertia  ,double Mass,double *Vec ,double *rotmat, int flag);
int Lwriteene(FILE *outfil,double ene);
int Lwritecrd(FILE *outfil,double *crd,int numatom);
int Lreadmdin(FILE *mdin, double dt, int numstep);
void LreadAmberCoord(FILE *input, double *crd, int numatom);
void LsetIni(double *omg,double *vel, double *Vec, double q[4],double qc[4][5],double corv[6][6], double *rotmat,double dt, int flag);

int main(int argc, char *argv[]) {
  int i,j,k,numstep,numatom,flag,n;

  double dt;
  double *crd,*vec,*crd_org;
  double *Inertia,Mass,*mas;
  double ene;
  double *Vec;
  double *vel,*omg,*acc,*domg,omg_temp[3];
  double *rotmat,rotmattemp[3][3];
  double prev[6][6],corv[6][6],qp[4][5],qc[4][5],op[3][6],oc[3][6],vp[3][6],vc[3][6];
  double q[4],r[4],dq[4],roted[4],lenq;

  char *infilname1,*infilname2,*infilname3,*outfilname1,*outfilname2;
  FILE *infil,*infil2,*infil3,*outfil1,*outfil2;

  int flag2;

  crd_org = (double *)emalloc(sizeof(double)*3);
  Inertia = (double *)emalloc(sizeof(double)*3*3);
  Vec     = (double *)emalloc(sizeof(double)*3);
  vel     = (double *)emalloc(sizeof(double)*3);
  omg     = (double *)emalloc(sizeof(double)*3);
  acc     = (double *)emalloc(sizeof(double)*3);
  domg    = (double *)emalloc(sizeof(double)*3);
  rotmat  = (double *)emalloc(sizeof(double)*3*3);

  if (argc < 6) {
    printf("USAGE: ./RBD1.exe flag(c/n) inputfile(crd) inputfile(parmtop) inputfile(mdin) outputfile(traj) outputfile(ene)\n");
    exit(1);
  }
  flag=(*++argv)[0];
  if (flag == 'c')
    flag=COM;
  else if (flag == 'n')
    flag=NOTCOM;
  else if (flag != 'c' && flag != 'n') {
    printf("flag error: must be c or n");
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
  for (i=0;i<3;++i)
    crd_org[i]=crd[i];
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j)
      crd[i*3+j]=crd[i*3+j]-crd_org[j];
  infil3 =efopen(infilname3,"r");
  fscanf(infil3,"%lf",&dt);
  fscanf(infil3,"%d",&numstep);
  fclose(infil);
  fclose(infil2);
  fclose(infil3);

  Lcvec(crd,numatom,vec);
  Mass=setcom(crd,mas,numatom,Vec);
  setInertia(crd,crd_org,mas,numatom,Inertia,flag);

  LsetIni(omg,vel,Vec,q,qc,corv,rotmat,dt,flag);
  for (i=0;i<3;++i) {
    oc[i][1]=dt*omg[i];
    vc[i][1]=dt*vel[i];
  }

  intg_set();

  outfil1=efopen(outfilname1,"w");
  outfil2=efopen(outfilname2,"w");

  for (i=0;i<numstep;++i){
    ene=Lcenergy(vel,omg,Inertia,Mass,Vec,rotmat,flag);

    qua_set_rotmat(q,rotmattemp);
    for (j=0;j<3;++j)
      for (k=0;k<3;++k)
	rotmat[j*3+k]=rotmattemp[j][k];
    Lseof(vel,omg,Inertia,Mass,Vec,rotmat,acc,domg,flag);
    for (j=0;j<3;++j)
      omg_temp[j]=omg[j];
    qua_trans_omgtodqua(omg_temp,q,dq);

    for (j=0;j<3;++j) {
      omg[j]+=dt*domg[j];
      vel[j]+=dt*acc[j];
    }
    
    for (j=0;j<4;++j)
      q[j]+=dt*dq[j];

    for (j=0;j<3;++j)
      crd[j]=corv[j][0];
    for (j=0;j<numatom;++j) {
      for (k=0;k<3;++k)
	r[k]=vec[j*3+k];
      r[3]=0.0;
      qua_rot(q,r,roted);
      for (k=0;k<3;++k)
	crd[j*3+k]=crd[k]+roted[k];
    }
    Lwriteene(outfil2,ene);
  }
  fclose(outfil2);

  free(crd);
  free(vec);
  free(mas);

  free(crd_org);
  free(Inertia);
  free(Vec);
  free(vel);    
  free(omg);    
  free(acc);    
  free(domg);   
  free(rotmat); 

  return 0;
}

void Lseof(double *vel,double *omg,double *Inertia,double Mass,double *Vec,double *rotmat,double *acc,double *domg,int flag ){
  int i;
  double *mat,*invmat;
  double *corif,*mcorif,*a;
  double *omgx,*omgxomg,*omgxomgxVec,*rotmatomgxomgxVec,*omgxomgxdomg,*rotomgxomgxdomg;

  mat    = (double *)ecalloc(6*6,sizeof(double));
  invmat = (double *)emalloc(sizeof(double)*6*6);  
  a      = (double *)ecalloc(6,sizeof(double));
  corif  = (double *)ecalloc(6,sizeof(double));
  mcorif = (double *)ecalloc(6,sizeof(double));

  omgx   = (double *)emalloc(sizeof(double)*3*3);
  omgxomg= (double *)emalloc(sizeof(double)*3*3);
  omgxomgxVec=(double *)emalloc(sizeof(double)*3);
  omgxomgxdomg=(double *)emalloc(sizeof(double)*3);
  rotmatomgxomgxVec=(double *)emalloc(sizeof(double)*3);
  rotomgxomgxdomg=(double *)ecalloc(3,sizeof(double)*3);

  if (flag==COM) {
    LsetmatCOM(Inertia,Vec,Mass,rotmat,mat);
    LsetcoriolisCOM(omg,Inertia,Vec,Mass,rotmat,corif);
  }
  else {
    Lsetmat(Inertia,Vec,Mass,rotmat,mat);
    Lsetcoriolis(omg,vel,Inertia,Vec,Mass,rotmat,corif);
  }
  invm(mat, invmat, 6);
  for (i=0;i<6;++i)
    mcorif[i]=-corif[i];
  mvmult(invmat,mcorif,a,6);
  for (i=0;i<3;++i) {
    domg[i]=a[i];
    acc[i]=a[i+3];
  }

  /*********************************************************************/
  /* v_product(omg,omgx);					       */
  /* mmult(omgx,omgx,omgxomg,3);				       */
  /* mvmult(omgxomg,Vec,omgxomgxVec,3);				       */
  /* mvmult(rotmat,omgxomgxVec,rotmatomgxomgxVec,3);		       */
  /* 								       */
  /* mvmult(omgxomg,domg,omgxomgxdomg,3);			       */
  /* mvmult(rotmat,omgxomgxdomg,rotomgxomgxdomg,3);		       */
  /* 								       */
  /* /\******************************************************\/	       */
  /* /\* for (i=0;i<3;++i) {			        *\/	       */
  /* /\*   acc[i]=-rotmatomgxomgxVec[i]+rotomgxomgxdomg[i]; *\/	       */
  /* /\* }						        *\/    */
  /* /\******************************************************\/	       */
  /*********************************************************************/


  free(invmat);
  free(corif);
  free(mcorif);
  free(a);
  free(mat);

  free(omgx);
  free(omgxomg);
  free(omgxomgxVec);
  free(omgxomgxdomg);
  free(rotmatomgxomgxVec);
  free(rotomgxomgxdomg);


}

double Lcenergy(double *vel ,double *omg ,double *Inertia ,double Mass,double *Vec ,double *rotmat, int flag) {
  double *omgx,*omgxvec,*rotomgxvec,*Ixomg;
  double ene;

  omgx        = (double *)emalloc(sizeof(double)*3*3);
  omgxvec     = (double *)emalloc(sizeof(double)*3);
  rotomgxvec  = (double *)emalloc(sizeof(double)*3);
  Ixomg       = (double *)emalloc(sizeof(double)*3*3);

  if (flag == NOTCOM) {
    v_product(omg,omgx);
    mvmult(omgx,Vec,omgxvec,3);
    mvmult(rotmat,omgxvec,rotomgxvec,3);
  }
  mvmult(Inertia,omg,Ixomg,3);

   if (flag == NOTCOM) 
    ene=0.5*Mass*inprod(vel,vel,3)+Mass*inprod(vel,rotomgxvec,3)+0.5*inprod(omg,Ixomg,3);
  else
    ene=0.5*Mass*inprod(vel,vel,3)+0.5*inprod(omg,Ixomg,3);

  free(omgx);
  free(omgxvec);
  free(rotomgxvec);
  free(Ixomg);

  return ene;
}

int Lwritecrd(FILE *out,double *crd,int numatom) {
  int i,j;

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      fprintf(out,"%12.8lf ",crd[i*numatom+j]);
    }
    fprintf(out,"\n");
  }
  fprintf(out,"\n");

  return 1;
}

int Lwriteene(FILE *out,double ene) {
  fprintf(out,"%12.8lf \n",ene);
  return 1;
}

int Lreadmdin(FILE *mdin, double dt, int numstep) {
  fscanf(mdin,"%lf",&dt);
  fscanf(mdin,"%d",&numstep);
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

void Lcvec(double *crd,int numatom,double *vec) {
  int i,j;

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j)
      vec[i*3+j]=crd[i*3+j]-crd[j];
}

void Lsetmat(double *Inertia, double *Vec, double Mass,double *rotmat, double *mat) {
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
  v_product(Vec,Vectil);
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

void Lsetcoriolis(double *omg, double *vel,double *Inertia, double *Vec, double Mass, double *rotmat, double *corif) {
  int i;
  double *omgtil,*omgInertia,*omgInertiaomg,*omgxomg,*omgxomgV,*rotomgxomgV;
  double *Vectil,*omgxVecx,*omgxVecxvel;

  omgtil             = (double *)emalloc(sizeof(double)*3*3);
  omgInertia         = (double *)emalloc(sizeof(double)*3*3);
  omgInertiaomg      = (double *)emalloc(sizeof(double)*3*3);
  omgxomg            = (double *)emalloc(sizeof(double)*3*3);
  omgxomgV           = (double *)emalloc(sizeof(double)*3);  
  rotomgxomgV        = (double *)emalloc(sizeof(double)*3);  

  v_product(omg,omgtil);
  mmult(omgtil,Inertia,omgInertia,3);
  mvmult(omgInertia,omg,omgInertiaomg,3);

  mmult(omgtil,omgtil,omgxomg,3);
  mvmult(omgxomg,Vec,omgxomgV,3);
  mvmult(rotmat,omgxomgV,rotomgxomgV,3);

  for (i=0;i<3;++i)
    corif[i]=omgInertiaomg[i];
  for (i=3;i<6;++i)
    corif[i]=Mass*rotomgxomgV[i-3];

  free(omgInertiaomg);
  free(omgInertia);
  free(omgtil);
  free(omgxomg);
  free(omgxomgV);
  free(rotomgxomgV);
}

void LsetmatCOM(double *Inertia, double *Vec, double Mass,double *rotmat, double *mat) {
  int i,j;

  for (i = 0; i < 36; ++i) {
    mat[i]=0.0;   	

  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      mat[i*6+j]=Inertia[i*3+j];
  mat[3*6+3]=Mass;mat[4*6+4]=Mass;mat[5*6+5]=Mass;

}

void LsetcoriolisCOM(double *omg, double *Inertia, double *Vec, double Mass, double *rotmat, double *corif) {
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

  /******************************************************************/
  /* for (i=0;i<numstep;++i){					    */
  /*   for (j=0;j<3;++j){					    */
  /*     Gear(1,prev[j],corv[j],dt,0.0);			    */
  /*     omg[j]=prev[j][1]/dt;					    */
  /*   }							    */
  /*   for (j=3;j<6;++j){					    */
  /*     Gear(1,prev[j],corv[j],dt,0.0);			    */
  /*     vel[j-3]=prev[j][1]/dt;				    */
  /*   }							    */
  /*   Lseof(vel,omg,Inertia,Mass,Vec,rotmat,acc,domg,flag);	    */
  /*   for (j=0;j<3;++j) {					    */
  /*     Gear(2,prev[j],corv[j],dt,domg[j]);			    */
  /*     omg[j]=corv[j][1]/dt;					    */
  /*   }							    */
  /*   for (j=3;j<6;++j){					    */
  /*     Gear(2,prev[j],corv[j],dt,acc[j-3]);			    */
  /*     vel[j-3]=corv[j][1]/dt;				    */
  /*   }							    */
  /*   for (j=0;j<3;++j)					    */
  /*     omg_temp[j]=omg[j];					    */
  /*   qua_trans_omgtodqua(omg_temp,q,dq);			    */
  /*   for (j=0;j<4;++j)					    */
  /*     q[j]+=dq[j]*dt;					    */
  /*   lenq=inprod(q,q,4);					    */
  /*   for (j=0;j<4;++j)					    */
  /*     q[j]=q[j]/lenq;					    */
  /*   qua_set_rotmat(q,rotmattemp);				    */
  /*   for (j=0;j<3;++j)					    */
  /*     for (k=0;k<3;++k)					    */
  /* 	rotmat[j*3+k]=rotmattemp[j][k];				    */
  /*   for (j=0;j<3;++j)					    */
  /*     crd[j]=corv[j][0];					    */
  /*   for (j=0;j<numatom;++j) {				    */
  /*     for (k=0;k<3;++k)					    */
  /* 	r[k]=vec[j*3+k];					    */
  /*     r[3]=0.0;						    */
  /*     qua_rot(q,r,roted);					    */
  /*     for (k=0;k<3;++k)					    */
  /* 	crd[j*3+k]=crd[k]+roted[k];				    */
  /*   }							    */
  /*   ene=Lcenergy(vel,omg,Inertia,Mass,Vec,rotmat,flag);	    */
  /*   Lwritecrd(outfil1,crd,numatom);				    */
  /*   Lwriteene(outfil2,ene);					    */
  /* }								    */
  /******************************************************************/

void LsetIni(double *omg,double *vel, double *Vec, double q[4],double qc[4][5],double corv[6][6], double *rotmat, double dt,int flag){
  int i,j;
  double *omgx,*omgxVec;

  omgx    = (double *)emalloc(sizeof(double)*3*3);
  omgxVec = (double *)emalloc(sizeof(double)*3);

  msetIni(rotmat,3);
  for (i=0;i<6;++i)
    for (j=0;j<6;++j)
      corv[i][j]=0.0;
  for (i=0;i<4;++i)
    for (j=0;j<5;++j)
      qc[i][j]=0.0;
  qc[0][0]=1.0;qc[1][0]=0.0;qc[2][0]=0.0;qc[3][0]=0.0;
  for (i=0;i<4;++i)
    q[i]=qc[i][0];
  for (i=0;i<3;++i) {
    omg[i]=0.5;
  }
  v_product(omg,omgx);
  mvmult(omgx,Vec,omgxVec,3);
  if (flag==NOTCOM) {
    for (i=0;i<3;++i) {
      vel[i]=-omgxVec[i];
    }
  }

  free(omgx);
  free(omgxVec);

}
