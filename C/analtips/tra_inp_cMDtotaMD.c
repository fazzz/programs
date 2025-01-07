#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "f2c.h"
#include "clapack.h"

#include "EF.h"
#include "PT.h"
#include "RBD.h"
#include "MB.h"
#include "LA.h"

int invm2(double *mat, double *invmat, int num);

int main(int argc, char *argv[]) {
  int i,j,k;
  int numatom,numdihed;
  int ininum,ternum,numtermbody;
  int *atom_dihed_pair;
  int flag;

  double dt;
  double *com1,*com2;
  double *Inertia,*InvInertia;
  double *mom;
  double *crd_org;
  double *coord,*coord2,*coordtb,*coordtb2,*velotb,*mass,*velo;
  double *dihed1,*dihed2,*dihed_all,*ddihed;
  double *velcom;
  double *coorddummy,*velotbdummy,*cv;
  double *og;

  double *mat,*Invmat;
  
  char *inputfilename1,*inputfilename2,*inputfilename3,*outputfilename;
  FILE *inputfile1,*inputfile2,*inputfile3,*outputfile;

  char *line;
  size_t len=0;
  
  if (argc < 5) {
    printf("USAGE: tra_cMDtotaMD flag(c or p or o or k) inputfilename1(rst) inputfilename2(cond) inputfilename3(parm) outputfilename(vel.in)\n");
    exit(1);
  }
  flag=(*++argv)[0];
  if (flag != 'c' && flag != 'p' && flag != 'o' && flag != 'k' ) {
    printf("flag error: must be c or p or o or k");
    exit(1);
  }
  inputfilename1 =  *++argv;
  inputfilename2 = *++argv;
  inputfilename3 = *++argv;
  outputfilename = *++argv;

  Inertia=(double *)emalloc(sizeof(double)*3*3);
  InvInertia=(double *)emalloc(sizeof(double)*3*3);

  velcom=(double *)ecalloc(sizeof(double),3);
  og=(double *)ecalloc(sizeof(double),3);
  mom=(double *)ecalloc(sizeof(double),3);

  com1=(double *)ecalloc(sizeof(double),3);
  com2=(double *)ecalloc(sizeof(double),3);

  coorddummy=(double *)ecalloc(sizeof(double),3);
  velotbdummy=(double *)ecalloc(sizeof(double),3);

  cv=(double *)ecalloc(sizeof(double),3);
  
  inputfile2=efopen(inputfilename2,"r");
  fscanf(inputfile2,"%d",&ininum);
  fscanf(inputfile2,"%d",&ternum);
  fscanf(inputfile2,"%lf",&dt);
  if (flag == 'c')
    fscanf(inputfile2,"%d",&numdihed);
  else
    numdihed=AP.NPHIH+AP.MPHIA;
  atom_dihed_pair=(int *)emalloc(sizeof(int)*4*numdihed);
  if (flag == 'c') {
    for (i=0;i<numdihed;++i) {
      fscanf(inputfile2,"%d",&atom_dihed_pair[i*4+0]);
      fscanf(inputfile2,"%d",&atom_dihed_pair[i*4+1]);
      fscanf(inputfile2,"%d",&atom_dihed_pair[i*4+2]);
      fscanf(inputfile2,"%d",&atom_dihed_pair[i*4+3]);
    }
  }
  numtermbody=ternum-ininum+1;
  if (numtermbody<0)
    printf("error: num of terminal is negative");
  coordtb =(double *)ecalloc(sizeof(double),numtermbody*3);
  coordtb2=(double *)ecalloc(sizeof(double),numtermbody*3);
  velotb  =(double *)ecalloc(sizeof(double),numtermbody*3);
  fclose(inputfile2);

  inputfile3=efopen(inputfilename3,"r");
  readParmtop(inputfile3);
  numatom=AP.NATOM;
  coord =(double *)ecalloc(sizeof(double),numatom*3);
  coord2=(double *)ecalloc(sizeof(double),numatom*3);
  velo  =(double *)ecalloc(sizeof(double),numatom*3);
  //atom_dihed_pair=(int *)emalloc(sizeof(int)*4*(AP.NPHIH+AP.MPHIA));
  mass=(double *)ecalloc(sizeof(double),numatom);
  for (i=0;i<numatom;++i)
    mass[i]=AP.AMASS[i];
  fclose(inputfile3);

  inputfile1=efopen(inputfilename1,"r");
  getline(&line,&len,inputfile1);
  getline(&line,&len,inputfile1);
  io_scanconf(inputfile1,numatom,coord ,'x');
  io_scanconf(inputfile1,numatom,velo,'x');
  fclose(inputfile1);
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j)
      coord2[i*3+j]=coord[i*3+j]+dt*20.455*velo[i*3+j];
  if (flag=='k')
    numdihed=pick_atom_dihed_pairs(atom_dihed_pair,'k');
  dihed1=(double *)ecalloc(sizeof(double),numdihed);
  dihed2=(double *)ecalloc(sizeof(double),numdihed);
  ddihed=(double *)ecalloc(sizeof(double),numdihed);
  pick_dihed_all(coord ,dihed1,numdihed,atom_dihed_pair);
  pick_dihed_all(coord2,dihed2,numdihed,atom_dihed_pair);
  for (i=0;i<numdihed;++i)
    ddihed[i]=(dihed2[i]-dihed1[i])/dt;

  /*********************************************/
  /* RBD_setcom(coord, mass,numtermbody,com1); */
  /* RBD_setcom(coord2,mass,numtermbody,com2); */
  /*********************************************/
  //  for (i=0;i<3;++i)
    //    velcom[i]=(com2[i]-com1[i])/dt;
  for (i=0;i<numtermbody;++i) {
    for (j=0;j<3;++j) {
      coordtb[i*3+j] =coord[i*3+j] -com1[j];
      coordtb2[i*3+j]=coord2[i*3+j]-com2[j];
    }
  }
  crd_org=(double *)ecalloc(sizeof(double),3);
  for (i=0;i<3;++i)
   crd_org[i]=0.0;
  outputfile=efopen(outputfilename,"w");
  fclose(outputfile);

  RBD_setInertia(coordtb,crd_org,mass,numtermbody,Inertia);
  for (i=0;i<numtermbody;++i)
    for (j=0;j<3;++j)
      velotb[i*3+j]=(coord2[i*3+j]-coord[i*3+j])/dt;
  for (i=0;i<3;++i)
      mom[i]=0.0;
  for (i=0;i<numtermbody;++i) {
    for (j=0;j<3;++j) {
      coorddummy[j]=coord[i*3+j];
      velotbdummy[j]=velotb[i*3+j];
    }
    outprod(coorddummy,velotbdummy,cv);
    for (j=0;j<3;++j)
      mom[j]+=mass[i]*cv[j];
  }
  invm2(Inertia,InvInertia,3);
  mvmult(InvInertia,mom,og,3);

  /****************************************/
  /*  for (i=0;i<3;++i)			  */
  /*   velotbdummy[i]=velo[i]-velo[3+i];  */
  /* Inertia[0]=0.0;			  */
  /* Inertia[4]=0.0;			  */
  /* Inertia[8]=0.0;			  */
  /* Inertia[1]=-(coord[3+2]-coord[2]);	  */
  /* Inertia[2]= (coord[3+1]-coord[1]);	  */
  /* Inertia[3]= (coord[3+2]-coord[2]);	  */
  /* Inertia[5]=-(coord[3+0]-coord[0]);	  */
  /* Inertia[6]=-(coord[3+1]-coord[1]);	  */
  /* Inertia[7]= (coord[3+0]-coord[0]);	  */
  /* 					  */
  /* invm2(Inertia,InvInertia,3);	  */
  /* mvmult(InvInertia,velotbdummy,og,3); */
  /****************************************/

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<3;++i)
    fprintf(outputfile,"%lf\n",og[i]);
  for (i=0;i<3;++i)
    fprintf(outputfile,"%lf\n",/*velcom[i]*/20.455*velo[i]);
  for (i=0;i<numdihed;++i)
    fprintf(outputfile,"%lf\n",ddihed[i]);
  fclose(outputfile);

  /*********************/
  /* free(Inertia);    */
  /* free(InvInertia); */
  /* free(coord);      */
  /* free(coord2);     */
  /* free(coordtb);    */
  /* free(coordtb2);   */
  /* free(velotb);     */
  /* free(mass);       */
  /* free(dihed1);     */
  /* free(dihed2);     */
  /* free(dihed_all);  */
  /* free(velcom);     */
  /* free(og);	       */
  /* free(mom);	       */
  /*********************/
    
  return 0;
}

int invm2(double *mat, double *invmat, const int num) {

  int i,j,k;
  double *mattemp;
  static long int m1,n1,lda,info,piv[500],lwork=500;
  static double work[500];

  m1 = num;
  n1 = num;
  lda=m1;
  
  mattemp=ecalloc(sizeof(double),m1*n1);
  memcpy(mattemp,mat,sizeof(double)*m1*n1);
  mtrans(mat,mattemp,m1);
  dgetrf_(&m1,&n1,mattemp,&lda,piv,&info);
  if (info!=0) return 0;
  dgetri_(&n1,mattemp,&lda,piv,work,&lwork,&info);
  if (info!=0) return 0;
  mtrans(mattemp,invmat,m1);
  memcpy(mattemp,mat,sizeof(double)*m1*n1);
  free(mattemp);
  

  return 1;
}


