
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "FFL.h"
#include "UMBP.h"
#include "PTL.h"

#define UNIT 4.184070*100.0

double ffL_calcffandforce_check(char *inputfilename,char *parmfilename,
				int numspatom,double dx,
				double f_es[3],double f_LJ[3],
				double f_14_es[3],double f_14_LJ[3],
				double f_d[3], double f_a[3], double f_b[3]
				) {
  int i,j,k,d;
  double *crd,*crddx,*crddy,*crddz,*refcrd;
  int numatom,numres;
  struct potential e,e_dx,e_dy,e_dz;
  struct force f,f_dx,f_dy,f_dz;

  double x[3];

  char *line;
  size_t len=0;

  FILE *inputfile,*outputfile,*parmfile;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  numres=AP.NRES;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

  crddx=(double *)gcemalloc(sizeof(double)*numatom*3);
  crddy=(double *)gcemalloc(sizeof(double)*numatom*3);
  crddz=(double *)gcemalloc(sizeof(double)*numatom*3);

  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfile,"%lf",&crd[i*3+j]);
  fclose(inputfile);

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      crddx[i*3+j]=crd[i*3+j];
      crddy[i*3+j]=crd[i*3+j];
      crddz[i*3+j]=crd[i*3+j];
    }
  }

  crddx[numspatom*3]+=dx;
  crddy[numspatom*3+1]+=dx;
  crddz[numspatom*3+2]+=dx;

  ffL_set_calcffandforce(&e,&f);
  ffL_set_calcffandforce(&e_dx,&f_dx);
  ffL_set_calcffandforce(&e_dy,&f_dy);
  ffL_set_calcffandforce(&e_dz,&f_dz);

  ffL_calcffandforce(crd,numatom,&e,&f);
  ffL_calcffandforce(crddx,numatom,&e_dx,&f_dx);
  ffL_calcffandforce(crddy,numatom,&e_dy,&f_dy);
  ffL_calcffandforce(crddz,numatom,&e_dz,&f_dz);

  f_es[0]=-0.5*(e_dx.p_e_t-e.p_e_t)/dx*4.184070*100.0;
  f_es[1]=-0.5*(e_dy.p_e_t-e.p_e_t)/dx*4.184070*100.0;
  f_es[2]=-0.5*(e_dz.p_e_t-e.p_e_t)/dx*4.184070*100.0;

  f_LJ[0]=-0.5*(e_dx.p_LJ_t-e.p_LJ_t)/dx*4.184070*100.0;
  f_LJ[1]=-0.5*(e_dy.p_LJ_t-e.p_LJ_t)/dx*4.184070*100.0;
  f_LJ[2]=-0.5*(e_dz.p_LJ_t-e.p_LJ_t)/dx*4.184070*100.0;

  f_14_es[0]=-0.5*(e_dx.p_e_14_t-e.p_e_14_t)/dx*4.184070*100.0;
  f_14_es[1]=-0.5*(e_dy.p_e_14_t-e.p_e_14_t)/dx*4.184070*100.0;
  f_14_es[2]=-0.5*(e_dz.p_e_14_t-e.p_e_14_t)/dx*4.184070*100.0;

  f_14_LJ[0]=-0.5*(e_dx.p_LJ_14_t-e.p_LJ_14_t)/dx*4.184070*100.0;
  f_14_LJ[1]=-0.5*(e_dy.p_LJ_14_t-e.p_LJ_14_t)/dx*4.184070*100.0;
  f_14_LJ[2]=-0.5*(e_dz.p_LJ_14_t-e.p_LJ_14_t)/dx*4.184070*100.0;

  f_d[0]=-(e_dx.p_d_t-e.p_d_t)/dx*4.184070*100.0;
  f_d[1]=-(e_dy.p_d_t-e.p_d_t)/dx*4.184070*100.0;
  f_d[2]=-(e_dz.p_d_t-e.p_d_t)/dx*4.184070*100.0;

  f_a[0]=-(e_dx.p_a_t-e.p_a_t)/dx*4.184070*100.0;
  f_a[1]=-(e_dy.p_a_t-e.p_a_t)/dx*4.184070*100.0;
  f_a[2]=-(e_dz.p_a_t-e.p_a_t)/dx*4.184070*100.0;

  f_b[0]=-(e_dx.p_b_t-e.p_b_t)/dx*4.184070*100.0;
  f_b[1]=-(e_dy.p_b_t-e.p_b_t)/dx*4.184070*100.0;
  f_b[2]=-(e_dz.p_b_t-e.p_b_t)/dx*4.184070*100.0;

  //////////////////////////////////////////////////////////////////////////////////////////////

  return 0.0;
}

double UMB_calcffandforce_check(double *crd,int numatom,int *pairp,int nump,double *fcp,double *dih_equ,int numspatom, double dx, double f_UMB[3]){
  int i,j,k,d;
  double *crddx,*crddy,*crddz;

  double *p,*p_dx,*p_dy,*p_dz;
  double p_t,p_t_dx,p_t_dy,p_t_dz;
  double **f,**f_dx,**f_dy,**f_dz;

  p=(int *)gcemalloc(sizeof(int)*nump);
  p_dx=(int *)gcemalloc(sizeof(int)*nump);
  p_dy=(int *)gcemalloc(sizeof(int)*nump);
  p_dz=(int *)gcemalloc(sizeof(int)*nump);

  f=(double **)gcemalloc(sizeof(double *)*numatom);
  f_dx=(double **)gcemalloc(sizeof(double *)*numatom);
  f_dy=(double **)gcemalloc(sizeof(double *)*numatom);
  f_dz=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) {
    f[i]=(double *)gcemalloc(sizeof(double)*3);
    f_dx[i]=(double *)gcemalloc(sizeof(double)*3);
    f_dy[i]=(double *)gcemalloc(sizeof(double)*3);
    f_dz[i]=(double *)gcemalloc(sizeof(double)*3);
  }

  crddx=(double *)gcemalloc(sizeof(double)*numatom*3);
  crddy=(double *)gcemalloc(sizeof(double)*numatom*3);
  crddz=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      crddx[i*3+j]=crd[i*3+j];
      crddy[i*3+j]=crd[i*3+j];
      crddz[i*3+j]=crd[i*3+j];
    }
  }

  crddx[numspatom*3]+=dx;
  crddy[numspatom*3+1]+=dx;
  crddz[numspatom*3+2]+=dx;

  p_t=UMB_calc_dihetype_ff(crd,numatom,pairp,nump,fcp,dih_equ,p,f);
  p_t_dx=UMB_calc_dihetype_ff(crddx,numatom,pairp,nump,fcp,dih_equ,p_dx,f_dx);
  p_t_dy=UMB_calc_dihetype_ff(crddy,numatom,pairp,nump,fcp,dih_equ,p_dy,f_dy);
  p_t_dz=UMB_calc_dihetype_ff(crddz,numatom,pairp,nump,fcp,dih_equ,p_dz,f_dz);

  f_UMB[0]=-(p_t_dx-p_t)/dx*4.184070*100.0;
  f_UMB[1]=-(p_t_dy-p_t)/dx*4.184070*100.0;
  f_UMB[2]=-(p_t_dz-p_t)/dx*4.184070*100.0;

}
