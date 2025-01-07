
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "FF.h"
#include "PT.h"
#include "EF.h"
#include "IO.h"

#define ON 1
#define OFF 0

double calc_FF_sub(char *parmtopname1,char *inputfilename,char *indexfilename,int numstep,int numnb,int num14,int flag,int amberflag,double *p_t);

int main(int argc, char *argv[]) {
  int flag;
  int i,j,k,amberflag;
  double *p_t1,*p_t2;
  int numstep,numnb,num14;

  char *inputfilename,*parmtopname1,*parmtopname2,*indexfilename;
  char *outputfilename;
  FILE *inputfile;
  FILE *outputfile;

  if (argc < 5) {
    printf("USAGE: ./%s [at] [54db]  numnb num14 numstep inputfilename(crd) parmtop1 parmtop2 indexfilename(index) outputfilename(w) \n",argv[0]);
    exit(1);
  }
  flag=(*++argv)[0];
  if (flag != 'a' && flag != 't') {
    printf("flag error: must be a or 4 or t");
    exit(1);
  }  
  if (flag=='a')
    amberflag=ON;
  else
    amberflag=OFF;
  flag=(*++argv)[0];
  if (flag != '5' && flag != '4' && flag != 'd' && flag != 'b') {
    printf("flag error: must be 5 or 4 or d or b");
    exit(1);
  }
  numnb=atoi(*++argv);
  num14=atoi(*++argv);
  numstep=atoi(*++argv);
  inputfilename = *++argv;
  parmtopname1 = *++argv;
  parmtopname2 = *++argv;
  indexfilename = *++argv;
  outputfilename = *++argv;

  p_t1=(double *)gcemalloc(sizeof(double)*numstep);
  p_t2=(double *)gcemalloc(sizeof(double)*numstep);

  calc_FF_sub(parmtopname1,inputfilename,indexfilename,numstep,numnb,num14,flag,amberflag,p_t1);
  calc_FF_sub(parmtopname2,inputfilename,indexfilename,numstep,numnb,num14,flag,amberflag,p_t2);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i)
    fprintf(outputfile,"%lf\n",p_t1[i]-p_t2[i]);
  fclose(outputfile);

  return 0;
}

double calc_FF_sub(char *parmtopname1,char *inputfilename,char *indexfilename,int numstep,int numnb,int num14,int flag,int amberflag,double *p_t){
  int i,j,k,f;
  double *ele,*ALJ,*BLJ;
  double *p_e,*p_LJ,*p_d,*p_a,*p_b;
  double p_e_t,p_LJ_t,p_e_14_t,p_LJ_14_t,p_d_t,p_a_t,p_b_t;
  double *f_e,*f_LJ,*n_d;
  double *p_e_14,*p_LJ_14;
  double *f_e_14,*f_LJ_14;
  double *crd;
  char *line,dummy;
  size_t len=0;

  int *indexnb,*index14;
  int numatom,numpara;
  FILE *inputfile,*parmtop1,*parmtop2,*indexfile;

  parmtop1=efopen(parmtopname1,"r");
  readParmtop(parmtop1);
  fclose(parmtop1);
  numatom=AP.NATOM;
  numpara=AP.NTYPES*(AP.NTYPES+1)/2;

  indexnb=(int *)emalloc(sizeof(int)*numnb*2);
  ele=(double *)emalloc(sizeof(double)*numatom);
  ALJ=(double *)emalloc(sizeof(double)*numatom*numatom);
  BLJ=(double *)emalloc(sizeof(double)*numatom*numatom);
  index14=(int *)emalloc(sizeof(int)*num14*2);

  ff_set_NB_PARM(ele,ALJ,BLJ,numatom);

  indexfile=efopen(indexfilename,"r");
  for (i=0;i<numnb;++i) {
    fscanf(indexfile,"%d",&f);indexnb[i*2]=f-1;
    fscanf(indexfile,"%s",&dummy);
    fscanf(indexfile,"%d",&f);indexnb[i*2+1]=f-1;
  }
  for (i=0;i<num14;++i) {
    fscanf(indexfile,"%d",&f);index14[i*2]=f-1;
    fscanf(indexfile,"%s",&dummy);
    fscanf(indexfile,"%d",&f);index14[i*2+1]=f-1;
  }
  fclose(indexfile);

  inputfile=efopen(inputfilename,"r");

  if (amberflag==ON)
    getline(&line,&len,inputfile);
  for (i=0;i<numstep;++i) {
    crd=(double *)ecalloc(sizeof(double),numatom*3);
    p_e=(double *)ecalloc(sizeof(double),numnb);
    p_LJ=(double *)ecalloc(sizeof(double),numnb);

    p_e_14=(double *)ecalloc(sizeof(double),numnb);
    p_LJ_14=(double *)ecalloc(sizeof(double),numnb);
    p_d=(double *)ecalloc(sizeof(double),(AP.NPHIH+AP.MPHIA));
    p_a=(double *)ecalloc(sizeof(double),(AP.NTHETH+AP.MTHETA));
    p_b=(double *)ecalloc(sizeof(double),(AP.NBONH+AP.MBONA));

    io_scanconf(inputfile,numatom,crd,'x');
    ff_calcFFNB(ele,ALJ,BLJ,p_e,p_LJ,f_e,f_LJ,numnb,indexnb,numatom,crd,2,0);
    ff_calcFFNB(ele,ALJ,BLJ,p_e_14,p_LJ_14,f_e_14,f_LJ_14,num14,index14,numatom,crd,2,0);
    ff_calcDIHE(p_d,n_d,crd,1,0,0);
    ff_calcANGLE(p_a,crd);
    ff_calcBOND(p_b,crd);

    p_t[i]=0.0;
    
    for (j=0;j<numnb;++j)
      p_t[i]+=p_e[j]+p_LJ[j];
    if (flag=='4' || flag=='d' || flag=='b')
      for (j=0;j<num14;++j)
	p_t[i]+=1.0/1.2*p_e_14[j]+0.5*p_LJ_14[j];
    if (flag=='d' || flag=='b')
      for (j=0;j<AP.NPHIH+AP.MPHIA;++j)
	p_t[i]+=p_d[j];
    if (flag=='b') {
      for (j=0;j<AP.NTHETH+AP.MTHETA;++j) {
	p_t[i]+=p_a[j];
	p_t[i]+=p_b[j];
      }
    }

    free(crd);
    free(p_e);
    free(p_LJ);
    free(p_e_14);
    free(p_LJ_14);
    free(p_d);
    free(p_a);
    free(p_b);
  }

  fclose(inputfile);

  indexnb;
  free(ele);
  free(ALJ);
  free(BLJ);
  free(index14);

  return 1.0;
}
