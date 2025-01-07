#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "GOLM_Clementi_set.h"
#include "GOLM_Clementi.h"
#include "GOLM_Clementi_check.h"

#include "PTL.h"
#include "EF.h"
#include "NC.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d;
  int numatom,numCAatom,numres;

  double f_natatt[3],f_repul[3],f_d[3],f_a[3],f_b[3];
  double dx=0.00001;
  int numspatom=10,nums=3,numa=4;

  double f_d1[4][3],f_d2[4][3];
  double f_natatt1[2][3],f_natatt2[2][3];

  int **ncmap;
  double dt=0.001;
  double *crd,*refcrd,*refcrdAA,*mass;
  struct potential_GOLM_Clementi e;

  double x[3];

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*reffilename,*outputfilename,*parmfilename,*mapfilename;
  FILE *inputfile,*reffile,*outputfile,*parmfile,*mapfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"dx",1,NULL,'x'},
    {"nums",1,NULL,'s'},
    {"numa",1,NULL,'a'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hx:s:a:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'x':
      dx=atof(optarg);
      break;
    case 's':
      numspatom=atoi(optarg);
      break;
    case 'a':
      numa=atoi(optarg);
      break;
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (argc < 5) {
    USAGE(progname);
    exit(1);
  }
  inputfilename     = *argv;
  reffilename       = *++argv;
  parmfilename      = *++argv;
  outputfilename    = *++argv;
  mapfilename       = *++argv;

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
  numres=AP.NRES;
  mass=(double *)gcemalloc(sizeof(double)*numCAatom);
  j=0;
  for (i=0;i<numatom;++i) {
    if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
      mass[j]=AP.AMASS[i];
      ++j;
    }
  }

  crd=(double *)gcemalloc(sizeof(double)*numCAatom*3);
  refcrd=(double *)gcemalloc(sizeof(double)*numCAatom*3);
  refcrdAA=(double *)gcemalloc(sizeof(double)*numatom*3);

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

  GOLM_Clementi_ff_set_calcff(&e,refcrd,refcrdAA,numCAatom,numatom);

  GOLM_Clementi_ff_calcff(crd,numCAatom,&e);
  GOLM_Clementi_calcff_check(inputfilename,reffilename,parmfilename,numspatom,dx,f_natatt,f_repul,f_d,f_a,f_b,f_d1,f_d2,nums,f_natatt1,f_natatt2,numa);

  printf("p_t=%e \n",e.p_t);
  printf("p_nat=%e \n",e.p_natatt_t);
  printf("p_rep=%e \n",e.p_repul_t);
  printf("p_dih=%e \n",e.p_d_t);
  printf("p_ang=%e \n",e.p_a_t);
  printf("p_bon=%e \n",e.p_b_t);

  printf("f_natt =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",e.f_natatt[numspatom][i]);
  }
  printf(")\n");

  printf("f_natt =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f_natatt[i]);
  }
  printf(")\n");

  printf("f_repl =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",e.f_repul[numspatom][i]);
  }
  printf(")\n");

  printf("f_repl =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f_repul[i]);
  }
  printf(")\n");

  printf("f_d =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",e.f_d[numspatom][i]);
  }
  printf(")\n");

  printf("f_d =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f_d[i]);
  }
  printf(")\n");

  printf("f_a =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",e.f_a[numspatom][i]);
  }
  printf(")\n");

  printf("f_a =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f_a[i]);
  }
  printf(")\n");

  printf("f_b =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",e.f_b[numspatom][i]);
  }
  printf(")\n");

  printf("f_b =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f_b[i]);
  }
  printf(")\n");

  /**********************************/
  /* for (i=0;i<4;++i) {	    */
  /*   printf("f_d(%d) =( ",i+1);   */
  /*   for (j=0;j<3;++j) {	    */
  /*     printf("%e ,",f_d1[i][j]); */
  /*   }			    */
  /*   printf(")\n");		    */
  /*   printf("f_d(%d) =( ",i+1);   */
  /*   for (j=0;j<3;++j) {	    */
  /*     printf("%e ,",f_d2[i][j]); */
  /*   }			    */
  /*   printf(")\n");		    */
  /* }				    */
  /**********************************/

  /***************************************/
  /* for (i=0;i<2;++i) {		 */
  /*   printf("f_n(%d) =( ",i+1);	 */
  /*   for (j=0;j<3;++j) {		 */
  /*     printf("%e ,",f_natatt1[i][j]); */
  /*   }				 */
  /*   printf(")\n");			 */
  /*   printf("fn(%d) =( ",i+1);	 */
  /*   for (j=0;j<3;++j) {		 */
  /*     printf("%e ,",f_natatt2[i][j]); */
  /*   }				 */
  /*   printf(")\n");			 */
  /* }					 */
  /***************************************/

  ncmap=GOLM_Clementi_make_native_contact(refcrdAA,6.5,&((e).num_natatt),numatom,numCAatom);
  mapfile=efopen(mapfilename,"w");
  for (i=0;i<numCAatom;++i) {
    for (j=0;j<numCAatom;++j) {
      fprintf(mapfile,"%2d ",ncmap[i][j]);
    }
  }
  fclose(mapfile);


  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename clustfilename parmfilename outputfilename\n",progname);
}


