
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "netcdf_mine.h"
#include "EF.h"
#include "IO.h"

#include "MBAR.h"

#define ON 0
#define OFF 1

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,t;
  int interval=1;
  int d1=0,d2=0,num;
  double f1,f2;
  double *pmf,pmf_min,pmf_max,index_max=0;

  int macminflag=OFF;
  int udflag=OFF;

  double din;

  int n_sim;
  int *n;

  double *fene;
  double ***enek;

  double criteria_BAR=1.0e-7;
  int MAXITE=1000;

  int n_total=0;
  double *W,*WT;

  int frame;
  double max,min;
  double width;
  double **oned_data,*hist;

  double *c1,*c2;

  char *progname;
  char *enekfilename,*eneklistfilename,*fenefilename;
  char *datalistfilename,*datafilename,*pmffilename;

  FILE *enekfile,*eneklistfile,*fenefile;
  FILE *datalistfile,*datafile,*pmffile;

  int ncid,fene_varid;
  size_t start[1],count[1];

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c=getopt(argc,argv,"hup:i:t:x:m:c:"))!=-1) {
    switch(c) {
    case 'i':
      interval=atoi(optarg);
      break;
    case 't':
      criteria_BAR=atof(optarg);
      break;
    case 'm':
      macminflag=ON;
      max=atof(optarg);
      break;
    case 'c':
      macminflag=ON;
      min=atof(optarg);
      break;
    case 'x':
      MAXITE=atoi(optarg);
      break;
    case 'p':
      interval=atoi(optarg);
      break;
    case 'u':
      udflag=ON;
      break;
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }


  argc-=optind;
  argv+=optind;

  if (argc < 6) {
    USAGE(progname);
    exit(1);
  }
  n_sim = atoi(*argv);
  width = atof(*++argv);
  fenefilename = *++argv;
  eneklistfilename = *++argv;
  datalistfilename = *++argv;
  pmffilename = *++argv;

  fene=(double *)gcemalloc(sizeof(double)*n_sim);
  start[0]=0;
  count[0]=n_sim;
  enc_open(fenefilename,NC_NOWRITE,&ncid);
  nc_inq_varid(ncid,"fene",&fene_varid);
  nc_get_vara_double(ncid,fene_varid,start,count,fene);
  encclose(ncid);

  n=(int *)gcemalloc(sizeof(int)*n_sim);
  enek=(double ***)gcemalloc(sizeof(double **)*n_sim);
  for (i=0;i<n_sim;++i) enek[i]=(double **)gcemalloc(sizeof(double *)*n_sim);

  for (i=0;i<n_sim;++i) for (j=0;j<n_sim;++j) enek[i][j]=(double *)gcemalloc(sizeof(double));

  oned_data=(double **)gcemalloc(sizeof(double *)*n_sim);
  for (i=0;i<n_sim;++i) oned_data[i]=(double *)gcemalloc(sizeof(double)*1);

  eneklistfile=efopen(eneklistfilename,"r");
  for (i=0;i<n_sim;++i) {
    for (j=0;j<n_sim;++j) {
      getline(&line,&len,eneklistfile);
      line[strlen(line)-1]='\0';
      enekfile=efopen(line,"r");
      getline(&line,&len,enekfile);
      k=0;
      num=0;
      d1 = 1;
      while ( d1 != -1  )  {
	d1=fscanf(enekfile,"%lf",&f1);
	d1=fscanf(enekfile,"%lf",&f1);
	if (k%interval == 0) {
	  enek[j][i]=(double *)gcerealloc(enek[j][i],sizeof(double)*(num+1));
	  enek[j][i][num]=f1;
	  ++num;
	}
	++k;
      } 
      fclose(enekfile);
      n[i]=num-1;
    }
  }
  fclose(eneklistfile);

  datalistfile=efopen(datalistfilename,"r");
  for (i=0;i<n_sim;++i) {
    getline(&line,&len,datalistfile);
    line[strlen(line)-1]='\0';
    datafile=efopen(line,"r");
    k=0;
    num=0;
    d1 = 1;
    d2 = 1;
    while ( d1 != -1 )  {
      d1=fscanf(datafile,"%lf",&f1);
      //      d2=fscanf(datafile,"%lf",&f2);
      if (k%interval == 0) {
	oned_data[i]=(double *)gcerealloc(oned_data[i],sizeof(double)*(num+1));
	oned_data[i][num]=f1;
	++num;
      }
      ++k;
    }
    fclose(datafile);
  }
  fclose(datalistfile);

  n_total=0;for (i=0;i<n_sim;++i) n_total+=n[i];

  //  if (macminflag==OFF)
  hist=MBAR_AVE_oned(fene,enek,n_sim,n,oned_data,width,&max,&min,&frame);
  //  else
  //    hist=MBAR_AVE_oned_wmaxmin(fene,enek,n_sim,n,oned_data,width,max,min,frame);

  pmf=(double *)gcemalloc(sizeof(double)*frame);
  /***************************************************************************************************/
  /* pmf_min=0.0;										     */
  /* for (i=0;i<frame[0]*frame[1];++i) {							     */
  /*   if ( pmf_min == 0.0 && hist[i] > 0.0 )							     */
  /* 	pmf_min=hist[i];									     */
  /*   if ( pmf_min > hist[i] && hist[i] > 0.0 )						     */
  /* 	pmf_min=hist[i];									     */
  /* }												     */
  /* 												     */
  /* for (i=0;i<frame[0]*frame[1];++i) {							     */
  /*   if (hist[i]!=0.0) {									     */
  /*     pmf[i]=-log(hist[i])+log(pmf_min);							     */
  /*   }											     */
  /* }												     */
  /* 												     */
  /* pmffile=efopen(pmffilename,"w");								     */
  /* fprintf(pmffile,"# histx histy \n");							     */
  /* for (i=0;i<=frame[0];++i) {								     */
  /*   for (j=0;j<frame[1];++j)									     */
  /*     fprintf(pmffile,"%e %e %e\n",width[0]*i+min[0],width[1]*j+min[1],pmf[i*frame[1]+j]);	     */
  /*   fprintf(pmffile,"\n");									     */
  /* }												     */
  /* fclose(pmffile);										     */
  /***************************************************************************************************/

  for (i=0;i</*=*/frame;++i) {
    if (hist[i]!=0.0) {
      //	pmf[i*frame[1]+j]=/*-*/log(hist[i*frame[1]+j]);
      pmf[i]=-log(hist[i]);
    }
  }

  pmf_min=pmf[0];
  pmf_max=pmf[0];
  index_max=0;
  for (i=0;i</*=*/frame;++i) {
    if ( pmf_min < pmf[i] && pmf[i]!=0 )
      pmf_min=pmf[i];
    if ( pmf_max > pmf[i] && pmf[i]!=0 ) {
	pmf_max=pmf[i];
	index_max=i;
    }
  }

  for (i=0;i</*=*/frame;++i) {
    if (hist[i]!=0.0) {
      if (udflag==OFF)
	pmf[i]-=pmf_min;
      else
	pmf[i]-=pmf_max;
    }
  }

  pmffile=efopen(pmffilename,"w");
  for (i=0;i<frame;++i) {
    if (pmf[i]!=0.0)
      fprintf(pmffile,"%10.4lf %10.4lf\n",min+width*i,pmf[i]);
    else {
      if ( udflag==OFF )
	fprintf(pmffile,"%10.4lf 0.0\n",min+width*i);
      else {
	if ( i!=index_max )
	  fprintf(pmffile,"%10.4lf ******\n",min+width*i);
	else
	  fprintf(pmffile,"%10.4lf 0.0\n",min+width*i);
      }
    }
  }
  fclose(pmffile);

  return 0;
}

void USAGE(char *progname) {
  printf("[-m fenefilename] -- w mbar ite\n");
  printf("[-t criteria_BAR ] -- criteria of MBAR iteration (default 10-7)\n");
  printf("[-x max_ite_step_of_BAR ] -- max num of MBAR iteration\n");
  printf("[-p] -- interval\n");
  printf("[-h] -- help\n");
  printf("USAGE: %s  n_sim width fenefilename eneklistfilename datalistfilename pmffilename \n", progname);
}
