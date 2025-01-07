
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "EF.h"
#include "IO.h"

#include "MBAR.h"

#define ON 1
#define OFF 0

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,t;
  int flag=OFF,iflag=OFF,jflag=OFF;
  int interval=1;
  int d=0,num;
  double sum=0.0;
  double f;

  double din;

  int n_sim,numk=0;
  int *n;

  double T=300,*T_sim;
  double k_B=1.98723e-3;
  double beta;

  double *fene;
  double ***enek,**expeneraw;

  int frame;
  double max,min;
  double width;
  double **od_data,*hist,*pmf;
  double pmf_min;
  double *betaene;

  int n_total=0;
  double *W,*WT;
  double *hist_error,*pmf_error;

  double A,*dA;

  double criteria_BAR=1.0e-7;
  int MAXITE=1000;

  char *progname;
  char *betaenefilename;
  char *enekfilename,*eneklistfilename,*fenefilename;
  char *datalistfilename,*datafilename;
  char *histfilename,*pmffilename;
  char *errorfilename,*perrorfilename;
  char *Tsimlistfilename;

  FILE *betaenefile;
  FILE *enekfile,*eneklistfile,*fenefile;
  FILE *datalistfile,*datafile;
  FILE *histfile,*pmffile;
  FILE *errorfile,*perrorfile;
  FILE *Tsimlistfile;

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c=getopt(argc,argv,"hijmt:p:"))!=-1) {
    switch(c) {
    case 'm':
      flag=ON;
      break;
    case 'i':
      iflag=ON;
      break;
    case 'j':
      jflag=ON;
      break;
    case 'p':
      interval=atoi(optarg);
      break;
    case 't':
      T=atof(optarg);
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

  if (argc < /*8*/9) {
    USAGE(progname);
    exit(1);
  }
  //  betaenefilename = *argv;
  numk = atoi(*/*++*/argv)-1;
  n_sim = atoi(*++argv);
  width = atof(*++argv);
  fenefilename = *++argv;
  eneklistfilename = *++argv;
  Tsimlistfilename = *++argv;
  datalistfilename = *++argv;
  histfilename = *++argv;
  pmffilename = *++argv;

  errorfilename = *++argv;
  perrorfilename = *++argv;

  fene=(double *)gcemalloc(sizeof(double)*n_sim);
  if (flag==OFF) {
    fenefile=efopen(fenefilename,"r");
    getline(&line,&len,fenefile);
    for (i=0;i<n_sim;++i) {
      fscanf(fenefile,"%d %lf",&d,&fene[i]);
      fene[i]=exp(fene[i]);
    }
    fclose(fenefile);
  }

  n=(int *)gcemalloc(sizeof(int)*n_sim);
  enek=(double ***)gcemalloc(sizeof(double **)*n_sim);
  for (i=0;i<n_sim;++i) {
    enek[i]=(double **)gcemalloc(sizeof(double *)*n_sim);
  }
  for (i=0;i<n_sim;++i) {
    for (j=0;j<n_sim;++j) {
      enek[i][j]=(double *)gcemalloc(sizeof(double ));
    }
  }
  od_data=(double **)gcemalloc(sizeof(double *)*n_sim);
  for (i=0;i<n_sim;++i) {
    od_data[i]=(double **)gcemalloc(sizeof(double ));
  }
  //  betaene=(double *)gcemalloc(sizeof(double));
  expeneraw=(double **)gcemalloc(sizeof(double *)*n_sim);
  T_sim=(double *)gcemalloc(sizeof(double)*n_sim);

  Tsimlistfile=efopen(Tsimlistfilename,"r");
  for (i=0;i<n_sim;++i) {
    fscanf(Tsimlistfile,"%lf",&T_sim[i]);
  }
  fclose(Tsimlistfile);

  eneklistfile=efopen(eneklistfilename,"r");
  for (i=0;i<n_sim;++i) {
    getline(&line,&len,eneklistfile);
    l=strlen(line);
    line[l-1]='\0';
    enekfile=efopen(line,"r");
    getline(&line,&len,enekfile);
    k=0;
    num=0;
    d = 1;
    while ( d != -1  )  {
      d=fscanf(enekfile,"%d",&l);
      d=fscanf(enekfile,"%lf",&f);
      if (k%interval == 0) {
	expeneraw[i]=(double *)gcerealloc(expeneraw[i],sizeof(double)*(/*k*/num+1));
	expeneraw[i][num]=f;
	++num;
      }
      ++k;
    } 
    fclose(enekfile);
    n[i]=num-1;
  }
  fclose(eneklistfile);

  for (i=0;i<n_sim;++i) {
    beta=1.0/k_B*(1.0/T_sim[i]-1.0/T);
    for (j=0;j<n_sim;++j) {
      num=0;
      for (k=0;k<n[j];++k) {
	f=expeneraw[j][k];
	enek[i][j]=(double *)gcerealloc(enek[i][j],sizeof(double)*(num+1));
	enek[i][j][num]=exp(-1.0*f*beta);
	++num;
      } 
    }
  }

  //  betaenefile=efopen(betaenefilename,"r");

  //  k=0;
  //  num=0;
  //  d = 1;
  //  while ( d != -1  )  {
  //    if (iflag==ON) d=fscanf(betaenefile,"%d",&l);
  //    d=fscanf(betaenefile,"%lf",&f);
  //    if (k%interval == 0) {
  //      betaene=(double *)gcerealloc(betaene,sizeof(double)*(num+1));
  //      betaene[num]=/*f*/0.0;
  //      //      betaene[num]=f*1.0/k_B*(1.0/T);
  //      ++num;
  //    }
  //    ++k;
  //  } 
  //  fclose(betaenefile);

  datalistfile=efopen(datalistfilename,"r");
  for (i=0;i<n_sim;++i) {
    getline(&line,&len,datalistfile);
    line[strlen(line)-1]='\0';
    datafile=efopen(line,"r");
    k=0;
    num=0;
    d = 1;
    while ( d != -1  )  {
      if (jflag==ON)
	d=fscanf(datafile,"%d",&l);
      d=fscanf(datafile,"%lf",&f);
      if (k%interval == 0) {
	od_data[i]=(double *)gcerealloc(od_data[i],sizeof(double)*(num+1));
	od_data[i][num]=f;
	++num;
      }
      ++k;
    } 
    fclose(datafile);
  }
  fclose(datalistfile);

  dA=gcemalloc(sizeof(double)*n_sim); 
  for (i=0;i<n_sim;++i) fene[i]=1.0;
  MBAR_A(fene, enek, n_sim, n, &A, dA);

  if (flag==ON) {
    MBAR_ite_high_speed_2(fene,enek,n_sim,n,criteria_BAR,MAXITE);
  }
  for (i=0;i<n_sim;++i) n_total+=n[i];
    
  W=gcemalloc(sizeof(double)*n_total*(n_sim+2)); 
  WT=(double *)gcemalloc(sizeof(double)*n_total*(n_sim+2));
    
  l=-1;
  for (i=0;i<n_sim;++i) {
    for (t=0;t<n[i];++t) {
      ++l;
      for (j=0;j<n_sim;++j) {
  	din=0.0;
  	for (k=0;k<n_sim;++k) {
  	  din+=n[k]*exp(fene[k]-enek[k][i][t]);
  	}
  	W[l*(n_sim+2)+j]=exp(fene[j]-enek[j][i][t])/din;
	WT[j*n_total+l] =W[l*(n_sim+2)+j];
      }
    }
  }

  hist=MBAR_AVE_oned_multi_temp_2(fene,enek,n_sim,n,od_data,width,&max,&min,&frame,numk);
  //  hist_error=MBAR_ACM2_histod(fene,enek,n_sim,n,od_data,width,&max,&min,&frame,W,WT);

  /**************************************/
  /* sum=0.0;			        */
  /* for (i=0;i<frame;++i) {	        */
  /*   sum+=hist[i];		        */
  /* }				        */
  /* for (i=0;i<frame;++i) {	        */
  /*   hist_error[i]=hist_error[i]/sum; */
  /* }				        */
  /**************************************/

  pmf=(double *)gcemalloc(sizeof(double)*frame);
  pmf_error=(double *)gcemalloc(sizeof(double)*frame);
  /*****************************************************/
  /* pmf_min=0.0;				       */
  /* for (i=0;i<=frame;++i) {			       */
  /*   if ( pmf_min == 0.0 && hist[i] > 0.0 )	       */
  /* 	pmf_min=hist[i];			       */
  /*   if ( pmf_min > pmf[i] && hist[i] > 0.0 )	       */
  /* 	pmf_min=hist[i];			       */
  /* }						       */
  /*****************************************************/

  pmf_min=0.0;
  for (i=0;i<=frame;++i) {
    if ( pmf_min == 0.0 && hist[i] > 0.0 )
	pmf_min=hist[i];
    if ( pmf_min < hist[i] && hist[i] > 0.0 )
	pmf_min=hist[i];
  }

  for (i=0;i<=frame;++i) {
    if (hist[i]!=0.0) {
      pmf[i]=-log(hist[i])+log(pmf_min);
    }
  }
  /*   if (hist_error[i]!=0.0) {	    */
  /*     pmf_error[i]=-log(hist_error[i]);  */
  /*   }				    */
  /* }					    */
  /******************************************/

  if (flag==ON) {
    fenefile=efopen(fenefilename,"w");
    fprintf(fenefile,"# fene\n");
    for (i=0;i<n_sim;++i)
      fprintf(fenefile,"%d %12.10lf\n",i,fene[i]);
    fclose(fenefile);
  }

  histfile=efopen(histfilename,"w");
  fprintf(histfile,"# hist \n");
  for (i=0;i<frame;++i)
    fprintf(histfile,"%10.4lf %10.4lf\n",min+width*i,hist[i]);
  fclose(histfile);

  pmffile=efopen(pmffilename,"w");
  fprintf(pmffile,"# pmf \n");
  for (i=0;i<frame;++i)
    fprintf(pmffile,"%10.4lf %10.4lf\n",min+width*i,pmf[i]);
  fclose(pmffile);

  /*******************************************************************/
  /* errorfile=efopen(errorfilename,"w");			     */
  /* fprintf(errorfile,"# hist \n");				     */
  /* for (i=0;i<frame;++i)					     */
  /*   fprintf(errorfile,"%10.4lf %10.4lf\n",width*i,hist_error[i]); */
  /* fclose(errorfile);						     */
  /* 								     */
  /* perrorfile=efopen(perrorfilename,"w");			     */
  /* fprintf(perrorfile,"# hist \n");				     */
  /* for (i=0;i<frame;++i)					     */
  /*   fprintf(perrorfile,"%10.4lf %10.4lf\n",width*i,pmf_error[i]); */
  /* fclose(errorfile);						     */
  /*******************************************************************/
  
  return 0;
}

void USAGE(char *progname) {
  printf("-p -- interval\n");
  printf("-h -- help\n");
  printf("USAGE: %s n_sim width fenefilename enelistfilename datalistfilename histfilename histfilename pmffilename errorfile perrorfilename \n", progname);
}

