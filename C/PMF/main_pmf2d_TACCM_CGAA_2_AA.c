
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "netcdf_mine.h"
#include "IO.h"

#include "MBAR.h"

#include "PMF.h"
#include "PT.h"
#include "EF.h"
#include "mymath.h"

double *io_scandcoldata2(FILE *inputfile,int numi,int numcol,int xcol,int ycol,int *numstep,double *data);

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,t;
  int interval=1;
  int d1=0,d2=0,num;
  double f1,f2,sum;
  double *pmfCGAA,pmf_min;

  double din;

  int n_sim;
  int *n;

  double *fene;
  double ***enek;

  double criteria_BAR=1.0e-7;
  int MAXITE=1000;

  int n_total=0;
  double *W,*WT;

  int *frame;
  double *max,*min;
  double *width;
  double ***twod_data,*hist_CGAA;

  double *c1,*c2;

  int numstep;
  double *data;
  int numi=1,numcol=2,xcol=1,ycol=2;

  double *pmfCG,*pmfAA;

  double pi;

  double T=0.0,beta=1.0;
  double k_B=1.98723e-3;

  double UNITT=418.4070;

  char *enekfilename,*eneklistfilename,*fenefilename;
  char *datalistfilename,*datafilename,*pmffilename;

  FILE *enekfile,*eneklistfile,*fenefile;
  FILE *datalistfile,*datafile,*pmffile;

  char *datafilename_CG,*pmffilename_AA;
  FILE *datafile_CG,*pmffile_AA;

  int ncid,fene_varid;
  size_t start[1],count[1];

  char *line;
  size_t len=0;

  int d;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *progname,*logfilename;

  FILE *logfile;

  int opt_idx=1;

  frame=(int *)gcemalloc(sizeof(int)*2);
  max=(double *)gcemalloc(sizeof(double)*2);
  min=(double *)gcemalloc(sizeof(double)*2);
  width=(double *)gcemalloc(sizeof(double)*2);

  struct option long_opt[] = {
    {"widthx",1,NULL,'w'},
    {"widthy",1,NULL,'z'},
    {"interval",1,NULL,'i'},
    {"temp",1,NULL,'t'},
    {"beta",1,NULL,'b'},
    {"numi",1,NULL,'I'},
    {"numcol",1,NULL,'c'},
    {"xcol",1,NULL,'x'},
    {"ycol",1,NULL,'y'},
    {"help",0,NULL,'H'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((d=getopt_long(argc,argv,"Hhi:a:j:B:w:z:t:b:I:c:x:y:",long_opt,&opt_idx))!=-1) {
    switch(d) {
    case 'w':
      width[0]=atof(optarg);
      break;
    case 'z':
      width[1]=atof(optarg);
      break;
    case 'i':
      interval=atoi(optarg);
      break;
    case 't':
      T=atof(optarg);
      break;
    case 'b':
      beta=atof(optarg);
      break;
    case 'I':
      numi=atoi(optarg);
      break;
    case 'x':
      xcol=atoi(optarg);
      break;
    case 'y':
      ycol=atoi(optarg);
      break;
    case 'H':
      USAGE(progname);
      exit(1);
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  progname=*argv;

  pi=acos(-1.0);

  argc-=optind;
  argv+=optind;

  if (argc < 6) {
    USAGE(progname);
    exit(1);
  }  
  n_sim = atoi(*argv);
  fenefilename = *++argv;
  eneklistfilename = *++argv;
  datalistfilename = *++argv;
  datafilename_CG   = *++argv;
  pmffilename_AA  = *++argv;

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

  twod_data=(double ***)gcemalloc(sizeof(double **)*2);
  for (i=0;i<2;++i)
    twod_data[i]=(double **)gcemalloc(sizeof(double *)*n_sim);
  for (i=0;i<2;++i) for (j=0;j<n_sim;++j) twod_data[i][j]=(double *)gcemalloc(sizeof(double)*1);

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
    while ( d1 != -1 && d2 != -1 )  {
      d1=fscanf(datafile,"%lf",&f1);
      d2=fscanf(datafile,"%lf",&f2);
      if (k%interval == 0) {
	twod_data[0][i]=(double *)gcerealloc(twod_data[0][i],sizeof(double)*(num+1));
	twod_data[1][i]=(double *)gcerealloc(twod_data[1][i],sizeof(double)*(num+1));
	twod_data[0][i][num]=f1;
	twod_data[1][i][num]=f2;
	++num;
      }
      ++k;
    }
    fclose(datafile);
  }
  fclose(datalistfile);

  n_total=0;for (i=0;i<n_sim;++i) n_total+=n[i];
  pmfCGAA=MBAR_AVE_twod(fene,enek,n_sim,n,twod_data,width,max,min,frame);

  numstep=0;
  data=(double *)gcemalloc(sizeof(double)*2);
  datafile_CG=efopen(datafilename_CG,"r");
  data=io_scandcoldata2(datafile_CG,numi,numcol,xcol,ycol,&numstep,data);
  fclose(datafile_CG);
  pmfCG=pmf_2dmap_wmaxmin(data,numstep,width[0],width[1],max[0],max[0],min[1],min[1],&frame[0],&frame[1]);

  pmfAA=(double *)gcemalloc(sizeof(double)*frame[0]*frame[1]);

  for (i=0;i<frame[0];++i) {
    for (j=0;j<frame[1];++j) {
      if ( pmfCG[i*frame[1]+j] != 0.0) pmfAA[i*frame[1]+j]=pmfCGAA[i*frame[1]+j]/pmfCG[i*frame[1]+j];
      else pmfAA[i*frame[1]+j]=0.0;
    }
  }

  sum=0.0;
  for (i=0;i<frame[0];++i)
    for (j=0;j<frame[1];++j)
      sum+=pmfAA[i*(frame[1])+j];

  for (i=0;i<frame[0];++i)
    for (j=0;j<frame[1];++j)
      pmfAA[i*(frame[1])+j]/=sum;

  pmf_min=0.0;
  for (i=0;i<=frame[0];++i) {
    for (j=0;j<frame[1];++j) {
      if ( pmf_min == 0.0 && pmfAA[i*(frame[1])+j] > 0.0 )
	pmf_min=pmfAA[i*(frame[1])+j];
      if ( pmf_min > pmfAA[i*(frame[1])+j] && pmfAA[i*(frame[1])+j] > 0.0 )
	pmf_min=pmfAA[i*(frame[1])+j];
    }
  }

  for (i=0;i<=frame[0];++i)
    for (j=0;j<frame[1];++j)
      if (pmfAA[i*(frame[1])+j]!=0.0)
	pmfAA[i*(frame[1])+j]=-log(pmfAA[i*(frame[1])+j])+log(pmf_min);

  if (T>0) beta=1.0/(k_B*T);
  for (i=0;i<frame[0];++i) for (j=0;j<frame[1];++j) pmfAA[i*frame[1]+j]=1.0/beta*pmfAA[i*frame[1]+j];

  pmffile_AA=efopen(pmffilename_AA,"w");
  for (i=0;i<frame[0];++i) {
    for (j=0;j<frame[1];++j)
      fprintf(pmffile_AA,"%e %e %e\n",width[0]*i+min[0],width[1]*j+min[1],pmfAA[i*frame[1]+j]);
    fprintf(pmffile_AA,"\n");
  }
  fclose(pmffile_AA);
  
  return 0;
}

double *io_scandcoldata2(FILE *inputfile,int numi,int numcol,int xcol,int ycol,int *numstep,double *data){
  int i,j,k;
  double f;

  char *line;
  size_t len=0;

  for (i=0;i<numi;++i)
    getline(&line,&len,inputfile);

  for (i=(*numstep);;++i) {
    for (j=0;j<numcol;++j) {
      if (fscanf(inputfile,"%lf",&f)!=-1) {
	if (j==xcol-1) {
	  data=(double *)gcerealloc(data,sizeof(double)*(i+1)*2);
	  data[i*2]=f;
	}
	else if (j==ycol-1) {
	  data[i*2+1]=f;
	}
      }
      else {
	*numstep=i-1;
	return data;
      }
    }
  }
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[--minx] \n");
  printf("[--maxx] \n");
  printf("[--miny] \n");
  printf("[--maxy] \n");
  printf("[--widthx] \n");
  printf("[--widthy] \n");
  printf("[--help] help \n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename_CGAA inputfilename_CG outputfilename_AA\n",progname);
}
