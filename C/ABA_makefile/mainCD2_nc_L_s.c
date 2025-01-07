#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "netcdf_mine.h"

#include "PDB.h"
#include "MB.h"
#include "PTL.h"
#include "EF.h"
#include "IO.h"

#define ON 1
#define OFF 0


int usage(void);

int main(int argc, char *argv[]) {
  int i,j,k,l,d,m,nt,dummy,num,initialstep=0;
  int numatom,numstep,*numdihed;
  int flag='P',flagcn='n',flagof='f';
  int **adpairs;
  int interval=1;

  double fc;
  double atom[4][3];
  double *protcoord;
  double theta,**old_theta;
  double pi;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double crd_nc[MAXATOM][3],*crd;
  struct my_netcdf_out_id_AMBER nc_id;

  char *inputfilename,*outputfilename,*parmfilename;
  char logfilename[100];
  char *progname;
  FILE *inputfile,*outputfile,*logfile,*parmfile;

  while((c=getopt(argc,argv,"hCPOKNcnmflD:t:i:p:o:d:k:a:"))!=-1) {
    switch(c) {
    case 'C':
      flag='c';
      break;
    case 'P':
      flag='P';
      nt=1;
      break;
    case 'O':
      flag='O';
      nt=2;
      break;
    case 'K':
      flag='K';
      nt=5;
      break;
    case 'c':
      flagcn=c;
      break;
    case 'n':
      flagcn=c;
      break;
    case 'v':
      interval=atoi(optarg);
      break;
    case 'k':
      initialstep=atoi(optarg);
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

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  inputfilename = *argv;
  parmfilename  =  *++argv;
  outputfilename = *++argv;

  pi=acos(-1.0);

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;

  numstep=mync_get_present_step_AMBER(inputfilename,&nc_id);

  adpairs=(int **)gcemalloc(sizeof(int *)*5);
  adpairs[0]=(int *)gcemalloc(sizeof(int)*4); // PHI,PSI
  adpairs[1]=(int *)gcemalloc(sizeof(int)*4); // OMEGA
  adpairs[2]=(int *)gcemalloc(sizeof(int)*4); // KI
  adpairs[3]=(int *)gcemalloc(sizeof(int)*8); // ACE
  adpairs[4]=(int *)gcemalloc(sizeof(int)*8); // NME
  numdihed=(int *)gcemalloc(sizeof(int)*5);  
  readdihedpairsL(adpairs,numdihed);
  
  old_theta=(int **)gcemalloc(sizeof(int *)*5);
  old_theta[0]=gcemalloc(sizeof(double)*numdihed[0]);
  old_theta[1]=gcemalloc(sizeof(double)*numdihed[1]);
  old_theta[2]=gcemalloc(sizeof(double)*numdihed[2]);
  old_theta[3]=gcemalloc(sizeof(double)*numdihed[3]);
  old_theta[4]=gcemalloc(sizeof(double)*numdihed[4]);

  outputfile=efopen(outputfilename,"w");

  for (i=initialstep;i<numstep;++i) {
    mync_open_inq_get_sh_AMBER(inputfilename,numatom,i,1,i+1,&nc_id,crd_nc);

    if ( (i%interval) == 0 ) {
      for (j=0;j<nt;++j) {
	for (k=0;k<numdihed[j];++k) {
	  for (m=0;m<4;++m) 
	    for (l=0;l<3;++l)
	      atom[m][l]=crd_nc[(adpairs[j][k*4+m])][l];
      
	  if (i==0) old_theta[j][k]=theta;
	  if (flagcn=='n') {
	    theta=pick_dihed(atom[0],atom[1],atom[2],atom[3],0,old_theta[j][k]);
	    if (theta > pi)  theta = theta -2.0*pi;
	  }
	  else if (flagcn=='c') {
	    theta=pick_dihed(atom[0],atom[1],atom[2],atom[3],1,old_theta[j][k]);
	  }
	  
	  fprintf(outputfile,"%e ",theta*180/acos(-1.0));
	  old_theta[j][k]=theta;
	}
      }
      fprintf(outputfile,"\n");
    }
  }
  fclose(outputfile);

  sprintf(logfilename,"%s_logCD.txt",outputfilename);
  logfile=efopen(logfilename,"w");
  for (i=0;i<nt;++i)
    fprintf(logfile,"%d ",numdihed[i]);

  fprintf(logfile,"\n");
  for (i=0;i<nt;++i) {
    for (j=0;j<numdihed[i];++j) {
      if (i==0) fprintf(logfile,"ph %d=",j+1);
      else if (i==1) fprintf(logfile,"omg%d=",j+1);
      else if (i==2) fprintf(logfile,"kai%d=",j+1);
      else if (i==3) fprintf(logfile,"ace%d=",j+1);
      else if (i==4) fprintf(logfile,"nme%d=",j+1);
      for (k=0;k<4;++k) 
	fprintf(logfile,"%4d ",adpairs[i][j*4+k]+1);
      for (k=0;k<4;++k) 
	fprintf(logfile,"%s ",AP.IGRAPH[adpairs[i][j*4+k]]);
      fprintf(logfile,"\n");
    }
  }
  fclose(logfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-P] phi psi \n");
  printf("[-O] phi psi omega \n");
  printf("[-K] phi psi omega kai \n");
  printf("[-c] continue mode ON \n");
  printf("[-n] continue mode OFF \n");
  printf("[-k interval] specify interval \n");
  printf("%s inputfilename parmfilename outputfilename \n",progname);
}

 
