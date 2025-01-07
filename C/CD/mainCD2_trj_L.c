#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

//#include "netcdf_mine.h"
#include "netcdf_mineL.h"

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
  int flag='P',flagcn='n',flagof='f',amberflag=OFF,flagna=OFF,Drestoutflag=OFF,crdflag=OFF,amberncflag=OFF,Radflag=OFF;
  int ambcrdflag=OFF,pdbflag=OFF;
  int flagp=3;
  int flagmo=ON;
  int **adpairs;
  int interval=1,pickinterval;

  double fc;
  double atom[4][3];
  double *protcoord;
  double theta,**old_theta;
  double pi;

  PDBF PD;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double crd_nc[MAXATOM][3],*crd;

  char *inputfilename,*outputfilename,*parmfilename,*condfilename;
  char *progname;
  FILE *inputfile,*outputfile,*logfile,*parmfile,*condfile;

  while((c=getopt(argc,argv,"hAIBCLPORKNcnmflD:t:i:p:o:d:k:v:a:"))!=-1) {
    switch(c) {
    case 'A':
      amberflag=ON;
      break;
    case 'I':
      crdflag=ON;
      break;
    case 'B':
      ambcrdflag=ON;
      break;
    case 'N':
      amberncflag=ON;
      break;
    case 'D':
      Drestoutflag=ON;
      fc=atof(optarg);
      break;
    case 'L':
      pdbflag=ON;
      numstep=1;
      interval=1;
      break;
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
    case 'R':
      Radflag=ON;
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
    case 'm':
      flagcn=c;
      break;
    case 'f':
      flagof='c';
      break;
    case 'l':
      flagof='x';
      break;
    case 't':
      flagmo=OFF;
      numstep=atoi(optarg);
      break;
    case 'a':
      flagna=ON;
      numatom=atoi(optarg);
      break;
    case 'v':
      interval=atoi(optarg);
      break;
    case 'k':
      //      pickinterval=atoi(optarg);
      initialstep=atoi(optarg);
      break;
    case 'i':
      inputfilename=optarg;
      break;
    case 'p':
      flagp=ON;
      parmfilename=optarg;
      break;
    case 'o':
      outputfilename=optarg;
    case 'd':
      condfilename=optarg;
      break;
    case 'h':
      USAGE(progname);
      exit(1);
    /*******************************/
    /* case 's':		   */
    /*   initialstep=atoi(optarg); */
    /*   break;			   */
    /*******************************/
    default:
      USAGE(progname);
      exit(1);
    }
  }

  if (Drestoutflag==ON && ambcrdflag==OFF) {
    printf("option error! -D must -B\n");
    exit(1);
  }

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  numstep = atoi(*argv);
  inputfilename  = *++argv;
  outputfilename = *++argv;

  /**********************************************************************************************/
  /* if (flagmo==ON && ambcrdflag==OFF && pdbflag==OFF && crdflag==OFF && amberncflag==OFF)     */
  /*   numstep=myncL_get_present_step_MCD(inputfilename,&nc_id_MCD);			        */
  /* else if (flagmo==ON && ambcrdflag==OFF && pdbflag==OFF && crdflag==OFF && amberncflag==ON) */
  /*   numstep=myncL_get_present_step_AMBER(inputfilename,&nc_id_MD);			        */
  /* else if (ambcrdflag==ON || pdbflag==ON || crdflag==ON)				        */
  /*   numstep=1;									        */
  /**********************************************************************************************/

  if (flagna==OFF && flag!='c') {
    parmfile=efopen(parmfilename,"r");
    readParmtopL(parmfile);
    fclose(parmfile);
    numatom=AP.NATOM;
  }
  //  crd_nc=(double **)gcemalloc(sizeof(double *)*numatom);
  //  for (i=0;i<numatom;++i) crd_nc[i]=(double *)gcemalloc(sizeof(double)*3);

  adpairs=(int **)gcemalloc(sizeof(int *)*5);
  if (flag!='c') {
    adpairs[0]=(int *)gcemalloc(sizeof(int)*4); // PHI,PSI
    adpairs[1]=(int *)gcemalloc(sizeof(int)*4); // OMEGA
    adpairs[2]=(int *)gcemalloc(sizeof(int)*4); // KI
    adpairs[3]=(int *)gcemalloc(sizeof(int)*8); // ACE
    adpairs[4]=(int *)gcemalloc(sizeof(int)*8); // NME
    numdihed=(int *)gcemalloc(sizeof(int)*5);  
    readdihedpairsL(adpairs,numdihed);
  }
  else {
    condfile=efopen(condfilename,"r");
    numdihed=(int *)gcemalloc(sizeof(int));  
    nt=1;
    fscanf(condfile,"%d",&num);
    numdihed[0]=num;
    adpairs[0]=(int *)gcemalloc(sizeof(int)*4*numdihed[0]);
    for (i=0;i<numdihed[0];++i)
      for (j=0;j<4;++j)
	fscanf(condfile,"%d",&adpairs[0][i*4+j]);
    fclose(condfile);
  }
  
  old_theta=(int **)gcemalloc(sizeof(int *)*5);
  old_theta[0]=gcemalloc(sizeof(double)*numdihed[0]);
  old_theta[1]=gcemalloc(sizeof(double)*numdihed[1]);
  old_theta[2]=gcemalloc(sizeof(double)*numdihed[2]);
  old_theta[3]=gcemalloc(sizeof(double)*numdihed[3]);
  old_theta[4]=gcemalloc(sizeof(double)*numdihed[4]);

  outputfile=efopen(outputfilename,"w");
  if (flagof=='c') {
    fprintf(outputfile,"numstep ");
    k=0;
    for (i=0;i<nt;++i) {
      for (j=0;j<numdihed[i];++j) {
	++k;
	fprintf(outputfile,"%d ",k);
      }
    }
    fprintf(outputfile,"\n");
  }

  pi=acos(-1.0);

  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);

  for (i=0;i<numstep;++i) {
    if (flagof=='c') fprintf(outputfile,"%d ",i);

    for (j=0;j<numatom;++j) {
      for (k=0;k<3;++k) {
	fscanf(inputfile,"%lf",&crd_nc[j][k]);
      }
    }
    
    if ( (i%interval) == 0 || i==initialstep || ambcrdflag==ON ) {
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
	  else if (flagcn=='c') theta=pick_dihed(atom[0],atom[1],atom[2],atom[3],1,old_theta[j][k]);
	  
	  if (Drestoutflag==ON) {
	    fprintf(outputfile,"\n%d %d %d %d ",adpairs[j][k*4]+1,adpairs[j][k*4+1]+1,adpairs[j][k*4+2]+1,adpairs[j][k*4+3]+1);
	  }
	  if (Radflag==OFF)
	    fprintf(outputfile,"%e ",theta*180/acos(-1.0));
	  else 
	    fprintf(outputfile,"%e ",theta);

	  if (numstep==1) {
	    if (flag=='P') {
	      if (k%2==1)
		fprintf(outputfile,"\n");
	    }
	    else 
	      fprintf(outputfile,"\n");
	  }
	  old_theta[j][k]=theta;
	}
      }
      fprintf(outputfile,"\n");
    }
  }
  fclose(inputfile);
  fclose(outputfile);

  logfile=efopen("log_CD.txt","w");
  for (i=0;i<nt;++i)
    fprintf(logfile,"%d ",numdihed[i]);
  fprintf(logfile,"\n");
  for (i=0;i<nt;++i) {
    for (j=0;j<numdihed[i];++j) {
      for (k=0;k<4;++k) 
	fprintf(logfile,"%4d ",adpairs[i][j*4+k]+1);
      for (k=0;k<4;++k) 
	fprintf(logfile,"%s ",AP.IGRAPH[adpairs[i][j*4+k]-1]);
      fprintf(logfile,"\n");
    }
  }
  fclose(logfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-A] amberflag \n");
  printf("[-I] crdflag \n");
  printf("[-B] amcrdflag \n");
  printf("[-N] amberncflag \n");
  printf("[-L] pdbfilflag \n");
  printf("[-D] Drestoutflag \n");
  printf("[-C] specify dihed in cond file \n");
  printf("[-P] phi psi \n");
  printf("[-O] phi psi omega \n");
  printf("[-K] phi psi omega kai \n");
  printf("[-c] continue mode ON \n");
  printf("[-n] continue mode OFF \n");
  printf("[-f] formated file \n");
  printf("[-l] unformated file \n");
  printf("[-t numstep] specify step \n");
  printf("[-a numatom] specify num of atom \n");
  printf("[-d condfile] specify cond file \n");
  printf("[-p topfilefile] specify top file \n");
  printf("[-v interval] specify interval \n");
  printf("[-k interval] specify initial \n");
  printf("%s inputfilename outputfilename \n",progname);
}

 
