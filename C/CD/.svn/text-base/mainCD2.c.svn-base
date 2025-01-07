#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "MB.h"
#include "PT.h"
#include "EF.h"
#include "IO.h"

#define ON 1
#define OFF 0

int usage(void);

int main(int argc, char *argv[]) {
  int i,j,k,l,m,nt,dummy,num,initialstep=0;
  int numatom,numstep,*numdihed;
  int flag='P',flagcn='n',flagof='f',amberflag=OFF,flagna=OFF;
  int **adpairs;
  int interval=1,pickinterval;

  double atom[4][3];
  double *protcoord;
  double theta,**old_theta;
  double pi;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*outputfilename,*parmfilename,*condfilename;
  FILE *inputfile,*outputfile,*logfile,*parmfile,*condfile;

  while((c=getopt(argc,argv,"ACPOKcnmflt:i:p:o:d:k:a:"))!=-1) {
    switch(c) {
    case 'A':
      amberflag=ON;
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
      parmfilename=optarg;
      break;
    case 'o':
      outputfilename=optarg;
    case 'd':
      condfilename=optarg;
      break;
    /*******************************/
    /* case 's':		   */
    /*   initialstep=atoi(optarg); */
    /*   break;			   */
    /*******************************/
    default:
      usage();
      exit(1);
    }
  }

  if (flagna==OFF) {
    parmfile=efopen(parmfilename,"r");
    readParmtop(parmfile);
    fclose(parmfile);
    numatom=AP.NATOM;
  }
  protcoord=(double *)gcemalloc(sizeof(double)*numatom*3);

  adpairs=(int **)gcemalloc(sizeof(int *)*5);
  if (flag!='c') {
    adpairs[0]=(int *)gcemalloc(sizeof(int)*4); // PHI,PSI
    adpairs[1]=(int *)gcemalloc(sizeof(int)*4); // OMEGA
    adpairs[2]=(int *)gcemalloc(sizeof(int)*4); // KI
    adpairs[3]=(int *)gcemalloc(sizeof(int)*8); // ACE
    adpairs[4]=(int *)gcemalloc(sizeof(int)*8); // NME
    numdihed=(int *)gcemalloc(sizeof(int)*5);  
    readdihedpairs(adpairs,numdihed);
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

  inputfile=efopen(inputfilename,"r");
  if (amberflag==ON)
    getline(&line,&len,inputfile);

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
  if (initialstep>0) {
    for (i=0;i<initialstep;++i) {
      if (flagof=='c') fprintf(outputfile,"%d ",i);
      io_scanconf(inputfile,numatom,protcoord,'x');
    }
  }

  for (i=initialstep;i<numstep;++i) {
    if (flagof=='c') fprintf(outputfile,"%d ",i);
    io_scanconf(inputfile,numatom,protcoord,'x');
    
    if ( (i%interval) == 0) {
      for (j=0;j<nt;++j) {
	for (k=0;k<numdihed[j];++k) {
	  for (m=0;m<4;++m) 
	    for (l=0;l<3;++l)
	      atom[m][l]=protcoord[(adpairs[j][k*4+m])*3+l];
      
	  if (i==0) old_theta[j][k]=theta;
	  if (flagcn=='n') {
	    theta=pick_dihed(atom[0],atom[1],atom[2],atom[3],0,old_theta[j][k]);
	    if (theta > pi)  theta = theta -2.0*pi;
	  }
	  else if (flagcn=='c') theta=pick_dihed(atom[0],atom[1],atom[2],atom[3],1,old_theta[j][k]);
	  
	  fprintf(outputfile,"%e ",theta*180/acos(-1.0));
	  old_theta[j][k]=theta;
	}
      }
      fprintf(outputfile,"\n");
    }
  }
  fclose(outputfile);
  fclose(inputfile);

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

int usage(void) {

  printf("USAGE:CD \n");
  printf("-A -C -P -O -K \n");
  printf("-c -n -m -f -l \n");
  printf("-i [crd] -p [top] -a [numatom] -t [numstep] -o [output] -d [cond] -k [111111 interval] -s [initialstep]  \n");

}

 
