#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "TOPO.h"
#include "LA.h"
#include "mymath.h"

#include "MB.h"
#include "PTL.h"
#include "EF.h"
#include "IO.h"

#define ON 1
#define OFF 0

double CDon(double atom[4][3], double pi, double *theta);

int usage(void);

int main(int argc, char *argv[]) {
  int i,j,k,l,d,dummy,initialstep=0;
  int numatom,numstep,numcrd=0;
  int absnt;

  int *numdihed,numdihedall=0;
  int flag='P',flagcn='n',flagof='f',amberflag=OFF,flagna=OFF,Drestoutflag=OFF,crdflag=OFF,amberncflag=OFF,Radflag=OFF;
  int TLflag=ON;
  int DEGflag=OFF;
  int ambcrdflag=OFF,pdbflag=OFF;
  int flagp=3;
  int flagmo=ON;
  int **adpairs;

  int m,num,nt,nd;
  double fc;

  int flagout;
  //  int **pairs;

  int interval=1;

  double *dihed_target,width;

  double atom[4][3];
  double *theta;
  double pi;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double *crd;

  char *inputfilename,*parmfilename,*condfilename;
  char *crdfilenamebase,crdfilename[100],*logfilename;
  char *progname;
  FILE *inputfile,*parmfile,*condfile,*crdfile,*logfile;

  while((c=getopt(argc,argv,"hAIBCLPOMRKNcnmflTED:t:i:p:o:d:k:v:a:"))!=-1) {
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
    case 'T':
      TLflag=OFF;
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
    case 'E':
      DEGflag=ON;
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
    case 'M':
      flag='M';
      nt=-2;
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
      //    case 'o':
      //      outputfilename=optarg;
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
  pi=acos(-1.0);

  argc-=optind;
  argv+=optind;

  if (argc < 5) {
    USAGE(progname);
    exit(1);
  }
  numstep         = atoi(*argv);
  inputfilename   = *++argv;
  width           = atof(*++argv);
  if (DEGflag==ON) width = width*pi/180.0;
  crdfilenamebase = *++argv;
  logfilename     = *++argv;
  
  if (flagna==OFF && flag!='c') {
    parmfile=efopen(parmfilename,"r");
    readParmtopL(parmfile);
    fclose(parmfile);
    numatom=AP.NATOM;
  }
  
  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

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
      for (j=0;j<4;++j) {
	fscanf(condfile,"%d",&adpairs[0][i*4+j]);
	adpairs[0][i*4+j]-=1;
      }
    fclose(condfile);
  }

  if (nt<0) {
    absnt=abs(nt)-1;
    numdihedall=numdihed[absnt];
  }
  else
    for (i=0;i<nt;++i) numdihedall+=numdihed[i];
  
  //pairs=(int **)gcemalloc(sizeof(int *)*numdihedall;
  dihed_target=(double *)gcemalloc(sizeof(double)*numdihedall);

  if (argc < 5+numdihedall) {
    USAGE(progname);
    exit(1);
  }
  for (i=0;i<numdihedall;++i) {
    dihed_target[i] = atof(*++argv);
    if (DEGflag==ON) {
      dihed_target[i] = dihed_target[i]*pi/180.0;
    }
    while (dihed_target[i] > pi) dihed_target[i]-=2.0*pi;
    while (dihed_target[i] < -1.0*pi) dihed_target[i]+=2.0*pi;
  }

  theta=(double *)gcemalloc(sizeof(double)*numdihedall);

  inputfile=efopen(inputfilename,"r");
  if (TLflag==ON)   getline(&line,&len,inputfile);

  if (initialstep>0) {
    for (i=0;i<initialstep;++i) {
      for (j=0;j<numatom;++j) {
	for (k=0;k<3;++k) {
	  fscanf(inputfile,"%lf",&crd[j*3+k]);
	}
      }
    }
  }

  logfile=efopen(logfilename,"w");

  for (i=initialstep;i<numstep;++i) {
    for (j=0;j<numatom;++j) {
      for (k=0;k<3;++k) {
	fscanf(inputfile,"%lf",&crd[j*3+k]);
      }
    }

    if ( (i%interval) == 0 || i==initialstep ) {
      flagout=ON;
      nd=0;
      if (nd < 0) {
	absnt=abs(nt)-1;
	for (k=0;k<numdihed[absnt];++k) {
	  for (m=0;m<4;++m) 
	    for (l=0;l<3;++l)
	      atom[m][l]=crd[(adpairs[absnt][k*4+m])*3+l];
	  
	  CDon(atom,pi,&theta[nd]);

	  if (theta[nd] < dihed_target[nd]-width || theta[nd] >= dihed_target[nd]+width) {
	    flagout=OFF;
	    break;
	  }
	  ++nd;
	}
	if (flagout==OFF) break;
      }
      else {
	for (j=0;j<nt;++j) {
	  for (k=0;k<numdihed[j];++k) {
	    for (m=0;m<4;++m) 
	      for (l=0;l<3;++l)
		atom[m][l]=crd[(adpairs[j][k*4+m])*3+l];
	  
	    CDon(atom,pi,&theta[nd]);

	    if (theta[nd] < dihed_target[nd]-width || theta[nd] >= dihed_target[nd]+width) {
	      flagout=OFF;
	      break;
	    }
	    ++nd;
	  }
	  if (flagout==OFF) break;
	}
      }

      if (flagout==ON) {
	nd=0;
	if (nt < 0) {
	  absnt=abs(nt)-1;
	  for (k=0;k<numdihed[absnt];++k) {
	    if (DEGflag==ON) {
	      fprintf(logfile,"%10.8lf ",theta[nd]*180.0/pi);
	    }
	    else {
	      fprintf(logfile,"%10.8lf ",theta[nd]);
	    }
	    ++nd;
	  }
	}
	else {	
	  for (j=0;j<nt;++j) {
	    for (k=0;k<numdihed[j];++k) {
	      if (DEGflag==ON) {
		fprintf(logfile,"%10.8lf ",theta[nd]*180.0/pi);
	      }
	      else {
		fprintf(logfile,"%10.8lf ",theta[nd]);
	      }
	      ++nd;
	    }
	  }
	  fprintf(logfile,"\n");
	}

	++numcrd;
	sprintf(crdfilename,"%s_%d",crdfilenamebase,numcrd);
	crdfile=efopen(crdfilename,"w");
	
	io_outputconf_Amberform(crdfile,numatom,crd);
	fclose(crdfile);
      }
    }
  }
  fclose(logfile);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-k interval] specify initial \n");
  printf("%s inputfilename outputfilename \n",progname);
}

double CDon(double atom[4][3], double pi, double *theta) {
  int i,j,k,l;
  int ii,jj,kk,ll;

  double m[3],n[3],m_n[3],n_n[3],lm,ln;
  double vij[3],vkj[3],vkl[3];
  double lkj;
  double vijvkj,vklvkj;

  double dihed;

  for (j=0;j<3;++j) {
    vij[j] = atom[1][j]-atom[0][j];
    vkj[j] = atom[1][j]-atom[2][j];
    vkl[j] = atom[3][j]-atom[2][j];
  }
  lkj=sqrt(inprod(vkj,vkj,3));
  
  outprod(vij,vkj,m);
  outprod(vkj,vkl,n);
  lm=sqrt(inprod(m,m,3));
  ln=sqrt(inprod(n,n,3));
  for (j=0;j<3;++j) {
    m_n[j]=m[j]/lm;
    n_n[j]=n[j]/ln;
  }
  
  dihed=inprod(m_n,n_n,3);
  if (dihed>=1.0)
    dihed=0.0;
  else if (dihed<=-1.0)
      dihed=pi;
  else
    dihed=acos(dihed);
  if (inprod(vij,n,3)>0) dihed=-dihed;
  if (dihed<-1.0*pi) dihed=2.0*pi+dihed;
  if (dihed>pi) dihed=-2.0*pi+dihed;
  
  *theta=dihed;
    
  return *theta;
}
