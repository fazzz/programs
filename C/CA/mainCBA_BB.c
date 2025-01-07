#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "netcdf_mineL.h"

#include "TOPO.h"
#include "PTL.h"
#include "EF.h"

#define ON 1
#define OFF 0

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d,m,nt,dummy,num;
  int numatom,numstep;
  int amberflag=OFF;

  int numB,numA;
  int **pairsA,**pairsB;
  int interval=1;

  double atom[3][3];

  double theta,length;
  double pi;

  char *line;
  size_t len1=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double crd_nc[MAXATOM][3],*crd;
  struct my_netcdf_out_id_MCD nc_id_MCD;
  struct my_netcdf_out_id_AMBER nc_id_MD;

  char *inputfilename,*outputfilenameA,*outputfilenameB,*parmtopfilename,*progname;
  char *logfilenameA="log_CBA_BB_A.txt",*logfilenameB="log_CBA_BB_B.txt";
  FILE *inputfile,*outputfileA,*outputfileB,*parmfile,*logfileA,*logfileB;

  pi=acos(-1.0);

  while((c=getopt(argc,argv,"hAl:L:"))!=-1) {
    switch(c) {
    case 'A':
      amberflag=ON;
      break;
    case 'l':
      logfilenameA=optarg;
      break;
    case 'L':
      logfilenameB=optarg;
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

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  inputfilename   = *argv;
  parmtopfilename = *++argv;
  outputfilenameB = *++argv;
  outputfilenameA = *++argv;

  if ( amberflag==OFF ) numstep=myncL_get_present_step_MCD(inputfilename,&nc_id_MCD);
  else  numstep=myncL_get_present_step_AMBER(inputfilename,&nc_id_MD);

  parmfile=efopen(parmtopfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;

  pairsA=(int **)gcemalloc(sizeof(int *)*1);
  pairsA[0]=(int *)gcemalloc(sizeof(int)*2);
  pairsB=(int **)gcemalloc(sizeof(int *)*1);
  pairsB[0]=(int *)gcemalloc(sizeof(int)*3);

  logfileB=efopen(logfilenameB,"w");
  numB=0;
  for (i=0;i<AP.MBONA;++i) {
    j=abs(AP.BA[i][0])/3;
    k=abs(AP.BA[i][1])/3;
    if ( strncmp(AP.ITREE[i],"M",1)==0 && strncmp(AP.ITREE[j],"M",1)==0)  {
      pairsB=(int **)gcerealloc(pairsB,sizeof(int *)*(numB+1));
      pairsB[numB]=(int *)gcemalloc(sizeof(int)*2);
      pairsB[numB][0]=j;
      pairsB[numB][1]=k;
      ++numB;
      fprintf(logfileB,"%4d-%4d %4s-%4s\n",j,k,AP.IGRAPH[j],AP.IGRAPH[j]);
    }
  }
  fclose(logfileB);

  logfileA=efopen(logfilenameA,"w");
  numA=0;
  for (i=0;i<AP.MTHETA;++i) {
    j=abs(AP.TA[i][0])/3;
    k=abs(AP.TA[i][1])/3;
    l=abs(AP.TA[i][2])/3;
    if ( strncmp(AP.ITREE[j],"M",1)==0 && strncmp(AP.ITREE[k],"M",1)==0 && strncmp(AP.ITREE[l],"M",1)==0 ) {
      pairsA=(int **)gcerealloc(pairsA,sizeof(int *)*(numA+1));
      pairsA[numA]=(int *)gcemalloc(sizeof(int)*3);
      pairsA[numA][0]=j;
      pairsA[numA][1]=k;
      pairsA[numA][2]=l;
      ++numA;
      fprintf(logfileA,"%4d-%4d-%4d %4s-%4s-%4s\n",j,k,l,AP.IGRAPH[j],AP.IGRAPH[j]);
    }
  }
  fclose(logfileA);
  
  outputfileB=efopen(outputfilenameB,"w");
  outputfileA=efopen(outputfilenameA,"w");
  for (i=0;i<numstep;++i) {
    if ( amberflag==OFF ) myncL_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id_MCD,crd_nc);
    else myncL_open_inq_get_sh_AMBER(inputfilename,numatom,i,1,i+1,&nc_id_MD,crd_nc);
    
    if ( (i%interval) == 0 ) {
      for (j=0;j<numA;++j) {
	for (k=0;k<3;++k) 
	  for (l=0;l<3;++l)
	    atom[k][l]=crd_nc[(pairsA[j][k])][l];
      
	theta=ang(atom[0],atom[1],atom[2]);
	if (theta < -1.0*pi) theta+=2.0*pi;
	else if (theta > pi) theta-=2.0*pi;

	fprintf(outputfileA,"%e ",theta*180.0/pi);
      }
      fprintf(outputfileA,"\n");
      
      for (j=0;j<numB;++j) {
	for (k=0;k<2;++k) 
	  for (l=0;l<3;++l)
	    atom[k][l]=crd_nc[(pairsB[j][k])][l];
      
	length=len(atom[0],atom[1]);
	fprintf(outputfileB,"%e ",length);
	  
      }
      fprintf(outputfileB,"\n");
    }
  }

  fclose(outputfileB);
  fclose(outputfileA);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-A] amberflag \n");
  printf("[-h ] help \n");
  printf("%s inputfilename parmfilename outputfilenameB outputfilenameA\n",progname);
}

 
