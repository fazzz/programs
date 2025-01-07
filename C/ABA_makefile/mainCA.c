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

int main(int argc, char *argv[]) {
  int i,j,k,l,m,nt,dummy;
  int numatom,numstep,*numdihed;
  int flag='P',flagcn='n',flagof='f',amberflag=OFF;
  int *adpairs;
  int interval=1,pickinterval;

  double atom[4][3];
  double *protcoord;
  double theta,*old_theta;
  double pi;

  char *line;
  size_t len=0;

  char *inputfilename,*outputfilename,*parmfilename,*condfilename;
  FILE *inputfile,*outputfile,*logfile,*parmfile,*condfile;

  if (argc < 9) {
    printf("USAGE: ./%s Amberflag[at] cnflag[cn] offlag[of] numsnap interval inputfilename1(trj) inputfilename2(parm) outputfilename(atrj) \n",argv[0]);
    exit(1);
  }
  amberflag=*(++argv)[0];
  if (amberflag=='a')
    amberflag=ON;
  else
    amberflag=OFF;
  flagcn=*(++argv)[0];
  flagof=*(++argv)[0];
  numstep=atoi(*++argv);
  interval=atoi(*++argv);
  inputfilename=*++argv;
  parmfilename=*++argv;
  outputfilename=*++argv;
 
  parmfile=efopen(parmfilename,"r");
  readParmtop(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;
  protcoord=(double *)gcemalloc(sizeof(double)*numatom*3);
  old_theta=(double *)gcemalloc(sizeof(double)*(AP.NTHETH+AP.MTHETA));

  if (flag=='c') {

  }

  inputfile=efopen(inputfilename,"r");
  if (amberflag==ON)
    getline(&line,&len,inputfile);
  outputfile=efopen(outputfilename,"w");
  if (flagof=='f') {
    fprintf(outputfile,"numstep ");
    for (i=0;i<(AP.NTHETH+AP.MTHETA);++i)
      fprintf(outputfile,"%d ",i+1);
    fprintf(outputfile,"\n");
  }

  pi=acos(-1.0);
  for (i=0;i<numstep;++i) {
    if (flagof=='f') fprintf(outputfile,"%d ",i);
    io_scanconf(inputfile,numatom,protcoord,'x');
    
    if ( (i%interval) == 0) {
      for (j=0;j<AP.NTHETH;++j) {
	for (k=0;k<3;++k) 
	  for (l=0;l<3;++l)
	    atom[k][l]=protcoord[(AP.TH[j][k])-3+l];
      
	if (i==0) old_theta[j]=theta;
	if (flagcn=='n') {
	  theta=pick_angle(atom[0],atom[1],atom[2],0,0.0);
	  if (theta > pi)  theta = theta -2.0*pi;
	}
	else if (flagcn=='c') theta=pick_angle(atom[0],atom[1],atom[2],1,old_theta[j]);
	  
	fprintf(outputfile,"%e ",theta*180/acos(-1.0));
	old_theta[j]=theta;
      }
    }

    for (j=0;j<AP.MTHETA;++j) {
      for (k=0;k<3;++k) 
	for (l=0;l<3;++l)
	  atom[k][l]=protcoord[(AP.TA[j][k])-3+l];
      
      if (i==0) old_theta[j]=theta;
      if (flagcn=='n') {
	theta=pick_angle(atom[0],atom[1],atom[2],0,0.0);
	if (theta > pi)  theta = theta -2.0*pi;
      }
      else if (flagcn=='c') theta=pick_angle(atom[0],atom[1],atom[2],1,old_theta[j]);
	  
      fprintf(outputfile,"%e ",theta*180/acos(-1.0));
      old_theta[j]=theta;
    }
    fprintf(outputfile,"\n");
  }

  fclose(outputfile);
  fclose(inputfile);

  logfile=efopen("log_CA.txt","w");
  fprintf(logfile,"%d \n",AP.NTHETH+AP.MTHETA);
  for (j=0;j<AP.NTHETH;++j) {
    for (k=0;k<3;++k) fprintf(logfile,"%4d ",AP.TH[j][k]/3+1);
    //    for (k=0;k<4;++k) fprintf(logfile,"%s ",AP.IGRAPH[AP.TH[j][k]/3]);
    fprintf(logfile,"\n");
  }
  for (j=0;j<AP.MTHETA;++j) {
    for (k=0;k<3;++k) fprintf(logfile,"%4d ",AP.TA[j][k]/3+1);
    //    for (k=0;k<4;++k) fprintf(logfile,"%s ",AP.IGRAPH[AP.TA[j][k]/3]);
    fprintf(logfile,"\n");
  }
  fclose(logfile);
  
  return 0;
}

