
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MB.h"
#include "PT.h"
#include "EF.h"
#include "IO.h"

#define ON 1
#define OFF 0

int main(int argc, char *argv[]) {
  int i,j,k,numatom,numstep,numdihed,flag,flag2,flagcn,flagof,amberflag,dihedtype;
  int *atom_dihed_pair;
  int interval;
  double atom_i[3],atom_j[3],atom_k[3],atom_l[3];
  double *protcoord;
  double theta,*old_theta;
  double pi;
  char *line;
  size_t len=0;


  char *inputfilename,*outputfilename,*parmfilename,*condfilename;
  FILE *inputfile,*outputfile,*logfile,*parmfile,*condfile;

  amberflag=OFF;
  
  if (argc < 8) {
    printf("USAGE: ./calcdihed.exe flag(c or p or o or k)  flag2(c or n) flag3(c or t) inputfilename(coord) parmfile inputfile(cond) outputfilename(dihed) option(flag i) \n");
    exit(1);
  }
  else {
    flag=(*++argv)[0];
    if (flag != 'c' && flag != 'p' && flag != 'o' && flag != 'k' ) {
      printf("flag error: must be c or p or o or k");
      exit(1);
    }
    flagcn=(*++argv)[0];
    if (flagcn != 'c' && flagcn != 'n' ) {
      printf("flag error: must be c or v");
      exit(1);
    }
    flagof=(*++argv)[0];
    if (flagof != 'c' && flagof != 't' && flagof != 'a' && flagof != 'b') {
      printf("flag error: must be c or t");
      exit(1);
    }
    if (flagof == 'a' ) {
      flagof='c';
      amberflag=ON;
    }
    else if (flagof == 'b') {
      flagof='t';
      amberflag=ON;
    }
    inputfilename  = *++argv;
    parmfilename   = *++argv;
    condfilename   = *++argv;
    outputfilename = *++argv;
  }

  parmfile=efopen(parmfilename,"r");
  readParmtop(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;
  protcoord=(double *)emalloc(sizeof(double)*numatom*3);

  condfile=efopen(condfilename,"r");
  fscanf(condfile,"%d",&numstep);
  if (flag == 'c') {
    fscanf(condfile,"%d",&numdihed);
  }
  else {
    numdihed=AP.NPHIH+AP.MPHIA;
  }
  atom_dihed_pair=(int *)emalloc(sizeof(int)*4*numdihed);
  if (flag == 'c') {
    for (i=0;i<numdihed;++i) {
      fscanf(condfile,"%d",&atom_dihed_pair[i*4+0]);
      fscanf(condfile,"%d",&atom_dihed_pair[i*4+1]);
      fscanf(condfile,"%d",&atom_dihed_pair[i*4+2]);
      fscanf(condfile,"%d",&atom_dihed_pair[i*4+3]);
    }
  }
  if (argc >= 9)
    fscanf(condfile,"%d",&interval);
  else
    interval=1;
  fclose(condfile);
 
  /*****************************************************************/
  /* if (flag != 'c') {						   */
  /*   atom_dihed_pair=emalloc(sizeof(int)*(AP.NPHIH+AP.MPHIA)*4); */
  /* }								   */
  /*****************************************************************/
  if (flag != 'c') {						    
  numdihed=0;
  for (i=0;i<AP.MPHIA;++i) {
    if (AP.PA[i][0] >=0 && AP.PA[i][1] >= 0 && AP.PA[i][2] >= 0 && AP.PA[i][3] >=0) {
      dihedtype=check_phi_psi_omega_kai1(abs(AP.PA[i][0])/3,abs(AP.PA[i][1])/3,abs(AP.PA[i][2])/3,abs(AP.PA[i][3])/3);
      if (dihedtype==phi || dihedtype==psi) {
	if (flag == 'p' || flag == 'o' || flag == 'k') {
	  for (j=0;j<4;++j) {
	    atom_dihed_pair[numdihed*4+j]=abs(AP.PA[i][j])/3+1;
	  }
	  ++numdihed;
	}	
      }
      else if (dihedtype==omega) {
	if (flag == 'o' || flag == 'k') {
	  for (j=0;j<4;++j) {
	    atom_dihed_pair[numdihed*4+j]=abs(AP.PA[i][j])/3+1;
	  }
	  ++numdihed;
	}	
      }
    }
  }
  for (i=0;i<AP.NPHIH;++i) {
    if (AP.PH[i][0] >=0 && AP.PH[i][1] >= 0 && AP.PH[i][2] >= 0 && AP.PH[i][3] >=0) {
      dihedtype=check_phi_psi_omega_kai1(abs(AP.PH[i][0])/3,abs(AP.PH[i][1])/3,abs(AP.PH[i][2])/3,abs(AP.PH[i][3])/3);
      if (dihedtype==phi || dihedtype==psi) {
	if (flag == 'p' || flag == 'o' || flag == 'k') {
	  for (j=0;j<4;++j) {
	    atom_dihed_pair[numdihed*4+j]=abs(AP.PH[i][j])/3+1;
	  }
	  ++numdihed;
	}	
      }
      else if (dihedtype==omega) {
	if (flag == 'o' || flag == 'k') {
	  for (j=0;j<4;++j) {
	    atom_dihed_pair[numdihed*4+j]=abs(AP.PH[i][j])/3+1;
	  }
	  ++numdihed;
	}	
      }
    }
  }
  if (flag == 'k') {
    for (i=0;i<AP.NPHIH;++i) {
      if (AP.PH[i][0] >=0 && AP.PH[i][1] >= 0 && AP.PH[i][2] >= 0 && AP.PH[i][3] >=0) {
	dihedtype=check_phi_psi_omega_kai1(abs(AP.PH[i][0])/3,abs(AP.PH[i][1])/3,abs(AP.PH[i][2])/3,abs(AP.PH[i][3])/3);
	if (dihedtype==kai1 ) {
	  for (j=0;j<4;++j) {
	    atom_dihed_pair[numdihed*4+j]=abs(AP.PH[i][j])/3+1;
	  }
	  ++numdihed;
	}	
      }
    }
  }
  }

  old_theta=emalloc(sizeof(double)*numdihed);
  inputfile=efopen(inputfilename,"r");
  if (amberflag==ON)
    getline(&line,&len,inputfile);

  outputfile=efopen(outputfilename,"w");
  if (flagof=='c') {
    fprintf(outputfile,"numstep ");
    for (i=0;i<numdihed;++i)
      fprintf(outputfile,"%d ",i);
    fprintf(outputfile,"\n");
  }

  pi=acos(-1.0);
  for (i=0;i<numstep;++i) {
    if (flagof=='c')
      fprintf(outputfile,"%d ",i);
    io_scanconf(inputfile,numatom,protcoord,'x');
    if ( (i%interval) == 0) {
      for (j=0;j<numdihed;++j) {
	for (k=0;k<3;++k) {
	  atom_i[k]=protcoord[(atom_dihed_pair[j*4]-1)*3+k];
	  atom_j[k]=protcoord[(atom_dihed_pair[j*4+1]-1)*3+k];
	  atom_k[k]=protcoord[(atom_dihed_pair[j*4+2]-1)*3+k];
	  atom_l[k]=protcoord[(atom_dihed_pair[j*4+3]-1)*3+k];
	}

	if (i==0) {
	  old_theta[j]=theta;
	}
	if (flagcn=='n') {
	  theta=pick_dihed(atom_i,atom_j,atom_k,atom_l,0,old_theta[j]);
	  if (theta > pi)
	    theta = theta -2.0*pi;
	}
	else if (flagcn=='c') {
	  theta=pick_dihed(atom_i,atom_j,atom_k,atom_l,1,old_theta[j]);
	}

	fprintf(outputfile,"%e ",theta*180/acos(-1.0));
	old_theta[j]=theta;
      }
      fprintf(outputfile,"\n");
    }
  }
  fclose(outputfile);

  logfile=efopen("log.txt","w");
  fprintf(logfile,"%d\n",numdihed);
  for (i=0;i<numdihed;++i) {
    for (j=0;j<4;++j) {
      fprintf(logfile,"%d ",atom_dihed_pair[i*4+j]);
    }
    for (j=0;j<4;++j) {
      fprintf(logfile,"%s ",AP.IGRAPH[atom_dihed_pair[i*4+j]-1]);
    }
    fprintf(logfile,"\n");
  }
  fclose(logfile);

  return 0;
}

