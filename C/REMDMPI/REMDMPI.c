
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "mpi.h"

#include "REMDMPI.h"

#include "PTL.h"
#include "EF.h"

#include "RAND.h"
#include "BOXMULL.h"
#include "MC.h"

#include "MDrun.h"

#include "netcdf_mineL.h"

#define BLANK 0
#define INPF  1
#define TEMP  2
#define TRJ   3
#define OUT   4

#define ON  0
#define OFF 1

int MPI_TREMD_pep_NHC_MP1998_Amber_AAFF(int myrank, int num_procs,int tag, MPI_Status* status,
					int numEX,  int numRE,
					double **crd,double **vel, double *mass, int numatom,
					double *zeta,double *V_zeta, double *Q,
					struct potential *e, struct force *f,
					double *T,double *NfKT, int numstep,int interval,
					double dt,double dt2,double wdt2[3],double wdt4[3], int nc,
					double *avePE, double *aveKE,double *aveT,
					double *varPE, double *varKE,double *varT, 
					double UNITT, double k_B, double tau, double pi,
					struct my_netcdf_out_id_MCD *nc_id_MCD,  FILE **outputfile ) {
  int i,j,k,c,l=0;
  int your_rank;

  double E_myrank,E_yours,beta_myrank,beta_yours,T_myrank,T_yours;
  double delta=0.0;
  double T_ym,T_my;
  double *crd_temp,*vel_temp;

  unsigned long init[4]={0x123, 0x234, 0x345, 0x456},length=4;
  init_by_array(init,length);

  crd_temp=(double *)gcemalloc(sizeof(double)*numatom*3);
  vel_temp=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numEX;++i) {
    E_myrank=runMD_NHC_MP1998_Amber_AAFF(crd[myrank],vel[myrank],mass,numatom,
					 &zeta[myrank],&V_zeta[myrank],Q[myrank],
					 e[myrank],f[myrank],
					 T[myrank],NfKT[myrank],numstep,interval,&l,
					 dt,dt2,wdt2,wdt4,nc,
					 &(avePE[myrank]),&(aveKE[myrank]),&(aveT[myrank]),
					 &(varPE[myrank]),&(varKE[myrank]),&(varT[myrank]),
					 UNITT,k_B,
					 nc_id_MCD[myrank],outputfile[myrank]);

    //    if(myrank==0) printf("Exchang now: %d !\n",i+1);

    if(myrank%2==0 ) {
      if ( i%2 == 0 ) your_rank=myrank+1;
      else {
	if ( myrank==0 )  your_rank=numRE-1;
	else  your_rank=myrank-1;
      }      

      MPI_Send((void *)&E_myrank, 1, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD);
    }
    else {
      if ( i%2 == 0 ) your_rank=myrank-1;
      else {
	if ( myrank==numRE-1 )  your_rank=0;
	else  your_rank=myrank+1;
      }      

      MPI_Recv((void *)&E_yours, 1, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD, status);

      T_yours=T[your_rank];
      T_myrank=T[myrank];

      beta_yours=1.0/(k_B*T_yours);
      beta_myrank=1.0/(k_B*T_myrank);
      delta=(beta_yours-beta_myrank)*(E_myrank-E_yours);

      if((c=Metropolis(delta))==1) {

	T_ym=sqrt(T_yours/T_myrank);
	T_my=sqrt(T_myrank/T_yours);

        for (j=0;j<numatom;++j) {
	  for (k=0;k<3;++k) {
	    vel[myrank][j*3+k]=T_ym*vel[myrank][j*3+k];
	    vel[your_rank][j*3+k]=T_my*vel[your_rank][j*3+k];
	  }
	}

        for (j=0;j<numatom;++j)  for (k=0;k<3;++k)   vel_temp[j*3+k]=vel[myrank][j*3+k];
        for (j=0;j<numatom;++j)  for (k=0;k<3;++k)   vel[myrank][j*3+k]=vel[your_rank][j*3+k];
        for (j=0;j<numatom;++j)  for (k=0;k<3;++k)   vel[your_rank][j*3+k]=vel_temp[j*3+k];
	
        for (j=0;j<numatom;++j)  for (k=0;k<3;++k)   crd_temp[j*3+k]=crd[myrank][j*3+k];
        for (j=0;j<numatom;++j)  for (k=0;k<3;++k)   crd[myrank][j*3+k]=crd[your_rank][j*3+k];
        for (j=0;j<numatom;++j)  for (k=0;k<3;++k)   crd[your_rank][j*3+k]=crd_temp[j*3+k];

	printf("%d -th EX acc between %d-%d : E(%d,T=%3.0lf)= %8.4lf E(%d,T=%3.0lf)= %8.4lf\n", 
	       i+1,myrank,your_rank,myrank,T[myrank],E_myrank,your_rank,T[your_rank],E_yours);
      }
      else  
	printf("%d -th EX rej between %d-%d : E(%d,T=%3.0lf)= %8.4lf E(%d,T=%3.0lf)=%8.4lf\n",
	       i+1,myrank,your_rank,myrank,T[myrank],E_myrank,your_rank,T[your_rank],E_yours);
    }

  }
  
  return 0;
}

void readInputs(FILE *inputfile, double **crd, double **vel, int numatom, double *T0) {
  int i,j,k;
  int c;
  int  d;
  double f1;
  char crdfilename[1000];
  FILE *crdfile;

  int State,Tflag=OFF;

  char *line;
  size_t len=0;

  i=0;  j=0;  f1=0.0;
  State=INPF;
  Tflag=OFF;

  while ((c=getc(inputfile))!=-1){
    if (State==INPF) {
      if (c==' ' || c=='\n') {
	if (j>0) {
	  crdfilename[j]=NULL;
	  crdfile=efopen(crdfilename,"r");
	  getline(&line,&len,crdfile);
	  fscanf(crdfile,"%d",&d);
	  for (k=0;k<numatom*3;++k) 
	    fscanf(crdfile,"%lf",&crd[i][k]);
	  fclose(crdfile);
	  State=TEMP;
	  j=0;
	}
      }
      else {
	crdfilename[j]=c;
	++j;
      }
    }
    else if (State==TEMP) {
      if (isdigit(c)) {
	d=(c-'0');
	f1=f1*10+(double)d;
	Tflag=ON;
      }
      else if (c==' ' || c=='\n') {
	if (Tflag==ON) {
	  T0[i]=(double)f1;
	  ++i;
	  f1=0.0;
	  State=INPF;
	  Tflag=OFF;
	}
      }
      else {
	printf("error\n");
	exit(1);
      }
    }
  }    
}

void readOutputs(FILE *outputfile, char *trjfilename[100], char *outputfilename[100]) {
  int i,j,k;
  int c;
  double f1;

  int State;

  i=0;  j=0;  f1=0.0;
  State=TRJ;

  while ((c=getc(outputfile))!=-1){
    if (State==TRJ) {
      if (c==' ' || c=='\n') {
	if (j>0) {
	  trjfilename[i][j]=NULL;
	  State=OUT;
	  j=0;
	}
      }
      else {
	trjfilename[i][j]=c;
	++j;
      }
    }
    else if (State==OUT) {
      if (c==' ' || c=='\n') {
	if (j>0) {
	  outputfilename[i][j]=NULL;
	  State=TRJ;
	  j=0;
	}
      }
      else {
	outputfilename[i][j]=c;
	++j;
      }
    }
  }    
}
