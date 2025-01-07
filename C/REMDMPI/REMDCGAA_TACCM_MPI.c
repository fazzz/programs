
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "mpi.h"

#include "REMDMPI.h"
#include "REMD_TACCM_MPI.h"
#include "REMDCGAA_TACCM_MPI.h"

#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008.h"

#include "PTL.h"
#include "EF.h"

#include "RAND.h"
#include "BOXMULL.h"
#include "MC.h"

#include "MDrun.h"
#include "TACCM_MDrun.h"
#include "TACCM_MDrun_GOLMAA_PROTEINS2008.h"
#include "TACCM.h"

#include "netcdf_mineL.h"

#define ON  0
#define OFF 1

int MPI_CGAATREM_TACCM_MD_pep_NHC_MP1998_Amber_AAFF(int myrank, int num_procs,int tag, MPI_Status* status,int numEX,
						    ////////////// AA ////////////////////////////
						    double *crdAA,double *velAA, double *mass, int numatom,
						    double zetaAA,double V_zetaAA, double QAA,
						    struct potential e, struct force f, double NfKTAA,
						    double *avePEAA, double *aveKEAA,double *aveTAA,
						    double *varPEAA, double *varKEAA,double *varTAA, 
						    struct my_netcdf_out_id_MCD nc_id_MCDAA,  FILE *outputfileAA, 
						    ////////////// CG ////////////////////////////
						    double *crdCG,double *velCG,
						    double zetaCG,double V_zetaCG, double QCG,
						    struct potential_GOLMAA_PROTEINS2008 e_GOLM, double NfKTCG,
						    double *avePECG, double *aveKECG,double *aveTCG,
						    double *varPECG, double *varKECG,double *varTCG, 
						    struct my_netcdf_out_id_MCD nc_id_MCDCG,  FILE *outputfileCG, 
						    ////////////// COMMON /////////////////////////
						    double dt,double dt2,double wdt2[3],double wdt4[3], int nc,
						    int numstep, int interval, double TAA, double TCG,
						    double UNITT, double k_B, double tau, double pi,
						    //////////////// TACCM ///////////////////////
						    double **Z,double **velZ,double massZ,
						    double *zetaZ,double *V_zetaZ,
						    double T0Z,double QZ,double NfKTZ,int numZ,
						    double Kapa,int **pairs,
						    double *avePEZ, double *aveKEZ,double *aveTZ,
						    double *varPEZ, double *varKEZ,double *varTZ, 
						    FILE **trjfileZ, FILE **trjfilTheta
						    //////////////// TACCM ///////////////////////
						    ) {
  int i,j,k,c,l=0;
  int count;

  double beta;
  double E_CG_2_send,E_CG_2_reciev,E_AA_1,E_CG_1,E_AA_2_send,E_AA_2_reciev;
  double *SendData,*ReceivData;
  double *crdCG_receiv,*Z_CG_receiv;
  double delta=0.0;
  double *theta,*frcZ;
  double *Z_temp,*velZ_temp;

  unsigned long init[4]={0x123, 0x234, 0x345, 0x456},length=4;
  init_by_array(init,length);

  Z_temp=(double *)gcemalloc(sizeof(double)*2);
  velZ_temp=(double *)gcemalloc(sizeof(double)*2);
  theta=(double *)gcemalloc(sizeof(double)*numZ);
  frcZ=(double *)gcemalloc(sizeof(double)*numZ);

  SendData=(double *)gcemalloc(sizeof(double)*(1+numatom*3+numZ));
  ReceivData=(double *)gcemalloc(sizeof(double)*(1+numatom*3+numZ));

  crdCG_receiv=(double *)gcemalloc(sizeof(double)*numatom*3);
  Z_CG_receiv=(double *)gcemalloc(sizeof(double)*numZ);

  for (i=0;i<numEX;++i) {
    if(myrank==0 ) {
      E_AA_1=runTACCM_MD_NHC_MP1998_Amber_AAFF(crdAA,velAA,mass,numatom,&zetaAA,&V_zetaAA,QAA,
					       e,f,TAA,NfKTAA,numstep,interval,&l,
					       dt,dt2,wdt2,wdt4,nc,
					       avePEAA,aveKEAA,aveTAA,varPEAA,varKEAA,varTAA,
					       UNITT,k_B,nc_id_MCDAA,outputfileAA,
					       Z[myrank],velZ[myrank],massZ,&zetaZ[myrank],&V_zetaZ[myrank],
					       T0Z,QZ,NfKTZ,numZ,Kapa,pairs,pi,
					       &(avePEZ[myrank]),&(aveKEZ[myrank]),&(aveTZ[myrank]),
					       &(varPEZ[myrank]),&(varKEZ[myrank]),&(varTZ[myrank]), 
					       trjfileZ[myrank],trjfilTheta[myrank]);

      zetaAA=0.0; V_zetaAA=0.0;
      zetaZ[myrank]=0.0; V_zetaZ[myrank]=0.0;

    }
    else if(myrank==1 ) {
      /*E_CG_2_send=*/runTACCM_MD_NHC_MP1998_GOLMAA_PROTEINS2008(crdCG,velCG,mass,numatom,
								 &zetaCG,&V_zetaCG,QCG,
								 e,e_GOLM,TCG,NfKTCG,numstep,interval,&l,
								 dt,dt2,wdt2,wdt4,nc,
								 &avePECG,&aveKECG,&aveTCG,&varPECG,&varKECG,&varTCG,
								 UNITT,k_B,nc_id_MCDCG,outputfileCG,
								 Z[myrank],velZ[myrank],massZ,
								 &zetaZ[myrank],&V_zetaZ[myrank],
								 T0Z,QZ,NfKTZ,numZ,Kapa,pairs,pi,
								 &(avePEZ[myrank]),&(aveKEZ[myrank]),
								 &(aveTZ[myrank]),
								 &(varPEZ[myrank]),&(varKEZ[myrank]),
								 &(varTZ[myrank]),
								 trjfileZ[myrank],trjfilTheta[myrank],&E_CG_2_send);
      zetaCG=0.0; V_zetaCG=0.0;
      zetaZ[myrank]=0.0; V_zetaZ[myrank]=0.0;

      SendData[0]=E_CG_2_send;
      for (j=0;j<numatom*3;++j)	SendData[j+1]=crdCG[j];
      for (j=0;j<numZ;++j)	SendData[j+1+numatom*3]=Z[myrank][j];

      MPI_Send(SendData, 1+numatom*3+numZ, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if( myrank==0 ) {
      MPI_Recv(ReceivData, 1+numatom*3+numZ, MPI_DOUBLE, 1, tag, MPI_COMM_WORLD, status);

      E_CG_2_reciev=ReceivData[0];
      for (j=0;j<numatom*3;++j)	crdCG_receiv[j]=ReceivData[j+1];
      for (j=0;j<numZ;++j)	Z_CG_receiv[j]=ReceivData[j+1+numatom*3];

      //      printf("crdCG_receiv[10]=%8.4lf Z_CG_receiv[10]=%8.4lf\n",crdCG_receiv[10],Z_CG_receiv[10]);

      TACCM_CTheta(crdAA,numatom,theta,numZ,pairs,pi);
      E_AA_2_reciev=TACCM_calc_eff_FF_Z(Z_CG_receiv,numZ,theta,Kapa,frcZ,pi);

      TACCM_CTheta(crdCG_receiv,numatom,theta,numZ,pairs,pi);
      E_CG_1=TACCM_calc_eff_FF_Z(Z[0],numZ,theta,Kapa,frcZ,pi);

      printf("E_AA_1=%8.4lf E_AA_2=%8.4lf\n",E_AA_1,E_AA_2_reciev);
      printf("E_CG_2=%8.4lf E_CG_1=%8.4lf\n",E_CG_2_reciev,E_CG_1);

      beta=1.0/(k_B*T0Z);
      delta=beta*((E_AA_2_reciev-E_AA_1)-(E_CG_2_reciev-E_CG_1));
      printf("delta=%8.4lf\n",delta);

      if((c=Metropolis(delta))==1) {
	for (j=0;j<numZ;++j) {
      	  velZ_temp[j]=velZ[1][j];
      	  velZ[1][j]=velZ[0][j];
      	  velZ[0][j]=velZ_temp[j];
          
      	  Z_temp[j]=Z[1][j];
      	  Z[1][j]=Z[0][j];
      	  Z[0][j]=Z_temp[j];
      	}
    
      	printf("%d -th EX acc\n", i+1);
      }
      else printf("%d -th EX rej\n", i+1);
      
      zetaAA=0.0; V_zetaAA=0.0;
      zetaZ[0]=0.0; V_zetaZ[0]=0.0;
      zetaCG=0.0; V_zetaCG=0.0;
      zetaZ[1]=0.0; V_zetaZ[1]=0.0;
      
    }
  }
  
  return 0;
}
