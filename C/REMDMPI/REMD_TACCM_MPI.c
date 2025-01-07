
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "mpi.h"

#include "REMDMPI.h"
#include "REMD_TACCM_MPI.h"

#include "PTL.h"
#include "EF.h"

#include "RAND.h"
#include "BOXMULL.h"
#include "MC.h"

#include "MDrun.h"
#include "TACCM_MDrun.h"

#include "netcdf_mineL.h"

#define ON  0
#define OFF 1

int MPI_TREM_TACCM_MD_pep_NHC_MP1998_Amber_AAFF(int myrank, int num_procs,int tag, MPI_Status* status,
						int numEX,  int numRE,
						double **crd,double **vel, double *mass, int numatom,
						double *zeta,double *V_zeta, double Q,
						struct potential *e, struct force *f,
						double T,double NfKT, int numstep,int interval,
						double dt,double dt2,double wdt2[3],double wdt4[3], int nc,
						double *avePE, double *aveKE,double *aveT,
						double *varPE, double *varKE,double *varT, 
						double UNITT, double k_B, double tau, double pi,
						struct my_netcdf_out_id_MCD *nc_id_MCD,  FILE **outputfile, 
						//////////////// TACCM ///////////////////////
						double **Z,double **velZ,double massZ,
						double *zetaZ,double *V_zetaZ,
						double *TZ,double *QZ,double *NfKTZ,int numZ,
						double Kapa,int **pairs,
						double *avePEZ, double *aveKEZ,double *aveTZ,
						double *varPEZ, double *varKEZ,double *varTZ, 
						FILE **trjfileZ, FILE **trjfilTheta
						//////////////// TACCM ///////////////////////
						) {
  int i,j,k,c,l=0;
  int your_rank;

  double E_myrank,E_yours,beta_myrank,beta_yours,T_myrank,T_yours;
  double delta=0.0;
  double T_ym,T_my;
  double *Z_temp,*velZ_temp;

  unsigned long init[4]={0x123, 0x234, 0x345, 0x456},length=4;
  init_by_array(init,length);

  Z_temp=(double *)gcemalloc(sizeof(double)*numRE);
  velZ_temp=(double *)gcemalloc(sizeof(double)*numRE);

  //  printf("yes this is line 65 in MPI%d\n",myrank);

  /******************************************************************************/
  /* for (j=0;j<numZ;++j) 						        */
  /*   printf("%d %d %d %d\n",pairs[j][0],pairs[j][1],pairs[j][2],pairs[j][3]); */
  /******************************************************************************/

  for (i=0;i<numEX;++i) {
    E_myrank=runTACCM_MD_NHC_MP1998_Amber_AAFF(crd[myrank],vel[myrank],mass,numatom,
					       &zeta[myrank],&V_zeta[myrank],Q,
					       e[myrank],f[myrank],
					       T,NfKT,numstep,interval,&l,
					       dt,dt2,wdt2,wdt4,nc,
					       &(avePE[myrank]),&(aveKE[myrank]),&(aveT[myrank]),
					       &(varPE[myrank]),&(varKE[myrank]),&(varT[myrank]),
					       UNITT,k_B,
					       nc_id_MCD[myrank],outputfile[myrank],
					       Z[myrank],velZ[myrank],massZ,
					       &zetaZ[myrank],&V_zetaZ[myrank],
					       TZ[myrank],QZ[myrank],NfKTZ[myrank],numZ,
					       Kapa,pairs,pi,
					       &(avePEZ[myrank]),&(aveKEZ[myrank]),&(aveTZ[myrank]),
					       &(varPEZ[myrank]),&(varKEZ[myrank]),&(varTZ[myrank]), 
					       trjfileZ[myrank], trjfilTheta[myrank]);

    //    if(myrank==0) printf("Exchang now: %d !\n",i+1);

    if(myrank%2==0 ) {
      if ( i%2 == 0 ) your_rank=myrank+1;
      else {
	if ( myrank==0 )  your_rank=numRE-1;
	else  your_rank=myrank-1;
      }      

      MPI_Send((void *)&E_myrank, 1, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD);

      zeta[myrank]=0.0;  V_zeta[myrank]=0.0;
      zetaZ[myrank]=0.0; V_zetaZ[myrank]=0.0;

    }
    else {
      if ( i%2 == 0 ) your_rank=myrank-1;
      else {
	if ( myrank==numRE-1 )  your_rank=0;
	else  your_rank=myrank+1;
      }      

      MPI_Recv((void *)&E_yours, 1, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD, status);

      T_yours=TZ[your_rank];
      T_myrank=TZ[myrank];

      beta_yours=1.0/(k_B*T_yours);
      beta_myrank=1.0/(k_B*T_myrank);
      delta=(beta_yours-beta_myrank)*(E_myrank-E_yours);

      if((c=Metropolis(delta))==1) {

	T_ym=sqrt(T_yours/T_myrank);
	T_my=sqrt(T_myrank/T_yours);

        for (j=0;j<numZ;++j) {
	  velZ[myrank][j]=T_ym*velZ[myrank][j];
	  velZ[your_rank][j]=T_my*velZ[your_rank][j];
	}

	printf("B Exc\n");
	for (j=0;j<numZ;++j) printf("Z[%d][%d]=%lf Z[%d][%d]=%lf \n",
				    myrank,j,Z[myrank][j],your_rank,j,Z[your_rank][j]);

        for (j=0;j<numZ;++j)  velZ_temp[j]=velZ[myrank][j];
        for (j=0;j<numZ;++j)  velZ[myrank][j]=velZ[your_rank][j];
        for (j=0;j<numZ;++j)  velZ[your_rank][j]=velZ_temp[j];
	
        for (j=0;j<numZ;++j)  Z_temp[j]=Z[myrank][j];
        for (j=0;j<numZ;++j)  Z[myrank][j]=Z[your_rank][j];
        for (j=0;j<numZ;++j)  Z[your_rank][j]=Z_temp[j];

	//	zeta[myrank]=0.0;  V_zeta[myrank]=0.0;
	//	zetaZ[myrank]=0.0; V_zetaZ[myrank]=0.0;

	printf("A Exc\n");
	for (j=0;j<numZ;++j) printf("Z[%d][%d]=%lf Z[%d][%d]=%lf \n",
				    myrank,j,Z[myrank][j],your_rank,j,Z[your_rank][j]);

	printf("%d -th EX acc between %d-%d : E(%d,T=%3.0lf)= %8.4lf E(%d,T=%3.0lf)= %8.4lf\n", 
	       i+1,myrank,your_rank,myrank,TZ[myrank],E_myrank,your_rank,TZ[your_rank],E_yours);
      }
      else  
	printf("%d -th EX rej between %d-%d : E(%d,T=%3.0lf)= %8.4lf E(%d,T=%3.0lf)=%8.4lf\n",
	       i+1,myrank,your_rank,myrank,TZ[myrank],E_myrank,your_rank,TZ[your_rank],E_yours);
    }
    zeta[myrank]=0.0;  V_zeta[myrank]=0.0; 
    zetaZ[myrank]=0.0; V_zetaZ[myrank]=0.0;
  }
  
  return 0;
}
