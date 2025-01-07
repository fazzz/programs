
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "MuSTARMD_ext_MPI.h"
#include "REMD_functions.h"

#define ON  0
#define OFF 1

double ***MPI_MuSTARMD_ext_pep_NHC_MP1998(int ND,int *myrank,
					  int num_procs, int tag, MPI_Status* status,
					  int *numRE, int numEX,
					  double **KZAlpha,
					  struct AAData_MuSTARMD *AlphaData,
					  struct TACCMData_MuSTARMD ZData,
					  struct AACGCommonData_MuSTARMD CData,
					  struct AmberParmL *ap_Alpha,
					  double *T0Alpha, double T0Z, int numstep, int interval, 
					  double dt,double dt2,double wdt2[3],double wdt4[3], int nc,
					  double UNITT, double k_B, double tau, double pi,
					  FILE* logfile ) {
  int i,j,k,c,l=0;
  int nd;
  int m,n;
  int your_rank;
  int count;

  double *betaAlpha/*ND*/,betaZ;

  double *EAlpha_m_Xi/*ND*/,*EAlpha_m_Xj/*ND*/,*EAlpha_n_Xi/*ND*/,*EAlpha_n_Xj/*ND*/;
  double EZ_m_Xi,EZ_m_Xj,EZ_n_Xi,EZ_n_Xj;

  double *KZAlpha_n/*ND*/;

  double **crdAlpha_receiv/*NDxnumatomx3*/;
  double *Z_receiv;
  double *crd_temp,*vel_temp;
  double *Z_temp,*vel_Z_temp;

  double *SData/*4*/,*RData/*4*/,*SKZ,*RKZ;
  double delta=0.0;

  int *exchange_data/*4*/,*exchange_data_receiv/*4*/;
  int **exchange_data_sum/*NDx(numRE*4)*/,**exchanged_data/*NDx(numRE*2)*/,**exchanged_data_receiv/*NDx(numRE*2)*/;

  double ***acc_ratio/*ND*numRE*numRE*/;

  unsigned long init[4]={0x123, 0x234, 0x345, 0x456},length=4;
  init_by_array(init,length);

  betaAlpha=(double *)gcemalloc(sizeof(double)*ND);
  for (i=0;i<ND;++i) betaAlpha[i]=1.0/(k_B*T0AA);
  betaZ=1.0/(k_B*T0Z);

  SData=(double *)gcemalloc(sizeof(double)*4);
  RData=(double *)gcemalloc(sizeof(double)*4);
  SKZ=(double *)gcemalloc(sizeof(double)*1);
  RKZ=(double *)gcemalloc(sizeof(double)*1);

  exchange_data=(double **)gcemalloc(sizeof(double *)*ND);
  exchange_data_receiv=(double **)gcemalloc(sizeof(double *)*ND);

  exchange_data_sum=(double **)gcemalloc(sizeof(double *)*ND);
  exchanged_data=(double **)gcemalloc(sizeof(double *)*ND);
  exchanged_data_receiv=(double **)gcemalloc(sizeof(double *)*ND);

  acc_ratio=(double ***)gcemalloc(sizeof(double **)*ND);

  for (nd=0;nd<ND;++nd) {
    exchange_data_sum=(int *)gcemalloc(sizeof(int)*numRE[nd]*4);
    exchanged_data=(int *)gcemalloc(sizeof(int)*numRE[nd]*2);
    exchanged_data_receiv=(int *)gcemalloc(sizeof(int)*numRE[nd]*2);

    acc_ratio[nd]=(double **)gcemalloc(sizeof(double *)*numRE[nd]);
    for (j=0;j<numRE[nd];++j) acc_ratio[nd]=(double *)gcemalloc(sizeof(double)*numRE[nd]);
    for (j=0;j<numRE[nd];++j) for (k=0;k<numRE[nd];++k) acc_ratio[nd][j][k]=0.0;
  }

  for (i=0;i<numEX;++i) {
    runMuSTARMD_AmberType_NHC_MP1998(ND,AlphaData,ZData,CData,
				     numstep,interval,&l,
				     dt,dt2,wdt2,wdt4,nc,
				     UNITT,k_B,pi,
				     &EAlpha_m_Xi,&EZm_Xi);

    fprintf(logfile,"%d-th: \n ",i+1);
    
    for (nd=0;nd<ND;++nd) {
      m=REMD_purmutation_func(myrank[nd]);

      fprintf(logfile,"m_%d= %d KZAlpha = %8.4lf\n",nd,m,ZData.KZAlpha[nd]);

      AlphaData[nd].zeta=0.0; AlphaData[nd].V_zeta=0.0;
      ZData.zetaZ=0.0; ZData.V_zetaZ=0.0;

      if(m%2==0 ) {
	if ( i%2 == 0 ) n=m+1;
	else 
	  if ( m==0 )  n=numRE[nd]-1;
	  else  n=m-1;
      }
      else {
	if ( i%2 == 0 ) n=m-1;
	else {
	  if ( m==numRE[nd]-1 )  n=0;
	  else  n=m+1;
	}      
      }
    
      your_rank=REMD_purmutation_inverse_func(n);

      SKZ=ZData.KZAlpha[nd];

      MPI_Send(SKZ, 1, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD);
      MPI_Recv(RKZ, 1, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD, status);

      KZAlpha_n=RKZ;

      MPI_Barrier(MPI_COMM_WORLD);

      CE_TACCMb_CGAA_test(AAData.crd,CGData.crd,ZData.Z,
			  (CData.numatom),ZData.numZ,
			  KZAAn,KZCGn,ZData.pairs,pi,
			  &EAAn_Xi,&ECGn_Xi,&EZn_Xi);

      SData[0]=EAlpha_m_Xi;
      SData[1]=EZm_Xi;

      SData[2]=EAlpha_n_Xi;
      SData[3]=EZn_Xi;

      MPI_Send(SData, 4, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD);
      MPI_Recv(RData, 4, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD, status);

      EAlpha_n_Xj=RData[0];
      EZn_Xj=RData[1];

      EAlpha_m_Xj=RData[2];
      EZm_Xj=RData[3];

      MPI_Barrier(MPI_COMM_WORLD);

      if (n>m) {

	delta=betaAlpha[nd]*((EAlpha_m_Xj-EAlpha_m_Xi)-(EAlpha_n_Xj-EAlpha_n_Xi))
  	     +betaZ*((EZm_Xj-EZm_Xi)-(EZn_Xj-EZn_Xi));

	MPI_Send(&delta, 1, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD);
      }
      else MPI_Recv(&delta, 1, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD, status);

      MPI_Barrier(MPI_COMM_WORLD);

      if((c=Metropolis(delta))==1) {
	if (myrank<your_rank)
	  printf("%d-th %d  Acc %d-%d : /_\\ = %8.4lf \n EAmi=%8.4lf EZmi=%8.4lf EAni=%8.4lf EZni=%8.4lf\n EAmj=%8.4lf EZmj=%8.4lf EAnj=%8.4lf EZnj=%8.4lf\n KZAm=%8.4lf KZAn=%8.4lf\n",
		 i+1,nd,myrank[nd],your_rank,myrank[nd],delta,
		 EAlpha_m_Xi,EZm_Xi,EAlpha_n_Xi,EZn_Xi,
		 EAlpha_m_Xj,EZm_Xj,EAlpha_n_Xj,EZn_Xj,
		 ZData.KZAlpha[nd],KZAlpha_n);
	ZData.KZAlpha[nd]=KZAlpha_n;

	if (myrank!=0) {
	  exchange_data[0]=myrank;
	  exchange_data[1]=your_rank;
	  exchange_data[2]=m;
	  exchange_data[3]=n;
	}
	else {
	  REMD_exchange_purmutation_funcs(myrank[nd],your_rank,m,n);
	  //	printf("176 m=%d n=%d\n",m,n);
	  acc_ratio[nd][m][n]+=1.0;
	}
      }
      else {
	if (myrank!=0) {
	  exchange_data[0]=-1;
	  exchange_data[1]=-1;
	  exchange_data[2]=-1;
	  exchange_data[3]=-1;
	}

	if (myrank<your_rank)
	  printf("%d-th %d Rej %d-%d : /_\\ = %8.4lf \n EAmi=%8.4lf EZmi=%8.4lf EAni=%8.4lf EZni=%8.4lf\n EAmj=%8.4lf EZmj=%8.4lf EAnj=%8.4lf EZnj=%8.4lf\n KZAm=%8.4lf KZAn=%8.4lf\n",
		 i+1,nd,myrank[nd],your_rank,myrank[nd],delta,
		 EAlpha_m_Xi,EZm_Xi,EAlpha_n_Xi,EZn_Xi,
		 EAlpha_m_Xj,EZm_Xj,EAlpha_n_Xj,EZn_Xj,
		 ZData.KZAlpha_nnd],KZAlpha_n);
      }
      
      if (myrank!=0) MPI_Send(&exchange_data, 4, MPI_INT, 0, tag, MPI_COMM_WORLD);

      MPI_Barrier(MPI_COMM_WORLD);

      if (myrank==0) {
	for (j=1;j<numRE[nd];++j) {
	  MPI_Recv(&exchange_data_receiv, 4, MPI_INT, j, tag, MPI_COMM_WORLD, status);
	  exchange_data_sum[nd][j*4+0]=exchange_data_receiv[0];
	  exchange_data_sum[nd][j*4+1]=exchange_data_receiv[1];
	  exchange_data_sum[nd][j*4+2]=exchange_data_receiv[2];
	  exchange_data_sum[nd][j*4+3]=exchange_data_receiv[3];
	  m=exchange_data_receiv[2];
	  n=exchange_data_receiv[3];
	  //	printf("%d:210 j=%d m=%d n=%d\n",myrank,j,m,n);
	  if ( m != -1 && n != -1 ) acc_ratio[nd][m][n]+=1.0;
	}
	for (j=1;j<numRE;++j) {
	  REMD_exchange_purmutation_funcs(exchange_data_sum[nd][j*4+0],exchange_data_sum[nd][j*4+1],
					  exchange_data_sum[nd][j*4+2],exchange_data_sum[nd][j*4+3]);
	}
	
	for (j=0;j<numRE[nd];++j) exchanged_data[nd][j]=index_parameters[j];
	for (j=0;j<numRE[nd];++j) exchanged_data[nd][j+numRE]=index_replicas[j];

	for (j=1;j<numRE[nd];++j) {
	  MPI_Send(exchanged_data[nd], numRE[nd]*2, MPI_INT, j, tag, MPI_COMM_WORLD);
	}
      }
      MPI_Barrier(MPI_COMM_WORLD);

      if (myrank!=0) {
	MPI_Recv(exchanged_data_receiv[nd], numRE[nd]*2, MPI_INT, 0, tag, MPI_COMM_WORLD, status);
	for (j=0;j<numRE[nd];++j) {
	  index_parameters[j]=exchanged_data_receiv[nd][j];
	  index_replicas[j]=exchanged_data_receiv[nd][numRE[nd]+j];
	}
      }
    
      MPI_Barrier(MPI_COMM_WORLD);
  }

  for (nd=0;nd<ND;++nd)
    for (i=0;i<numRE[nd];++i) 
      for (j=0;j<numRE[nd];++j) 
	acc_ratio[nd][i][j]=acc_ratio[nd][i][j]/numEX*2.0;
  
  return acc_ratio;
}

void MuSTARMD_ext_readInputs(FILE *inputfile,int ND, 
			     int numatom,int *numRE,int *myrank,
			     double **crdAlpha,double **velAlpha,
			     double **KZAlpha){
  int i,j,k;
  int c;
  int d;

  int nd=0;

  int AlphaINPflag=0;
  int AlphaKZ=1;

  double f1,fdummy;
  char crdfilename[1000];
  FILE *crdfile;

  int State,Tflag=OFF;

  char *line;
  size_t len=0;

  i=0;  j=0;  f1=0.0;
  State=AAINPF;
  Tflag=OFF;

  while ((c=getc(inputfile))!=-1){
    if (State==INPflag) {
      if (c==' ' || c=='\n') {
	if (j>0) {
	  crdfilename[j]='\0';
	  crdfile=efopen(crdfilename,"r");
	  getline(&line,&len,crdfile);
	  fscanf(crdfile,"%d",&d);
	  if (i==myrank)  for (k=0;k<numatom*3;++k) fscanf(crdfile,"%lf",&crdAlpha[nd][k]);
	  else for (k=0;k<numatom*3;++k) fscanf(crdfile,"%lf",&fdummy);
	  fclose(crdfile);
	  State=KZflag;
	  j=0;
	}
      }
      else {
	crdfilename[j]=c;
	++j;
      }
    }
    else {
      if (isdigit(c)) {
	d=(c-'0');
	f1=f1*10+(double)d;
	Tflag=ON;
      }
      else if (c==' ' || c=='\n') {
	if (Tflag==ON) {
	  KZAlpha[nd][i]=(double)f1;
	  ++i;
	  if (i==numRE) break;
	  f1=0.0;
	  State=INPFlag;
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
