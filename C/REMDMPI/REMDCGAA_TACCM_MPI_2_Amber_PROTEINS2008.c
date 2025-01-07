
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "REMDCGAA_TACCM_MPI_2_Amber_PROTEINS2008.h"
#include "REMD_functions.h"

#define ON  0
#define OFF 1

double **MPI_CGAATREM_TACCM_MD_pep_NHC_MP1998_Amber_PROTEINS2008(int myrank,int num_procs,
								 int tag, MPI_Status* status,
								 int numRE, int numEX, double *KZAA, double *KZCG,
								 struct AADataforREMD_Amber AAData,
								 struct CGDataforREMD_PROTEINS2008 CGData,
								 struct TACCMDataforREMD_A_P2008 ZData,
								 struct AACGCommonDataforREMD_A_P2008 CData,
								 double T0AA,double T0CG, double T0Z, 
								 int numstep, int interval, 
								 double dt,double dt2,
								 double wdt2[3],double wdt4[3], int nc,
								 double UNITT, double k_B, double tau, double pi,
								 FILE* logfile ) {
  int i,j,k,c,l=0;
  int m,n;
  int your_rank;
  int count;

  double betaAA,betaCG,betaZ;

  double EAAm_Xi,EAAm_Xj,EAAn_Xi,EAAn_Xj;
  double ECGm_Xi,ECGm_Xj,ECGn_Xi,ECGn_Xj;
  double EZm_Xi,EZm_Xj,EZn_Xi,EZn_Xj;

  double KZAAn,KZCGn;

  double *crdAA_receiv,*crdCG_receiv;
  double *Z_receiv;
  double *crd_temp,*vel_temp;
  double *Z_temp,*vel_Z_temp;

  double *SData,*RData,*SKZ,*RKZ;
  double delta=0.0;

  int exchange_data[4],exchange_data_receiv[4],*exchange_data_sum,*exchanged_data,*exchanged_data_receiv;
  double **acc_ratio;

  unsigned long init[4]={0x123, 0x234, 0x345, 0x456},length=4;
  init_by_array(init,length);

  betaAA=1.0/(k_B*T0AA);
  betaCG=1.0/(k_B*T0CG);
  betaZ=1.0/(k_B*T0Z);

  SData=(double *)gcemalloc(sizeof(double)*6/*8*/);
  RData=(double *)gcemalloc(sizeof(double)*6/*8*/);
  SKZ=(double *)gcemalloc(sizeof(double)*2);
  RKZ=(double *)gcemalloc(sizeof(double)*2);

  exchange_data_sum=(int *)gcemalloc(sizeof(int)*numRE*4);
  exchanged_data=(int *)gcemalloc(sizeof(int)*numRE*2);
  exchanged_data_receiv=(int *)gcemalloc(sizeof(int)*numRE*2);

  acc_ratio=(double **)gcemalloc(sizeof(double *)*numRE);
  for (i=0;i<numRE;++i) acc_ratio[i]=(double *)gcemalloc(sizeof(double)*numRE);
  for (i=0;i<numRE;++i) for (j=0;j<numRE;++j) acc_ratio[i][j]=0.0;

  for (i=0;i<numEX;++i) {
    runTACCM_CGAA_MD_NHC_MP1998_Amber_PROTEINS2008(// AA /////////////////////////////////////////////////////////
						   AAData.crd,AAData.vel,
						   &(AAData.zeta),&(AAData.V_zeta),AAData.Q,
						   AAData.e,AAData.f,AAData.T,AAData.NfKT,
						   AAData.avePE,AAData.aveKE,AAData.aveT,
						   AAData.varPE,AAData.varKE,AAData.varT,
						   AAData.nc_id_MCD,AAData.outputfile,
						   // CG /////////////////////////////////////////////////////////
						   CGData.crd,CGData.vel,
						   &(CGData.zeta),&(CGData.V_zeta),CGData.Q,
						   CGData.e,
						   CGData.T,CGData.NfKT,
						   CGData.avePE,CGData.aveKE,CGData.aveT,
						   CGData.varPE,CGData.varKE,CGData.varT,
						   CGData.nc_id_MCD,CGData.outputfile,
						   // Z  /////////////////////////////////////////////////////////
						   ZData.Z,ZData.velZ,ZData.massZ,
						   &(ZData.zetaZ),&(ZData.V_zetaZ),
						   ZData.QZ,ZData.NfKTZ,ZData.T,
						   ZData.numZ,ZData.KZAA,ZData.KZCG,ZData.pairs,
						   ZData.avePEZ,ZData.aveKEZ,ZData.aveTZ,
						   ZData.varPEZ,ZData.varKEZ,ZData.varTZ, 
						   ZData.trjfileZ,ZData.trjfilThetaAA,ZData.trjfilThetaCG,
						   // CM  /////////////////////////////////////////////////////////
						   CData.mass,(CData.numatom),(CData.numheavyatom),
						   numstep,interval,&l,
						   dt,dt2,wdt2,wdt4,nc,UNITT,k_B,pi,
						   &EAAm_Xi,&ECGm_Xi,&EZm_Xi);

    m=REMD_purmutation_func(myrank);
    fprintf(logfile,"%d-th: m= %d KZAA= %8.4lf KZCG= %8.4lf\n",i+1,m,ZData.KZAA,ZData.KZCG);

    AAData.zeta=0.0; AAData.V_zeta=0.0;
    CGData.zeta=0.0; CGData.V_zeta=0.0;
    ZData.zetaZ=0.0; ZData.V_zetaZ=0.0;

    if(m%2==0 ) {
      if ( i%2 == 0 ) n=m+1;
      else {
	if ( m==0 )  n=numRE-1;
	else  n=m-1;
      }
    }      
    else {
      if ( i%2 == 0 ) n=m-1;
      else {
	if ( m==numRE-1 )  n=0;
	else  n=m+1;
      }      
    }

    your_rank=REMD_purmutation_inverse_func(n);

    SKZ[0]=ZData.KZAA;
    SKZ[1]=ZData.KZCG;
    //    printf("107\n");
    MPI_Send(SKZ, 2, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD);
    MPI_Recv(RKZ, 2, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD, status);
    //    printf("110\n");
    KZAAn=RKZ[0];
    KZCGn=RKZ[1];

    MPI_Barrier(MPI_COMM_WORLD);

    CE_TACCM_CGAA(AAData.crd,CGData.crd,ZData.Z,
		  (CData.numatom),ZData.numZ,
		  KZAAn,KZCGn,ZData.pairs,pi,
		  &EAAn_Xi,&ECGn_Xi,&EZn_Xi);

    SData[0]=EAAm_Xi;
    SData[1]=ECGm_Xi;
    SData[2]=EZm_Xi;

    SData[3]=EAAn_Xi;
    SData[4]=ECGn_Xi;
    SData[5]=EZn_Xi;

    //    SData[6]=ZData.KZAA;
    //    SData[7]=ZData.KZCG;
    //    printf("148\n");
    MPI_Send(SData, /*8*/6, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD);
    MPI_Recv(RData, /*8*/6, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD, status);
    //    printf("151\n");
    EAAn_Xj=RData[0];
    ECGn_Xj=RData[1];
    EZn_Xj=RData[2];

    EAAm_Xj=RData[3];
    ECGm_Xj=RData[4];
    EZm_Xj=RData[5];

    //    KZAAn=RData[6];
    //    KZCGn=RData[7];

    MPI_Barrier(MPI_COMM_WORLD);

    if (n>m) {
      //      delta=betaAA*((EAAm_Xj-EAAm_Xi)-(EAAn_Xi-EAAn_Xj))
      //      	+betaCG*((ECGm_Xj-ECGm_Xi)-(ECGn_Xi-ECGn_Xj))
      //      	+betaZ*((EZm_Xj-EZm_Xi)-(EZn_Xi-EZn_Xj));

      delta=betaAA*((EAAm_Xj-EAAm_Xi)-(EAAn_Xj-EAAn_Xi))
      	+betaCG*((ECGm_Xj-ECGm_Xi)-(ECGn_Xj-ECGn_Xi))
      	+betaZ*((EZm_Xj-EZm_Xi)-(EZn_Xj-EZn_Xi));

      /**************************************************************/
      /* delta=-betaAA*((EAAm_Xj-EAAm_Xi)-(EAAn_Xi-EAAn_Xj))	    */
      /* 	-betaCG*((ECGm_Xj-ECGm_Xi)-(ECGn_Xi-ECGn_Xj))	    */
      /* 	-betaZ*((EZm_Xj-EZm_Xi)-(EZn_Xi-EZn_Xj));	    */
      /**************************************************************/

      MPI_Send(&delta, 1, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD);
    }
    else MPI_Recv(&delta, 1, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD, status);

    MPI_Barrier(MPI_COMM_WORLD);

    if((c=Metropolis(delta))==1) {
      if (myrank<your_rank)
	printf("%d-th Acc %d-%d : /_\\ = %8.4lf \n EAAmi=%8.4lf ECGmi=%8.4lf EZmi=%8.4lf EAAni=%8.4lf ECGni=%8.4lf EZni=%8.4lf\n EAAmj=%8.4lf ECGmj=%8.4lf EZmj=%8.4lf EAAnj=%8.4lf ECGnj=%8.4lf EZnj=%8.4lf\n KZAAm=%8.4lf KZAAn=%8.4lf KZCGm=%8.4lf KZCGn=%8.4lf\n",
	       i+1,myrank,your_rank,myrank,delta,
	       EAAm_Xi,ECGm_Xi,EZm_Xi,EAAn_Xi,ECGn_Xi,EZn_Xi,
	       EAAm_Xj,ECGm_Xj,EZm_Xj,EAAn_Xj,ECGn_Xj,EZn_Xj,
	       ZData.KZAA,KZAAn,ZData.KZCG,KZCGn);
      ZData.KZAA=KZAAn;
      ZData.KZCG=KZCGn;

      if (myrank!=0) {
	exchange_data[0]=myrank;
	exchange_data[1]=your_rank;
	exchange_data[2]=m;
	exchange_data[3]=n;
      }
      else {
	REMD_exchange_purmutation_funcs(myrank,your_rank,m,n);
	//	printf("176 m=%d n=%d\n",m,n);
	acc_ratio[m][n]+=1.0;
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
	printf("%d-th Rej %d-%d : /_\\ = %8.4lf \n EAAmi=%8.4lf ECGmi=%8.4lf EZmi=%8.4lf EAAni=%8.4lf ECGni=%8.4lf EZni=%8.4lf\n EAAmj=%8.4lf ECGmj=%8.4lf EZmj=%8.4lf EAAnj=%8.4lf ECGnj=%8.4lf EZnj=%8.4lf\n KZAAm=%8.4lf KZAAn=%8.4lf KZCGm=%8.4lf KZCGn=%8.4lf\n",
	       i+1,myrank,your_rank,myrank,delta,
	       EAAm_Xi,ECGm_Xi,EZm_Xi,EAAn_Xi,ECGn_Xi,EZn_Xi,
	       EAAm_Xj,ECGm_Xj,EZm_Xj,EAAn_Xj,ECGn_Xj,EZn_Xj,
	       ZData.KZAA,KZAAn,ZData.KZCG,KZCGn);
    }

    if (myrank!=0) MPI_Send(&exchange_data, 4, MPI_INT, 0, tag, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    if (myrank==0) {
      for (j=1;j<numRE;++j) {
	MPI_Recv(&exchange_data_receiv, 4, MPI_INT, j, tag, MPI_COMM_WORLD, status);
	exchange_data_sum[j*4+0]=exchange_data_receiv[0];
	exchange_data_sum[j*4+1]=exchange_data_receiv[1];
	exchange_data_sum[j*4+2]=exchange_data_receiv[2];
	exchange_data_sum[j*4+3]=exchange_data_receiv[3];
	m=exchange_data_receiv[2];
	n=exchange_data_receiv[3];
	//	printf("%d:210 j=%d m=%d n=%d\n",myrank,j,m,n);
	if ( m != -1 && n != -1 ) acc_ratio[m][n]+=1.0;
      }
      for (j=1;j<numRE;++j) {
	REMD_exchange_purmutation_funcs(exchange_data_sum[j*4+0],exchange_data_sum[j*4+1],
					exchange_data_sum[j*4+2],exchange_data_sum[j*4+3]);
      }

      for (j=0;j<numRE;++j) exchanged_data[j]=index_parameters[j];
      for (j=0;j<numRE;++j) exchanged_data[j+numRE]=index_replicas[j];

      for (j=1;j<numRE;++j) {
	MPI_Send(exchanged_data, numRE*2, MPI_INT, j, tag, MPI_COMM_WORLD);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (myrank!=0) {
      MPI_Recv(exchanged_data_receiv, numRE*2, MPI_INT, 0, tag, MPI_COMM_WORLD, status);
      for (j=0;j<numRE;++j) {
	index_parameters[j]=exchanged_data_receiv[j];
	index_replicas[j]=exchanged_data_receiv[numRE+j];
      }
    }
    
    /******************************************************************************************************/
    /* if (myrank==1) {											  */
    /*   printf("%d: 244 \n",myrank);									  */
    /*   for (j=0;j<numRE;++j) {									  */
    /* 	printf("i=%d-m=%d i=%d-m=%d ",j,REMD_purmutation_func(j),j,index_parameters[j]);		  */
    /* 	printf("m=%d-i=%d m=%d-i=%d\n",j,REMD_purmutation_inverse_func(j),j,index_replicas[j]);		  */
    /*   }												  */
    /*   printf("%d: 249 \n",myrank);									  */
    /* }												  */
    /******************************************************************************************************/
    MPI_Barrier(MPI_COMM_WORLD);
  }

  for (i=0;i<numRE;++i) for (j=0;j<numRE;++j) acc_ratio[i][j]=acc_ratio[i][j]/numEX*2.0;
  
  return acc_ratio;
}

double **MPI_CGFGTREM_TACCM_ABAbMD_NH_new_Amber_PROTEINS2008(int myrank,int num_procs,
							     int tag, MPI_Status* status,
							     int numRE, int numEX, double *KZAA, double *KZCG,
							     struct AADataforREMD_Amber AAData,
							     struct CGDataforREMD_PROTEINS2008 CGData,
							     struct TACCMDataforREMD_A_P2008 ZData,
							     struct AACGCommonDataforREMD_A_P2008 CData,
							     double T0AA,double T0CG, double T0Z, 
							     int numstep, int interval, 
							     double dt,double tau, double tau2,
							     double UNITT, double k_B, double pi, FILE* logfile ) {
  int i,j,k,c,l=0;
  int m,n;
  int your_rank;
  int count;

  double betaAA,betaCG,betaZ;

  double EAAm_Xi,EAAm_Xj,EAAn_Xi,EAAn_Xj;
  double ECGm_Xi,ECGm_Xj,ECGn_Xi,ECGn_Xj;
  double EZm_Xi,EZm_Xj,EZn_Xi,EZn_Xj;

  double KZAAn,KZCGn;

  double *crdAA_receiv,*crdCG_receiv;
  double *Z_receiv;
  double *crd_temp,*vel_temp;
  double *Z_temp,*vel_Z_temp;

  double *SData,*RData,*SKZ,*RKZ;
  double delta=0.0;

  int exchange_data[4],exchange_data_receiv[4],*exchange_data_sum,*exchanged_data,*exchanged_data_receiv;
  double **acc_ratio;

  unsigned long init[4]={0x123, 0x234, 0x345, 0x456},length=4;
  init_by_array(init,length);

  betaAA=1.0/(k_B*T0AA);
  betaCG=1.0/(k_B*T0CG);
  betaZ=1.0/(k_B*T0Z);

  SData=(double *)gcemalloc(sizeof(double)*6/*8*/);
  RData=(double *)gcemalloc(sizeof(double)*6/*8*/);
  SKZ=(double *)gcemalloc(sizeof(double)*2);
  RKZ=(double *)gcemalloc(sizeof(double)*2);

  exchange_data_sum=(int *)gcemalloc(sizeof(int)*numRE*4);
  exchanged_data=(int *)gcemalloc(sizeof(int)*numRE*2);
  exchanged_data_receiv=(int *)gcemalloc(sizeof(int)*numRE*2);

  acc_ratio=(double **)gcemalloc(sizeof(double *)*numRE);
  for (i=0;i<numRE;++i) acc_ratio[i]=(double *)gcemalloc(sizeof(double)*numRE);
  for (i=0;i<numRE;++i) for (j=0;j<numRE;++j) acc_ratio[i][j]=0.0;

  //  printf("myrank=%d:340\n",myrank);

  for (i=0;i<numEX;++i) {
    runTACCM_CGFG_ABAbMD_NH_new_Amber_PROTEINS2008(// AA /////////////////////////////////////////////////////////
						   AAData.crd,AAData.q,AAData.qvel,
						   AAData.predict,AAData.correct,
						   AAData.s,AAData.s_vel,AAData.predict_s,AAData.correct_s,
						   AAData.gzi,AAData.gzi_vel,AAData.predict_gzi,AAData.correct_gzi,
						   CData.numclutparent,CData.terminal,CData.origin,
						   AAData.vel_Term,AAData.predict_Term,AAData.predict_Term2,
						   AAData.correct_Term,AAData.correct_Term2,
						   AAData.clt,AAData.Q,
						   AAData.e,AAData.f,AAData.T,T0AA,AAData.KEo,
						   AAData.avePE,AAData.aveKE,AAData.aveT,
						   AAData.varPE,AAData.varKE,AAData.varT,
						   AAData.nc_id_MCD,AAData.outputfile,
						   // CG /////////////////////////////////////////////////////////
						   CGData.crd,CGData.q,CGData.qvel,
						   CGData.predict,CGData.correct,
						   CGData.s,CGData.s_vel,CGData.predict_s,CGData.correct_s,
						   CGData.gzi,CGData.gzi_vel,CGData.predict_gzi,CGData.correct_gzi,
						   CData.numclutparent,CData.terminal,CData.origin,
						   CGData.vel_Term,CGData.predict_Term,CGData.predict_Term2,
						   CGData.correct_Term,CGData.correct_Term2,
						   CGData.clt,CGData.Q,
						   CGData.e,CGData.T,T0CG,CGData.KEo,
						   CGData.avePE,CGData.aveKE,CGData.aveT,
						   CGData.varPE,CGData.varKE,CGData.varT,
						   CGData.nc_id_MCD,CGData.outputfile,
						   // Z  /////////////////////////////////////////////////////////
						   ZData.Z,ZData.velZ,ZData.massZ,
						   ZData.predict,ZData.correct,
						   ZData.sZ,ZData.s_velZ,ZData.gziZ,ZData.gzi_velZ,
						   ZData.predict_gziZ,ZData.correct_gziZ,
						   ZData.predict_sZ,ZData.correct_sZ,
						   ZData.QZ,ZData.T,T0Z,ZData.KEo,
						   ZData.numZ,ZData.KZAA,ZData.KZCG,ZData.pairs,
						   ZData.avePEZ,ZData.aveKEZ,ZData.aveTZ,
						   ZData.varPEZ,ZData.varKEZ,ZData.varTZ, 
						   ZData.trjfileZ,ZData.trjfilThetaAA,ZData.trjfilThetaCG,
						   // CM  /////////////////////////////////////////////////////////
						   CData.mass,(CData.numatom),(CData.numclut),(CData.DOF),
						   numstep,interval,&l,dt,tau,tau2,UNITT,k_B,pi,
						   &EAAm_Xi,&ECGm_Xi,&EZm_Xi);

    m=REMD_purmutation_func(myrank);
    fprintf(logfile,"%d-th: m= %d KZAA= %8.4lf KZCG= %8.4lf\n",i+1,m,ZData.KZAA,ZData.KZCG);

    //    printf("myrank=%d:388\n",myrank);

    //    AAData.zeta=0.0; AAData.V_zeta=0.0;
    //    CGData.zeta=0.0; CGData.V_zeta=0.0;
    //    ZData.zetaZ=0.0; ZData.V_zetaZ=0.0;

    if(m%2==0 ) {
      if ( i%2 == 0 ) n=m+1;
      else {
	if ( m==0 )  n=numRE-1;
	else  n=m-1;
      }
    }      
    else {
      if ( i%2 == 0 ) n=m-1;
      else {
	if ( m==numRE-1 )  n=0;
	else  n=m+1;
      }      
    }

    your_rank=REMD_purmutation_inverse_func(n);

    SKZ[0]=ZData.KZAA;
    SKZ[1]=ZData.KZCG;
    //    printf("107\n");
    MPI_Send(SKZ, 2, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD);
    MPI_Recv(RKZ, 2, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD, status);
    //    printf("110\n");
    KZAAn=RKZ[0];
    KZCGn=RKZ[1];

    MPI_Barrier(MPI_COMM_WORLD);

    CE_TACCM_CGAA(AAData.crd,CGData.crd,ZData.Z,
		  (CData.numatom),ZData.numZ,
		  KZAAn,KZCGn,ZData.pairs,pi,
		  &EAAn_Xi,&ECGn_Xi,&EZn_Xi);

    SData[0]=EAAm_Xi;
    SData[1]=ECGm_Xi;
    SData[2]=EZm_Xi;

    SData[3]=EAAn_Xi;
    SData[4]=ECGn_Xi;
    SData[5]=EZn_Xi;

    //    SData[6]=ZData.KZAA;
    //    SData[7]=ZData.KZCG;
    //    printf("148\n");
    MPI_Send(SData, /*8*/6, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD);
    MPI_Recv(RData, /*8*/6, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD, status);
    //    printf("151\n");
    EAAn_Xj=RData[0];
    ECGn_Xj=RData[1];
    EZn_Xj=RData[2];

    EAAm_Xj=RData[3];
    ECGm_Xj=RData[4];
    EZm_Xj=RData[5];

    //    KZAAn=RData[6];
    //    KZCGn=RData[7];

    MPI_Barrier(MPI_COMM_WORLD);

    if (n>m) {
      //      delta=betaAA*((EAAm_Xj-EAAm_Xi)-(EAAn_Xi-EAAn_Xj))
      //      	+betaCG*((ECGm_Xj-ECGm_Xi)-(ECGn_Xi-ECGn_Xj))
      //      	+betaZ*((EZm_Xj-EZm_Xi)-(EZn_Xi-EZn_Xj));

      delta=betaAA*((EAAm_Xj-EAAm_Xi)-(EAAn_Xj-EAAn_Xi))
      	+betaCG*((ECGm_Xj-ECGm_Xi)-(ECGn_Xj-ECGn_Xi))
      	+betaZ*((EZm_Xj-EZm_Xi)-(EZn_Xj-EZn_Xi));

      /**************************************************************/
      /* delta=-betaAA*((EAAm_Xj-EAAm_Xi)-(EAAn_Xi-EAAn_Xj))	    */
      /* 	-betaCG*((ECGm_Xj-ECGm_Xi)-(ECGn_Xi-ECGn_Xj))	    */
      /* 	-betaZ*((EZm_Xj-EZm_Xi)-(EZn_Xi-EZn_Xj));	    */
      /**************************************************************/

      MPI_Send(&delta, 1, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD);
    }
    else MPI_Recv(&delta, 1, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD, status);

    MPI_Barrier(MPI_COMM_WORLD);

    if((c=Metropolis(delta))==1) {
      if (myrank<your_rank)
	printf("%d-th Acc %d-%d : /_\\ = %8.4lf \n EAAmi=%8.4lf ECGmi=%8.4lf EZmi=%8.4lf EAAni=%8.4lf ECGni=%8.4lf EZni=%8.4lf\n EAAmj=%8.4lf ECGmj=%8.4lf EZmj=%8.4lf EAAnj=%8.4lf ECGnj=%8.4lf EZnj=%8.4lf\n KZAAm=%8.4lf KZAAn=%8.4lf KZCGm=%8.4lf KZCGn=%8.4lf\n",
	       i+1,myrank,your_rank,myrank,delta,
	       EAAm_Xi,ECGm_Xi,EZm_Xi,EAAn_Xi,ECGn_Xi,EZn_Xi,
	       EAAm_Xj,ECGm_Xj,EZm_Xj,EAAn_Xj,ECGn_Xj,EZn_Xj,
	       ZData.KZAA,KZAAn,ZData.KZCG,KZCGn);
      ZData.KZAA=KZAAn;
      ZData.KZCG=KZCGn;

      if (myrank!=0) {
	exchange_data[0]=myrank;
	exchange_data[1]=your_rank;
	exchange_data[2]=m;
	exchange_data[3]=n;
      }
      else {
	REMD_exchange_purmutation_funcs(myrank,your_rank,m,n);
	//	printf("176 m=%d n=%d\n",m,n);
	acc_ratio[m][n]+=1.0;
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
	printf("%d-th Rej %d-%d : /_\\ = %8.4lf \n EAAmi=%8.4lf ECGmi=%8.4lf EZmi=%8.4lf EAAni=%8.4lf ECGni=%8.4lf EZni=%8.4lf\n EAAmj=%8.4lf ECGmj=%8.4lf EZmj=%8.4lf EAAnj=%8.4lf ECGnj=%8.4lf EZnj=%8.4lf\n KZAAm=%8.4lf KZAAn=%8.4lf KZCGm=%8.4lf KZCGn=%8.4lf\n",
	       i+1,myrank,your_rank,myrank,delta,
	       EAAm_Xi,ECGm_Xi,EZm_Xi,EAAn_Xi,ECGn_Xi,EZn_Xi,
	       EAAm_Xj,ECGm_Xj,EZm_Xj,EAAn_Xj,ECGn_Xj,EZn_Xj,
	       ZData.KZAA,KZAAn,ZData.KZCG,KZCGn);
    }

    if (myrank!=0) MPI_Send(&exchange_data, 4, MPI_INT, 0, tag, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    if (myrank==0) {
      for (j=1;j<numRE;++j) {
	MPI_Recv(&exchange_data_receiv, 4, MPI_INT, j, tag, MPI_COMM_WORLD, status);
	exchange_data_sum[j*4+0]=exchange_data_receiv[0];
	exchange_data_sum[j*4+1]=exchange_data_receiv[1];
	exchange_data_sum[j*4+2]=exchange_data_receiv[2];
	exchange_data_sum[j*4+3]=exchange_data_receiv[3];
	m=exchange_data_receiv[2];
	n=exchange_data_receiv[3];
	//	printf("%d:210 j=%d m=%d n=%d\n",myrank,j,m,n);
	if ( m != -1 && n != -1 ) acc_ratio[m][n]+=1.0;
      }
      for (j=1;j<numRE;++j) {
	REMD_exchange_purmutation_funcs(exchange_data_sum[j*4+0],exchange_data_sum[j*4+1],
					exchange_data_sum[j*4+2],exchange_data_sum[j*4+3]);
      }

      for (j=0;j<numRE;++j) exchanged_data[j]=index_parameters[j];
      for (j=0;j<numRE;++j) exchanged_data[j+numRE]=index_replicas[j];

      for (j=1;j<numRE;++j) {
	MPI_Send(exchanged_data, numRE*2, MPI_INT, j, tag, MPI_COMM_WORLD);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (myrank!=0) {
      MPI_Recv(exchanged_data_receiv, numRE*2, MPI_INT, 0, tag, MPI_COMM_WORLD, status);
      for (j=0;j<numRE;++j) {
	index_parameters[j]=exchanged_data_receiv[j];
	index_replicas[j]=exchanged_data_receiv[numRE+j];
      }
    }
    
    /******************************************************************************************************/
    /* if (myrank==1) {											  */
    /*   printf("%d: 244 \n",myrank);									  */
    /*   for (j=0;j<numRE;++j) {									  */
    /* 	printf("i=%d-m=%d i=%d-m=%d ",j,REMD_purmutation_func(j),j,index_parameters[j]);		  */
    /* 	printf("m=%d-i=%d m=%d-i=%d\n",j,REMD_purmutation_inverse_func(j),j,index_replicas[j]);		  */
    /*   }												  */
    /*   printf("%d: 249 \n",myrank);									  */
    /* }												  */
    /******************************************************************************************************/
    MPI_Barrier(MPI_COMM_WORLD);
  }

  for (i=0;i<numRE;++i) for (j=0;j<numRE;++j) acc_ratio[i][j]=acc_ratio[i][j]/numEX*2.0;
  
  return acc_ratio;
}

double **MPI_CGAATREM_TACCM_MD_pep_NHC_MP1998_2PROTEINS2008(int myrank,int num_procs,
							    int tag, MPI_Status* status,
							    int numRE, int numEX, double *KZCG1, double *KZCG2,
							    struct CGDataforREMD_PROTEINS2008 CG1Data,
							    struct CGDataforREMD_PROTEINS2008 CG2Data,
							    struct TACCMDataforREMD_A_P2008 ZData,
							    struct AACGCommonDataforREMD_A_P2008 CData,
							    double T0CG1,double T0CG2, double T0Z, 
							    int numstep, int interval, 
							    double dt,double dt2,
							    double wdt2[3],double wdt4[3], int nc,
							    double UNITT, double k_B, double tau, double pi,
							    FILE* logfile ) {
  int i,j,k,c,l=0;
  int m,n;
  int your_rank;
  int count;

  double betaCG1,betaCG2,betaZ;

  double ECG1m_Xi,ECG1m_Xj,ECG1n_Xi,ECG1n_Xj;
  double ECG2m_Xi,ECG2m_Xj,ECG2n_Xi,ECG2n_Xj;
  double EZm_Xi,EZm_Xj,EZn_Xi,EZn_Xj;

  double KZCG1n,KZCG2n;

  double *crdCG1_receiv,*crdCG2_receiv;
  double *Z_receiv;
  double *crd_temp,*vel_temp;
  double *Z_temp,*vel_Z_temp;

  double *SData,*RData,*SKZ,*RKZ;
  double delta=0.0;

  int exchange_data[4],exchange_data_receiv[4],*exchange_data_sum,*exchanged_data,*exchanged_data_receiv;
  double **acc_ratio;

  unsigned long init[4]={0x123, 0x234, 0x345, 0x456},length=4;
  init_by_array(init,length);

  betaCG1=1.0/(k_B*T0CG1);
  betaCG2=1.0/(k_B*T0CG2);
  betaZ=1.0/(k_B*T0Z);

  SData=(double *)gcemalloc(sizeof(double)*6/*8*/);
  RData=(double *)gcemalloc(sizeof(double)*6/*8*/);
  SKZ=(double *)gcemalloc(sizeof(double)*2);
  RKZ=(double *)gcemalloc(sizeof(double)*2);

  exchange_data_sum=(int *)gcemalloc(sizeof(int)*numRE*4);
  exchanged_data=(int *)gcemalloc(sizeof(int)*numRE*2);
  exchanged_data_receiv=(int *)gcemalloc(sizeof(int)*numRE*2);

  acc_ratio=(double **)gcemalloc(sizeof(double *)*numRE);
  for (i=0;i<numRE;++i) acc_ratio[i]=(double *)gcemalloc(sizeof(double)*numRE);
  for (i=0;i<numRE;++i) for (j=0;j<numRE;++j) acc_ratio[i][j]=0.0;

  for (i=0;i<numEX;++i) {
    runTACCM_CGAA_MD_NHC_MP1998_2PROTEINS2008(// CG1 /////////////////////////////////////////////////////////
					      CG1Data.crd,CG1Data.vel,
					      &(CG1Data.zeta),&(CG1Data.V_zeta),CG1Data.Q,
					      CG1Data.e,
					      CG1Data.T,CG1Data.NfKT,
					      CG1Data.avePE,CG1Data.aveKE,CG1Data.aveT,
					      CG1Data.varPE,CG1Data.varKE,CG1Data.varT,
					      CG1Data.nc_id_MCD,CG1Data.outputfile,
					      // CG2 /////////////////////////////////////////////////////////
					      CG2Data.crd,CG2Data.vel,
					      &(CG2Data.zeta),&(CG2Data.V_zeta),CG2Data.Q,
					      CG2Data.e,
					      CG2Data.T,CG2Data.NfKT,
					      CG2Data.avePE,CG2Data.aveKE,CG2Data.aveT,
					      CG2Data.varPE,CG2Data.varKE,CG2Data.varT,
					      CG2Data.nc_id_MCD,CG2Data.outputfile,
					      // Z  /////////////////////////////////////////////////////////
					      ZData.Z,ZData.velZ,ZData.massZ,
					      &(ZData.zetaZ),&(ZData.V_zetaZ),
					      ZData.QZ,ZData.NfKTZ,ZData.T,
					      ZData.numZ,ZData.KZAA,ZData.KZCG,ZData.pairs,
					      ZData.avePEZ,ZData.aveKEZ,ZData.aveTZ,
					      ZData.varPEZ,ZData.varKEZ,ZData.varTZ, 
					      ZData.trjfileZ,ZData.trjfilThetaAA,ZData.trjfilThetaCG,
					      // CM  /////////////////////////////////////////////////////////
					      CData.mass,(CData.numatom),(CData.numheavyatom),
					      numstep,interval,&l,
					      dt,dt2,wdt2,wdt4,nc,UNITT,k_B,pi,
					      &ECG1m_Xi,&ECG2m_Xi,&EZm_Xi);

    m=REMD_purmutation_func(myrank);
    fprintf(logfile,"%d-th: m= %d KZCG1= %8.4lf KZCG2= %8.4lf\n",i+1,m,ZData.KZAA,ZData.KZCG);

    CG1Data.zeta=0.0; CG1Data.V_zeta=0.0;
    CG2Data.zeta=0.0; CG2Data.V_zeta=0.0;
    ZData.zetaZ=0.0; ZData.V_zetaZ=0.0;

    if(m%2==0 ) {
      if ( i%2 == 0 ) n=m+1;
      else {
	if ( m==0 )  n=numRE-1;
	else  n=m-1;
      }
    }      
    else {
      if ( i%2 == 0 ) n=m-1;
      else {
	if ( m==numRE-1 )  n=0;
	else  n=m+1;
      }      
    }

    your_rank=REMD_purmutation_inverse_func(n);

    SKZ[0]=ZData.KZAA;
    SKZ[1]=ZData.KZCG;
    //    printf("107\n");
    MPI_Send(SKZ, 2, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD);
    MPI_Recv(RKZ, 2, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD, status);
    //    printf("110\n");
    KZCG1n=RKZ[0];
    KZCG2n=RKZ[1];

    MPI_Barrier(MPI_COMM_WORLD);

    CE_TACCM_CGAA(CG1Data.crd,CG2Data.crd,ZData.Z,
		  (CData.numatom),ZData.numZ,
		  KZCG1n,KZCG2n,ZData.pairs,pi,
		  &ECG1n_Xi,&ECG2n_Xi,&EZn_Xi);

    SData[0]=ECG1m_Xi;
    SData[1]=ECG2m_Xi;
    SData[2]=EZm_Xi;

    SData[3]=ECG1n_Xi;
    SData[4]=ECG2n_Xi;
    SData[5]=EZn_Xi;

    //    SData[6]=ZData.KZCG1;
    //    SData[7]=ZData.KZCG2;
    //    printf("148\n");
    MPI_Send(SData, /*8*/6, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD);
    MPI_Recv(RData, /*8*/6, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD, status);
    //    printf("151\n");
    ECG1n_Xj=RData[0];
    ECG2n_Xj=RData[1];
    EZn_Xj=RData[2];

    ECG1m_Xj=RData[3];
    ECG2m_Xj=RData[4];
    EZm_Xj=RData[5];

    //    KZCG1n=RData[6];
    //    KZCG2n=RData[7];

    MPI_Barrier(MPI_COMM_WORLD);

    if (n>m) {
      //      delta=betaCG1*((ECG1m_Xj-ECG1m_Xi)-(ECG1n_Xi-ECG1n_Xj))
      //      	+betaCG2*((ECG2m_Xj-ECG2m_Xi)-(ECG2n_Xi-ECG2n_Xj))
      //      	+betaZ*((EZm_Xj-EZm_Xi)-(EZn_Xi-EZn_Xj));

      delta=betaCG1*((ECG1m_Xj-ECG1m_Xi)-(ECG1n_Xj-ECG1n_Xi))
      	+betaCG2*((ECG2m_Xj-ECG2m_Xi)-(ECG2n_Xj-ECG2n_Xi))
      	+betaZ*((EZm_Xj-EZm_Xi)-(EZn_Xj-EZn_Xi));

      /**************************************************************/
      /* delta=-betaCG1*((ECG1m_Xj-ECG1m_Xi)-(ECG1n_Xi-ECG1n_Xj))	    */
      /* 	-betaCG2*((ECG2m_Xj-ECG2m_Xi)-(ECG2n_Xi-ECG2n_Xj))	    */
      /* 	-betaZ*((EZm_Xj-EZm_Xi)-(EZn_Xi-EZn_Xj));	    */
      /**************************************************************/

      MPI_Send(&delta, 1, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD);
    }
    else MPI_Recv(&delta, 1, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD, status);

    MPI_Barrier(MPI_COMM_WORLD);

    if((c=Metropolis(delta))==1) {
      if (myrank<your_rank)
	printf("%d-th Acc %d-%d : /_\\ = %8.4lf \n ECG1mi=%8.4lf ECG2mi=%8.4lf EZmi=%8.4lf ECG1ni=%8.4lf ECG2ni=%8.4lf EZni=%8.4lf\n ECG1mj=%8.4lf ECG2mj=%8.4lf EZmj=%8.4lf ECG1nj=%8.4lf ECG2nj=%8.4lf EZnj=%8.4lf\n KZCG1m=%8.4lf KZCG1n=%8.4lf KZCG2m=%8.4lf KZCG2n=%8.4lf\n",
	       i+1,myrank,your_rank,myrank,delta,
	       ECG1m_Xi,ECG2m_Xi,EZm_Xi,ECG1n_Xi,ECG2n_Xi,EZn_Xi,
	       ECG1m_Xj,ECG2m_Xj,EZm_Xj,ECG1n_Xj,ECG2n_Xj,EZn_Xj,
	       ZData.KZAA,KZCG1n,ZData.KZCG,KZCG2n);
      ZData.KZAA=KZCG1n;
      ZData.KZCG=KZCG2n;

      if (myrank!=0) {
	exchange_data[0]=myrank;
	exchange_data[1]=your_rank;
	exchange_data[2]=m;
	exchange_data[3]=n;
      }
      else {
	REMD_exchange_purmutation_funcs(myrank,your_rank,m,n);
	//	printf("176 m=%d n=%d\n",m,n);
	acc_ratio[m][n]+=1.0;
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
	printf("%d-th Rej %d-%d : /_\\ = %8.4lf \n ECG1mi=%8.4lf ECG2mi=%8.4lf EZmi=%8.4lf ECG1ni=%8.4lf ECG2ni=%8.4lf EZni=%8.4lf\n ECG1mj=%8.4lf ECG2mj=%8.4lf EZmj=%8.4lf ECG1nj=%8.4lf ECG2nj=%8.4lf EZnj=%8.4lf\n KZCG1m=%8.4lf KZCG1n=%8.4lf KZCG2m=%8.4lf KZCG2n=%8.4lf\n",
	       i+1,myrank,your_rank,myrank,delta,
	       ECG1m_Xi,ECG2m_Xi,EZm_Xi,ECG1n_Xi,ECG2n_Xi,EZn_Xi,
	       ECG1m_Xj,ECG2m_Xj,EZm_Xj,ECG1n_Xj,ECG2n_Xj,EZn_Xj,
	       ZData.KZAA,KZCG1n,ZData.KZCG,KZCG2n);
    }

    if (myrank!=0) MPI_Send(&exchange_data, 4, MPI_INT, 0, tag, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    if (myrank==0) {
      for (j=1;j<numRE;++j) {
	MPI_Recv(&exchange_data_receiv, 4, MPI_INT, j, tag, MPI_COMM_WORLD, status);
	exchange_data_sum[j*4+0]=exchange_data_receiv[0];
	exchange_data_sum[j*4+1]=exchange_data_receiv[1];
	exchange_data_sum[j*4+2]=exchange_data_receiv[2];
	exchange_data_sum[j*4+3]=exchange_data_receiv[3];
	m=exchange_data_receiv[2];
	n=exchange_data_receiv[3];
	//	printf("%d:210 j=%d m=%d n=%d\n",myrank,j,m,n);
	if ( m != -1 && n != -1 ) acc_ratio[m][n]+=1.0;
      }
      for (j=1;j<numRE;++j) {
	REMD_exchange_purmutation_funcs(exchange_data_sum[j*4+0],exchange_data_sum[j*4+1],
					exchange_data_sum[j*4+2],exchange_data_sum[j*4+3]);
      }

      for (j=0;j<numRE;++j) exchanged_data[j]=index_parameters[j];
      for (j=0;j<numRE;++j) exchanged_data[j+numRE]=index_replicas[j];

      for (j=1;j<numRE;++j) {
	MPI_Send(exchanged_data, numRE*2, MPI_INT, j, tag, MPI_COMM_WORLD);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (myrank!=0) {
      MPI_Recv(exchanged_data_receiv, numRE*2, MPI_INT, 0, tag, MPI_COMM_WORLD, status);
      for (j=0;j<numRE;++j) {
	index_parameters[j]=exchanged_data_receiv[j];
	index_replicas[j]=exchanged_data_receiv[numRE+j];
      }
    }
    
    /******************************************************************************************************/
    /* if (myrank==1) {											  */
    /*   printf("%d: 244 \n",myrank);									  */
    /*   for (j=0;j<numRE;++j) {									  */
    /* 	printf("i=%d-m=%d i=%d-m=%d ",j,REMD_purmutation_func(j),j,index_parameters[j]);		  */
    /* 	printf("m=%d-i=%d m=%d-i=%d\n",j,REMD_purmutation_inverse_func(j),j,index_replicas[j]);		  */
    /*   }												  */
    /*   printf("%d: 249 \n",myrank);									  */
    /* }												  */
    /******************************************************************************************************/
    MPI_Barrier(MPI_COMM_WORLD);
  }

  for (i=0;i<numRE;++i) for (j=0;j<numRE;++j) acc_ratio[i][j]=acc_ratio[i][j]/numEX*2.0;
  
  return acc_ratio;
}

void  CGAAREMDreadInputs_test(FILE *inputfile,int numatom,int numRE,int myrank,
			      double *crdAA,double *velAA, double *crdCG,double *velCG,
			      double *KZAA, double *KZCG){
  int i,j,k;
  int c;
  int  d;
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
    if (State==AAINPF) {
      if (c==' ' || c=='\n') {
	if (j>0) {
	  crdfilename[j]='\0';
	  crdfile=efopen(crdfilename,"r");
	  getline(&line,&len,crdfile);
	  fscanf(crdfile,"%d",&d);
	  if (i==myrank)  for (k=0;k<numatom*3;++k) fscanf(crdfile,"%lf",&crdAA[k]);
	  else for (k=0;k<numatom*3;++k) fscanf(crdfile,"%lf",&fdummy);
	  fclose(crdfile);
	  State=CGINPF;
	  j=0;
	}
      }
      else {
	crdfilename[j]=c;
	++j;
      }
    }
    else if (State==CGINPF) {
      if (c==' ' || c=='\n') {
	if (j>0) {
	  crdfilename[j]='\0';
	  crdfile=efopen(crdfilename,"r");
	  getline(&line,&len,crdfile);
	  fscanf(crdfile,"%d",&d);
	  if (i==myrank) for (k=0;k<numatom*3;++k) fscanf(crdfile,"%lf",&crdCG[k]);
	  else for (k=0;k<numatom*3;++k) fscanf(crdfile,"%lf",&fdummy);
	  fclose(crdfile);
	  State=AAKZ;
	  j=0;
	}
      }
      else {
	crdfilename[j]=c;
	++j;
      }
    }
    else if (State==AAKZ) {
      if (isdigit(c)) {
	d=(c-'0');
	f1=f1*10+(double)d;
	Tflag=ON;
      }
      else if (c==' ' || c=='\n') {
	if (Tflag==ON) {
	  KZAA[i]=(double)f1;
	  f1=0.0;
	  State=CGKZ;
	  Tflag=OFF;
	}
      }
      else {
	printf("error\n");
	exit(1);
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
	  KZCG[i]=(double)f1;
	  ++i;
	  if (i==numRE) break;
	  f1=0.0;
	  State=AAINPF;
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
