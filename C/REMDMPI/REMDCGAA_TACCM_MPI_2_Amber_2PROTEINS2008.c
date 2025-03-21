
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "REMDCGAA_TACCM_MPI_2_Amber_PROTEINS2008.h"
#include "REMDCGAA_TACCM_MPI_2_Amber_2PROTEINS2008.h"
#include "REMDCGAA_TACCM_MPI_2_Amber_2PROTEINS2008.h"
#include "REMD_functions.h"

#define ON  0
#define OFF 1

double **MPI_CGFGTREM_TACCM_ABAbMD_NH_new_Amber_2PROTEINS2008(int myrank,int num_procs,
							     int tag, MPI_Status* status,
							     int numRE, int numEX, double *KZAA, double **KZCG,
							     struct AADataforREMD_Amber AAData,
							     struct CGDataforREMD_PROTEINS2008 CGData[2],
							     struct TACCMDataforREMD_A_2P2008 ZData,
							     struct AACGCommonDataforREMD_A_P2008 CData,
							     double T0AA,double T0CG[2], double T0Z, 
							     int numstep, int interval, 
							     double dt,double tau, double tau2,
							     double UNITT, double k_B, double pi, FILE* logfile ) {
  int i,j,k,c,l=0;
  int m,n;
  int your_rank;
  int count;

  double betaAA,betaCG[2],betaZ;

  double EAAm_Xi,EAAm_Xj,EAAn_Xi,EAAn_Xj;
  double ECGm_Xi[2],ECGm_Xj[2],ECGn_Xi[2],ECGn_Xj[2];
  double EZm_Xi,EZm_Xj,EZn_Xi,EZn_Xj;

  double KZAAn,KZCGn[2];

  double *crdAA_receiv,**crdCG_receiv;
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
  for (i=0;i<2;++i)  betaCG[i]=1.0/(k_B*T0CG[i]);
  betaZ=1.0/(k_B*T0Z);

  SData=(double *)gcemalloc(sizeof(double)/**6*/*8);
  RData=(double *)gcemalloc(sizeof(double)/**6*/*8);
  SKZ=(double *)gcemalloc(sizeof(double)*/*2*/3);
  RKZ=(double *)gcemalloc(sizeof(double)*/*2*/3);

  exchange_data_sum=(int *)gcemalloc(sizeof(int)*numRE*4);
  exchanged_data=(int *)gcemalloc(sizeof(int)*numRE*2);
  exchanged_data_receiv=(int *)gcemalloc(sizeof(int)*numRE*2);

  acc_ratio=(double **)gcemalloc(sizeof(double *)*numRE);
  for (i=0;i<numRE;++i) acc_ratio[i]=(double *)gcemalloc(sizeof(double)*numRE);
  for (i=0;i<numRE;++i) for (j=0;j<numRE;++j) acc_ratio[i][j]=0.0;

  //  printf("myrank=%d:340\n",myrank);

  for (i=0;i<numEX;++i) {
    runTACCM_2CG1FG_ABAbMD_NH_new_Amber_PROTEINS2008(// AA /////////////////////////////////////////////////////////
						     AAData.crd,AAData.q,AAData.qvel,
						     AAData.predict,AAData.correct,
						     AAData.s,AAData.s_vel,AAData.predict_s,AAData.correct_s,
						     AAData.gzi,AAData.gzi_vel,
						     AAData.predict_gzi,AAData.correct_gzi,
						     CData.numclutparent,CData.terminal,CData.origin,
						     AAData.vel_Term,AAData.predict_Term,AAData.predict_Term2,
						     AAData.correct_Term,AAData.correct_Term2,
						     AAData.clt,AAData.Q,
						     AAData.e,AAData.f,AAData.T,T0AA,AAData.KEo,
						     AAData.avePE,AAData.aveKE,AAData.aveT,
						     AAData.varPE,AAData.varKE,AAData.varT,
						     AAData.nc_id_MCD,AAData.outputfile,
						     // CG1 ////////////////////////////////////////////////////////
						     CGData[0].crd,CGData[0].q,CGData[0].qvel,
						     CGData[0].predict,CGData[0].correct,
						     CGData[0].s,CGData[0].s_vel,
						     CGData[0].predict_s,CGData[0].correct_s,
						     CGData[0].gzi,CGData[0].gzi_vel,
						     CGData[0].predict_gzi,CGData[0].correct_gzi,
						     CData.numclutparent,CData.terminal,CData.origin,
						     CGData[0].vel_Term,
						     CGData[0].predict_Term,CGData[0].predict_Term2,
						     CGData[0].correct_Term,CGData[0].correct_Term2,
						     CGData[0].clt,CGData[0].Q,
						     CGData[0].e,CGData[0].T,T0CG[0],CGData[0].KEo,
						     CGData[0].avePE,CGData[0].aveKE,CGData[0].aveT,
						     CGData[0].varPE,CGData[0].varKE,CGData[0].varT,
						     CGData[0].nc_id_MCD,CGData[0].outputfile,
						     // CG2 ////////////////////////////////////////////////////////
						     CGData[1].crd,CGData[1].q,CGData[1].qvel,
						     CGData[1].predict,CGData[1].correct,
						     CGData[1].s,CGData[1].s_vel,
						     CGData[1].predict_s,CGData[1].correct_s,
						     CGData[1].gzi,CGData[1].gzi_vel,
						     CGData[1].predict_gzi,CGData[1].correct_gzi,
						     CData.numclutparent,CData.terminal,CData.origin,
						     CGData[1].vel_Term,
						     CGData[1].predict_Term,CGData[1].predict_Term2,
						     CGData[1].correct_Term,CGData[1].correct_Term2,
						     CGData[1].clt,CGData[1].Q,
						     CGData[1].e,CGData[1].T,T0CG[1],CGData[1].KEo,
						     CGData[1].avePE,CGData[1].aveKE,CGData[1].aveT,
						     CGData[1].varPE,CGData[1].varKE,CGData[1].varT,
						     CGData[1].nc_id_MCD,CGData[1].outputfile,
						     // Z  /////////////////////////////////////////////////////////
						     ZData.Z,ZData.velZ,ZData.massZ,
						     ZData.predict,ZData.correct,
						     ZData.sZ,ZData.s_velZ,ZData.gziZ,ZData.gzi_velZ,
						     ZData.predict_gziZ,ZData.correct_gziZ,
						     ZData.predict_sZ,ZData.correct_sZ,
						     ZData.QZ,ZData.T,T0Z,ZData.KEo,
						     ZData.numZ,ZData.KZAA,ZData.KZCG[0],ZData.KZCG[1],ZData.pairs,
						     ZData.avePEZ,ZData.aveKEZ,ZData.aveTZ,
						     ZData.varPEZ,ZData.varKEZ,ZData.varTZ, 
						     ZData.trjfileZ,ZData.trjfilThetaAA,
						     ZData.trjfilThetaCG1,ZData.trjfilThetaCG2,
						     // CM  ///////////////////////////////////////////////////////
						     CData.mass,(CData.numatom),(CData.numclut),(CData.DOF),
						     numstep,interval,&l,dt,tau,tau2,UNITT,k_B,pi,
						     &EAAm_Xi,&ECGm_Xi[0],&ECGm_Xi[1],&EZm_Xi);

    m=REMD_purmutation_func(myrank);
    fprintf(logfile,"%d-th: m= %d KZAA= %8.4lf KZCG1= %8.4lf KZCG2= %8.4lf\n",
	    i+1,m,ZData.KZAA,ZData.KZCG[0],ZData.KZCG[1]);

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
    SKZ[1]=ZData.KZCG[0];
    SKZ[2]=ZData.KZCG[1];
    //    printf("107\n");
    MPI_Send(SKZ, /*2*/3, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD);
    MPI_Recv(RKZ, /*2*/3, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD, status);
    //    printf("110\n");
    KZAAn=RKZ[0];
    KZCGn[0]=RKZ[1];
    KZCGn[1]=RKZ[2];

    MPI_Barrier(MPI_COMM_WORLD);

    CE_TACCM_2CG1AA(AAData.crd,CGData[0].crd,CGData[1].crd,ZData.Z,
		    (CData.numatom),ZData.numZ,
		    KZAAn,KZCGn[0],KZCGn[1],ZData.pairs,pi,
		    &EAAn_Xi,&ECGn_Xi[0],&ECGn_Xi[1],&EZn_Xi);

    SData[0]=EAAm_Xi;
    SData[1]=ECGm_Xi[0];
    SData[2]=ECGm_Xi[1];
    SData[3]=EZm_Xi;

    SData[4]=EAAn_Xi;
    SData[5]=ECGn_Xi[0];
    SData[6]=ECGn_Xi[1];
    SData[7]=EZn_Xi;

    //    SData[6]=ZData.KZAA;
    //    SData[7]=ZData.KZCG;
    //    printf("148\n");
    MPI_Send(SData, 8/*6*/, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD);
    MPI_Recv(RData, 8/*6*/, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD, status);
    //    printf("151\n");
    EAAn_Xj=RData[0];
    ECGn_Xj[0]=RData[1];
    ECGn_Xj[1]=RData[2];
    EZn_Xj=RData[3];

    EAAm_Xj=RData[4];
    ECGm_Xj[0]=RData[5];
    ECGm_Xj[1]=RData[6];
    EZm_Xj=RData[7];

    //    KZAAn=RData[6];
    //    KZCGn=RData[7];

    MPI_Barrier(MPI_COMM_WORLD);

    if (n>m) {
      //      delta=betaAA*((EAAm_Xj-EAAm_Xi)-(EAAn_Xi-EAAn_Xj))
      //      	+betaCG*((ECGm_Xj-ECGm_Xi)-(ECGn_Xi-ECGn_Xj))
      //      	+betaZ*((EZm_Xj-EZm_Xi)-(EZn_Xi-EZn_Xj));

      delta=betaAA*((EAAm_Xj-EAAm_Xi)-(EAAn_Xj-EAAn_Xi))
      	+betaCG[0]*((ECGm_Xj[0]-ECGm_Xi[0])-(ECGn_Xj[0]-ECGn_Xi[0]))
      	+betaCG[1]*((ECGm_Xj[1]-ECGm_Xi[1])-(ECGn_Xj[1]-ECGn_Xi[1]))
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
	printf("%d-th Acc %d-%d : /_\\ = %8.4lf \n EAAmi=%8.4lf ECG1mi=%8.4lf EC2Gmi=%8.4lf EZmi=%8.4lf EAAni=%8.4lf ECG1ni=%8.4lf ECG2ni=%8.4lf EZni=%8.4lf\n EAAmj=%8.4lf ECG1mj=%8.4lf ECG2mj=%8.4lf EZmj=%8.4lf EAAnj=%8.4lf ECG1nj=%8.4lf ECG2nj=%8.4lf EZnj=%8.4lf\n KZAAm=%8.4lf KZAAn=%8.4lf KZCG1m=%8.4lf KZCG2m=%8.4lf KZCG1n=%8.4lf KZCG2n=%8.4lf\n",
	       i+1,myrank,your_rank,myrank,delta,
	       EAAm_Xi,ECGm_Xi[0],ECGm_Xi[1],EZm_Xi,EAAn_Xi,ECGn_Xi[0],ECGn_Xi[1],EZn_Xi,
	       EAAm_Xj,ECGm_Xj[0],ECGm_Xj[1],EZm_Xj,EAAn_Xj,ECGn_Xj[0],ECGn_Xj[1],EZn_Xj,
	       ZData.KZAA,KZAAn,ZData.KZCG[0],ZData.KZCG[1],KZCGn[0],KZCGn[1]);
      ZData.KZAA=KZAAn;
      ZData.KZCG[0]=KZCGn[0];
      ZData.KZCG[1]=KZCGn[1];

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
	printf("%d-th Rej %d-%d : /_\\ = %8.4lf \n EAAmi=%8.4lf ECG1mi=%8.4lf ECG2mi=%8.4lf EZmi=%8.4lf EAAni=%8.4lf ECG1ni=%8.4lf ECG2ni=%8.4lf EZni=%8.4lf\n EAAmj=%8.4lf ECG1mj=%8.4lf ECG2mj=%8.4lf EZmj=%8.4lf EAAnj=%8.4lf ECG1nj=%8.4lf ECG2nj=%8.4lf EZnj=%8.4lf\n KZAAm=%8.4lf KZAAn=%8.4lf KZCG1m=%8.4lf KZCG2m=%8.4lf KZCG1n=%8.4lf KZCG2n=%8.4lf\n",
	       i+1,myrank,your_rank,myrank,delta,
	       EAAm_Xi,ECGm_Xi[0],ECGm_Xi[1],EZm_Xi,EAAn_Xi,ECGn_Xi[0],ECGn_Xi[1],EZn_Xi,
	       EAAm_Xj,ECGm_Xj[0],ECGm_Xj[1],EZm_Xj,EAAn_Xj,ECGn_Xj[0],ECGn_Xj[1],EZn_Xj,
	       ZData.KZAA,KZAAn,ZData.KZCG[0],ZData.KZCG[1],KZCGn[0],KZCGn[1]);
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

void  CGAAREMDreadInputs_Amber_PROTEINS2008_Amber_hybrid_1FG2CG(FILE *inputfile,int numatom,int numRE,int myrank,
								double *crdAA,double *velAA, 
								double *crdCG1,double *velCG1,
								double *crdCG2,double *velCG2,
								double *KZAA, double *KZCG1, double *KZCG2){
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
  State=AA1INPF;
  Tflag=OFF;

  while ((c=getc(inputfile))!=-1){
    if (State==AA1INPF) {
      if (c==' ' || c=='\n') {
	if (j>0) {
	  crdfilename[j]='\0';
	  crdfile=efopen(crdfilename,"r");
	  getline(&line,&len,crdfile);
	  fscanf(crdfile,"%d",&d);
	  if (i==myrank)  for (k=0;k<numatom*3;++k) fscanf(crdfile,"%lf",&crdAA[k]);
	  else for (k=0;k<numatom*3;++k) fscanf(crdfile,"%lf",&fdummy);
	  fclose(crdfile);
	  State=CG1INPF;
	  j=0;
	}
      }
      else {
	crdfilename[j]=c;
	++j;
      }
    }
    else if (State==CG1INPF) {
      if (c==' ' || c=='\n') {
	if (j>0) {
	  crdfilename[j]='\0';
	  crdfile=efopen(crdfilename,"r");
	  getline(&line,&len,crdfile);
	  fscanf(crdfile,"%d",&d);
	  if (i==myrank) for (k=0;k<numatom*3;++k) fscanf(crdfile,"%lf",&crdCG1[k]);
	  else for (k=0;k<numatom*3;++k) fscanf(crdfile,"%lf",&fdummy);
	  fclose(crdfile);
	  //	  State=AAKZ;
	  State=CG2INPF;
	  j=0;
	}
      }
      else {
	crdfilename[j]=c;
	++j;
      }
    }
    else if (State==CG2INPF) {
      if (c==' ' || c=='\n') {
	if (j>0) {
	  crdfilename[j]='\0';
	  crdfile=efopen(crdfilename,"r");
	  getline(&line,&len,crdfile);
	  fscanf(crdfile,"%d",&d);
	  if (i==myrank) for (k=0;k<numatom*3;++k) fscanf(crdfile,"%lf",&crdCG2[k]);
	  else for (k=0;k<numatom*3;++k) fscanf(crdfile,"%lf",&fdummy);
	  fclose(crdfile);
	  State=AA1KZ;
	  j=0;
	}
      }
      else {
	crdfilename[j]=c;
	++j;
      }
    }
    else if (State==AA1KZ) {
      if (isdigit(c)) {
	d=(c-'0');
	f1=f1*10+(double)d;
	Tflag=ON;
      }
      else if (c==' ' || c=='\n') {
	if (Tflag==ON) {
	  KZAA[i]=(double)f1;
	  f1=0.0;
	  State=CG1KZ;
	  Tflag=OFF;
	}
      }
      else {
	printf("error\n");
	exit(1);
      }
    }
    else if (State==CG1KZ) {
      if (isdigit(c)) {
	d=(c-'0');
	f1=f1*10+(double)d;
	Tflag=ON;
      }
      else if (c==' ' || c=='\n') {
	if (Tflag==ON) {
	  KZCG1[i]=(double)f1;
	  //	  ++i;
	  //	  if (i==numRE) break;
	  f1=0.0;
	  State=CG2KZ;
	  Tflag=OFF;
	}
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
	  KZCG2[i]=(double)f1;
	  ++i;
	  if (i==numRE) break;
	  f1=0.0;
	  State=AA1INPF;
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

double CE_TACCM_2CG1AA(double *crdAA,double *crdCG1,double *crdCG2,double *Z, int numatom,int numZ,
		       double KZAA,double KZCG1,double KZCG2,int **pairs,double pi,
		       double *EAA,double *ECG1,double *ECG2,double *EZ){
  int i,j;
  double *thetaAA,*thetaCG1,*thetaCG2;
  double delta;

  thetaAA=(double *)gcemalloc(sizeof(double)*numZ);
  thetaCG1=(double *)gcemalloc(sizeof(double)*numZ);
  thetaCG2=(double *)gcemalloc(sizeof(double)*numZ);

  TACCM_CTheta(crdAA,numatom,thetaAA,numZ,pairs,pi);
  TACCM_CTheta(crdCG1,numatom,thetaCG1,numZ,pairs,pi);
  TACCM_CTheta(crdCG2,numatom,thetaCG2,numZ,pairs,pi);

  *EAA=0.0;
  for (i=0;i<numZ;++i) {
    if ((delta=Z[i]-thetaAA[i])>pi) delta-=2.0*pi;
    else if ((delta=Z[i]-thetaAA[i])<-1.0*pi) delta+=2.0*pi;
    *EAA+=0.5*KZAA*delta*delta;
  }

  *ECG1=0.0;
  for (i=0;i<numZ;++i) {
    if ((delta=Z[i]-thetaCG1[i])>pi) delta-=2.0*pi;
    else if ((delta=Z[i]-thetaCG1[i])<-1.0*pi) delta+=2.0*pi;
    *ECG1+=0.5*KZCG1*delta*delta;
  }

  *ECG2=0.0;
  for (i=0;i<numZ;++i) {
    if ((delta=Z[i]-thetaCG2[i])>pi) delta-=2.0*pi;
    else if ((delta=Z[i]-thetaCG2[i])<-1.0*pi) delta+=2.0*pi;
    *ECG2+=0.5*KZCG2*delta*delta;
  }

  *EZ=(*EAA)+(*ECG1)+(*ECG2);
}

double **MPI_CGAATREM_TACCM_MD_pep_NHC_MP1998_Amber_2PROTEINS2008(int myrank,int num_procs,
								  int tag, MPI_Status* status,
								  int numRE, int numEX, double *KZAA, double **KZCG,
								  struct AADataforREMD_Amber AAData,
								  struct CGDataforREMD_PROTEINS2008 CGData[2],
								  struct TACCMDataforREMD_A_2P2008 ZData,
								  struct AACGCommonDataforREMD_A_P2008 CData,
								  double T0AA,double T0CG[2], double T0Z, 
								  int numstep, int interval, 
								  double dt,double dt2,
								  double wdt2[3],double wdt4[3], int nc,
								  double UNITT, double k_B, double tau, double pi,
								  FILE* logfile ) {
  int i,j,k,c,l=0;
  int m,n;
  int your_rank;
  int count;

  double betaAA,betaCG[2],betaZ;

  double EAAm_Xi,EAAm_Xj,EAAn_Xi,EAAn_Xj;
  double ECGm_Xi[2],ECGm_Xj[2],ECGn_Xi[2],ECGn_Xj[2];
  double EZm_Xi,EZm_Xj,EZn_Xi,EZn_Xj;

  double KZAAn,KZCGn[2];

  double *crdAA_receiv,**crdCG_receiv;
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
  betaCG[0]=1.0/(k_B*T0CG[0]);
  betaCG[1]=1.0/(k_B*T0CG[1]);
  betaZ=1.0/(k_B*T0Z);

  SData=(double *)gcemalloc(sizeof(double)/**6*/*8);
  RData=(double *)gcemalloc(sizeof(double)/**6*/*8);
  SKZ=(double *)gcemalloc(sizeof(double)*/*2*/3);
  RKZ=(double *)gcemalloc(sizeof(double)*/*2*/3);

  exchange_data_sum=(int *)gcemalloc(sizeof(int)*numRE*4);
  exchanged_data=(int *)gcemalloc(sizeof(int)*numRE*2);
  exchanged_data_receiv=(int *)gcemalloc(sizeof(int)*numRE*2);

  acc_ratio=(double **)gcemalloc(sizeof(double *)*numRE);
  for (i=0;i<numRE;++i) acc_ratio[i]=(double *)gcemalloc(sizeof(double)*numRE);
  for (i=0;i<numRE;++i) for (j=0;j<numRE;++j) acc_ratio[i][j]=0.0;

  for (i=0;i<numEX;++i) {
    runTACCM_2CG1FG_MD_NHC_MP1998_Amber_PROTEINS2008(// AA /////////////////////////////////////////////////////////
						     AAData.crd,AAData.vel,
						     &(AAData.zeta),&(AAData.V_zeta),AAData.Q,
						     AAData.e,AAData.f,AAData.T,AAData.NfKT,
						     AAData.avePE,AAData.aveKE,AAData.aveT,
						     AAData.varPE,AAData.varKE,AAData.varT,
						     AAData.nc_id_MCD,AAData.outputfile,
						     // CG1 ////////////////////////////////////////////////////////
						     CGData[0].crd,CGData[0].vel,
						     &(CGData[0].zeta),&(CGData[0].V_zeta),CGData[0].Q,
						     CGData[0].e,
						     CGData[0].T,CGData[0].NfKT,
						     CGData[0].avePE,CGData[0].aveKE,CGData[0].aveT,
						     CGData[0].varPE,CGData[0].varKE,CGData[0].varT,
						     CGData[0].nc_id_MCD,CGData[0].outputfile,
						     // CG2 ////////////////////////////////////////////////////////
						     CGData[1].crd,CGData[1].vel,
						     &(CGData[1].zeta),&(CGData[1].V_zeta),CGData[1].Q,
						     CGData[1].e,
						     CGData[1].T,CGData[1].NfKT,
						     CGData[1].avePE,CGData[1].aveKE,CGData[1].aveT,
						     CGData[1].varPE,CGData[1].varKE,CGData[1].varT,
						     CGData[1].nc_id_MCD,CGData[1].outputfile,
						     // Z  /////////////////////////////////////////////////////////
						     ZData.Z,ZData.velZ,ZData.massZ,
						     &(ZData.zetaZ),&(ZData.V_zetaZ),
						     ZData.QZ,ZData.NfKTZ,ZData.T,
						     ZData.numZ,ZData.KZAA,ZData.KZCG[0],ZData.KZCG[1],ZData.pairs,
						     ZData.avePEZ,ZData.aveKEZ,ZData.aveTZ,
						     ZData.varPEZ,ZData.varKEZ,ZData.varTZ, 
						     ZData.trjfileZ,ZData.trjfilThetaAA,
						     ZData.trjfilThetaCG1,ZData.trjfilThetaCG2,
						     // CM  ////////////////////////////////////////////////////////
						     CData.mass,(CData.numatom),(CData.numheavyatom),
						     numstep,interval,&l,
						     dt,dt2,wdt2,wdt4,nc,UNITT,k_B,pi,
						     &EAAm_Xi,&ECGm_Xi[0],&ECGm_Xi[1],&EZm_Xi);

    m=REMD_purmutation_func(myrank);
    fprintf(logfile,"%d-th: m= %d KZAA= %8.4lf KZCG1= %8.4lf KZCG2= %8.4lf\n",
	    i+1,m,ZData.KZAA,ZData.KZCG[0],ZData.KZCG[1]);

    AAData.zeta=0.0; AAData.V_zeta=0.0;
    CGData[0].zeta=0.0; CGData[0].V_zeta=0.0;
    CGData[1].zeta=0.0; CGData[1].V_zeta=0.0;
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
    SKZ[1]=ZData.KZCG[0];
    SKZ[2]=ZData.KZCG[1];
    //    printf("107\n");
    MPI_Send(SKZ, /*2*/3, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD);
    MPI_Recv(RKZ, /*2*/3, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD, status);
    //    printf("110\n");
    KZAAn=RKZ[0];
    KZCGn[0]=RKZ[1];
    KZCGn[1]=RKZ[2];

    MPI_Barrier(MPI_COMM_WORLD);

    CE_TACCM_2CG1AA(AAData.crd,CGData[0].crd,CGData[1].crd,ZData.Z,
		    (CData.numatom),ZData.numZ,
		    KZAAn,KZCGn[0],KZCGn[1],ZData.pairs,pi,
		    &EAAn_Xi,&ECGn_Xi[0],&ECGn_Xi[1],&EZn_Xi);

    SData[0]=EAAm_Xi;
    SData[1]=ECGm_Xi[0];
    SData[2]=ECGm_Xi[1];
    SData[3]=EZm_Xi;

    SData[4]=EAAn_Xi;
    SData[5]=ECGn_Xi[0];
    SData[6]=ECGn_Xi[1];
    SData[7]=EZn_Xi;

    //    SData[6]=ZData.KZAA;
    //    SData[7]=ZData.KZCG;
    //    printf("148\n");
    MPI_Send(SData, 8/*6*/, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD);
    MPI_Recv(RData, 8/*6*/, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD, status);
    //    printf("151\n");
    EAAn_Xj=RData[0];
    ECGn_Xj[0]=RData[1];
    ECGn_Xj[1]=RData[2];
    EZn_Xj=RData[3];

    EAAm_Xj=RData[4];
    ECGm_Xj[0]=RData[5];
    ECGm_Xj[1]=RData[6];
    EZm_Xj=RData[7];

    //    KZAAn=RData[6];
    //    KZCGn=RData[7];

    MPI_Barrier(MPI_COMM_WORLD);

    if (n>m) {
      //      delta=betaAA*((EAAm_Xj-EAAm_Xi)-(EAAn_Xi-EAAn_Xj))
      //      	+betaCG*((ECGm_Xj-ECGm_Xi)-(ECGn_Xi-ECGn_Xj))
      //      	+betaZ*((EZm_Xj-EZm_Xi)-(EZn_Xi-EZn_Xj));

      delta=betaAA*((EAAm_Xj-EAAm_Xi)-(EAAn_Xj-EAAn_Xi))
      	+betaCG[0]*((ECGm_Xj[0]-ECGm_Xi[0])-(ECGn_Xj[0]-ECGn_Xi[0]))
      	+betaCG[1]*((ECGm_Xj[1]-ECGm_Xi[1])-(ECGn_Xj[1]-ECGn_Xi[1]))
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
	printf("%d-th Acc %d-%d : /_\\ = %8.4lf \n EAAmi=%8.4lf ECG1mi=%8.4lf EC2Gmi=%8.4lf EZmi=%8.4lf EAAni=%8.4lf ECG1ni=%8.4lf ECG2ni=%8.4lf EZni=%8.4lf\n EAAmj=%8.4lf ECG1mj=%8.4lf ECG2mj=%8.4lf EZmj=%8.4lf EAAnj=%8.4lf ECG1nj=%8.4lf ECG2nj=%8.4lf EZnj=%8.4lf\n KZAAm=%8.4lf KZAAn=%8.4lf KZCG1m=%8.4lf KZCG2m=%8.4lf KZCG1n=%8.4lf KZCG2n=%8.4lf\n",
	       i+1,myrank,your_rank,myrank,delta,
	       EAAm_Xi,ECGm_Xi[0],ECGm_Xi[1],EZm_Xi,EAAn_Xi,ECGn_Xi[0],ECGn_Xi[1],EZn_Xi,
	       EAAm_Xj,ECGm_Xj[0],ECGm_Xj[1],EZm_Xj,EAAn_Xj,ECGn_Xj[0],ECGn_Xj[1],EZn_Xj,
	       ZData.KZAA,KZAAn,ZData.KZCG[0],ZData.KZCG[1],KZCGn[0],KZCGn[1]);
      ZData.KZAA=KZAAn;
      ZData.KZCG[0]=KZCGn[0];
      ZData.KZCG[1]=KZCGn[1];

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
	printf("%d-th Rej %d-%d : /_\\ = %8.4lf \n EAAmi=%8.4lf ECG1mi=%8.4lf ECG2mi=%8.4lf EZmi=%8.4lf EAAni=%8.4lf ECG1ni=%8.4lf ECG2ni=%8.4lf EZni=%8.4lf\n EAAmj=%8.4lf ECG1mj=%8.4lf ECG2mj=%8.4lf EZmj=%8.4lf EAAnj=%8.4lf ECG1nj=%8.4lf ECG2nj=%8.4lf EZnj=%8.4lf\n KZAAm=%8.4lf KZAAn=%8.4lf KZCG1m=%8.4lf KZCG2m=%8.4lf KZCG1n=%8.4lf KZCG2n=%8.4lf\n",
	       i+1,myrank,your_rank,myrank,delta,
	       EAAm_Xi,ECGm_Xi[0],ECGm_Xi[1],EZm_Xi,EAAn_Xi,ECGn_Xi[0],ECGn_Xi[1],EZn_Xi,
	       EAAm_Xj,ECGm_Xj[0],ECGm_Xj[1],EZm_Xj,EAAn_Xj,ECGn_Xj[0],ECGn_Xj[1],EZn_Xj,
	       ZData.KZAA,KZAAn,ZData.KZCG[0],ZData.KZCG[1],KZCGn[0],KZCGn[1]);
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
