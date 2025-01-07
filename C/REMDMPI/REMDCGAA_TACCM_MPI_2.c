
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "REMDCGAA_TACCM_MPI_2.h"

#define ON  0
#define OFF 1

int MPI_CGAATREM_TACCM_MD_pep_NHC_MP1998(int myrank, int num_procs,int tag, MPI_Status* status,
					 int numRE, int numEX,
					 struct AADataforREMD AAData,struct CGDataforREMD CGData,
					 struct TACCMDataforREMD ZData, struct AACGCommonDataforREMD CData,
					 double T0AA,double T0CG, double T0Z, int numstep, int interval, 
					 double dt,double dt2,double wdt2[3],double wdt4[3], int nc,
					 double UNITT, double k_B, double tau, double pi ) {
  int i,j,k,c,l=0;
  int your_rank;
  int count;

  double betaAA,betaCG,betaZ;

  double EAAm_Xi,EAAm_Xj,EAAn_Xi,EAAn_Xj;
  double ECGm_Xi,ECGm_Xj,ECGn_Xi,ECGn_Xj;
  double EZm_Xi,EZm_Xj,EZn_Xi,EZn_Xj;

  double *crdAA_receiv,*crdCG_receiv;
  double *Z_receiv;
  double *crd_temp,*vel_temp;
  double *Z_temp,*vel_Z_temp;

  double *SData,*RData;
  double delta=0.0;

  unsigned long init[4]={0x123, 0x234, 0x345, 0x456},length=4;
  init_by_array(init,length);

  betaAA=1.0/(k_B*T0AA);
  betaCG=1.0/(k_B*T0CG);
  betaZ=1.0/(k_B*T0Z);

  SData=(double *)gcemalloc(sizeof(double)*(3+(CData.numatom)*3*2+ZData.numZ));
  RData=(double *)gcemalloc(sizeof(double)*(3+(CData.numatom)*3*2+ZData.numZ));

  crdAA_receiv=(double *)gcemalloc(sizeof(double)*(CData.numatom)*3);
  crdCG_receiv=(double *)gcemalloc(sizeof(double)*(CData.numatom)*3);
  Z_receiv=(double *)gcemalloc(sizeof(double)*ZData.numZ);

  crd_temp=(double *)gcemalloc(sizeof(double)*(CData.numatom)*3);
  vel_temp=(double *)gcemalloc(sizeof(double)*(CData.numatom)*3);
  Z_temp=(double *)gcemalloc(sizeof(double)*(ZData.numZ));
  vel_Z_temp=(double *)gcemalloc(sizeof(double)*(ZData.numZ));

  //  printf("yes 59 in MPI\n");

  for (i=0;i<numEX;++i) {
    runTACCM_CGAA_MD_NHC_MP1998(// AA /////////////////////////////////////////////////////////
				AAData.crd[myrank],AAData.vel[myrank],
				&(AAData.zeta[myrank]),&(AAData.V_zeta[myrank]),AAData.Q,
				AAData.e[myrank],AAData.f[myrank],AAData.T[myrank],AAData.NfKT,
				&(AAData.avePE[myrank]),&(AAData.aveKE[myrank]),&(AAData.aveT[myrank]),
				&(AAData.varPE[myrank]),&(AAData.varKE[myrank]),&(AAData.varT[myrank]),
				AAData.nc_id_MCD[myrank],AAData.outputfile[myrank],
				// CG /////////////////////////////////////////////////////////
				CGData.crd[myrank],CGData.vel[myrank],
				&(CGData.zeta[myrank]),&(CGData.V_zeta[myrank]),CGData.Q,
				CGData.e_GOLM[myrank],CGData.T[myrank],CGData.NfKT,
				&(CGData.avePE[myrank]),&(CGData.aveKE[myrank]),&(CGData.aveT[myrank]),
				&(CGData.varPE[myrank]),&(CGData.varKE[myrank]),&(CGData.varT[myrank]),
				CGData.nc_id_MCD[myrank],CGData.outputfile[myrank],
				// Z  /////////////////////////////////////////////////////////
				ZData.Z[myrank],ZData.velZ[myrank],ZData.massZ,
				&(ZData.zetaZ[myrank]),&(ZData.V_zetaZ[myrank]),
				ZData.QZ,ZData.NfKTZ,ZData.T[myrank],
				ZData.numZ,ZData.KZAA[myrank],ZData.KZCG[myrank],ZData.pairs,
				&(ZData.avePEZ[myrank]),&(ZData.aveKEZ[myrank]),&(ZData.aveTZ[myrank]),
				&(ZData.varPEZ[myrank]),&(ZData.varKEZ[myrank]),&(ZData.varTZ[myrank]), 
				ZData.trjfileZ[myrank],ZData.trjfilThetaAA[myrank],ZData.trjfilThetaCG[myrank],
				// CM  /////////////////////////////////////////////////////////
				CData.mass,(CData.numatom),numstep,interval,&l,
				dt,dt2,wdt2,wdt4,nc,UNITT,k_B,pi,
				&EAAm_Xi,&ECGm_Xi,&EZm_Xi);
    //    printf("yes 88 in MPI\n");
    //    printf("EAAm_Xi=%8.4lf ECGm_Xi=%8.4lf EZm_Xi=%8.4lf\n",EAAm_Xi,ECGm_Xi,EZm_Xi);

    AAData.zeta[myrank]=0.0; AAData.V_zeta[myrank]=0.0;
    CGData.zeta[myrank]=0.0; CGData.V_zeta[myrank]=0.0;
    ZData.zetaZ[myrank]=0.0; ZData.V_zetaZ[myrank]=0.0;

    SData[0]=EAAm_Xi;
    SData[1]=ECGm_Xi;
    SData[2]=EZm_Xi;
    k=3;
    for (j=0;j<(CData.numatom)*3;++j)  {
      SData[k]=AAData.crd[myrank][j];      ++k;
    }
    for (j=0;j<(CData.numatom)*3;++j)  {
      SData[k]=CGData.crd[myrank][j];      ++k;
    }
    for (j=0;j<ZData.numZ;++j)	{
      SData[k]=ZData.Z[myrank][j];      ++k;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if(myrank%2==0 ) {
      if ( i%2 == 0 ) your_rank=myrank+1;
      else {
	if ( myrank==0 )  your_rank=numRE-1;
	else  your_rank=myrank-1;
      }      

      MPI_Send(SData, 3+(CData.numatom)*3*2+ZData.numZ, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD);
    }
    else {
      if ( i%2 == 0 ) your_rank=myrank-1;
      else {
	if ( myrank==numRE-1 )  your_rank=0;
	else  your_rank=myrank+1;
      }      

      MPI_Recv(RData, 3+(CData.numatom)*3*2+ZData.numZ, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD, status);

      EAAn_Xj=RData[0]; ECGn_Xj=RData[1]; EZn_Xj=RData[2];
      for (j=0;j<(CData.numatom)*3;++j) {
	crdAA_receiv[j]=RData[3+j];
	crdCG_receiv[j]=RData[3+(CData.numatom)*3+j];
      }
      for (j=0;j<ZData.numZ;++j) Z_receiv[j]=RData[3+(CData.numatom)*3*2+j];

      CE_TACCM_CGAA(crdAA_receiv,crdCG_receiv,Z_receiv,
      		    (CData.numatom),ZData.numZ,ZData.KZAA[myrank],ZData.KZCG[myrank],ZData.pairs,
      		    &EAAm_Xj,&ECGm_Xj,&EZm_Xj);

      CE_TACCM_CGAA(AAData.crd[myrank],CGData.crd[myrank],ZData.Z[myrank],
		    CData.numatom,ZData.numZ,ZData.KZAA[your_rank],ZData.KZCG[your_rank],ZData.pairs,
      		    &EAAn_Xi,&ECGn_Xi,&EZn_Xi);

      delta=betaAA*((EAAm_Xj-EAAm_Xi)-(EAAn_Xi-EAAn_Xj))
      	   +betaCG*((ECGm_Xj-ECGm_Xi)-(ECGn_Xi-ECGn_Xj))
      	   +betaZ*((EZm_Xj-EZm_Xi)-(EZn_Xi-EZn_Xj));

    printf("EAAm_Xi=%8.4lf ECGm_Xi=%8.4lf EZm_Xi=%8.4lf\n EAAm_Xj=%8.4lf ECGm_Xj=%8.4lf EZm_Xj=%8.4lf\n EAAn_Xi=%8.4lf ECGn_Xi=%8.4lf EZn_Xi=%8.4lf\n EAAn_Xj=%8.4lf ECGn_Xj=%8.4lf EZn_Xj=%8.4lf\n",
	   EAAm_Xi,ECGm_Xi,EZm_Xi,EAAm_Xj,ECGm_Xj,EZm_Xj,EAAn_Xi,ECGn_Xi,EZn_Xi,EAAn_Xj,ECGn_Xj,EZn_Xj);

      if((c=Metropolis(delta))==1) {
        for (j=0;j<(CData.numatom)*3;++j)  {
	  vel_temp[j]=AAData.vel[myrank][j];
	  AAData.vel[myrank][j]=AAData.vel[your_rank][j];
	  AAData.vel[your_rank][j]=vel_temp[j];
	  crd_temp[j]=AAData.crd[myrank][j];
	  AAData.crd[myrank][j]=AAData.crd[your_rank][j];
	  AAData.crd[your_rank][j]=crd_temp[j];

	  vel_temp[j]=CGData.vel[myrank][j];
	  CGData.vel[myrank][j]=CGData.vel[your_rank][j];
	  CGData.vel[your_rank][j]=vel_temp[j];	
	  crd_temp[j]=CGData.crd[myrank][j];
	  CGData.crd[myrank][j]=CGData.crd[your_rank][j];
	  CGData.crd[your_rank][j]=crd_temp[j];
	}

        for (j=0;j<ZData.numZ;++j)  {
	  vel_Z_temp[j]=ZData.velZ[myrank][j];
	  ZData.velZ[myrank][j]=ZData.velZ[your_rank][j];
	  ZData.velZ[your_rank][j]=vel_Z_temp[j];

	  Z_temp[j]=ZData.Z[myrank][j];
	  ZData.Z[myrank][j]=ZData.Z[your_rank][j];
	  ZData.Z[your_rank][j]=Z_temp[j];
	}

	printf("%d -th EX acc between %d-%d : /_\\ = %8.4lf \n",i+1,myrank,your_rank,myrank,delta);
      }
      else  printf("%d -th EX rej between %d-%d : /_\\ = %8.4lf \n",i+1,myrank,your_rank,myrank,delta);
    }
  }
  
  return 0;
}

void  CGAAREMDreadInputs(FILE *inputfile,int numatom,int numRE,
			 double **crdAA,double **velAA, double **crdCG,double **velCG,
			 double *KZAA, double *KZCG){
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
  State=AAINPF;
  Tflag=OFF;

  while ((c=getc(inputfile))!=-1){
    if (State==AAINPF) {
      if (c==' ' || c=='\n') {
	if (j>0) {
	  crdfilename[j]=NULL;
	  crdfile=efopen(crdfilename,"r");
	  getline(&line,&len,crdfile);
	  fscanf(crdfile,"%d",&d);
	  for (k=0;k<numatom*3;++k) fscanf(crdfile,"%lf",&crdAA[i][k]);
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
	  crdfilename[j]=NULL;
	  crdfile=efopen(crdfilename,"r");
	  getline(&line,&len,crdfile);
	  fscanf(crdfile,"%d",&d);
	  for (k=0;k<numatom*3;++k) fscanf(crdfile,"%lf",&crdCG[i][k]);
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
