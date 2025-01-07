
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "REMDCGAA_TACCM_MPI_2_Amber-CAGo-MB.h"
#include "REMD_functions.h"

#include "TACCM_CGAAMDrun_Amber_CAGo_MB.h"

#include "TACCM.h"

#include "TACCM_CGAAMDrun_test.h"

#include "MD_NHC_MP1996.h"
#include "MDrun.h"
#include "MD.h"

#include "FVDIHED.h"

#define ON  0
#define OFF 1

void  CGAAREMDreadInputs_Amber_CAGo_MB(FILE *inputfile,int numatom,int numRE,int myrank,
				       double *crdAA,double *velAA, double *crdCG,double *velCG,
				       double *KZAA, double *KZCG, struct AmberParmL ap_CG) {
  int i,j,k,l;
  int numCA;
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
	  if (i==myrank) {
	    numCA=0;
	    for (k=0;k<numatom;++k) {
	      if (strncmp(ap_CG.IGRAPH[k],"CA",2)==0) {
		for (l=0;l<3;++l) {
		  fscanf(crdfile,"%lf",&crdCG[numCA*3+l]);
		}
		++numCA;
	      }
	      else {
		for (l=0;l<3;++l) {
		  fscanf(crdfile,"%lf",&fdummy);
		}
	      }
	    }
	  }
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

////////////////////////////////////////////////////////////////////////////////////////////////

/*******************************************************************************************************************************/
/* double runTACCM_CGAA_MD_NHC_MP1998_Amber_CAGo_MB(// AA /////////////////////////////////////////////////////////	       */
/* 						 double *crdAA,double *velAA, 						       */
/* 						 double *zetaAA,double *V_zetaAA, double QAA,				       */
/* 						 struct potential e, struct force f, double TAA, double NfKTAA,		       */
/* 						 double *avePEAA, double *aveKEAA,double *aveTAA,			       */
/* 						 double *varPEAA, double *varKEAA,double *varTAA,			       */
/* 						 struct my_netcdf_out_id_MCD nc_id_MCDAA,  FILE *outputfileAA,		       */
/* 						 // CG /////////////////////////////////////////////////////////	       */
/* 						 double *crdCG,double *velCG, 						       */
/* 						 double *zetaCG,double *V_zetaCG, double QCG,				       */
/* 						 struct potential_GOLM_Clementi_MB e_CG,				       */
/* 						 double de, double d2,							       */
/* 						 double TCG, double NfKTCG,						       */
/* 						 double *avePECG, double *aveKECG,double *aveTCG,			       */
/* 						 double *varPECG, double *varKECG,double *varTCG,			       */
/* 						 struct my_netcdf_out_id_MCD nc_id_MCDCG,  FILE *outputfileCG,		       */
/* 						 // Z  /////////////////////////////////////////////////////////	       */
/* 						 double *Z,double *velZ,double massZ,					       */
/* 						 double *zetaZ,double *V_zetaZ,						       */
/* 						 double QZ,double NfKTZ,double TZ,					       */
/* 						 double KZAA, double KZCG,						       */
/* 						 int numZ_dih,int **pairs_dihe_AA,int **pairs_dihe_CG,			       */
/* 						 int numZ_ang,int **pairs_angl_AA,int **pairs_angl_CG,			       */
/* 						 int numZ_bon,int **pairs_bond_AA,int **pairs_bond_CG,			       */
/* 						 double *avePEZ, double *aveKEZ,double *aveTZ,				       */
/* 						 double *varPEZ, double *varKEZ,double *varTZ, 				       */
/* 						 FILE *trjfileZ, FILE *trjfilThetaAA, FILE *trjfilThetaCG,		       */
/* 						 // CM  /////////////////////////////////////////////////////////	       */
/* 						 double *mass,double *massCA,						       */
/* 						 int numatom, int numCAatom, 						       */
/* 						 int numstep, int interval,int *l,					       */
/* 						 double dt,double dt2,double wdt2[3] ,double wdt4[3] ,int nc,		       */
/* 						 double UNITT, double k_B,double pi,					       */
/* 						 double *PEZAA, double *PEZCG, double *PEZ) {				       */
/*   int i,j,k;														       */
/* 															       */
/*   double PEAA=0.0,KEAA=0.0,EtAA,PEvAA,KEvAA;										       */
/*   double PECG=0.0,KECG=0.0,EtCG,PEvCG,KEvCG;										       */
/*   double KEZ,KEvZ,PEvZ,EtZ;												       */
/* 															       */
/*   double *thetaAA,*thetaCG,*frcZ,*fAA,*fCG;										       */
/*   double summass,COM[3],crd_nc[MAXATOM][3];										       */
/* 															       */
/*   struct force f_CG;													       */
/*   struct potential e_CG2;												       */
/* 															       */
/*   //  printf("nb=%5.3d na=%5.3d nd=%5.3d natom=%5.3d\n");								       */
/*   //  printf("here is line 70 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");						       */
/*   //  *aveKEZ=0.0;													       */
/*   //  *varPEZ=0.0;													       */
/*   //  *varKEZ=0.0;													       */
/* 															       */
/*   //  *aveTAA=0.0;													       */
/*   //  *varTAA=0.0;													       */
/* 															       */
/*   //  *aveTCG=0.0;													       */
/*   //  *varTCG=0.0;													       */
/* 															       */
/*   //  *aveTZ=0.0;													       */
/*   //  *varTZ=0.0;													       */
/* 															       */
/*   ffL_calcffandforce_AA(crdAA,numatom,&e,&f);									       */
/*   GOLM_Clementi_MB_ff_calcff(crdCG,numCAatom,de,d2,&e_CG);								       */
/* 															       */
/*   fAA=(double *)gcemalloc(sizeof(double)*(numZ_dih+numZ_ang+numZ_bon));						       */
/*   fCG=(double *)gcemalloc(sizeof(double)*(numZ_dih+numZ_ang+numZ_bon));						       */
/*   frcZ=(double *)gcemalloc(sizeof(double)*(numZ_dih+numZ_ang+numZ_bon));						       */
/* 															       */
/*   thetaAA=(double *)gcemalloc(sizeof(double)*(numZ_dih+numZ_ang+numZ_bon));						       */
/*   thetaCG=(double *)gcemalloc(sizeof(double)*(numZ_dih+numZ_ang+numZ_bon));						       */
/* 															       */
/*   TACCM_CTheta_Amber_CAGo_MB(crdAA,numatom,thetaAA,									       */
/* 			     numZ_dih,pairs_dihe_AA,									       */
/* 			     numZ_ang,pairs_angl_AA,									       */
/* 			     numZ_bon,pairs_bond_AA,									       */
/* 			     pi);											       */
/*   TACCM_CTheta_Amber_CAGo_MB(crdCG,numCAatom,thetaCG,								       */
/* 			     numZ_dih,pairs_dihe_CG,									       */
/* 			     numZ_ang,pairs_angl_CG,									       */
/* 			     numZ_bon,pairs_bond_CG,									       */
/* 			     pi);											       */
/* 															       */
/*   TACCM_calc_eff_FF_Z(Z,numZ_dih+numZ_ang+numZ_bon,									       */
/* 		      thetaAA,KZAA,fAA,pi);										       */
/*   TACCM_calc_eff_FF_Z(Z,numZ_dih+numZ_ang+numZ_bon,									       */
/* 		      thetaCG,KZCG,fCG,pi);										       */
/* 															       */
/*   for (i=0;i<numZ_dih+numZ_ang+numZ_bon;++i) frcZ[i]=fAA[i]+fCG[i];							       */
/* 															       */
/*   //  printf("here is line 110 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");						       */
/* 															       */
/*   for (i=0;i<numstep;++i) {												       */
/*     //    printf("here is line 113 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");						       */
/*     //    printf("\n numstep=%5.3d\n",i);										       */
/*     //    printf("nb=%5.3d na=%5.3d nd=%5.3d \n",numZ_bon,numZ_ang,numZ_dih);					       */
/*     //    printf("natom=%5.3d \n",numatom);										       */
/* 															       */
/*     TACCM_CTheta_Amber_CAGo_MB(crdAA,numatom,thetaAA,								       */
/*     			       numZ_dih,pairs_dihe_AA,									       */
/*     			       numZ_ang,pairs_angl_AA,									       */
/*     			       numZ_bon,pairs_bond_AA,									       */
/*     			       pi);											       */
/* 															       */
/*     //    printf("here is line 119 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");						       */
/* 															       */
/*     TACCM_CTheta_Amber_CAGo_MB(crdCG,numCAatom,thetaCG,								       */
/*     			       numZ_dih,pairs_dihe_CG,									       */
/*     			       numZ_ang,pairs_angl_CG,									       */
/*     			       numZ_bon,pairs_bond_CG,									       */
/*     			       pi);											       */
/* 															       */
/*     //    printf("here is line 129 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");						       */
/* 															       */
/*     PEAA=TACCM_MD_Propagetor_NH_MP1998_AA_Amber_CAGo_MB(crdAA,velAA,mass,zetaAA,V_zetaAA,				       */
/*     							QAA,NfKTAA,numatom,&KEAA,&KEvAA,&PEvAA,				       */
/*     							dt,dt2,nc,wdt4,wdt2,						       */
/*     							&e,&f,Z,numZ_dih,numZ_ang,numZ_bon,				       */
/*     							thetaAA,KZAA,							       */
/*     							pairs_dihe_AA,pairs_angl_AA,pairs_bond_AA,			       */
/*     							PEZAA,pi);							       */
/* 															       */
/*     //    printf("here is line 138 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");						       */
/* 															       */
/*     PECG=TACCM_MD_Propagetor_NH_MP1998_CG_Amber_CAGo_MB(crdCG,velCG,massCA,zetaCG,V_zetaCG,				       */
/*     							QCG,NfKTCG,numCAatom,&KECG,&KEvCG,&PEvCG,			       */
/*     							dt,dt2,nc,wdt4,wdt2,						       */
/*     							&e_CG,de,d2,Z,numZ_dih,numZ_ang,numZ_bon,			       */
/*     							thetaCG,KZCG,							       */
/*     							pairs_dihe_CG,pairs_angl_CG,pairs_bond_CG,			       */
/*     							PEZCG,pi);							       */
/* 															       */
/*     //    printf("here is line 147 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");						       */
/* 															       */
/*     TACCM_CGAA_MD_Propagetor_NH_MP1998_Z_2(Z,velZ,massZ,thetaAA,thetaCG,zetaZ,V_zetaZ,				       */
/*     					   QZ,NfKTZ,numZ_dih+numZ_ang+numZ_bon,						       */
/*     					   &KEZ,&KEvZ,&PEvZ,								       */
/*     					   dt,dt2,nc,wdt4,wdt2,KZAA,KZCG,						       */
/*     					   PEZAA,PEZCG,PEZ,frcZ,pi);							       */
/* 															       */
/*     //    printf("here is line 156 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");						       */
/* 															       */
/*     if (i%interval==0) {												       */
/*     															       */
/*       //      printf("here is line 165 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");					       */
/*     															       */
/*       KEAA=KEAA/UNITT;     TAA=KEAA/((3*numatom)*k_B)*2.0;								       */
/*       PEvAA=PEvAA/UNITT;   KEvAA=KEvAA/UNITT;									       */
/*     															       */
/*       KECG=KECG/UNITT;     TCG=KECG/((3*numCAatom)*k_B)*2.0;								       */
/*       PEvCG=PEvCG/UNITT;  KEvCG=KEvCG/UNITT;										       */
/*     															       */
/*       PEAA=e.p_d_t+e.p_a_t+e.p_b_t;											       */
/*       EtAA=PEAA+KEAA+PEvAA+KEvAA;											       */
/*       fprintf(outputfileAA,"%d %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e "							       */
/*     	                    ,i+1,PEAA,KEAA,KEvAA,PEvAA,EtAA,TAA);							       */
/*     															       */
/*       PECG=e_CG.p_MB;												       */
/*       EtCG=PECG+KECG+PEvCG+KEvCG;											       */
/*       fprintf(outputfileCG,"%d %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e "							       */
/*     	                  ,i+1,PECG,KECG,KEvCG,PEvCG,EtCG,TCG);								       */
/*     															       */
/*       //      *avePEAA=(i*(*avePEAA)+PEAA)/(i+1); *varPEAA=(i*(*varPEAA)+PEAA*PEAA)/(i+1);				       */
/*       //      *aveKEAA=(i*(*aveKEAA)+KEAA)/(i+1); *varKEAA=(i*(*varKEAA)+KEAA*KEAA)/(i+1);				       */
/*       //      *aveTAA=(i*(*aveTAA)+TAA)/(i+1);  *varTAA=(i*(*varTAA)+TAA*TAA)/(i+1);					       */
/*     															       */
/*       //      *avePECG=(i*(*avePECG)+PECG)/(i+1); *varPECG=(i*(*varPECG)+PECG*PECG)/(i+1);				       */
/*       //      *aveKECG=(i*(*aveKECG)+KECG)/(i+1); *varKECG=(i*(*varKECG)+KECG*KECG)/(i+1);				       */
/*       //      *aveTCG=(i*(*aveTCG)+TCG)/(i+1);  *varTCG=(i*(*varTCG)+TCG*TCG)/(i+1);					       */
/*     															       */
/*       summass=0.0; for (j=0;j<numatom;++j) summass+=mass[j];								       */
/*       for (j=0;j<3;++j) COM[j]=0.0;											       */
/*       for (j=0;j<numatom;++j)  for (k=0;k<3;++k) COM[k]+=mass[j]*crdAA[j*3+k]/summass;				       */
/*       for (j=0;j<numatom;++j)  for (k=0;k<3;++k) crdAA[j*3+k]-=COM[k];						       */
/*       for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crdAA[j*3+k];						       */
/*       myncL_put_crd_MCD(nc_id_MCDAA,*l,crd_nc);									       */
/*     															       */
/*       summass=0.0; for (j=0;j<numCAatom;++j) summass+=massCA[j];							       */
/*       for (j=0;j<3;++j) COM[j]=0.0;											       */
/*       for (j=0;j<numCAatom;++j)  for (k=0;k<3;++k) COM[k]+=massCA[j]*crdCG[j*3+k]/summass;				       */
/*       for (j=0;j<numCAatom;++j)  for (k=0;k<3;++k) crdCG[j*3+k]-=COM[k];						       */
/*       for (j=0;j<numCAatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crdCG[j*3+k];						       */
/*       for (j=0;j<numCAatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crdCG[j*3+k];						       */
/*       myncL_put_crd_MCD(nc_id_MCDCG,*l,crd_nc);									       */
/*       ++(*l);													       */
/*     															       */
/*       ///////////////// TACCM //////////////////////									       */
/*       KEZ=KEZ/UNITT;      TZ=KEZ/((numZ_dih+numZ_ang+numZ_bon)*k_B)*2.0;						       */
/*     															       */
/*       PEvZ=PEvZ/UNITT;      KEvZ=KEvZ/UNITT;										       */
/*     															       */
/*       EtZ=*PEZ+KEZ+PEvZ+KEvZ;											       */
/*       fprintf(outputfileAA,"%d %e %e %e %e %e %e %e\n",i+1,*PEZ,KEZ,KEvZ,PEvZ,EtZ,TZ);				       */
/*     															       */
/*       for (j=0;j<(numZ_dih+numZ_ang+numZ_bon);++j) fprintf(trjfileZ,"%e ",Z[j]);					       */
/*       fprintf(trjfileZ,"\n");											       */
/*       for (j=0;j<(numZ_dih+numZ_ang+numZ_bon);++j) fprintf(trjfilThetaAA,"%e ",thetaAA[j]);				       */
/*       fprintf(trjfilThetaAA,"\n");											       */
/*       for (j=0;j<(numZ_dih+numZ_ang+numZ_bon);++j) fprintf(trjfilThetaCG,"%e ",thetaCG[j]);				       */
/*       fprintf(trjfilThetaCG,"\n");											       */
/*       ///////////////// TACCM //////////////////////									       */
/*     }														       */
/*   }															       */
/*   //  printf("here is line 246 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");						       */
/* 															       */
/*   return *PEZ;													       */
/* 															       */
/* 															       */
/* }															       */
/* 															       */
/* double TACCM_MD_Propagetor_NH_MP1998_AA_Amber_CAGo_MB(double *crd,double *vel,double *mass,				       */
/* 						      double *zeta,double *V_zeta,double Q,				       */
/* 						      double NfKT,int numatom,double *KE,double *KEv,double *PEv,	       */
/* 						      double dt,double dt2,int nc,double wdt4[3],double wdt2[3],	       */
/* 						      struct potential *e, struct force *f,				       */
/* 						      double *Z, int numZ_dih,int numZ_ang, int numZ_bon,		       */
/* 						      double *theta,double Kapa,					       */
/* 						      int **pairs_dihe_AA, int **pairs_angl_AA, int **pairs_bond_AA,	       */
/* 						      double *PEZ,double pi) {						       */
/*   int i,j,k;														       */
/*   double *frc;													       */
/*   double **frcZ;													       */
/* 															       */
/*   frc=(double *)gcemalloc(sizeof(double)*numatom*3);									       */
/*   frcZ=(double **)gcemalloc(sizeof(double *)*numatom);								       */
/*   for (i=0;i<numatom;++i) frcZ[i]=(double *)gcemalloc(sizeof(double)*3);						       */
/* 															       */
/*   TACCM_CTheta_Amber_CAGo_MB(crd,numatom,theta,									       */
/* 		       numZ_dih,pairs_dihe_AA,										       */
/* 		       numZ_ang,pairs_angl_AA,										       */
/* 		       numZ_bon,pairs_bond_AA,										       */
/* 		       pi);												       */
/*   *PEZ=TACCM_calc_eff_FF_MD_Amber_CAGo_MB(crd,numatom,theta,Z,							       */
/* 					  numZ_dih,numZ_ang,numZ_bon,							       */
/* 					  Kapa,frcZ,									       */
/* 					  pairs_dihe_AA,pairs_angl_AA,pairs_bond_AA,pi);				       */
/* 															       */
/*   for (i=0;i<numatom;++i)												       */
/*     for (j=0;j<3;++j) 												       */
/*       frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+frcZ[i][j];					       */
/* 															       */
/*   MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);				       */
/* 															       */
/*   //////////////////////////////////////////////////////////////////////////////////////				       */
/*   for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];					       */
/*   for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];						       */
/* 															       */
/*   ffL_calcffandforce_AA(crd,numatom,e,f);										       */
/*   TACCM_CTheta_Amber_CAGo_MB(crd,numatom,theta,									       */
/* 		       numZ_dih,pairs_dihe_AA,										       */
/* 		       numZ_ang,pairs_angl_AA,										       */
/* 		       numZ_bon,pairs_bond_AA,										       */
/* 		       pi);												       */
/*   *PEZ=TACCM_calc_eff_FF_MD_Amber_CAGo_MB(crd,numatom,theta,Z,							       */
/* 					  numZ_dih,numZ_ang,numZ_bon,							       */
/* 					  Kapa,frcZ,									       */
/* 					  pairs_dihe_AA,pairs_angl_AA,pairs_bond_AA,pi);				       */
/* 															       */
/*   for (i=0;i<numatom;++i)												       */
/*     for (j=0;j<3;++j) 												       */
/*       frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+frcZ[i][j];					       */
/*   for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];					       */
/*   //////////////////////////////////////////////////////////////////////////////////////				       */
/* 															       */
/*   MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);				       */
/* 															       */
/*   *KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) *KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];				       */
/* 															       */
/*   *KEv=0.5*Q*(*V_zeta)*(*V_zeta);											       */
/*   *PEv=NfKT*(*zeta);													       */
/* 															       */
/*   return 0.0;													       */
/* }															       */
/* 															       */
/* double TACCM_MD_Propagetor_NH_MP1998_CG_Amber_CAGo_MB(double *crd,double *vel,double *mass,				       */
/* 						      double *zeta,double *V_zeta,double Q,				       */
/* 						      double NfKT,int numCAatom,double *KE,double *KEv,double *PEv,	       */
/* 						      double dt,double dt2,int nc,double wdt4[3],double wdt2[3],	       */
/* 						      struct potential_GOLM_Clementi_MB *e_CG,				       */
/* 						      double de, double d2,						       */
/* 						      double *Z, int numZ_dih,int numZ_ang, int numZ_bon,		       */
/* 						      double *theta,double Kapa,					       */
/* 						      int **pairs_dihe_CG, int **pairs_angl_CG, int **pairs_bond_CG,	       */
/* 						      double *PEZ,double pi) {						       */
/*   int i,j,k;														       */
/*   double *frc;													       */
/*   double **frcZ;													       */
/* 															       */
/*   frc=(double *)gcemalloc(sizeof(double)*numCAatom*3);								       */
/*   frcZ=(double **)gcemalloc(sizeof(double *)*numCAatom);								       */
/*   for (i=0;i<numCAatom;++i) frcZ[i]=(double *)gcemalloc(sizeof(double)*3);						       */
/* 															       */
/*   TACCM_CTheta_Amber_CAGo_MB(crd,numCAatom,theta,									       */
/* 		       numZ_dih,pairs_dihe_CG,										       */
/* 		       numZ_ang,pairs_angl_CG,										       */
/* 		       numZ_bon,pairs_bond_CG,										       */
/* 		       pi);												       */
/*   *PEZ=TACCM_calc_eff_FF_MD_Amber_CAGo_MB(crd,numCAatom,theta,Z,							       */
/* 					  numZ_dih,numZ_ang,numZ_bon,							       */
/* 					  Kapa,frcZ,									       */
/* 					  pairs_dihe_CG,pairs_angl_CG,pairs_bond_CG,pi);				       */
/* 															       */
/*   for (i=0;i<numCAatom;++i)												       */
/*     for (j=0;j<3;++j) 												       */
/*       frc[i*3+j]=(*e_CG).f_MB[i][j]+frcZ[i][j];									       */
/* 															       */
/*   MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numCAatom,nc,wdt4,wdt2);				       */
/* 															       */
/*   //////////////////////////////////////////////////////////////////////////////////////				       */
/*   for (i=0;i<numCAatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];					       */
/*   for (i=0;i<numCAatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];						       */
/* 															       */
/*   GOLM_Clementi_MB_ff_calcff(crd,numCAatom,de,d2,e_CG);								       */
/*   TACCM_CTheta_Amber_CAGo_MB(crd,numCAatom,theta,									       */
/* 		       numZ_dih,pairs_dihe_CG,										       */
/* 		       numZ_ang,pairs_angl_CG,										       */
/* 		       numZ_bon,pairs_bond_CG,										       */
/* 		       pi);												       */
/*   *PEZ=TACCM_calc_eff_FF_MD_Amber_CAGo_MB(crd,numCAatom,theta,Z,							       */
/* 					  numZ_dih,numZ_ang,numZ_bon,							       */
/* 					  Kapa,frcZ,									       */
/* 					  pairs_dihe_CG,pairs_angl_CG,pairs_bond_CG,pi);				       */
/* 															       */
/*   for (i=0;i<numCAatom;++i)												       */
/*     for (j=0;j<3;++j) 												       */
/*       frc[i*3+j]=(*e_CG).f_MB[i][j]+frcZ[i][j];									       */
/*   for (i=0;i<numCAatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];					       */
/*   for (i=0;i<numCAatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];					       */
/*   //////////////////////////////////////////////////////////////////////////////////////				       */
/* 															       */
/*   MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numCAatom,nc,wdt4,wdt2);				       */
/* 															       */
/*   *KE=0.0; for (i=0;i<numCAatom;++i) for (j=0;j<3;++j) *KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];			       */
/* 															       */
/*   *KEv=0.5*Q*(*V_zeta)*(*V_zeta);											       */
/*   *PEv=NfKT*(*zeta);													       */
/* 															       */
/*   return 0.0;													       */
/* }															       */
/* 															       */
/* double TACCM_CTheta_Amber_CAGo_MB(double *crd,int numatom,double *theta, 						       */
/* 				  int numdihe, int **pairs_dih_AA,							       */
/* 				  int numangl, int **pairs_ang_AA,							       */
/* 				  int numbond, int **pairs_bon_AA, 							       */
/* 				  double pi){										       */
/*   int i,j,k,l;													       */
/*   int ii,jj,kk,ll;													       */
/* 															       */
/*   double lenij,lenkj;												       */
/*   double cosijk,angijk;												       */
/* 															       */
/*   double m[3],n[3],m_n[3],n_n[3],lm,ln;										       */
/*   double vij[3],vkj[3],vkl[3];											       */
/*   double lkj;													       */
/*   double vijvkj,vklvkj;												       */
/* 															       */
/*   double atom[4][3];													       */
/*   double angl,dihed;													       */
/* 															       */
/*   //  printf("here is line 381 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");						       */
/* 															       */
/*   /\*****************************************************************\/						       */
/*   /\* for (i=0;i<numbond;++i) {					   *\/						       */
/*   /\*   ii=pairs_bon_AA[i][0]/\\*-1*\\/;				   *\/						       */
/*   /\*   jj=pairs_bon_AA[i][1]/\\*-1*\\/;				   *\/						       */
/*   /\*   for (j=0;j<3;++j) {					   *\/							       */
/*   /\*     atom[0][j]=crd[ii*3+j];				   *\/							       */
/*   /\*     atom[1][j]=crd[jj*3+j];				   *\/							       */
/*   /\*   }							   *\/							       */
/*   /\* 								   *\/						       */
/*   /\*   lenij=0.0;						   *\/							       */
/*   /\*   for (j=0;j<3;++j) {					   *\/							       */
/*   /\*     lenij += (atom[0][j]-atom[1][j])*(atom[0][j]-atom[1][j]); *\/						       */
/*   /\*   }							   *\/							       */
/*   /\*   lenij=sqrt(lenij);					   *\/							       */
/*   /\*   theta[i]=lenij;						   *\/						       */
/*   /\* }								   *\/						       */
/*   /\*****************************************************************\/						       */
/* 															       */
/*   //  printf("here is line 397 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");						       */
/*   //  printf("numangl=%5.3d\n",numangl);										       */
/*   ////////////////////////////////////////////////////////////							       */
/* 															       */
/*   /\**************************************\/										       */
/*   /\* for (i=0;i<numangl;++i) {	        *\/									       */
/*   /\*   ii=pairs_ang_AA[i][0]/\\*-1*\\/;   *\/									       */
/*   /\*   jj=pairs_ang_AA[i][1]/\\*-1*\\/;   *\/									       */
/*   /\*   kk=pairs_ang_AA[i][2]/\\*-1*\\/;   *\/									       */
/*   /\*   for (j=0;j<3;++j) {	        *\/										       */
/*   /\*     atom[0][j]=crd[ii*3+j];        *\/										       */
/*   /\*     atom[1][j]=crd[jj*3+j];        *\/										       */
/*   /\*     atom[2][j]=crd[kk*3+j];        *\/										       */
/*   /\*   }			        *\/										       */
/*   /\* 				        *\/									       */
/*   /\*   lenij = len(atom[0],atom[1]);    *\/										       */
/*   /\*   lenkj = len(atom[2],atom[1]);    *\/										       */
/*   /\*   for (j=0;j<3;++j) {	        *\/										       */
/*   /\*     vij[j]=atom[1][j]-atom[0][j];  *\/										       */
/*   /\*     vkj[j]=atom[1][j]-atom[2][j];  *\/										       */
/*   /\*   }			        *\/										       */
/*   /\*   cosijk=inprod(vij,vkj,3);        *\/										       */
/*   /\*   cosijk=cosijk/lenij/lenkj;       *\/										       */
/*   /\*   if (cosijk>=1.0) angijk=0.0;     *\/										       */
/*   /\*   else if (cosijk<=0.0) angijk=pi; *\/										       */
/*   /\*   else  angijk = acos(cosijk);     *\/										       */
/*   /\* 				        *\/									       */
/*   /\*   theta[i+numbond]=lenij;	        *\/									       */
/*   /\* }				        *\/									       */
/*   /\**************************************\/										       */
/* 															       */
/*   //  printf("here is line 426 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");						       */
/*   ////////////////////////////////////////////////////////////							       */
/* 															       */
/*   /\********************************************\/									       */
/*   /\* for (i=0;i<numdihe;++i) {		      *\/								       */
/*   /\*   ii=pairs_dih_AA[i][0]/\\*-1*\\/;	      *\/								       */
/*   /\*   jj=pairs_dih_AA[i][1]/\\*-1*\\/;	      *\/								       */
/*   /\*   kk=pairs_dih_AA[i][2]/\\*-1*\\/;	      *\/								       */
/*   /\*   ll=pairs_dih_AA[i][3]/\\*-1*\\/;	      *\/								       */
/*   /\* 					      *\/								       */
/*   /\*   for (j=0;j<3;++j) {		      *\/									       */
/*   /\*     atom[0][j]=crd[ii*3+j];	      *\/									       */
/*   /\*     atom[1][j]=crd[jj*3+j];	      *\/									       */
/*   /\*     atom[2][j]=crd[kk*3+j];	      *\/									       */
/*   /\*     atom[3][j]=crd[ll*3+j];	      *\/									       */
/*   /\*   }				      *\/									       */
/*   /\* 					      *\/								       */
/*   /\*   for (j=0;j<3;++j) {		      *\/									       */
/*   /\*     vij[j] = atom[1][j]-atom[0][j];      *\/									       */
/*   /\*     vkj[j] = atom[1][j]-atom[2][j];      *\/									       */
/*   /\*     vkl[j] = atom[3][j]-atom[2][j];      *\/									       */
/*   /\*   }				      *\/									       */
/*   /\*   lkj=sqrt(inprod(vkj,vkj,3));	      *\/									       */
/*   /\* 					      *\/								       */
/*   /\*   outprod(vij,vkj,m);		      *\/									       */
/*   /\*   outprod(vkj,vkl,n);		      *\/									       */
/*   /\*   lm=sqrt(inprod(m,m,3));		      *\/								       */
/*   /\*   ln=sqrt(inprod(n,n,3));		      *\/								       */
/*   /\*   for (j=0;j<3;++j) {		      *\/									       */
/*   /\*     m_n[j]=m[j]/lm;		      *\/									       */
/*   /\*     n_n[j]=n[j]/ln;		      *\/									       */
/*   /\*   }				      *\/									       */
/*   /\* 					      *\/								       */
/*   /\*   dihed=inprod(m_n,n_n,3);		      *\/								       */
/*   /\*   if (dihed>=1.0)			      *\/								       */
/*   /\*     dihed=0.0;			      *\/									       */
/*   /\*   else if (dihed<=-1.0)		      *\/								       */
/*   /\*     dihed=pi;			      *\/									       */
/*   /\*   else				      *\/									       */
/*   /\*     dihed=acos(dihed);		      *\/									       */
/*   /\*   if (inprod(vij,n,3)>0) dihed=-dihed;   *\/									       */
/*   /\*   if (dihed<-1.0*pi) dihed=2.0*pi+dihed; *\/									       */
/*   /\*   if (dihed>pi) dihed=-2.0*pi+dihed;     *\/									       */
/*   /\* 					      *\/								       */
/*   /\*   theta[i+numbond+numangl]=dihed;	      *\/								       */
/*   /\* }					      *\/								       */
/*   /\********************************************\/									       */
/* 															       */
/*   //  printf("here is line 471 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");						       */
/* 															       */
/*   return 0.0;													       */
/* }															       */
/* 															       */
/* double TACCM_calc_eff_FF_MD_Amber_CAGo_MB(double *crd,int numatom, double *theta, double *Z,  			       */
/* 					  int numZ_dih,int numZ_ang,int numZ_bon,					       */
/* 					  double Kapa, double **frcZ,							       */
/* 					  int **pairs_dih,int **pairs_ang,int **pairs_bon,				       */
/* 					  double pi){									       */
/*   int i,j,k,l;													       */
/* 															       */
/*   int ii,jj,kk,ll;													       */
/* 															       */
/*   double f;														       */
/* 															       */
/*   double lenij,lenkj;												       */
/*   double vij[3],vkj[3];												       */
/*   double cosijk,angijk;												       */
/*   double f1,f2;													       */
/* 															       */
/*   double atom[4][3];													       */
/* 															       */
/*   double *dvdpsi;													       */
/*   double delta;													       */
/* 															       */
/*   double PE=0.0;													       */
/* 															       */
/*   int **pairs_temp;													       */
/* 															       */
/*   /////////////////////////////////////////////////////////////////////////////////////				       */
/*   dvdpsi=(double *)gcemalloc(sizeof(double)*numZ_dih);								       */
/*   pairs_temp=(int **)gcemalloc(sizeof(int *)*numZ_dih);								       */
/*   for (i=0;i<numZ_dih;++i) pairs_temp[i]=(int *)gcemalloc(sizeof(int)*4);						       */
/* 															       */
/*   for (i=0;i<numZ_dih;++i) for (j=0;j<4;++j) pairs_temp[i][j]=pairs_dih[i][j]/\*-1*\/;				       */
/* 															       */
/*   for (i=0;i<numZ_dih;++i) {												       */
/*     if ((delta=Z[i]-theta[i+numZ_bon+numZ_ang])>pi) delta-=2.0*pi;							       */
/*     else if ((delta=Z[i]-theta[i+numZ_bon+numZ_ang])<-1.0*pi) delta+=2.0*pi;						       */
/*     dvdpsi[i]=-Kapa*delta;												       */
/*     PE+=0.5*Kapa*delta*delta;											       */
/*   }															       */
/* 															       */
/*   FVDIHED_force_dihed(crd,numatom,frcZ,pairs_temp,dvdpsi,numZ_dih);							       */
/*   /////////////////////////////////////////////////////////////////////////////////////////				       */
/*   //  printf("here is line 512 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");						       */
/* 															       */
/*   for (i=0;i<numZ_ang;++i) {												       */
/*     ii=pairs_ang[i][0]/\*-1*\/;											       */
/*     jj=pairs_ang[i][1]/\*-1*\/;											       */
/*     kk=pairs_ang[i][2]/\*-1*\/;											       */
/* 															       */
/*     for (j=0;j<3;++j) {												       */
/*       atom[0][j]=crd[ii*3+j];											       */
/*       atom[1][j]=crd[jj*3+j];											       */
/*       atom[2][j]=crd[kk*3+j];											       */
/*     }														       */
/* 															       */
/*     lenij = len(atom[0],atom[1]);											       */
/*     lenkj = len(atom[2],atom[1]);											       */
/*     for (j=0;j<3;++j) {												       */
/*       vij[j]=atom[1][j]-atom[0][j];											       */
/*       vkj[j]=atom[1][j]-atom[2][j];											       */
/*     }														       */
/*     cosijk=inprod(vij,vkj,3);											       */
/*     cosijk=cosijk/lenij/lenkj;											       */
/*     if (cosijk>=1.0) angijk=0.0;											       */
/*     else if (cosijk<=0.0) angijk=pi;											       */
/*     else  angijk = acos(cosijk);											       */
/* 															       */
/*     angijk = ang(atom[0],atom[1],atom[2]);										       */
/* 															       */
/*     PE += Kapa*(angijk-theta[i+numZ_bon])*(angijk-theta[i+numZ_bon]);						       */
/* 															       */
/*     for (j=0;j<3;++j) {												       */
/*       f1 = -2.0*Kapa*(angijk-theta[i+numZ_bon])/(lenij*sin(angijk))*(vkj[j]/lenkj-cosijk*vij[j]/lenij)		       */
/* 	*4.184070*100.0;												       */
/*       f2 = -2.0*Kapa*(angijk-theta[i+numZ_bon])/(lenkj*sin(angijk))*(vij[j]/lenij-cosijk*vkj[j]/lenkj)		       */
/* 	*4.184070*100.0;												       */
/* 															       */
/*       frcZ[ii][j] += f1;												       */
/*       frcZ[jj][j] += f2;												       */
/*       frcZ[kk][j] += -f1-f2;												       */
/*     }														       */
/*   }															       */
/*   /////////////////////////////////////////////////////////////////////////////////////////				       */
/*   //  printf("here is line 552 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");						       */
/* 															       */
/*   for (i=0;i<numZ_bon;++i) {												       */
/*     ii=pairs_bon[i][0]/\*-1*\/;											       */
/*     jj=pairs_bon[i][1]/\*-1*\/;											       */
/* 															       */
/*     for (j=0;j<3;++j) {												       */
/*       atom[0][j]=crd[ii*3+j];											       */
/*       atom[1][j]=crd[jj*3+j];											       */
/*     }														       */
/*   															       */
/*     lenij = len(atom[0],atom[1]);											       */
/*     PE+=Kapa*(lenij-theta[i])*(lenij-theta[i]);									       */
/* 															       */
/*     for (j=0;j<3;++j) {												       */
/*       f = 2.0*Kapa*(lenij-theta[i])*(atom[1][j]-atom[0][j])/lenij*4.184070*100.0;					       */
/*       frcZ[ii][j] += f;												       */
/*       frcZ[jj][j] += -f;												       */
/*     }														       */
/*   }															       */
/*   //  printf("here is line 571 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");						       */
/* 															       */
/*   return PE;														       */
/* }															       */
/* 															       */
/* double CE_TACCMb_CGAA_Amber_CAGo_MB(double *crdAA,double *crdCG,double *Z, int numatom,int numCAatom,		       */
/* 				    int numZ_dih, int **pairs_dihed_AA, int **pairs_dihed_CG,				       */
/* 				    int numZ_ang, int **pairs_angle_AA, int **pairs_angle_CG,				       */
/* 				    int numZ_bon, int **pairs_bond_AA,  int **pairs_bond_CG,				       */
/* 				    double KZAA,double KZCG, double pi,							       */
/* 				    double *EAA,double *ECG,double *EZ) {						       */
/*   int i,j;														       */
/*   int numZ;														       */
/*   double *thetaAA,*thetaCG;												       */
/*   double delta;													       */
/* 															       */
/*   thetaAA=(double *)gcemalloc(sizeof(double)*numZ);									       */
/*   thetaCG=(double *)gcemalloc(sizeof(double)*numZ);									       */
/*   //  printf("here is line 614 on TACCM_CGAAMDrun_Amber_CAGo_MB.c\n");						       */
/*   //  TACCM_CTheta_Amber_CAGo_MB(crdAA,numatom,thetaAA,								       */
/*   //			     numZ_dih,pairs_dihed_AA,									       */
/*   //			     numZ_ang,pairs_angle_AA,									       */
/*   //			     numZ_bon,pairs_bond_AA,pi);								       */
/*   //  printf("here is line 620 on TACCM_CGAAMDrun_Amber_CAGo_MB.c\n");						       */
/*   //  TACCM_CTheta_Amber_CAGo_MB(crdCG,numCAatom,thetaCG,								       */
/*   //  		       numZ_dih,pairs_dihed_CG,									       */
/*   //  		       numZ_ang,pairs_angle_CG,									       */
/*   //  		       numZ_bon,pairs_bond_CG,pi);								       */
/*   //  printf("here is line 623 on TACCM_CGAAMDrun_Amber_CAGo_MB.c\n");						       */
/*   *EAA=0.0;														       */
/*   for (i=0;i<numZ_bon;++i) {												       */
/*     delta=Z[i]-thetaAA[i];												       */
/*     *EAA+=0.5*KZAA*delta*delta;											       */
/*   }															       */
/*   for (i=0;i<numZ_ang;++i) {												       */
/*     if ((delta=Z[i+numZ_bon]-thetaAA[i+numZ_bon])>pi) delta-=2.0*pi;							       */
/*     else if ((delta=Z[i+numZ_bon]-thetaAA[i+numZ_bon])<-1.0*pi) delta+=2.0*pi;					       */
/*     *EAA+=0.5*KZAA*delta*delta;											       */
/*   }															       */
/*   for (i=0;i<numZ_dih;++i) {												       */
/*     if ((delta=Z[i+numZ_bon+numZ_ang]-thetaAA[i+numZ_bon+numZ_ang])>pi) delta-=2.0*pi;				       */
/*     else if ((delta=Z[i+numZ_bon+numZ_ang]-thetaAA[i+numZ_bon+numZ_ang])<-1.0*pi) delta+=2.0*pi;			       */
/*     *EAA+=0.5*KZAA*delta*delta;											       */
/*   }															       */
/* 															       */
/*   *ECG=0.0;														       */
/*   for (i=0;i<numZ_bon;++i) {												       */
/*     delta=Z[i]-thetaCG[i];												       */
/*     *ECG+=0.5*KZCG*delta*delta;											       */
/*   }															       */
/*   for (i=0;i<numZ_ang;++i) {												       */
/*     if ((delta=Z[i+numZ_bon]-thetaCG[i+numZ_bon])>pi) delta-=2.0*pi;							       */
/*     else if ((delta=Z[i+numZ_bon]-thetaCG[i+numZ_bon])<-1.0*pi) delta+=2.0*pi;					       */
/*     *ECG+=0.5*KZCG*delta*delta;											       */
/*   }															       */
/*   for (i=0;i<numZ_dih;++i) {												       */
/*     if ((delta=Z[i+numZ_bon+numZ_ang]-thetaCG[i+numZ_bon+numZ_ang])>pi) delta-=2.0*pi;				       */
/*     else if ((delta=Z[i+numZ_bon+numZ_ang]-thetaCG[i+numZ_bon+numZ_ang])<-1.0*pi) delta+=2.0*pi;			       */
/*     *ECG+=0.5*KZCG*delta*delta;											       */
/*   }															       */
/* 															       */
/*   *EZ=(*EAA)+(*ECG);													       */
/*   //  printf("here is line 657 on TACCM_CGAAMDrun_Amber_CAGo_MB.c\n");						       */
/* }															       */
/* 															       */
/* /\*															       */
/*   ffL_calcffandforce_AA(crdAA,numatom,&e,&f);									       */
/*   ffL_calcffandforce_AA(crdCG,numatom,&e_CG2,&f_CG);									       */
/* 															       */
/*   fAA=(double *)gcemalloc(sizeof(double)*numZ_dih);									       */
/*   fCG=(double *)gcemalloc(sizeof(double)*numZ_dih);									       */
/*   frcZ=(double *)gcemalloc(sizeof(double)*numZ_dih);									       */
/* 															       */
/*   thetaAA=(double *)gcemalloc(sizeof(double)*numZ_dih);								       */
/*   thetaCG=(double *)gcemalloc(sizeof(double)*numZ_dih);								       */
/* 															       */
/*   TACCM_CTheta(crdAA,numatom,thetaAA,numZ_dih,pairs_dihe_AA,pi);							       */
/*   TACCM_CTheta(crdCG,numCAatom,thetaCG,numZ_dih,pairs_dihe_CG,pi);							       */
/*   TACCM_calc_eff_FF_Z(Z,numZ_dih,thetaAA,KZAA,fAA,pi);								       */
/*   TACCM_calc_eff_FF_Z(Z,numZ_dih,thetaCG,KZCG,fCG,pi);								       */
/*   for (i=0;i<numZ_dih;++i) frcZ[i]=fAA[i]+fCG[i];									       */
/* 															       */
/*   for (i=0;i<numstep;++i) {												       */
/*     printf("here is line 505 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");						       */
/*     TACCM_CTheta(crdAA,numatom,thetaAA,numZ_dih,pairs_dihe_AA,pi);							       */
/*     TACCM_CTheta(crdCG,numatom,thetaCG,numZ_dih,pairs_dihe_CG,pi);							       */
/* 															       */
/*     PEAA=TACCM_MD_Propagetor_NH_MP1998_AA_test(crdAA,velAA,mass,zetaAA,V_zetaAA,					       */
/* 					       QAA,NfKTAA,numatom,&KEAA,&KEvAA,&PEvAA,					       */
/* 					       dt,dt2,nc,wdt4,wdt2,							       */
/* 					       &e,&f,Z,numZ_dih,thetaAA,KZAA,pairs_dihe_AA,PEZAA,pi);			       */
/* 															       */
/*     PECG=TACCM_MD_Propagetor_NH_MP1998_AA_test(crdCG,velCG,mass,zetaCG,V_zetaCG,					       */
/* 					       QCG,NfKTCG,numatom,&KECG,&KEvCG,&PEvCG,					       */
/* 					       dt,dt2,nc,wdt4,wdt2,							       */
/* 					       &e_CG2,&f_CG,Z,numZ_dih,thetaCG,KZCG,pairs_dihe_CG,PEZCG,pi);		       */
/* 															       */
/*     TACCM_CGAA_MD_Propagetor_NH_MP1998_Z_2(Z,velZ,massZ,thetaAA,thetaCG,zetaZ,V_zetaZ,				       */
/* 					   QZ,NfKTZ,numZ_dih,&KEZ,&KEvZ,&PEvZ,						       */
/* 					   dt,dt2,nc,wdt4,wdt2,KZAA,KZCG,						       */
/* 					   PEZAA,PEZCG,PEZ,frcZ,pi);							       */
/*       														       */
/*     if (i%interval==0) {												       */
/*       KEAA=KEAA/UNITT;     TAA=KEAA/((3*numatom)*k_B)*2.0;								       */
/*       PEvAA=PEvAA/UNITT;   KEvAA=KEvAA/UNITT;									       */
/* 															       */
/*       KECG=KECG/UNITT;     TCG=KECG/((3*numatom)*k_B)*2.0;								       */
/*       PEvCG=PEvCG/UNITT;  KEvCG=KEvCG/UNITT;										       */
/* 															       */
/*       PEAA=e.p_d_t+e.p_a_t+e.p_b_t;											       */
/*       EtAA=PEAA+KEAA+PEvAA+KEvAA;											       */
/*       fprintf(outputfileAA,"%d %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e "							       */
/* 	                    ,i+1,PEAA,KEAA,KEvAA,PEvAA,EtAA,TAA);							       */
/* 															       */
/*       PECG=e_CG2.p_d_t+e_CG2.p_a_t+e_CG2.p_b_t;									       */
/*       EtCG=PECG+KECG+PEvCG+KEvCG;											       */
/*       fprintf(outputfileCG,"%d %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e "							       */
/* 	                  ,i+1,PECG,KECG,KEvCG,PEvCG,EtCG,TCG);								       */
/* 															       */
/*       //      *avePEAA=(i*(*avePEAA)+PEAA)/(i+1); *varPEAA=(i*(*varPEAA)+PEAA*PEAA)/(i+1);				       */
/*       //      *aveKEAA=(i*(*aveKEAA)+KEAA)/(i+1); *varKEAA=(i*(*varKEAA)+KEAA*KEAA)/(i+1);				       */
/*       //      *aveTAA=(i*(*aveTAA)+TAA)/(i+1);  *varTAA=(i*(*varTAA)+TAA*TAA)/(i+1);					       */
/* 															       */
/*       //      *avePECG=(i*(*avePECG)+PECG)/(i+1); *varPECG=(i*(*varPECG)+PECG*PECG)/(i+1);				       */
/*       //      *aveKECG=(i*(*aveKECG)+KECG)/(i+1); *varKECG=(i*(*varKECG)+KECG*KECG)/(i+1);				       */
/*       //      *aveTCG=(i*(*aveTCG)+TCG)/(i+1);  *varTCG=(i*(*varTCG)+TCG*TCG)/(i+1);					       */
/* 															       */
/*       summass=0.0; for (j=0;j<numatom;++j) summass+=mass[j];								       */
/*       for (j=0;j<3;++j) COM[j]=0.0;											       */
/*       for (j=0;j<numatom;++j)  for (k=0;k<3;++k) COM[k]+=mass[j]*crdAA[j*3+k]/summass;				       */
/*       for (j=0;j<numatom;++j)  for (k=0;k<3;++k) crdAA[j*3+k]-=COM[k];						       */
/*       for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crdAA[j*3+k];						       */
/*       myncL_put_crd_ene_MCD(nc_id_MCDAA,*l,crd_nc,e,0.0);								       */
/* 															       */
/*       for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crdCG[j*3+k];						       */
/*       myncL_put_crd_ene_MCD(nc_id_MCDCG,*l,crd_nc,e,0.0);								       */
/*       ++(*l);													       */
/* 															       */
/*       ///////////////// TACCM //////////////////////									       */
/*       KEZ=KEZ/UNITT;      TZ=KEZ/(numZ_dih*k_B)*2.0;									       */
/* 															       */
/*       PEvZ=PEvZ/UNITT;      KEvZ=KEvZ/UNITT;										       */
/* 															       */
/*       EtZ=*PEZ+KEZ+PEvZ+KEvZ;											       */
/*       fprintf(outputfileAA,"%d %e %e %e %e %e %e %e\n",i+1,*PEZ,KEZ,KEvZ,PEvZ,EtZ,TZ);				       */
/* 															       */
/*       for (j=0;j<numZ_dih;++j) fprintf(trjfileZ,"%e ",Z[j]);								       */
/*       fprintf(trjfileZ,"\n");											       */
/*       for (j=0;j<numZ_dih;++j) fprintf(trjfilThetaAA,"%e ",thetaAA[j]);						       */
/*       fprintf(trjfilThetaAA,"\n");      										       */
/*       for (j=0;j<numZ_dih;++j) fprintf(trjfilThetaCG,"%e ",thetaCG[j]);						       */
/*       fprintf(trjfilThetaCG,"\n");      										       */
/*       ///////////////// TACCM //////////////////////									       */
/*     }														       */
/*   }															       */
/* *\/															       */
/*******************************************************************************************************************************/
/*  */

