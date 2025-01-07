#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "EF.h"

#include "REMDCGAA_TACCM_calc_uene_Amber_PROTEINS2008.h"
#include "REMDCGAA_TACCM_MPI_2_Amber_2PROTEINS2008.h"

#include "REMDCGAA_TACCM_MPI_2_Amber_PROTEINS2008.h"

#define FG  0
#define CG1 1
#define CG2 2

#define ON 1
#define OFF 0

double CE_TACCM_2CG1AA_oc(double *crdAA,double *crdCG1,double *crdCG2,double *Z, int numatom,int numZ,
			  double KZAA,double KZCG1,double KZCG2,int **pairs,double pi,
			  double *EAA,double *ECG1,double *ECG2,double *EZ);

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,m,n,d;
  int index;
  int *numbering;

  int numstep=1000,interval=100;
  int numEX=1,numRE=2;

  int equflag=OFF,objmodeflag=FG;
  int numstepequ=0,intervalequ=1;

  struct AADataforREMD_Amber AAdata;
  struct CGDataforREMD_PROTEINS2008 CGdata[2];
  struct AACGCommonDataforREMD_A_P2008 Cdata;
  struct TACCMDataforREMD_A_2P2008 Zdata;
  double *KZAAs,*KZCG1s,*KZCG2s;
  double KZAAobj=0.0,KZCG1obj=0.0,KZCG2obj=0.0;

  struct potential e;
  struct force f;

  double a,b;

  double EAA,ECG1,ECG2,EZ;
  double EAAdiff,ECG1diff,ECG2diff,EZdiff;

  int addflag=OFF,nobetaflag=OFF,AMBERMODEflag=OFF;

  double T=300,T_sim=300;
  double k_B=1.98723e-3;
  double beta=1.0;

  double p_t=0.0;

  int *indexTACCM,**pairsZ;

  int **sereis;

  double *crd_AA,*crd_CG1,*crd_CG2,*Z;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double pi;

  char **trjfileAAname,**trjfileCG1name,**trjfileCG2name,**trjfilenameZ,*refcrdfilename1,*refcrdfilename2;
  char *inputfilename,*parmfilename,*TACCMfilename,*parametertrjfilename;
  FILE *inputfile,*parmfile,*TACCMfile,*parametertrjfile,*refcrdfile1,*refcrdfile2;

  char *outputfilenamebase,outputfilename[2000];

  FILE **trjfileZ,**outputfile;

  double *refcrd1,*refcrd2,ep=0.3,criteria=6.5;
  int NCmode=3,nibnum=3,numnb,num14;

  double crd_nc_AA[MAXATOM][3],crd_nc_CG1[MAXATOM][3],crd_nc_CG2[MAXATOM][3];
  struct my_netcdf_out_id_MCD *nc_id_MD_AA,*nc_id_MD_CG1,*nc_id_MD_CG2;

  char *progname;

  int opt_idx=1;

  pi=acos(-1.0);

  struct option long_opt[] = {
    {"tempobj",1,NULL,'t'}, {"tempsim",1,NULL,'s'},
    {"numRE",1,NULL,'N'}, {"numEX",1,NULL,'e'},
    {"ep",1,NULL,'p'}, {"cutoff",1,NULL,'c'},
    {"nums",1,NULL,'S'}, {"int",1,NULL,'i'},
    {"equ",1,NULL,'E'}, {"intequ",1,NULL,'I'},
    {"nobeta",0,NULL,'n'},  {"a",0,NULL,'a'},
    {"AA",0,NULL,'0'}, {"CG1",0,NULL,'1'}, {"CG2",0,NULL,'2'},
    {"KZAAo",1,NULL,'A'}, {"KZCG1o",1,NULL,'C'}, {"KZCG2o",1,NULL,'3'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"h012naoPt:s:p:c:b:e:N:I:E:S:i:A:C:3:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 't':
      T=atof(optarg); break;
    case 's':
      T_sim=atof(optarg); break;
    case 'a':
      addflag=ON; break;
    case 'n':
      nobetaflag=ON; break;
    case 'e':
      numEX=atoi(optarg); break;
    case 'N':
      numRE=atoi(optarg); break;
    case 'S':
      numstep=atoi(optarg); break;
    case 'i':
      interval=atoi(optarg); break;
    case 'E':
      equflag=ON;
      numstepequ=atoi(optarg);  break;
    case 'I':
      intervalequ=atoi(optarg);  break;
    case '0':
      objmodeflag=FG;  break;
    case '1':
      objmodeflag=CG1;  break;
    case '2':
      objmodeflag=CG2;  break;
    case 'p':
      ep=atof(optarg); break;
    case 'c':
      criteria=atof(optarg); break;
    case 'A':
      KZAAobj=atof(optarg);  break;
    case 'C':
      KZCG1obj=atof(optarg);  break;
    case '3':
      KZCG2obj=atof(optarg);  break;
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (argc < 6) {
    USAGE(progname);
    exit(1);
  }
  inputfilename        = *argv;
  parametertrjfilename = *++argv;
  parmfilename         = *++argv;
  refcrdfilename1      = *++argv;
  refcrdfilename2      = *++argv;
  TACCMfilename        = *++argv;
  outputfilenamebase   = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 

  Cdata.numatom=AP.NATOM;
  j=0;  for (i=0;i<Cdata.numatom;++i)  if (strncmp(AP.IGRAPH[i],"H",1)==0)  ++j;
  Cdata.numheavyatom=Cdata.numatom-j;  Cdata.numres=AP.NRES;

  TACCMfile=efopen(TACCMfilename,"r");
  fscanf(TACCMfile,"%d",&Zdata.numZ);
  Zdata.pairs=(int **)gcemalloc(sizeof(int *)*Zdata.numZ);
  for (i=0;i<Zdata.numZ;++i) Zdata.pairs[i]=(int *)gcemalloc(sizeof(int)*5);
  for (i=0;i<Zdata.numZ;++i) {
    for (j=0;j<4;++j) fscanf(TACCMfile,"%d",&(Zdata.pairs[i][j]));
    fscanf(TACCMfile,"%d",&(Zdata.pairs[i][j]));
  }
  fclose(TACCMfile);

  crd_AA=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
  crd_CG1=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
  crd_CG2=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
  Z=(double *)gcemalloc(sizeof(double)*Zdata.numZ);

  trjfileAAname=(char **)gcemalloc(sizeof(char *)*numRE);
  trjfileCG1name=(char **)gcemalloc(sizeof(char *)*numRE);
  trjfileCG2name=(char **)gcemalloc(sizeof(char *)*numRE);
  trjfilenameZ=(char **)gcemalloc(sizeof(char *)*numRE);
  for (i=0;i<numRE;++i) {
    trjfileAAname[i]=(char *)gcemalloc(sizeof(char)*1000);
    trjfileCG1name[i]=(char *)gcemalloc(sizeof(char)*1000);
    trjfileCG2name[i]=(char *)gcemalloc(sizeof(char)*1000);
    trjfilenameZ[i]=(char *)gcemalloc(sizeof(char)*1000);
  }
  nc_id_MD_AA=(struct my_netcdf_out_id_MCD *)gcemalloc(sizeof(struct my_netcdf_out_id_MCD)*numRE);
  nc_id_MD_CG1=(struct my_netcdf_out_id_MCD *)gcemalloc(sizeof(struct my_netcdf_out_id_MCD)*numRE);
  nc_id_MD_CG2=(struct my_netcdf_out_id_MCD *)gcemalloc(sizeof(struct my_netcdf_out_id_MCD)*numRE);
  trjfileZ=(FILE **)gcemalloc(sizeof(trjfileZ)*numRE);
  KZAAs=(double *)gcemalloc(sizeof(double)*numRE);
  KZCG1s=(double *)gcemalloc(sizeof(double)*numRE);
  KZCG2s=(double *)gcemalloc(sizeof(double)*numRE);

  inputfile=efopen(inputfilename,"r");
  CGAAREMDreadInputs_calc_uene_1FG2CG(inputfile,Cdata.numatom,numRE,KZAAs,KZCG1s,KZCG2s,
				      trjfileAAname,trjfileCG1name,trjfileCG2name,
				      nc_id_MD_AA,nc_id_MD_CG1,nc_id_MD_CG2,trjfileZ);
  fclose(inputfile);

  sereis=(int **)gcemalloc(sizeof(int *)*numRE);
  for (i=0;i<numRE;++i) sereis[i]=(int *)gcemalloc(sizeof(int)*numEX);
  numbering=(int *)gcemalloc(sizeof(int)*numRE*numRE);
  for (i=0;i<numRE*numRE;++i) numbering[i]=0;

  refcrd1=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
  refcrdfile1=efopen(refcrdfilename1,"r");
  //  getline(&line,&len,refcrdfile);
  fscanf(refcrdfile1,"%d",&d);
  for (i=0;i<Cdata.numatom;++i) for (j=0;j<3;++j) fscanf(refcrdfile1,"%lf",&refcrd1[i*3+j]);
  fclose(refcrdfile1);

  refcrd2=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
  refcrdfile2=efopen(refcrdfilename2,"r");
  //  getline(&line,&len,refcrdfile);
  fscanf(refcrdfile2,"%d",&d);
  for (i=0;i<Cdata.numatom;++i) for (j=0;j<3;++j) fscanf(refcrdfile2,"%lf",&refcrd2[i*3+j]);
  fclose(refcrdfile2);

  ffL_set_calcffandforce(&(AAdata.e),&(AAdata.f));

  ffL_set_calcffandforce(&e,&f);
  ffL_set_non_bonding_index_1(&numnb,&num14);
  e.parm.numnb=numnb;
  e.parm.num14=num14;
  e.parm.indexnb=(int *)gcemalloc(sizeof(int)*numnb*2);
  e.parm.index14=(int *)gcemalloc(sizeof(int)*num14*2);
  ffL_set_non_bonding_index_2(e.parm.indexnb,e.parm.index14);

  GOLMAA_PROTEINS2008_ff_set_calcff_b(&(CGdata[0].e),refcrd1,Cdata.numatom,Cdata.numres,
  				      /*AAData.*/e.parm.indexnb,/*AAData.*/e.parm.numnb,ep,nibnum,criteria);

  GOLMAA_PROTEINS2008_ff_set_calcff_b(&(CGdata[1].e),refcrd2,Cdata.numatom,Cdata.numres,
  				      /*AAData.*/e.parm.indexnb,/*AAData.*/e.parm.numnb,ep,nibnum,criteria);

  parametertrjfile=efopen(parametertrjfilename,"r");
  for (i=0;i<numRE;++i)
    for (j=0;j<numEX;++j) 
      fscanf(parametertrjfile,"%d",&(sereis[i][j]));
  fclose(parametertrjfile);

  if ( nobetaflag==OFF ) beta=1.0/k_B*(1.0/T_sim-1.0/T);

  outputfile=(FILE **)gcemalloc(sizeof(FILE *)*numRE*numRE);
  for (i=0;i<numRE;++i) {
    for (j=0;j<numRE;++j) {
      sprintf(outputfilename,"%s_bo_%d_wp_%d\0",outputfilenamebase,i+1,j+1);
      if ( addflag==OFF ) outputfile[i*numRE+j]=efopen(outputfilename,"w");
      else if ( addflag==ON )  outputfile[i*numRE+j]=efopen(outputfilename,"a");
    }
  }

  numstep=numstep/interval;
  numstepequ=numstepequ/intervalequ;
  for (i=0;i<numRE;++i) {
    l=0;
    if (equflag==ON) {
      for (j=0;j<numstepequ;++j) {
	mync_open_inq_get_sh_MCD(trjfileAAname[i],Cdata.numatom,l,1,l+1,&(nc_id_MD_AA[i]),crd_nc_AA);
	mync_open_inq_get_sh_MCD(trjfileCG1name[i],Cdata.numatom,l,1,l+1,&(nc_id_MD_CG1[i]),crd_nc_CG1);
	mync_open_inq_get_sh_MCD(trjfileCG2name[i],Cdata.numatom,l,1,l+1,&(nc_id_MD_CG2[i]),crd_nc_CG2);
	for (k=0;k<Zdata.numZ;++k) fscanf(trjfileZ[i],"%lf",&Z[k]);
	++l;
      }
    }

    for (j=0;j<numEX;++j) {
      index=sereis[i][j];

      for (k=0;k<numstep;++k) {
	mync_open_inq_get_sh_MCD(trjfileAAname[i],Cdata.numatom,l,1,l+1,&(nc_id_MD_AA[i]),crd_nc_AA);
	mync_open_inq_get_sh_MCD(trjfileCG1name[i],Cdata.numatom,l,1,l+1,&(nc_id_MD_CG1[i]),crd_nc_CG1);
	mync_open_inq_get_sh_MCD(trjfileCG2name[i],Cdata.numatom,l,1,l+1,&(nc_id_MD_CG2[i]),crd_nc_CG2);
	++l;

	for (m=0;m<Zdata.numZ;++m) fscanf(trjfileZ[i],"%lf",&Z[m]);
	for (m=0;m<Cdata.numatom;++m) {
	  for (n=0;n<3;++n) {
	    crd_AA[m*3+n]=crd_nc_AA[m][n];  crd_CG1[m*3+n]=crd_nc_CG1[m][n];  crd_CG2[m*3+n]=crd_nc_CG2[m][n];
	  }
	}

	for (m=0;m<numRE;++m) {
	  CE_TACCM_2CG1AA_oc(crd_AA,crd_CG1,crd_CG2,Z,
			     Cdata.numatom,Zdata.numZ,
			     KZAAs[m],KZCG1s[m],KZCG2s[m],Zdata.pairs,pi,
			     &EAA,&ECG1,&ECG2,&EZ);
	  if (objmodeflag==CG1) {
	    ffL_calcffandforce(crd_AA,Cdata.numatom,&(AAdata.e),&(AAdata.f),AAdata.e);
	    GOLMAA_PROTEINS2008_ff_calcff_wobaimp_b(crd_CG2,Cdata.numatom,&(CGdata[1].e));
	  }
	  if (objmodeflag==CG2) {
	    ffL_calcffandforce(crd_AA,Cdata.numatom,&(AAdata.e),&(AAdata.f),AAdata.e);
	    GOLMAA_PROTEINS2008_ff_calcff_wobaimp_b(crd_CG1,Cdata.numatom,&(CGdata[0].e));
	  }
	  if (objmodeflag==FG) {
	  /*****************************************************************************************/
          /* GOLMAA_PROTEINS2008_ff_calcff_wobaimp_b(crd_CG2,Cdata.numatom,&(CGdata[0].e));	     */
	  /* a=CGdata[0].e.p_t;								     */
	  /* GOLMAA_PROTEINS2008_ff_calcff_wobaimp_b(crd_CG1,Cdata.numatom,&(CGdata[0].e));	     */
	  /* GOLMAA_PROTEINS2008_ff_calcff_wobaimp_b(crd_CG1,Cdata.numatom,&(CGdata[1].e));	     */
	  /* b=CGdata[1].e.p_t;								     */
	  /* GOLMAA_PROTEINS2008_ff_calcff_wobaimp_b(crd_CG2,Cdata.numatom,&(CGdata[1].e));	     */
          /*****************************************************************************************/
	    GOLMAA_PROTEINS2008_ff_calcff_wobaimp_b(crd_CG1,Cdata.numatom,&(CGdata[0].e));
	    GOLMAA_PROTEINS2008_ff_calcff_wobaimp_b(crd_CG2,Cdata.numatom,&(CGdata[1].e));
	  }

	  /*****************************************************************************************/
          /* ffL_calcffandforce(crd_AA,Cdata.numatom,&(AAdata.e),&(AAdata.f),AAdata.e);		   */
	  /* GOLMAA_PROTEINS2008_ff_calcff_wobaimp_b(crd_CG1,Cdata.numatom,&(CGdata[0].e));	   */
	  /* GOLMAA_PROTEINS2008_ff_calcff_wobaimp_b(crd_CG2,Cdata.numatom,&(CGdata[1].e));	   */
          /*****************************************************************************************/


	  CE_TACCM_2CG1AA_oc(crd_AA,crd_CG1,crd_CG2,Z,
			     Cdata.numatom,Zdata.numZ,
			     KZAAobj,KZCG1obj,KZCG2obj,Zdata.pairs,pi,
			     &EAAdiff,&ECG1diff,&ECG2diff,&EZdiff);

	  if (objmodeflag==CG1) {
	    EZ=((ECG1diff-ECG1)+EAA+AAdata.e.p_e_t+AAdata.e.p_LJ_t+AAdata.e.p_e_14_t+AAdata.e.p_LJ_14_t+AAdata.e.p_d_t+AAdata.e.p_a_t+AAdata.e.p_b_t)+ECG2+CGdata[1].e.p_t;
	  }
	  if (objmodeflag==CG2) {
	    EZ=((ECG2diff-ECG2)+EAA+AAdata.e.p_e_t+AAdata.e.p_LJ_t+AAdata.e.p_e_14_t+AAdata.e.p_LJ_14_t+AAdata.e.p_d_t+AAdata.e.p_a_t+AAdata.e.p_b_t)+ECG1+CGdata[0].e.p_t;
	  }
	  if (objmodeflag==FG) {
	    EZ=(EAAdiff-EAA)+ECG1+CGdata[0].e.p_t+ECG2+CGdata[1].e.p_t;
	    //	    EZ=EAA+ECG1+CGdata[0].e.p_t+ECG2+CGdata[1].e.p_t;
	  }

	  numbering[index*numRE+m]+=1;          
          if (EZ>/*500*/1000) // 12-06-12
	    EZ=0.0;

	  fprintf(outputfile[index*numRE+m],"%d %e \n",numbering[index*numRE+m],EZ*beta);
	}
      }
    }
  }
  for (i=0;i<numRE;++i) for (j=0;j<numRE;++j) fclose(outputfile[i*numRE+j]);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename parametertrjfilename parmfilename TACCMfilename outputfilename\n",progname);
}

double CE_TACCM_2CG1AA_oc(double *crdAA,double *crdCG1,double *crdCG2,double *Z, int numatom,int numZ,
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
