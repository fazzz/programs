
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "EF.h"

#include "REMDCGAA_TACCM_MPI_2_testb.h"
#include "REMDCGAA_TACCM_MPI_2_testc.h"

//#include "REMDCGAA_TACCM_calc_uene_Amber_PROTEINS2008.h"

//#include "REMDCGAA_TACCM_MPI_2_Amber_PROTEINS2008.h"

#define FG 1
#define CG 0

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,m,n,d;
  int index;
  int *numbering;

  int numstep=1000,interval=100;
  int numEX=1,numRE=2;

  int equflag=OFF,objmodeflag=FG;
  int numstepequ=0,intervalequ=1;

  struct AADataforREMD_testb AAdata, CGdata;
  struct AACGCommonDataforREMD_testb Cdata;
  struct TACCMDataforREMD_testb Zdata;
  struct AmberParmL ap_FG,ap_CG;

  double *KZAAs,*KZCGs;
  double KZAAobj=0.0,KZCGobj=0.0;

  double EAA,ECG,EZ;
  double EAAdiff,ECGdiff,EZdiff,Ediff;

  int addflag=OFF,nobetaflag=OFF,AMBERMODEflag=OFF;

  double T=300,T_sim=300;
  double k_B=1.98723e-3;
  double beta=1.0;

  double p_t=0.0;

  int *indexTACCM,**pairsZ;

  int **sereis;

  double *crd_AA,*crd_CG,*Z;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double pi;

  char **trjfileAAname,**trjfileCGname,**trjfilenameZ;
  char *inputfilename,*parmtopfilename_FG,*parmtopfilename_CG,*TACCMfilename,*parametertrjfilename;
  FILE *inputfile,*parmtopfile_FG,*parmtopfile_CG,*TACCMfile,*parametertrjfile;

  char *outputfilenamebase,outputfilename[2000];

  FILE **trjfileZ,**outputfile;

  double crd_nc_AA[MAXATOM][3],crd_nc_CG[MAXATOM][3];
  struct my_netcdf_out_id_MCD *nc_id_MD_AA,*nc_id_MD_CG;

  char *progname;

  int opt_idx=1;

  pi=acos(-1.0);

  struct option long_opt[] = {
    {"tempobj",1,NULL,'t'}, {"tempsim",1,NULL,'s'},
    {"numRE",1,NULL,'N'}, {"numEX",1,NULL,'e'},
    {"nums",1,NULL,'S'}, {"int",1,NULL,'i'},
    {"equ",1,NULL,'E'}, {"intequ",1,NULL,'I'},
    {"nobeta",0,NULL,'n'},  {"a",0,NULL,'a'},
    {"AA",0,NULL,'o'}, {"CG",0,NULL,'p'}, {"KZAAo",1,NULL,'A'}, {"KZCGo",1,NULL,'C'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hnat:s:b:e:N:I:E:S:i:o:A:C:",long_opt,&opt_idx))!=-1) {
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
    case 'o':
      objmodeflag=FG;  break;
    case 'p':
      objmodeflag=CG;  break;
    case 'A':
      KZAAobj=atof(optarg);  break;
    case 'C':
      KZCGobj=atof(optarg);  break;
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
  parmtopfilename_FG   = *++argv;
  parmtopfilename_CG   = *++argv;
  TACCMfilename        = *++argv;
  outputfilenamebase   = *++argv;

  parmtopfile_FG=efopen(parmtopfilename_FG,"r");
  readParmtopLb(parmtopfile_FG,&ap_FG);
  fclose(parmtopfile_FG); 

  parmtopfile_CG=efopen(parmtopfilename_CG,"r");
  readParmtopLb(parmtopfile_CG,&ap_CG);
  fclose(parmtopfile_CG); 

  Cdata.numatom=ap_FG.NATOM;
  j=0;  for (i=0;i<Cdata.numatom;++i)  if (strncmp(ap_FG.IGRAPH[i],"H",1)==0)  ++j;
  Cdata.numheavyatom=Cdata.numatom-j;  Cdata.numres=ap_FG.NRES;

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
  crd_CG=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
  Z=(double *)gcemalloc(sizeof(double)*Zdata.numZ);

  trjfileAAname=(char **)gcemalloc(sizeof(char *)*numRE);
  trjfileCGname=(char **)gcemalloc(sizeof(char *)*numRE);
  trjfilenameZ=(char **)gcemalloc(sizeof(char *)*numRE);
  for (i=0;i<numRE;++i) {
    trjfileAAname[i]=(char *)gcemalloc(sizeof(char)*1000);
    trjfileCGname[i]=(char *)gcemalloc(sizeof(char)*1000);
    trjfilenameZ[i]=(char *)gcemalloc(sizeof(char)*1000);
  }
  nc_id_MD_AA=(struct my_netcdf_out_id_MCD *)gcemalloc(sizeof(struct my_netcdf_out_id_MCD)*numRE);
  nc_id_MD_CG=(struct my_netcdf_out_id_MCD *)gcemalloc(sizeof(struct my_netcdf_out_id_MCD)*numRE);
  trjfileZ=(FILE **)gcemalloc(sizeof(trjfileZ)*numRE);
  KZAAs=(double *)gcemalloc(sizeof(double)*numRE);
  KZCGs=(double *)gcemalloc(sizeof(double)*numRE);

  inputfile=efopen(inputfilename,"r");
  CGAAREMDreadInputs_calc_uene(inputfile,Cdata.numatom,numRE,KZAAs,KZCGs,
			       trjfileAAname,trjfileCGname,
			       nc_id_MD_AA,nc_id_MD_CG,trjfileZ);
  fclose(inputfile);

  sereis=(int **)gcemalloc(sizeof(int *)*numRE);
  for (i=0;i<numRE;++i) sereis[i]=(int *)gcemalloc(sizeof(int)*numEX);
  numbering=(int *)gcemalloc(sizeof(int)*numRE*numRE);
  for (i=0;i<numRE*numRE;++i) numbering[i]=0;

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

  ffLc_set_calcffandforce(&(AAdata.e),&(AAdata.f),ap_FG);
  ffLc_set_calcffandforce((&CGdata.e),&(CGdata.f),ap_CG);

  numstep=numstep/interval;
  numstepequ=numstepequ/intervalequ;
  for (i=0;i<numRE;++i) {
    l=0;
    if (equflag==ON) {
      for (j=0;j<numstepequ;++j) {
	mync_open_inq_get_sh_MCD(trjfileAAname[i],Cdata.numatom,l,1,l+1,&(nc_id_MD_AA[i]),crd_nc_AA);
	mync_open_inq_get_sh_MCD(trjfileCGname[i],Cdata.numatom,l,1,l+1,&(nc_id_MD_CG[i]),crd_nc_CG);
	for (k=0;k<Zdata.numZ;++k) fscanf(trjfileZ[i],"%lf",&Z[k]);
	++l;
      }
    }

    for (j=0;j<numEX;++j) {
      index=sereis[i][j];

      for (k=0;k<numstep;++k) {
	mync_open_inq_get_sh_MCD(trjfileAAname[i],Cdata.numatom,l,1,l+1,&(nc_id_MD_AA[i]),crd_nc_AA);
	mync_open_inq_get_sh_MCD(trjfileCGname[i],Cdata.numatom,l,1,l+1,&(nc_id_MD_CG[i]),crd_nc_CG);
	++l;

	for (m=0;m<Zdata.numZ;++m) fscanf(trjfileZ[i],"%lf",&Z[m]);
	for (m=0;m<Cdata.numatom;++m) {
	  for (n=0;n<3;++n) {
	    crd_AA[m*3+n]=crd_nc_AA[m][n];  crd_CG[m*3+n]=crd_nc_CG[m][n];
	  }
	}

	for (m=0;m<numRE;++m) {
	  CE_TACCM_CGAA(crd_AA,crd_CG,Z,Cdata.numatom,Zdata.numZ,KZAAs[m],KZCGs[m],Zdata.pairs,pi,&EAA,&ECG,&EZ);
	  if (objmodeflag==CG)
	    ffLc_calcffandforce(crd_AA,Cdata.numatom,&(AAdata.e),&(AAdata.f),ap_FG);
	  if (objmodeflag==FG)
	    ffLc_calcffandforce(crd_CG,Cdata.numatom,&(CGdata.e),&(CGdata.f),ap_CG);

	  CE_TACCM_CGAA(crd_AA,crd_CG,Z,Cdata.numatom,Zdata.numZ,KZAAobj,KZCGobj,Zdata.pairs,pi,
			&EAAdiff,&ECGdiff,&EZdiff);
	  if (objmodeflag==CG) {
	    EZ=(ECGdiff-ECG)-1.0*(EAA+/*+AAdata.e.p_t*/AAdata.e.p_e_t+AAdata.e.p_LJ_t+AAdata.e.p_e_14_t+AAdata.e.p_LJ_14_t+AAdata.e.p_d_t+AAdata.e.p_a_t+AAdata.e.p_b_t);
	    //	    EZ=-1.0*(ECGdiff-ECG)+(EAA+AAdata.e.p_d_t+AAdata.e.p_a_t+AAdata.e.p_b_t);
	    //	    EZ=-1.0*((ECGdiff-ECG)+EAA+AAdata.e.p_e_t+AAdata.e.p_LJ_t+AAdata.e.p_e_14_t+AAdata.e.p_LJ_14_t+AAdata.e.p_d_t+AAdata.e.p_a_t+AAdata.e.p_b_t);
	  }
	  if (objmodeflag==FG) {
	    //	    EZ=(EAAdiff-EAA)+ECG+CGdata.e.p_e_t+CGdata.e.p_LJ_t+CGdata.e.p_e_14_t+CGdata.e.p_LJ_14_t+CGdata.e.p_d_t+CGdata.e.p_a_t+CGdata.e.p_b_t;
	    EZ=(EAAdiff-EAA)-1.0*(ECG+CGdata.e.p_e_t+CGdata.e.p_LJ_t+CGdata.e.p_e_14_t+CGdata.e.p_LJ_14_t+CGdata.e.p_d_t+CGdata.e.p_a_t+CGdata.e.p_b_t);
	    //	    EZ=-1.0*(EAAdiff-EAA)+(ECG+CGdata.e.p_d_t+CGdata.e.p_a_t+CGdata.e.p_b_t);
	  }

	  numbering[index*numRE+m]+=1;
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
