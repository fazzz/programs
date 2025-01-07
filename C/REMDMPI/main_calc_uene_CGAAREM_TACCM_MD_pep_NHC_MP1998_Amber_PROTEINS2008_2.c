#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "EF.h"

#include "REMDCGAA_TACCM_calc_uene_Amber_PROTEINS2008.h"

#include "REMDCGAA_TACCM_MPI_2_Amber_PROTEINS2008.h"

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

  struct AADataforREMD_Amber AAdata;
  struct CGDataforREMD_PROTEINS2008 CGdata;
  struct AACGCommonDataforREMD_A_P2008 Cdata;
  struct TACCMDataforREMD_A_P2008 Zdata;
  double *KZAAs,*KZCGs;
  double KZAAobj=0.0,KZCGobj=0.0;

  struct potential e;
  struct force f;

  double EAA,ECG,EZ;
  double EAAdiff,ECGdiff,EZdiff;

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

  char **trjfileAAname,**trjfileCGname,**trjfilenameZ,*refcrdfilename;
  char *inputfilename,*parmfilename,*TACCMfilename,*parametertrjfilename;
  FILE *inputfile,*parmfile,*TACCMfile,*parametertrjfile,*refcrdfile;

  char *outputfilenamebase,outputfilename[2000];

  FILE **trjfileZ,**outputfile;

  double *refcrd,ep=0.3,criteria=6.5;
  int NCmode=3,nibnum=3,numnb,num14;

  double crd_nc_AA[MAXATOM][3],crd_nc_CG[MAXATOM][3];
  struct my_netcdf_out_id_MCD *nc_id_MD_AA,*nc_id_MD_CG;

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
    {"AA",0,NULL,'o'}, {"CG",0,NULL,'P'}, {"KZAAo",1,NULL,'A'}, {"KZCGo",1,NULL,'C'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hnaoPt:s:p:c:b:e:N:I:E:S:i:A:C:",long_opt,&opt_idx))!=-1) {
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
    case 'P':
      objmodeflag=CG;  break;
    case 'p':
      ep=atof(optarg); break;
    case 'c':
      criteria=atof(optarg); break;
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

  if (argc < 5) {
    USAGE(progname);
    exit(1);
  }
  inputfilename        = *argv;
  parametertrjfilename = *++argv;
  parmfilename         = *++argv;
  refcrdfilename       = *++argv;
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

  refcrd=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
  refcrdfile=efopen(refcrdfilename,"r");
  //  getline(&line,&len,refcrdfile);
  fscanf(refcrdfile,"%d",&d);
  for (i=0;i<Cdata.numatom;++i) for (j=0;j<3;++j) fscanf(refcrdfile,"%lf",&refcrd[i*3+j]);
  fclose(refcrdfile);

  ffL_set_calcffandforce(&(AAdata.e),&(AAdata.f));

  ffL_set_calcffandforce(&e,&f);
  ffL_set_non_bonding_index_1(&numnb,&num14);
  e.parm.numnb=numnb;
  e.parm.num14=num14;
  e.parm.indexnb=(int *)gcemalloc(sizeof(int)*numnb*2);
  e.parm.index14=(int *)gcemalloc(sizeof(int)*num14*2);
  ffL_set_non_bonding_index_2(e.parm.indexnb,e.parm.index14);

  GOLMAA_PROTEINS2008_ff_set_calcff_b(&(CGdata.e),refcrd,Cdata.numatom,Cdata.numres,
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
	    ffL_calcffandforce(crd_AA,Cdata.numatom,&(AAdata.e),&(AAdata.f),AAdata.e);
	  if (objmodeflag==FG)
	    GOLMAA_PROTEINS2008_ff_calcff_wobaimp_b(crd_CG,Cdata.numatom,&(CGdata.e));

	  CE_TACCM_CGAA(crd_AA,crd_CG,Z,Cdata.numatom,Zdata.numZ,KZAAobj,KZCGobj,Zdata.pairs,pi,
			&EAAdiff,&ECGdiff,&EZdiff);
	  if (objmodeflag==CG) {
	    EZ=-1.0*((ECGdiff-ECG)+EAA+AAdata.e.p_e_t+AAdata.e.p_LJ_t+AAdata.e.p_e_14_t+AAdata.e.p_LJ_14_t+AAdata.e.p_d_t+AAdata.e.p_a_t+AAdata.e.p_b_t);
	  }
	  if (objmodeflag==FG) {
	    EZ=(EAAdiff-EAA)+ECG+CGdata.e.p_t;
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

