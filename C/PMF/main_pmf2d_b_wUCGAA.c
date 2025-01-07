
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

//#include "PMF.h"
#include "PTLb.h"
#include "EF.h"
#include "FFLc.h"


//#include "REMDCGAA_TACCM_calc_uene_Amber_PROTEINS2008.h"

//#include "REMDCGAA_TACCM_MPI_2_Amber_PROTEINS2008.h"

#include "netcdf_mineL.h"

#include "TACCM_CGAAMDrun.h"

#define ON 1
#define OFF 0

double *io_scandcoldata2(FILE *inputfile,int numi,int numcol,int xcol,int ycol,int *numstep,double *data);

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;

  int minmaxflag=OFF;

  int numstep;
  double *data;
  int numi=0,numcol=2,xcol=1,ycol=2;
  double maxx=1.0,maxy=1.0,minx=0.0,miny=0.0;
  int framex,framey;
  double width,pi;
  double *pmf_AA_CG,*pmf_AA,*pmf_CG,pmf_min_AA_CG=0.0,pmf_min_AA=0.0,pmf_min_CG=0.0,sum_AA_CG,sum_AA,sum_CG;

  struct AmberParmL ap_FG,ap_CG;

  struct potential e_FG,e_CG;
  struct force f_FG,f_CG;
  
  int **pairs;
  double *crd_AA,*crd_CG,*Z;
  double KZAA,KZCG;
  double EAA,ECG,EZ;
  double expBEAA,expBECG;

  int numatom,numZ;

  double T=300.0,beta=1.0;
  double k_B=1.98723e-3;
  double UNITT=418.4070;
  
  char *TACCMfilename,*parmtopfilename_FG,*parmtopfilename_CG;
  char *trjfilename_AA,*trjfilename_CG,*trjfilenameZ;
  char *inputfilename,*outputfilename_AA_CG,*outputfilename_AA,*outputfilename_CG;

  FILE *TACCMfile,*parmtopfile_FG,*parmtopfile_CG;
  FILE *trjfileAA,*trjfileCG,*trjfileZ;
  FILE *inputfile1,*outputfile_AA_CG,*outputfile_AA,*outputfile_CG;
  FILE *logfile;

  char *line;
  size_t len=0;

  int d;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double crd_nc_AA[MAXATOM][3],crd_nc_CG[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MD_AA,nc_id_MD_CG;

  char *progname,*logfilename;

  int opt_idx=1;

  struct option long_opt[] = {
    {"width",1,NULL,'w'},
    {"temp",1,NULL,'t'},
    {"beta",1,NULL,'b'},
    {"numi",1,NULL,'i'},
    {"numcol",1,NULL,'c'},
    {"xcol",1,NULL,'x'},
    {"ycol",1,NULL,'y'},
    {"maxx",1,NULL,'A'},
    {"minx",1,NULL,'I'},
    {"maxy",1,NULL,'@'},
    {"miny",1,NULL,'j'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((d=getopt_long(argc,argv,"hw:i:t:b:c:x:y:A:I:@:j:",long_opt,&opt_idx))!=-1) {
    switch(d) {
    case 'w':
      width=atof(optarg);
      break;
    case 't':
      T=atof(optarg);
      break;
    case 'b':
      beta=atof(optarg);
      break;
    case 'i':
      numi=atoi(optarg);
      break;
    case 'x':
      xcol=atoi(optarg);
      break;
    case 'y':
      ycol=atoi(optarg);
      break;
    case 'A':
      minmaxflag=ON;
      maxx=atof(optarg);
      break;
    case 'I':
      minmaxflag=ON;
      minx=atof(optarg);
      break;
    case '@':
      minmaxflag=ON;
      maxy=atof(optarg);
      break;
    case 'j':
      minmaxflag=ON;
      miny=atof(optarg);
      break;
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

  pi=acos(-1.0);

  if (argc < 12) {
    USAGE(progname);
    exit(1);
  }  
  KZAA = atof(*argv);
  KZCG = atof(*++argv);
  inputfilename     = *++argv;
  trjfilename_AA    = *++argv;
  trjfilename_CG    = *++argv;
  trjfilenameZ      = *++argv;
  parmtopfilename_FG    = *++argv;
  parmtopfilename_CG    = *++argv;
  TACCMfilename     = *++argv;
  outputfilename_AA_CG = *++argv;
  outputfilename_AA = *++argv;
  outputfilename_CG = *++argv;

  numstep=0;
  data=(double *)gcemalloc(sizeof(double)*2);
  inputfile1=efopen(inputfilename,"r");
  data=io_scandcoldata2(inputfile1,numi,numcol,xcol,ycol,&numstep,data);
  fclose(inputfile1);

  TACCMfile=efopen(TACCMfilename,"r");
  fscanf(TACCMfile,"%d",&numZ);
  pairs=(int **)gcemalloc(sizeof(int *)*numZ);
  for (i=0;i<numZ;++i) pairs[i]=(int *)gcemalloc(sizeof(int)*5);
  for (i=0;i<numZ;++i) {
    for (j=0;j<4;++j) fscanf(TACCMfile,"%d",&(pairs[i][j]));
    fscanf(TACCMfile,"%d",&(pairs[i][j]));
  }
  fclose(TACCMfile);

  parmtopfile_FG=efopen(parmtopfilename_FG,"r");
  readParmtopLb(parmtopfile_FG,&ap_FG);
  fclose(parmtopfile_FG); 

  parmtopfile_CG=efopen(parmtopfilename_CG,"r");
  readParmtopLb(parmtopfile_CG,&ap_CG);
  fclose(parmtopfile_CG); 

  for (i=0;i<ap_FG.NATOM;++i) ap_FG.CHRG[i]=0.0;
  for (i=0;i<(ap_FG.NTYPES)*(ap_FG.NTYPES+1)/2;++i){
    ap_FG.CN1[i]=0.0; ap_FG.CN2[i]=0.0;
  }

  for (i=0;i<ap_CG.NATOM;++i) ap_CG.CHRG[i]=0.0;
  for (i=0;i<(ap_CG.NTYPES)*(ap_CG.NTYPES+1)/2;++i){
    ap_CG.CN1[i]=0.0; ap_CG.CN2[i]=0.0;
  }

  numatom=ap_FG.NATOM;

  ffLc_set_calcffandforce(&(e_FG),&(f_FG),ap_FG);
  ffLc_set_calcffandforce(&(e_CG),&(f_CG),ap_CG);

  crd_AA=(double *)gcemalloc(sizeof(double)*numatom*3);
  crd_CG=(double *)gcemalloc(sizeof(double)*numatom*3);
  Z=(double *)gcemalloc(sizeof(double)*numZ);

  if (minmaxflag==OFF) {
    maxx=data[0];minx=data[0];
    maxy=data[1];miny=data[1];
    for (i=0;i<numstep;++i) {
      if (maxx < data[i*2]) maxx=data[i*2];
      if (minx > data[i*2]) minx=data[i*2];
      if (maxy < data[i*2+1]) maxy=data[i*2+1];
      if (miny > data[i*2+1]) miny=data[i*2+1];
    }
  }

  framex=(int)((maxx-minx)/width)+1;
  framey=(int)((maxy-miny)/width)+1;

  pmf_AA_CG=(double *)gcemalloc(sizeof(double)*(framex)*(framey));
  pmf_AA=(double *)gcemalloc(sizeof(double)*(framex)*(framey));
  pmf_CG=(double *)gcemalloc(sizeof(double)*(framex)*(framey));
  for (i=0;i<framex;++i) {
    for (j=0;j<framey;++j) {
      pmf_AA_CG[i*(framey)+j]=0.0;
      pmf_AA[i*(framey)+j]=0.0;
      pmf_CG[i*(framey)+j]=0.0;
    }
  }

  beta=1.0/(k_B*T);

  l=0;

  mync_open_inq_get_sh_MCD(trjfilename_AA,numatom,l,1,l+1,&(nc_id_MD_AA),crd_nc_AA);
  mync_open_inq_get_sh_MCD(trjfilename_CG,numatom,l,1,l+1,&(nc_id_MD_CG),crd_nc_CG);
  trjfileZ=efopen(trjfilenameZ,"r");

  logfile=efopen("logfile_pmf2d_b_wUCGAA.txt","w");
  for (i=0;i<numstep;++i) {
    mync_open_inq_get_sh_MCD(trjfilename_AA,numatom,l,1,l+1,&(nc_id_MD_AA),crd_nc_AA);
    mync_open_inq_get_sh_MCD(trjfilename_CG,numatom,l,1,l+1,&(nc_id_MD_CG),crd_nc_CG);
    ++l;

    for (j=0;j<numZ;++j) fscanf(trjfileZ,"%lf",&Z[j]);
    for (j=0;j<numatom;++j) {
      for (k=0;k<3;++k) {
	crd_AA[j*3+k]=crd_nc_AA[j][k];
	crd_CG[j*3+k]=crd_nc_CG[j][k];
      }
    }

    CE_TACCM_CGAA(crd_AA,crd_CG,Z,numatom,numZ,KZAA,KZCG,pairs,pi,&EAA,&ECG,&EZ);
    ffLc_calcffandforce_AA(crd_AA,numatom,&e_FG,&f_FG,ap_FG);
    ffLc_calcffandforce_AA(crd_CG,numatom,&e_CG,&f_CG,ap_CG);

    if (data[2*i] <= maxx && data[2*i+1] <= maxy) {
      ECG+=e_CG.p_d_t+e_CG.p_a_t+e_CG.p_b_t;
      EAA+=e_FG.p_d_t+e_FG.p_a_t+e_FG.p_b_t;

      expBECG=exp(beta*ECG);
      expBEAA=exp(beta*EAA);
      fprintf(logfile,"%8.4lf %8.4lf\n",EAA,ECG);
      pmf_AA_CG[((int)((data[i*2]-minx)/width))*(framey)+((int)((data[i*2+1]-miny)/width))]+=1.0;
      pmf_AA[((int)((data[i*2]-minx)/width))*(framey)+((int)((data[i*2+1]-miny)/width))]+=expBECG;
      pmf_CG[((int)((data[i*2]-minx)/width))*(framey)+((int)((data[i*2+1]-miny)/width))]+=expBEAA;
    }
  }
  fclose(trjfileZ);
  fclose(logfile);

  sum_AA_CG=0.0; sum_AA=0.0; sum_CG=0.0;
  for (i=0;i<framex;++i){ 
    for (j=0;j<framey;++j) { 
      sum_AA_CG+=pmf_AA_CG[i*(framey)+j];
      sum_AA+=pmf_AA[i*(framey)+j]; sum_CG+=pmf_CG[i*(framey)+j]; 
    }
  }
  for (i=0;i<framex;++i){ 
    for (j=0;j<framey;++j){
      pmf_AA_CG[i*(framey)+j]=pmf_AA_CG[i*(framey)+j]/sum_AA_CG;
      pmf_AA[i*(framey)+j]=pmf_AA[i*(framey)+j]/sum_AA; pmf_CG[i*(framey)+j]=pmf_CG[i*(framey)+j]/sum_CG;
    }
  }

  for (i=0;i<=framex;++i) {
    for (j=0;j<framey;++j) {
      if ( pmf_min_AA_CG == 0.0 && pmf_AA_CG[i*(framey)+j] > 0.0 ) pmf_min_AA_CG=pmf_AA_CG[i*(framey)+j];
      if ( pmf_min_AA_CG > pmf_AA_CG[i*(framey)+j] && pmf_AA_CG[i*(framey)+j] > 0.0 ) 
	pmf_min_AA_CG=pmf_AA_CG[i*(framey)+j];
      if ( pmf_min_AA == 0.0 && pmf_AA[i*(framey)+j] > 0.0 ) pmf_min_AA=pmf_AA[i*(framey)+j];
      if ( pmf_min_AA > pmf_AA[i*(framey)+j] && pmf_AA[i*(framey)+j] > 0.0 ) pmf_min_AA=pmf_AA[i*(framey)+j];
      if ( pmf_min_CG == 0.0 && pmf_CG[i*(framey)+j] > 0.0 ) pmf_min_CG=pmf_CG[i*(framey)+j];
      if ( pmf_min_CG > pmf_CG[i*(framey)+j] && pmf_CG[i*(framey)+j] > 0.0 ) pmf_min_CG=pmf_CG[i*(framey)+j];
    }
  }
  
  for (i=0;i<=framex;++i) {
    for (j=0;j<framey;++j) {
      if (pmf_AA_CG[i*(framey)+j]!=0.0) pmf_AA_CG[i*(framey)+j]=-log(pmf_AA_CG[i*(framey)+j])+log(pmf_min_AA_CG);
      if (pmf_AA[i*(framey)+j]!=0.0) pmf_AA[i*(framey)+j]=-log(pmf_AA[i*(framey)+j])+log(pmf_min_AA);
      if (pmf_CG[i*(framey)+j]!=0.0) pmf_CG[i*(framey)+j]=-log(pmf_CG[i*(framey)+j])+log(pmf_min_CG);
    }
  }

  outputfile_AA_CG=efopen(outputfilename_AA_CG,"w");
  outputfile_AA=efopen(outputfilename_AA,"w");
  outputfile_CG=efopen(outputfilename_CG,"w");
  for (i=0;i<framex;++i) {
    for (j=0;j<framey;++j) {
      fprintf(outputfile_AA_CG,"%e %e %e\n",width*i+minx,width*j+miny,1.0/beta*pmf_AA_CG[i*framey+j]);
      fprintf(outputfile_AA,"%e %e %e\n",width*i+minx,width*j+miny,1.0/beta*pmf_AA[i*framey+j]);
      fprintf(outputfile_CG,"%e %e %e\n",width*i+minx,width*j+miny,1.0/beta*pmf_CG[i*framey+j]);
    }
    fprintf(outputfile_AA_CG,"\n");
    fprintf(outputfile_AA,"\n");
    fprintf(outputfile_CG,"\n");
  }
  fclose(outputfile_AA_CG);
  fclose(outputfile_AA);
  fclose(outputfile_CG);
  
  return 0;
}

double *io_scandcoldata2(FILE *inputfile,int numi,int numcol,int xcol,int ycol,int *numstep,double *data){
  int i,j,k;
  double f;

  char *line;
  size_t len=0;

  for (i=0;i<numi;++i)
    getline(&line,&len,inputfile);

  for (i=(*numstep);;++i) {
    for (j=0;j<numcol;++j) {
      if (fscanf(inputfile,"%lf",&f)!=-1) {
	if (j==xcol-1) {
	  data=(double *)gcerealloc(data,sizeof(double)*(i+1)*2);
	  data[i*2]=f;
	}
	else if (j==ycol-1) {
	  data[i*2+1]=f;
	}
      }
      else {
	*numstep=i-1;
	return data;
      }
    }
  }
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilenam outputfilename\n",progname);
}
