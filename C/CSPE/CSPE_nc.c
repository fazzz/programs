
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <getopt.h>

#include "fftw3.h"

#include "PTL.h"
#include "EF.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

int USAGE(char *progname);

#define kbkcl 1.98723e-3
#define kbuap 1.98723e-3*4.18407*100.0

void FTd(int numdihed, int nstep,double *dtrj, double *spe);
void FT(int natom, int nstep,double *trj, double *spe);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int flag='c',flag2='a',AMBncflag=ON;
  int numatom,numres,numdihed,numstep;
  int num_of_split=1,num_step_ini=0,num_step_fin;
  double c=2.999792e-2;
  double kb=1.98723e-3*4.18407*100.0;
  double deltat,pi,kbT,temp,sumspe,dof;
  double *spe,*trj,*speave;

  double *crd;
  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;
  struct my_netcdf_out_id_AMBER nc_id_MD;

  char *inputfilename,*parmtopfilename,*outputfilename,*progname;
  FILE *inputfile,*parmfile,*outputfile;
  FILE *log;

  int c2;
  extern char *optarg;
  extern int optind,opterr,optopt;

  int opt_idx=1;

  struct option long_opt[] = {
    {"c",0,NULL,'c'},
    {"v",0,NULL,'v'},
    {"CA",0,NULL,'C'},
    {"AA",0,NULL,'A'},
    {"DI",0,NULL,'D'},
    {"AMBERON",0,NULL,'O'},
    {"AMBEROFF",0,NULL,'B'},
    {"dt",1,NULL,'d'},
    {"num_split",1,NULL,'s'},
    {"num_ini",1,NULL,'i'},
    {"temp",1,NULL,'T'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c2=getopt(argc,argv,"hcvCADOBd:s:i:f:T:"))!=-1) {
    switch(c2) {
    case 'c':
      flag='c';
      break;
    case 'v':
      flag='v';
      break;
    case 'C':
      flag2='c';
      break;
    case 'A':
      flag2='a';
      break;
    case 'D':
      flag2='d';
      break;
    case 'O':
      AMBncflag=ON;
      break;
    case 'B':
      AMBncflag=OFF;
      break;
    case 'd':
      deltat=atof(optarg);
      break;
    case 's':
      num_of_split=atoi(optarg);
      break;
    case 'i':
      num_step_ini=atoi(optarg);
      break;
    case 'f':
      num_step_fin=atoi(optarg);
      break;
    case 'T':
      temp=atof(optarg);
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
  
  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  inputfilename   = *argv;
  parmtopfilename = *++argv;
  outputfilename  = *++argv;

  if (AMBncflag==OFF)  num_step_fin=myncL_get_present_step_MCD(inputfilename,&nc_id_MCD);
  else num_step_fin=myncL_get_present_step_AMBER(inputfilename,&nc_id_MD);

  numstep=(num_step_fin-num_step_ini)/num_of_split;

  if (flag=='v')   kbT=kbkcl*temp;
  else  kbT=kb*temp;

  parmfile=efopen(parmtopfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;
  numres=AP.NRES;
  
  speave=(double *)gcemalloc(sizeof(double)*numstep);
  l=num_step_ini;
  for (i=0;i<num_of_split;++i) {
    spe=(double *)gcemalloc(sizeof(double)*numstep);

    if (flag2=='a') {
      trj=(double *)gcemalloc(sizeof(double)*numstep*numatom*3);

      for (j=0;j<numstep;++j) {
	if (AMBncflag==OFF)  myncL_open_inq_get_sh_MCD(inputfilename,numatom,l+i,1,l+i+1,&nc_id_MCD,crd_nc);
	else myncL_open_inq_get_sh_AMBER(inputfilename,numatom,l+i,1,l+i+1,&nc_id_MD,crd_nc);

	for (k=0;k<numatom;++k) for (l=0;l<3;++l) trj[j*numatom*3+k*3+l] = crd_nc[k][l]*sqrt(AP.AMASS[k]);
      }

      FT(numatom,numstep,trj,spe);
    }
    else {
      trj=(double *)gcemalloc(sizeof(double)*numstep*numres*3);

      for (j=0;j<numstep;++j) {
	if (AMBncflag==OFF)  myncL_open_inq_get_sh_MCD(inputfilename,numatom,l+i,1,l+i+1,&nc_id_MCD,crd_nc);
	else myncL_open_inq_get_sh_AMBER(inputfilename,numatom,l+i,1,l+i+1,&nc_id_MD,crd_nc);

	for (k=0;k<numatom;++k) 
	  if (strncmp(AP.IGRAPH[j],"CA",2) == 0)
	    for (l=0;l<3;++l) trj[j*numatom*3+k*3+l] = crd_nc[k][l]*sqrt(AP.AMASS[k]);
      }

      FT(numres,numstep,trj,spe);
    }
    for (j=0;j<numstep;++j)
      speave[j]+=spe[j];
  }
  for (i=0;i<numstep;++i)
    speave[i]=speave[i]/num_of_split;

  pi=acos(-1.0);
    
  dof=0.0;
  log=efopen("log_CSPE.txt","w");
  for (i=0;i<(int)(numstep/2);++i) {
    if (flag=='v') {
      dof+=speave[i]/kbT*2;
    }
    else {
      dof+=(2.0*pi*i/numstep/deltat)*(2.0*pi*i/numstep/deltat)*speave[i]/kbT*2;
    }
  } 
  fprintf(log,"dof=%lf\n",dof);
  fclose(log);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
    if (flag=='v') {
      fprintf(outputfile,"%e %e \n",(double)i/numstep/deltat/c,speave[i]/kbT*2);
    }
    else {
      fprintf(outputfile,"%e %e \n",(double)i/numstep/deltat/c,(2.0*pi*i/numstep/deltat)*(2.0*pi*i/numstep/deltat)*speave[i]/kbT*2);
    }  
  } 
  fclose(outputfile);
  free(speave);

  
  return 0;
}

void FT(int natom, int nstep,double *trj, double *spe){
  int i,j,k,l;
  fftw_complex *in,*out;
  fftw_plan p;
  FILE *outputfile;

  double sumspe,K=0.0,sum=0.0;

  for (i=0;i<nstep;++i)
    spe[i]=0.0;

  for (i=0;i<natom;++i) {
    for (j=0;j<3;++j) {
      in  = fftw_malloc(sizeof(fftw_complex)*nstep);
      out = fftw_malloc(sizeof(fftw_complex)*nstep);
      p   = fftw_plan_dft_1d(nstep,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
      
      for (k=0;k<nstep;++k) {
	in[k][0]=trj[k*natom*3+i*3+j];
	in[k][1]=0.0;
	out[k][0]=0.0;
	out[k][1]=0.0;
      }
      fftw_execute(p);
      
      for (k=0;k<nstep;++k)
	spe[k] += out[k][0]/nstep*out[k][0]/nstep+out[k][1]/nstep*out[k][1]/nstep;
            
      fftw_destroy_plan(p);
      fftw_free(in); fftw_free(out);
    }
  }    

  for (k=0;k<nstep*natom*3;++k)
    K+=trj[k]*trj[k]/nstep;
  for (k=0;k<nstep;++k)
    sum += spe[k];


}

void FTd(int numdihed, int nstep,double *dtrj, double *spe){
  int i,j,k;
  fftw_complex *in,*out;
  fftw_plan p;
  FILE *outputfile;

  for (i=0;i<nstep;++i)
    spe[i]=0.0;

  for (i=0;i<numdihed;++i) {
    in  = fftw_malloc(sizeof(fftw_complex)*nstep);
    out = fftw_malloc(sizeof(fftw_complex)*nstep);
    p   = fftw_plan_dft_1d(nstep,in,out,FFTW_FORWARD,FFTW_ESTIMATE);

    for (j=0;j<nstep;++j)
      for (k=0;k<2;++k)
	in[j][k]=0.0;

    for (j=0;j<nstep;++j)
	in[j][0]=dtrj[j*numdihed+i];
    
    fftw_execute(p);
    
    for (j=0;j<nstep;++j)
      spe[j] += out[j][0]*out[j][0]+out[j][1]*out[j][1];
    
    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);
  }  
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-c] \n");
  printf("[-v] \n");
  printf("[--CA] \n");
  printf("[--AA] \n");
  printf("[--AMBERON] \n");
  printf("[--AMBEROFF] \n");
  printf("[--dt] dt \n");
  printf("[--num_of_split] num_of_split \n");
  printf("[--num_ini] num_ini \n");
  printf("[--temp] temp \n");
  printf("%s inputfilename parmfilename outputfilename \n",progname);
}
