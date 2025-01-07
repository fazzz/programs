#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "EF.h"

#include "REMDCGAA_TACCM_calc_uene_Amber_PROTEINS2008.h"

#define AAINPF 0
#define CGINPF 1
#define ZINPF  3
#define AAKZ   4
#define CGKZ   5

#define AA1INPF 0
#define CG1INPF 1
#define CG2INPF 2
#define Z1INPF  3
#define AA1KZ   4
#define CG1KZ   5
#define CG2KZ   6

#define WINPF   7

#define ON 0
#define OFF 1

void  CGAAREMDreadInputs_calc_uene(FILE *inputfile,int numatom,int numRE,
				   double *KZAAs,double *KZCGs,
				   char **trjfilenameAA, char **trjfilenameCG,
				   struct my_netcdf_out_id_MCD *nc_id_MD_AA,
				   struct my_netcdf_out_id_MCD *nc_id_MD_CG,
				   FILE **trjfileZ){
  int i,j,k;
  int c,d;
  int numstep;

  double f1,fdummy;

  char /*trjfilenameAA[1000],trjfilenameCG[1000],*/trjfilenameZ[1000];

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
	  trjfilenameAA[i][j]='\0';
	  numstep=mync_get_present_step_MCD(trjfilenameAA[i],&(nc_id_MD_AA[i]));
	  State=CGINPF;
	  j=0;
	}
      }
      else {
	trjfilenameAA[i][j]=c;
	++j;
      }
    }
    else if (State==CGINPF) {
      if (c==' ' || c=='\n') {
	if (j>0) {
	  trjfilenameCG[i][j]='\0';
	  numstep=mync_get_present_step_MCD(trjfilenameCG[i],&(nc_id_MD_CG[i]));
	  State=ZINPF;
	  j=0;
	}
      }
      else {
	trjfilenameCG[i][j]=c;
	++j;
      }
    }
    else if (State==ZINPF) {
      if (c==' ' || c=='\n') {
	if (j>0) {
	  trjfilenameZ[j]='\0';
	  trjfileZ[i]=efopen(trjfilenameZ,"r");
	  State=AAKZ;
	  j=0;
	}
      }
      else {
	trjfilenameZ[j]=c;
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
	  KZAAs[i]=(double)f1;
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
	  KZCGs[i]=(double)f1;
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

void  CGAAREMDreadInputs_pickup_trj(FILE *inputfile, int numRE,
				    char **trjfilename, FILE **trjfile ){
  int i,j,k;
  int c,d;
  int numstep;

  double f1,fdummy;

  int State,Tflag=OFF;

  char *line;
  size_t len=0;

  i=0;  j=0;  f1=0.0;

  while ((c=getc(inputfile))!=-1){
    if (c==' ' || c=='\n') {
      if (j>0) {
	trjfilename[i][j]='\0';
	trjfile[i]=efopen(trjfilename[i],"r");
	++i;
	if (i==numRE) break;
	j=0;
      }
    }
    else {
      trjfilename[i][j]=c;
      ++j;
    }
  }
}

void  CGAAREMDreadInputs_calc_uene_1FG2CG(FILE *inputfile,int numatom,int numRE,
					  double *KZAAs,double *KZCG1s,double *KZCG2s,
					  char **trjfilenameAA, char **trjfilenameCG1, char **trjfilenameCG2,
					  struct my_netcdf_out_id_MCD *nc_id_MD_AA,
					  struct my_netcdf_out_id_MCD *nc_id_MD_CG1,
					  struct my_netcdf_out_id_MCD *nc_id_MD_CG2,
					  FILE **trjfileZ){
  int i,j,k;
  int c,d;
  int numstep;

  double f1,fdummy;

  char /*trjfilenameAA[1000],trjfilenameCG[1000],*/trjfilenameZ[1000];

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
	  trjfilenameAA[i][j]='\0';
	  numstep=mync_get_present_step_MCD(trjfilenameAA[i],&(nc_id_MD_AA[i]));
	  State=CG1INPF;
	  j=0;
	}
      }
      else {
	trjfilenameAA[i][j]=c;
	++j;
      }
    }
    else if (State==CG1INPF) {
      if (c==' ' || c=='\n') {
	if (j>0) {
	  trjfilenameCG1[i][j]='\0';
	  numstep=mync_get_present_step_MCD(trjfilenameCG1[i],&(nc_id_MD_CG1[i]));
	  State=CG2INPF;
	  j=0;
	}
      }
      else {
	trjfilenameCG1[i][j]=c;
	++j;
      }
    }
    else if (State==CG2INPF) {
      if (c==' ' || c=='\n') {
	if (j>0) {
	  trjfilenameCG2[i][j]='\0';
	  numstep=mync_get_present_step_MCD(trjfilenameCG2[i],&(nc_id_MD_CG2[i]));
	  State=ZINPF;
	  j=0;
	}
      }
      else {
	trjfilenameCG2[i][j]=c;
	++j;
      }
    }
    else if (State==ZINPF) {
      if (c==' ' || c=='\n') {
	if (j>0) {
	  trjfilenameZ[j]='\0';
	  trjfileZ[i]=efopen(trjfilenameZ,"r");
	  State=AA1KZ;
	  j=0;
	}
      }
      else {
	trjfilenameZ[j]=c;
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
	  KZAAs[i]=(double)f1;
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
	  KZCG1s[i]=(double)f1;
	  f1=0.0;
	  State=CG2KZ;
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
	  KZCG2s[i]=(double)f1;
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

void  CGAAREMDreadInputs_calc_uene_InV2InW(FILE *inputfile,int numatom,int numRE,
					   double *KZAAs,double *KZCGs,
					   char **trjfilenameAA, char **trjfilenameCG,
					   struct my_netcdf_out_id_MCD *nc_id_MD_AA,
					   struct my_netcdf_out_id_MCD *nc_id_MD_CG,
					   FILE **trjfileZ,FILE **inputwfile){
  int i,j,k;
  int c,d;
  int numstep;

  double f1,fdummy;

  char /*trjfilenameAA[1000],trjfilenameCG[1000],*/trjfilenameZ[1000],inputwfilename[1000];

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
	  trjfilenameAA[i][j]='\0';
	  numstep=mync_get_present_step_MCD(trjfilenameAA[i],&(nc_id_MD_AA[i]));
	  State=CGINPF;
	  j=0;
	}
      }
      else {
	trjfilenameAA[i][j]=c;
	++j;
      }
    }
    else if (State==CGINPF) {
      if (c==' ' || c=='\n') {
	if (j>0) {
	  trjfilenameCG[i][j]='\0';
	  numstep=mync_get_present_step_MCD(trjfilenameCG[i],&(nc_id_MD_CG[i]));
	  State=ZINPF;
	  j=0;
	}
      }
      else {
	trjfilenameCG[i][j]=c;
	++j;
      }
    }
    else if (State==ZINPF) {
      if (c==' ' || c=='\n') {
	if (j>0) {
	  trjfilenameZ[j]='\0';
	  trjfileZ[i]=efopen(trjfilenameZ,"r");
	  State=AAKZ;
	  j=0;
	}
      }
      else {
	trjfilenameZ[j]=c;
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
	  KZAAs[i]=(double)f1;
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
    else  if (State==CGKZ) {
      if (isdigit(c)) {
	d=(c-'0');
	f1=f1*10+(double)d;
	Tflag=ON;
      }
      else if (c==' ' || c=='\n') {
	if (Tflag==ON) {
	  KZCGs[i]=(double)f1;
	  f1=0.0;
	  State=WINPF;
	  Tflag=OFF;
	}
      }
    }
    else {
      if (c==' ' || c=='\n') {
	if (j>0) {
	  inputwfilename[j]='\0';
	  inputwfile[i]=efopen(inputwfilename,"r");
	  State=AA1INPF;
	  ++i;
	  if (i==numRE) break;
	  j=0;
	}
      }
      else {
	inputwfilename[j]=c;
	++j;
      }
    }    
  }
}

void  CGAAREMDreadInputs_calc_uene_1FG2CG_InV2InW(FILE *inputfile,int numatom,int numRE,
						  double *KZAAs,double *KZCG1s,double *KZCG2s,
						  char **trjfilenameAA, char **trjfilenameCG1, char **trjfilenameCG2,
						  struct my_netcdf_out_id_MCD *nc_id_MD_AA,
						  struct my_netcdf_out_id_MCD *nc_id_MD_CG1,
						  struct my_netcdf_out_id_MCD *nc_id_MD_CG2,
						  FILE **trjfileZ,FILE **inputwfile){
  int i,j,k;
  int c,d;
  int numstep;

  double f1,fdummy;

  char /*trjfilenameAA[1000],trjfilenameCG[1000],*/trjfilenameZ[1000],inputwfilename[1000];

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
	  trjfilenameAA[i][j]='\0';
	  numstep=mync_get_present_step_MCD(trjfilenameAA[i],&(nc_id_MD_AA[i]));
	  State=CG1INPF;
	  j=0;
	}
      }
      else {
	trjfilenameAA[i][j]=c;
	++j;
      }
    }
    else if (State==CG1INPF) {
      if (c==' ' || c=='\n') {
	if (j>0) {
	  trjfilenameCG1[i][j]='\0';
	  numstep=mync_get_present_step_MCD(trjfilenameCG1[i],&(nc_id_MD_CG1[i]));
	  State=CG2INPF;
	  j=0;
	}
      }
      else {
	trjfilenameCG1[i][j]=c;
	++j;
      }
    }
    else if (State==CG2INPF) {
      if (c==' ' || c=='\n') {
	if (j>0) {
	  trjfilenameCG2[i][j]='\0';
	  numstep=mync_get_present_step_MCD(trjfilenameCG2[i],&(nc_id_MD_CG2[i]));
	  State=ZINPF;
	  j=0;
	}
      }
      else {
	trjfilenameCG2[i][j]=c;
	++j;
      }
    }
    else if (State==ZINPF) {
      if (c==' ' || c=='\n') {
	if (j>0) {
	  trjfilenameZ[j]='\0';
	  trjfileZ[i]=efopen(trjfilenameZ,"r");
	  State=AA1KZ;
	  j=0;
	}
      }
      else {
	trjfilenameZ[j]=c;
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
	  KZAAs[i]=(double)f1;
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
	  KZCG1s[i]=(double)f1;
	  f1=0.0;
	  State=CG2KZ;
	  Tflag=OFF;
	}
      }
      else {
	printf("error\n");
	exit(1);
      }
    }
    else if (State==CG2KZ)  {
      if (isdigit(c)) {
	d=(c-'0');
	f1=f1*10+(double)d;
	Tflag=ON;
      }
      else if (c==' ' || c=='\n') {
	if (Tflag==ON) {
	  KZCG2s[i]=(double)f1;
	  f1=0.0;
	  State=WINPF;
	  Tflag=OFF;
	}
      }
      else {
	printf("error\n");
	exit(1);
      }
    }
    else {
      if (c==' ' || c=='\n') {
	if (j>0) {
	  inputwfilename[j]='\0';
	  inputwfile[i]=efopen(inputwfilename,"r");
	  State=AA1INPF;
	  ++i;
	  if (i==numRE) break;
	  j=0;
	}
      }
      else {
	inputwfilename[j]=c;
	++j;
      }
    }    
  }    
}
