
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "NC.h"
#include "NCcount.h"
#include "PTL.h"
#include "EF.h"
#include "TOPO.h"
#include "MB.h"

#include "RADG.h"

#include "netcdf_mine.h"

#define AA 0
#define CA 1
#define HA 2

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d;
  int na,nr,resi,resj;
  int numstep;
  double Q,aveQ=0.0,varQ=0.0;
  double RADG,aRADG=0.0,vRADG=0.0;

  double *crd,*crdref,vec[3];
  int numatom,numres;

  double **crdRADG,*mass;
  int numatomRADG;

  int AMBERMODEflag=OFF,MODE=HA,AAMODE=OFF,crdMODE=OFF;
  int nibnum=3;

  double criteria=criteria_NC;
  int numncaa,numncres;
  int  *index_natatt,**ncmap,**ncmapres;
  double atom1[3],atom2[3];

  double nc_ratio=NCratio_default,nc_ext=NCext_default;

  double *cradii_natatt,*cradii_natatt_res;
  double length;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;
  struct my_netcdf_out_id_AMBER nc_id;

  char *line;
  size_t le=0;

  char *progname;
  char *inputfilename,*parmtopname,*crdreffilename;
  char *outputfilename;

  FILE *crdfile,*parmtop,*crdreffile;
  FILE *outputfile;

  int opt_idx=1;

  struct option long_opt[] = {
    {"Amber",0,NULL,'A'},
    {"crd",0,NULL,'d'},
    {"wratio",0,NULL,'1'},
    {"wext",0,NULL,'2'},
    {"AA",0,NULL,'x'},
    {"nibnum",1,NULL,'b'},
    {"criteria",1,NULL,'c'},
    {"NCratio",1,NULL,'N'},
    {"NCext",1,NULL,'C'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hdAx12c:b:N:C:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'A':
      AMBERMODEflag=ON;
      break;
    case 'd':
      crdMODE=ON;
      break;
    case 'x':
      AAMODE=ON;
      break;
    case 'b':
      nibnum=atoi(optarg);
      break;
    case 'c':
      criteria=atof(optarg);
      break;
    case 'N':
      nc_ratio=atof(optarg);
      break;
    case 'C':
      nc_ext=atof(optarg);
      break;
    default:
      USAGE(progname);
      exit(1);
    case 'h':
      USAGE(progname);
      exit(1);
    }
  }

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  inputfilename = *argv;
  parmtopname = *++argv;
  crdreffilename = *++argv;
  outputfilename = *++argv;

  parmtop=efopen(parmtopname,"r");
  readParmtopL(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;
  numres=AP.NRES;

  numatomRADG=AP.NATOM;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  crdref=(double *)gcemalloc(sizeof(double)*numatom*3);

  crdreffile=efopen(crdreffilename,"r");
  io_scanconf_Amber_ini(crdreffile,numatom,crdref);
  fclose(crdreffile);

  if (crdMODE==ON) numstep=1;
  else {
    if (AMBERMODEflag==ON) numstep=mync_get_present_step_AMBER(inputfilename,&nc_id);
    else numstep=mync_get_present_step_MCD(inputfilename,&nc_id_MCD);
  }

  ncmapres=(int **)gcemalloc(sizeof(int *)*numres);
  for (i=0;i<numres;++i) ncmapres[i]=(int *)gcemalloc(sizeof(int)*numres);
  ncmap=(int **)gcemalloc(sizeof(int *)*numatom);
  for (i=0;i<numatom;++i) ncmap[i]=(int *)gcemalloc(sizeof(int)*numatom);

  index_natatt=make_native_contact_list_aa_wnibnum(&numncaa,&numncres,crdref,numatom,numres,criteria,ncmap,ncmapres,EXC,nibnum);
  cradii_natatt=(double *)gcemalloc(sizeof(double)*numncaa);
  cradii_natatt_res=(double *)gcemalloc(sizeof(double)*numncres);

  na=0;
  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      if (ncmap[i][j]==0 ) {
	length = 0.0;
	for(k=0;k<3;++k){
	  vec[k] = crdref[i*3+k]-crdref[j*3+k];
	  length += vec[k]*vec[k];
	}
	length = sqrt(length);
	cradii_natatt[na]=length;
	++na;
      }
    }
  }

  nr=0;
  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      resi=PTL_resnum(i)-1;
      resj=PTL_resnum(j)-1;
      if ( ncmapres[resi][resj]==0 && strcmp(AP.IGRAPH[i],"CA",2)==0 && strcmp(AP.IGRAPH[j],"CA",2)==0) {
	length = 0.0;
	for(k=0;k<3;++k){
	  vec[k] = crdref[i*3+k]-crdref[j*3+k];
	  length += vec[k]*vec[k];
	}
	length = sqrt(length);
	cradii_natatt_res[nr]=length;
	++nr;
      }
    }
  }

  if (MODE==CA) {
    j=0;
    for (i=0;i<numatomRADG;++i) {
      if (strncmp(AP.IGRAPH[i],"CA",2) == 0) {
	++j;
      }
    }
    numatomRADG=j;
  }
  else if (MODE==HA) {
    j=0;
    for (i=0;i<numatomRADG;++i) {
      if (strncmp(AP.IGRAPH[i],"H",1) != 0) {
	++j;
      }
    }
    numatomRADG=j;
  }

  if (MODE==CA) {
    crdRADG=(double **)gcemalloc(sizeof(double *)*numatomRADG);
    for (i=0;i<numatomRADG;++i) {
      crdRADG[i]=(double *)gcemalloc(sizeof(double)*3);
    }
    mass=(double *)gcemalloc(sizeof(double)*numatomRADG);
    j=0;
    for (i=0;i<AP.NATOM;++i) {
      if (strncmp(AP.IGRAPH[i],"CA",2) == 0) {
	mass[j]=AP.AMASS[i];
	++j;
      }
    }
  }
  else if (MODE==HA) {
    crdRADG=(double **)gcemalloc(sizeof(double *)*numatomRADG);
    for (i=0;i<numatomRADG;++i) {
      crdRADG[i]=(double *)gcemalloc(sizeof(double)*3);
    }
    mass=(double *)gcemalloc(sizeof(double)*numatomRADG);
    j=0;
    for (i=0;i<AP.NATOM;++i) {
      if (strncmp(AP.IGRAPH[i],"H",1) != 0) {
	mass[j]=AP.AMASS[i];
	++j;
      }
    }
  }
  else {
    crdRADG=(double **)gcemalloc(sizeof(double *)*numatomRADG);
    for (i=0;i<numatomRADG;++i) {
      crdRADG[i]=(double *)gcemalloc(sizeof(double)*3);
    }
    mass=(double *)gcemalloc(sizeof(double)*numatomRADG);
    for (i=0;i<numatomRADG;++i) {
      mass[i]=AP.AMASS[i];
    }
  }

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
    if (AMBERMODEflag==ON) mync_open_inq_get_sh_AMBER(inputfilename,numatom,i,1,i+1,&nc_id,crd_nc);
    else mync_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id_MCD,crd_nc);
    for (j=0;j<numatom;++j)for (k=0;k<3;++k) crd[j*3+k]=crd_nc[j][k];
    if (MODE==AA) {
      for (j=0;j<numatom;++j)for (k=0;k<3;++k) crdRADG[j][k]=crd_nc[j][k];
    }
    else if (MODE==HA) {
      l=0;
      for (j=0;j<numatom;++j) {
	if (strncmp(AP.IGRAPH[j],"H",1) != 0) {
	  for (k=0;k<3;++k) crdRADG[l][k]=crd_nc[j][k];
	  ++l;
	}
      }
    }
    else if (MODE==CA) {
      l=0;
      for (j=0;j<numatom;++j) {
	if (strncmp(AP.IGRAPH[j],"CA",2) == 0) {
	  for (k=0;k<3;++k) crdRADG[l][k]=crd_nc[j][k];
	  ++l;
	}
      }
    }

    if (AAMODE==OFF) {
      if (MODE==1) 
	Q=NCcount_native_contact_wratio(numncres,crd,numatom,numres,ncmap,ncmapres,cradii_natatt_res,nc_ratio);
      else 
	Q=NCcount_native_contact_wext(numncres,crd,numatom,numres,ncmap,ncmapres,cradii_natatt_res,nc_ext);
    }
    else {
      if (MODE==1) 
	Q=NCcount_native_contact_AA_wratio(numncres,crd,numatom,numres,ncmap,ncmapres,cradii_natatt,nc_ratio);
      else 
	Q=NCcount_native_contact_AA_wext(numncres,crd,numatom,numres,ncmap,ncmapres,cradii_natatt,nc_ext);
    }

    
    aveQ=(i*aveQ+Q)/(i+1);
    varQ=(i*varQ+Q*Q)/(i+1);

    if (MODE==AA)
      RADG=RADG_calc_radg(crdRADG,mass,numatomRADG);
    else if (MODE=CA)
      RADG=RADG_calc_radg(crdRADG,mass,numatomRADG);
    else
      RADG=RADG_calc_radg(crdRADG,mass,numatomRADG);
    
    aRADG=(i*aRADG+RADG)/(i+1);
    vRADG=(i*vRADG+RADG*RADG)/(i+1);

    fprintf(outputfile,"%e %e\n",RADG,Q);
  }
  fclose(outputfile);

  varQ=sqrt(varQ-aveQ*aveQ);

  printf("%s \n",inputfilename);
  printf("%s \n",outputfilename);
  printf("ave= %e var= %e\n",aveQ,varQ);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("[--Amber] \n");
  printf("[--wratio] \n");
  printf("[--wext] \n");
  printf("[--AA] \n");
  printf("[--nibnum] \n");
  printf("[--criteria] \n");
  printf("[--NCratio] \n");
  printf("[--NCext] \n");
  printf("%s [-h] inputfilename parmfilename crdreffilename outputfilename \n",progname);
}
