
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

#include "netcdf_mine.h"

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d;
  int na,nr,resi,resj;
  int numstep;
  double Qa_h,aveQa_h=0.0,varQa_h=0.0;
  double Qh_a,aveQh_a=0.0,varQh_a=0.0;

  double *crd,*crdrefa,*crdrefh,vec[3];
  int numatom,numres;

  int AMBERMODEflag=OFF,MODE=1,AAMODE=OFF,crdMODE=OFF;
  int nibnum=3;

  double criteria=criteria_NC;
  int numncaa_a,numncres_a;
  int numncaa_h,numncres_h;
  int  *index_natatt_a,**ncmap_a,**ncmapres_a;
  int  *index_natatt_h,**ncmap_h,**ncmapres_h;
  double atom1[3],atom2[3];

  double nc_ratio=NCratio_default,nc_ext=NCext_default;

  int **ncmap_a_h,**ncmap_h_a,**ncmapres_a_h,**ncmapres_h_a;
  int numncaa_a_h=0,numncaa_h_a=0;
  int numncres_a_h=0,numncres_h_a=0;
  double *cradii_natatt_a_h,*cradii_natatt_res_a_h;
  double *cradii_natatt_h_a,*cradii_natatt_res_h_a;
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
  char *inputfilename,*parmtopname,*crdrefafilename,*crdrefhfilename;
  char *outputfilename;

  FILE *crdfile,*parmtop,*crdrefafile,*crdrefhfile;
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

  if (argc < 5) {
    USAGE(progname);
    exit(1);
  }
  inputfilename = *argv;
  parmtopname = *++argv;
  crdrefafilename = *++argv;
  crdrefhfilename = *++argv;
  outputfilename = *++argv;

  parmtop=efopen(parmtopname,"r");
  readParmtopL(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;
  numres=AP.NRES;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  crdrefa=(double *)gcemalloc(sizeof(double)*numatom*3);
  crdrefh=(double *)gcemalloc(sizeof(double)*numatom*3);

  crdrefafile=efopen(crdrefafilename,"r");
  io_scanconf_Amber_ini(crdrefafile,numatom,crdrefa);
  fclose(crdrefafile);

  crdrefhfile=efopen(crdrefhfilename,"r");
  io_scanconf_Amber_ini(crdrefhfile,numatom,crdrefh);
  fclose(crdrefhfile);

  if (crdMODE==ON) numstep=1;
  else {
    if (AMBERMODEflag==ON) numstep=mync_get_present_step_AMBER(inputfilename,&nc_id);
    else numstep=mync_get_present_step_MCD(inputfilename,&nc_id_MCD);
  }

  ncmapres_a=(int **)gcemalloc(sizeof(int *)*numres);
  for (i=0;i<numres;++i) ncmapres_a[i]=(int *)gcemalloc(sizeof(int)*numres);
  ncmap_a=(int **)gcemalloc(sizeof(int *)*numatom);
  for (i=0;i<numatom;++i) ncmap_a[i]=(int *)gcemalloc(sizeof(int)*numatom);

  ncmapres_h=(int **)gcemalloc(sizeof(int *)*numres);
  for (i=0;i<numres;++i) ncmapres_h[i]=(int *)gcemalloc(sizeof(int)*numres);
  ncmap_h=(int **)gcemalloc(sizeof(int *)*numatom);
  for (i=0;i<numatom;++i) ncmap_h[i]=(int *)gcemalloc(sizeof(int)*numatom);

  ncmap_a_h=(int **)gcemalloc(sizeof(int *)*numatom);
  for (i=0;i<numatom;++i) ncmap_a_h[i]=(int *)gcemalloc(sizeof(int)*numatom);
  ncmap_h_a=(int **)gcemalloc(sizeof(int *)*numatom);
  for (i=0;i<numatom;++i) ncmap_h_a[i]=(int *)gcemalloc(sizeof(int)*numatom);

  ncmapres_a_h=(int **)gcemalloc(sizeof(int *)*numres);
  for (i=0;i<numres;++i) ncmapres_a_h[i]=(int *)gcemalloc(sizeof(int)*numres);
  ncmapres_h_a=(int **)gcemalloc(sizeof(int *)*numres);
  for (i=0;i<numres;++i) ncmapres_h_a[i]=(int *)gcemalloc(sizeof(int)*numres);


  index_natatt_a=make_native_contact_list_aa_wnibnum(&numncaa_a,&numncres_a,crdrefa,numatom,numres,criteria,ncmap_a,ncmapres_a,EXC,nibnum);

  index_natatt_h=make_native_contact_list_aa_wnibnum(&numncaa_h,&numncres_h,crdrefh,numatom,numres,criteria,ncmap_h,ncmapres_h,EXC,nibnum);

  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
	ncmap_h_a[i][j]=-2;
	ncmap_a_h[i][j]=-2;
    }
  }

  for (i=0;i<numres;++i) {
    for (j=i+1;j<numres;++j) {
	ncmapres_a_h[i][j]=-2;
	ncmapres_h_a[i][j]=-2;
    }
  }

  numncaa_h_a=0;
  numncaa_a_h=0;
  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      if (ncmap_a[i][j]!=0 && ncmap_h[i][j]==0) {
	ncmap_h_a[i][j]=0;
	++numncaa_h_a;
      }
      else if (ncmap_a[i][j]==0 && ncmap_h[i][j]!=0) {
	ncmap_a_h[i][j]=0;
	++numncaa_a_h;
      }
    }
  }

  numncres_h_a=0;
  numncres_a_h=0;
  for (i=0;i<numres;++i) {
    for (j=i+1;j<numres;++j) {
      if (ncmapres_a[i][j]!=0 && ncmapres_h[i][j]==0) {
	ncmapres_h_a[i][j]=0;
	++numncres_h_a;
      }
      else if (ncmapres_a[i][j]==0 && ncmapres_h[i][j]!=0) {
	ncmapres_a_h[i][j]=0;
	++numncres_a_h;
      }
    }
  }

  cradii_natatt_a_h=(double *)gcemalloc(sizeof(double)*numncaa_a_h);
  cradii_natatt_res_a_h=(double *)gcemalloc(sizeof(double)*numncres_a_h);
  cradii_natatt_h_a=(double *)gcemalloc(sizeof(double)*numncaa_h_a);
  cradii_natatt_res_h_a=(double *)gcemalloc(sizeof(double)*numncres_h_a);

  na=0;
  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      if (ncmap_a_h[i][j]==0 ) {
	length = 0.0;
	for(k=0;k<3;++k){
	  vec[k] = crdrefa[i*3+k]-crdrefa[j*3+k];
	  length += vec[k]*vec[k];
	}
	length = sqrt(length);
	cradii_natatt_a_h[na]=length;
	++na;
      }
    }
  }
  numncaa_a_h=na;

  na=0;
  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      if (ncmap_h_a[i][j]==0 ) {
	length = 0.0;
	for(k=0;k<3;++k){
	  vec[k] = crdrefh[i*3+k]-crdrefh[j*3+k];
	  length += vec[k]*vec[k];
	}
	length = sqrt(length);
	cradii_natatt_h_a[na]=length;
	++na;
      }
    }
  }
  numncaa_h_a=na;

  nr=0;
  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      resi=PTL_resnum(i)-1;
      resj=PTL_resnum(j)-1;
      if ( ncmapres_a_h[resi][resj]==0 && strcmp(AP.IGRAPH[i],"CA",2)==0 && strcmp(AP.IGRAPH[j],"CA",2)==0) {
	length = 0.0;
	for(k=0;k<3;++k){
	  vec[k] = crdrefa[i*3+k]-crdrefa[j*3+k];
	  length += vec[k]*vec[k];
	}
	length = sqrt(length);
	cradii_natatt_res_a_h[nr]=length;
	++nr;
      }
    }
  }
  numncres_a_h=nr;

  nr=0;
  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      resi=PTL_resnum(i)-1;
      resj=PTL_resnum(j)-1;
      if ( ncmapres_h_a[resi][resj]==0 && strcmp(AP.IGRAPH[i],"CA",2)==0 && strcmp(AP.IGRAPH[j],"CA",2)==0) {
	length = 0.0;
	for(k=0;k<3;++k){
	  vec[k] = crdrefh[i*3+k]-crdrefh[j*3+k];
	  length += vec[k]*vec[k];
	}
	length = sqrt(length);
	cradii_natatt_res_h_a[nr]=length;
	++nr;
      }
    }
  }
  numncres_h_a=nr;

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
    if (crdMODE==OFF) {
      if (AMBERMODEflag==ON) mync_open_inq_get_sh_AMBER(inputfilename,numatom,i,1,i+1,&nc_id,crd_nc);
      else mync_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id_MCD,crd_nc);
      for (j=0;j<numatom;++j)for (k=0;k<3;++k)crd[j*3+k]=crd_nc[j][k];
    }
    else {
      crdfile=efopen(inputfilename,"r");
      getline(&line,&le,crdfile);
      fscanf(crdfile,"%d",&d);
      for (j=0;j<numatom;++j) for (k=0;k<3;++k) fscanf(crdfile,"%lf",&crd[j*3+k]);
      fclose(crdfile);
    }

    if (AAMODE==OFF) {
      if (MODE==1) {
	Qa_h=NCcount_native_contact_wratio(numncres_a_h,crd,numatom,numres,
					   ncmap_a_h,ncmapres_a_h,
					   cradii_natatt_res_a_h,
					   nc_ratio);
	Qh_a=NCcount_native_contact_wratio(numncres_h_a,crd,numatom,numres,
					   ncmap_h_a,ncmapres_h_a,
					   cradii_natatt_res_h_a,
					   nc_ratio);
      }
      else {
	Qa_h=NCcount_native_contact_wext(numncres_a_h,crd,numatom,numres,
					 ncmap_a_h,ncmapres_a_h,
					 cradii_natatt_res_a_h,
					 nc_ext);
	Qh_a=NCcount_native_contact_wext(numncres_h_a,crd,numatom,numres,
					 ncmap_h_a,ncmapres_h_a,
					 cradii_natatt_res_h_a,
					 nc_ext);
      }
    }
    else {
      if (MODE==1) {
	Qa_h=NCcount_native_contact_AA_wratio(numncres_a_h,crd,numatom,numres,
					      ncmap_a_h,ncmapres_a_h,
					      cradii_natatt_a_h,nc_ratio);
	Qh_a=NCcount_native_contact_AA_wratio(numncres_h_a,crd,numatom,numres,
					      ncmap_h_a,ncmapres_h_a,
					      cradii_natatt_h_a,nc_ratio);
      }
      else {
	Qa_h=NCcount_native_contact_AA_wext(numncres_a_h,crd,numatom,numres,
					    ncmap_a_h,ncmapres_a_h,
					    cradii_natatt_a_h,nc_ext);
	Qh_a=NCcount_native_contact_AA_wext(numncres_h_a,crd,numatom,numres,
					    ncmap_h_a,ncmapres_h_a,
					    cradii_natatt_h_a,nc_ext);
      }
    }
    
    aveQa_h=(i*aveQa_h+Qa_h)/(i+1);
    varQa_h=(i*varQa_h+Qa_h*Qa_h)/(i+1);
    aveQh_a=(i*aveQh_a+Qh_a)/(i+1);
    varQh_a=(i*varQh_a+Qh_a*Qh_a)/(i+1);
    fprintf(outputfile,"%e %e\n",Qh_a,Qa_h);
  }
  fclose(outputfile);

  varQa_h=sqrt(varQa_h-aveQa_h*aveQa_h);
  varQh_a=sqrt(varQh_a-aveQh_a*aveQh_a);

  printf("%s \n",inputfilename);
  printf("%s \n",outputfilename);
  printf("ave_a_h= %e var_a_h= %e ave_h_a= %e var_h_a= %e\n",aveQa_h,varQa_h,aveQh_a,varQh_a);

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
  printf("%s [-h] inputfilename parmfilename crdrefafilename outputfilename \n",progname);
}
