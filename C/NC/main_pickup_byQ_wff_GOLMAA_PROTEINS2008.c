
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "PTL.h"

#include "NC.h"
#include "NCcount.h"
#include "EF.h"
#include "TOPO.h"
#include "MB.h"

#include "PDBL.h"

#include "netcdf_mine.h"

#include "FFL.h"
#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008.h"

#define Netcdf 0
#define PDBf 1

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,dummy;
  int na,nr,resi,resj;
  int numstep;
  double Q;
  double Qmin=0.1,Qmax=0.9;

  double *crd,*crdref,vec[3],*mass,summass;
  double COM[3];
  int numatom,numres;

  int AMBERMODEflag=OFF,MODE=1,AAMODE=OFF,crdMODE=OFF;
  int nibnum=3;

  double criteria=criteria_NC;
  int numncaa,numncres;
  int  *index_natatt,**ncmap,**ncmapres;
  double atom1[3],atom2[3];

  double ep=0.3;
  int NCmode=3;
  double d=1.0,de=1.0,d2;

  double nc_ratio=NCratio_default,nc_ext=NCext_default;

  double *cradii_natatt,*cradii_natatt_res;
  double length;

  PDBLF PDBL;

  int outMODE=NC,outtypeMode=HV;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;
  struct my_netcdf_out_id_AMBER nc_id;

  struct potential e;
  struct force f;
  struct potential_GOLMAA_PROTEINS2008 e_GOLM;

  char *line;
  size_t le=0;

  char *progname;
  char *inputfilename,*parmtopname,*crdreffilename;
  char *outputfilename,*outputfilename2,outputfilenamepdb[1000];

  FILE *crdfile,*parmtop,*crdreffile;
  FILE *outputfile,*outputfile2;

  int opt_idx=1;

  struct option long_opt[] = {
    {"Amber",0,NULL,'A'},
    {"crd",0,NULL,'d'},
    {"Qmin",1,NULL,'i'},
    {"Qmax",1,NULL,'a'},
    {"wratio",0,NULL,'1'},
    {"wext",0,NULL,'2'},
    {"AA",0,NULL,'x'},
    {"nibnum",1,NULL,'b'},
    {"criteria",1,NULL,'c'},
    {"NCratio",1,NULL,'N'},
    {"NCext",1,NULL,'C'},
    {"nc",0,NULL,'n'},
    {"pdb",0,NULL,'p'},
    {"AA",0,NULL,'l'},
    {"HV",0,NULL,'v'},
    {"CA",0,NULL,'k'},
    {"ep",1,NULL,'j'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hdAxnplvk12c:b:N:C:i:a:j:e:",long_opt,&opt_idx))!=-1) {
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
    case 'n':
      outMODE=Netcdf;
      break;
    case 'p':
      outMODE=PDBf;
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
    case 'i':
      Qmin=atof(optarg);
      break;
    case 'a':
      Qmax=atof(optarg);
      break;
    case 'l':
      outtypeMode=AA;
      break;
    case 'v':
      outtypeMode=HV;
      break;
    case 'k':
      outtypeMode=CA;
      break;
    case 'j':
      ep=atof(optarg);
      break;
    case 'q':
      de=atof(optarg);
      break;
    case 'g':
      d=atof(optarg);
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
  crdreffilename = *++argv;
  outputfilename = *++argv;
  outputfilename2= *++argv;

  parmtop=efopen(parmtopname,"r");
  readParmtopL(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;
  numres=AP.NRES;
  mass=(double *)gcemalloc(sizeof(double)*numatom);
  for (i=0;i<numatom;++i) mass[i]=AP.AMASS[i];

  PDBL.numatom=numatom;
  PDBL.PDBLa=(PDBLA *)gcemalloc(sizeof(PDBLA)*numatom);
  readPDBLdatafromParmtop(PDBL);

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

  ffL_set_calcffandforce(&e,&f);
  GOLMAA_PROTEINS2008_ff_set_calcff(&e_GOLM,crdref,numatom,numres,e.parm.indexnb,e.parm.numnb,ep,nibnum,criteria);
  GOLMAA_PROTEINS2008_ff_calcff(crdref,numatom,&e_GOLM);

  summass=0.0; for (j=0;j<numatom;++j) summass+=mass[j];
  for (j=0;j<3;++j) COM[j]=0.0;
  for (j=0;j<numatom;++j) 
    for (k=0;k<3;++k) 
      COM[k]+=mass[j]*crdref[j*3+k]/summass;
  for (j=0;j<numatom;++j) 
    for (k=0;k<3;++k) 
      crdref[j*3+k]-=COM[k];

  GOLMAA_PROTEINS2008_ff_calcff(crdref,numatom,&e_GOLM);

  if (outMODE==Netcdf)
    myncL_create_def_AMBER(outputfilename,numatom,&nc_id);

  outputfile2=efopen(outputfilename2,"w");
  l=0;
  for (i=0;i<numstep;++i) {
    if (crdMODE==OFF) {
      if (AMBERMODEflag==ON) mync_open_inq_get_sh_AMBER(inputfilename,numatom,i,1,i+1,&nc_id,crd_nc);
      else mync_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id_MCD,crd_nc);
      for (j=0;j<numatom;++j)for (k=0;k<3;++k)crd[j*3+k]=crd_nc[j][k];
    }
    else {
      crdfile=efopen(inputfilename,"r");
      getline(&line,&le,crdfile);
      fscanf(crdfile,"%d",&dummy);
      for (j=0;j<numatom;++j) for (k=0;k<3;++k) fscanf(crdfile,"%lf",&crd[j*3+k]);
      fclose(crdfile);
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

    if ( Q >= Qmin && Q <= Qmax) {
      GOLMAA_PROTEINS2008_ff_calcff_wobaimp(crd,numatom,&e_GOLM);
      fprintf(outputfile2,"%3d p_t=%8.3e p_n=%8.3e p_r=%8.3e p_d=%8.3e\n"
	      ,i,e_GOLM.p_t,e_GOLM.p_natatt_t,e_GOLM.p_repul_t,e_GOLM.p_d_t);

      if (outMODE==Netcdf) {
	for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crd[j*3+k];
	myncL_put_crd_AMBER(nc_id,l,crd_nc);
	++l;
      }
      else if (outMODE==PDBf) {
	++l;
	for (j=0;j<numatom;++j) {
	  for (k=0;k < 3; ++k)	{ 
	    PDBL.PDBLa[j].coord[k]=crd[j*3+k];
	  }
	}
	sprintf(outputfilenamepdb,"%s_%d.pdb",outputfilename,l);
	outputfile=efopen(outputfilenamepdb,"w");
	fprintf(outputfile,"MODEL\n");
	//	writPDB/*_wopt*/(outputfile,PDB/*,HV*/);
	writPDBL_wopt(outputfile,PDBL,outtypeMode);
	fprintf(outputfile,"ENDMOD\n");
	fclose(outputfile);
      }
    }
  }

  fclose(outputfile2);
  printf("%d \n",l);

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
