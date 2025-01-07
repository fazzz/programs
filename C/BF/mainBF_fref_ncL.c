
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "bestfit.h"
#include "PTL.h"
#include "EF.h"

#include "netcdf_mineL.h"

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,s,a;
  int numatom,numCAatom,numstep,numinterval=1;
  int MODE=AA,NCMODE=MD,INMODE=AA;

  double *mass,massCA,*crd,*refcrd,rmsd,*rmsd_trj,f[3];

  double crd_nc[MAXATOM][3];

  char *inputfilename,*parmtopfilename,*refcrdfilename;
  char *outputfilename,*outputtrjfilename;
  char *progname;

  FILE *inputfile, *parmtop,*refcrdfile;
  FILE *outputfile;

  struct my_netcdf_out_id_MCD nc_id_MCD;
  struct my_netcdf_out_id_AMBER nc_id_MD;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  int opt_idx=1;

  struct option long_opt[] = {
    {"Amber",0,NULL,'A'},
    {"CA",0,NULL,'C'},
    {"H",0,NULL,'H'},
    {"INCA",0,NULL,'I'},
    {"interval",1,NULL,'i'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"ACHIhi:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'A':
      NCMODE=AA; break;
    case 'C':
      MODE=CA; break;
    case 'H':
      MODE=HV; break;
    case 'I':
      INMODE=CA; break;
    case 'i':
      numinterval=atoi(optarg); break;
    case 'h':
      USAGE(progname);   break;
    default:
      USAGE(progname);  exit(1);
    }
  }

  progname=argv[0];
  argc-=optind;
  argv+=optind;

  if (argc < 5) {
    USAGE(progname);
    exit(1);
  }
  inputfilename      =  *argv;
  parmtopfilename    =  *++argv;
  refcrdfilename     =  *++argv;
  outputfilename     =  *++argv;
  outputtrjfilename  =  *++argv;

  parmtop=efopen(parmtopfilename,"r");
  readParmtopL(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;

  if (INMODE==CA) {
    j=0;
    for (i=0;i<numatom;++i) {
      if (strncmp(AP.IGRAPH[i],"CA",2)==0){
	++j; massCA=AP.AMASS[i];
      }
    }
    numatom=j;
  }

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  refcrd=(double *)gcemalloc(sizeof(double)*numatom*3);

  refcrdfile=efopen(refcrdfilename,"r");
  getline(&line,&len,refcrdfile);
  fscanf(refcrdfile,"%d",&a);
  k=0;
  for (i=0;i<AP.NATOM;++i) {
    for (j=0;j<3;++j) fscanf(refcrdfile,"%lf",&f[j]);
    if (INMODE==CA) {
      if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
	for (j=0;j<3;++j) refcrd[k*3+j]=f[j];
	++k;
      }
    }
    else { 
      for (j=0;j<3;++j) refcrd[i*3+j]=f[j];
    }
  }
  fclose(refcrdfile);

  mass=(double *)gcemalloc(sizeof(double)*numatom);
  if (INMODE==CA) for (i=0;i<numatom;++i) mass[i] = massCA;
  else for (i=0;i<numatom;++i) mass[i] = AP.AMASS[i];

  if (NCMODE==AMBER) numstep=myncL_get_present_step_MCD(inputfilename,&nc_id_MD);
  else numstep=myncL_get_present_step_MCD(inputfilename,&nc_id_MCD);

  /************************************************************************************************************/
  /* outputfile=efopen(outputfilename,"w");								      */
  /* 													      */
  /* for (i=0;i<numstep;++i) {										      */
  /*   if (NCMODE==AMBER) myncL_open_inq_get_sh_AMBER(inputfilename,numatom,i,1,i+1,&nc_id_MD,crd_nc);	      */
  /*   else myncL_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id_MCD,crd_nc);			      */
  /* 													      */
  /*   l=0;												      */
  /*   for (j=0;j<numatom;++j) {									      */
  /*     if (MODE==CA) {										      */
  /* 	if (strncmp(AP.IGRAPH[i],"CA",2)==0) {								      */
  /* 	  for (k=0;k<3;++k) crd[l*3+k]=crd_nc[j][k];							      */
  /* 	}												      */
  /* 	++l;												      */
  /*     }												      */
  /*     else {												      */
  /* 	for (k=0;k<3;++k) crd[j*3+k]=crd_nc[j][k];							      */
  /*     }												      */
  /*   }												      */
  /* 													      */
  /*   rmsd=bestfit_nc(refcrd,crd,mass,numatom);							      */
  /*   fprintf(outputfile,"%12.8lf\n",rmsd);								      */
  /* }													      */
  /* 													      */
  /* fclose(outputfile);										      */
  /************************************************************************************************************/

  rmsd_trj =(double *)gcemalloc(sizeof(double)*numstep);

  bf_trajectry_ncbL(numatom,numstep,numinterval,mass,refcrd,rmsd_trj,
		    inputfilename,outputtrjfilename,MODE,"x",NCMODE,INMODE);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) fprintf(outputfile,"%12.8lf\n",rmsd_trj[i]);
  fclose(outputfile);

   return 0;
}

void USAGE(char *progname) {
  printf("-C [alpha carbon] \n");
  printf("-H [exclude hydrogen] \n");
  printf("-K [Amber type] \n");
  printf("-h [help] \n");
  printf("%s inputfilename parmtopfilename outputfilenamebase\n",progname);
}
