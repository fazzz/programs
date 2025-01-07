#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "GOLMAA.h"
#include "GOLMAA_set.h"
#include "GOLMAA_check.h"

#include "PTL.h"
#include "EF.h"
#include "PDB.h"
#include "NC.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,ii,jj;
  int numatom,numres;
  int pdbflag=OFF;
  int flagd=ON,flagn=ON,flagr=ON;

  int flagfrccheck=OFF,flagnatcheck=OFF,flagadd=OFF;

  int **nb_matrix;
  double constant=1.0;

  int numspatom=11;
  double dx=0.00001;
  double *f,f_n[3],f_r[3];
  
  double *crd,*refcrd;
  struct potential_GOLMAA e;
  double R_C_D=1.0;

  PDBF PDB,PDBref;

  /****************************/
  /* int num_NCcrd,num_NCref; */
  /* int *ncmapcrd,*ncmapref; */
  /****************************/
  int numnc,*indexncb;
  int **ncmap;
  double Q_NC;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*reffilename,*parmfilename,*outputfilename;
  char *progname;
  FILE *inputfile,*parmfile,*reffile,*outputfile;

  int opt_idx=1;

  struct option long_opt[] = {
    {"pdb",0,NULL,'p'},
    {"f",0,NULL,'f'},
    {"nc",0,NULL,'c'},
    {"a",0,NULL,'a'},
    {"d",0,NULL,'d'},
    {"n",0,NULL,'n'},
    {"r",0,NULL,'r'},
    {"s",1,NULL,'s'},
    {"dx",1,NULL,'x'},
    {"rate",1,NULL,'|'},
    {"cons",1,NULL,'>'},
    {"H",0,NULL,'H'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"pfcadnrHs:|:x:>:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'p':
      pdbflag=ON;
      break;
    case 'f':
      flagfrccheck=ON;
      break;
    case 'c':
      flagnatcheck=ON;
      break;
    case 'a':
      flagadd=ON;
      break;
    case 'd':
      flagd=OFF;
      break;
    case 'x':
      dx=atof(optarg);
      break;
    case 's':
      numspatom=atoi(optarg);
      break;
    case '|':
      R_C_D=atof(optarg);
      break;
    case '>':
      constant=atof(optarg);
      break;
    case 'H':
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
  reffilename = *++argv;
  parmfilename = *++argv;
  outputfilename = *++argv;

  parmfile=efopen(parmfilename,"r");  
  readParmtopL(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;
  numres=AP.NRES;
  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  refcrd=(double *)gcemalloc(sizeof(double)*numatom*3);
  if (pdbflag==ON) {
    PDB.PDBa=(PDBA *)gcemalloc(sizeof(PDBA)*numatom);
  }

  inputfile=efopen(inputfilename,"r");
  if (pdbflag==OFF) {
    io_scanconf_Amber_ini(inputfile,numatom,crd);
    fclose(inputfile);
  }
  else {
    readPDB(inputfile,PDB,numatom);
    for (i=0;i<numatom;++i)
      for (j=0;j<3;++j)
	crd[i*3+j]=PDB.PDBa[i].coord[j];
  }

  reffile=efopen(reffilename,"r");
  if (pdbflag==OFF) {
    io_scanconf_Amber_ini(reffile,numatom,refcrd);
    fclose(reffile);
  }
  else {
    readPDB(reffile,PDB,numatom);
    for (i=0;i<numatom;++i)
      for (j=0;j<3;++j)
	refcrd[i*3+j]=PDB.PDBa[i].coord[j];
  }

  /***********************************************************************************/
  /* ncmapcrd=(int *)gcemalloc(sizeof(int)*numatom*numatom);			     */
  /* ncmapref=(int *)gcemalloc(sizeof(int)*numatom*numatom);			     */
  /* make_native_contact_list_aa(&num_NCref,refcrd,numatom,criteria_NC,ncmapref,ON); */
  /* make_native_contact_list_aa(&num_NCcrd,crd,numatom,criteria_NC,ncmapcrd,ON);    */
  /* Q_NC=((double)num_NCcrd)/((double)num_NCref);				     */
  /***********************************************************************************/

  ncmap=(int **)gcemalloc(sizeof(int *)*numres);
  for (i=0;i<numres;++i) ncmap[i]=(int *)gcemalloc(sizeof(int)*numres);
  indexncb=make_native_contact_list(&numnc,refcrd,numatom,numres,criteria_NC,ncmap/*,EXC*/);
  Q_NC=count_native_contact(numnc,crd,numatom,numres,indexncb,criteria_NC,EXC);

  nb_matrix=(int **)gcemalloc(sizeof(int *)*numatom);
  GOLMAAff_set_calcff(&e,refcrd,numatom,nb_matrix,R_C_D,constant);

  if (flagadd==OFF)
    outputfile=efopen(outputfilename,"w");
  else
    outputfile=efopen(outputfilename,"a");

  if (flagnatcheck==ON) {
    for (i=0;i<e.num_natatt;++i) {
      for (j=0;j<1;++j) {
	fprintf(outputfile,"%c ",AP.IGRAPH[e.index_natatt[i*2]][j]);
      }
      ii=PTL_resnum(e.index_natatt[i*2]);
      fprintf(outputfile,"(%-3d_%-3d) - ",e.index_natatt[i*2]+1,ii);
      for (j=0;j<1;++j) {
	fprintf(outputfile,"%c ",AP.IGRAPH[e.index_natatt[i*2+1]][j]);
      }
      jj=PTL_resnum(e.index_natatt[i*2+1]);
      fprintf(outputfile,"(%-3d_%-3d)\n",e.index_natatt[i*2+1]+1,jj);
    }
  }

  GOLMAAff_calcff(crd,numatom,&e,flagd,flagn,flagr,nb_matrix);
  /***************************************************************************************/
  /* fprintf(outputfile,"p_t=%8.3lf \n",e.p_t);						 */
  /* fprintf(outputfile,"p_b=%8.3lf p_a=%8.3lf  p_d=%8.3lf \n",e.p_b_t,e.p_a_t,e.p_d_t); */
  /* fprintf(outputfile,"p_n=%8.3lf p_r=%8.3lf \n",e.p_natatt_t,e.p_repul_t);		 */
  /***************************************************************************************/

  fprintf(outputfile,"%8.3lf %8.3lf %8.3lf %8.3lf %8.3lf \n",Q_NC,e.p_t,e.p_d_t,e.p_natatt_t,e.p_repul_t);

  if (flagfrccheck==ON) {
    fprintf(outputfile,"f_x=%8.3lf f_y=%8.3lf f_z=%8.3lf \n",e.f_t[numspatom][0],e.f_t[numspatom][1],e.f_t[numspatom][2]);
    fprintf(outputfile,"f_x=%8.3lf f_y=%8.3lf f_z=%8.3lf \n",e.f_natatt[numspatom][0],e.f_natatt[numspatom][1],e.f_natatt[numspatom][2]);
    fprintf(outputfile,"f_x=%8.3lf f_y=%8.3lf f_z=%8.3lf \n",e.f_repul[numspatom][0],e.f_repul[numspatom][1],e.f_repul[numspatom][2]);
  }

  printf("%8.3lf %8.3lf %8.3lf %8.3lf %8.3lf \n",Q_NC,e.p_t,e.p_d_t,e.p_natatt_t,e.p_repul_t);

  if (flagfrccheck==ON) {
    f=GOLMAAff_calcff_check(crd,numatom,&e,numspatom,dx,f_n,f_r,nb_matrix);

    printf("f_t_x=%8.3lf f_t_y=%8.3lf f_t_z=%8.3lf \n",e.f_t[numspatom][0],e.f_t[numspatom][1],e.f_t[numspatom][2]);
    printf("f_t_x=%8.3lf f_t_y=%8.3lf f_t_z=%8.3lf \n",f[0],f[1],f[2]);

    printf("f_n_x=%8.3lf f_n_y=%8.3lf f_n_z=%8.3lf \n",e.f_natatt[numspatom][0],e.f_natatt[numspatom][1],e.f_natatt[numspatom][2]);
    printf("f_n_x=%8.3lf f_n_y=%8.3lf f_n_z=%8.3lf \n",f_n[0],f_n[1],f_n[2]);
    
    printf("f_r_x=%8.3lf f_r_y=%8.3lf f_r_z=%8.3lf \n",e.f_repul[numspatom][0],e.f_repul[numspatom][1],e.f_repul[numspatom][2]);
    printf("f_r_x=%8.3lf f_r_y=%8.3lf f_r_z=%8.3lf \n",f_r[0],f_r[1],f_r[2]);
  }

  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-help]\n");
  printf("%s inputfilename reffilename parmfilename\n",progname);
}

