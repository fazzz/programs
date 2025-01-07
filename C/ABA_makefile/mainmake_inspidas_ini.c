#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "netcdf_mine.h"

#include "IO.h"
#include "ECEPE.h"
#include "PROTOPO.h"
#include "PDB.h"
#include "TOPO.h"
#include "PT.h"
#include "EF.h"

#define ON 1
#define OFF 0

int usage(void);

int main(int argc, char *argv[]) {
  int i,j,k,l,d,m,nt,dummy,num,n;
  int flagcn='n',flagrst=OFF;
  int numatom,numstep;
  int interval=1;

  int aromflag=ON;

  double atom[4][3];
  double *protcoord;
  double theta,*old_theta;
  double pi;

  double *co;
  double *dihed,*delta_dihed;
  double dihed_check;
  double *co_dummy,*dihed_dummy,*ene;

  int *bp;
  int **bp_f;
  int *numb;

  int ntotaldih,ndihinres,numres; 
  int *numreslist;

  struct ECEPE_parms ECEPE_p;
  struct pnb nb_p;
  struct ECEPE_pote p;
  struct ECEPE_force f;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double crd_nc[MAXATOM][3],*crd;

  char *line;
  size_t len=0;

  char *name_atom_list;

  char *inputfilename,*outputfilename;
  char *preofilename,*bd8filename,*coofilename_for_sflag,*angfilename,*coofilename;
  char *progname,*parmtopname,*name_atom1,*name_atom;
  FILE *inputfile,*outputfile,*parmtop;
  FILE *preofile,*bd8file,*coofile_for_sflag,*angfile,*coofile;

  while((c=getopt(argc,argv,"h"))!=-1) {
    switch(c) {
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

  if (argc < 7) {
    USAGE(progname);
    exit(1);
  }
  inputfilename  = *argv;
  preofilename   = *++argv;
  bd8filename   = *++argv;
  angfilename = *++argv;
  coofilename = *++argv;
  parmtopname = *++argv;
  outputfilename = *++argv;

  read_ECEPE_parm(preofilename,bd8filename,&ECEPE_p,&nb_p);
  numatom=ECEPE_p.NUMATM;

  parmtop=efopen(parmtopname,"r");
  readParmtop(parmtop);
  fclose(parmtop);

  pi=acos(-1.0);

  bp=(int *)gcemalloc(sizeof(int)*(ECEPE_p.NUMATM-1)*2);
  bp_f=(int **)gcemalloc(sizeof(int *)*ECEPE_p.NUMATM);

  numb=(int *)gcemalloc(sizeof(int)*ECEPE_p.NUMATM);
  name_atom_list=(char *)gcemalloc(sizeof(char)*ECEPE_p.NUMATM*4);
  for (i=0;i<ECEPE_p.NUMATM;++i) 
    for (j=0;j<4;++j)
      name_atom_list[i*4+j]=ECEPE_p.atom[i].name_atom[j];
  make_bp(ECEPE_p.NUMATM,bp,bp_f,numb,name_atom_list,aromflag,OFF);

  make_dihed_pairs_list(ECEPE_p,bp_f,numb);

  co=(double *)gcemalloc(sizeof(double)*(ECEPE_p.NUMATM)*3);
  dihed=(double *)gcemalloc(sizeof(double)*(ECEPE_p.NUMVAR));
  delta_dihed=(double *)gcemalloc(sizeof(double)*(ECEPE_p.NUMVAR));
  for (i=0;i<ECEPE_p.NUMVAR;++i) delta_dihed[i]=0.0;

  co_dummy=(double *)gcemalloc(sizeof(double)*40*3);
  dihed_dummy=(double *)gcemalloc(sizeof(double)*10);
  ene=(double *)gcemalloc(sizeof(double)*6);

  coofile=efopen(coofilename,"r");
  read_ECEPE_coo(coofile,co,dihed_dummy,ECEPE_p.NUMATM);
  fclose(coofile);

  ntotaldih=0;
  ndihinres=1;
  numres=1;
  numreslist=(int *)gcemalloc(sizeof(int)*ECEPE_p.NUMRES);
  for (i=0;i<ECEPE_p.NUMVAR;++i) {
    if (ECEPE_p.dihed[i].indexv1>numres) {
      numres=ECEPE_p.dihed[i].indexv1;
      numreslist[numres-2]=ndihinres;
      ndihinres=1;
    }
    if (ECEPE_p.dihed[i].indexv2>ndihinres) {
      ndihinres=ECEPE_p.dihed[i].indexv2;
    }
  }

  angfile=efopen(angfilename,"r");
  k=0;
  for (i=0;i<ECEPE_p.NUMRES;++i) {
    for (j=0;j<numreslist[i];++j) {
      fscanf(angfile,"%lf",&dihed_dummy[k]);
      ++k;
    }
    for (j=0;j<10-numreslist[i];++j) {
      fscanf(angfile,"%lf",&f);
    }
  }
  fclose(angfile);

  for (i=0;i<ECEPE_p.NUMVAR;++i) {
    dihed_dummy[i]=dihed_dummy[i]*pi/180.0;
    if (dihed_dummy[i]<-pi)
      dihed_dummy[i]+=2.0*pi;
    else if (dihed_dummy[i]>pi)
      dihed_dummy[i]-=2.0*pi;
  }

  calc_TORS_for_get_sabun(ECEPE_p.NUMVAR,co,ECEPE_p,dihed_dummy,delta_dihed);
  
  inputfile=efopen(inputfilename,"r");
  crd=(double *)gcemalloc(sizeof(double )*numatom*3);
  
  name_atom1=(char *)gcemalloc(sizeof(char)*ECEPE_p.NUMATM*4);

  for (i=0;j<4;++j) {
    for (j=0;j<3;++j) {
      crd[i*4+j]=0.0;
    }
  }
  for (j=0;j<3;++j) {
    fscanf(inputfile,"%lf",&crd[4*3+j]); // N
  }
  for (j=0;j<3;++j) {
    fscanf(inputfile,"%lf",&crd[5*3+j]); // H
  }
  //  name_atom1[4]
  for (i=6;i<numatom;++i) {
    for (j=0;j<3;++j) {
      fscanf(inputfile,"%lf",&crd[i*3+j]); // CA~
    }
    for (j=0;j<4;++j) {
      name_atom1[i*4+j]=AP.IGRAPH[i-2][j];
    }
  }

  name_atom=(char *)gcemalloc(sizeof(char)*ECEPE_p.NUMATM*4);

  for (j=0;j<ECEPE_p.NUMATM;++j) {
    for (k=0;k<3;++k) {
      crd_nc[j][k]=crd[(ECEPE_p.atom[j].katom-1)*3+k];
    }
    for (k=0;k<4;++k) {
      name_atom[j*4+k]=name_atom1[(ECEPE_p.atom[j].katom-1)*4+k];
    }
  }
  for (k=0;k<ECEPE_p.NUMVAR;++k) {
     for (m=0;m<4;++m) {
      for (l=0;l<3;++l) {
	atom[m][l]=crd_nc[(ECEPE_p.dihed[k].dpairs[m])][l];
      }
      for (l=0;l<4;++l) {
      	printf("%c",name_atom[(ECEPE_p.dihed[k].dpairs[m])*4+l]);
      }
    }    
    printf("\n");
    for (m=0;m<4;++m) {
      for (l=0;l<4;++l) {
	printf("%c",ECEPE_p.atom[(ECEPE_p.dihed[k].dpairs[m])].name_atom[l]);
      }
      printf(" ");
    }
    printf("\n");
    theta=dih(atom[0],atom[1],atom[2],atom[3]);
    theta+=delta_dihed[k];
    if (theta > pi)  theta-=2.0*pi;
    else if (theta < -pi)  theta+=2.0*pi;
    dihed[k]=theta*180/pi;
  }
  fclose(inputfile);

  outputfile=efopen(outputfilename,"w");
  k=0;
  for (i=0;i<ECEPE_p.NUMRES;++i) {
    for (j=0;j<numreslist[i];++j) {
      fprintf(outputfile,"%8.3lf",dihed[k]);
      ++k;
    }
    for (j=0;j<10-numreslist[i];++j) {
      fprintf(inputfile,"   0.000");
    }
    fprintf(inputfile,"\n");
  }
  fclose(outputfile);

  
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h ] help \n");
  printf("%s inputfilename preofilename bd8filename coofilename_for_sflag outputfilename \n",progname);
}


 


