#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "ECEPE.h"
#include "EF.h"

#define numatommasstype 5

char *atomdata[100]= {
  "C",
  "H",
  "N",
  "O",
  "S"
};
  
double massdata[100]= {
  1.20100000E+01,
  1.00800000E+00,
  1.40100000E+01,
  1.60000000E+01,
  3.20600000E+01
};

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k;
  int num;

  double *mass;

  char *progname;

  char *preofilename;
  char *bd8filename;
  char *massfilename;
  FILE *massfile;

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double *f;
  gchar *atomnameaskey;
  gdouble *atommass;
  GHashTable *masslist;

  char *nameatom;
  struct ECEPE_parms ECEPE_p;
  struct pnb nb_p;

  masslist=g_hash_table_new_full(g_str_hash,g_str_equal,g_free,g_free);

  for (i=0;i<numatommasstype;++i) {
    atomnameaskey=g_strdup(atomdata[i]);
    atommass=g_new(gdouble,1);atommass=&massdata[i];
    g_hash_table_insert(masslist,atomnameaskey,atommass);
    if ((f=g_hash_table_lookup(masslist,atomdata[i]))!=NULL) {
      //      printf("%s - %lf\n",atomdata[i],*f);
      ;
    }
  }

  progname=argv[0];
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

  argc-=optind;
  argv+=optind;

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  preofilename  = *argv;
  bd8filename   = *++argv;
  massfilename  = *++argv;

  read_ECEPE_parm(preofilename,bd8filename,&ECEPE_p,&nb_p);

  //  mass=(double *)gcemalloc(sizeof(double)*ECEPE_p.NUMATM);

  nameatom=(char *)gcemalloc(sizeof(char)*1);

  massfile=efopen(massfilename,"w");
  for (i=0;i<ECEPE_p.NUMATM;++i) {
    for (j=0;j<ECEPE_p.NUMATM;++j) {
	if (i==ECEPE_p.atom[j].katom-1) {
	  num=j;
	  break;
	}
    }
    nameatom[0]=ECEPE_p.atom[num].name_atom[0];
    if ((mass=g_hash_table_lookup(masslist,nameatom))==NULL) {
      printf("There is no key %s in hashtable",nameatom);
      exit(1);
    }
    printf("%lf  \n",*mass);
    fprintf(massfile,"%lf  \n",*mass);
  }
  fclose(massfile);

  /******************************************/
  /* for (i=0;i<ECEPE_p.NUMATM;++i) {	    */
  /*   fprintf(massfile,"%lf  \n",mass[i]); */
  /* }					    */
  /* printf("\n");			    */
  /******************************************/
  
  return 0;
}

void USAGE(char *progname) {
  printf("-h -- help\n");
  printf("USAGE: %s profilename bd8filename massfilename\n", progname);
}

