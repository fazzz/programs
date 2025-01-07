#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "ECEPE.h"
#include "EF.h"

#include <glib.h>

#define ON  1
#define OFF 0

void USAGE(char *progname);

#define numatommasstype 5

char *atommasstypelist[2000]= {
  "C",
  "N",
  "H",
  "O",
  "S"
};

char *atommasslist[2000]={
  "12.01     ",
  "14.0      ",
  "1.0       ",
  "16.0      ",
  "32.060    "
};

int main(int argc, char *argv[]) {
  int i,j;
  int numatom;
  char *mass;
  char *named[4];
  char *name_atom_list;

  gchar *atomnameaskey;
  gchar *atommass;
  gint *atonum;
  GHashTable *masslist;

  struct ECEPE_parms ECEPE_p;
  struct pnb *nb_p;

  char *progname;
  char *preofilename,*bd8filename;
  char *outputfilename;
  FILE *outputfile;

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

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
  outputfilename = *++argv;

  read_ECEPE_parm_wobd8(preofilename/*,bd8filename*/,&ECEPE_p/*,&nb_p*/);
  numatom=ECEPE_p.NUMATM;

  masslist=g_hash_table_new_full(g_str_hash,g_str_equal,g_free,g_free);

  for (i=0;i<numatommasstype;++i) {
    atomnameaskey=g_strdup(atommasstypelist[i]);
    atommass=g_new(gchar,1);
    atommass=atommasslist[i];
    g_hash_table_insert(masslist,atomnameaskey,atommass);
    if ((mass=g_hash_table_lookup(masslist,atommasstypelist[i]))==NULL) {
      printf("There is no key %s in hashtable",atommasstypelist[i]);
      exit(1);
    } 
  }

  name_atom_list=(char *)gcemalloc(sizeof(char)*numatom*4);
  for (i=0;i<numatom;++i)
    for (j=0;j<4;++j)
      name_atom_list[(ECEPE_p.atom[i].katom-1)*4+j]=ECEPE_p.atom[i].name_atom[j];

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numatom;++i) {    
    for (j=0;j<4;++j) named[j]=name_atom_list[i*4+j];
    if ((mass=g_hash_table_lookup(masslist,named))==NULL) {
      printf("There is no key %s in hashtable",named);
      exit(1);
    } 
    fprintf(outputfile,"%s ",mass);
    if ((i+1)%10==0) fprintf(outputfile,"\n");
  }
  fclose(outputfile);

}

void USAGE(char *progname) {
  printf("-h -- help\n");
  printf("USAGE: %s profilename outputfilename\n", progname);
}
