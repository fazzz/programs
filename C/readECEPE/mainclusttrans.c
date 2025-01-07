#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "ECEPE.h"

#include "EF.h"

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,nb;

  int numclust;
  int numnewclust;
  int sum;


  int *origin_atom_a,*terminal;
  int *num_atom_clust,*num_branch,**terminal_atom_a,*hingmat;
  int *nNumClutOfParent,**nNumClutOfChild,*IndexOfABICycle,MAP;

  char *progname;
  char *clustfilename,*outputfilename;
  FILE *clustfile,*outputfile;

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
  numnewclust = atoi(*argv);
  clustfilename  = *++argv;
  outputfilename = *++argv;

  clustfile=efopen(clustfilename,"r");
  
  fscanf(clustfile,"%d",&numclust);

  origin_atom_a=(int *)gcemalloc(sizeof(int)*numclust);
  terminal=(int *)gcemalloc(sizeof(int)*numclust);
  num_atom_clust=(int *)gcemalloc(sizeof(int)*numclust);
  num_branch=(int *)gcemalloc(sizeof(int)*numclust);
  terminal_atom_a=(int **)gcemalloc(sizeof(int *)*numclust);
  hingmat=(int *)gcemalloc(sizeof(int)*numclust);
  for (i=0;i<numclust;++i) terminal_atom_a[i]=(int *)gcemalloc(sizeof(int)*2);
  nNumClutOfParent=(int *)gcemalloc(sizeof(int)*numclust);
  nNumClutOfChild=(int **)gcemalloc(sizeof(int *)*numclust);
  for (i=0;i<numclust;++i) nNumClutOfChild[i]=(int *)gcemalloc(sizeof(int)*2);
  IndexOfABICycle=(int *)gcemalloc(sizeof(int)*numclust);

  for(k=0;k<numclust;++k)
    fscanf(clustfile,"%d",&origin_atom_a[k]);

  for(k=0;k<numclust;++k)
    fscanf(clustfile,"%d",&terminal[k]);

  for(k=0;k<numclust;++k)
    fscanf(clustfile,"%d",&num_atom_clust[k]);

  for(k=0;k<numclust;++k)
    fscanf(clustfile,"%d",&num_branch[k]);

  for(k=0;k<numclust;++k)
    fscanf(clustfile,"%d",&hingmat[k]);

  for(k=0;k<numclust;++k)
    for(nb=0;nb<num_branch[k];++nb)
      fscanf(clustfile,"%d",&terminal_atom_a[k][nb]);


  for (i=0;i<numclust;++i)
    fscanf(clustfile, "%d", &nNumClutOfParent[i]);

  for (i=0;i<numclust;++i)
    for(nb=0;nb<num_branch[i];++nb)
      fscanf(clustfile, "%d",&nNumClutOfChild[i][nb]);

  for (i=0;i<numclust;++i)
    fscanf(clustfile, "%d", &IndexOfABICycle[i]);

  fscanf(clustfile, "%d",&MAP);
  fclose(clustfile);

  outputfile=efopen(outputfilename,"w");

  fprintf(outputfile,"%d\n",numnewclust);
  
  for(k=0;k<numnewclust;++k) {
    fprintf(outputfile,"%d ",origin_atom_a[k]);
    if ((k+1)%10==0) fprintf(outputfile,"\n");
  }
  for(k=0;k<numnewclust-1;++k) {
    fprintf(outputfile,"%d ",terminal[k]);
    if ((k+1)%10==0) fprintf(outputfile,"\n");
  }
  fprintf(outputfile,"0\n");
  for(k=0;k<numnewclust-1;++k) {
    fprintf(outputfile,"%d ",num_atom_clust[k]);
    if ((k+1)%10==0) fprintf(outputfile,"\n");
  }
  sum=0;
  for (k=numnewclust-1;k<numclust;++k) {
    sum+=num_atom_clust[k];
  }
  fprintf(outputfile,"%d\n",sum);

  for(k=0;k<numnewclust;++k) {
    fprintf(outputfile,"%d ",num_branch[k]);
    if ((k+1)%10==0) fprintf(outputfile,"\n");
  }
  for(k=0;k<numnewclust;++k) {
    fprintf(outputfile,"%d ",hingmat[k]);
    if ((k+1)%10==0) fprintf(outputfile,"\n");
  }
  for(k=0;k<numnewclust;++k) {
    for(nb=0;nb<num_branch[k];++nb) {
      fprintf(outputfile,"%d ",terminal_atom_a[k][nb]);
    }
    if ((k+1)%10==0) fprintf(outputfile,"\n");
  }
  for (i=0;i<numnewclust;++i) {
    fprintf(outputfile, "%d ", nNumClutOfParent[i]);
    if ((k+1)%10==0) fprintf(outputfile,"\n");
  }

  for (i=0;i<numnewclust;++i) {
    for(nb=0;nb<num_branch[i];++nb) {
      fprintf(outputfile, "%d ",nNumClutOfChild[i][nb]);
    }
    if ((i+1)%10==0) fprintf(outputfile,"\n");
  }
  for (i=0;i<numnewclust;++i) {
    fprintf(outputfile, "%d ", IndexOfABICycle[i]);
    if ((i+1)%10==0) fprintf(outputfile,"\n");
  }

  fclose(outputfile);


  return 0;
}

void USAGE(char *progname) {
  printf("-h -- help\n");
  printf("USAGE: %s profilename bd8filename \n", progname);
}



