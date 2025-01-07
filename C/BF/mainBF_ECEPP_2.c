
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "QUA.h"
#include "bestfit.h"
#include "PT.h"
#include "EF.h"
#include "IO.h"

#include "netcdf_mine.h"

#include "f2c.h"
#include "clapack.h"

#define MAXNUNITERATION 20

void USAGE(void);

int main(int argc, char *argv[]) {
  int i,j,k,s,numatom,numstep,numite;
  int flagv='c',flagamber='a',flago='x',flagc='a';
  int MAX_numiteration=MAXNUNITERATION;
  double epsilon=0.0001;
  double *mass,*coord_ref,*rmsd_trj;
  double rmsd_ave=0.0,rmsd_ave_old=0.0;
  char *pname;
  char *inputfilename,*inputfilename2, *inputfilename3,outputfilename[100],outputfilename2[100];
  char *outputfilenamebase,outputfilenamermsd[100];

  FILE *inputfile, *inputfile2,*inputfile3;
  FILE *outputfile,*outputfile2,*outputfilermsd,*log;

  char *line;
  size_t len=0;
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  while((c=getopt(argc,argv,"vache:n:"))!=-1) {
    switch(c) {
    case 'v':
      flagv='b';
      break;
    case 'c':
      flagc='c';
      break;
    case 'a':
      flagamber='a';
      break;
    case 'e':
      epsilon=atof(optarg);
      break;
    case 'n':
      //      numite=atoi(optarg);
      MAX_numiteration=atoi(optarg);
      break;
    case 'h':
      USAGE();
      printf("USAGE:%s numstep numatom input(traj) (input(velo)) input(ParmTop) outputbasename\n",argv[0]);
      break;
    default:
      USAGE();
      printf("USAGE:%s numstep numatom input(traj) (input(velo)) input(ParmTop) outputbasename\n",argv[0]);
      exit(1);
    }
  }

  pname=argv[0];
  argc-=optind;
  argv+=optind;

  if (flagv == 'b') {
    if (argc < 6) {
      printf("USAGE:%s numstep input(traj) input(velo) input(ParmTop) outputbasename\n",pname);
      exit(1);
    }
  }
  else {
    if (argc < 5) {
      printf("USAGE:%s numstep input(traj) input(ParmTop) outputbasename\n",pname);
      exit(1);
    }
  }
  numstep = atoi(*argv);
  numatom = atoi(*++argv);
  inputfilename      =  *++argv;
  if (flagv=='b')
    inputfilename2   =  *++argv;
  inputfilename3     =  *++argv;
  outputfilenamebase =  *++argv;

  inputfile3=efopen(inputfilename3,"r");
  readParmtop(inputfile3);
  fclose(inputfile3);
  if (flagc=='c') {
    //    numatom = AP.NRES;
    mass=(double *)gcemalloc(sizeof(double)*numatom);
    io_scancamass(inputfile,numatom,mass);
  }
  else {
    //    numatom = AP.NATOM;
    mass=(double *)gcemalloc(sizeof(double)*numatom);
    for (i=0;i<numatom;++i)
      mass[i] = AP.AMASS[i];
  }

  coord_ref=(double *)gcemalloc(sizeof(double)*numatom*3);
  rmsd_trj =(double *)gcemalloc(sizeof(double)*numstep);


  for (i=0;i<MAX_numiteration;++i) {
    if (i>0) {
      sprintf(inputfilename,"%s_cyc=%d_bf.trj",outputfilenamebase,i);
      if (flagv == 'b')
	sprintf(inputfilename2,"%s_cyc=%d_bf.vel",outputfilenamebase,i);
      /**************************************************/
      /* inputfile=efopen(inputfilename,"r");	        */
      /* if (flagv == 'b')			        */
      /* 	inputfile2=fopen(inputfilename2,"r");   */
      /**************************************************/
    }

    inputfile=efopen(inputfilename,"r");
    if (flagamber=='a' && i==0)
      getline(&line,&len,inputfile );
    ave_coord(coord_ref,inputfile,numatom,numstep,flagc);
    fclose(inputfile);

    /*******************************************/
    /* if (i>0) 			       */
    /*   inputfile=efopen(outputfilename,"r"); */
    /* else 				       */
    /*   inputfile=efopen(inputfilename,"r");  */
    /*******************************************/
    if (flago=='o')
      sprintf(outputfilename,"%s_bf.trj",outputfilenamebase);
    else
      sprintf(outputfilename,"%s_cyc=%d_bf.trj",outputfilenamebase,i+1);
    /******************************************/
    /* outputfile=efopen(outputfilename,"w"); */
    /******************************************/
    if (flagv == 'b'){
      if (flago=='o')
	sprintf(outputfilename2,"%s_bf.vel",outputfilenamebase);
      else
	sprintf(outputfilename2,"%s_cyc=%d_bf.vel",outputfilenamebase,i+1);
      /********************************************/
      /* outputfile2=efopen(outputfilename2,"w"); */
      /********************************************/
    }

    if (i>0)
      flagamber='x';
    bf_trajectry(numatom,numstep,mass,coord_ref,rmsd_trj,inputfilename,inputfilename2,outputfilename,outputfilename2,flagv,flagamber,'o',flagc);
    if (flago=='o') {
      sprintf(outputfilenamermsd,"%s_rmsd.txt",outputfilenamebase);
      outputfilermsd=efopen(outputfilenamermsd,"w");
      fprintf(outputfilermsd,"step    rmsd\n");
      for (j=0;j<numstep;++j)
	fprintf(outputfilermsd,"%d %12.8lf\n",j+1,rmsd_trj[j]);
      fclose(outputfilermsd);
      break;
    }
    /**************************/
    /* if (flagv=='b')	      */
    /*   fclose(inputfile2);  */
    /* fclose(outputfile);    */
    /* if (flagv=='b')	      */
    /*   fclose(outputfile2); */
    /**************************/

    flago='x';
    rmsd_ave=0.0;
    for (j=0;j<numstep;++j)
      rmsd_ave=(j*rmsd_ave+rmsd_trj[j])/(j+1);
    printf("cyc: %d  rmsd: %10.6lf\n",i+1,rmsd_ave);
    if (abs(rmsd_ave-rmsd_ave_old) < epsilon && i > 0)
      flago='o';
    rmsd_ave_old=rmsd_ave;
  }
   return 0;
}

void USAGE(void) {
  printf("-v[BF for velocity] -a[Amber mode] -c [alpha carbon] -h [help] -e [epsilon criteria of delta rmsd] -n [num of iteration (<10)]");
}
