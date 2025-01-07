
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "PCA.h"
#include "IO.h"
#include "PT.h"
#include "EF.h"

int main(int argc, char *argv[]) {
  int i,j;
  int numatom,numstep,numtype,typeflag;
  char atomtype[20][4];
  double *traj,*cov,*eigenvalue;
  
  char *inputfilename,*inputfilename2,*inputfilename3,*outputfilename,*outputfilename2;
  FILE *inputfile,*inputfile2,*inputfile3,*outputfile,*outputfile2;
  
  if (argc < 6) {
    printf("USAGE: %s inputfilename(trj) inputfilename2(cond) inputfilename3(parm)  outputfilename(trj) outputfilename2(cov)\n",argv[0]);
    exit(1);
  }
  inputfilename =  *++argv;
  inputfilename2 = *++argv;
  inputfilename3 = *++argv;
  outputfilename = *++argv;
  outputfilename2= *++argv;
  
  inputfile2=efopen(inputfilename2,"r");
  fscanf(inputfile2,"%d",&numstep);
  fscanf(inputfile2,"%d",&typeflag);
  if (typeflag==1) {
    fscanf(inputfile2,"%d",&numtype);
    for (i=0;i<numtype;++i)
      fscanf(inputfile2,"%4s",&atomtype[i]);
  }
  fclose(inputfile2);

  inputfile3=efopen(inputfilename3,"r");
  readParmtop(inputfile3);
  fclose(inputfile3);
  if (typeflag==1)
    numatom=io_count_atomtype(atomtype,numtype);
  else
    numatom=AP.NATOM;
  traj       = (double *)egccalloc(sizeof(double),numatom*3*numstep);
  cov        = (double *)egccalloc(sizeof(double),numatom*3*numatom*3);
  eigenvalue = (double *)egccalloc(sizeof(double),numatom*3);

  inputfile=efopen(inputfilename,"r");
  io_scan_atomtype_traj_aw(inputfile,numstep,traj,typeflag,atomtype,numtype,numatom);
  fclose(inputfile);

  pca_norm(traj,numstep,numatom);
  pca_covm(traj,numstep,numatom,cov);
  pca_diag(cov,eigenvalue,numatom);
  pca_proj(traj,cov,numstep,numatom);

  outputfile=efopen(outputfilename,"w");
  io_outputtimeseries_f(outputfile,numstep,numatom*3,traj);
  fclose(outputfile);

  outputfile2=efopen(outputfilename2,"w");
  fprintf(outputfile,"n eigenvalue\n");
  io_outputdata_f(outputfile2,numatom*3,eigenvalue);
  fclose(outputfile2);

  return 0;
}








/***************************************************************************************************************/
/* void scantraj(FILE *inputfile, /\*double *traj,*\/int time, int numatom);				       */
/* int scanmass(FILE *inputfile,int numatom,double mass[MAXNUMATOM]);					       */
/* 													       */
/* int main(int argc, char *argv[]) {									       */
/*   int i,j,k;												       */
/*   int minnumoption=4;										       */
/*   char *inputfilename0,*inputfilename1,*inputfilename2, *outputfilename;				       */
/*   char *USAGE="USAGE: ./calc_PCA condfile trajfile parmfile outputfile\n";				       */
/*   char *ERROR="does not exist in this directry\n";							       */
/*   FILE *inputfile0,*inputfile1,*inputfile2,*outputfile;						       */
/* 													       */
/*   int numatom,time;											       */
/*   double mass[MAXNUMATOM];										       */
/*   //  double *traj/\*numatom*time*3*\/;								       */
/*   //  double *cov/\*numatom*time*3*numatom*time*3*\/;						       */
/*   //  double *crd_ave;/\*numatom*3*\/								       */
/*   double covR[3*MAXNUMATOM*3*MAXNUMATOM];								       */
/*   double w[3*MAXNUMATOM];										       */
/*   double sum,cont;											       */
/* 													       */
/* 													       */
/*   if (argc < minnumoption){										       */
/*     printf("%s",USAGE);										       */
/*     exit(1);												       */
/*   }													       */
/* 													       */
/*   inputfilename0 = *++argv;										       */
/*   inputfilename1 = *++argv;										       */
/*   inputfilename2 = *++argv;										       */
/*   outputfilename = *++argv;										       */
/* 													       */
/*   if ((inputfile0=fopen(inputfilename0,"r"))==NULL){							       */
/*     printf("%s %s",inputfilename0, ERROR);								       */
/*     exit(1);												       */
/*   }													       */
/*   fscanf(inputfile0,"%d",&time);									       */
/*   fscanf(inputfile0,"%d",&numatom);									       */
/*   fclose(inputfile0);										       */
/* 													       */
/*   if ((inputfile1=fopen(inputfilename1,"r"))==NULL){							       */
/*     printf("%s %s",inputfilename1, ERROR);								       */
/*     exit(1);												       */
/*   }													       */
/*   scantraj(inputfile1,/\*traj,*\/time,numatom);							       */
/*   fclose(inputfile1);										       */
/* 													       */
/*   if ((inputfile2=fopen(inputfilename2,"r"))==NULL){							       */
/*     printf("%s %s",inputfilename2, ERROR);								       */
/*     exit(1);												       */
/*   }													       */
/*   readParmtop(inputfile2);										       */
/*   for (i=0;i<numatom;++i){										       */
/*     mass[i]=AP.AMASS[i];										       */
/*   }													       */
/*   //  scanmass(inputfile2,numatom,mass);								       */
/*   fclose(inputfile2);										       */
/* 													       */
/* 													       */
/*   //  cov=malloc(sizeof(double)*numatom*3*numatom*3);						       */
/*   crd_ave=malloc(sizeof(double)*numatom*3);								       */
/*   calc_PCA(traj,/\*cov,*\/crd_ave,time,numatom,mass,covR,w);						       */
/*   proj_PCA(traj,/\*traj_trans,*\/crd_ave,time,numatom,mass,covR,w);					       */
/* 													       */
/*   if ((outputfile=fopen("eigen_value","w"))==NULL){							       */
/*     printf("eigen_value %s",ERROR);									       */
/*     exit(1);												       */
/*   }													       */
/*   if ((outputfile=fopen("axis","w"))==NULL){								       */
/*     printf("eigen_value %s",ERROR);									       */
/*     exit(1);												       */
/*   }													       */
/*   sum=0.0;												       */
/*   for (i=0;i<numatom*3;++i)										       */
/*     sum+=w[i];											       */
/* 													       */
/*   cont=0.0;												       */
/*   for (i=numatom*3-1;i>=0;--i){									       */
/*     cont+=w[i];											       */
/*     fprintf(outputfile,"%4d -th  %12.8lf %12.8lf \n",numatom*3-i, w[i], cont/sum);			       */
/*   }													       */
/* 													       */
/*   for (j=0;j<numatom*3;++j){										       */
/*     fprintf(outputfile,"%8.3lf",crd_ave[j]);								       */
/*   }													       */
/* 													       */
/*   for (j=0;j<numatom*3;++j){										       */
/*     fprintf(outputfile,"%8.3lf",covR[(numatom)*3-1+j*numatom*3]);					       */
/*   }													       */
/*  fprintf(outputfile,"\n");										       */
/*  for (j=0;j<numatom*3;++j){										       */
/*    fprintf(outputfile,"%8.3lf",covR[(numatom)*3-2+j*numatom*3]);					       */
/*   }													       */
/*  fprintf(outputfile,"\n");										       */
/*   fclose(outputfile);										       */
/* 													       */
/*   if ((outputfile=fopen(outputfilename,"w"))==NULL){							       */
/*     printf("%s %s",outputfilename,ERROR);								       */
/*     exit(1);												       */
/*   }													       */
/* 													       */
/* 													       */
/*   if ((outputfile=fopen(outputfilename,"w"))==NULL){							       */
/*     printf("%s %s",outputfilename,ERROR);								       */
/*     exit(1);												       */
/*   }													       */
/* 													       */
/*   for (i=0;i<time;++i) {										       */
/*     fprintf(outputfile,"%lf %lf\n",traj_trans[i*2/\*numatom*3*\/],traj_trans[i*2/\*numatom*3*\/+1]);	       */
/*   }													       */
/* 													       */
/*   fclose(outputfile);										       */
/* 													       */
/*   if ((outputfile=fopen("traj.txt","w"))==NULL){							       */
/*     printf("%s %s",outputfilename,ERROR);								       */
/*     exit(1);												       */
/*   }													       */
/* 													       */
/*   for (i=0;i<time;++i){										       */
/*     for (j=0;j<numatom;++j){										       */
/*       for (k=0;k<3;++k){										       */
/* 	fprintf(outputfile,"%lf ",traj[i*numatom*3+j*3+k]);						       */
/*       }												       */
/*       fprintf(outputfile,"\n");									       */
/*     }												       */
/*   }													       */
/* 													       */
/* 													       */
/*   return 0;												       */
/* 													       */
/* }													       */
/* 													       */
/* void scantraj(FILE *inputfile, /\*double *traj,*\/int time, int numatom)				       */
/* {													       */
/*   int i,j,k;												       */
/* 													       */
/*   traj=malloc(sizeof(double)*time*numatom*3);							       */
/* 													       */
/*   for (i=0;i<time;++i)										       */
/*     for (j=0;j<numatom;++j)										       */
/*       for (k=0;k<3;++k)										       */
/* 	fscanf(inputfile,"%lf", &traj[i*numatom*3+j*3+k]);						       */
/* 													       */
/* }													       */
/* 													       */
/* int scanmass(FILE *inputfile,int numatom,double mass[MAXNUMATOM])					       */
/* {													       */
/*   int i;												       */
/*   for (i=0;i<numatom;++i)										       */
/*     fscanf(inputfile,"%lf",&mass[i]);								       */
/* 													       */
/* }													       */
/***************************************************************************************************************/
