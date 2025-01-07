
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SSL.h"
#include "EF.h"

int scnadata(FILE *inputfile,int nums,int numv,double *data,int flag);
int scnaslidedata(FILE *inputfile,int nums, int numslide,int numv,double *data,int flag);
int setmat(double *Lambdap,double *Lambda,int numv);
int output(char *outputfilenamebase,char *outputfilenamebase2,double *Lambda, double *Sigma,double *COVM,int numv,double rou,int numite,int outflag);
int outputKLdiv(FILE *outputfile,double *KLdiv,int numv, int numite);
int outputKLdivtopten(FILE *outputfile,int topten[10], double vtopten[10], int num, int numite, double sumKLdiv);

int main(int argc, char *argv[]) {
  int i,j,k;
  int nums,numv,numite,numstotal,numslide,topten[10],inflag,num,outflag,outflag2;
  double *data,*datanorm,*COVM,*Lambda,*Sigma,*Lambdap,*Sigmap,*Lambdaini,*Sigmaini,rou,*KLdiv,vtopten[10],sumKLdiv;

  char *inputfilename,*inputfilename2;
  char *outputfilenamebase,*outputfilenamebase2,*outputfilename3,*outputfilename4,*outputfilename5,*outputfilename6;
  FILE *inputfile,*inputfile2;
  FILE *outputfile3,*outputfile4,*outputfile5,*outputfile6;

  if (argc < 3) {
    printf("USAGE: ./ssl.exe inflag outflag inputfilename(data) inputfilename2(cond) outputfilenamebase outputfilebase2 outputfile3 outputfile4\n");
    exit(1);
  }
  inflag=(*++argv)[0];
  if (inflag != 'c' && inflag != 'v') {
    printf("inflag error: must be c or v");
    exit(1);
  }
  outflag2=(*++argv)[0];
  if (outflag2 != 'e' && outflag2 != 'b') {
    printf("outflag error: must be e or b");
    exit(1);
  }

  inputfilename   = *++argv;
  inputfilename2  = *++argv;
  outputfilenamebase  = *++argv;
  outputfilenamebase2 = *++argv;
  outputfilename3 = *++argv;
  outputfilename4 = *++argv;
  if (outflag2 == 'b') {
    outputfilename5 = *++argv;
    outputfilename6 = *++argv;
  }

  inputfile2=efopen(inputfilename2,"r");
  fscanf(inputfile2,"%d",&numstotal);
  fscanf(inputfile2,"%d",&nums);
  fscanf(inputfile2,"%d",&numslide);
  fscanf(inputfile2,"%d",&numv);
  fscanf(inputfile2,"%d",&num);
  fscanf(inputfile2,"%lf",&rou);
  fscanf(inputfile2,"%d",&outflag);
  fclose(inputfile2);

  inputfile   =  efopen(inputfilename,"r");
  if (outflag2 == 'b') {
    outputfile5 = efopen(outputfilename5,"w");
    outputfile6 = efopen(outputfilename6,"w");
  }
  outputfile3 = efopen(outputfilename3,"w");
  outputfile4 = efopen(outputfilename4,"w");

  data      = (double *)ecalloc(numv*nums,sizeof(double));
  datanorm  = (double *)ecalloc(numv*nums,sizeof(double));
  COVM      = (double *)ecalloc(numv*numv,sizeof(double));
  Lambda    = (double *)ecalloc(numv*numv,sizeof(double));
  Sigma     = (double *)ecalloc(numv*numv,sizeof(double));
  Lambdap   = (double *)ecalloc(numv*numv,sizeof(double));
  Sigmap    = (double *)ecalloc(numv*numv,sizeof(double));
  Lambdaini = (double *)ecalloc(numv*numv,sizeof(double));
  Sigmaini  = (double *)ecalloc(numv*numv,sizeof(double));
  KLdiv     = (double *)ecalloc(numv,sizeof(double));

  scnadata(inputfile,nums,numv,data,inflag);
  ssl_normalize(data,nums,numv,datanorm);
  ssl_covm(datanorm,nums,numv,COVM);
  ssl_gralassomain(COVM, rou, numv, Lambdaini, Sigmaini);
  output(outputfilenamebase,outputfilenamebase2,Lambdaini,Sigmaini,COVM,numv,rou,0,outflag);

  numite=(numstotal-nums)/numslide;
  for (i=0;i<numite;++i) {
    if (outflag2=='b') {
      setmat(Lambdap,Lambda,numv);
      setmat(Sigmap,Sigma,numv);
    }
    scnaslidedata(inputfile,nums,numslide,numv,data,inflag);
    ssl_normalize(data,nums,numv,datanorm);
    ssl_covm(datanorm,nums,numv,COVM);
    ssl_gralassomain(COVM, rou, numv, Lambda, Sigma);
    output(outputfilenamebase,outputfilenamebase2,Lambda,Sigma,COVM,numv,rou,i+1,outflag);
    if (outflag2=='b') {
      sumKLdiv=ssl_KLdiv(Lambda,Sigma,Lambdap,Sigmap,numv,KLdiv,topten,vtopten);
      outputKLdivtopten(outputfile5,topten,vtopten,num,i+1,sumKLdiv);
      outputKLdiv(outputfile6,KLdiv,numv,i+1);
    }
    sumKLdiv=ssl_KLdiv(Lambda,Sigma,Lambdaini,Sigmaini,numv,KLdiv,topten,vtopten);
    outputKLdivtopten(outputfile3,topten,vtopten,num,i+1,sumKLdiv);
    outputKLdiv(outputfile4,KLdiv,numv,i+1);
  }

  fclose(inputfile);
  fclose(outputfile3);
  fclose(outputfile4);
  if (outflag2 == 'b') {
    fclose(outputfile5);
    fclose(outputfile6);
  }
  free(data);
  free(COVM);
  free(Lambda);    free(Sigma);
  free(Lambdap);   free(Sigmap);
  free(Lambdaini); free(Sigmaini);

  return 0;
}

int scnadata(FILE *inputfile,int nums,int numv,double *data,int flag){
  int i,j;
  double f,d;

  for (i=0;i<nums;++i) {
    if (flag=='c')
      fscanf(inputfile,"%d",&d);
    for (j=0;j<numv;++j) {
      fscanf(inputfile,"%lf",&f);
      data[i*numv+j]=f;
    }
  }

  return 0;
}

int scnaslidedata(FILE *inputfile,int nums, int numslide,int numv,double *data,int flag) {
  int i,j,d;
  double f;
  double *datadummy;

  datadummy=(double *)ecalloc(numv*(nums-numslide),sizeof(double));

  for (i=numslide;i<nums;++i)
    for (j=0;j<numv;++j)
      datadummy[(i-numslide)*numv+j]=data[i*numv+j];

  for (i=0;i<nums-numslide;++i)
    for (j=0;j<numv;++j)
      data[i*numv+j]=datadummy[i*numv+j];

  for (i=0;i<numslide;++i) {
    if (flag=='c')
      fscanf(inputfile,"%d",&d);
    for (j=0;j<numv;++j) {
      fscanf(inputfile,"%lf",&f);
      data[(i+nums-numslide)*numv+j]=f;
    }
  }

  return 0;
}

 int setmat(double *Lambdap,double *Lambda,int numv) {
  int i,j;

  for (i=0;i<numv;++i)
    for (j=0;j<numv;++j)
      Lambdap[i*numv+j]=Lambda[i*numv+j];

}

int output(char *outputfilenamebase,char *outputfilenamebase2,double *Lambda, double *Sigma,double *COVM,int numv,double rou,int numite,int outflag) {
  int i,j,numspars;
  double sparsity;
  char outputfilename[100],outputfilename2[100];
  FILE *outputfile,*outputfile2;

  sprintf(outputfilename,"%s_%d.txt",outputfilenamebase,numite);
  sprintf(outputfilename2,"%s_%d.txt",outputfilenamebase2,numite);

  outputfile  =  efopen(outputfilename,"w");
  outputfile2 = efopen(outputfilename2,"w");

  if (outflag==0) {
    numspars=0;
    for (i=0;i<numv;++i) {
      for (j=0;j<numv;++j) {
	if (j<i) {
	  fprintf(outputfile,"0 ");
	}
	else {
	  if (Lambda[i*numv+j]!=0)
	    fprintf(outputfile,"1 ");
	  else {
	    fprintf(outputfile,"0 ");
	    ++numspars;
	  }
	}
      }
      fprintf(outputfile,"\n");    
    }
    sparsity = 2.0*numspars/(numv*(numv-1));
    fprintf(outputfile,"%lf \n",sparsity);
    fprintf(outputfile,"\n");
    fclose(outputfile);
  }
  else {
    for (i=0;i<numv;++i) {
      for (j=0;j<numv;++j) {
	if (j<i)
	  fprintf(outputfile,"0.0 ");
	else
	  fprintf(outputfile,"%lf ",Lambda[i*numv+j]);
      }
      fprintf(outputfile,"\n");    
    }
    fprintf(outputfile,"\n");
    for (i=0;i<numv;++i) {
      for (j=0;j<numv;++j) {
	if (j<i)
	  fprintf(outputfile,"0.0 ");
	else
	  fprintf(outputfile,"%lf ",Sigma[i*numv+j]);
      }
      fprintf(outputfile,"\n");    
    }
    fprintf(outputfile,"\n test");
    fclose(outputfile);
  }

  for (i=0;i<numv;++i) {
    for (j=0;j<numv;++j) {
      if (j<i) {
	fprintf(outputfile2,"0 ");
      }
      else {
	fprintf(outputfile2,"%e ",COVM[i*numv+j]);
      }
    }
    fprintf(outputfile2,"\n");    
  }
  fprintf(outputfile2,"\n");
  fclose(outputfile2);

}

int outputKLdiv(FILE *outputfile,double *KLdiv,int numv, int numite){
  int i;

  if (numite==0) {
    fprintf(outputfile,"step ");
    for (i=0;i<numv;++i) {
      fprintf(outputfile,"%d ",i);
    }
  }
  fprintf(outputfile,"%6d ",numite);
  for (i=0;i<numv;++i) {
    fprintf(outputfile,"%e ",KLdiv[i]);
  }
  fprintf(outputfile,"\n");    

}

int outputKLdivtopten(FILE *outputfile,int topten[10], double vtopten[10], int num, int numite, double sumKLdiv){
  int i;

  fprintf(outputfile,"%d to %d \n",numite-1,numite);
  for (i=0;i<num;++i)
    fprintf(outputfile,"%4d ",topten[i]);
  fprintf(outputfile,"\n     ");
  for (i=0;i<num;++i)
    fprintf(outputfile,"%e ",vtopten[i]);
  fprintf(outputfile,"\n");
  fprintf(outputfile,"%4.2lf \n",sumKLdiv);

}

