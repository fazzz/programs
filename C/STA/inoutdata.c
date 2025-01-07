
#include <stdio.h>
#include <stdlib.h>

#include "EF.h"
#include "PT.h"

int io_scanconf(FILE *inputfile,int numatom,double *data,int flag){
  int i,j,d;
  double f;
  
  for (i=0;i<numatom;++i) {
    if (flag=='c')
      fscanf(inputfile,"%d",&d);
    for (j=0;j<3;++j)
      fscanf(inputfile,"%lf",&data[i*3+j]);
  }

  return 0;
}

int io_scan_aw(FILE *inputfile,int numstep,double *data) {
  int i,j,k;
  int numatom;
  double d;

  numatom=AP.NATOM;
  for (i=0;i<numstep;++i)  {
    for (j=0;j<numatom;++j)  {
      for (k=0;k<3;++k) {
	fscanf(inputfile,"%lf",&d);
	data[i*numatom*3+j*3+k] = d*sqrt(AP.AMASS[j]);
      }
    }
  }

  return 1;
}

int io_count_atomtype(char atomtype[20][4],int numtype) {
  int i,j;
  int numatom,numatomtype=0;

  numatom=AP.NATOM;
  for (i=0;i<numatom;++i)
    for (j=0;j<numtype;++j)
      if (strncmp(AP.IGRAPH[i],atomtype[j],3) == 0)
	++numatomtype;

  if (numatomtype==0) {
    printf("error:There are no %s\n",atomtype);
    exit(1);
  }

  return numatomtype;
}

int io_scan_atomtype_traj_aw(FILE *inputfile,int numstep,double *data, int typeflag,char atomtype[20][4], int numtype,int numatomtype) {
  int i,j,k,l;
  int numatom;
  double d;

  numatom=AP.NATOM;

  for (i=0;i<numstep;++i)  {
    for (j=0;j<numatom;++j)  {
      for (k=0;k<3;++k) {
	fscanf(inputfile,"%lf",&d);
	if (typeflag==1) {
	  for (l=0;l<numtype;++l)
	    if (strncmp(AP.IGRAPH[j],atomtype[l],3) == 0) {
	      data[i*numatomtype*3+j*3+k] = d*sqrt(AP.AMASS[j]);
	  }
	}
	else {
	    data[i*numatomtype*3+j*3+k] = d*sqrt(AP.AMASS[j]);
	}
      }
    }
  }

  return numatomtype;
}

int io_scan_data(FILE *inputfile,int numdimension,double *data,int flag){
  int i,d;

  if (flag=='c')
    fscanf(inputfile,"%d",&d);
  for (i=0;i<numdimension;++i) {
    fscanf(inputfile,"%lf",&data[i]);
  }

  return 0;
}

int io_scantimeseries(FILE *inputfile,int numstep,int numdimension,double *data,int flag){
  int i,j,d;
  
  for (i=0;i<numstep;++i) {
    if (flag=='c')
      fscanf(inputfile,"%d",&d);
    for (j=0;j<numdimension;++j) {
      fscanf(inputfile,"%lf",&data[i*numdimension+j]);
    }
  }
  
  return 0;
}

int io_scanslidetimeseries(FILE *inputfile,int numstep, int numslide,int numdimension,double *data,int flag) {
  int i,j,d;
  double f;
  double *datadummy;
  
  datadummy=(double *)ecalloc(numdimension*(numstep-numslide),sizeof(double));
  
  for (i=numslide;i<numstep;++i)
    for (j=0;j<numdimension;++j)
      datadummy[(i-numslide)*numdimension+j]=data[i*numdimension+j];
  
  for (i=0;i<numstep-numslide;++i)
    for (j=0;j<numstep;++j)
      data[i*numstep+j]=datadummy[i*numstep+j];

  for (i=0;i<numslide;++i) {
    if (flag=='c')
      fscanf(inputfile,"%d",&d);
    for (j=0;j<numstep;++j) {
      fscanf(inputfile,"%lf",&f);
      data[(i+numstep-numslide)*numdimension+j]=f;
    }
  }
  
  return 0;
}

int io_scanmatrix(FILE *inputfile,int numrow,int numcolum,double *data){
  int i,j;
  
  for (i=0;i<numrow;++i) {
    for (j=0;j<numcolum;++j) {
      fscanf(inputfile,"%lf",&data[i*numcolum+j]);
    }
  }
    
  return 0;
}
  
int io_outputconf(FILE *outputfile,int numatom,double *data,int flag){
  int i,j,d;
  double f;
  
  for (i=0;i<numatom;++i) {
    if (flag=='c')
      fprintf(outputfile,"%d ",i);
    for (j=0;j<3;++j)
      fprintf(outputfile,"%10.6lf",data[i*3+j]);
    fprintf(outputfile,"\n");
  }
  fprintf(outputfile,"\n");

  return 0;
}

int io_outputdata(FILE *outputfile, int numdimension, double *data) {
  int i;
  
  for (i=0;i<numdimension;++i)
    fprintf(outputfile,"%e ",data[i]);
  fprintf(outputfile,"\n");

  return 0;
}

int io_outputdata_f(FILE *outputfile, int numdimension, double *data) {
  int i;
  
  for (i=0;i<numdimension;++i)
    fprintf(outputfile,"%d %e \n",i+1,data[i]);

  return 0;
}


int io_outputtimeseries(FILE *outputfile,int numstep,int numdimension,double *data){
  int i,j,d;
  
  for (i=0;i<numstep;++i) {
    for (j=0;j<numdimension;++j) {
      fprintf(outputfile,"%e",data[i*numdimension+j]);
    }
    fprintf(outputfile,"\n");
  }
  
  return 0;
}

int io_outputtimeseries_f(FILE *outputfile,int numstep,int numdimension,double *data){
  int i,j,d;
  
  for (i=0;i<numstep;++i) {
    fprintf(outputfile,"%d ",i+1);
    for (j=0;j<numdimension;++j) {
      fprintf(outputfile," %e",data[i*numdimension+j]);
    }
    fprintf(outputfile,"\n");
  }
  
  return 0;
}


int io_outputmatrix(FILE *inputfile,int numrow,int numcolum,double *data){
  int i,j;
  
  for (i=0;i<numrow;++i) {
    for (j=0;j<numcolum;++j) {
      fprintf(inputfile,"%10.6lf ",data[i*numcolum+j]);
    }
    fprintf(inputfile,"\n");
  }
  fprintf(inputfile,"\n");
    
  return 0;
}
