#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double *data;

double calc_rmsd(int num_initial, int num_final,int numvalue, int nv);
double calc_least_sqrt_devi(int num_initial, int num_final,int numvalue, int nv);
int scan_ene_pro(FILE *inputfile,int numstep, int numvalue);

int main(int argc,char* argv[]){
  int numstep,numvalue,nv,num_initial,num_final;
  double rmsd,drift;
  double dt;
  char *inputfilename,*inputfilename2,*inputfilename3,*outputfilename,outputfilename2[100];
  FILE *inputfile,*inputfile2;
  FILE *outputfile;

  if (argc < 3) {
    printf("USAGE: ./calc_ene_conservation.exe inputfilename(enepro) inputfilename2(cond)  outputfilename\n");
    exit(1);
  }

  inputfilename =  *++argv;
  inputfilename2 = *++argv;
  outputfilename = *++argv;


  if((inputfile2=fopen(inputfilename2,"r"))==NULL)  {
    printf("There is not %s\n",inputfilename2);
    exit(1);
  }
  fscanf(inputfile2,"%d",&numstep);
  fscanf(inputfile2,"%d",&numvalue);
  fscanf(inputfile2,"%d",&nv);
  fscanf(inputfile2,"%d",&num_initial);
  fscanf(inputfile2,"%d",&num_final);
  fscanf(inputfile2,"%lf",&dt);
  fclose(inputfile2);

  if((inputfile=fopen(inputfilename,"r"))==NULL)  {
    printf("There is not %s\n",inputfilename);
    exit(1);
  }
  data = (double *)malloc(sizeof(double)*numstep*numvalue);
  scan_ene_pro(inputfile,numstep,numvalue);
  fclose(inputfile);

  rmsd=calc_rmsd(num_initial,num_final,numvalue,nv);
  drift=calc_least_sqrt_devi(num_initial,num_final,numvalue,nv);

  printf("%lf %e %e \n",dt,rmsd,drift);
  free(data);

  return 1;

}

int scan_ene_pro(FILE *inputfile,int numstep, int numvalue) {
  int i,j;
  double d;

  for (i=0;i<numstep;++i)  {
    for (j=0;j<numvalue;++j)  {
      fscanf(inputfile,"%lf",&d);
      data[i*numvalue+j]=d;
    }
  }
  
  return 0;
}

double calc_rmsd(int num_initial, int num_final,int numvalue, int nv) {
  int i,j,numstep;
  double val_initial;
  double ave=0.0;
  double var=0.0;
  double rmsd;

  numstep=num_final-num_initial;
  val_initial = data[num_initial];
  for (i=num_initial;i<num_final;++i) {
    for (j=0;j<numvalue;++j) {
      if (j==nv) {
	ave += data[i*numvalue+j];
	var += data[i*numvalue+j]*data[i*numvalue+j];
      }
    }
  }

  ave = ave/numstep;
  var = var/numstep;
  var = var - ave*ave;

  if (var != 0.0 && val_initial != 0.0) {
    rmsd = sqrt(var)/abs(val_initial);
  }
  else {
    rmsd = 0.0;
  }

  return rmsd;
}

double calc_least_sqrt_devi(int num_initial, int num_final,int numvalue, int nv) {
  int i,j,numstep;
  double sum_x=0.0,sum_x_2=0.0,sum_xy=0.0,sum_y=0.0;
  double a,b;

  numstep=num_final-num_initial;

  for (i=num_initial;i<num_final;++i) {
    for (j=0;j<numvalue;++j) {
      if (j==nv) {
	sum_x   += i;
	sum_xy  += data[i*numvalue+j]*i;
	sum_x_2 += i*i;
	sum_y   += data[i*numvalue+j];
      }
    }
  }

  a=1.0/(sum_x_2*numstep-sum_x*sum_x)*(numstep*sum_xy-sum_x*sum_y);
  b=1.0/(sum_x_2*numstep-sum_x*sum_x)*(-sum_x*sum_xy-sum_x_2*sum_y);

  return a;
}
