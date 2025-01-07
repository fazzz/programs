#include <stdio.h>
#include <math.h>

//#include "ABA.h" // 2014-06-18
#include "ABAb.h"  // 2014-06-18

#include "EF.h" // 2-14-06-17

#define ON 1
#define OFF 0

void ABAs_restat_read(char *inputvelofilename,int numclut,
		      double *correct,double **correctt_Term,double **correctt_Term2,double correct_s[6]
		      ,int MODE,int TERMMODE) {
  int i,j;
  FILE *inputvelofile;

  inputvelofile=efopen(inputvelofilename,"r");

  if (TERMMODE==ON) {
    for (i=0;i<6;++i) 
      for (j=0;j<6;++j)
	fscanf(inputvelofile,"%lf",&correctt_Term[i][j]);
    for (i=0;i<6;++i)  
      for (j=0;j<6;++j)
	fscanf(inputvelofile,"%lf",&correctt_Term2[i][j]);
  }

  for (i=0;i<numclut;++i)
    for (j=0;j<6;++j)
      fscanf(inputvelofile,"%lf",&correct[i*6+j]);

  if (MODE==NVT)  for (i=0;i<6;++i)  fscanf(inputvelofile,"%lf",&correct_s[i]);

  fclose(inputvelofile);
}

void ABAs_restat_write_vel(char *restatvelofilename,int numclut,
			   double *correct,double **correctt_Term,double **correctt_Term2,double correct_s[6]
			   ,int MODE,int TERMMODE) {
  int i,j,l=0;
  FILE *restatvelofile;

  restatvelofile=efopen(restatvelofilename,"w");

  if (TERMMODE==ON) {
    for (i=0;i<6;++i)  {
      for (j=0;j<6;++j) {
	fprintf(restatvelofile,"%10.8lf ",correctt_Term[i][j]);
	++l;
	if (l%10==0) fprintf(restatvelofile,"\n");
      }
    }
    for (i=0;i<6;++i)  {
      for (j=0;j<6;++j) {
	fprintf(restatvelofile,"%10.8lf ",correctt_Term2[i][j]);
	++l;
	if (l%10==0) fprintf(restatvelofile,"\n");
      }
    }
  }

  for (i=0;i<numclut;++i) {
    for (j=0;j<6;++j) {
      fprintf(restatvelofile,"%10.8lf ",correct[i*6+j]);
      if (l%10==0) fprintf(restatvelofile,"\n");
    }
  }

  if (MODE==NVT) {
    for (i=0;i<6;++i) {
      fprintf(restatvelofile,"%10.8lf\n",correct_s[i]);
      if (l%10==0) fprintf(restatvelofile,"\n");
    }
  }

  fclose(restatvelofile);
}

void ABAs_restat_read_new(char *inputvelofilename,int numclut,
			  double *correct,double **correctt_Term,double **correctt_Term2,double correct_gzi[5],
			  int MODE,int TERMMODE) {
  int i,j;
  FILE *inputvelofile;

  inputvelofile=efopen(inputvelofilename,"r");

  if (TERMMODE==ON) {
    for (i=0;i<6;++i) 
      for (j=0;j<6;++j)
	fscanf(inputvelofile,"%lf",&correctt_Term[i][j]);
    for (i=0;i<6;++i)  
      for (j=0;j<6;++j)
	fscanf(inputvelofile,"%lf",&correctt_Term2[i][j]);
  }

  for (i=0;i<numclut;++i)
    for (j=0;j<6;++j)
      fscanf(inputvelofile,"%lf",&correct[i*6+j]);

  if (MODE==NVT)  for (i=0;i<5;++i)  fscanf(inputvelofile,"%lf",&correct_gzi[i]);

  fclose(inputvelofile);
}

void ABAs_restat_write_vel_new(char *restatvelofilename,int numclut,
			       double *correct,double **correctt_Term,double **correctt_Term2,double correct_gzi[5],
			       int MODE,int TERMMODE) {
  int i,j,l=0;
  FILE *restatvelofile;

  restatvelofile=efopen(restatvelofilename,"w");

  if (TERMMODE==ON) {
    for (i=0;i<6;++i)  {
      for (j=0;j<6;++j) {
	fprintf(restatvelofile,"%10.8lf ",correctt_Term[i][j]);
	++l;
	if (l%10==0) fprintf(restatvelofile,"\n");
      }
    }
    for (i=0;i<6;++i)  {
      for (j=0;j<6;++j) {
	fprintf(restatvelofile,"%10.8lf ",correctt_Term2[i][j]);
	++l;
	if (l%10==0) fprintf(restatvelofile,"\n");
      }
    }
  }

  for (i=0;i<numclut;++i) {
    for (j=0;j<6;++j) {
      fprintf(restatvelofile,"%10.8lf ",correct[i*6+j]);
      if (l%10==0) fprintf(restatvelofile,"\n");
    }
  }

  if (MODE==NVT) {
    for (i=0;i<5;++i) {
      fprintf(restatvelofile,"%10.8lf\n",correct_gzi[i]);
      if (l%10==0) fprintf(restatvelofile,"\n");
    }
  }

  fclose(restatvelofile);
}

void ABAs_restat_read_new_mvV(char *inputvelofilename,int numclut,
			      double *qvel_b1,double *qvel_b2,
			      double *vel_Term_b1,double *vel_Term_b2,
			      double zeta,
			      int MODE,int TERMMODE) {
  int i,j;
  FILE *inputvelofile;

  inputvelofile=efopen(inputvelofilename,"r");

  if (TERMMODE==ON) {
    for (i=0;i<6;++i) fscanf(inputvelofile,"%lf",&vel_Term_b1[i]);
    for (i=0;i<6;++i) fscanf(inputvelofile,"%lf",&vel_Term_b2[i]);
  }

  for (i=0;i<numclut;++i)
    fscanf(inputvelofile,"%lf",&qvel_b1[i]);

  for (i=0;i<numclut;++i)
    fscanf(inputvelofile,"%lf",&qvel_b2[i]);

  if (MODE==NVT) fscanf(inputvelofile,"%lf",&zeta);

  fclose(inputvelofile);
}

void ABAs_restat_write_vel_new_mvV(char *restatvelofilename,int numclut,
				   double *qvel_b1,double *qvel_b2,
				   double *vel_Term_b1,double *vel_Term_b2,
				   double zeta,
				   int MODE,int TERMMODE) {
  int i,j,l=0;
  FILE *restatvelofile;

  restatvelofile=efopen(restatvelofilename,"w");

  if (TERMMODE==ON) {
    for (i=0;i<6;++i)  {
	fprintf(restatvelofile,"%10.8lf ",vel_Term_b1[i]);
	++l;
	if (l%10==0) fprintf(restatvelofile,"\n");
    }
    for (i=0;i<6;++i)  {
    	fprintf(restatvelofile,"%10.8lf ",vel_Term_b2[i]);
	++l;
	if (l%10==0) fprintf(restatvelofile,"\n");
    }
  }

  for (i=0;i<numclut;++i) {
    fprintf(restatvelofile,"%10.8lf ",qvel_b1[i]);
    if (l%10==0) fprintf(restatvelofile,"\n");
  }

  for (i=0;i<numclut;++i) {
    fprintf(restatvelofile,"%10.8lf ",qvel_b2[i]);
    if (l%10==0) fprintf(restatvelofile,"\n");
  }

  if (MODE==NVT) {
    fprintf(restatvelofile,"%10.8lf\n",zeta);
    if (l%10==0) fprintf(restatvelofile,"\n");
  }

  fclose(restatvelofile);
}

void ABAs_restat_read_chain(char *inputvelofilename,int numclut,
			    double *correct,double **correctt_Term,double **correctt_Term2,double **correct_gzi,
			    int M,
			    int MODE,int TERMMODE) {
  int i,j;
  FILE *inputvelofile;

  inputvelofile=efopen(inputvelofilename,"r");

  if (TERMMODE==ON) {
    for (i=0;i<6;++i) 
      for (j=0;j<6;++j)
	fscanf(inputvelofile,"%lf",&correctt_Term[i][j]);
    for (i=0;i<6;++i)  
      for (j=0;j<6;++j)
	fscanf(inputvelofile,"%lf",&correctt_Term2[i][j]);
  }

  for (i=0;i<numclut;++i)
    for (j=0;j<6;++j)
      fscanf(inputvelofile,"%lf",&correct[i*6+j]);

  if (MODE==NVT)  for (i=0;i<M;++i) for (j=0;j<6;++j)  fscanf(inputvelofile,"%lf",&correct_gzi[i][j]);

  fclose(inputvelofile);
}

void ABAs_restat_write_vel_chain(char *restatvelofilename,int numclut,
				 double *correct,double **correctt_Term,double **correctt_Term2,double **correct_gzi,
				 int M,
				 int MODE,int TERMMODE) {
  int i,j,l=0;
  FILE *restatvelofile;

  restatvelofile=efopen(restatvelofilename,"w");

  if (TERMMODE==ON) {
    for (i=0;i<6;++i)  {
      for (j=0;j<6;++j) {
	fprintf(restatvelofile,"%10.8lf ",correctt_Term[i][j]);
	++l;
	if (l%10==0) fprintf(restatvelofile,"\n");
      }
    }
    for (i=0;i<6;++i)  {
      for (j=0;j<6;++j) {
	fprintf(restatvelofile,"%10.8lf ",correctt_Term2[i][j]);
	++l;
	if (l%10==0) fprintf(restatvelofile,"\n");
      }
    }
  }

  for (i=0;i<numclut;++i) {
    for (j=0;j<6;++j) {
      fprintf(restatvelofile,"%10.8lf ",correct[i*6+j]);
      if (l%10==0) fprintf(restatvelofile,"\n");
    }
  }

  if (MODE==NVT) {
    for (i=0;i<M;++i) {
      for (j=0;j<6;++j) {
	fprintf(restatvelofile,"%10.8lf\n",correct_gzi[i][j]);
      }
      if (l%10==0) fprintf(restatvelofile,"\n");
    }
  }

  fclose(restatvelofile);
}
