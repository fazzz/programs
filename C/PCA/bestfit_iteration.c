//#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "const.h"
#include "bestfit.h"
#include "quaternion.h"

void calc_ave_crd(double *traj,double coord_ref[MAXNUMATOM][3],int time,int numatom);

void bestfit_iteration(double *traj,double *crd_ave,
		       int numiteration,int time,int numatom,double mass[MAXNUMATOM])
{
  int i,j,k,l;
  double coord_ref[MAXNUMATOM][3],coord_tag[MAXNUMATOM][3],velo_tag[MAXNUMATOM][3];
  double rmsd[MAXTIME];
  FILE *output;

  for (i=0;i<numiteration;++i){
    calc_ave_crd(traj,coord_ref,time,numatom);
    for (j=0;j<time;++j){
      for (k=0;k<numatom;++k)
	for (l=0;l<3;++l)
	  coord_tag[k][l]=traj[j*numatom*3+k*3+l];
      rmsd[i]+=bestfit(coord_ref,coord_tag,velo_tag,mass,numatom);
      for (k=0;k<numatom;++k)
	for (l=0;l<3;++l)
	  traj[j*numatom*3+k*3+l]=coord_tag[k][l];
    }
  }

  for (k=0;k<numatom;++k)
	for (l=0;l<3;++l)
	  crd_ave[k*3+l] = coord_ref[k][l];

  if ((output=fopen("rmsd.txt","w"))==NULL){
    exit(1);
  }
  for (i=0;i<10;++i){
    fprintf(output,"%d %lf\n", i, rmsd[i]/time);
  }
  fclose(output);


}

void calc_ave_crd(double *traj,double coord_ref[MAXNUMATOM][3],int time,int numatom)
{
  int i,j,k;

  for (i=0;i<time;++i)
    for (j=0;j<numatom;++j)
      for (k=0;k<3;++k)
	coord_ref[j][k]=(i*coord_ref[j][k]+traj[i*numatom*3+j*3+k])/(i+1);
}
