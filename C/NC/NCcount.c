
#include <stdio.h>
#include <math.h>

#include "PTL.h"
#include "TOPO.h"

#include "NC.h"
#include "NCcount.h"

double NCcount_native_contact_wratio(int numnc, double *crd, 
				     int numatom,int numres,
				     int **ncmap,int **ncmap_res,
				     double *cradii_natatt_CA, double nc_ratio) {
  int i,j,k,na;
  int resi,resj;
  int n;
  double vec[3];
  double len=0.0;
  double Q=0.0;

  na=0;
  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      resi=PTL_resnum(i)-1;
      resj=PTL_resnum(j)-1;
      if ( ncmap_res[resi][resj]==0  && strcmp(AP.IGRAPH[i],"CA",2)==0 && strcmp(AP.IGRAPH[j],"CA",2)==0) {
	len = 0.0;
	for(k=0;k<3;++k){
	  vec[k] = crd[i*3+k]-crd[j*3+k];
	  len += vec[k]*vec[k];
	}
	len = sqrt(len);
	if (len <= cradii_natatt_CA[na]*nc_ratio ) {
	  Q+=1.0;
	}
	++na;
      }
    }
  }

  return (double)Q/numnc;
}

double NCcount_native_contact_wext(int numnc, double *crd, 
				   int numatom,int numres,
				   int **ncmap,int **ncmap_res,
				   double *cradii_natatt_CA, double nc_ext) {
  int i,j,k,na;
  int resi,resj;
  int n;
  double vec[3];
  double len=0.0;
  double Q=0.0;

  na=0;
  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      resi=PTL_resnum(i)-1;
      resj=PTL_resnum(j)-1;
      if (ncmap_res[resi][resj]==0 && strcmp(AP.IGRAPH[i],"CA",2)==0 &&  strcmp(AP.IGRAPH[j],"CA",2)==0) {
	len = 0.0;
	for(k=0;k<3;++k){
	  vec[k] = crd[i*3+k]-crd[j*3+k];
	  len += vec[k]*vec[k];
	}
	len = sqrt(len);
	if (len <= cradii_natatt_CA[na]+nc_ext ) {
	  Q+=1.0;
	}
	++na;
      }
    }
  }

  return (double)Q/numnc;
}

double NCcount_native_contact_AA_wratio(int numnc, double *crd, 
					int numatom,int numres,
					int **ncmap,int **ncmap_res,
					double *cradii_natatt, double nc_ratio) {
  int i,j,k,na;
  int resi,resj;
  int n;
  double vec[3];
  double len=0.0;
  double QAA=0.0;

  na=0;
  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      resi=PTL_resnum(i)-1;
      resj=PTL_resnum(j)-1;
      if (ncmap[i][j]==0) {
	len = 0.0;
	for(k=0;k<3;++k){
	  vec[k] = crd[i*3+k]-crd[j*3+k];
	  len += vec[k]*vec[k];
	}
	len = sqrt(len);
	if (len <= cradii_natatt[na]*nc_ratio && ncmap_res[resi][resj]==0 ) {
	  ncmap_res[resi][resj]=1;
	  QAA+=1.0;
	}
	++na;
      }
    }
  }

  return (double)QAA/numnc;
}

double NCcount_native_contact_AA_wext(int numnc, double *crd, 
				      int numatom,int numres,
				      int **ncmap,int **ncmap_res,
				      double *cradii_natatt, double nc_ext) {
  int i,j,k,na;
  int resi,resj;
  int n;
  double vec[3];
  double len=0.0;
  double QAA=0.0;

  na=0;
  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      resi=PTL_resnum(i)-1;
      resj=PTL_resnum(j)-1;
      if (ncmap[i][j]==0) {
	len = 0.0;
	for(k=0;k<3;++k){
	  vec[k] = crd[i*3+k]-crd[j*3+k];
	  len += vec[k]*vec[k];
	}
	len = sqrt(len);
	if (len <= cradii_natatt[na]+nc_ext && ncmap_res[resi][resj]==0 ) {
	  ncmap_res[resi][resj]=1;
	  QAA+=1.0;
	}
	++na;
      }
    }
  }

  return (double)QAA/numnc;
}

