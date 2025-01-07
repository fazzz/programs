
#include <stdio.h>

#include "PT.h"
#include "TOPO.h"

#include "NC_check.h"

int *make_native_contact_list(int *numnc, double *cord, 
			      int numatom, int numres,
			      double criteria, int *ncmap
			      ) {
  int i,j,k,ii,jj;
  int flag,cflag;
  int *indexncb;
  double length,atom[2][3];

  indexncb=(int *)gcemalloc(sizeof(int)*2);
  *numnc=0;

  for (i=0;i<numres*numres;++i) ncmap[i]=-2;

  for (i=0;i<numres;++i) {
    for (j=i+4;j<numres;++j) {
      flag=NN;
      for (ii=AP.IPRES[i]-1;ii<AP.IPRES[i+1]-1;++ii) {
	for (jj=AP.IPRES[j]-1;jj<AP.IPRES[j+1]-1;++jj) {
	  if (strncmp(AP.IGRAPH[ii],"H",1)!=0 && strncmp(AP.IGRAPH[jj],"H",1)!=0 ) {
	    for (k=0;k<3;++k) {
	      atom[0][k]=cord[ii*3+k];
	      atom[1][k]=cord[jj*3+k];
	    }
	    length=len(atom[0],atom[1]);
	    if (length <= criteria) {
	      flag=NC;
	      break;
	    }
	  }
	}
	if (flag==NC)
	  break;
      }
      if (flag==NC) {
	indexncb=(int *)gcerealloc(indexncb,sizeof(int)*(*numnc+1)*2);
	indexncb[(*numnc)*2]=i+1;
	indexncb[(*numnc)*2+1]=j+1;
	++(*numnc);
	ncmap[i*numres+j]=0;
      }
      else if (flag == NN) {
	ncmap[i*numres+j]=1;
      }
    }
  }
  return indexncb;
}

double count_native_contact(int numnc, double *cord, 
			    int numatom, int numres,
			    int *indexncb, double criteria, int HMODE) {
  int i,j,k,ii,jj;
  int resi,resj;
  int flag,cflag,q=0;
  double length,atom[2][3];

  for (i=0;i<numnc;++i) {
    resi=indexncb[i*2]-1;
    resj=indexncb[i*2+1]-1;
    flag=NN;
    for (ii=AP.IPRES[resi]/*iniatomnumofres(resi)*/-1;PT_resnum(ii)-1==resi;++ii) {
      for (jj=AP.IPRES[resj]/*iniatomnumofres(resj)*/-1;PT_resnum(jj)-1==resj;++jj) {
	if (HMODE!=INC) {
	  if (atomnamencmp(ii,"H",1)!=0 && atomnamencmp(jj,"H",1)!=0) {
	    cflag=ON;
	  }
	  else {
	    cflag=OFF;
	  }
	}
	else {
	  cflag=ON;
	}
	if (cflag==ON) {
	  for (k=0;k<3;++k) {
	    atom[0][k]=cord[ii*3+k];
	    atom[1][k]=cord[jj*3+k];
	  }
	  length=len(atom[0],atom[1]);
	  if (len(atom[0],atom[1]) <= criteria) {
	    flag=NC;
	  }
	}
      }
    }

    if (flag==NC) {
      ++q;
    }
  }

  return (double)q/numnc;
}

int *make_native_contact_list_aa(int *numnc, double *cord, 
				 int numatom, double criteria,
				 int *ncmap, int HMODE) {
  int i,j,k;
  int resi,resj,cflag;
  int *indexncb;
  double length,atom[2][3];

  indexncb=(int *)gcemalloc(sizeof(int)*2);

  *numnc=0;

  for (i=0;i<numatom*numatom;++i) ncmap[i]=-2;

  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      if (HMODE!=INC) {
	//	if (atomnamencmp(i,"H",1)!=0 && atomnamencmp(j,"H",1)!=0) {
	if (strncmp(AP.IGRAPH[i],"H",1)!=0 && strncmp(AP.IGRAPH[j],"H",1)!=0) {
	  cflag=ON;
	}
	else {
	  cflag=OFF;
	}
      }
      else {
	cflag=ON;
      }
      if (cflag==ON) {
	for (k=0;k<3;++k) {
	  atom[0][k]=cord[i*3+k];
	  atom[1][k]=cord[j*3+k];
	}
	length = len(atom[0],atom[1]);
	resi=PT_resnum(i);
	resj=PT_resnum(j);

	if ( resi <  resj-3 ) {
	  if (length <= criteria) {
	    indexncb=(int *)gcerealloc(indexncb,sizeof(int)*((*numnc)+1)*3);
	    indexncb[(*numnc)*2]=i;
	    indexncb[(*numnc)*2+1]=j;
	    ++(*numnc);
	    ncmap[i*numatom+j]=0;
	  }
	  else {
	    ncmap[i*numatom+j]=1;
	  }
      }
	else {
	  ncmap[i*numatom+j]=-1;
	}
      }
    }
  }
  return indexncb;
}

int chaeck_within_3_neibor_res(int num1, int num2) {
  int numres1,numres2;

  numres1=PT_resnum(num1);
  numres2=PT_resnum(num2);

  if (numres2 < numres1+3)
    return YES;
  else 
    return NO;
}

