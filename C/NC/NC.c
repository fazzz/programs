
#include <stdio.h>

#include "PTL.h"
#include "TOPO.h"

#include "NC.h"

int *make_native_contact_list(int *numnc, double *cord, 
			      int numatom, int numres,
			      double criteria, int **ncmap
			      ) {
  int i,j,k,ii,jj;
  int flag,cflag;
  int *indexncb;
  double length,atom[2][3];

  indexncb=(int *)gcemalloc(sizeof(int)*2);
  *numnc=0;

  for (i=0;i<numres;++i) 
    for (j=0;j<numres;++j) 
      ncmap[i][j]=-2;


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
	ncmap[i][j]=0;
      }
      else if (flag == NN) {
	ncmap[i][j]=1;
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
    for (ii=AP.IPRES[resi]/*iniatomnumofres(resi)*/-1;PTL_resnum(ii)-1==resi;++ii) {
      for (jj=AP.IPRES[resj]/*iniatomnumofres(resj)*/-1;PTL_resnum(jj)-1==resj;++jj) {
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

double count_native_contact_aa(int numnc, double *crd, 
			       int numatom,int numres,
			       int **ncmap, double *cradii_natatt) {
  int i,j,k,na;
  int n;
  double vec[3];
  double len=0.0;
  double QAA=0.0;

  int **ncmapres;

  ncmapres=(int **)gcemalloc(sizeof(int *)*numres);
  for (i=0;i<numres;++i)  ncmapres[i]=(int *)gcemalloc(sizeof(int)*numres);

  for (i=0;i<numres;++i)  for (j=0;j<numres;++j) ncmapres[i][j]=-1;

  na=0;
  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      if (ncmap[i][j]==0) {
	len = 0.0;
	for(k=0;k<3;++k){
	  vec[k] = crd[i*3+k]-crd[j*3+k];
	  len += vec[k]*vec[k];
	}
	len = sqrt(len);
	if (len <= 1.2*cradii_natatt[na]) {
	  if (ncmapres[PTL_resnum(i)][PTL_resnum(j)]==-1) {
	    ncmapres[PTL_resnum(i)][PTL_resnum(j)]=0;
	    ++QAA;
	  }
	}
	++na;
      }
    }
  }

  return (double)QAA/numnc;
}

double count_native_contact_hybrid(int numncaa, double *crd, 
				   int numatom,int numres,
				   int **ncmap, double *cradii_natatt) {
  int i,j,k,na;
  int n;
  double vec[3];
  double len=0.0;
  double QAA=0.0;

  na=0;
  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      if (ncmap[i][j]==0) {
	len = 0.0;
	for(k=0;k<3;++k){
	  vec[k] = crd[i*3+k]-crd[j*3+k];
	  len += vec[k]*vec[k];
	}
	len = sqrt(len);
	if (len <= cradii_natatt[na]+1.0) {
	  ++QAA;
	}
	++na;
      }
    }
  }

  return (double)QAA/numncaa;
}

double count_native_contact_ca(int numnc,double *crd, 
			       int numres,
			       int **ncmapres,double *cradii_natatt) {
  int i,j,k,na;
  int ii,jj;
  double vec[3];
  double len=0.0;
  double QAA=0.0;

  na=0;
  for (i=0;i<numres;++i) {
    for (j=i+1;j<numres;++j) {
      if (ncmapres[i][j]==0 ) {
	ii=PTL_canum_fromresnum(i);
	jj=PTL_canum_fromresnum(j);
	len = 0.0;
	for(k=0;k<3;++k){
	  vec[k] = crd[ii*3+k]-crd[jj*3+k];
	  len += vec[k]*vec[k];
	}
	len = sqrt(len);
	if (len <= 1.2*cradii_natatt[na]) ++QAA;
	++na;
      }
    }
  }

  return (double)QAA/numnc;
}


int *make_native_contact_list_aa(int *numnc, double *cord, 
				 int numatom, double criteria,
				 int **ncmap, int HMODE) {
  int i,j,k;
  int resi,resj,cflag;
  int *indexncb;
  double length,atom[2][3];

  indexncb=(int *)gcemalloc(sizeof(int)*2);

  *numnc=0;

  for (i=0;i<numatom;++i) for (j=0;j<numatom;++j) ncmap[i][j]=-2;

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
	resi=PTL_resnum(i);
	resj=PTL_resnum(j);

	if ( resi <  resj-3 ) {
	  if (length <= criteria) {
	    indexncb=(int *)gcerealloc(indexncb,sizeof(int)*((*numnc)+1)*3);
	    indexncb[(*numnc)*2]=i;
	    indexncb[(*numnc)*2+1]=j;
	    ++(*numnc);
	    ncmap[i][j]=0;
	  }
	  else {
	    ncmap[i][j]=1;
	  }
      }
	else {
	  ncmap[i][j]=-1;
	}
      }
    }
  }
  return indexncb;
}

int *make_native_contact_list_aa_2(int *numncaa, int *numncres,double *cord, 
				   int numatom, int numres, double criteria,
				   int **ncmap, int HMODE) {
  int i,j,k,n;
  int resi,resj,cflag;
  int *indexncb;
  double length,atom[2][3];
  
  int **ncmapres;

  indexncb=(int *)gcemalloc(sizeof(int)*2);

  *numncaa=0;
  *numncres=0;

  for (i=0;i<numatom;++i) 
    for (j=0;j<numatom;++j) ncmap[i][j]=-2;

  ncmapres=(int **)gcemalloc(sizeof(int *)*numres);
  for (i=0;i<numres;++i) ncmapres[i]=(int *)gcemalloc(sizeof(int)*numres);
  for (i=0;i<numres;++i)  for (j=0;j<numres;++j) ncmapres[i][j]=-1;

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
	resi=PTL_resnum(i);
	resj=PTL_resnum(j);

	if ( resi <  resj-3 ) {
	  if (length <= criteria) {
	    indexncb=(int *)gcerealloc(indexncb,sizeof(int)*((*numncaa)+1)*3);
	    indexncb[(*numncaa)*2]=i;
	    indexncb[(*numncaa)*2+1]=j;
	    ++(*numncaa);
	    ncmap[i][j]=0;
	    ///////////////////////////////////////////////////
	    n=PTL_resnum(i)*numres+PTL_resnum(j);
	    if (ncmapres[PTL_resnum(i)][PTL_resnum(j)]==-1) {
	      ncmapres[PTL_resnum(i)][PTL_resnum(j)]=0;
	      ++(*numncres);
	    }
	    ///////////////////////////////////////////////////
	  }
	  else {
	    ncmap[i][j]=1;
	  }
	}
	else {
	  ncmap[i][j]=-1;
	}
      }
    }
  }
  return indexncb;
}

int *make_native_contact_list_aa_3(int *numncaa, int *numncres,double *cord, 
				   int numatom, int numres, double criteria,
				   int **ncmap,  int **ncmapres,int HMODE) {
  int i,j,k,n;
  int resi,resj,cflag;
  int *indexncb;
  double length,atom[2][3];
  
  indexncb=(int *)gcemalloc(sizeof(int)*2);

  *numncaa=0;
  *numncres=0;

  for (i=0;i<numatom;++i) for (j=0;j<numatom;++j) ncmap[i][j]=-2;

  for (i=0;i<numres;++i)  for (j=0;j<numres;++j) ncmapres[i][j]=-1;

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
	resi=PTL_resnum(i);
	resj=PTL_resnum(j);

	if ( resi <  resj-3 ) {
	  if (length <= criteria) {
	    indexncb=(int *)gcerealloc(indexncb,sizeof(int)*((*numncaa)+1)*3);
	    indexncb[(*numncaa)*2]=i;
	    indexncb[(*numncaa)*2+1]=j;
	    ++(*numncaa);
	    ncmap[i][j]=0;
	    ///////////////////////////////////////////////////
	    n=PTL_resnum(i)*numres+PTL_resnum(j);
	    if (ncmapres[PTL_resnum(i)][PTL_resnum(j)]==-1) {
	      if (strncmp(AP.LABERES[PTL_resnum(i)],"ACE",3) != 0 && strncmp(AP.LABERES[PTL_resnum(i)],"NME",3) != 0 
		  && strncmp(AP.LABERES[PTL_resnum(j)],"ACE",3) != 0 && strncmp(AP.LABERES[PTL_resnum(j)],"NME",3) != 0) {
		ncmapres[PTL_resnum(i)][PTL_resnum(j)]=0;
		++(*numncres);
		//		printf("%d-%d %d\n",resi,resj,n);
	      }
	    }
	    ///////////////////////////////////////////////////
	  }
	  else {
	    ncmap[i][j]=1;
	  }
	}
	else {
	  ncmap[i][j]=-1;
	}
      }
    }
  }
  return indexncb;
}

int *make_native_contact_list_aa_3_nadjacent(int *numncaa, int *numncres,double *cord, 
					     int numatom, int numres, double criteria,
					     int **ncmap,  int **ncmapres,int HMODE) {
  int i,j,k,n;
  int resi,resj,cflag;
  int *indexncb;
  double length,atom[2][3];
  
  indexncb=(int *)gcemalloc(sizeof(int)*2);

  *numncaa=0;
  *numncres=0;

  for (i=0;i<numatom;++i) for (j=0;j<numatom;++j) ncmap[i][j]=-2;

  for (i=0;i<numres;++i)  for (j=0;j<numres;++j) ncmapres[i][j]=-1;

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
	resi=PTL_resnum(i);
	resj=PTL_resnum(j);

	if ( resi <  resj-1/*3*/ ) {
	  if (length <= criteria) {
	    indexncb=(int *)gcerealloc(indexncb,sizeof(int)*((*numncaa)+1)*2/*3*/);
	    indexncb[(*numncaa)*2+0]=i;
	    indexncb[(*numncaa)*2+1]=j;
	    ++(*numncaa);
	    ncmap[i][j]=0;
	    ///////////////////////////////////////////////////
	    n=PTL_resnum(i)*numres+PTL_resnum(j);
	    if (ncmapres[PTL_resnum(i)][PTL_resnum(j)]==-1) {
	      if (strncmp(AP.LABERES[PTL_resnum(i)],"ACE",3) != 0 && strncmp(AP.LABERES[PTL_resnum(i)],"NME",3) != 0 
		  && strncmp(AP.LABERES[PTL_resnum(j)],"ACE",3) != 0 && strncmp(AP.LABERES[PTL_resnum(j)],"NME",3) != 0) {
		ncmapres[PTL_resnum(i)][PTL_resnum(j)]=0;
		++(*numncres);
		//		printf("%d-%d %d\n",resi,resj,n);
	      }
	    }
	    ///////////////////////////////////////////////////
	  }
	  else {
	    ncmap[i][j]=1;
	  }
	}
	else {
	  ncmap[i][j]=-1;
	}
      }
    }
  }
  return indexncb;
}

int *make_native_contact_list_aa_wnibnum(int *numncaa, int *numncres,double *cord, 
					 int numatom, int numres, double criteria,
					 int **ncmap,  int **ncmapres,int HMODE, int nibnum) {
  int i,j,k,n;
  int resi,resj,cflag;
  int *indexncb;
  double length,atom[2][3];
  
  indexncb=(int *)gcemalloc(sizeof(int)*2);

  *numncaa=0;
  *numncres=0;

  for (i=0;i<numatom;++i) for (j=0;j<numatom;++j) ncmap[i][j]=-2;

  for (i=0;i<numres;++i)  for (j=0;j<numres;++j) ncmapres[i][j]=-1;

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
	resi=PTL_resnum(i)-1;
	resj=PTL_resnum(j)-1;

	if ( resi <  resj-nibnum/*3*/ ) {
	  if (length <= criteria) {
	    indexncb=(int *)gcerealloc(indexncb,sizeof(int)*((*numncaa)+1)*2/*3*/);
	    indexncb[(*numncaa)*2+0]=i;
	    indexncb[(*numncaa)*2+1]=j;
	    ++(*numncaa);
	    ncmap[i][j]=0;
	    ///////////////////////////////////////////////////
	    if (ncmapres[resi][resj]==-1) {
	      if (   strncmp(AP.LABERES[resi],"ACE",3) != 0 && strncmp(AP.LABERES[resi],"NME",3) != 0 
		  && strncmp(AP.LABERES[resj],"ACE",3) != 0 && strncmp(AP.LABERES[resj],"NME",3) != 0) {
		ncmapres[resi][resj]=0;
		++(*numncres);
	      }
	    }
	    ///////////////////////////////////////////////////
	  }
	  else {
	    ncmap[i][j]=1;
	  }
	}
	else {
	  ncmap[i][j]=-1;
	}
      }
    }
  }
  return indexncb;
}

int *make_native_contact_list_aa_3_nadjacent_2(int *numncaa, int *numncres,double *cord, 
					       int numatom, int numres, double criteria,
					       int **ncmap,  int **ncmapres,int HMODE) {
  int i,j,k,n;
  int resi,resj,cflag;
  int *indexncb;
  double length,atom[2][3];
  
  indexncb=(int *)gcemalloc(sizeof(int)*2);

  *numncaa=0;
  *numncres=0;

  for (i=0;i<numatom;++i) for (j=0;j<numatom;++j) ncmap[i][j]=-2;

  for (i=0;i<numres;++i)  for (j=0;j<numres;++j) ncmapres[i][j]=-1;

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
	resi=PTL_resnum(i);
	resj=PTL_resnum(j);

	if ( resi <  resj-3 ) {
	  if (length <= criteria) {
	    indexncb=(int *)gcerealloc(indexncb,sizeof(int)*((*numncaa)+1)*2/*3*/);
	    indexncb[(*numncaa)*2+0]=i;
	    indexncb[(*numncaa)*2+1]=j;
	    ++(*numncaa);
	    ncmap[i][j]=0;
	    ///////////////////////////////////////////////////
	    n=PTL_resnum(i)*numres+PTL_resnum(j);
	    if (ncmapres[PTL_resnum(i)][PTL_resnum(j)]==-1) {
	      if (strncmp(AP.LABERES[PTL_resnum(i)],"ACE",3) != 0 && strncmp(AP.LABERES[PTL_resnum(i)],"NME",3) != 0 
		  && strncmp(AP.LABERES[PTL_resnum(j)],"ACE",3) != 0 && strncmp(AP.LABERES[PTL_resnum(j)],"NME",3) != 0) {
		ncmapres[PTL_resnum(i)][PTL_resnum(j)]=0;
		++(*numncres);
		//		printf("%d-%d %d\n",resi,resj,n);
	      }
	    }
	    ///////////////////////////////////////////////////
	  }
	  else {
	    ncmap[i][j]=1;
	  }
	}
	else {
	  ncmap[i][j]=-1;
	}
      }
    }
  }
  return indexncb;
}

int chaeck_within_3_neibor_res(int num1, int num2) {
  int numres1,numres2;

  numres1=PTL_resnum(num1);
  numres2=PTL_resnum(num2);

  if (numres2 < numres1+3)
    return YES;
  else 
    return NO;
}

