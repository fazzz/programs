

#include "PROTOPO.h"

char *atomtplist[2000]= {
  "CA ","N  ",
  "CB ","CA ",
  "CG ","CB ",
  "CG1","CB ",
  "CG2","CB ",
  "CG3","CB ",
  "CD ","CG ",
  "CD1","CG ",
  "CD2","CG ",
  "CE ","CD ",
  "CE1","ND1",
  "CE2","CD1",
  "CE3","CD2",
  "CZ ","CE ",
  "CZ2","CE2",
  "CZ3","CE3",
  "CH2","CZ2",
  "HA ","CA ",
  "HB ","CB ",
  "HD ","CD ",
  "HD1","CD1",
  "HD2","CD2",
  "HD3","CD3",
  "HG ","CG ",
  "HG1","CG1",
  "HG2","CG2",
  "HG3","CG3",
  "HE ","CE ",
  "HE1","CE1",
  "HE2","CE2",
  "HE3","CE3",
  "HZ ","CZ ",
  "HZ2","CZ2",
  "HZ3","CZ3",
  "HH2","CH2",
  "HN ","N  ",
  "HH ","OH ",
  "H3 ","C  ",
  "N  ","C  ",
  "NZ ","CE ",
  "NE ","CD ",
  "NE1","CD1",
  "NE2","CD2",
  "ND ","CG ",
  "ND1","CG ",
  "ND2","CG ",
  "O  ","C  ",
  "OG ","CB ",
  "OG1","CB ",
  "OD ","CG ",
  "OD1","CG ",
  "OD2","CG ",
  "OE ","CD ",
  "OE1","CD ",
  "OE2","CD ",
  "OH ","CZ ",
  "SG ","CB ",
  "SD ","CG ",
  "C  ","C  "
};

/****************************/
/* char *atomtplist[600]= { */
/*   "CA ","N  ",	    */
/*   "CB ","CA ",	    */
/*   "CG ","CB ",	    */
/*   "CG1","CB ",	    */
/*   "CG2","CB ",	    */
/*   "CG3","CB ",	    */
/*   "CD ","CG ",	    */
/*   "CD1","CG ",	    */
/*   "CD2","CG ",	    */
/*   "CE ","CD ",	    */
/*   "CE1","CD ",	    */
/*   "CE2","CD2",	    */
/*   "CE3","CD2",	    */
/*   "CZ ","CE ",	    */
/*   "CZ1","CE2",	    */
/*   "CZ2","CE3",	    */
/*   "CH ","CZ ",	    */
/*   "CH2","CZ2",	    */
/*   "HA ","CA ",	    */
/*   "HB ","CB ",	    */
/*   "HD ","CD ",	    */
/*   "HD1","CD1",	    */
/*   "HD2","CD2",	    */
/*   "HD3","CD3",	    */
/*   "HG ","CG ",	    */
/*   "HG1","CG1",	    */
/*   "HG2","CG2",	    */
/*   "HG3","CG3",	    */
/*   "HE ","CE ",	    */
/*   "HE1","CE1",	    */
/*   "HE2","CE2",	    */
/*   "HZ ","CZ ",	    */
/*   "HZ1","CZ1",	    */
/*   "HZ2","CZ2",	    */
/*   "HH2","CH2",	    */
/*   "HN ","N  ",	    */
/*   "HH ","OH ",	    */
/*   "H3 ","C  ",	    */
/*   "N  ","C  ",	    */
/*   "NZ ","CE ",	    */
/*   "NE ","CD ",	    */
/*   "NE2","CD2",	    */
/*   "ND ","CG ",	    */
/*   "ND1","CG ",	    */
/*   "ND2","CG ",	    */
/*   "O  ","C  ",	    */
/*   "OG ","CB ",	    */
/*   "OG1","CB ",	    */
/*   "OD ","CG ",	    */
/*   "OD1","CG ",	    */
/*   "OD2","CG ",	    */
/*   "OE ","CD ",	    */
/*   "OE1","CD ",	    */
/*   "OE2","CD ",	    */
/*   "OH ","CZ ",	    */
/*   "SG ","CB ",	    */
/*   "SD ","CG ",	    */
/*   "C  ","C  "	    */
/* };			    */
/****************************/

char *atomtplist_c2[100]= {
  "HG ","OG ",
  "C  ","CA ",
  "CD1","CG1",
  "CE ","SD ",
  "CE1","CD ",
  "CE2","CD2",
  "CZ ","NEX",
  "HG1","OG1",
  "HD1","ND1",
  "HD2","ND2",
  "HE1","NE1",
  "HE2","NE2",
  "HZ ","NZ ",
  "CE3","CD2",
  "N  ","HN "
};

char *atomtplist_c3[numatomtp*2]= {
  "C   ","H3 1",
  "CE1","CD1",
  "HD2 ","OD2 ",
  "HE2 ","OE2 ",
  "CZ  ","NEX "
};

char *atomtplist_c4[numatomtp*2]= {
  "CZ  ","CE1 "
};

char *atomtplist_c5[numatomtp*2]= {
  "CZ  ","CE2 "
};


int make_bp(int numatom,int *bp,int **bp_f , int *numb, char *name_atom_list, int aromaflag) {
  int i,j,k;
  int numbd=0;

  int flagimso=0;
  int flagimdo=0;

  gchar *atomnameaskey;
  gchar *atomname;
  char named[4];
  char *name;
  gint *atonum;
  GHashTable *bplist,*bplist_c2,*bplist_c3,*bplist_c4,*bplist_c5;

  int ap[2],apflag;

  bplist=g_hash_table_new_full(g_str_hash,g_str_equal,g_free,g_free);
  bplist_c2=g_hash_table_new_full(g_str_hash,g_str_equal,g_free,g_free);
  bplist_c3=g_hash_table_new_full(g_str_hash,g_str_equal,g_free,g_free);
  bplist_c4=g_hash_table_new_full(g_str_hash,g_str_equal,g_free,g_free);
  bplist_c5=g_hash_table_new_full(g_str_hash,g_str_equal,g_free,g_free);

  for (i=0;i<numatomtp2;++i) {
    atomnameaskey=g_strdup(atomtplist[i*2]);
    atomname=g_new(gchar,1);atomname=atomtplist[i*2+1];
    g_hash_table_insert(bplist,atomnameaskey,atomname);
    if ((name=g_hash_table_lookup(bplist,atomtplist[i*2]))!=NULL) {
      ;
    }
  }
  for (i=0;i<numatomtp3;++i) {
    atomnameaskey=g_strdup(atomtplist_c2[i*2]);
    atomname=g_new(gchar,1);atomname=atomtplist_c2[i*2+1];
    g_hash_table_insert(bplist_c2,atomnameaskey,atomname);
    if ((name=g_hash_table_lookup(bplist_c2,atomtplist_c2[i*2]))!=NULL) {
      ;
    }
  }
  for (i=0;i<numatomtp4;++i) {
    atomnameaskey=g_strdup(atomtplist_c3[i*2]);
    atomname=g_new(gchar,1);atomname=atomtplist_c3[i*2+1];
    g_hash_table_insert(bplist_c3,atomnameaskey,atomname);
    if ((name=g_hash_table_lookup(bplist_c3,atomtplist_c3[i*2]))!=NULL) {
      ;
    }
  }
  atomnameaskey=g_strdup(atomtplist_c4[0]);
  atomname=g_new(gchar,1);atomname=atomtplist_c4[1];
  g_hash_table_insert(bplist_c4,atomnameaskey,atomname);

  atomnameaskey=g_strdup(atomtplist_c5[0]);
  atomname=g_new(gchar,1);atomname=atomtplist_c5[1];
  g_hash_table_insert(bplist_c5,atomnameaskey,atomname);

  for (i=1;i<numatom;++i) {
    for (j=0;j<4;++j)
      named[j]=name_atom_list[i*4+j];
    sea_p_by_name(named,i,bp,ap,apflag,bplist,bplist_c2,bplist_c3,bplist_c4,bplist_c5,numatom,name_atom_list);
  }

  order_bp(bp, bp_f,numatom,numb);

  //  close_ring(numatom,bp_f,numb,name_atom_list);

  flagimso=ch_imso_topo(bp_f,numb,numatom,name_atom_list);

  flagimdo=ch_imdo_topo(bp_f,numb,numatom,name_atom_list);

  //  flagimdo=1;
  if (flagimso==1) aromaflag=OFF;
  if (flagimdo==1) aromaflag=OFF;

  //  ch_PRO_topo(bp_f,numb,numatom,name_atom_list);

  if (aromaflag==ON)
    ch_aromatic_topo(bp_f,numb,numatom,name_atom_list);


  //  g_hash_table_destroy(bplist);
}

int close_ring(int numatom,int **bp_f , int *numb, char *name_atom_list) {
  int i,j,k;
  int flag=OFF;
  int totalnumb=0;
  char named[4],named2[4];

  for (i=0;i<numatom;++i) totalnumb+=numb[i];

  for (i=0;i<numatom;++i) {
    for (j=0;j<4;++j) named[j]=name_atom_list[i*4+j];
    if ((strncmp(named,"C",1)==0 && numb[i]<3) || (strncmp(named,"N",1)==0 && numb[i]<2) ) {
      flag=OFF;
      for (j=i-1;j>=0;--j) {
	for (k=0;k<4;++k) named2[k]=name_atom_list[j*4+k];
	if ((strncmp(named2,"C",1)==0 && numb[j]<3) || (strncmp(named2,"N",1)==0 && numb[j]<2)) {
	  bp_f[i][numb[i]]=j;
	  bp_f[j][numb[j]]=i;
	  ++totalnumb;
	  numb[i]++;
	  numb[j]++;
	  flag=ON;
	  break;
	}
	if (flag==OFF) {
	  for (j=i+1;j<numatom;++j) {
	    for (k=0;k<4;++k) named2[k]=name_atom_list[j*4+k];
	    if ((strncmp(named2,"C",1)==0 && numb[j]<3) || (strncmp(named2,"N",1)==0 && numb[j]<2) || (strncmp(named2,"NE1",3)==0 && numb[j]<3) ) {
	      bp_f[i][numb[i]]=j;
	      bp_f[j][numb[j]]=i;
	      ++totalnumb;
	      numb[i]++;
	      numb[j]++;
	      flag=ON;
	      break;
	    }
	  }
	}
      }
    }
  }

}


void sea_p_by_name(char name[4], int nb, int *bp,int ap[2], int apflag,GHashTable *bplist,GHashTable *bplist_c2,GHashTable *bplist_c3,GHashTable *bplist_c4,GHashTable *bplist_c5,int numatom, char *name_atom_list) {
  int i,j,k;
  int *num_p;
  char *name_p;
  char *name_d,*name_d2,*name_d3;
  char named[4];

  int flag;
  int flag2;  

  gint *atomnum;
  gchar *atomnameaskey;

  name_d=(char *)gcemalloc(sizeof(char)*4);
  name_d2=(char *)gcemalloc(sizeof(char)*2);
  name_d3=(char *)gcemalloc(sizeof(char)*3);

  for (i=0;i<4;++i) name_d[i]=name[i];
  for (i=0;i<2;++i) name_d2[i]=name[i];
  for (i=0;i<3;++i) name_d3[i]=name[i];

  flag=OFF;
  flag2=OFF;
  apflag=OFF;
  for (i=nb-1;i>=0;--i) {
    for (j=0;j<4;++j)
      named[j]=name_atom_list[i*4+j];
    if ((name_p=g_hash_table_lookup(bplist,name_d3))==NULL) {
      printf("There is no key %s in hashtable",name_d3);
      exit(1);
    }
    if (strncmp(named,name_p,3)==0) {
      bp[(nb-1)*2]=i;
      bp[(nb-1)*2+1]=nb;
      flag=ON;
      break;
    }
    if (flag==OFF) {
      if ((name_p=g_hash_table_lookup(bplist_c2,name_d3))!=NULL) {
	if (strncmp(named,name_p,3)==0) {
	  bp[(nb-1)*2]=i;
	  bp[(nb-1)*2+1]=nb;
	  flag=ON;
	  break;
	}
      }
    }
    if (flag==OFF) {
      if ((name_p=g_hash_table_lookup(bplist_c3,name_d))!=NULL) {
	if (strncmp(named,name_p,4)==0) {
	  bp[(nb-1)*2]=i;
	  bp[(nb-1)*2+1]=nb;
	  flag=ON;
	  break;
	}
      }
    }
    if (flag==OFF) {
      if (flag2==OFF) {
	if ((name_p=g_hash_table_lookup(bplist_c4,name_d))!=NULL) {
	  if (strncmp(named,name_p,4)==0) {
	    bp[(nb-1)*2]=i;
	    bp[(nb-1)*2+1]=nb;
	    flag2=ON;
	  }
	}
      }
      if (flag2==ON) {
	for (j=nb-1;j>=0;--j) {
	  if ((name_p=g_hash_table_lookup(bplist_c5,name_d))!=NULL) {
	    if (strncmp(named,name_p,4)==0) {
	      ap[0]=j;
	      ap[1]=nb;
	      flag2=OFF;
	      flag=ON;
	      apflag=ON;
	      break;
	    }
	  }
	}
	if (flag==OFF) {
	  for (j=nb+1;j<numatom;++j) {
	    if ((name_p=g_hash_table_lookup(bplist_c5,name_d))!=NULL) {
	      if (strncmp(named,name_p,4)==0) {
		ap[0]=j;
		ap[1]=nb;
		flag2=OFF;
		flag=ON;
		apflag=ON;
		break;
	      }
	    }
	    if (flag==ON)
	      break;
	  }
	}
      }
    } 
  }
  
  if (flag==OFF) {
    for (i=nb+1;i<numatom;++i) {
      for (j=0;j<4;++j)
	named[j]=name_atom_list[i*4+j];
      if ((name_p=g_hash_table_lookup(bplist,name_d3))==NULL) {
	//	printf("There is no key %s in hashtable",name_d3);
	//	exit(1);
	;
      }
      if (strncmp(named,name_p,3/*2*//*3*/)==0) {
	bp[(nb-1)*2]=nb;
	bp[(nb-1)*2+1]=i;
	flag=ON;
	break;
      }
      if (flag==OFF) {
	if ((name_p=g_hash_table_lookup(bplist_c2,name_d3))==NULL) {
	  //	  printf("There is no key %s in hashtable",name_d2);
	  //  exit(1);
	  ;
	}
	else {
	  if (strncmp(named,name_p,3/*2*//*4*/)==0) {
	    bp[(nb-1)*2]=nb;
	    bp[(nb-1)*2+1]=i;
	    flag=ON;
	    break;
	  }
	}
      }
      if (flag==OFF) {
	if ((name_p=g_hash_table_lookup(bplist_c3,name_d))==NULL) {
	  //	  printf("There is no key %s in hashtable",name_d);
	  //	  exit(1);
	}
	else {
	  if (strncmp(named,name_p,/*2*/4)==0) {
	    bp[(nb-1)*2]=nb;
	    bp[(nb-1)*2+1]=i;
	    flag=ON;
	    break;
	  }
	}
      }
    }
    if (flag==OFF) {
      if (flag2==OFF) {
	if ((name_p=g_hash_table_lookup(bplist_c4,name_d))!=NULL) {
	  if (strncmp(named,name_p,4)==0) {
	    bp[(nb-1)*2]=i;
	    bp[(nb-1)*2+1]=nb;
	    flag2=ON;
	  }
	}
      }
      if (flag2==ON) {
	for (j=nb-1;j>=0;--j) {
	  if ((name_p=g_hash_table_lookup(bplist_c5,name_d))!=NULL) {
	    if (strncmp(named,name_p,4)==0) {
	      ap[(nb-1)*2]=j;
	      ap[(nb-1)*2+1]=nb;
	      flag2=OFF;
	      flag=ON;
	      apflag=ON;
	      break;
	    }
	  }
	}
	if (flag==OFF) {
	  for (i=nb+1;i<numatom;++i) {
	    if ((name_p=g_hash_table_lookup(bplist_c5,name_d))!=NULL) {
	      if (strncmp(named,name_p,4)==0) {
		ap[0]=j;
		ap[1]=nb;
		flag2=OFF;
		flag=ON;
		apflag=ON;
		break;
	      }
	    }
	    if (flag==ON)
	      break;
	  }
	}
      }
    }
  }
}

int order_bp(int *bp, int **bp_f, int numatom, int *numb) {
  int i,j,k;

  for (i=0;i<numatom;++i) numb[i]=0;

  for (i=0;i<numatom;++i) {
    for (j=0;j<numatom-1;++j) {
      if (bp[j*2]==i /*&& bp[j*2] < bp[j*2+1]*/ ) {
	numb[i]+=1;
	bp_f[i]=(int *)gcerealloc(bp_f[i],sizeof(int)*numb[i]);
	bp_f[i][numb[i]-1]=bp[j*2+1];
	/***********************************************************/
        /* if (i<bp_f[i][numb[i]-1])				   */
	/*   printf("%d - %d \n",i+1,bp_f[i][numb[i]-1]+1);	   */
        /***********************************************************/
      }
      else if (bp[j*2+1]==i /*&& bp[j*2+1] < bp[j*2]*/ ) {
	numb[i]+=1;
	bp_f[i]=(int *)gcerealloc(bp_f[i],sizeof(int)*numb[i]);
	bp_f[i][numb[i]-1]=bp[j*2];
	/***********************************************************/
        /* if (i<bp_f[i][numb[i]-1])				   */
	/*   printf("%d - %d \n",i+1,bp_f[i][numb[i]-1]+1);	   */
        /***********************************************************/
      }
    }
  }

  /******************************/
  /* for (i=0;i<numatom;++i) {  */
  /*   bp_f[i][numb[i]-1]='\0'; */
  /* }			        */
  /******************************/

}


/****************************************************************************/
/* int make_bp(char *atom[4],int numatom, int **bpairs, int *numb) {	    */
/*   int i,j;								    */
/* 									    */
/*   for (i=0;i<numatom;++i) {						    */
/*     for (j=1;j<numatom;++j) {					    */
/*       if (i!=j) {							    */
/* 	if (canbeconnect(atom[i],atom[j])==YES) {			    */
/* 	  bpairs[i][numb[i]]=j;						    */
/* 	  numb[i]++;							    */
/* 	  bpairs[numb[i]][numb[j]]=i;					    */
/* 	  numb[j]++;							    */
/* 	  break;							    */
/* 	}								    */
/*       }								    */
/*     }								    */
/*   }									    */
/* }									    */
/* 									    */
/* int canbeconnect(char atomA[4],char atomB[4]) {			    */
/*   int i,j,k;								    */
/* 									    */
/* 									    */
/* 									    */
/* }									    */
/****************************************************************************/


int make_nb_matrix(int **bpairs, int *numb, int depth,
		   int **matrix, int numatom){
  int i,j,k,n;

  for (i=0;i<numatom;++i) for (j=0;j<numatom;++j) matrix[i][j]=-1;

  for (i=0;i<numatom;++i) matrix[i][i]=0;

  for (i=0;i<numatom;++i) {
    for (j=0;j<numb[i];++j) {
      matrix[i][bpairs[i][j]]=1;
      matrix[bpairs[i][j]][i]=1;
    }
  }

  for (n=0;n<depth-1;++n)
    for (i=0;i<numatom;++i)
      for (j=0;j<numatom;++j)
	if (matrix[i][j]==n+1)
	  for (k=0;k<numatom;++k)
	    if (matrix[k][j]==1)
	      if (matrix[i][k]==-1) 
		matrix[i][k]=n+2;

}

int set_nb_pairs(int **matrix, int numatom, 
		 int **pair1_5, int **pair1_4, 
		 int *num1_5, int *num14){
  int i,j;

  for (i=0;i<numatom;++i) {
    num1_5[i]=0;
    num14[i]=0;
  }


  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      if (matrix[i][j]==3) {
	pair1_4[i]=(int *)gcerealloc(pair1_4[i],sizeof(int)*(num14[i]+1));
	pair1_4[i][num14[i]]=j;
	num14[i]=num14[i]+1;
      }
      if (matrix[i][j]==-1) {
	pair1_5[i]=(int *)gcerealloc(pair1_5[i],sizeof(int)*(num1_5[i]+1));
	pair1_5[i][num1_5[i]]=j;
	num1_5[i]=num1_5[i]+1;
      }
    }
  }
}

int ch_aromatic_topo(int **bp_f, int *numb,
		     int numatom,
		     char *name_atom_list){
  int i,j,k,n,m;
  int flag=0;
  char named[4];
  int ncZ,ncE1,ncE2,ncD1,ncD2,ncG,nhE1,nhE2,nhD1,nhD2,nhZ;
  int maxnum1=-1,maxnum2;
  int countcheck[11];

  for (i=0;i<numatom;++i) {
    for (j=0;j<4;++j) named[j]=name_atom_list[i*4+j];
    maxnum2=maxnum1;
    if (strncmp(named,"CZ",2)==0) {
      ncZ=i;
      flag=0;
      for (j=0;j<11;++j) countcheck[j]=0;
      for (j=1;flag<10;++j) {
	n=i+j;
	m=i-j;
	if ( n > maxnum2  ) {
	  for (k=0;k<4;++k) named[k]=name_atom_list[n*4+k];
	  if (strncmp(named,"CE1",3)==0 && countcheck[0]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    ncE1=n;
	    countcheck[0]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"CE2",3)==0 && countcheck[1]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    ncE2=n;
	    countcheck[1]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"CD2",3)==0 && countcheck[2]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    ncD2=n;
	    countcheck[2]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"CD1",3)==0 && countcheck[3]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    ncD1=n;
	    countcheck[3]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"HE2",3)==0 && countcheck[4]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    nhE2=n;
	    countcheck[4]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"HE1",3)==0 && countcheck[5]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    nhE1=n;
	    countcheck[5]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"HD2",3)==0 && countcheck[6]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    nhD2=n;
	    countcheck[6]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"HD1",3)==0 && countcheck[7]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    nhD1=n;
	    countcheck[7]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"CG",2)==0 && countcheck[8]==0) {
	    ncG=n;
	    countcheck[8]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"HZ",2)==0 && countcheck[9]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    nhZ=n;
	    countcheck[9]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"OH",2)==0 && countcheck[10]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    nhZ=n;
	    countcheck[10]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	}
	if ( m > maxnum2 ) {
	  //////////////////////////////////////////////////
	  for (k=0;k<4;++k) named[k]=name_atom_list[m*4+k];
	  if (strncmp(named,"CE1",3)==0 && countcheck[0]==0) {
	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	    numb[m]=0;
	    ncE1=m;
	    countcheck[0]=1;
	    if (maxnum1<m) maxnum1=m;
	    flag++;
	  }
	  else if (strncmp(named,"CE2",3)==0 && countcheck[1]==0) {
	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	    numb[m]=0;
	    ncE2=m;
	    countcheck[1]=1;
	    if (maxnum1<m) maxnum1=m;
	    flag++;
	  }
	  else if (strncmp(named,"CD2",3)==0 && countcheck[2]==0) {
	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	    numb[m]=0;
	    ncD2=m;
	    countcheck[2]=1;
	    if (maxnum1<m) maxnum1=m;
	    flag++;
	  }
	  else if (strncmp(named,"CD1",3)==0 && countcheck[3]==0) {
	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	    numb[m]=0;
	    ncD1=m;
	    countcheck[3]=1;
	    if (maxnum1<m) maxnum1=m;
	    flag++;
	  }
	  else if (strncmp(named,"HE2",3)==0 && countcheck[4]==0) {
	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	    numb[m]=0;
	    nhE2=m;
	    countcheck[4]=1;
	    if (maxnum1<m) maxnum1=m;
	    flag++;
	  }
	  else if (strncmp(named,"HE1",3)==0 && countcheck[5]==0) {
	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	    numb[m]=0;
	    nhE1=m;
	    countcheck[5]=1;
	    if (maxnum1<m) maxnum1=m;
	    flag++;
	  }
	  else if (strncmp(named,"HD2",3)==0 && countcheck[6]==0) {
	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	    numb[m]=0;
	    nhD2=m;
	    countcheck[6]=1;
	    if (maxnum1<m) maxnum1=m;
	    flag++;
	  }
	  else if (strncmp(named,"HD1",3)==0 && countcheck[7]==0) {
	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	    numb[m]=0;
	    nhD1=m;
	    countcheck[7]=1;
	    if (maxnum1<m) maxnum1=m;
	    flag++;
	  }
	  else if (strncmp(named,"CG",2)==0 && countcheck[8]==0) {
	    ncG=m;
	    countcheck[8]=1;
	    if (maxnum1<m) maxnum1=m;
	    flag++;
	  }
	  else if (strncmp(named,"HZ",2)==0 && countcheck[9]==0) {
	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	    numb[m]=0;
	    nhZ=m;
	    countcheck[9]=1;
	    if (maxnum1<m) maxnum1=m;
	    flag++;
	  }
	  else if (strncmp(named,"OH",2)==0 && countcheck[10]==0) {
	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	    numb[m]=0;
	    nhZ=m;
	    countcheck[10]=1;
	    if (maxnum1<m) maxnum1=m;
	    flag++;
	  }
	}
      }
      /////////////////////////////////////////////////
      for (j=0;j<numb[ncZ];++j) {
	if (bp_f[ncZ][j]==ncE1) {
	  for (k=j;k<numb[ncZ]-1;++k) {
	    bp_f[ncZ][k]=bp_f[ncZ][k+1];
	  }
	}
	else if (bp_f[ncZ][j]==ncE1) {
	  for (k=j;k<numb[ncZ]-1;++k) {
	    bp_f[ncZ][k]=bp_f[ncZ][k+1];
	  }
	}
      }
      /////////////////////////////////////////////////
      numb[ncZ]-=2;
      bp_f[ncZ][numb[ncZ]]=ncG;
      numb[ncZ]+=1;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=ncZ;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      bp_f[ncE2][numb[ncE2]]=ncG;
      numb[ncE2]=1;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=ncE2;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      bp_f[ncE1][numb[ncE1]]=ncG;
      numb[ncE1]=1;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=ncE1;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      bp_f[ncD2][numb[ncD2]]=ncG;
      numb[ncD2]=1;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=ncD2;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      bp_f[ncD1][numb[ncD1]]=ncG;
      numb[ncD1]=1;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=ncD1;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      bp_f[nhE2][numb[nhE2]]=ncG;
      numb[nhE2]=1;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=nhE2;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      bp_f[nhE1][numb[nhE1]]=ncG;
      numb[nhE1]=1;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=nhE1;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      bp_f[nhD2][numb[nhD2]]=ncG;
      numb[nhD2]=1;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=nhD2;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      bp_f[nhD1][numb[nhD1]]=ncG;
      numb[nhD1]=1;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=nhD1;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      bp_f[nhZ][numb[nhZ]]=ncG;
      numb[nhZ]=1;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=nhZ;
      numb[ncG]+=1;
    }
  }
}

int ch_imso_topo(int **bp_f, int *numb,
		 int numatom,
		 char *name_atom_list){
  int i,j,k,n,m;
  int flag=0,flag2=0;
  char named[4];
  char name1[4],name2[4],name3[4];
  int ncE1,nnE2,nhE1,nhD1,nhD2,ncG,nnD1,ncD2;
  int maxnum1=-1,maxnum2;
  int countcheck[11];

  for (i=0;i<numatom;++i) {
    for (j=0;j<4;++j) named[j]=name_atom_list[i*4+j];
    maxnum2=maxnum1;
    if (strncmp(named,"CE1",3)==0) {
      for (j=0;j<4;++j) { name1[j]=name_atom_list[(bp_f[i][0])*4+j]; name2[j]=name_atom_list[(bp_f[i][1])*4+j]; name3[j]=name_atom_list[(bp_f[i][2])*4+j];}
      if (    (strncmp(name1,"ND1",3)==0 && strncmp(name2,"NE2",3)==0 && strncmp(name3,"HE1",3)==0)
	    ||(strncmp(name1,"ND1",3)==0 && strncmp(name3,"NE2",3)==0 && strncmp(name2,"HE1",3)==0)
	    ||(strncmp(name2,"ND1",3)==0 && strncmp(name1,"NE2",3)==0 && strncmp(name3,"HE1",3)==0)
	    ||(strncmp(name2,"ND1",3)==0 && strncmp(name3,"NE2",3)==0 && strncmp(name1,"HE1",3)==0)
	    ||(strncmp(name3,"ND1",3)==0 && strncmp(name1,"NE2",3)==0 && strncmp(name2,"HE1",3)==0)
	    ||(strncmp(name3,"ND1",3)==0 && strncmp(name2,"NE2",3)==0 && strncmp(name1,"HE1",3)==0)
	      ) {	
	ncE1=i;
	flag2=1;
	for (j=0;j<7;++j) countcheck[j]=0;
	for (j=1;flag<7;++j) {
	  n=i+j;
	  m=i-j;
	  if ( n > maxnum2  ) {
	    for (k=0;k<4;++k) named[k]=name_atom_list[n*4+k];
	    if (strncmp(named,"NE2",3)==0 && countcheck[0]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      nnE2=n;
	      countcheck[0]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"HE1",3)==0 && countcheck[1]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      nhE1=n;
	      countcheck[1]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"HD1",3)==0 && countcheck[2]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      nhD1=n;
	      countcheck[2]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"HD2",3)==0 && countcheck[3]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      nhD2=n;
	      countcheck[3]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"HE1",3)==0 && countcheck[4]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      nhE1=n;
	      countcheck[4]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"ND1",3)==0 && countcheck[5]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      nnD1=n;
	      countcheck[5]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"CD2",3)==0 && countcheck[6]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      ncD2=n;
	      countcheck[6]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"CG",2)==0 && countcheck[7]==0) {
	      ncG=n;
	      countcheck[7]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	  }
	  if ( m > maxnum2 ) {
	    //////////////////////////////////////////////////
	    for (k=0;k<4;++k) named[k]=name_atom_list[m*4+k];
	    if (strncmp(named,"NE2",3)==0 && countcheck[0]==0) {
	      for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	      numb[m]=0;
	      nnE2=m;
	      countcheck[0]=1;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"HE1",3)==0 && countcheck[1]==0) {
	      for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	      numb[m]=0;
	      nhE1=m;
	      countcheck[1]=1;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"HD1",3)==0 && countcheck[2]==0) {
	      for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	      numb[m]=0;
	      nhD1=m;
	      countcheck[2]=1;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"HD2",3)==0 && countcheck[3]==0) {
	      for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	      numb[m]=0;
	      nhD2=m;
	      countcheck[3]=1;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"HE1",3)==0 && countcheck[4]==0) {
	      for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	      numb[m]=0;
	      nhE1=m;
	      countcheck[4]=1;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"ND1",3)==0 && countcheck[5]==0) {
	      for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	      numb[m]=0;
	      nnD1=m;
	      countcheck[5]=1;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"CD2",3)==0 && countcheck[6]==0) {
	      for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	      numb[m]=0;
	      ncD2=m;
	      countcheck[6]=1;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"CG",2)==0 && countcheck[7]==0) {
	      ncG=m;
	      countcheck[7]=1;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	  }
	}
	/////////////////////////////////////////////////
	numb[ncE1]=1;
	bp_f[ncE1]=(int *)gcerealloc(bp_f[ncE1],sizeof(int)*(numb[ncE1]));
	bp_f[ncE1][0]=ncG;
	bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
	bp_f[ncG][numb[ncG]]=ncE1;
	numb[ncG]+=1;
	/////////////////////////////////////////////////
	numb[nnE2]=1;
	bp_f[nnE2]=(int *)gcerealloc(bp_f[nnE2],sizeof(int)*(numb[nnE2]));
	bp_f[nnE2][0]=ncG;
	bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
	bp_f[ncG][numb[ncG]]=nnE2;
	numb[ncG]+=1;
	/////////////////////////////////////////////////
	numb[nhD1]=1;
	bp_f[nhD1]=(int *)gcerealloc(bp_f[nhD1],sizeof(int)*(numb[nhD1]));
	bp_f[nhD1][0]=ncG;
	bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
	bp_f[ncG][numb[ncG]]=nhD1;
	numb[ncG]+=1;
	/////////////////////////////////////////////////
	numb[nhE1]=1;
	bp_f[nhE1]=(int *)gcerealloc(bp_f[nhE1],sizeof(int)*(numb[nhE1]));
	bp_f[nhE1][0]=ncG;
	bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
	bp_f[ncG][numb[ncG]]=nhE1;
	numb[ncG]+=1;
	/////////////////////////////////////////////////
	numb[nhD2]=1;
	bp_f[nhD2]=(int *)gcerealloc(bp_f[nhD2],sizeof(int)*(numb[nhD2]));
	bp_f[nhD2][0]=ncG;
	bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
	bp_f[ncG][numb[ncG]]=nhD2;
	numb[ncG]+=1;
	/////////////////////////////////////////////////
	numb[nnD1]=1;
	bp_f[nnD1]=(int *)gcerealloc(bp_f[nnD1],sizeof(int)*(numb[nnD1]));
	bp_f[nnD1][0]=ncG;
	bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
	bp_f[ncG][numb[ncG]]=nnD1;
	numb[ncG]+=1;
	/////////////////////////////////////////////////
	numb[ncD2]=1;
	bp_f[ncD2]=(int *)gcerealloc(bp_f[ncD2],sizeof(int)*(numb[ncD2]));
	bp_f[ncD2][0]=ncG;
	bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
	bp_f[ncG][numb[ncG]]=ncD2;
	numb[ncG]+=1;
	/////////////////////////////////////////////////
      }
    }
  }

  return flag2;
}

int ch_imdo_topo(int **bp_f, int *numb,
		 int numatom,
		 char *name_atom_list){
  int i,j,k,n,m;
  int flag=0,flag2=0;
  char named[4];
  char name1[4],name2[4],name3[4];
  int ncG,ncD1,ncD2,ncE2,ncE3,ncZ2,ncZ3,ncH;
  int nnE1;
  int nhD1,nhE1,nhE3,nhZ2,nhZ3,nhH;
  int maxnum1=-1,maxnum2;
  int countcheck[15];

  for (i=0;i<numatom;++i) {
    for (j=0;j<4;++j) named[j]=name_atom_list[i*4+j];
    maxnum2=maxnum1;
    if (strncmp(named,"CH2",3)==0) {
      flag2=1;
      ncH=i;
      flag=1;
      for (j=0;j<15;++j) countcheck[j]=0;
      for (j=1;flag<15;++j) {
	n=i+j;
	m=i-j;
	if ( n > maxnum2  ) {
	  for (k=0;k<4;++k) named[k]=name_atom_list[n*4+k];
	  if (strncmp(named,"CZ3",3)==0 && countcheck[0]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    ncZ3=n;
	    countcheck[0]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"CZ2",3)==0 && countcheck[1]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    ncZ2=n;
	    countcheck[1]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"CE3",3)==0 && countcheck[2]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    ncE3=n;
	    countcheck[2]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"CE2",3)==0 && countcheck[3]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    ncE2=n;
	    countcheck[3]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"CD2",3)==0 && countcheck[4]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    ncD2=n;
	    countcheck[4]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"CD1",3)==0 && countcheck[5]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    ncD1=n;
	    countcheck[5]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"CG",2)==0 && countcheck[6]==0) {
	    ncG=n;
	    countcheck[6]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  if (strncmp(named,"NE1",3)==0 && countcheck[7]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    nnE1=n;
	    countcheck[7]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"HH2",3)==0 && countcheck[8]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    nhH=n;
	    countcheck[8]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"HZ3",3)==0 && countcheck[9]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    nhZ3=n;
	    countcheck[9]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"HZ2",3)==0 && countcheck[10]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    nhZ2=n;
	    countcheck[10]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"HE3",3)==0 && countcheck[11]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    nhE3=n;
	    countcheck[11]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"HE1",3)==0 && countcheck[12]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    nhE1=n;
	    countcheck[12]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"HD1",3)==0 && countcheck[13]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    nhD1=n;
	    countcheck[13]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	}
	if ( m > maxnum2 ) {
	  //////////////////////////////////////////////////
	  for (k=0;k<4;++k) named[k]=name_atom_list[m*4+k];
	  if (strncmp(named,"CZ3",3)==0 && countcheck[0]==0) {
	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	    numb[m]=0;
	    ncZ3=m;
	    countcheck[0]=1;
	    if (maxnum1<m) maxnum1=m;
	    flag++;
	  }
	  else if (strncmp(named,"CZ2",3)==0 && countcheck[1]==0) {
	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	    numb[m]=0;
	    ncZ2=m;
	    countcheck[1]=1;
	    if (maxnum1<m) maxnum1=m;
	    flag++;
	  }
	else if (strncmp(named,"CE3",3)==0 && countcheck[2]==0) {
	  for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	  numb[m]=0;
	  countcheck[2]=1;
	  if (maxnum1<m) maxnum1=m;
	  ncE3=m;
	  flag++;
	}
	else if (strncmp(named,"CE2",3)==0 && countcheck[3]==0) {
	  for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	  numb[m]=0;
	  ncE2=m;
	  countcheck[3]=1;
	  if (maxnum1<m) maxnum1=m;
	  flag++;
	}
	else if (strncmp(named,"CD2",3)==0 && countcheck[4]==0) {
	  for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	  numb[m]=0;
	  ncD2=m;
	  countcheck[4]=1;
	  if (maxnum1<m) maxnum1=m;
	  flag++;
	}
	else if (strncmp(named,"CD1",3)==0 && countcheck[5]==0) {
	  for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	  numb[m]=0;
	  ncD1=m;
	  countcheck[5]=1;
	  if (maxnum1<m) maxnum1=m;
	  flag++;
	}
	else if (strncmp(named,"CG",2)==0 && countcheck[6]==0) {
	  ncG=m;
	  countcheck[6]=1;
	  if (maxnum1<m) maxnum1=m;
	  flag++;
	}
	if (strncmp(named,"NE1",3)==0 && countcheck[7]==0) {
	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	    numb[m]=0;
	    nnE1=m;
	    countcheck[7]=1;
	    if (maxnum1<m) maxnum1=m;
	    flag++;
	}
	else if (strncmp(named,"HH1",3)==0 && countcheck[8]==0) {
	  for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	  numb[m]=0;
	  nhH=m;
	  countcheck[8]=1;
	  if (maxnum1<m) maxnum1=m;
	  flag++;
	}
	else if (strncmp(named,"HZ3",3)==0 && countcheck[9]==0) {
	  for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	  numb[m]=0;
	  nhZ3=m;
	  countcheck[9]=1;
	  if (maxnum1<m) maxnum1=m;
	  flag++;
	}
	else if (strncmp(named,"HZ2",3)==0 && countcheck[10]==0) {
	  for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	  numb[m]=0;
	  nhZ2=m;
	  countcheck[10]=1;
	  if (maxnum1<m) maxnum1=m;
	  flag++;
	}
	else if (strncmp(named,"HE3",3)==0 && countcheck[11]==0) {
	  for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	  numb[m]=0;
	  nhE3=m;
	  countcheck[11]=1;
	  if (maxnum1<m) maxnum1=m;
	  flag++;
	}
	else if (strncmp(named,"HE1",3)==0 && countcheck[12]==0) {
	  for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	  numb[m]=0;
	  nhE1=m;
	  countcheck[12]=1;
	  if (maxnum1<m) maxnum1=m;
	  flag++;
	}
	else if (strncmp(named,"HD1",3)==0 && countcheck[13]==0) {

	  for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	  numb[m]=0;
	  nhD1=m;
	  countcheck[13]=1;
	  if (maxnum1<m) maxnum1=m;
	  flag++;
	}
      }
      /////////////////////////////////////////////////
      numb[ncH]=1;
      bp_f[ncH]=(int *)gcerealloc(bp_f[ncH],sizeof(int)*(numb[ncH]));
      bp_f[ncH][0]=ncG;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=ncH;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      numb[ncZ3]=1;
      bp_f[ncZ3]=(int *)gcerealloc(bp_f[ncZ3],sizeof(int)*(numb[ncZ3]));
      bp_f[ncZ3][0]=ncG;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=ncZ3;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      numb[ncZ2]=1;
      bp_f[ncZ2]=(int *)gcerealloc(bp_f[ncZ2],sizeof(int)*(numb[ncZ2]));
      bp_f[ncZ2][0]=ncG;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=ncZ2;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      numb[ncE3]=1;
      bp_f[ncE3]=(int *)gcerealloc(bp_f[ncE3],sizeof(int)*(numb[ncE3]));
      bp_f[ncE3][0]=ncG;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=ncE3;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      numb[ncE2]=1;
      bp_f[ncE2]=(int *)gcerealloc(bp_f[ncE2],sizeof(int)*(numb[ncE2]));
      bp_f[ncE2][0]=ncG;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=ncE2;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      numb[ncD2]=1;
      bp_f[ncD2]=(int *)gcerealloc(bp_f[ncD2],sizeof(int)*(numb[ncD2]));
      bp_f[ncD2][0]=ncG;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=ncD2;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      numb[ncD1]=1;
      bp_f[ncD1]=(int *)gcerealloc(bp_f[ncD1],sizeof(int)*(numb[ncD1]));
      bp_f[ncD1][0]=ncG;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=ncD1;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      numb[nhD1]=1;
      bp_f[nhD1]=(int *)gcerealloc(bp_f[nhD1],sizeof(int)*(numb[nhD1]));
      bp_f[nhD1][0]=ncG;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=nhD1;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      numb[nhE1]=1;
      bp_f[nhE1]=(int *)gcerealloc(bp_f[nhE1],sizeof(int)*(numb[nhE1]));
      bp_f[nhE1][0]=ncG;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=nhE1;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      numb[nhD1]=1;
      bp_f[nhD1]=(int *)gcerealloc(bp_f[nhD1],sizeof(int)*(numb[nhD1]));
      bp_f[nhD1][0]=ncG;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=nhD1;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      numb[nnE1]=1;
      bp_f[nnE1]=(int *)gcerealloc(bp_f[nnE1],sizeof(int)*(numb[nnE1]));
      bp_f[nnE1][0]=ncG;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=nnE1;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      numb[nhH]=1;
      bp_f[nhH]=(int *)gcerealloc(bp_f[nhH],sizeof(int)*(numb[nhH]));
      bp_f[nhH][0]=ncG;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=nhH;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      numb[nhZ3]=1;
      bp_f[nhZ3]=(int *)gcerealloc(bp_f[nhZ3],sizeof(int)*(numb[nhZ3]));
      bp_f[nhZ3][0]=ncG;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=nhZ3;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      numb[nhZ2]=1;
      bp_f[nhZ2]=(int *)gcerealloc(bp_f[nhZ2],sizeof(int)*(numb[nhZ2]));
      bp_f[nhZ2][0]=ncG;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=nhZ2;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      numb[nhE3]=1;
      bp_f[nhE3]=(int *)gcerealloc(bp_f[nhE3],sizeof(int)*(numb[nhE3]));
      bp_f[nhE3][0]=ncG;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=nhE3;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      numb[nhE1]=1;
      bp_f[nhE1]=(int *)gcerealloc(bp_f[nhE1],sizeof(int)*(numb[nhE1]));
      bp_f[nhE1][0]=ncG;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=nhE1;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      numb[nhD1]=1;
      bp_f[nhD1]=(int *)gcerealloc(bp_f[nhD1],sizeof(int)*(numb[nhD1]));
      bp_f[nhD1][0]=ncG;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=nhD1;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
    }
  }
  return flag2;
}

int ch_PRO_topo(int **bp_f, int *numb,
		int numatom,
		char *name_atom_list){
  int i,j,k,n,m;
  int flag=0,flag2=0;
  char named[4],named2[4];
  char name1[4],name2[4],name3[4];
  int ncA,ncB,ncG,ncD;
  int nhB1,nhB2,nhG1,nhG2,nhD1,nhD2;
  int nn,nnext;
  int maxnum1=-1,maxnum2;
  int countcheck[12];

  for (i=0;i<numatom;++i) {
    for (j=0;j<4;++j) named[j]=name_atom_list[i*4+j];
    maxnum2=maxnum1;
    if (strncmp(named,"CD",2)==0) {
      ncD=i;
      if (numb[ncD]==3) {
	flag2=1;	
	for (j=0;j<12;++j) countcheck[j]=0;

	flag=1;
	for (j=1;flag<11;++j) {
	  n=i+j;
	  m=i-j;
	  if ( n > maxnum2  ) {
	    for (k=0;k<4;++k) named[k]=name_atom_list[n*4+k];
	    if (strncmp(named,"HD 1",4)==0 && countcheck[0]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      nhD1=n;
	      countcheck[0]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"HD 2",4)==0 && countcheck[1]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      nhD2=n;
	      countcheck[1]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"CG",2)==0 && countcheck[2]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      ncG=n;
	      countcheck[2]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"HG 1",4)==0 && countcheck[3]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      nhG1=n;	
	      countcheck[3]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"HG 2",4)==0 && countcheck[4]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      nhG2=n;
	      countcheck[4]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"CB",2)==0 &&  countcheck[5]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      ncB=n;
	      countcheck[5]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"HB 1",4)==0 && countcheck[6]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      nhB1=n;
	      countcheck[6]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"HB 2",4)==0 && countcheck[7]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      nhB2=n;
	      countcheck[7]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"CA",2)==0 && countcheck[8]==0) {
	      /************************************************/
	      /* for (k=0;k<numb[n];++k) bp_f[n][k]=0;	    */
	      /* numb[n]=0;				    */
	      /************************************************/
	      ncA=n;
	      countcheck[8]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"N",1)==0 && numb[n]==2 && countcheck[9]==0) {
	      nn=n;
	      for (k=0;k<4;++k) named2[k]=name_atom_list[(bp_f[nn][0])*4+k];
	      if (strncmp(named2,"CA",2)==0)
		nnext=bp_f[nn][1];
	      else
		nnext=bp_f[nn][0];
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      if (maxnum1<n) maxnum1=n;
	      countcheck[9]=1;
	      flag++;
	    }
	  }
	  if ( m > maxnum2 ) {
	    /////////////////////////////////////////////////
	    for (k=0;k<4;++k) named[k]=name_atom_list[m*4+k];
	    if (strncmp(named,"HD 1",4)==0 && countcheck[0]==0) {
	      for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	      numb[m]=0;
	      countcheck[0]=1;
	      nhD1=m;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"HD 2",4)==0 && countcheck[1]==0) {
	      for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	      numb[m]=0;
	      countcheck[1]=1;
	      nhD2=m;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"CG",2)==0 && countcheck[2]==0) {
	      for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	      numb[m]=0;
	      countcheck[2]=1;
	      ncG=m;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"HG 1",4)==0 && countcheck[3]==0) {
	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	    numb[m]=0;
	    countcheck[3]=1;
	    nhG1=m;
	    if (maxnum1<m) maxnum1=m;
	    flag++;
	    }
	    else if (strncmp(named,"HG 2",4)==0 && countcheck[4]==0) {
	      for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	      numb[m]=0;
	      countcheck[4]=1;
	      nhG2=m;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"CB",2)==0 && countcheck[5]==0 ) {
	      for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	      numb[m]=0;
	      countcheck[5]=1;
	      ncB=m;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"HB 1",4)==0 && countcheck[6]==0 ) {
	      for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	      numb[m]=0;
	      countcheck[6]=1;
	      nhB1=m;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"HB 2",4)==0 && countcheck[7]==0) {
	      for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	      numb[m]=0;
	      countcheck[7]=1;
	      nhB2=m;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"CA",2)==0 && countcheck[8]==0) {
	      /************************************************/
	      /* for (k=0;k<numb[m];++k) bp_f[m][k]=0;	    */
	      /* numb[m]=0;				    */
	      /************************************************/
	      countcheck[8]=1;
	      ncA=m;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"N",1)==0 && numb[m]==2 && countcheck[9]==0) {
	      nn=m;
	      countcheck[9]=1;
	      for (k=0;k<4;++k) named2[k]=name_atom_list[(bp_f[nn][0])*4+k];
	      if (strncmp(named2,"CA",2)==0)
		nnext=bp_f[nn][1];
	      else
		nnext=bp_f[nn][0];
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	  }
	}
	//////////////////////////////////////////////////
	numb[nn]=1;
	bp_f[nn]=(int *)gcerealloc(bp_f[nn],sizeof(int)*(numb[nn]));
	bp_f[nn][0]=ncA;
	/****************************************************************************/
        /* bp_f[ncA]=(int *)gcerealloc(bp_f[ncA],sizeof(int)*(numb[ncA]+1));	    */
	/* bp_f[ncA][numb[ncA]]=nn;						    */
	/* numb[ncA]+=1;							    */
        /****************************************************************************/
	/////////////////////////////////////////////////
	numb[ncD]=1;
	bp_f[ncD]=(int *)gcerealloc(bp_f[ncD],sizeof(int)*(numb[ncD]));
	bp_f[ncD][0]=ncA;
	bp_f[ncA]=(int *)gcerealloc(bp_f[ncA],sizeof(int)*(numb[ncA]+1));
	bp_f[ncA][numb[ncA]]=ncD;
	numb[ncA]+=1;
	/////////////////////////////////////////////////
	numb[ncG]=1;
	bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]));
	bp_f[ncG][0]=ncA;
	bp_f[ncA]=(int *)gcerealloc(bp_f[ncA],sizeof(int)*(numb[ncA]+1));
	bp_f[ncA][numb[ncA]]=ncG;
	numb[ncA]+=1;
	/////////////////////////////////////////////////
	numb[ncB]=1;
	bp_f[ncB]=(int *)gcerealloc(bp_f[ncB],sizeof(int)*(numb[ncB]));
	bp_f[ncB][0]=ncA;
	/****************************************************************************/
        /* bp_f[ncA]=(int *)gcerealloc(bp_f[ncA],sizeof(int)*(numb[ncA]+1));	    */
	/* bp_f[ncA][numb[ncA]]=ncB;						    */
	/* numb[ncA]+=1;							    */
        /****************************************************************************/
	/////////////////////////////////////////////////
	numb[nhD1]=1;
	bp_f[nhD1]=(int *)gcerealloc(bp_f[nhD1],sizeof(int)*(numb[nhD1]));
	bp_f[nhD1][0]=ncA;
	bp_f[ncA]=(int *)gcerealloc(bp_f[ncA],sizeof(int)*(numb[ncA]+1));
	bp_f[ncA][numb[ncA]]=nhD1;
	numb[ncA]+=1;
	/////////////////////////////////////////////////
	numb[nhD2]=1;
	bp_f[nhD2]=(int *)gcerealloc(bp_f[nhD2],sizeof(int)*(numb[nhD2]));
	bp_f[nhD2][0]=ncA;
	bp_f[ncA]=(int *)gcerealloc(bp_f[ncA],sizeof(int)*(numb[ncA]+1));
	bp_f[ncA][numb[ncA]]=nhD2;
	numb[ncA]+=1;
	/////////////////////////////////////////////////
	numb[nhG1]=1;
	bp_f[nhG1]=(int *)gcerealloc(bp_f[nhG1],sizeof(int)*(numb[nhG1]));
	bp_f[nhG1][0]=ncA;
	bp_f[ncA]=(int *)gcerealloc(bp_f[ncA],sizeof(int)*(numb[ncA]+1));
	bp_f[ncA][numb[ncA]]=nhG1;
	numb[ncA]+=1;
	/////////////////////////////////////////////////
	numb[nhG2]=1;
	bp_f[nhG2]=(int *)gcerealloc(bp_f[nhG2],sizeof(int)*(numb[nhG2]));
	bp_f[nhG2][0]=ncA;
	bp_f[ncA]=(int *)gcerealloc(bp_f[ncA],sizeof(int)*(numb[ncA]+1));
	bp_f[ncA][numb[ncA]]=nhG2;
	numb[ncA]+=1;
	/////////////////////////////////////////////////
	numb[nhB1]=1;
	bp_f[nhB1]=(int *)gcerealloc(bp_f[nhB1],sizeof(int)*(numb[nhB1]));
	bp_f[nhB1][0]=ncA;
	bp_f[ncA]=(int *)gcerealloc(bp_f[ncA],sizeof(int)*(numb[ncA]+1));
	bp_f[ncA][numb[ncA]]=nhB1;
	numb[ncA]+=1;
	/////////////////////////////////////////////////
	numb[nhB2]=1;
	bp_f[nhB2]=(int *)gcerealloc(bp_f[nhB2],sizeof(int)*(numb[nhB2]));
	bp_f[nhB2][0]=ncA;
	bp_f[ncA]=(int *)gcerealloc(bp_f[ncA],sizeof(int)*(numb[ncA]+1));
	bp_f[ncA][numb[ncA]]=nhB2;
	numb[ncA]+=1;
	/////////////////////////////////////////////////
	for (j=0;j<numb[nnext];++j) {
	  for (k=0;k<4;++k) named[k]=name_atom_list[bp_f[nnext][j]*4+k];
	  if (strncmp(named,"N",1)==0) {
	    bp_f[nnext][j]=ncA;
	    break;
	  }
	}
	bp_f[ncA]=(int *)gcerealloc(bp_f[ncA],sizeof(int)*(numb[ncA]+1));
	bp_f[ncA][numb[ncA]]=nnext;
	numb[ncA]+=1;
      }
    }
  }
}



















	      apflag=ON;
	      break;
	    }
	  }
	}
	if (flag==OFF) {
	  for (i=nb+1;i<numatom;++i) {
	    if ((name_p=g_hash_table_lookup(bplist_c5,name_d))!=NULL) {
	      if (strncmp(named,name_p,4)==0) {
		ap[0]=j;
		ap[1]=nb;
		flag2=OFF;
		flag=ON;
		apflag=ON;
		break;
	      }
	    }
	    if (flag==ON)
	      break;
	  }
	}
      }
    }
  }
}

int order_bp(int *bp, int **bp_f, int numatom, int *numb) {
  int i,j,k;

  for (i=0;i<numatom;++i) numb[i]=0;

  for (i=0;i<numatom;++i) {
    for (j=0;j<numatom-1;++j) {
      if (bp[j*2]==i /*&& bp[j*2] < bp[j*2+1]*/ ) {
	numb[i]+=1;
	bp_f[i]=(int *)gcerealloc(bp_f[i],sizeof(int)*numb[i]);
	bp_f[i][numb[i]-1]=bp[j*2+1];
	/***********************************************************/
        /* if (i<bp_f[i][numb[i]-1])				   */
	/*   printf("%d - %d \n",i+1,bp_f[i][numb[i]-1]+1);	   */
        /***********************************************************/
      }
      else if (bp[j*2+1]==i /*&& bp[j*2+1] < bp[j*2]*/ ) {
	numb[i]+=1;
	bp_f[i]=(int *)gcerealloc(bp_f[i],sizeof(int)*numb[i]);
	bp_f[i][numb[i]-1]=bp[j*2];
	/***********************************************************/
        /* if (i<bp_f[i][numb[i]-1])				   */
	/*   printf("%d - %d \n",i+1,bp_f[i][numb[i]-1]+1);	   */
        /***********************************************************/
      }
    }
  }

  /******************************/
  /* for (i=0;i<numatom;++i) {  */
  /*   bp_f[i][numb[i]-1]='\0'; */
  /* }			        */
  /******************************/

}


/****************************************************************************/
/* int make_bp(char *atom[4],int numatom, int **bpairs, int *numb) {	    */
/*   int i,j;								    */
/* 									    */
/*   for (i=0;i<numatom;++i) {						    */
/*     for (j=1;j<numatom;++j) {					    */
/*       if (i!=j) {							    */
/* 	if (canbeconnect(atom[i],atom[j])==YES) {			    */
/* 	  bpairs[i][numb[i]]=j;						    */
/* 	  numb[i]++;							    */
/* 	  bpairs[numb[i]][numb[j]]=i;					    */
/* 	  numb[j]++;							    */
/* 	  break;							    */
/* 	}								    */
/*       }								    */
/*     }								    */
/*   }									    */
/* }									    */
/* 									    */
/* int canbeconnect(char atomA[4],char atomB[4]) {			    */
/*   int i,j,k;								    */
/* 									    */
/* 									    */
/* 									    */
/* }									    */
/****************************************************************************/


int make_nb_matrix(int **bpairs, int *numb, int depth,
		   int **matrix, int numatom){
  int i,j,k,n;

  for (i=0;i<numatom;++i) for (j=0;j<numatom;++j) matrix[i][j]=-1;

  for (i=0;i<numatom;++i) matrix[i][i]=0;

  for (i=0;i<numatom;++i) {
    for (j=0;j<numb[i];++j) {
      matrix[i][bpairs[i][j]]=1;
      matrix[bpairs[i][j]][i]=1;
    }
  }

  for (n=0;n<depth-1;++n)
    for (i=0;i<numatom;++i)
      for (j=0;j<numatom;++j)
	if (matrix[i][j]==n+1)
	  for (k=0;k<numatom;++k)
	    if (matrix[k][j]==1)
	      if (matrix[i][k]==-1) 
		matrix[i][k]=n+2;

}

int set_nb_pairs(int **matrix, int numatom, 
		 int **pair1_5, int **pair1_4, 
		 int *num1_5, int *num14){
  int i,j;

  for (i=0;i<numatom;++i) {
    num1_5[i]=0;
    num14[i]=0;
  }


  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      if (matrix[i][j]==3) {
	pair1_4[i]=(int *)gcerealloc(pair1_4[i],sizeof(int)*(num14[i]+1));
	pair1_4[i][num14[i]]=j;
	num14[i]++;
      }
      if (matrix[i][j]==-1) {
	pair1_5[i]=(int *)gcerealloc(pair1_5[i],sizeof(int)*(num1_5[i]+1));
	pair1_5[i][num1_5[i]]=j;
	num1_5[i]++;
      }
    }
  }
}

int ch_aromatic_topo(int **bp_f, int *numb,
		     int numatom,
		     char *name_atom_list){
  int i,j,k,n,m;
  int flag=0;
  char named[4];
  int ncZ,ncE1,ncE2,ncD1,ncD2,ncG,nhE1,nhE2,nhD1,nhD2,nhZ;
  int maxnum1=-1,maxnum2;
  int countcheck[11];

  for (i=0;i<numatom;++i) {
    for (j=0;j<4;++j) named[j]=name_atom_list[i*4+j];
    maxnum2=maxnum1;
    if (strncmp(named,"CZ",2)==0) {
      ncZ=i;
      flag=0;
      for (j=0;j<11;++j) countcheck[j]=0;
      for (j=1;flag<10;++j) {
	n=i+j;
	m=i-j;
	if ( n > maxnum2  ) {
	  for (k=0;k<4;++k) named[k]=name_atom_list[n*4+k];
	  if (strncmp(named,"CE1",3)==0 && countcheck[0]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    ncE1=n;
	    countcheck[0]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"CE2",3)==0 && countcheck[1]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    ncE2=n;
	    countcheck[1]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"CD2",3)==0 && countcheck[2]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    ncD2=n;
	    countcheck[2]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"CD1",3)==0 && countcheck[3]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    ncD1=n;
	    countcheck[3]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"HE2",3)==0 && countcheck[4]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    nhE2=n;
	    countcheck[4]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"HE1",3)==0 && countcheck[5]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    nhE1=n;
	    countcheck[5]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"HD2",3)==0 && countcheck[6]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    nhD2=n;
	    countcheck[6]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"HD1",3)==0 && countcheck[7]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    nhD1=n;
	    countcheck[7]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"CG",2)==0 && countcheck[8]==0) {
	    ncG=n;
	    countcheck[8]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"HZ",2)==0 && countcheck[9]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    nhZ=n;
	    countcheck[9]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	  else if (strncmp(named,"OH",2)==0 && countcheck[10]==0) {
	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	    numb[n]=0;
	    nhZ=n;
	    countcheck[10]=1;
	    if (maxnum1<n) maxnum1=n;
	    flag++;
	  }
	}
	if ( m > maxnum2 ) {
	  //////////////////////////////////////////////////
	  for (k=0;k<4;++k) named[k]=name_atom_list[m*4+k];
	  if (strncmp(named,"CE1",3)==0 && countcheck[0]==0) {
	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	    numb[m]=0;
	    ncE1=m;
	    countcheck[0]=1;
	    if (maxnum1<m) maxnum1=m;
	    flag++;
	  }
	  else if (strncmp(named,"CE2",3)==0 && countcheck[1]==0) {
	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	    numb[m]=0;
	    ncE2=m;
	    countcheck[1]=1;
	    if (maxnum1<m) maxnum1=m;
	    flag++;
	  }
	  else if (strncmp(named,"CD2",3)==0 && countcheck[2]==0) {
	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	    numb[m]=0;
	    ncD2=m;
	    countcheck[2]=1;
	    if (maxnum1<m) maxnum1=m;
	    flag++;
	  }
	  else if (strncmp(named,"CD1",3)==0 && countcheck[3]==0) {
	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	    numb[m]=0;
	    ncD1=m;
	    countcheck[3]=1;
	    if (maxnum1<m) maxnum1=m;
	    flag++;
	  }
	  else if (strncmp(named,"HE2",3)==0 && countcheck[4]==0) {
	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	    numb[m]=0;
	    nhE2=m;
	    countcheck[4]=1;
	    if (maxnum1<m) maxnum1=m;
	    flag++;
	  }
	  else if (strncmp(named,"HE1",3)==0 && countcheck[5]==0) {
	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	    numb[m]=0;
	    nhE1=m;
	    countcheck[5]=1;
	    if (maxnum1<m) maxnum1=m;
	    flag++;
	  }
	  else if (strncmp(named,"HD2",3)==0 && countcheck[6]==0) {
	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	    numb[m]=0;
	    nhD2=m;
	    countcheck[6]=1;
	    if (maxnum1<m) maxnum1=m;
	    flag++;
	  }
	  else if (strncmp(named,"HD1",3)==0 && countcheck[7]==0) {
	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	    numb[m]=0;
	    nhD1=m;
	    countcheck[7]=1;
	    if (maxnum1<m) maxnum1=m;
	    flag++;
	  }
	  else if (strncmp(named,"CG",2)==0 && countcheck[8]==0) {
	    ncG=m;
	    countcheck[8]=1;
	    if (maxnum1<m) maxnum1=m;
	    flag++;
	  }
	  else if (strncmp(named,"HZ",2)==0 && countcheck[9]==0) {
	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	    numb[m]=0;
	    nhZ=m;
	    countcheck[9]=1;
	    if (maxnum1<m) maxnum1=m;
	    flag++;
	  }
	  else if (strncmp(named,"OH",2)==0 && countcheck[10]==0) {
	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	    numb[m]=0;
	    nhZ=m;
	    countcheck[10]=1;
	    if (maxnum1<m) maxnum1=m;
	    flag++;
	  }
	}
      }
      /////////////////////////////////////////////////
      for (j=0;j<numb[ncZ];++j) {
	if (bp_f[ncZ][j]==ncE1) {
	  for (k=j;k<numb[ncZ]-1;++k) {
	    bp_f[ncZ][k]=bp_f[ncZ][k+1];
	  }
	}
	else if (bp_f[ncZ][j]==ncE1) {
	  for (k=j;k<numb[ncZ]-1;++k) {
	    bp_f[ncZ][k]=bp_f[ncZ][k+1];
	  }
	}
      }
      /////////////////////////////////////////////////
      numb[ncZ]-=2;
      bp_f[ncZ][numb[ncZ]]=ncG;
      numb[ncZ]+=1;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=ncZ;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      bp_f[ncE2][numb[ncE2]]=ncG;
      numb[ncE2]=1;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=ncE2;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      bp_f[ncE1][numb[ncE1]]=ncG;
      numb[ncE1]=1;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=ncE1;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      bp_f[ncD2][numb[ncD2]]=ncG;
      numb[ncD2]=1;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=ncD2;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      bp_f[ncD1][numb[ncD1]]=ncG;
      numb[ncD1]=1;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=ncD1;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      bp_f[nhE2][numb[nhE2]]=ncG;
      numb[nhE2]=1;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=nhE2;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      bp_f[nhE1][numb[nhE1]]=ncG;
      numb[nhE1]=1;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=nhE1;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      bp_f[nhD2][numb[nhD2]]=ncG;
      numb[nhD2]=1;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=nhD2;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      bp_f[nhD1][numb[nhD1]]=ncG;
      numb[nhD1]=1;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=nhD1;
      numb[ncG]+=1;
      /////////////////////////////////////////////////
      bp_f[nhZ][numb[nhZ]]=ncG;
      numb[nhZ]=1;
      bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
      bp_f[ncG][numb[ncG]]=nhZ;
      numb[ncG]+=1;
    }
  }
}

int ch_imso_topo(int **bp_f, int *numb,
		 int numatom,
		 char *name_atom_list){
  int i,j,k,n,m;
  int flag=0,flag2=0;
  char named[4];
  char name1[4],name2[4],name3[4];
  int ncE1,nnE2,nhE1,nhD1,nhD2,ncG,nnD1,ncD2;
  int maxnum1=-1,maxnum2;
  int countcheck[7];

  for (i=0;i<numatom;++i) {
    for (j=0;j<4;++j) named[j]=name_atom_list[i*4+j];
    maxnum2=maxnum1;
    if (strncmp(named,"CE1",3)==0) {
      for (j=0;j<4;++j) { name1[j]=name_atom_list[(bp_f[i][0])*4+j]; name2[j]=name_atom_list[(bp_f[i][1])*4+j]; name3[j]=name_atom_list[(bp_f[i][2])*4+j];}
      if (    (strncmp(name1,"ND1",3)==0 && strncmp(name2,"NE2",3)==0 && strncmp(name3,"HE1",3)==0)
	    ||(strncmp(name1,"ND1",3)==0 && strncmp(name3,"NE2",3)==0 && strncmp(name2,"HE1",3)==0)
	    ||(strncmp(name2,"ND1",3)==0 && strncmp(name1,"NE2",3)==0 && strncmp(name3,"HE1",3)==0)
	    ||(strncmp(name2,"ND1",3)==0 && strncmp(name3,"NE2",3)==0 && strncmp(name1,"HE1",3)==0)
	    ||(strncmp(name3,"ND1",3)==0 && strncmp(name1,"NE2",3)==0 && strncmp(name2,"HE1",3)==0)
	    ||(strncmp(name3,"ND1",3)==0 && strncmp(name2,"NE2",3)==0 && strncmp(name1,"HE1",3)==0)
	      ) {	
	ncE1=i;
	flag2=1;
	for (j=0;j<7;++j) countcheck[j]=0;
	for (j=1;flag<7;++j) {
	  n=i+j;
	  m=i-j;
	  if ( n > maxnum2  ) {
	    for (k=0;k<4;++k) named[k]=name_atom_list[n*4+k];
	    if (strncmp(named,"NE2",3)==0 && countcheck[0]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      nnE2=n;
	      countcheck[0]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"HE1",3)==0 && countcheck[1]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      nhE1=n;
	      countcheck[1]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"HD1",3)==0 && countcheck[2]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      nhD1=n;
	      countcheck[2]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"HD2",3)==0 && countcheck[3]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      nhD2=n;
	      countcheck[3]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"HE1",3)==0 && countcheck[4]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      nhE1=n;
	      countcheck[4]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"ND1",3)==0 && countcheck[5]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      nnD1=n;
	      countcheck[5]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"CD2",3)==0 && countcheck[6]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      ncD2=n;
	      countcheck[6]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"CG",2)==0 && countcheck[7]==0) {
	      ncG=n;
	      countcheck[7]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	  }
	  if ( m > maxnum2 ) {
	    //////////////////////////////////////////////////
	    for (k=0;k<4;++k) named[k]=name_atom_list[m*4+k];
	    if (strncmp(named,"NE2",3)==0 && countcheck[0]==0) {
	      for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	      numb[m]=0;
	      nnE2=m;
	      countcheck[0]=1;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"HE1",3)==0 && countcheck[1]==0) {
	      for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	      numb[m]=0;
	      nhE1=m;
	      countcheck[1]=1;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"HD1",3)==0 && countcheck[2]==0) {
	      for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	      numb[m]=0;
	      nhD1=m;
	      countcheck[2]=1;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"HD2",3)==0 && countcheck[3]==0) {
	      for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	      numb[m]=0;
	      nhD2=m;
	      countcheck[3]=1;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"HE1",3)==0 && countcheck[4]==0) {
	      for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	      numb[m]=0;
	      nhE1=m;
	      countcheck[4]=1;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"ND1",3)==0 && countcheck[5]==0) {
	      for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	      numb[m]=0;
	      nnD1=m;
	      countcheck[5]=1;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"CD2",3)==0 && countcheck[6]==0) {
	      for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	      numb[m]=0;
	      ncD2=m;
	      countcheck[6]=1;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"CG",2)==0 && countcheck[7]==0) {
	      ncG=m;
	      countcheck[7]=1;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	  }
	}
	/////////////////////////////////////////////////
	numb[ncE1]=1;
	bp_f[ncE1]=(int *)gcerealloc(bp_f[ncE1],sizeof(int)*(numb[ncE1]));
	bp_f[ncE1][0]=ncG;
	bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
	bp_f[ncG][numb[ncG]]=ncE1;
	numb[ncG]+=1;
	/////////////////////////////////////////////////
	numb[nnE2]=1;
	bp_f[nnE2]=(int *)gcerealloc(bp_f[nnE2],sizeof(int)*(numb[nnE2]));
	bp_f[nnE2][0]=ncG;
	bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
	bp_f[ncG][numb[ncG]]=nnE2;
	numb[ncG]+=1;
	/////////////////////////////////////////////////
	numb[nhD1]=1;
	bp_f[nhD1]=(int *)gcerealloc(bp_f[nhD1],sizeof(int)*(numb[nhD1]));
	bp_f[nhD1][0]=ncG;
	bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
	bp_f[ncG][numb[ncG]]=nhD1;
	numb[ncG]+=1;
	/////////////////////////////////////////////////
	numb[nhE1]=1;
	bp_f[nhE1]=(int *)gcerealloc(bp_f[nhE1],sizeof(int)*(numb[nhE1]));
	bp_f[nhE1][0]=ncG;
	bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
	bp_f[ncG][numb[ncG]]=nhE1;
	numb[ncG]+=1;
	/////////////////////////////////////////////////
	numb[nhD2]=1;
	bp_f[nhD2]=(int *)gcerealloc(bp_f[nhD2],sizeof(int)*(numb[nhD2]));
	bp_f[nhD2][0]=ncG;
	bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
	bp_f[ncG][numb[ncG]]=nhD2;
	numb[ncG]+=1;
	/////////////////////////////////////////////////
	numb[nnD1]=1;
	bp_f[nnD1]=(int *)gcerealloc(bp_f[nnD1],sizeof(int)*(numb[nnD1]));
	bp_f[nnD1][0]=ncG;
	bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
	bp_f[ncG][numb[ncG]]=nnD1;
	numb[ncG]+=1;
	/////////////////////////////////////////////////
	numb[ncD2]=1;
	bp_f[ncD2]=(int *)gcerealloc(bp_f[ncD2],sizeof(int)*(numb[ncD2]));
	bp_f[ncD2][0]=ncG;
	bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));
	bp_f[ncG][numb[ncG]]=ncD2;
	numb[ncG]+=1;
	/////////////////////////////////////////////////
      }
    }
  }

  return flag2;
}

/***********************************************************************************/
/* int ch_imdo_topo(int **bp_f, int *numb,					   */
/* 		 int numatom,							   */
/* 		 char *name_atom_list){						   */
/*   int i,j,k,n,m;								   */
/*   int flag=0,flag2=0;							   */
/*   char named[4];								   */
/*   char name1[4],name2[4],name3[4];						   */
/*   int ncG,ncD1,ncD2,ncE2,ncE3,ncZ2,ncZ3,ncH;					   */
/*   int nnE1;									   */
/*   int nhD1,nhE1,nhE3,nhZ2,nhZ3,nhH;						   */
/*   int maxnum1=-1,maxnum2;							   */
/*   int countcheck[15];							   */
/* 										   */
/*   for (i=0;i<numatom;++i) {							   */
/*     for (j=0;j<4;++j) named[j]=name_atom_list[i*4+j];			   */
/*     maxnum2=maxnum1;								   */
/*     if (strncmp(named,"CH2",3)==0) {						   */
/*       flag2=1;								   */
/*       ncH=i;									   */
/*       flag=1;								   */
/*       for (j=0;j<15;++j) countcheck[j]=0;					   */
/*       for (j=1;flag<15;++j) {						   */
/* 	n=i+j;									   */
/* 	m=i-j;									   */
/* 	if ( n > maxnum2  ) {							   */
/* 	  for (k=0;k<4;++k) named[k]=name_atom_list[n*4+k];			   */
/* 	  if (strncmp(named,"CZ3",3)==0 && countcheck[0]==0) {			   */
/* 	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;				   */
/* 	    numb[n]=0;								   */
/* 	    ncZ3=n;								   */
/* 	    countcheck[0]=1;							   */
/* 	    if (maxnum1<n) maxnum1=n;						   */
/* 	    flag++;								   */
/* 	  }									   */
/* 	  else if (strncmp(named,"CZ2",3)==0 && countcheck[1]==0) {		   */
/* 	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;				   */
/* 	    numb[n]=0;								   */
/* 	    ncZ2=n;								   */
/* 	    countcheck[1]=1;							   */
/* 	    if (maxnum1<n) maxnum1=n;						   */
/* 	    flag++;								   */
/* 	  }									   */
/* 	  else if (strncmp(named,"CE3",3)==0 && countcheck[2]==0) {		   */
/* 	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;				   */
/* 	    numb[n]=0;								   */
/* 	    ncE3=n;								   */
/* 	    countcheck[2]=1;							   */
/* 	    if (maxnum1<n) maxnum1=n;						   */
/* 	    flag++;								   */
/* 	  }									   */
/* 	  else if (strncmp(named,"CE2",3)==0 && countcheck[3]==0) {		   */
/* 	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;				   */
/* 	    numb[n]=0;								   */
/* 	    ncE2=n;								   */
/* 	    countcheck[3]=1;							   */
/* 	    if (maxnum1<n) maxnum1=n;						   */
/* 	    flag++;								   */
/* 	  }									   */
/* 	  else if (strncmp(named,"CD2",3)==0 && countcheck[4]==0) {		   */
/* 	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;				   */
/* 	    numb[n]=0;								   */
/* 	    ncD2=n;								   */
/* 	    countcheck[4]=1;							   */
/* 	    if (maxnum1<n) maxnum1=n;						   */
/* 	    flag++;								   */
/* 	  }									   */
/* 	  else if (strncmp(named,"CD1",3)==0 && countcheck[5]==0) {		   */
/* 	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;				   */
/* 	    numb[n]=0;								   */
/* 	    ncD1=n;								   */
/* 	    countcheck[5]=1;							   */
/* 	    if (maxnum1<n) maxnum1=n;						   */
/* 	    flag++;								   */
/* 	  }									   */
/* 	  else if (strncmp(named,"CG",2)==0 && countcheck[6]==0) {		   */
/* 	    ncG=n;								   */
/* 	    countcheck[6]=1;							   */
/* 	    if (maxnum1<n) maxnum1=n;						   */
/* 	    flag++;								   */
/* 	  }									   */
/* 	  if (strncmp(named,"NE1",3)==0 && countcheck[7]==0) {			   */
/* 	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;				   */
/* 	    numb[n]=0;								   */
/* 	    nnE1=n;								   */
/* 	    countcheck[7]=1;							   */
/* 	    if (maxnum1<n) maxnum1=n;						   */
/* 	    flag++;								   */
/* 	  }									   */
/* 	  else if (strncmp(named,"HH2",3)==0 && countcheck[8]==0) {		   */
/* 	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;				   */
/* 	    numb[n]=0;								   */
/* 	    nhH=n;								   */
/* 	    countcheck[8]=1;							   */
/* 	    if (maxnum1<n) maxnum1=n;						   */
/* 	    flag++;								   */
/* 	  }									   */
/* 	  else if (strncmp(named,"HZ3",3)==0 && countcheck[9]==0) {		   */
/* 	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;				   */
/* 	    numb[n]=0;								   */
/* 	    nhZ3=n;								   */
/* 	    countcheck[9]=1;							   */
/* 	    if (maxnum1<n) maxnum1=n;						   */
/* 	    flag++;								   */
/* 	  }									   */
/* 	  else if (strncmp(named,"HZ2",3)==0 && countcheck[10]==0) {		   */
/* 	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;				   */
/* 	    numb[n]=0;								   */
/* 	    nhZ2=n;								   */
/* 	    countcheck[10]=1;							   */
/* 	    if (maxnum1<n) maxnum1=n;						   */
/* 	    flag++;								   */
/* 	  }									   */
/* 	  else if (strncmp(named,"HE3",3)==0 && countcheck[11]==0) {		   */
/* 	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;				   */
/* 	    numb[n]=0;								   */
/* 	    nhE3=n;								   */
/* 	    countcheck[11]=1;							   */
/* 	    if (maxnum1<n) maxnum1=n;						   */
/* 	    flag++;								   */
/* 	  }									   */
/* 	  else if (strncmp(named,"HE1",3)==0 && countcheck[12]==0) {		   */
/* 	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;				   */
/* 	    numb[n]=0;								   */
/* 	    nhE1=n;								   */
/* 	    countcheck[12]=1;							   */
/* 	    if (maxnum1<n) maxnum1=n;						   */
/* 	    flag++;								   */
/* 	  }									   */
/* 	  else if (strncmp(named,"HD1",3)==0 && countcheck[13]==0) {		   */
/* 	    for (k=0;k<numb[n];++k) bp_f[n][k]=0;				   */
/* 	    numb[n]=0;								   */
/* 	    nhD1=n;								   */
/* 	    countcheck[13]=1;							   */
/* 	    if (maxnum1<n) maxnum1=n;						   */
/* 	    flag++;								   */
/* 	  }									   */
/* 	}									   */
/* 	if ( m > maxnum2 ) {							   */
/* 	  //////////////////////////////////////////////////			   */
/* 	  for (k=0;k<4;++k) named[k]=name_atom_list[m*4+k];			   */
/* 	  if (strncmp(named,"CZ3",3)==0 && countcheck[0]==0) {			   */
/* 	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;				   */
/* 	    numb[m]=0;								   */
/* 	    ncZ3=m;								   */
/* 	    countcheck[0]=1;							   */
/* 	    if (maxnum1<m) maxnum1=m;						   */
/* 	    flag++;								   */
/* 	  }									   */
/* 	  else if (strncmp(named,"CZ2",3)==0 && countcheck[1]==0) {		   */
/* 	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;				   */
/* 	    numb[m]=0;								   */
/* 	    ncZ2=m;								   */
/* 	    countcheck[1]=1;							   */
/* 	    if (maxnum1<m) maxnum1=m;						   */
/* 	    flag++;								   */
/* 	  }									   */
/* 	  else if (strncmp(named,"CE3",3)==0 && countcheck[2]==0) {		   */
/* 	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;				   */
/* 	    numb[m]=0;								   */
/* 	    countcheck[2]=1;							   */
/* 	    if (maxnum1<m) maxnum1=m;						   */
/* 	    ncE3=m;								   */
/* 	    flag++;								   */
/* 	  }									   */
/* 	  else if (strncmp(named,"CE2",3)==0 && countcheck[3]==0) {		   */
/* 	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;				   */
/* 	    numb[m]=0;								   */
/* 	    ncE2=m;								   */
/* 	    countcheck[3]=1;							   */
/* 	    if (maxnum1<m) maxnum1=m;						   */
/* 	    flag++;								   */
/* 	  }									   */
/* 	  else if (strncmp(named,"CD2",3)==0 && countcheck[4]==0) {		   */
/* 	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;				   */
/* 	    numb[m]=0;								   */
/* 	    ncD2=m;								   */
/* 	    countcheck[4]=1;							   */
/* 	    if (maxnum1<m) maxnum1=m;						   */
/* 	    flag++;								   */
/* 	  }									   */
/* 	  else if (strncmp(named,"CD1",3)==0 && countcheck[5]==0) {		   */
/* 	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;				   */
/* 	    numb[m]=0;								   */
/* 	    ncD1=m;								   */
/* 	    countcheck[5]=1;							   */
/* 	    if (maxnum1<m) maxnum1=m;						   */
/* 	    flag++;								   */
/* 	  }									   */
/* 	  else if (strncmp(named,"CG",2)==0 && countcheck[6]==0) {		   */
/* 	    ncG=m;								   */
/* 	    countcheck[6]=1;							   */
/* 	    if (maxnum1<m) maxnum1=m;						   */
/* 	    flag++;								   */
/* 	  }									   */
/* 	  if (strncmp(named,"NE1",3)==0 && countcheck[7]==0) {			   */
/* 	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;				   */
/* 	    numb[m]=0;								   */
/* 	    nnE1=m;								   */
/* 	    countcheck[7]=1;							   */
/* 	    if (maxnum1<m) maxnum1=m;						   */
/* 	    flag++;								   */
/* 	  }									   */
/* 	  else if (strncmp(named,"HH1",3)==0 && countcheck[8]==0) {		   */
/* 	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;				   */
/* 	    numb[m]=0;								   */
/* 	    nhH=m;								   */
/* 	    countcheck[8]=1;							   */
/* 	    if (maxnum1<m) maxnum1=m;						   */
/* 	    flag++;								   */
/* 	  }									   */
/* 	  else if (strncmp(named,"HZ3",3)==0 && countcheck[9]==0) {		   */
/* 	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;				   */
/* 	    numb[m]=0;								   */
/* 	    nhZ3=m;								   */
/* 	    countcheck[9]=1;							   */
/* 	    if (maxnum1<m) maxnum1=m;						   */
/* 	    flag++;								   */
/* 	  }									   */
/* 	  else if (strncmp(named,"HZ2",3)==0 && countcheck[10]==0) {		   */
/* 	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;				   */
/* 	    numb[m]=0;								   */
/* 	    nhZ2=m;								   */
/* 	    countcheck[10]=1;							   */
/* 	    if (maxnum1<m) maxnum1=m;						   */
/* 	    flag++;								   */
/* 	  }									   */
/* 	  else if (strncmp(named,"HE3",3)==0 && countcheck[11]==0) {		   */
/* 	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;				   */
/* 	    numb[m]=0;								   */
/* 	    nhE3=m;								   */
/* 	    countcheck[11]=1;							   */
/* 	    if (maxnum1<m) maxnum1=m;						   */
/* 	    flag++;								   */
/* 	  }									   */
/* 	  else if (strncmp(named,"HE1",3)==0 && countcheck[12]==0) {		   */
/* 	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;				   */
/* 	    numb[m]=0;								   */
/* 	    nhE1=m;								   */
/* 	    countcheck[12]=1;							   */
/* 	    if (maxnum1<m) maxnum1=m;						   */
/* 	    flag++;								   */
/* 	  }									   */
/* 	  else if (strncmp(named,"HD1",3)==0 && countcheck[13]==0) {		   */
/* 	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;				   */
/* 	    numb[m]=0;								   */
/* 	    nhD1=m;								   */
/* 	    countcheck[13]=1;							   */
/* 	    if (maxnum1<m) maxnum1=m;						   */
/* 	    flag++;								   */
/* 	  }									   */
/* 	}									   */
/*       }									   */
/*       /////////////////////////////////////////////////			   */
/*       numb[ncH]=1;								   */
/*       bp_f[ncH]=(int *)gcerealloc(bp_f[ncH],sizeof(int)*(numb[ncH]));	   */
/*       bp_f[ncH][0]=ncG;							   */
/*       bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));	   */
/*       bp_f[ncG][numb[ncG]]=ncH;						   */
/*       numb[ncG]+=1;								   */
/*       /////////////////////////////////////////////////			   */
/*       numb[ncZ3]=1;								   */
/*       bp_f[ncZ3]=(int *)gcerealloc(bp_f[ncZ3],sizeof(int)*(numb[ncZ3]));	   */
/*       bp_f[ncZ3][0]=ncG;							   */
/*       bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));	   */
/*       bp_f[ncG][numb[ncG]]=ncZ3;						   */
/*       numb[ncG]+=1;								   */
/*       /////////////////////////////////////////////////			   */
/*       numb[ncZ2]=1;								   */
/*       bp_f[ncZ2]=(int *)gcerealloc(bp_f[ncZ2],sizeof(int)*(numb[ncZ2]));	   */
/*       bp_f[ncZ2][0]=ncG;							   */
/*       bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));	   */
/*       bp_f[ncG][numb[ncG]]=ncZ2;						   */
/*       numb[ncG]+=1;								   */
/*       /////////////////////////////////////////////////			   */
/*       numb[ncE3]=1;								   */
/*       bp_f[ncE3]=(int *)gcerealloc(bp_f[ncE3],sizeof(int)*(numb[ncE3]));	   */
/*       bp_f[ncE3][0]=ncG;							   */
/*       bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));	   */
/*       bp_f[ncG][numb[ncG]]=ncE3;						   */
/*       numb[ncG]+=1;								   */
/*       /////////////////////////////////////////////////			   */
/*       numb[ncE2]=1;								   */
/*       bp_f[ncE2]=(int *)gcerealloc(bp_f[ncE2],sizeof(int)*(numb[ncE2]));	   */
/*       bp_f[ncE2][0]=ncG;							   */
/*       bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));	   */
/*       bp_f[ncG][numb[ncG]]=ncE2;						   */
/*       numb[ncG]+=1;								   */
/*       /////////////////////////////////////////////////			   */
/*       numb[ncD2]=1;								   */
/*       bp_f[ncD2]=(int *)gcerealloc(bp_f[ncD2],sizeof(int)*(numb[ncD2]));	   */
/*       bp_f[ncD2][0]=ncG;							   */
/*       bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));	   */
/*       bp_f[ncG][numb[ncG]]=ncD2;						   */
/*       numb[ncG]+=1;								   */
/*       /////////////////////////////////////////////////			   */
/*       numb[ncD1]=1;								   */
/*       bp_f[ncD1]=(int *)gcerealloc(bp_f[ncD1],sizeof(int)*(numb[ncD1]));	   */
/*       bp_f[ncD1][0]=ncG;							   */
/*       bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));	   */
/*       bp_f[ncG][numb[ncG]]=ncD1;						   */
/*       numb[ncG]+=1;								   */
/*       /////////////////////////////////////////////////			   */
/*       numb[nhD1]=1;								   */
/*       bp_f[nhD1]=(int *)gcerealloc(bp_f[nhD1],sizeof(int)*(numb[nhD1]));	   */
/*       bp_f[nhD1][0]=ncG;							   */
/*       bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));	   */
/*       bp_f[ncG][numb[ncG]]=nhD1;						   */
/*       numb[ncG]+=1;								   */
/*       /////////////////////////////////////////////////			   */
/*       numb[nhE1]=1;								   */
/*       bp_f[nhE1]=(int *)gcerealloc(bp_f[nhE1],sizeof(int)*(numb[nhE1]));	   */
/*       bp_f[nhE1][0]=ncG;							   */
/*       bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));	   */
/*       bp_f[ncG][numb[ncG]]=nhE1;						   */
/*       numb[ncG]+=1;								   */
/*       /////////////////////////////////////////////////			   */
/*       numb[nhD1]=1;								   */
/*       bp_f[nhD1]=(int *)gcerealloc(bp_f[nhD1],sizeof(int)*(numb[nhD1]));	   */
/*       bp_f[nhD1][0]=ncG;							   */
/*       bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));	   */
/*       bp_f[ncG][numb[ncG]]=nhD1;						   */
/*       numb[ncG]+=1;								   */
/*       /////////////////////////////////////////////////			   */
/*       numb[nnE1]=1;								   */
/*       bp_f[nnE1]=(int *)gcerealloc(bp_f[nnE1],sizeof(int)*(numb[nnE1]));	   */
/*       bp_f[nnE1][0]=ncG;							   */
/*       bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));	   */
/*       bp_f[ncG][numb[ncG]]=nnE1;						   */
/*       numb[ncG]+=1;								   */
/*       /////////////////////////////////////////////////			   */
/*       numb[nhH]=1;								   */
/*       bp_f[nhH]=(int *)gcerealloc(bp_f[nhH],sizeof(int)*(numb[nhH]));	   */
/*       bp_f[nhH][0]=ncG;							   */
/*       bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));	   */
/*       bp_f[ncG][numb[ncG]]=nhH;						   */
/*       numb[ncG]+=1;								   */
/*       /////////////////////////////////////////////////			   */
/*       numb[nhZ3]=1;								   */
/*       bp_f[nhZ3]=(int *)gcerealloc(bp_f[nhZ3],sizeof(int)*(numb[nhZ3]));	   */
/*       bp_f[nhZ3][0]=ncG;							   */
/*       bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));	   */
/*       bp_f[ncG][numb[ncG]]=nhZ3;						   */
/*       numb[ncG]+=1;								   */
/*       /////////////////////////////////////////////////			   */
/*       numb[nhZ2]=1;								   */
/*       bp_f[nhZ2]=(int *)gcerealloc(bp_f[nhZ2],sizeof(int)*(numb[nhZ2]));	   */
/*       bp_f[nhZ2][0]=ncG;							   */
/*       bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));	   */
/*       bp_f[ncG][numb[ncG]]=nhZ2;						   */
/*       numb[ncG]+=1;								   */
/*       /////////////////////////////////////////////////			   */
/*       numb[nhE3]=1;								   */
/*       bp_f[nhE3]=(int *)gcerealloc(bp_f[nhE3],sizeof(int)*(numb[nhE3]));	   */
/*       bp_f[nhE3][0]=ncG;							   */
/*       bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));	   */
/*       bp_f[ncG][numb[ncG]]=nhE3;						   */
/*       numb[ncG]+=1;								   */
/*       /////////////////////////////////////////////////			   */
/*       numb[nhE1]=1;								   */
/*       bp_f[nhE1]=(int *)gcerealloc(bp_f[nhE1],sizeof(int)*(numb[nhE1]));	   */
/*       bp_f[nhE1][0]=ncG;							   */
/*       bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));	   */
/*       bp_f[ncG][numb[ncG]]=nhE1;						   */
/*       numb[ncG]+=1;								   */
/*       /////////////////////////////////////////////////			   */
/*       numb[nhD1]=1;								   */
/*       bp_f[nhD1]=(int *)gcerealloc(bp_f[nhD1],sizeof(int)*(numb[nhD1]));	   */
/*       bp_f[nhD1][0]=ncG;							   */
/*       bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]+1));	   */
/*       bp_f[ncG][numb[ncG]]=nhD1;						   */
/*       numb[ncG]+=1;								   */
/*       /////////////////////////////////////////////////			   */
/*     }									   */
/*   }										   */
/*   return flag2;								   */
/* }										   */
/***********************************************************************************/
/*  */


int ch_PRO_topo(int **bp_f, int *numb,
		int numatom,
		char *name_atom_list){
  int i,j,k,n,m;
  int flag=0,flag2=0;
  char named[4],named2[4];
  char name1[4],name2[4],name3[4];
  int ncA,ncB,ncG,ncD;
  int nhB1,nhB2,nhG1,nhG2,nhD1,nhD2;
  int nn,nnext;
  int maxnum1=-1,maxnum2;
  int countcheck[12];

  for (i=0;i<numatom;++i) {
    for (j=0;j<4;++j) named[j]=name_atom_list[i*4+j];
    maxnum2=maxnum1;
    if (strncmp(named,"CD",2)==0) {
      ncD=i;
      if (numb[ncD]==3) {
	flag2=1;	
	for (j=0;j<12;++j) countcheck[j]=0;

	flag=1;
	for (j=1;flag<11;++j) {
	  n=i+j;
	  m=i-j;
	  if ( n > maxnum2  ) {
	    for (k=0;k<4;++k) named[k]=name_atom_list[n*4+k];
	    if (strncmp(named,"HD 1",4)==0 && countcheck[0]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      nhD1=n;
	      countcheck[0]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"HD 2",4)==0 && countcheck[1]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      nhD2=n;
	      countcheck[1]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"CG",2)==0 && countcheck[2]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      ncG=n;
	      countcheck[2]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"HG 1",4)==0 && countcheck[3]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      nhG1=n;	
	      countcheck[3]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"HG 2",4)==0 && countcheck[4]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      nhG2=n;
	      countcheck[4]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"CB",2)==0 &&  countcheck[5]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      ncB=n;
	      countcheck[5]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"HB 1",4)==0 && countcheck[6]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      nhB1=n;
	      countcheck[6]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"HB 2",4)==0 && countcheck[7]==0) {
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      nhB2=n;
	      countcheck[7]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"CA",2)==0 && countcheck[8]==0) {
	      /************************************************/
	      /* for (k=0;k<numb[n];++k) bp_f[n][k]=0;	    */
	      /* numb[n]=0;				    */
	      /************************************************/
	      ncA=n;
	      countcheck[8]=1;
	      if (maxnum1<n) maxnum1=n;
	      flag++;
	    }
	    else if (strncmp(named,"N",1)==0 && numb[n]==2 && countcheck[9]==0) {
	      nn=n;
	      for (k=0;k<4;++k) named2[k]=name_atom_list[(bp_f[nn][0])*4+k];
	      if (strncmp(named2,"CA",2)==0)
		nnext=bp_f[nn][1];
	      else
		nnext=bp_f[nn][0];
	      for (k=0;k<numb[n];++k) bp_f[n][k]=0;
	      numb[n]=0;
	      if (maxnum1<n) maxnum1=n;
	      countcheck[9]=1;
	      flag++;
	    }
	  }
	  if ( m > maxnum2 ) {
	    /////////////////////////////////////////////////
	    for (k=0;k<4;++k) named[k]=name_atom_list[m*4+k];
	    if (strncmp(named,"HD 1",4)==0 && countcheck[0]==0) {
	      for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	      numb[m]=0;
	      countcheck[0]=1;
	      nhD1=m;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"HD 2",4)==0 && countcheck[1]==0) {
	      for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	      numb[m]=0;
	      countcheck[1]=1;
	      nhD2=m;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"CG",2)==0 && countcheck[2]==0) {
	      for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	      numb[m]=0;
	      countcheck[2]=1;
	      ncG=m;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"HG 1",4)==0 && countcheck[3]==0) {
	    for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	    numb[m]=0;
	    countcheck[3]=1;
	    nhG1=m;
	    if (maxnum1<m) maxnum1=m;
	    flag++;
	    }
	    else if (strncmp(named,"HG 2",4)==0 && countcheck[4]==0) {
	      for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	      numb[m]=0;
	      countcheck[4]=1;
	      nhG2=m;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"CB",2)==0 && countcheck[5]==0 ) {
	      for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	      numb[m]=0;
	      countcheck[5]=1;
	      ncB=m;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"HB 1",4)==0 && countcheck[6]==0 ) {
	      for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	      numb[m]=0;
	      countcheck[6]=1;
	      nhB1=m;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"HB 2",4)==0 && countcheck[7]==0) {
	      for (k=0;k<numb[m];++k) bp_f[m][k]=0;
	      numb[m]=0;
	      countcheck[7]=1;
	      nhB2=m;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"CA",2)==0 && countcheck[8]==0) {
	      /************************************************/
	      /* for (k=0;k<numb[m];++k) bp_f[m][k]=0;	    */
	      /* numb[m]=0;				    */
	      /************************************************/
	      countcheck[8]=1;
	      ncA=m;
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	    else if (strncmp(named,"N",1)==0 && numb[m]==2 && countcheck[9]==0) {
	      nn=m;
	      countcheck[9]=1;
	      for (k=0;k<4;++k) named2[k]=name_atom_list[(bp_f[nn][0])*4+k];
	      if (strncmp(named2,"CA",2)==0)
		nnext=bp_f[nn][1];
	      else
		nnext=bp_f[nn][0];
	      if (maxnum1<m) maxnum1=m;
	      flag++;
	    }
	  }
	}
	//////////////////////////////////////////////////
	numb[nn]=1;
	bp_f[nn]=(int *)gcerealloc(bp_f[nn],sizeof(int)*(numb[nn]));
	bp_f[nn][0]=ncA;
	/****************************************************************************/
        /* bp_f[ncA]=(int *)gcerealloc(bp_f[ncA],sizeof(int)*(numb[ncA]+1));	    */
	/* bp_f[ncA][numb[ncA]]=nn;						    */
	/* numb[ncA]+=1;							    */
        /****************************************************************************/
	/////////////////////////////////////////////////
	numb[ncD]=1;
	bp_f[ncD]=(int *)gcerealloc(bp_f[ncD],sizeof(int)*(numb[ncD]));
	bp_f[ncD][0]=ncA;
	bp_f[ncA]=(int *)gcerealloc(bp_f[ncA],sizeof(int)*(numb[ncA]+1));
	bp_f[ncA][numb[ncA]]=ncD;
	numb[ncA]+=1;
	/////////////////////////////////////////////////
	numb[ncG]=1;
	bp_f[ncG]=(int *)gcerealloc(bp_f[ncG],sizeof(int)*(numb[ncG]));
	bp_f[ncG][0]=ncA;
	bp_f[ncA]=(int *)gcerealloc(bp_f[ncA],sizeof(int)*(numb[ncA]+1));
	bp_f[ncA][numb[ncA]]=ncG;
	numb[ncA]+=1;
	/////////////////////////////////////////////////
	numb[ncB]=1;
	bp_f[ncB]=(int *)gcerealloc(bp_f[ncB],sizeof(int)*(numb[ncB]));
	bp_f[ncB][0]=ncA;
	/****************************************************************************/
        /* bp_f[ncA]=(int *)gcerealloc(bp_f[ncA],sizeof(int)*(numb[ncA]+1));	    */
	/* bp_f[ncA][numb[ncA]]=ncB;						    */
	/* numb[ncA]+=1;							    */
        /****************************************************************************/
	/////////////////////////////////////////////////
	numb[nhD1]=1;
	bp_f[nhD1]=(int *)gcerealloc(bp_f[nhD1],sizeof(int)*(numb[nhD1]));
	bp_f[nhD1][0]=ncA;
	bp_f[ncA]=(int *)gcerealloc(bp_f[ncA],sizeof(int)*(numb[ncA]+1));
	bp_f[ncA][numb[ncA]]=nhD1;
	numb[ncA]+=1;
	/////////////////////////////////////////////////
	numb[nhD2]=1;
	bp_f[nhD2]=(int *)gcerealloc(bp_f[nhD2],sizeof(int)*(numb[nhD2]));
	bp_f[nhD2][0]=ncA;
	bp_f[ncA]=(int *)gcerealloc(bp_f[ncA],sizeof(int)*(numb[ncA]+1));
	bp_f[ncA][numb[ncA]]=nhD2;
	numb[ncA]+=1;
	/////////////////////////////////////////////////
	numb[nhG1]=1;
	bp_f[nhG1]=(int *)gcerealloc(bp_f[nhG1],sizeof(int)*(numb[nhG1]));
	bp_f[nhG1][0]=ncA;
	bp_f[ncA]=(int *)gcerealloc(bp_f[ncA],sizeof(int)*(numb[ncA]+1));
	bp_f[ncA][numb[ncA]]=nhG1;
	numb[ncA]+=1;
	/////////////////////////////////////////////////
	numb[nhG2]=1;
	bp_f[nhG2]=(int *)gcerealloc(bp_f[nhG2],sizeof(int)*(numb[nhG2]));
	bp_f[nhG2][0]=ncA;
	bp_f[ncA]=(int *)gcerealloc(bp_f[ncA],sizeof(int)*(numb[ncA]+1));
	bp_f[ncA][numb[ncA]]=nhG2;
	numb[ncA]+=1;
	/////////////////////////////////////////////////
	numb[nhB1]=1;
	bp_f[nhB1]=(int *)gcerealloc(bp_f[nhB1],sizeof(int)*(numb[nhB1]));
	bp_f[nhB1][0]=ncA;
	bp_f[ncA]=(int *)gcerealloc(bp_f[ncA],sizeof(int)*(numb[ncA]+1));
	bp_f[ncA][numb[ncA]]=nhB1;
	numb[ncA]+=1;
	/////////////////////////////////////////////////
	numb[nhB2]=1;
	bp_f[nhB2]=(int *)gcerealloc(bp_f[nhB2],sizeof(int)*(numb[nhB2]));
	bp_f[nhB2][0]=ncA;
	bp_f[ncA]=(int *)gcerealloc(bp_f[ncA],sizeof(int)*(numb[ncA]+1));
	bp_f[ncA][numb[ncA]]=nhB2;
	numb[ncA]+=1;
	/////////////////////////////////////////////////
	for (j=0;j<numb[nnext];++j) {
	  for (k=0;k<4;++k) named[k]=name_atom_list[bp_f[nnext][j]*4+k];
	  if (strncmp(named,"N",1)==0) {
	    bp_f[nnext][j]=ncA;
	    break;
	  }
	}
	bp_f[ncA]=(int *)gcerealloc(bp_f[ncA],sizeof(int)*(numb[ncA]+1));
	bp_f[ncA][numb[ncA]]=nnext;
	numb[ncA]+=1;
      }
    }
  }
}
