#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "PROTOPO.h"

#define fconst 4.184070*100.0
//#define fconst 1.0

#include "glib.h"

//#include "ABA.h"        // 2015-02-16
//#include "gener.h"      // 2015-02-16
//#include "MD.h"         // 2015-02-16
//#include "force.h"      // 2015-02-16
//#include "physics.h"    // 2015-02-16
//#include "math.h"       // 2015-02-16

//#include "TOPO.h"      // 2015-02-16
#include "TOPO.h"      // 2015-02-16
#include "EF.h"

#include "ECEPE.h"

// void print_iter(gpointer key, gpointer val, gpointer usr_data);

#define ON 1
#define OFF 0

char *atomtypelist[300]= {
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
  "CE1","CD ",
  "CE2","CD1",
  "CZ ","CE ",
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
  "HZ ","CZ ",
  "HN ","N  ",
  "HH ","OH ",
  "HH1","NH1",
  "HH2","NH2",
  "H3 ","C  ",
  "N  ","C  ",
  "NZ ","CE ",
  "NE ","CD ",
  "NE2","CD ",
  "ND ","CG ",
  "ND1","CG ",
  "ND2","CG ",
  "NH1","CZ ",
  "NH2","CZ ",
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

char *atomtypelist_c2[100]= {
  "HG ","OG ",
  "C  ","CA ",
  "CD1","CG1",
  "CE ","SD ",
  "CE1","CD1",
  "CZ ","NEX",
  "HG1","OG1",
  "HD2","ND2",
  "HE2","NE2",
  "HZ ","NZ "
};

char *atomtypelist_c3[numatomtype*2]= {
  "C   ","H3 1",
  "HD2 ","OD2 ",
  "HE2 ","OE2 ",
  "CZ  ","NEX "
};

char *atomtypelist_c4[numatomtype*2]= {
  "CZ  ","CE1 "
};

char *atomtypelist_c5[numatomtype*2]= {
  "CZ  ","CE2 "
};

void read_ECEPE_parm(char *preofilename, char *bd8filename, struct ECEPE_parms *ECEPE_p, struct pnb *nb_p) {
  int i,j;
  int d;
  double f;
  char *line,dummy;
  size_t len=0;

  FILE *preo,*bd8;

  struct ECEPE_atom *atom_dummy;

  preo=efopen(preofilename,"r");
  bd8=efopen(bd8filename,"r");

  getline(&line,&len,preo);

  fscanf(preo,"%d",&(*ECEPE_p).NUMATM);
  fscanf(preo,"%d",&(*ECEPE_p).NUMVAR);
  fscanf(preo,"%d",&(*ECEPE_p).NUMRES);
  fscanf(preo,"%d",&(*ECEPE_p).NUMINT);
  fscanf(preo,"%d",&(*ECEPE_p).NUMS);

  (*ECEPE_p).dihed=(struct ECEPE_dihe *)calloc((*ECEPE_p).NUMVAR,sizeof(struct ECEPE_dihe)); // 2015-02-16
  (*ECEPE_p).atom=(struct ECEPE_atom *)calloc((*ECEPE_p).NUMATM,sizeof(struct ECEPE_atom));  // 2015-02-16
  atom_dummy=(struct ECEPE_atom *)calloc((*ECEPE_p).NUMATM,sizeof(struct ECEPE_atom));       // 2015-02-16

  for (i=0;i<(*ECEPE_p).NUMVAR;++i) {
     fscanf(preo,"%lf",&(*ECEPE_p).dihed[i].angle);
     if((*ECEPE_p).dihed[i].angle==-75.000) 
       i-=1;
   }
   /*****************************/
   /* getline(&line,&len,preo); */
   /* getline(&line,&len,preo); */
   /* getline(&line,&len,preo); */
   /*****************************/

  for (i=0;i<5;++i)
    fscanf(preo,"%d",&d);

  for (i=0;i<(*ECEPE_p).NUMVAR;++i) {
    fscanf(preo,"%d",&(*ECEPE_p).dihed[i].indexv1);
    fscanf(preo,"%d",&(*ECEPE_p).dihed[i].indexv2);
    fscanf(preo,"%d",&(*ECEPE_p).dihed[i].ibnd1);
    fscanf(preo,"%d",&(*ECEPE_p).dihed[i].ibnd2);
    fscanf(preo,"%d",&(*ECEPE_p).dihed[i].ifront);
    fscanf(preo,"%d",&(*ECEPE_p).dihed[i].ibchar1);
    fscanf(preo,"%d",&(*ECEPE_p).dihed[i].ibchar2);
    fscanf(preo,"%d",&(*ECEPE_p).dihed[i].ibchar3);
    fscanf(preo,"%lf",&(*ECEPE_p).dihed[i].A);
    fscanf(preo,"%d",&(*ECEPE_p).dihed[i].NB);
    fscanf(preo,"%d",&(*ECEPE_p).dihed[i].NS);
    fscanf(preo,"%d",&(*ECEPE_p).dihed[i].IFTOR);
    fscanf(preo,"%d",&d);
    fscanf(preo,"%d",&d);
    fscanf(preo,"%d",&d);
    fscanf(preo,"%d",&d);
    fscanf(preo,"%d",&d);
  }

  for (i=0;i<(*ECEPE_p).NUMATM;++i) {
    for (j=0;j<3;++j)
      fscanf(preo,"%lf",&atom_dummy[i].refcoord[j]);
    fscanf(preo,"%lf",&atom_dummy[i].charge);
    fscanf(preo,"%d",&atom_dummy[i].nbtype);
    fscanf(preo,"%d",&atom_dummy[i].kunit);
    fscanf(preo,"%d",&atom_dummy[i].katom);
    fscanf(preo,"%d",&atom_dummy[i].jatom);
    for (j=0;j<2;++j)
      getc(preo);
    for (j=0;j<4;++j)
      atom_dummy[i].name_atom[j]=getc(preo);
    fscanf(preo,"%3s",&atom_dummy[i].name_res);
    fscanf(preo,"%d",&d);
  }

  for (i=0;i<(*ECEPE_p).NUMATM;++i) {
    for (j=0;j<3;++j)
      (*ECEPE_p).atom[atom_dummy[i].katom-1].refcoord[j]=atom_dummy[atom_dummy[i].jatom-1].refcoord[j];
    (*ECEPE_p).atom[atom_dummy[i].katom-1].charge=atom_dummy[atom_dummy[i].jatom-1].charge;
    (*ECEPE_p).atom[atom_dummy[i].katom-1].nbtype=atom_dummy[atom_dummy[i].jatom-1].nbtype;
    (*ECEPE_p).atom[atom_dummy[i].katom-1].kunit=atom_dummy[atom_dummy[i].jatom-1].kunit;
    for (j=0;j<4;++j)
      (*ECEPE_p).atom[atom_dummy[i].katom-1].name_atom[j]=atom_dummy[atom_dummy[i].jatom-1].name_atom[j];
    for (j=0;j<4;++j)
      (*ECEPE_p).atom[atom_dummy[i].katom-1].name_res[j]=atom_dummy[atom_dummy[i].jatom-1].name_res[j];
    (*ECEPE_p).atom[atom_dummy[i].katom-1].katom=atom_dummy[atom_dummy[i].katom-1].katom;
    (*ECEPE_p).atom[atom_dummy[i].jatom-1].jatom=atom_dummy[atom_dummy[i].jatom-1].jatom;
  }

  //  (*nb_p).Acff=(double *)gcemalloc(sizeof(double)*numatomtype*numatomtype); // 2015-02-16
  //  (*nb_p).Bcff=(double *)gcemalloc(sizeof(double)*numatomtype*numatomtype); // 2015-02-16
  (*nb_p).Acff=(double *)calloc(numatomtype*numatomtype,sizeof(double)); // 2015-02-16
  (*nb_p).Bcff=(double *)calloc(numatomtype*numatomtype,sizeof(double)); // 2015-02-16

  for (i=0;i<numatomtype;++i) {
    for (j=0;j<numatomtype;++j) {
      fscanf(bd8,"%lf",&(*nb_p).Acff[i*numatomtype+j]);
    }
    for (j=0;j<numatomtype;++j) {
      fscanf(bd8,"%lf",&(*nb_p).Bcff[i*numatomtype+j]);
    }
  }

  fclose(preo);
  fclose(bd8);

  //  free((*ECEPE_p).dihed); // 2015-02-16
  //  free((*ECEPE_p).atom);  // 2015-02-16
  free(atom_dummy);       // 2015-02-16
}

void read_ECEPE_detail_coo (FILE *file, double *co, int numatom) {
  int i,j,na;
  char *line,x[20],y[20],z[20];
  size_t len=0;
  
  na=0;
  for (i=0;na<numatom;++i) {
    getline(&line,&len,file);
    if (i>18) {
      for (j=0;j<14;++j) {
	x[j]=line[0+j];
	y[j]=line[14+j];
	z[j]=line[28+j];
	co[na*3]=atof(x);
	co[na*3+1]=atof(y);
	co[na*3+2]=atof(z);
      }
      ++na;
    }
  }
}

int read_ECEPE_detail_coo_cyc (FILE *file, double *co, double *dihed,double *ene) {
  int i,c;
  double d;
  
  for (i=0;i<6;++i) {
    c= fscanf(file,"%lf",&ene[i]);
    if(c==-1)
      return -1;
  }
  for (i=0;i<40*3;++i) {
    c= fscanf(file,"%lf",&d);
    if(c!=-1)
      co[i]=d;
    else
      return -1;
  }
  for (i=0;i<10;++i) {
    c= fscanf(file,"%lf",&d);
    if(c!=-1)
      dihed[i]=d;
    else
      return -1;
  }

  return 1;
}

int read_ECEPE_detail_coo_cyc_for_protein(FILE *file, double *co, double *dihed,double *ene) {
  int i,c;
  double d;
  
  for (i=0;i<6;++i) {
    c= fscanf(file,"%lf",&ene[i]);
    if(c==-1)
      return -1;
  }
  for (i=0;i<400*3;++i) {
    c= fscanf(file,"%lf",&d);
    if(c!=-1)
      co[i]=d;
    else
      return -1;
  }
  for (i=0;i<100;++i) {
    c= fscanf(file,"%lf",&d);
    if(c!=-1)
      dihed[i]=d;
    else
      return -1;
  }

  return 1;
}


double calc_ff_ECEPE(struct ECEPE_pote *p, struct ECEPE_force *f, struct ECEPE_parms ECEPE_p, struct pnb nb_p, int *pairs, int numint){
  /****************************************************************************/
  /* int i,j,k;								      */
  /* int numdih,numatom;						      */
  /* double *A,*dih;							      */
  /* int *ns,*nb,*nbtype;						      */
  /* int *iftors;							      */
  /* double *charge;							      */
  /* double p_nb,*f_nb,p_es,*f_es,t,*n,*co;				      */
  /* 									      */
  /* numdih=ECEPE_p.NUMVAR;						      */
  /* numatom=ECEPE_p.NUMATM;						      */
  /* 									      */
  /* co=(double *)gcemalloc(sizeof(double)*numatom*3);			      */
  /* A=(double *)gcemalloc(sizeof(double)*numdih);			      */
  /* ns=(int *)gcemalloc(sizeof(int)*numdih);				      */
  /* nb=(int *)gcemalloc(sizeof(int)*numdih);				      */
  /* dih=(double *)gcemalloc(sizeof(double)*numdih);			      */
  /* iftors=(int *)gcemalloc(sizeof(int)*numdih);			      */
  /* charge=(double *)gcemalloc(sizeof(double)*numatom);		      */
  /* nbtype=(int *)gcemalloc(sizeof(int)*numatom);			      */
  /* 									      */
  /* f_nb=(double *)gcemalloc(sizeof(double)*numatom*3);		      */
  /* f_es=(double *)gcemalloc(sizeof(double)*numatom*3);		      */
  /* n=(double *)gcemalloc(sizeof(double)*numdih);			      */
  /* 									      */
  /* for (i=0;i<numdih;++i) {						      */
  /*   A[i]=ECEPE_p.dihed[i].A;						      */
  /*   ns[i]=ECEPE_p.dihed[i].NS;					      */
  /*   nb[i]=ECEPE_p.dihed[i].NB;					      */
  /*   dih[i]=ECEPE_p.dihed[i].angle;					      */
  /*   iftors[i]=ECEPE_p.dihed[i].IFTOR;				      */
  /* }									      */
  /* 									      */
  /* for (i=0;i<numatom;++i) {						      */
  /*   charge[i]=ECEPE_p.atom[i].charge;				      */
  /*   nbtype[i]=ECEPE_p.atom[i].nbtype;				      */
  /*   for (j=0;j<3;++j)						      */
  /*     co[i*3+j]=ECEPE_p.atom[i].refcoord[j];				      */
  /* }									      */
  /* 									      */
  /* //  set_pairs_ECEPE();						      */
  /* //  calc_TORS_ECEPE(&t,n,A,ns,nb,dih,iftors,numdih);		      */
  /* calc_TORS_ECEPE2(co,ECEPE_p, delta_dihed);				      */
  /* calc_NB_ECEPE(&p_nb,f_nb,&p_es,f_es,co,pairs,numint,charge,nbtype,nb_p); */
  /* 									      */
  /* (*p).p_t=t+p_nb+p_es;						      */
  /* (*p).p_tors=t;							      */
  /* (*p).p_nb=p_nb;							      */
  /* (*p).p_es=p_es;							      */
  /****************************************************************************/

}

double calc_ff_ECEPE_for_db(struct ECEPE_pote *p, struct ECEPE_force *f, struct ECEPE_parms ECEPE_p, struct pnb nb_p, int **pairs_1_5,int **pairs_1_4, int *num1_5, int *num1_4){
  int i,j,k;
  int numdih,numatom;
  double *A,*dih;
  int *ns,*nb,*nbtype;
  int *iftors;
  double *charge;
  double p_nb,*f_nb,p_es,*f_es,t,*n,*co;

  numdih=ECEPE_p.NUMVAR;
  numatom=ECEPE_p.NUMATM;

  co=(double *)calloc(numatom*3,sizeof(double));        // 2015-02-16
  A=(double *)calloc(numdih,sizeof(double));            // 2015-02-16
  ns=(int *)calloc(numdih,sizeof(int));                 // 2015-02-16
  nb=(int *)calloc(numdih,sizeof(int));                 // 2015-02-16
  dih=(double *)calloc(numdih,sizeof(double));          // 2015-02-16
  iftors=(int *)calloc(numdih,sizeof(int));             // 2015-02-16
  charge=(double *)calloc(numatom,sizeof(double));      // 2015-02-16
  nbtype=(int *)calloc(numatom,sizeof(int));            // 2015-02-16

  f_nb=(double *)calloc(numatom*3,sizeof(double));      // 2015-02-16
  f_es=(double *)calloc(numatom*3,sizeof(double));      // 2015-02-16
  n=(double *)calloc(numdih,sizeof(double));            // 2015-02-16

  for (i=0;i<numdih;++i) {
    A[i]=ECEPE_p.dihed[i].A;
    ns[i]=ECEPE_p.dihed[i].NS;
    nb[i]=ECEPE_p.dihed[i].NB;
    dih[i]=ECEPE_p.dihed[i].angle;
    iftors[i]=ECEPE_p.dihed[i].IFTOR;
  }

  for (i=0;i<numatom;++i) {
    charge[i]=ECEPE_p.atom[i].charge;
    nbtype[i]=ECEPE_p.atom[i].nbtype;
    for (j=0;j<3;++j)
      co[i*3+j]=ECEPE_p.atom[i].refcoord[j];
  }

  //  set_pairs_ECEPE();
  calc_TORS_ECEPE(&t,n,A,ns,nb,dih,iftors,numdih);
  calc_NB_ECEPE_for_db(&p_nb,f_nb,&p_es,f_es,co,pairs_1_5,pairs_1_4,num1_5,num1_4,numatom,charge,nbtype,nb_p);

  (*p).p_t=t+p_nb+p_es;
  (*p).p_tors=t;
  (*p).p_nb=p_nb;
  (*p).p_es=p_es;


  free(co);        // 2015-02-16
  free(A);         // 2015-02-16
  free(ns);        // 2015-02-16
  free(nb);        // 2015-02-16
  free(dih);       // 2015-02-16
  free(iftors);    // 2015-02-16
  free(charge);    // 2015-02-16
  free(nbtype);    // 2015-02-16

  free(f_nb);      // 2015-02-16
  free(f_es);      // 2015-02-16
  free(n);         // 2015-02-16

}


double calc_TORS_ECEPE(double *t,double *n,double *A,int *ns,int *nb,double *dih, int *iftors,int numdih ){
  int i;
  
  *t=0.0;

  for (i=0;i<numdih;++i) {
    if (iftors[i]!=0) {
      *t+=A[i]*(1.0+(double)ns[i]*cos((double)nb[i]*dih[i]));
      n[i]=-A[i]*((double)nb[i]*(double)ns[i]*sin((double)nb[i]*dih[i]))*fconst;
    }
  }
}

double calc_TORS_ECEPE2(double *crd, struct ECEPE_parms p, double *delta_dihed, double *p_d, double *Q,int numclut,int *origin){ // 2015-02-16
  int i,j,flag;
  int numdih,nNumClut;
  double dihed;
  double atom_i[3],atom_j[3],atom_k[3],atom_l[3];
  double pi; //0811
  int ni,nj,nk,nl;

  double *p_d_each; // 2015-02-16

  p_d_each=(double *)calloc(numclut,sizeof(double)); // 2015-02-16

  pi=acos(-1.0);
  /*************************************************/
  /* double *A;					   */
  /* int *ns,*nb,*iftors;			   */
  /* 						   */
  /* numdih=p.NUMVAR;				   */
  /* A=(double *)gcemalloc(sizeof(double)*numdih); */
  /* ns=(int *)gcemalloc(sizeof(int)*numdih);	   */
  /* nb=(int *)gcemalloc(sizeof(int)*numdih);	   */
  /* iftors=(int *)gcemalloc(sizeof(int)*numdih);  */
  /* 						   */
  /* for (i=0;i<numdih;++i) {			   */
  /*   A[i]=p.dihed[i].A;			   */
  /*   ns[i]=p.dihed[i].NS;			   */
  /*   nb[i]=p.dihed[i].NB;			   */
  /*   iftors[i]=p.dihed[i].IFTOR;		   */
  /* }						   */
  /*************************************************/
 
  numdih=p.NUMVAR;
  for (i=0;i<numdih;++i) {
    //    ni=p.atom[p.dihed[i].dpairs[0]].katom-1;
    //    nj=p.atom[p.dihed[i].dpairs[1]].katom-1;
    //    nk=p.atom[p.dihed[i].dpairs[2]].katom-1;
    //    nl=p.atom[p.dihed[i].dpairs[3]].katom-1;
    for (j=0;j<3;++j) {
      atom_i[j]=crd[(p.dihed[i].dpairs[0])*3+j]; // 2015-02-16
      atom_j[j]=crd[(p.dihed[i].dpairs[1])*3+j]; // 2015-02-16
      atom_k[j]=crd[(p.dihed[i].dpairs[2])*3+j]; // 2015-02-16
      atom_l[j]=crd[(p.dihed[i].dpairs[3])*3+j]; // 2015-02-16
      //      atom_i[j]=prot.coord[p.dihed[i].dpairs[0]][j];  // 2015-02-16
      //      atom_j[j]=prot.coord[p.dihed[i].dpairs[1]][j];  // 2015-02-16
      //      atom_k[j]=prot.coord[p.dihed[i].dpairs[2]][j];  // 2015-02-16
      //      atom_l[j]=prot.coord[p.dihed[i].dpairs[3]][j];  // 2015-02-16
      //      atom_i[j]=prot.coord[ni][j];
      //      atom_j[j]=prot.coord[nj][j];
      //      atom_k[j]=prot.coord[nk][j];
      //      atom_l[j]=prot.coord[nl][j];
    }

    dihed=dih(atom_i,atom_j,atom_k,atom_l);
    if (dihed>pi) dihed-=2.0*pi;       // 0811
    else if (dihed<-1.0*pi) dihed+=2.0*pi;  // 0811

    dihed+=delta_dihed[i];
    flag=OFF;
    for (j=0;j</*prot.DOF*/numclut;++j) {
      //      if (clust[j].origin_atom_a-1==p.dihed[i].dpairs[1]) {
      //      if (clust[j].origin_atom_a-1==p.dihed[i].dpairs[2]) { // 2015-02-16
      if (origin[j]-1==p.dihed[i].dpairs[2]) {         // 2015-02-16
	nNumClut=j;
	flag=ON;
	break;
      }
    }
    if (p.dihed[i].IFTOR!=0 && flag == ON) {
      //      potential_pro.p_dihedc[nNumClut]=p.dihed[i].A*(1.0+(double)p.dihed[i].NS*cos((double)p.dihed[i].NB*dihed));
      *p_d+=p.dihed[i].A*(1.0+(double)p.dihed[i].NS*cos((double)p.dihed[i].NB*dihed)); // 2015-02-16

      p_d_each[nNumClut]=p.dihed[i].A*(1.0+(double)p.dihed[i].NS*cos((double)p.dihed[i].NB*dihed)); // 2015-02-16

      //      clust[nNumClut].f_c.f_dihed=-p.dihed[i].A*((double)p.dihed[i].NB*(double)p.dihed[i].NS*sin((double)p.dihed[i].NB*dihed))*fconst; // 2015-02-16
      Q[nNumClut]=-p.dihed[i].A*((double)p.dihed[i].NB*(double)p.dihed[i].NS*sin((double)p.dihed[i].NB*dihed))*fconst; // 2015-02-16
    }
  }

  free(p_d_each); // 2015-02-16

}

double calc_NB_ECEPE(double *p_nb, double *f_nb,double *p_es,double *f_es,double *co,int *pairs,int numint,double *charge, int *nbtype,struct pnb nb_p) {
  int i,j,k;
  int indxa,indxb,indxpb;
  double f;
  double vec[3];
  double len2,ro2,ro4,ro6=1.0,ro12;

  *p_es=0.0;
  *p_nb=0.0;
  for (i=0;i<numint;++i) {
    indxa=abs(pairs[i*2])-1;
    indxb=abs(pairs[i*2+1])-1;
    len2 = 0.0;
    for(j=0;j<3;++j){
      vec[j] = co[indxb*3+j]-co[indxa*3+j];
      len2 += vec[j]*vec[j];
    }

    f=charge[indxa]*charge[indxb]/sqrt(len2);
    *p_es+=charge[indxa]*charge[indxb]/sqrt(len2);
    indxpb=(nbtype[indxa]-1)*numatomtype+(nbtype[indxb]-1);
    ro2=nb_p.Bcff[indxpb]/len2;

    //    printf("%d\n",indxpb);
    // printf("ro2=%lf\n",ro2);
    // printf("len2=%lf\n",len2);
    //    printf("i=%d indxa=%d indxb=%d\n",i,indxa,indxb);
    ro6=1.0;
    for (j=0;j<3;++j)  ro6 = ro6*ro2;
    /**********************************************************/

    /*   if (nb_p[indxpb] >0)				      */
    /* 	if (pairs[i*2]<0) /\*for 14 int*\/		      */
    /* 	  *p_nb += nb_p[indxpb]*ro6*(0.5*ro6 - 2.0);	      */
    /* 	else						      */
    /* 	  *p_nb += nb_p[indxpb]*ro6*(ro6 - 2.0);		      */
    /*   else { 					      */
    /* 	ro4=ro2*ro2;					      */
    /* 	*p_nb -= nb_p[indxpb]*ro6*(ro6 - 2.0*ro4);	      */
    /*   }						      */
    /* else						      */
    /**********************************************************/
    if (pairs[i*2] <0) { /*for 14 int*/
      //      printf("%lf\n",nb_p.Bcff[indxpb]*ro6*(ro6 - 2.0));
      if (nb_p.Acff[indxpb] > 0.0) {
	f=nb_p.Acff[indxpb]*ro6*(0.5*ro6 - 2.0);
	//	printf("%10.4e\n",f);
	*p_nb += nb_p.Acff[indxpb]*ro6*(0.5*ro6 - 2.0);
      }
      else {
	ro4=ro2*ro2;
	f=-nb_p.Acff[indxpb]*ro6*(ro6 - 2.0*ro4);
	//	printf("%10.4e\n",f);
	*p_nb -= nb_p.Acff[indxpb]*ro6*(ro6 - 2.0*ro4);
      }
    }
    else { 
      if (nb_p.Acff[indxpb] > 0.0) {
	f=nb_p.Acff[indxpb]*ro6*(ro6 - 2.0);
	//	printf("%10.4e\n",f);
	*p_nb +=  nb_p.Acff[indxpb]*ro6*(ro6 - 2.0);
      }
      else {
	ro4=ro2*ro2;
	f=-nb_p.Acff[indxpb]*ro6*(ro6 - 2.0*ro4);
	//	printf("%10.4e\n",f);
	*p_nb -= nb_p.Acff[indxpb]*ro6*(ro6 - 2.0*ro4);
      }
      //      printf("%lf\n",nb_p.Bcff[indxpb]*ro6*(ro6 - 2.0*ro4));
    }
    f=charge[indxa]*charge[indxb]/sqrt(len2);
    /***************************/
    /* if (pairs[i*2] >0)      */
    /*   printf("%10.4e\n",f); */
    /***************************/
  /*****/
  /* } */
  /*****/
  }
}

double calc_NB_ECEPE_for_db(double *p_nb, double *f_nb,double *p_es,double *f_es,double *co,int **pairs_1_5,int **pairs_1_4,int *numnb,int *num14,int numatom,double *charge, int *nbtype,struct pnb nb_p) {
  int i,j,k,l;
  int indxa,indxb,indxpb;
  double f;
  double vec[3];
  double len2,ro2,ro4,ro6=1.0,ro12;

  double *p_nb_each;
  double *p_es_each;

  int nNumClut,num_a_prot,alpha,tnum_atom_clust;

  FILE *debug;

  debug=efopen("ene1_5_2_tamd.txt","w");

  //  p_nb_each=(double *)gcemalloc(sizeof(double)*numatom); // 2015-02-16
  //  p_es_each=(double *)gcemalloc(sizeof(double)*numatom); // 2015-02-16
  p_nb_each=(double *)calloc(numatom,sizeof(double));       // 2015-02-16
  p_es_each=(double *)calloc(numatom,sizeof(double));	    // 2015-02-16
  for (i=0;i<numatom;++i) {
    p_nb_each[i]=0.0;
    p_es_each[i]=0.0;
    for (j=0;j<3;++j) {
      f_nb[i*3+j]=0.0;
      f_es[i*3+j]=0.0;
    }
  }

  *p_es=0.0;
  *p_nb=0.0;
  for (i=0;i<numatom;++i) {
    if (flagnb14==ON) {
      for (j=0;j<num14[i];++j) {
	indxa=i;
	indxb=pairs_1_4[i][j];
	len2 = 0.0;
	for(k=0;k<3;++k){
	  vec[k] = co[indxb*3+k]-co[indxa*3+k];
	  len2 += vec[k]*vec[k];
	}
	indxpb=(nbtype[indxa]-1)*numatomtype+(nbtype[indxb]-1);
	ro2=nb_p.Bcff[indxpb]/len2;
      
	ro6=1.0;
	for (k=0;k<3;++k)  ro6 = ro6*ro2;
	if (nb_p.Acff[indxpb] > 0.0) {
	  p_nb_each[i]+=nb_p.Acff[indxpb]*ro6*(0.5*ro6 - 2.0);
	  //	  printf ("1-4NB:%d %d %8.6lf \n",indxa,indxb,nb_p.Acff[indxpb]*ro6*(0.5*ro6 - 2.0)/*f*/);
	  *p_nb += nb_p.Acff[indxpb]*ro6*(0.5*ro6 - 2.0);
	  for (l=0;l<3;++l) {
	    f_nb[indxa*3+l] += -nb_p.Acff[indxpb]*ro6*(0.5*ro6*12.0 - 2.0*6.0)/len2*vec[l]*fconst;
	    f_nb[indxb*3+l] -= -nb_p.Acff[indxpb]*ro6*(0.5*ro6*12.0 - 2.0*6.0)/len2*vec[l]*fconst;
	  }
	}
	else {
	  ro4=ro2*ro2;
	  p_nb_each[i]=-nb_p.Acff[indxpb]*ro6*(ro6 - 2.0*ro4);
	  //	  printf ("1-4NB:%d %d %8.6lf \n",indxa,indxb,/*f*/-nb_p.Acff[indxpb]*ro6*(ro6 - 2.0*ro4));
	  *p_nb -= nb_p.Acff[indxpb]*ro6*(ro6 - 2.0*ro4);
	  for (l=0;l<3;++l) {
	    f_nb[indxa*3+l] += nb_p.Acff[indxpb]*ro6*(ro6*12.0 - 2.0*ro4*10.0)/len2*vec[l]*fconst;
	    f_nb[indxb*3+l] -= nb_p.Acff[indxpb]*ro6*(ro6*12.0 - 2.0*ro4*10.0)/len2*vec[l]*fconst;
	  }
	}
      }
    }

    if (flages14==ON) {
      for (j=0;j<num14[i];++j) {
	indxa=i;
	indxb=pairs_1_4[i][j];
	len2 = 0.0;
	for(k=0;k<3;++k){
	  vec[k] = co[indxb*3+k]-co[indxa*3+k];
	  len2 += vec[k]*vec[k];
	}
      
	f=charge[indxa]*charge[indxb]/sqrt(len2);
	//      printf ("%d %d %10.4e 1-4es\n",indxa,indxb,f);
	p_es_each[i]+=charge[indxa]*charge[indxb]/sqrt(len2);
	*p_es+=charge[indxa]*charge[indxb]/sqrt(len2);
	for (l=0;l<3;++l) {
	  f_es[indxa*3+l] += -charge[indxa]*charge[indxb]/sqrt(len2)/len2*vec[l]*fconst;
	  f_es[indxb*3+l] -= -charge[indxa]*charge[indxb]/sqrt(len2)/len2*vec[l]*fconst;
	}
      }
    }

    if (flagnb15==ON) {
      for (j=0;j<numnb[i];++j) {
	indxa=i;
	indxb=pairs_1_5[i][j];
	len2 = 0.0;
	for(k=0;k<3;++k){
	  vec[k] = co[indxb*3+k]-co[indxa*3+k];
	  len2 += vec[k]*vec[k];
	}
	indxpb=(nbtype[indxa]-1)*numatomtype+(nbtype[indxb]-1);
	ro2=nb_p.Bcff[indxpb]/len2;
	
	ro6=1.0;
	for (k=0;k<3;++k)  ro6 = ro6*ro2;
	if (nb_p.Acff[indxpb] > 0.0) {
	  f=nb_p.Acff[indxpb]*ro6*(ro6 - 2.0)*100;
	  //	printf ("%d %d %10.4e 1-5nb\n",indxa,indxb,f);
	  fprintf (debug,"1-5NB:%d %d %8.6lf \n",indxa,indxb,nb_p.Acff[indxpb]*ro6*(ro6 - 2.0));
	  p_nb_each[i] += nb_p.Acff[indxpb]*ro6*(ro6 - 2.0);
	  *p_nb +=  nb_p.Acff[indxpb]*ro6*(ro6 - 2.0);
	  for (l=0;l<3;++l) {
	    f_nb[indxa*3+l] += -nb_p.Acff[indxpb]*ro6*(ro6*12.0 - 2.0*6.0)/len2*vec[l]*fconst;
	    f_nb[indxb*3+l] -= -nb_p.Acff[indxpb]*ro6*(ro6*12.0 - 2.0*6.0)/len2*vec[l]*fconst;
	  }
	}
	else {
	  ro4=ro2*ro2;
	  f=-nb_p.Acff[indxpb]*ro6*(ro6 - 2.0*ro4);
	  //	printf ("%d %d %10.4e 1-5nb\n",indxa,indxb,f);
	  fprintf(debug,"1-5NB:%d %d %8.6lf \n",indxa,indxb,nb_p.Acff[indxpb]*ro6*(ro6 - 2.0*ro4));
	  p_nb_each[i] += -nb_p.Acff[indxpb]*ro6*(ro6 - 2.0*ro4);
	  *p_nb -= nb_p.Acff[indxpb]*ro6*(ro6 - 2.0*ro4);
	  for (l=0;l<3;++l) {
	    f_nb[indxa*3+l] += nb_p.Acff[indxpb]*ro6*(ro6*12.0 - 2.0*ro4*10.0)/len2*vec[l]*fconst;
	    f_nb[indxb*3+l] -= nb_p.Acff[indxpb]*ro6*(ro6*12.0 - 2.0*ro4*10.0)/len2*vec[l]*fconst;
	  }
	}
      }
    }

    if (flages15==ON) {
      for (j=0;j<numnb[i];++j) {
	indxa=i;
	indxb=pairs_1_5[i][j];
	len2 = 0.0;
	for(k=0;k<3;++k){
	  vec[k] = co[indxb*3+k]-co[indxa*3+k];
	  len2 += vec[k]*vec[k];
	}
	indxa=i;
	indxb=pairs_1_5[i][j];
	f=charge[indxa]*charge[indxb]/sqrt(len2);
	p_es_each[i] += charge[indxa]*charge[indxb]/sqrt(len2);
	*p_es+=charge[indxa]*charge[indxb]/sqrt(len2);
	for (l=0;l<3;++l) {
	  f_es[indxa*3+l] += -charge[indxa]*charge[indxb]/sqrt(len2)/len2*vec[l]*fconst;
	  f_es[indxb*3+l] -= -charge[indxa]*charge[indxb]/sqrt(len2)/len2*vec[l]*fconst;
	  //      printf ("%d %d %10.4e 1-5es\n",indxa,indxb,f);
	}
      }
    }
  }

  // 2015-02-16
  /***********************************************************************************/
  /* num_a_prot=0;								     */
  /* for(nNumClut = 0;nNumClut < prot.DOF; ++nNumClut) {			     */
  /*   tnum_atom_clust = clust[nNumClut].num_atom_clust;			     */
  /*   for(i=0;i < tnum_atom_clust;++i) {					     */
  /*     potential_pro.p_L_J[num_a_prot] = p_nb_each[num_a_prot]*2.0;		     */
  /*     potential_pro.p_elesta[num_a_prot] = p_es_each[num_a_prot]*2.0;	     */
  /*     for(alpha=0; alpha<3; ++alpha) {					     */
  /* 	clust[nNumClut].f_c.f_L_J[i][alpha] = f_nb[num_a_prot*3+alpha];		     */
  /* 	clust[nNumClut].f_c.f_elesta[i][alpha] = f_es[num_a_prot*3+alpha];	     */
  /*     }									     */
  /*     ++num_a_prot;								     */
  /*   }									     */
  /* }										     */
  /* 										     */
  /* f=0;									     */
  /* for (i=0;i<numatom;++i)							     */
  /*   f+=p_nb_each[i];								     */
  /* f=0;									     */
  /* for (i=0;i<numatom;++i)							     */
  /*   f+=p_es_each[i];								     */
  /* f=0.0;									     */
  /* for (i=0;i<numatom;++i)							     */
  /*   f+=potential_pro.p_L_J[i];						     */
  /* f=0.0;									     */
  /* for (i=0;i<numatom;++i)							     */
  /*   f+=potential_pro.p_elesta[i];						     */
  /***********************************************************************************/
  // 2015-02-16

  fclose(debug);
}

int set_pairs_ECEPE(int *pairs,int numint,int *nbtype) {

}

void read_ECEPE_coo (FILE *file, double *co, double *dihed, int numatom) {
  int i,j,k,nd,na;
  int flag;
  char *line,x[6],y[6],z[6],d[3],khi[4];
  size_t len=0;

  nd=0;
  na=0;
  i=0;
  for (;;) {
    getline(&line,&len,file);
    if (line[0]=='4') {
      flag=ON;
      for (;nd<3;) {
	for (j=0;j<8;++j) {
	  d[j]=line[48+nd*8+j];
	}
	dihed[nd]=atof(d);
	++nd;
      }
    }
    if (flag==ON) {
      for (i=0;i<4;++i)
	khi[i]=line[8+i];
      if (strncmp(khi,"EKHI",4)==0) {
	for (i=0;i<7;++i)
	  d[i]=line[17+i];
	dihed[nd]=atof(d);
	++nd;
      }
      for (i=0;i<4;++i)
	khi[i]=line[40+i];
      if (strncmp(khi,"EKHI",4)==0) {
	for (i=0;i<7;++i)
	  d[i]=line[49+i];
	dihed[nd]=atof(d);
	++nd;
      }
      if (line[0]=='5') {
	//      printf("%s\n",line);
	for (j=0;j<6;++j) {
	  x[j]=line[50+j];
	  y[j]=line[58+j];
	  z[j]=line[66+j];
	  co[na*3]=atof(x);
	  co[na*3+1]=atof(y);
	  co[na*3+2]=atof(z);
	}
	++na;
      }
    }
    if (strncmp(line,"END",3)==0) {
      break;
    }
  }

  /************************************/
  /* for (j=0;j<i;++j) {	      */
  /*   for (k=0;k<3;++k)	      */
  /*     printf("%6.4lf ",co[j*3+k]); */
  /*   printf("\n");		      */
  /* }				      */
  /************************************/

}

//int make_int_pair_list(int **bp,int *numb,int numatom,int numbond, int **pair1_5, int **pair1_4, int *num14,int *num1_5) {
int make_int_pair_list(int **bp,int *numb,int numatom, int **pair1_5, int **pair1_4, int *num14,int *num1_5) {
  int i,j,k,l,n,m;
  int temp;
  int numatom_1,numatom_2,numatom_3,numatom_4;
  int *numex15;
  int **pair_ex_1_5;
  
  pair_ex_1_5=(int **)calloc(numatom,sizeof(int *)); // 2015-02-16
  //  pair_ex_1_5=(int **)gcemalloc(sizeof(int *)*numatom); // 2015-02-16
  //  num14=(int *)gcemalloc(sizeof(int)*numatom);
  //  num1_5=(int *)gcemalloc(sizeof(int)*numatom);
  numex15=(int *)calloc(numatom,sizeof(int));        // 2015-02-16
  //  numex15=(int *)gcemalloc(sizeof(int)*numatom); // 2015-02-16
  
  for (i=0;i<numatom;++i) {
    num14[i]=0;
    numex15[i]=0;
    for (j=0;j<numb[i];++j) {
      numatom_1=bp[i][j];
      ++numex15[i];
      //      pair_ex_1_5[i]=(int *)gcerealloc(pair_ex_1_5[i],sizeof(int)*numex15[i]); // 2015-02-16
      pair_ex_1_5[i]=(int *)realloc(pair_ex_1_5[i],sizeof(int)*numex15[i]); // 2015-02-16
      pair_ex_1_5[i][numex15[i]-1]=numatom_1;
      for (k=0;k<numb[numatom_1];++k) {
	numatom_2=bp[numatom_1][k];
	if (numatom_2!=i) {
	  ++numex15[i];
	  //	  pair_ex_1_5[i]=(int *)gcerealloc(pair_ex_1_5[i],sizeof(int)*numex15[i]); // 2015-02-16
	  pair_ex_1_5[i]=(int *)realloc(pair_ex_1_5[i],sizeof(int)*numex15[i]); // 2015-02-16
	  pair_ex_1_5[i][numex15[i]-1]=numatom_2;
	  for (l=0;l<numb[numatom_2];++l) {
	    numatom_3=bp[numatom_2][l];
	    if (numatom_3!=numatom_1) {
	      ++numex15[i];
	      //	      pair_ex_1_5[i]=(int *)gcerealloc(pair_ex_1_5[i],sizeof(int)*numex15[i]); // 2015-02-16
	      pair_ex_1_5[i]=(int *)realloc(pair_ex_1_5[i],sizeof(int)*numex15[i]);  // 2015-02-16
	      pair_ex_1_5[i][numex15[i]-1]=numatom_3;
	      if (numatom_3 > i) {
		++num14[i];
		//		pair1_4[i]=(int *)gcerealloc(pair1_4[i],sizeof(int)*num14[i]); // 2015-02-16
		pair1_4[i]=(int *)realloc(pair1_4[i],sizeof(int)*num14[i]); // 2015-02-16
		pair1_4[i][num14[i]-1]=numatom_3;
	      }
	    }
	  }
	}
      }
    }
  }

  for (i=0;i<numatom;++i) {
    num1_5[i]=1;
    for (j=i+1;j<numatom;++j) {
      if (lookup_ex_list(pair_ex_1_5,i,j,numex15)==1) {
	//  	pair1_5[i]=(int *)gcerealloc(pair1_5[i],sizeof(int)*num1_5[i]); // 2015-02-16
  	pair1_5[i]=(int *)realloc(pair1_5[i],sizeof(int)*num1_5[i]); // 2015-02-16
  	pair1_5[i][num1_5[i]-1]=j;
  	++num1_5[i];
      }
    }
  }

  for (i=0;i<numatom;++i) {
    for (j=1;j<num14[i];++j) {
      for (k=j;k>0;--k) {
	if (pair1_4[i][k]<pair1_4[i][k-1]) {
	  temp=pair1_4[i][k];
	  pair1_4[i][k]=pair1_4[i][k-1];
	  pair1_4[i][k-1]=temp;
	}
      }
    }
  }

  for (i=0;i<numatom;++i) {
    for (j=1;j<num1_5[i]-1;++j) {
      for (k=j;k>0;--k) {
	if (pair1_5[i][k]<pair1_5[i][k-1]) {
	  temp=pair1_5[i][k];
	  pair1_5[i][k]=pair1_5[i][k-1];
	  pair1_5[i][k-1]=temp;
	}
      }     
    }
  }

  /***********************************/
  /* for (i=0;i<numatom;++i) {	     */
  /*   pair1_5[i][num1_5[i]-1]='\0'; */
  /* }				     */
  /* for (i=0;i<numatom;++i) {	     */
  /*   pair1_4[i][num14[i]-1]='\0';  */
  /* }				     */
  /***********************************/

  /*****************************************************/
  /* for (i=0;i<numatom;++i) {			       */
  /*   for (j=0;j<num14[i];++j) {		       */
  /*     printf("%d to %d\n",i+1,pair1_4[i][j]+1);     */
  /*   }					       */
  /* }						       */
  /* printf("\n");				       */
  /* 						       */
  /* for (i=0;i<numatom;++i) {			       */
  /*   for (j=0;j<numex15[i];++j) {		       */
  /*     printf("%d to %d\n",i+1,pair_ex_1_5[i][j]+1); */
  /*   }					       */
  /* }						       */
  /* printf("\n");				       */
  /*****************************************************/

  free(pair_ex_1_5); // 2015-02-16
  free(numex15);     // 2015-02-16

}

int lookup_ex_list(int **pair_ex_1_5,int numi, int numj, int *numex15) {
  int i,j;

  for (i=0;i<numex15[numi];++i) {
    if (pair_ex_1_5[numi][i]==numj)
      return 0;
  }

  return 1;
}



int make_bd_pair_list(struct ECEPE_parms pa,int *bp,int *bp_f , int *numb) {
  int i,j,k;
  int numbd=0;

  gchar *atomnameaskey;
  gchar *atomname;
  char *name;
  gint *atonum;
  GHashTable *bplist,*bplist_c2,*bplist_c3,*bplist_c4,*bplist_c5;

  int ap[2],apflag;

  bplist=g_hash_table_new_full(g_str_hash,g_str_equal,g_free,g_free);
  bplist_c2=g_hash_table_new_full(g_str_hash,g_str_equal,g_free,g_free);
  bplist_c3=g_hash_table_new_full(g_str_hash,g_str_equal,g_free,g_free);
  bplist_c4=g_hash_table_new_full(g_str_hash,g_str_equal,g_free,g_free);
  bplist_c5=g_hash_table_new_full(g_str_hash,g_str_equal,g_free,g_free);

  for (i=0;i<numatomtype2;++i) {
    atomnameaskey=g_strdup(atomtypelist[i*2]);
    atomname=g_new(gchar,1);atomname=atomtypelist[i*2+1];
    g_hash_table_insert(bplist,atomnameaskey,atomname);
    if ((name=g_hash_table_lookup(bplist,atomtypelist[i*2]))!=NULL) {
      //      printf("%s - %s\n",atomtypelist[i*2],name);
      ;
    }
  }
  for (i=0;i<numatomtype3;++i) {
    atomnameaskey=g_strdup(atomtypelist_c2[i*2]);
    atomname=g_new(gchar,1);atomname=atomtypelist_c2[i*2+1];
    g_hash_table_insert(bplist_c2,atomnameaskey,atomname);
    if ((name=g_hash_table_lookup(bplist_c2,atomtypelist_c2[i*2]))!=NULL) {
      //      printf("%s - %s\n",atomtypelist_c2[i*2],name);
      ;
    }
  }
  for (i=0;i<numatomtype4;++i) {
    atomnameaskey=g_strdup(atomtypelist_c3[i*2]);
    atomname=g_new(gchar,1);atomname=atomtypelist_c3[i*2+1];
    g_hash_table_insert(bplist_c3,atomnameaskey,atomname);
    if ((name=g_hash_table_lookup(bplist_c3,atomtypelist_c3[i*2]))!=NULL) {
      //      printf("%s - %s\n",atomtypelist_c3[i*2],name);
      ;
    }
  }
  atomnameaskey=g_strdup(atomtypelist_c4[0]);
  atomname=g_new(gchar,1);atomname=atomtypelist_c4[1];
  g_hash_table_insert(bplist_c4,atomnameaskey,atomname);

  atomnameaskey=g_strdup(atomtypelist_c5[0]);
  atomname=g_new(gchar,1);atomname=atomtypelist_c5[1];
  g_hash_table_insert(bplist_c5,atomnameaskey,atomname);

  for (i=1;i<pa.NUMATM;++i) {
    sea_pair_by_name(pa.atom[i].name_atom,i,bp,ap,apflag,bplist,bplist_c2,bplist_c3,bplist_c4,bplist_c5,pa);
  }

  /*****************************************************/
  /* if (apflag==ON) {				       */
  /*   bp=(int *)gcerealloc(bp,sizeof(int)*numatom*2); */
  /*   bp[numatom*2]=ap[0];			       */
  /*   bp[numatom*2+1]=ap[1];			       */
  /* }						       */
  /*****************************************************/

  order_bd_pair_list(bp, bp_f, pa.NUMATM,numb);

  //  g_hash_table_destroy(bplist);
}

void make_dihed_pairs_list( struct ECEPE_parms p, int **bpl, int *numb ) {
  int i,j,k,l;
  int a1,a2,a3,a4;

  for (i=0;i<p.NUMVAR;++i) {
    a2=p.atom[p.dihed[i].ibnd1-1].katom-1;
    a3=p.atom[p.dihed[i].ibnd2-1].katom-1;
    for (j=0;j<numb[a2];++j) {
      if ((k=bpl[a2][j])!=a3) {
	a1=k;
	break;
      }
    }
    for (j=0;j<numb[a3];++j) {
      if ((k=bpl[a3][j])!=a2) {
	a4=k;
	break;
      }
    }
    p.dihed[i].dpairs[0]=a1;
    p.dihed[i].dpairs[1]=a2;
    p.dihed[i].dpairs[2]=a3;
    p.dihed[i].dpairs[3]=a4;
  }
}

void sea_pair_by_name(char *name, int nb, int *bp,int ap[2], int apflag,GHashTable *bplist,GHashTable *bplist_c2,GHashTable *bplist_c3,GHashTable *bplist_c4,GHashTable *bplist_c5,struct ECEPE_parms pa) {
  int i,j,k;
  int *num_p;
  char *name_p;
  char *name_d,*name_d2,*name_d3;

  int flag;
  int flag2;  

  gint *atomnum;
  gchar *atomnameaskey;

  //  name_d=(char *)gcemalloc(sizeof(char)*4);  // 2015-02-16
  //  name_d2=(char *)gcemalloc(sizeof(char)*2); // 2015-02-16
  //  name_d3=(char *)gcemalloc(sizeof(char)*3); // 2015-02-16
  name_d=(char *)calloc(4,sizeof(char));   // 2015-02-16
  name_d2=(char *)calloc(2,sizeof(char));  // 2015-02-16
  name_d3=(char *)calloc(3,sizeof(char));  // 2015-02-16

  for (i=0;i<4;++i) name_d[i]=name[i];
  for (i=0;i<2;++i) name_d2[i]=name[i];
  for (i=0;i<3;++i) name_d3[i]=name[i];

  flag=OFF;
  flag2=OFF;
  apflag=OFF;
  for (i=nb-1;i>=0;--i) {
    if ((name_p=g_hash_table_lookup(bplist,name_d3))==NULL) {
      printf("There is no key %s in hashtable",name_d3);
      exit(1);
    }
    if (strncmp(pa.atom[i].name_atom,name_p,3)==0) {
      bp[(nb-1)*2]=i;
      bp[(nb-1)*2+1]=nb;
      flag=ON;
      break;
    }
    if (flag==OFF) {
      if ((name_p=g_hash_table_lookup(bplist_c2,name_d3))!=NULL) {
	if (strncmp(pa.atom[i].name_atom,name_p,3)==0) {
	  bp[(nb-1)*2]=i;
	  bp[(nb-1)*2+1]=nb;
	  flag=ON;
	  break;
	}
      }
    }
    if (flag==OFF) {
      if ((name_p=g_hash_table_lookup(bplist_c3,name_d))!=NULL) {
	if (strncmp(pa.atom[i].name_atom,name_p,4)==0) {
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
	  if (strncmp(pa.atom[i].name_atom,name_p,4)==0) {
	    bp[(nb-1)*2]=i;
	    bp[(nb-1)*2+1]=nb;
	    flag2=ON;
	  }
	}
      }
      if (flag2==ON) {
	for (j=nb-1;j>=0;--j) {
	  if ((name_p=g_hash_table_lookup(bplist_c5,name_d))!=NULL) {
	    if (strncmp(pa.atom[j].name_atom,name_p,4)==0) {
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
	  for (j=nb+1;j<pa.NUMATM;++j) {
	    if ((name_p=g_hash_table_lookup(bplist_c5,name_d))!=NULL) {
	      if (strncmp(pa.atom[j].name_atom,name_p,4)==0) {
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
    for (i=nb+1;i<pa.NUMATM;++i) {
      if ((name_p=g_hash_table_lookup(bplist,name_d3))==NULL) {
	printf("There is no key %s in hashtable",name_d3);
	exit(1);
      }
      if (strncmp(pa.atom[i].name_atom,name_p,3/*2*//*3*/)==0) {
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
	  if (strncmp(pa.atom[i].name_atom,name_p,3/*2*//*4*/)==0) {
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
	  if (strncmp(pa.atom[i].name_atom,name_p,/*2*/4)==0) {
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
	  if (strncmp(pa.atom[i].name_atom,name_p,4)==0) {
	    bp[(nb-1)*2]=i;
	    bp[(nb-1)*2+1]=nb;
	    flag2=ON;
	  }
	}
      }
      if (flag2==ON) {
	for (j=nb-1;j>=0;--j) {
	  if ((name_p=g_hash_table_lookup(bplist_c5,name_d))!=NULL) {
	    if (strncmp(pa.atom[j].name_atom,name_p,4)==0) {
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
	  for (i=nb+1;i<pa.NUMATM;++i) {
	    if ((name_p=g_hash_table_lookup(bplist_c5,name_d))!=NULL) {
	      if (strncmp(pa.atom[j].name_atom,name_p,4)==0) {
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

  if (flag==OFF)
    printf("error: there is no candidate to bond %s\n",name);
  else
    printf("%s %d to %s %d\n",name,nb,name_p,i);

  free(name_d);   // 2015-02-16
  free(name_d2);  // 2015-02-16
  free(name_d3);  // 2015-02-16

}

int order_bd_pair_list(int *bp, int **bp_f, int numatom, int *numb) {
  int i,j,k;

  for (i=0;i<numatom;++i) numb[i]=0;

  for (i=0;i<numatom;++i) {
    for (j=0;j<numatom-1;++j) {
      if (bp[j*2]==i /*&& bp[j*2] < bp[j*2+1]*/ ) {
	numb[i]+=1;
	//	bp_f[i]=(int *)gcerealloc(bp_f[i],sizeof(int)*numb[i]); // 2015-02-16
	bp_f[i]=(int *)realloc(bp_f[i],sizeof(int)*numb[i]);            // 2015-02-16
	bp_f[i][numb[i]-1]=bp[j*2+1];
	/***********************************************************/
        /* if (i<bp_f[i][numb[i]-1])				   */
	/*   printf("%d - %d \n",i+1,bp_f[i][numb[i]-1]+1);	   */
        /***********************************************************/
      }
      else if (bp[j*2+1]==i /*&& bp[j*2+1] < bp[j*2]*/ ) {
	numb[i]+=1;
	//	bp_f[i]=(int *)gcerealloc(bp_f[i],sizeof(int)*numb[i]);  // 2015-02-16
	bp_f[i]=(int *)realloc(bp_f[i],sizeof(int)*numb[i]);             // 2015-02-16
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

/********************************************************************/
/* void print_iter(gpointer key, gpointer val, gpointer usr_data) { */
/*   g_debug("%s:%d",(gchar *)key,*(gint *)val);		    */
/* }								    */
/********************************************************************/

//double *pick_ECEPP_data(char *preofilename,char *bd8filename,char *coofilename,char *coofilename_for_sflag,int **pair1_5, int **pair1_4, int *num1_4,int *num1_5,struct ECEPE_parms *ECEPE_p, struct pnb *nb_p/*, double *delta_dihed*/,char *angfilename, int flagang, int flagpairs, char *pairsfilename,double *delta_dihed) { // 2015-02-16
double *pick_ECEPP_data(char *preofilename,char *bd8filename,char *coofilename,char *coofilename_for_sflag,int **pair1_5, int **pair1_4, int *num1_4,int *num1_5,struct ECEPE_parms *ECEPE_p, struct pnb *nb_p/*, double *delta_dihed*/,char *angfilename, int flagang, int flagpairs, char *pairsfilename,double *delta_dihed,double *dihed,double *co,char *InpfilCOORD, char *InpfilTOP/*, int oliflag*/) { // 2015-02-16
  int i,j,k,l,n1_4,n1_5;
  int *bp,**bp_f;
  int *numb;
  int numnb,num14;

  //  char *line="line"; // 2015-02-16
  //  size_t len=0;      // 2015-02-16

  //  double *co; // 2015-02-16

  double *co_dummy,*dihed_dummy,*ene,*co2;
  double pi;

  FILE *coofile,*coofile_for_sflag;

  FILE *angfile; // 0811
  int *numdiofres; // 0811
  int num; // 0811
  double fd; // 0811

  FILE *input;

  // double *delta_dihed; 0811
  double x;

  int **matrix;

  char *name_atom_list;

  double *dihed_dummy_dummy;
  int ntotaldih,ndihinres,numres;

  read_ECEPE_parm(preofilename,bd8filename,ECEPE_p,nb_p);
  //  bp=(int *)gcemalloc(sizeof(int)*((*ECEPE_p).NUMATM-1)*2); // 2015-02-16
  //  bp_f=(int **)gcemalloc(sizeof(int *)*(*ECEPE_p).NUMATM);  // 2015-02-16
  bp=(int *)calloc(((*ECEPE_p).NUMATM-1)*2,sizeof(int));     // 2015-02-16
  bp_f=(int **)calloc((*ECEPE_p).NUMATM,sizeof(int *));      // 2015-02-16

  //  numb=(int *)gcemalloc(sizeof(int)*(*ECEPE_p).NUMATM);  // 2015-02-16
  numb=(int *)calloc((*ECEPE_p).NUMATM,sizeof(int));         // 2015-02-16
  /*****************************************************************************************/
  /* make_bd_pair_list((*ECEPE_p),bp,bp_f,numb);					   */
  /* 											   */
  //  delta_dihed=(double *)gcemalloc(sizeof(double)*((*ECEPE_p).NUMVAR)); //0811
  /* 											   */
  pi=acos(-1.0);
  /* 											   */
  /* /\******************************************************\/				   */
  /* /\* printf("bond list\n");			        *\/				   */
  /* /\* for (i=0;i<(*ECEPE_p).NUMATM;++i) {	        *\/				   */
  /* /\*   for (j=0;j<numb[i];++j)			        *\/			   */
  /* /\*     if (i<bp_f[i][j])			        *\/				   */
  /* /\* 	printf("%d - %d\n",i+1,bp_f[i][j]+1);	        *\/			   */
  /* /\* }						        *\/			   */
  /* /\******************************************************\/				   */
  /* 											   */
  /* 											   */
  /* make_dihed_pairs_list((*ECEPE_p),bp_f,numb);					   */
  /* 											   */
  /* make_int_pair_list(bp_f,numb,(*ECEPE_p).NUMATM,pair1_5,pair1_4,num1_4,num1_5);	   */
  /* 											   */
  /* numnb=0;										   */
  /* num14=0;										   */
  /* for (i=0;i<(*ECEPE_p).NUMATM;++i)							   */
  /*   if ( num1_5[i]-1 > 0 )								   */
  /*     numnb+=num1_5[i]-1;								   */
  /* for (i=0;i<(*ECEPE_p).NUMATM;++i)							   */
  /*   if ( num1_4[i]-1 > 0)								   */
  /*     num14+=num1_4[i]-1;								   */
  /* 											   */
  /* for (i=0;i<(*ECEPE_p).NUMATM;++i) {						   */
  /*   num1_5[i]-=1;									   */
  /* }											   */
  /*****************************************************************************************/

  /////////////////////////////////////////
  //  name_atom_list=(char *)gcemalloc(sizeof(char)*(*ECEPE_p).NUMATM*4);  // 2015-02-16
  name_atom_list=(char *)calloc((*ECEPE_p).NUMATM*4,sizeof(char));         // 2015-02-16
  for (i=0;i<(*ECEPE_p).NUMATM;++i) for (j=0;j<4;++j)
	name_atom_list[i*4+j]=(*ECEPE_p).atom[i].name_atom[j];
  make_bp((*ECEPE_p).NUMATM,bp,bp_f,numb,name_atom_list,1,flagpairs);

  /*************************************************************/  // 2015-02-16
  /* for (i=0;i<(*ECEPE_p).NUMATM;++i) {		       */  // 2015-02-16
  /*   for (j=0;j<4;++j) printf("%c",name_atom_list[i*4+j]);   */  // 2015-02-16
  /*   printf("(%3d)-",i);				       */  // 2015-02-16
  /*   for (j=0;j<numb[i];++j) {			       */  // 2015-02-16
  /*     n=bp_f[i][j];					       */  // 2015-02-16
  /*     for (k=0;k<4;++k) printf("%c",name_atom_list[n*4+k]); */  // 2015-02-16
  /*     printf("(%3d) ",n);				       */  // 2015-02-16
  /*   }						       */  // 2015-02-16
  /*   printf("\n");					       */  // 2015-02-16
  /* }							       */  // 2015-02-16
  /*************************************************************/  // 2015-02-16


  make_dihed_pairs_list((*ECEPE_p),bp_f,numb);

  //  matrix=(int **)gcemalloc(sizeof(int *)*(*ECEPE_p).NUMATM); // 2015-02-16
  matrix=(int **)calloc((*ECEPE_p).NUMATM,sizeof(int *));  // 2015-02-16
  for (i=0;i<(*ECEPE_p).NUMATM;++i) {
    //    matrix[i]=(int *)gcemalloc(sizeof(int)*(*ECEPE_p).NUMATM); // 2015-02-16
    matrix[i]=(int *)calloc((*ECEPE_p).NUMATM,sizeof(int));       // 2015-02-16
  }

  make_nb_matrix(bp_f,numb,3,matrix,(*ECEPE_p).NUMATM);

  for (i=0;i<(*ECEPE_p).NUMATM;++i) {
    num1_5[i]=0;
    num1_4[i]=0;
  }

  set_nb_pairs(matrix,(*ECEPE_p).NUMATM,pair1_5,pair1_4,num1_5,num1_4);

  printf("nb list\n");
  for (i=0;i<(*ECEPE_p).NUMATM;++i) {
    printf("%d - ",i+1);
    for (j=0;j<num1_5[i];++j)
      printf("%d ",pair1_5[i][j]+1);
    printf("\n");
  }
  printf("14 list\n");
  for (i=0;i<(*ECEPE_p).NUMATM;++i) {
    printf("%d - ",i+1);
    for (j=0;j<num1_4[i];++j)
      printf("%d ",pair1_4[i][j]+1);
    printf("\n");
  }

  if ( flagpairs==ON ) {
    read_ECEPE_pairs(pairsfilename,pair1_5,pair1_4,num1_5,num1_4,(*ECEPE_p).NUMATM);
  }



  /**************************************************/
  /* for (i=0;i<(*ECEPE_p).NUMATM;++i) {	    */
  /*   n1_4=0;					    */
  /*   n1_5=0;					    */
  /*   for (j=i+1;j<(*ECEPE_p).NUMATM;++j) {	    */
  /*     if (matrix[i][j]==3) {			    */
  /* 	n1_4++;					    */
  /*     }					    */
  /*     if (matrix[i][j]==-1) {		    */
  /* 	n1_5++;					    */
  /*     }					    */
  /*   }					    */
  /*   num1_4[i]=n1_4;				    */
  /*   num1_5[i]=n1_5;				    */
  /* }						    */
  /**************************************************/


  for (i=0;i<(*ECEPE_p).NUMATM;++i) {
    for (j=0;j<4;++j) printf("%c",name_atom_list[i*4+j]);
    printf("--");
    for (j=0;j<num1_4[i];++j) 
      for (k=0;k<4;++k) printf("%c",name_atom_list[(pair1_4[i][j])*4+k]);
    printf("\n");
  }

  /*****************************************************************************/
  /* for (i=0;i<(*ECEPE_p).NUMATM;++i) {				       */
  /*   for (j=0;j<4;++j) printf("%c",name_atom_list[i*4+j]);		       */
  /*   printf("--");							       */
  /*   for (j=0;j<num1_5[i];++j) 					       */
  /*     for (k=0;k<4;++k) printf("%c",name_atom_list[(pairex1_5[i][j])*4+k]); */
  /*   printf("\n");							       */
  /* }									       */
  /*****************************************************************************/

  for (i=0;i<(*ECEPE_p).NUMATM;++i) {
    for (j=0;j<4;++j) printf("%c",name_atom_list[i*4+j]);
    printf("(%3d)--(%3d)\n",i,num1_5[i]);
  }

  numnb=0;
  num14=0;
  for (i=0;i<(*ECEPE_p).NUMATM;++i) numnb+=num1_5[i];
  for (i=0;i<(*ECEPE_p).NUMATM;++i) num14+=num1_4[i];

  /////////////////////////////////////////


  printf("nb list\n");
  for (i=0;i<(*ECEPE_p).NUMATM;++i) {
    printf("%d - ",i+1);
    for (j=0;j<num1_5[i];++j)
      printf("%d ",pair1_5[i][j]+1);
    printf("\n");
  }
  printf("14 list\n");
  for (i=0;i<(*ECEPE_p).NUMATM;++i) {
    printf("%d - ",i+1);
    for (j=0;j<num1_4[i];++j)
      printf("%d ",pair1_4[i][j]+1);
    printf("\n");
  }

  //  exit(1);

  //  co=(double *)gcemalloc(sizeof(double)*(*ECEPE_p).NUMATM*3); // 2015-02-16
  //  co=(double *)calloc((*ECEPE_p).NUMATM*3,sizeof(double));        // 2015-02-16

  coofile=efopen(coofilename,"r");
  read_ECEPE_coo(coofile,co,dihed,(*ECEPE_p).NUMATM);
  fclose(coofile);

  //////////////////////////////////////////////////////////////////////
  if ((*ECEPE_p).NUMATM>40) {
    //    co_dummy=(double *)gcemalloc(sizeof(double)*400*3); // 2015-02-16
    co_dummy=(double *)calloc(400*3,sizeof(double));          // 2015-02-16
  }
  else {
    //    co_dummy=(double *)gcemalloc(sizeof(double)*40*3);  // 2015-02-16
    co_dummy=(double *)calloc(40*3,sizeof(double));           // 2015-02-16
  }
  //  co2=(double *)gcemalloc(sizeof(double)*(*ECEPE_p).NUMATM*3); // 2015-02-16
  co2=(double *)calloc((*ECEPE_p).NUMATM*3,sizeof(double));        // 2015-02-16
  if ((*ECEPE_p).NUMVAR>10/*ECEPE_p->NUMVAR>10*/) {
    //    dihed_dummy=(double *)gcemalloc(sizeof(double)*100);    // 2015-02-16
    dihed_dummy=(double *)calloc(100,sizeof(double));             // 2015-02-16
  }
  else {
    //    dihed_dummy=(double *)gcemalloc(sizeof(double)*10);     // 2015-02-16
    dihed_dummy=(double *)calloc(10,sizeof(double));              // 2015-02-16
  }
  //  ene=(double *)gcemalloc(sizeof(double)*6);                  // 2015-02-16
  ene=(double *)calloc(6,sizeof(double));                         // 2015-02-16
  //  delta_dihed=(double *)gcemalloc(sizeof(double)*(ECEPE_p->NUMVAR));
  for (i=0;i<ECEPE_p->NUMVAR;++i) delta_dihed[i]=0.0;

  if (flagang==OFF) { // 0811
    coofile_for_sflag=efopen(coofilename_for_sflag,"r");
    if ((*ECEPE_p).NUMATM>40 ||(*ECEPE_p).NUMVAR>10 )
      read_ECEPE_detail_coo_cyc_for_protein(coofile_for_sflag,co_dummy,dihed_dummy,ene);
    else
      read_ECEPE_detail_coo_cyc(coofile_for_sflag,co_dummy,dihed_dummy,ene);
    for (i=0;i<ECEPE_p->NUMVAR;++i) {
      dihed_dummy[i]=dihed_dummy[i]*pi/180.0;
      if (dihed_dummy[i]<-pi)
	dihed_dummy[i]+=2.0*pi;
      else if (dihed_dummy[i]>pi)
	dihed_dummy[i]-=2.0*pi;
    }
    k=0;
    for (i=0;i<ECEPE_p->NUMATM;++i) {
      for (j=0;j<3;++j) {
	co2[i*3+j]=co_dummy[k];
	++k;
      }
    }    
    calc_TORS_for_get_sabun(ECEPE_p->NUMVAR,co2,*ECEPE_p,dihed_dummy,delta_dihed);
    fclose(coofile_for_sflag);
  } // 0811
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  // 0811
  if (flagang==ON) {
    //    numdiofres=(int *)gcemalloc(sizeof(int)*ECEPE_p->NUMRES); // 2015-02-16
    numdiofres=(int *)calloc(ECEPE_p->NUMRES,sizeof(int)); // 2015-02-16
    for (i=0;i<ECEPE_p->NUMRES;++i) numdiofres[i]=0;
    num=1;
    for (i=0;i<ECEPE_p->NUMVAR;++i) {
      if (num==(*ECEPE_p).dihed[i].indexv1) {
	numdiofres[num-1]+=1;
      }
      else {
	++num;
	numdiofres[num-1]+=1;
      }
    }
    angfile=efopen(angfilename,"r");
    k=0;
    for (i=0;i<ECEPE_p->NUMRES;++i) {
      for (j=0;j<numdiofres[i];++j) {
	fscanf(angfile,"%lf",&dihed_dummy[k]);
	++k;
      }
      for (j=0;j<10-numdiofres[i];++j)
	fscanf(angfile,"%lf",&fd);
    }
    fclose(angfile);
    ntotaldih=0;
    ndihinres=1;
    numres=1;
    //    dihed_dummy_dummy=(double *)gcemalloc(sizeof(double)*ECEPE_p->NUMVAR); // 2015-02-16
    dihed_dummy_dummy=(double *)calloc(ECEPE_p->NUMVAR,sizeof(double));       // 2015-02-16
    for (i=0;i<ECEPE_p->NUMVAR;++i) {
      if (ECEPE_p->dihed[i].indexv1>numres) {
	numres=ECEPE_p->dihed[i].indexv1;
	ntotaldih+=ndihinres;
	ndihinres=1;
      }
      if (ECEPE_p->dihed[i].indexv2>ndihinres) {
	ndihinres=ECEPE_p->dihed[i].indexv2;
      }
      dihed_dummy_dummy[i]=dihed_dummy[ntotaldih+ECEPE_p->dihed[i].indexv2-1];
    }
    for (i=0;i<ECEPE_p->NUMVAR;++i) {
      dihed_dummy[i]=dihed_dummy_dummy[i];
    }

    for (i=0;i<ECEPE_p->NUMVAR;++i) {
      dihed_dummy[i]=dihed_dummy[i]/180.0*pi;
      if (dihed_dummy[i]<-pi) dihed_dummy[i]+=2.0*pi;
      else if (dihed_dummy[i]>pi) dihed_dummy[i]-=2.0*pi;
    }

    for (i=0;i<ECEPE_p->NUMATM;++i)
      for (j=0;j<3;++j) co_dummy[(ECEPE_p->atom[i].katom-1)*3+j]=co[i*3+j];
    for (i=0;i<ECEPE_p->NUMATM;++i)
      for (j=0;j<3;++j)	co[i*3+j]=co_dummy[i*3+j];

    calc_TORS_for_get_sabun(ECEPE_p->NUMVAR,co,*ECEPE_p,dihed_dummy,delta_dihed);
  }
  // 0811
  //////////////////////////////////////////////////////////////////////



  /************************************************************/ // 2015-02-16
  /* for (i=0;i<(*ECEPE_p).NUMATM;++i) {		      */ // 2015-02-16
  /*   for (j=0;j<3;++j) {				      */ // 2015-02-16
  /*     prot.coord[(*ECEPE_p).atom[i].katom-1][j]=co[i*3+j]; */ // 2015-02-16
  /*   }						      */ // 2015-02-16
  /* }							      */ // 2015-02-16
  /************************************************************/ // 2015-02-16

  //////////////////////////////////////////
  if ((input=fopen(InpfilCOORD,"r")) == NULL) {
    printf("error %s cannot open\n", InpfilTOP);
    exit(1);
  }
  //  if (oliflag==ON) getline(&line,&len,input); // 2015-02-16
  fscanf(input,"%d",&i);  
  for(i=0;i<(*ECEPE_p).NUMATM;++i) {
    for(j=0;j<3;++j) {
      fscanf(input,"%lf",&x);
      //      prot.coord[i][j]=x;  // 2015-02-16
      co[i*3+j]=x;                 // 2015-02-16
    }
  }
  fclose(input);
  ///////////////////////////////////////////

  //  refcoordscan_ECEPE(prot.coord); // 2015-02-16
  //  refcoordscan_ECEPE(co); // 2015-02-16 // temporary comment out for debug

  free(numdiofres); // 2015-02-16


  //  return delta_dihed;
}

//int refcoordscan_ECEPE(double crd[MAXA][3]) { // 2015-02-16
int refcoordscan_ECEPE(double crd[MAXA][3],int numatom, int numdihed_rest, double **dihed_rest) { // 2015-02-16
  double x;
  int i,j,nNumDihed,nNumDihed_now;
  int flag = 0;
  int num = 0;
  double coord_rest[MAXA][3];

  double *theta_ref; // 2015-02-16

  FILE *log_rest;

  theta_ref=(double *)calloc(numdihed_rest,sizeof(double)); // 2015-02-16

  //  for(i=0;i<prot.num_atom;++i) { // 2015-02-16
  for(i=0;i<numatom;++i) {  // 2015-02-16
    for(j=0;j<3;++j) {
      coord_rest[i][j]=crd[i][j];
    }
    ++num;
  }


  //  for (nNumDihed=0;nNumDihed<prot.nNumDihed_rest;++nNumDihed) { // 2015-02-16
  for (nNumDihed=0;nNumDihed<numdihed_rest;++nNumDihed) { // 2015-02-16
    //    theta_ref[nNumDihed]/*rad*/ = pick_dihed_one_clust/*2*/(dihed_rest[nNumDihed][0]-1,dihed_rest[nNumDihed][1]-1,dihed_rest[nNumDihed][2]-1,dihed_rest[nNumDihed][3]-1,nNumDihed); // 2015-02-16 temporary comment out for debug
  }

  if ((log_rest=fopen("log_rest.txt","a"))==NULL) {
    printf("error: cannot open log_rest.txt\n");
    exit(1);
  }

  fprintf(log_rest,"log_rest.txt","a");
  //  for (i=0;i<prot.nNumDihed_rest;++i) // 2015-02-16
  for (i=0;i<numdihed_rest;++i)     // 2015-02-16
    fprintf(log_rest,"%e \n",theta_ref[i]);
  fclose(log_rest);

  free(theta_ref); // 2015-02-16
  
}

//void pick_mass_data(char *massfilename, struct ECEPE_parms ECEPE_p) { // 2015-02-16
void pick_mass_data(char *massfilename, struct ECEPE_parms ECEPE_p,int numclut, double *mass) { // 2015-02-16
  int i,j,num_a_prot;
  double *f;

  FILE *massfile;

  massfile=efopen(massfilename,"r");

  //  f=(double *)gcemalloc(sizeof(double)*ECEPE_p.NUMATM); // 2015-02-16
  f=(double *)calloc(ECEPE_p.NUMATM,sizeof(double));        // 2015-02-16

  for (i=0;i<ECEPE_p.NUMATM;++i) {
    fscanf(massfile,"%lf",&f[i]);
  }

  /***********************************************/ // 2015-02-16
  /* num_a_prot=0;				 */ // 2015-02-16
  /* for(i=0; i<prot.DOF; ++i){			 */ // 2015-02-16
  /*   for(j=0; j<clust[i].num_atom_clust; ++j){ */ // 2015-02-16
  /*     clust[i].mass_clust[j] = f[num_a_prot]; */ // 2015-02-16
  /*     ++num_a_prot;				 */ // 2015-02-16
  /*   }					 */ // 2015-02-16
  /* }						 */ // 2015-02-16
  /***********************************************/ // 2015-02-16

  for(i=0; i<ECEPE_p.NUMATM; ++i){  // 2015-02-16
    mass[i] = f[i];                 // 2015-02-16
  }                                 // 2015-02-16

  fclose(massfile);

  free(f); // 2015-02-16
}

double calc_TORS_for_get_sabun(int numdih, double *co, struct ECEPE_parms p ,double *dihed_inspidas, double *delta_dihed){
  int i,j;
  double dihed;
  double atom_i[3],atom_j[3],atom_k[3],atom_l[3];
  double pi;

  pi=acos(-1.0);
  
  for (i=0;i<numdih;++i) {
    for (j=0;j<3;++j) {
      atom_i[j]=co[p.dihed[i].dpairs[0]*3+j];
      atom_j[j]=co[p.dihed[i].dpairs[1]*3+j];
      atom_k[j]=co[p.dihed[i].dpairs[2]*3+j];
      atom_l[j]=co[p.dihed[i].dpairs[3]*3+j];
    }
    dihed=dih(atom_i,atom_j,atom_k,atom_l);
    if (dihed<-pi)
      dihed+=2.0*pi;
    else if (dihed>pi)
      dihed-=2.0*pi;
    delta_dihed[i]=dihed_inspidas[i]-dihed;
  }
}

void read_ECEPE_pairs(char *pairsfilename,int **pairs_1_5,int **pairs_1_4, int *num1_5, int *num1_4,int numatom) {
  int i,j;
  int d1,d2;
  int numtotal1_4,numtotal1_5;
  int na,nb;
  FILE *pairsfile;

  pairsfile=efopen(pairsfilename,"r");
  fscanf(pairsfile,"%d",&numtotal1_4);
  fscanf(pairsfile,"%d",&numtotal1_5);
  na=1;
  nb=0;
  for (i=0;i<numatom;++i) {
    num1_4[i]=0;
    num1_5[i]=0;
  }

  for (i=0;i<numtotal1_4;++i) {
    fscanf(pairsfile,"%d",&d1);
    fscanf(pairsfile,"%d",&d2);
    if (d1==na) {
      //      pairs_1_4[na-1]=(int *)gcerealloc(pairs_1_4[na-1],sizeof(int)*(nb+1)); // 2015-02-16
      pairs_1_4[na-1]=(int *)realloc(pairs_1_4[na-1],sizeof(int)*(nb+1)); // 2015-02-16
      pairs_1_4[na-1][nb]=d2-1;
      ++nb;
    }
    else {
      num1_4[na-1]=nb;
      na=d1;
      nb=1;
      //      pairs_1_4[na-1]=(int *)gcerealloc(pairs_1_4[na-1],sizeof(int)*1); // 2015-02-16
      pairs_1_4[na-1]=(int *)realloc(pairs_1_4[na-1],sizeof(int)*1); // 2015-02-16
      pairs_1_4[na-1][nb-1]=d2-1;
    }
  }
  num1_4[na-1]=nb;
  na=1;
  nb=0;
  for (i=0;i<numtotal1_5;++i) {
    fscanf(pairsfile,"%d",&d1);
    fscanf(pairsfile,"%d",&d2);
    if (d1==na) {
      //      pairs_1_5[na-1]=(int *)gcerealloc(pairs_1_5[na-1],sizeof(int)*(nb+1)); // 2015-02-16
      pairs_1_5[na-1]=(int *)realloc(pairs_1_5[na-1],sizeof(int)*(nb+1)); // 2015-02-16
      pairs_1_5[na-1][nb]=d2-1;
      ++nb;
    }
    else {
      num1_5[na-1]=nb;
      na=d1;
      nb=1;
      //      pairs_1_5[na-1]=(int *)gcerealloc(pairs_1_5[na-1],sizeof(int)*1); // 2015-02-16
      pairs_1_5[na-1]=(int *)realloc(pairs_1_5[na-1],sizeof(int)*1); // 2015-02-16
      pairs_1_5[na-1][nb-1]=d2-1;
    }
  }
  num1_5[na-1]=nb;
  fclose(pairsfile);
}

