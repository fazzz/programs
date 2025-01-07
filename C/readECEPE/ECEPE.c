#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//#include "glib.h"

#include "TOPO.h"

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
  "H3 ","C  ",
  "N  ","C  ",
  "NZ ","CE ",
  "NE ","CD ",
  "NE2","CD ",
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
  int c;
  int d;
  double f;
  char *line,dummy;
  size_t len=0;

  FILE *preo,*bd8;

  preo=efopen(preofilename,"r");
  bd8=efopen(bd8filename,"r");

  getline(&line,&len,preo);

  fscanf(preo,"%d",&(*ECEPE_p).NUMATM);
  fscanf(preo,"%d",&(*ECEPE_p).NUMVAR);
  fscanf(preo,"%d",&(*ECEPE_p).NUMRES);
  fscanf(preo,"%d",&(*ECEPE_p).NUMINT);
  fscanf(preo,"%d",&(*ECEPE_p).NUMS);

   (*ECEPE_p).dihed=(struct ECEPE_dihe *)gcemalloc(sizeof(struct ECEPE_dihe)*(*ECEPE_p).NUMVAR);
   (*ECEPE_p).atom=(struct ECEPE_atom *)gcemalloc(sizeof(struct ECEPE_atom)*(*ECEPE_p).NUMATM);

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
    /***************************************************/
    /* fscanf(preo,"%d",&(*ECEPE_p).dihed[i].indexv1); */
    /* fscanf(preo,"%d",&(*ECEPE_p).dihed[i].indexv2); */
    /* fscanf(preo,"%d",&(*ECEPE_p).dihed[i].ibnd1);   */
    /* fscanf(preo,"%d",&(*ECEPE_p).dihed[i].ibnd2);   */
    /* fscanf(preo,"%d",&(*ECEPE_p).dihed[i].ifront);  */
    /* fscanf(preo,"%d",&(*ECEPE_p).dihed[i].ibchar1); */
    /* fscanf(preo,"%d",&(*ECEPE_p).dihed[i].ibchar2); */
    /* fscanf(preo,"%d",&(*ECEPE_p).dihed[i].ibchar3); */
    /* fscanf(preo,"%lf",&(*ECEPE_p).dihed[i].A);      */
    /* fscanf(preo,"%d",&(*ECEPE_p).dihed[i].NB);      */
    /* fscanf(preo,"%d",&(*ECEPE_p).dihed[i].NS);      */
    /* fscanf(preo,"%d",&(*ECEPE_p).dihed[i].IFTOR);   */
    /* fscanf(preo,"%d",&d);			       */
    /* fscanf(preo,"%d",&d);			       */
    /* fscanf(preo,"%d",&d);			       */
    /* fscanf(preo,"%d",&d);			       */
    /* fscanf(preo,"%d",&d);			       */
    /***************************************************/
    (*ECEPE_p).dihed[i].indexv1=0;
    (*ECEPE_p).dihed[i].indexv2=0;
    (*ECEPE_p).dihed[i].ibnd1=0; 
    (*ECEPE_p).dihed[i].ibnd2=0; 
    (*ECEPE_p).dihed[i].ifront=0;
    (*ECEPE_p).dihed[i].ibchar1=0;
    (*ECEPE_p).dihed[i].ibchar2=0;
    (*ECEPE_p).dihed[i].ibchar3=0;
    c=getc(preo);
    for (j=0;j<4;++j) {
      c=getc(preo);
      if (c >= '0' && c <= '9') 
	(*ECEPE_p).dihed[i].indexv1=(*ECEPE_p).dihed[i].indexv1*10+c-'0';
    }
    for (j=0;j<2;++j) {
      c=getc(preo);
      if (c >= '0' && c <= '9') (*ECEPE_p).dihed[i].indexv2=(*ECEPE_p).dihed[i].indexv2*10+c-'0';
    }
    for (j=0;j<7;++j) {
      c=getc(preo);
      if (c >= '0' && c <= '9') (*ECEPE_p).dihed[i].ibnd1=(*ECEPE_p).dihed[i].ibnd1*10+c-'0';
    }
    for (j=0;j<5;++j) {
      c=getc(preo);
      if (c >= '0' && c <= '9') (*ECEPE_p).dihed[i].ibnd2=(*ECEPE_p).dihed[i].ibnd2*10+c-'0';
    }
    for (j=0;j<7;++j) {
      c=getc(preo);
      if (c >= '0' && c <= '9') (*ECEPE_p).dihed[i].ifront=(*ECEPE_p).dihed[i].ifront*10+c-'0';
    }
    for (j=0;j<5;++j) {
      c=getc(preo);
      if (c >= '0' && c <= '9') (*ECEPE_p).dihed[i].ibchar1=(*ECEPE_p).dihed[i].ibchar1*10+c-'0';
    }
    for (j=0;j<2;++j) {
      c=getc(preo);
      if (c >= '0' && c <= '9') (*ECEPE_p).dihed[i].ibchar2=(*ECEPE_p).dihed[i].ibchar2*10+c-'0';
    }
    for (j=0;j<2;++j) {
      c=getc(preo);
      if (c >= '0' && c <= '9') (*ECEPE_p).dihed[i].ibchar3=(*ECEPE_p).dihed[i].ibchar3*10+c-'0';
    }


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
      fscanf(preo,"%lf",&(*ECEPE_p).atom[i].refcoord[j]);
    fscanf(preo,"%lf",&(*ECEPE_p).atom[i].charge);
    fscanf(preo,"%d",&(*ECEPE_p).atom[i].nbtype);
    fscanf(preo,"%d",&(*ECEPE_p).atom[i].kunit);
    fscanf(preo,"%d",&(*ECEPE_p).atom[i].katom);
    fscanf(preo,"%d",&(*ECEPE_p).atom[i].jatom);
    for (j=0;j<2;++j)
      getc(preo);
    for (j=0;j<4;++j)
      (*ECEPE_p).atom[i].name_atom[j]=getc(preo);
    fscanf(preo,"%3s",&(*ECEPE_p).atom[i].name_res);
    fscanf(preo,"%d",&d);
  }

  (*nb_p).Acff=(double *)gcemalloc(sizeof(double)*numatomtype*numatomtype);
  (*nb_p).Bcff=(double *)gcemalloc(sizeof(double)*numatomtype*numatomtype);

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

}

void read_ECEPE_parm_wtransindex(char *preofilename, char *bd8filename, struct ECEPE_parms *ECEPE_p, struct pnb *nb_p) {
  int i,j;
  int c;
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

   (*ECEPE_p).dihed=(struct ECEPE_dihe *)gcemalloc(sizeof(struct ECEPE_dihe)*(*ECEPE_p).NUMVAR);
   (*ECEPE_p).atom=(struct ECEPE_atom *)gcemalloc(sizeof(struct ECEPE_atom)*(*ECEPE_p).NUMATM);
   atom_dummy=(struct ECEPE_atom *)gcemalloc(sizeof(struct ECEPE_atom)*(*ECEPE_p).NUMATM);

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
    /***************************************************/
    /* fscanf(preo,"%d",&(*ECEPE_p).dihed[i].indexv1); */
    /* fscanf(preo,"%d",&(*ECEPE_p).dihed[i].indexv2); */
    /* fscanf(preo,"%d",&(*ECEPE_p).dihed[i].ibnd1);   */
    /* fscanf(preo,"%d",&(*ECEPE_p).dihed[i].ibnd2);   */
    /* fscanf(preo,"%d",&(*ECEPE_p).dihed[i].ifront);  */
    /* fscanf(preo,"%d",&(*ECEPE_p).dihed[i].ibchar1); */
    /* fscanf(preo,"%d",&(*ECEPE_p).dihed[i].ibchar2); */
    /* fscanf(preo,"%d",&(*ECEPE_p).dihed[i].ibchar3); */
    /***************************************************/

    (*ECEPE_p).dihed[i].indexv1=0;
    (*ECEPE_p).dihed[i].indexv2=0;
    (*ECEPE_p).dihed[i].ibnd1=0; 
    (*ECEPE_p).dihed[i].ibnd2=0; 
    (*ECEPE_p).dihed[i].ifront=0;
    (*ECEPE_p).dihed[i].ibchar1=0;
    (*ECEPE_p).dihed[i].ibchar2=0;
    (*ECEPE_p).dihed[i].ibchar3=0;
    c=getc(preo);
    for (j=0;j<4;++j) {
      c=getc(preo);
      if (c >= '0' && c <= '9') 
	(*ECEPE_p).dihed[i].indexv1=(*ECEPE_p).dihed[i].indexv1*10+c-'0';
    }
    for (j=0;j<2;++j) {
      c=getc(preo);
      if (c >= '0' && c <= '9') (*ECEPE_p).dihed[i].indexv2=(*ECEPE_p).dihed[i].indexv2*10+c-'0';
    }
    for (j=0;j<7;++j) {
      c=getc(preo);
      if (c >= '0' && c <= '9') (*ECEPE_p).dihed[i].ibnd1=(*ECEPE_p).dihed[i].ibnd1*10+c-'0';
    }
    for (j=0;j<5;++j) {
      c=getc(preo);
      if (c >= '0' && c <= '9') (*ECEPE_p).dihed[i].ibnd2=(*ECEPE_p).dihed[i].ibnd2*10+c-'0';
    }
    for (j=0;j<7;++j) {
      c=getc(preo);
      if (c >= '0' && c <= '9') (*ECEPE_p).dihed[i].ifront=(*ECEPE_p).dihed[i].ifront*10+c-'0';
    }
    for (j=0;j<5;++j) {
      c=getc(preo);
      if (c >= '0' && c <= '9') (*ECEPE_p).dihed[i].ibchar1=(*ECEPE_p).dihed[i].ibchar1*10+c-'0';
    }
    for (j=0;j<2;++j) {
      c=getc(preo);
      if (c >= '0' && c <= '9') (*ECEPE_p).dihed[i].ibchar2=(*ECEPE_p).dihed[i].ibchar2*10+c-'0';
    }
    for (j=0;j<2;++j) {
      c=getc(preo);
      if (c >= '0' && c <= '9') (*ECEPE_p).dihed[i].ibchar3=(*ECEPE_p).dihed[i].ibchar3*10+c-'0';
    }

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

  (*nb_p).Acff=(double *)gcemalloc(sizeof(double)*numatomtype*numatomtype);
  (*nb_p).Bcff=(double *)gcemalloc(sizeof(double)*numatomtype*numatomtype);

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

}


void read_ECEPE_parm_wobd8(char *preofilename, struct ECEPE_parms *ECEPE_p) {
  int i,j;
  int d;
  double f;
  char *line,dummy;
  size_t len=0;

  FILE *preo;

  preo=efopen(preofilename,"r");

  getline(&line,&len,preo);

  fscanf(preo,"%d",&(*ECEPE_p).NUMATM);
  fscanf(preo,"%d",&(*ECEPE_p).NUMVAR);
  fscanf(preo,"%d",&(*ECEPE_p).NUMRES);
  fscanf(preo,"%d",&(*ECEPE_p).NUMINT);
  fscanf(preo,"%d",&(*ECEPE_p).NUMS);

   (*ECEPE_p).dihed=(struct ECEPE_dihe *)gcemalloc(sizeof(struct ECEPE_dihe)*(*ECEPE_p).NUMVAR);
   (*ECEPE_p).atom=(struct ECEPE_atom *)gcemalloc(sizeof(struct ECEPE_atom)*(*ECEPE_p).NUMATM);

   for (i=0;i<(*ECEPE_p).NUMVAR;++i) {
     fscanf(preo,"%lf",&(*ECEPE_p).dihed[i].angle);
     if((*ECEPE_p).dihed[i].angle==-75.000) 
       i-=1;
   }

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
      fscanf(preo,"%lf",&(*ECEPE_p).atom[i].refcoord[j]);
    fscanf(preo,"%lf",&(*ECEPE_p).atom[i].charge);
    fscanf(preo,"%d",&(*ECEPE_p).atom[i].nbtype);
    fscanf(preo,"%d",&(*ECEPE_p).atom[i].kunit);
    fscanf(preo,"%d",&(*ECEPE_p).atom[i].katom);
    fscanf(preo,"%d",&(*ECEPE_p).atom[i].jatom);
    for (j=0;j<2;++j)
      getc(preo);
    for (j=0;j<4;++j)
      (*ECEPE_p).atom[i].name_atom[j]=getc(preo);
    fscanf(preo,"%3s",&(*ECEPE_p).atom[i].name_res);
    fscanf(preo,"%d",&d);
  }

  fclose(preo);

}

double calc_ff_ECEPE(struct ECEPE_pote *p, struct ECEPE_force *f, struct ECEPE_parms ECEPE_p, struct pnb nb_p, int *pairs, int numint){
  int i,j,k;
  int numdih,numatom;
  double *A,*dihed;
  int *ns,*nb,*nbtype;
  int *iftors;
  double *charge;
  double p_nb,*f_nb,p_es,*f_es,t,*n,*co;

  numdih=ECEPE_p.NUMVAR;
  numatom=ECEPE_p.NUMATM;

  co=(double *)gcemalloc(sizeof(double)*numatom*3);
  A=(double *)gcemalloc(sizeof(double)*numdih);
  ns=(int *)gcemalloc(sizeof(int)*numdih);
  nb=(int *)gcemalloc(sizeof(int)*numdih);
  dihed=(double *)gcemalloc(sizeof(double)*numdih);
  iftors=(int *)gcemalloc(sizeof(int)*numdih);
  charge=(double *)gcemalloc(sizeof(double)*numatom);
  nbtype=(int *)gcemalloc(sizeof(int)*numatom);

  f_nb=(double *)gcemalloc(sizeof(double)*numatom*3);
  f_es=(double *)gcemalloc(sizeof(double)*numatom*3);
  n=(double *)gcemalloc(sizeof(double)*numdih);

  for (i=0;i<numdih;++i) {
    A[i]=ECEPE_p.dihed[i].A;
    ns[i]=ECEPE_p.dihed[i].NS;
    nb[i]=ECEPE_p.dihed[i].NB;
    dihed[i]=ECEPE_p.dihed[i].angle;
    iftors[i]=ECEPE_p.dihed[i].IFTOR;
  }

  for (i=0;i<numatom;++i) {
    charge[i]=ECEPE_p.atom[i].charge;
    nbtype[i]=ECEPE_p.atom[i].nbtype;
    for (j=0;j<3;++j)
      co[i*3+j]=ECEPE_p.atom[i].refcoord[j];
  }

  //  set_pairs_ECEPE();
  //  calc_TORS_ECEPE(&t,n,A,ns,nb,dihed,iftors,numdih);
  calc_TORS_ECEPE2(&t,n,A,ns,nb,iftors,numdih,co,ECEPE_p);
  calc_NB_ECEPE(&p_nb,f_nb,&p_es,f_es,co,pairs,numint,charge,nbtype,nb_p);

  (*p).p_t=t+p_nb+p_es;
  (*p).p_tors=t;
  (*p).p_nb=p_nb;
  (*p).p_es=p_es;

}

double calc_ff_ECEPE_for_db(struct ECEPE_pote *p, struct ECEPE_force *f, struct ECEPE_parms ECEPE_p, struct pnb nb_p, int **pairs_1_5,int **pairs_1_4, int *num1_5, int *num1_4){
  int i,j,k;
  int numdih,numatom;
  double *A,*dihed;
  int *ns,*nb,*nbtype;
  int *iftors;
  double *charge;
  double p_nb,*f_nb,p_es,*f_es,t,*n,*co;

  numdih=ECEPE_p.NUMVAR;
  numatom=ECEPE_p.NUMATM;

  co=(double *)gcemalloc(sizeof(double)*numatom*3);
  A=(double *)gcemalloc(sizeof(double)*numdih);
  ns=(int *)gcemalloc(sizeof(int)*numdih);
  nb=(int *)gcemalloc(sizeof(int)*numdih);
  dihed=(double *)gcemalloc(sizeof(double)*numdih);
  iftors=(int *)gcemalloc(sizeof(int)*numdih);
  charge=(double *)gcemalloc(sizeof(double)*numatom);
  nbtype=(int *)gcemalloc(sizeof(int)*numatom);

  f_nb=(double *)gcemalloc(sizeof(double)*numatom*3);
  f_es=(double *)gcemalloc(sizeof(double)*numatom*3);
  n=(double *)gcemalloc(sizeof(double)*numdih);

  for (i=0;i<numdih;++i) {
    A[i]=ECEPE_p.dihed[i].A;
    ns[i]=ECEPE_p.dihed[i].NS;
    nb[i]=ECEPE_p.dihed[i].NB;
    dihed[i]=ECEPE_p.dihed[i].angle;
    iftors[i]=ECEPE_p.dihed[i].IFTOR;
  }

  for (i=0;i<numatom;++i) {
    charge[i]=ECEPE_p.atom[i].charge;
    nbtype[i]=ECEPE_p.atom[i].nbtype;
    for (j=0;j<3;++j)
      co[i*3+j]=ECEPE_p.atom[i].refcoord[j];
  }

  //  set_pairs_ECEPE();
  //  calc_TORS_ECEPE(&t,n,A,ns,nb,dihed,iftors,numdih);
  calc_TORS_ECEPE2(&t,n,A,ns,nb,iftors,numdih,co,ECEPE_p);
  calc_NB_ECEPE_for_db(&p_nb,f_nb,&p_es,f_es,co,pairs_1_5,pairs_1_4,num1_5,num1_4,numatom,charge,nbtype,nb_p);

  (*p).p_t=t+p_nb+p_es;
  (*p).p_tors=t;
  (*p).p_nb=p_nb;
  (*p).p_es=p_es;

}

double calc_ff_ECEPE_for_db_cyc(struct ECEPE_pote *p, struct ECEPE_force *f, struct ECEPE_parms ECEPE_p, struct pnb nb_p, int **pairs_1_5,int **pairs_1_4, int *num1_5, int *num1_4, FILE *file, double *delta_dihed){
  int i,j,k;
  int numdih,numatom;
  double *A,*dihed;
  int *ns,*nb,*nbtype;
  int *iftors;
  double *charge;
  double p_nb,*f_nb,p_es,*f_es,t,*n,*co;

  numdih=ECEPE_p.NUMVAR;
  numatom=ECEPE_p.NUMATM;

  co=(double *)gcemalloc(sizeof(double)*numatom*3);
  A=(double *)gcemalloc(sizeof(double)*numdih);
  ns=(int *)gcemalloc(sizeof(int)*numdih);
  nb=(int *)gcemalloc(sizeof(int)*numdih);
  dihed=(double *)gcemalloc(sizeof(double)*numdih);
  iftors=(int *)gcemalloc(sizeof(int)*numdih);
  charge=(double *)gcemalloc(sizeof(double)*numatom);
  nbtype=(int *)gcemalloc(sizeof(int)*numatom);

  f_nb=(double *)gcemalloc(sizeof(double)*numatom*3);
  f_es=(double *)gcemalloc(sizeof(double)*numatom*3);
  n=(double *)gcemalloc(sizeof(double)*numdih);

  for (i=0;i<numdih;++i) {
    A[i]=ECEPE_p.dihed[i].A;
    ns[i]=ECEPE_p.dihed[i].NS;
    nb[i]=ECEPE_p.dihed[i].NB;
    dihed[i]=ECEPE_p.dihed[i].angle;
    iftors[i]=ECEPE_p.dihed[i].IFTOR;
  }

  for (i=0;i<numatom;++i) {
    charge[i]=ECEPE_p.atom[i].charge;
    nbtype[i]=ECEPE_p.atom[i].nbtype;
    for (j=0;j<3;++j)
      co[i*3+j]=ECEPE_p.atom[i].refcoord[j];
  }

  //  set_pairs_ECEPE();
  //  calc_TORS_ECEPE(&t,n,A,ns,nb,dihed,iftors,numdih);
  calc_TORS_ECEPE2_for_check(&t,n,A,ns,nb,iftors,numdih,co,ECEPE_p,file,delta_dihed);
  calc_NB_ECEPE_for_db(&p_nb,f_nb,&p_es,f_es,co,pairs_1_5,pairs_1_4,num1_5,num1_4,numatom,charge,nbtype,nb_p);

  (*p).p_t=t+p_nb+p_es;
  (*p).p_tors=t;
  (*p).p_nb=p_nb;
  (*p).p_es=p_es;

}


/**************************************************************************************************************************************************************************************************************************/
/* double calc_ff_ECEPE_for_db_cyc(struct ECEPE_pote *p, struct ECEPE_force *f, struct ECEPE_parms ECEPE_p, struct pnb nb_p, int **pairs_1_5,int **pairs_1_4, int *num1_5, int *num1_4, FILE *file, double *delta_dihed){ */
/*   int i,j,k;																										  */
/*   int numdih,numatom;																								  */
/*   double *A,*dihed;																									  */
/*   int *ns,*nb,*nbtype;																								  */
/*   int *iftors;																									  */
/*   double *charge;																									  */
/*   double p_nb,*f_nb,p_es,*f_es,t,*n,*co;																						  */
/* 																											  */
/*   numdih=ECEPE_p.NUMVAR;																								  */
/*   numatom=ECEPE_p.NUMATM;																								  */
/* 																											  */
/*   co=(double *)gcemalloc(sizeof(double)*numatom*3);																					  */
/*   A=(double *)gcemalloc(sizeof(double)*numdih);																					  */
/*   ns=(int *)gcemalloc(sizeof(int)*numdih);																						  */
/*   nb=(int *)gcemalloc(sizeof(int)*numdih);																						  */
/*   dihed=(double *)gcemalloc(sizeof(double)*numdih);																					  */
/*   iftors=(int *)gcemalloc(sizeof(int)*numdih);																					  */
/*   charge=(double *)gcemalloc(sizeof(double)*numatom);																				  */
/*   nbtype=(int *)gcemalloc(sizeof(int)*numatom);																					  */
/* 																											  */
/*   f_nb=(double *)gcemalloc(sizeof(double)*numatom*3);																				  */
/*   f_es=(double *)gcemalloc(sizeof(double)*numatom*3);																				  */
/*   n=(double *)gcemalloc(sizeof(double)*numdih);																					  */
/* 																											  */
/*   for (i=0;i<numdih;++i) {																								  */
/*     A[i]=ECEPE_p.dihed[i].A;																								  */
/*     ns[i]=ECEPE_p.dihed[i].NS;																							  */
/*     nb[i]=ECEPE_p.dihed[i].NB;																							  */
/*     dihed[i]=ECEPE_p.dihed[i].angle;																							  */
/*     iftors[i]=ECEPE_p.dihed[i].IFTOR;																						  */
/*   }																											  */
/* 																											  */
/*   for (i=0;i<numatom;++i) {																								  */
/*     charge[i]=ECEPE_p.atom[i].charge;																						  */
/*     nbtype[i]=ECEPE_p.atom[i].nbtype;																						  */
/*     for (j=0;j<3;++j)																								  */
/*       co[i*3+j]=ECEPE_p.atom[i].refcoord[j];																						  */
/*   }																											  */
/* 																											  */
/*   //  set_pairs_ECEPE();																								  */
/*   //  calc_TORS_ECEPE(&t,n,A,ns,nb,dihed,iftors,numdih);																				  */
/*   calc_TORS_ECEPE2_for_check(&t,n,A,ns,nb,iftors,numdih,co,ECEPE_p,file,delta_dihed);																  */
/*   calc_NB_ECEPE_for_db(&p_nb,f_nb,&p_es,f_es,co,pairs_1_5,pairs_1_4,num1_5,num1_4,numatom,charge,nbtype,nb_p);													  */
/* 																											  */
/*   (*p).p_t=t+p_nb+p_es;																								  */
/*   (*p).p_tors=t;																									  */
/*   (*p).p_nb=p_nb;																									  */
/*   (*p).p_es=p_es;																									  */
/* 																											  */
/* }																											  */
/**************************************************************************************************************************************************************************************************************************/


double calc_TORS_ECEPE(double *t,double *n,double *A,int *ns,int *nb,double *dihed, int *iftors,int numdih ){
  int i;
  
  *t=0.0;

  for (i=0;i<numdih;++i) {
    if (iftors[i]!=0) {
      *t+=A[i]*(1.0+(double)ns[i]*cos((double)nb[i]*dihed[i]));
      n[i]=-A[i]*((double)nb[i]*(double)ns[i]*sin((double)nb[i]*dihed[i]))*4.184070*100;
    }
  }
}

double calc_TORS_ECEPE2(double *t,double *n,double *A,int *ns,int *nb, int *iftors,int numdih, double *co, struct ECEPE_parms p ){
  int i,j;
  double dihed;
  double atom_i[3],atom_j[3],atom_k[3],atom_l[3];
  
  *t=0.0;

  for (i=0;i<numdih;++i) {
    for (j=0;j<3;++j) {
      atom_i[j]=co[p.dihed[i].dpairs[0]*3+j];
      atom_j[j]=co[p.dihed[i].dpairs[1]*3+j];
      atom_k[j]=co[p.dihed[i].dpairs[2]*3+j];
      atom_l[j]=co[p.dihed[i].dpairs[3]*3+j];
    }
    dihed=dih(atom_i,atom_j,atom_k,atom_l);
    if (iftors[i]!=0) {
      *t+=A[i]*(1.0+(double)ns[i]*cos((double)nb[i]*dihed));
      n[i]=-A[i]*((double)nb[i]*(double)ns[i]*sin((double)nb[i]*dihed))*4.184070*100;
    }
  }
}

double calc_TORS_ECEPE2_for_check(double *t,double *n,double *A,int *ns,int *nb, int *iftors,int numdih, double *co, struct ECEPE_parms p ,FILE *file, double *delta_dihed){
  int i,j;
  double dihed;
  double atom_i[3],atom_j[3],atom_k[3],atom_l[3];
  double pi;

  pi=acos(-1.0);
  
  *t=0.0;

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
    dihed+=delta_dihed[i];
    fprintf(file,"%16.10lf ",dihed*180.0/pi);
    if (iftors[i]!=0) {
      *t+=A[i]*(1.0+(double)ns[i]*cos((double)nb[i]*dihed));
      n[i]=-A[i]*((double)nb[i]*(double)ns[i]*sin((double)nb[i]*dihed))*4.184070*100;
    }
  }

  fprintf(file,"\n");
}

double calc_TORS_ECEPE2_for_check_out(int numdih, double *co,struct ECEPE_parms p ,FILE *file, double *delta_dihed){
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
    dihed+=delta_dihed[i];
    if (dihed<-pi)
      dihed+=2.0*pi;
    else if (dihed>pi)
      dihed-=2.0*pi;
    fprintf(file,"%d %16.10lf \n",i+1,dihed*180.0/pi);
  }

  fprintf(file,"\n");
}

double chang_index_dihed_pairs(int numdih, struct ECEPE_parms p ){
  int i,j;
  int **dummy;

  dummy=(int**)gcemalloc(sizeof(int*)*numdih);
  for (i=0;i<numdih;++i)
    dummy[i]=(int*)gcemalloc(sizeof(int)*4);

  for (i=0;i<numdih;++i) {
    for (j=0;j<4;++j) {
      dummy[(p.dihed[i].indexv1-1)+p.dihed[i].indexv2-1][j]=p.dihed[i].dpairs[j];
    }
  }

  for (i=0;i<numdih;++i) {
    for (j=0;j<4;++j) {
      p.dihed[i].dpairs[j]=dummy[i][j];
    }
  }

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


double calc_TORS_for_check(double *co, int atom1, int atom2, int atom3, int atom4){
  int i,j;
  double dihed;
  double atom_i[3],atom_j[3],atom_k[3],atom_l[3];
  
  for (j=0;j<3;++j) {
    atom_i[j]=co[atom1*3+j];
    atom_j[j]=co[atom2*3+j];
    atom_k[j]=co[atom3*3+j];
    atom_l[j]=co[atom4*3+j];
  }
  dihed=dih(atom_i,atom_j,atom_k,atom_l);

  return dihed;
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

    // printf("%d\n",indxpb);
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

  FILE *ele_sta_ASN;

  ele_sta_ASN=efopen("debug.txt","w");

  *p_es=0.0;
  *p_nb=0.0;
  for (i=0;i<numatom;++i) {
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
	f=nb_p.Acff[indxpb]*ro6*(0.5*ro6 - 2.0);
	//	printf ("%d %d %10.4e 1-4nb\n",indxa,indxb,f);
	*p_nb += nb_p.Acff[indxpb]*ro6*(0.5*ro6 - 2.0);
	for (l=0;l<3;++l)
	  f_nb[i*3+l] += -nb_p.Acff[indxpb]*ro6*(0.5*ro6*24.0 - 2.0*12.0)/len2*vec[j]*4.184070*100;
      }
      else {
	ro4=ro2*ro2;
	f=-nb_p.Acff[indxpb]*ro6*(ro6 - 2.0*ro4);
	//	printf ("%d %d %10.4e 1-4nb\n",indxa,indxb,f);
	*p_nb -= nb_p.Acff[indxpb]*ro6*(ro6 - 2.0*ro4);
	for (l=0;l<3;++l)
	  f_nb[i*3+l] += -nb_p.Acff[indxpb]*ro6*(ro6*24.0 - 2.0*ro4*20.0)/len2*vec[j]*4.184070*100;
      }
    }

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
      fprintf(ele_sta_ASN,"%d %d %10.4e 1-4es \nqa=%10.4e qb=%10.4e len2=%10.4e\n",indxa+1,indxb+1,f,charge[indxa],charge[indxb],len2);
      *p_es+=charge[indxa]*charge[indxb]/sqrt(len2);
      for (l=0;l<3;++l)
	f_es[i*3+l] += charge[indxa]*charge[indxb]/len2*vec[j]*4.184070*100.0;
    }

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
	f=nb_p.Acff[indxpb]*ro6*(ro6 - 2.0);
	//	printf ("%d %d %10.4e 1-5nb\n",indxa,indxb,f);
	*p_nb +=  nb_p.Acff[indxpb]*ro6*(ro6 - 2.0);
	for (l=0;l<3;++l)
	  f_nb[i*3+l] += -nb_p.Acff[indxpb]*ro6*(ro6*24.0 - 2.0*12.0)/len2*vec[j]*4.184070*100;
      }
      else {
	ro4=ro2*ro2;
	f=-nb_p.Acff[indxpb]*ro6*(ro6 - 2.0*ro4);
	//	printf ("%d %d %10.4e 1-5nb\n",indxa,indxb,f);
	*p_nb -= nb_p.Acff[indxpb]*ro6*(ro6 - 2.0*ro4);
	for (l=0;l<3;++l)
	  f_nb[i*3+l] += -nb_p.Acff[indxpb]*ro6*(ro6*24.0 - 2.0*ro4*20.0)/len2*vec[j]*4.184070*100;
      }
    }

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
      fprintf(ele_sta_ASN,"%d %d %10.4e 1-5es \nqa=%10.4e qb=%10.4e len2=%10.4e\n",indxa+1,indxb+1,f,charge[indxa],charge[indxb],len2);
      *p_es+=charge[indxa]*charge[indxb]/sqrt(len2);
      for (l=0;l<3;++l)
	f_es[i*3+l] += charge[indxa]*charge[indxb]/len2*vec[j]*4.184070*100.0;
      //      printf ("%d %d %10.4e 1-5es\n",indxa,indxb,f);
    }
  }

  fclose(ele_sta_ASN);

}

int set_pairs_ECEPE(int *pairs,int numint,int *nbtype) {

}

void read_ECEPE_coo (FILE *file, double *co, double *dihed, int numatom) {
  int i,j,k,nd,na;
  int flag;
  double diheddummy;
  char *line,x[6],y[6],z[6],d[3],khi[4];
  size_t len=0;

  nd=0;
  na=0;
  i=0;
  for (;;) {
    getline(&line,&len,file);
    if (line[0]=='4') {
      flag=ON;
      for (;nd<2;) {
	for (j=0;j<8;++j) {
	  d[j]=line[48+nd*8+j];
	}
	dihed[nd]=atof(d);
	++nd;
      }
      for (j=0;j<8;++j) {
	d[j]=line[48+nd*8+j];
      }
      diheddummy=atof(d);
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
  dihed[nd]=diheddummy;
  ++nd;

  /************************************/
  /* for (j=0;j<i;++j) {	      */
  /*   for (k=0;k<3;++k)	      */
  /*     printf("%6.4lf ",co[j*3+k]); */
  /*   printf("\n");		      */
  /* }				      */
  /************************************/

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

int make_int_pair_list(int **bp,int *numb,int numatom,int numbond, int **pair1_5, int **pair1_4, int *num14,int *num1_5) {
  int i,j,k,l,n,m;
  int temp;
  int numatom_1,numatom_2,numatom_3,numatom_4;
  int *numex15;
  int **pair_ex_1_5;
  
  pair_ex_1_5=(int **)gcemalloc(sizeof(int *)*numatom);
  //  num14=(int *)gcemalloc(sizeof(int)*numatom);
  //  num1_5=(int *)gcemalloc(sizeof(int)*numatom);
  numex15=(int *)gcemalloc(sizeof(int)*numatom);
  
  for (i=0;i<numatom;++i) {
    num14[i]=0;
    numex15[i]=0;
    for (j=0;j<numb[i];++j) {
      numatom_1=bp[i][j];
      ++numex15[i];
      pair_ex_1_5[i]=(int *)gcerealloc(pair_ex_1_5[i],sizeof(int)*numex15[i]);
      pair_ex_1_5[i][numex15[i]-1]=numatom_1;
      for (k=0;k<numb[numatom_1];++k) {
	numatom_2=bp[numatom_1][k];
	if (numatom_2!=i) {
	  ++numex15[i];
	  pair_ex_1_5[i]=(int *)gcerealloc(pair_ex_1_5[i],sizeof(int)*numex15[i]);
	  pair_ex_1_5[i][numex15[i]-1]=numatom_2;
	  for (l=0;l<numb[numatom_2];++l) {
	    numatom_3=bp[numatom_2][l];
	    if (numatom_3!=numatom_1) {
	      ++numex15[i];
	      pair_ex_1_5[i]=(int *)gcerealloc(pair_ex_1_5[i],sizeof(int)*numex15[i]);
	      pair_ex_1_5[i][numex15[i]-1]=numatom_3;
	      if (numatom_3 > i) {
		++num14[i];
		pair1_4[i]=(int *)gcerealloc(pair1_4[i],sizeof(int)*num14[i]);
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
  	pair1_5[i]=(int *)gcerealloc(pair1_5[i],sizeof(int)*num1_5[i]);
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
    a2=p.dihed[i].ibnd1-1;
    a3=p.dihed[i].ibnd2-1;
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

void make_dihed_pairs_list_v2( struct ECEPE_parms p, int **bpl, int *numb ) {
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

  /***************************************************************/
  /* if (flag==OFF)						 */
  /*   printf("error: there is no candidate to bond %s\n",name); */
  /* else							 */
  /*   printf("%s %d to %s %d\n",name,nb,name_p,i);		 */
  /***************************************************************/

}

int order_bd_pair_list(int *bp, int **bp_f, int numatom, int *numb) {
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

/********************************************************************/
/* void print_iter(gpointer key, gpointer val, gpointer usr_data) { */
/*   g_debug("%s:%d",(gchar *)key,*(gint *)val);		    */
/* }								    */
/********************************************************************/

int get_sabun_from_delta_dihedfile(FILE *deltadihedfile,int numvar,struct ECEPE_parms ECEPE_p,double *delta_dihed){
  int i,j,k;
  char *sequence;

  /********************************************************************************/
  /* delta_dihedlist=g_hash_table_new_full(g_str_hash,g_str_equal,g_free,g_free); */
  /* 										  */
  /* for (i=0;i<NUMTOTALRES;++i) {						  */
  /*   atomnameaskey=g_strdup(atommasstypelist[i]);				  */
  /*   atommass=g_new(gchar,1);							  */
  /*   atommass=atommasslist[i];						  */
  /*   g_hash_table_insert(masslist,atomnameaskey,atommass);			  */
  /*   if ((mass=g_hash_table_lookup(masslist,atommasstypelist[i]))==NULL) {	  */
  /*     printf("There is no key %s in hashtable",atommasstypelist[i]);		  */
  /*     exit(1);								  */
  /*   } 									  */
  /* }										  */
  /* 										  */
  /* delta_dihed[0]=-60.0/180*pi;						  */
  /* 										  */
  /* sequence=get_sequence(ECEPE_p);						  */
  /* 										  */
  /* k=1;									  */
  /* for (i=0;i<ECEPE_p.NUMRES;++i) {						  */
  /*   resname=sequence[i];							  */
  /*   numatomonres=;								  */
  /*   for (j=0;j<numatomonres;++j) {						  */
  /*     delta_dihed[k]=resang[][j];						  */
  /*     ++k;									  */
  /*   }									  */
  /* }										  */
  /* 										  */
  /* delta_dihed[0]=-60.0/180*pi;						  */
  /********************************************************************************/
}
