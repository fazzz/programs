#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MolBas.h"
#include "EF.h"
#include "ParmTop.h"

double pick_bond_leng(  double atom_i[3],double atom_j[3]);
int scantraj(double *coord_tag,int numatom, FILE *inputfile);

int main(int argc, char *argv[]) {
  int i,j,k,numbond;
  int *atom_bond_pair,numatom,numstep;
  double atom_i[3],atom_j[3];
  double *protcoord;
  double bl;

  char *inputfilename,*outputfilename,*parmfilename,*condfilename;
  FILE *inputfile, *outputfile,*parmfile,*condfile;
  
  if (argc < 4) {
    printf("USAGE: ./CBL.exe inputfilename(coord) parmfile inputfile(cond) outputfilename\n");
    exit(1);
  }
  inputfilename =  *++argv;
  parmfilename = *++argv;
  condfilename = *++argv;
  outputfilename = *++argv;

  condfile=efopen(condfilename,"r");
  fscanf(condfile,"%d",&numstep);
  fscanf(condfile,"%d",&numbond);
  atom_bond_pair=emalloc(sizeof(int)*2*numbond);
  for (j=0;j<numbond;++j) {
    fscanf(condfile,"%d",&i);
    atom_bond_pair[j*2]=i-1;
    fscanf(condfile,"%d",&i);
    atom_bond_pair[j*4+1]=i-1;
  }
  fclose(condfile);

  parmfile=fopen(parmfilename,"r");
  readParmtop(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;
  protcoord=malloc(sizeof(double)*numatom*3);

  outputfile=efopen(outputfilename,"w");
  inputfile=efopen(inputfilename,"r");
  for (i=0;i<numstep;++i) {
    scantraj(protcoord,numatom,inputfile);
    for (j=0;j<numbond;++j) {
      for (k=0;k<3;++k) {
	atom_i[k]=protcoord[(atom_bond_pair[j*2])*3+k];
	atom_j[k]=protcoord[(atom_bond_pair[j*2+1])*3+k];
      }
      bl=pick_bond_leng(atom_i,atom_j);
      fprintf(outputfile,"%12.8lf ",bl);
    }
    fprintf(outputfile,"\n");
  }
  fclose(inputfile);
  fclose(outputfile);

}

int scantraj(double *coord_tag,int numatom, FILE *inputfile) {
  int i,j;
  double d;

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      if (fscanf(inputfile,"%lf",&d)== EOF){
	break;
      }
      coord_tag[i*3+j]=d;
    }
  }

  return 1;
}

double pick_bond_leng(  double atom_i[3],double atom_j[3]) {
  int alpha;
  double len=0.0;

  for (alpha=0;alpha<3;++alpha) {
    len += (atom_j[alpha]-atom_i[alpha])*(atom_j[alpha]-atom_i[alpha]);
  }
  len=sqrt(len);
  return len;
}
