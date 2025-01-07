#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "quaternion.h"
#include "bestfit.h"

#include "f2c.h"
#include "clapack.h"

int  scancoord(double coord[MAXNUMATOM][3],char *filename);
void scanmass(double mass[MAXNUMATOM],char *filename, int numatom);
int  scantraj(double coord_tag[MAXNUMATOM][3],int numatom, FILE *inputfile);
void test_velo(double coord_tag[MAXNUMATOM][3],double coord_tag_old[MAXNUMATOM][3],double velo_tag[MAXNUMATOM][3],int numatom, double deltat);

int main(int argc, char *argv[]){
  int i,j,k,s,numatom;
  double coord_ref[MAXNUMATOM][3],coord_tag[MAXNUMATOM][3],velo_tag[MAXNUMATOM][3],coord_bestfit[3],coord_tag_old[MAXNUMATOM][3],test[3];
  double mass[MAXNUMATOM];
  double rmsd=0.0;
  double omega[3],angMom[3];
  //  long int info,lwork=16,n=4;
  char *inputfilename,*inputfilename2, *outputfilename, *outputfilename2;
  FILE *inputfile, *inputfile2,*outputfile,*outputfile2;

  static long int m=3,n=3,lda=3,info,piv[3],lwork=3;
  static double work[3];
  double Inertia[9], Inertia_dummy[3][3];
  char *line;
  size_t len=0;

  numatom=scancoord(coord_ref,"crd_ref.txt");
  scanmass(mass,"mass.txt",numatom);

  if (argc < 5)
  {
    printf("USAGE:./bestfit.exe input(traj) input(velo) output(traj) output(velo)\n");
    exit(1);
  }
  else
  {
    inputfilename = *++argv;
    inputfilename2 = *++argv;
    outputfilename = *++argv;
    outputfilename2 = *++argv;
  }

  if((inputfile=fopen(inputfilename,"r"))==NULL)
  {
    printf("There is not %s\n",inputfilename);
    exit(1);
  }

  if((inputfile2=fopen(inputfilename2,"r"))==NULL)
  {
    printf("There is not %s\n",inputfilename2);
    exit(1);
  }

  if((outputfile=fopen(outputfilename,"w"))==NULL)
  {
     printf("There is not %s\n",outputfilename);
     exit(1);
  }

  if((outputfile2=fopen(outputfilename2,"w"))==NULL)
  {
    printf("There is not %s\n",outputfilename2);
    exit(1);
  }

  for (i=0;i<3;++i)
    test[i] = 0.0;

  getline(&line,&len,inputfile);
  getline(&line,&len,inputfile2);

  s=0;
  for (;;)
  {
    ++s;

    if(scantraj(coord_tag,numatom,inputfile)==0)
    {
      break;
    }
    scantraj(velo_tag,numatom,inputfile2);
  
    rmsd=bestfit(coord_ref,coord_tag,velo_tag,mass,numatom);

    //    test_velo(coord_tag,coord_tag_old,velo_tag,numatom,1.0);

    for (i=0;i<3;++i)
      angMom[i] = 0.0;
    for (i=0;i<numatom;++i){
      angMom[0]+= mass[i]*(coord_tag[i][2]*velo_tag[i][1]-coord_tag[i][1]*velo_tag[i][2]);
      angMom[1]+= mass[i]*(coord_tag[i][0]*velo_tag[i][2]-coord_tag[i][2]*velo_tag[i][0]);
      angMom[2]+= mass[i]*(coord_tag[i][1]*velo_tag[i][0]-coord_tag[i][0]*velo_tag[i][1]);
    }

    for (i=0;i<9;++i)
      Inertia[i]=0.0;
    for (i=0;i<numatom;++i){
      Inertia[0]+=mass[i]*coord_tag[i][1]*coord_tag[i][1]+mass[i]*coord_tag[i][2]*coord_tag[i][2];
      Inertia[1]-=mass[i]*coord_tag[i][0]*coord_tag[i][1];
      Inertia[2]-=mass[i]*coord_tag[i][0]*coord_tag[i][2];

      Inertia[3]-=mass[i]*coord_tag[i][1]*coord_tag[i][0];
      Inertia[4]+=mass[i]*coord_tag[i][0]*coord_tag[i][0]+mass[i]*coord_tag[i][2]*coord_tag[i][2];
      Inertia[5]-=mass[i]*coord_tag[i][1]*coord_tag[i][2];

      Inertia[6]-=mass[i]*coord_tag[i][2]*coord_tag[i][0];
      Inertia[7]-=mass[i]*coord_tag[i][2]*coord_tag[i][1];
      Inertia[8]+=mass[i]*coord_tag[i][0]*coord_tag[i][0]+mass[i]*coord_tag[i][1]*coord_tag[i][1];
	}

    dgetrf_(&m,&n,Inertia,&lda,piv,&info);

    dgetri_(&n,Inertia,&lda,piv,work,&lwork,&info);

    k=0;
    for (i=0;i<3;++i){
      for (j=0;j<3;++j){
	Inertia_dummy[i][j]=Inertia[k];
	++k;
      }
    }

    for (i=0;i<3;++i)
      omega[i]=0.0;
    for (i=0;i<3;++i){
      for (j=0;j<3;++j){
	omega[i]+=Inertia_dummy[i][j]*angMom[j];
      }
    }

  for (i=0;i<numatom;++i)
  {
    velo_tag[i][0] -= coord_tag[i][1]*omega[2]-coord_tag[i][2]*omega[1];
    velo_tag[i][1] -= coord_tag[i][2]*omega[0]-coord_tag[i][0]*omega[2];
    velo_tag[i][2] -= coord_tag[i][0]*omega[1]-coord_tag[i][1]*omega[0];
  }

  for (i=0;i<3;++i)
    angMom[i] = 0.0;
  for (i=0;i<numatom;++i){
    angMom[0]+= mass[i]*(coord_tag[i][2]*velo_tag[i][1]-coord_tag[i][1]*velo_tag[i][2]);
    angMom[1]+= mass[i]*(coord_tag[i][0]*velo_tag[i][2]-coord_tag[i][2]*velo_tag[i][0]);
    angMom[2]+= mass[i]*(coord_tag[i][1]*velo_tag[i][0]-coord_tag[i][0]*velo_tag[i][1]);
  }

    for (i=0;i<numatom;++i)
      for (j=0;j<3;++j)
    coord_tag_old[i][j]=coord_tag[i][j];

    for (i=0;i<numatom;++i)
     fprintf(outputfile,"%12.8lf %12.8lf %12.8lf\n",coord_tag[i][0],coord_tag[i][1],coord_tag[i][2]);

    fprintf(outputfile,"\n\n");

    for (i=0;i<numatom;++i)
     fprintf(outputfile2,"%12.8lf %12.8lf %12.8lf\n",velo_tag[i][0],velo_tag[i][1],velo_tag[i][2]);

    fprintf(outputfile2,"\n\n");

    //    printf("%d %lf %lf\n",s,rmsd,rmsd2);
    printf("%d %lf \n",s,rmsd);
  }

  fprintf(outputfile,"111111");
  fclose(outputfile);
  fclose(outputfile2);

  fclose(inputfile);
   return 0;
}

int scancoord(double coord[MAXNUMATOM][3],char *filename)
{
  int i,j,numatom;
  double d;
  FILE *file;

  if ((file=fopen(filename,"r"))== NULL )
    exit(1);

  fscanf(file,"%d",&i);
  numatom=i;

  for (i=0;i<numatom;++i)
  {
    for (j=0;j<3;++j)
    {
      fscanf(file,"%lf",&d);
      coord[i][j]=d;
    }
  }

  fclose(file);

  return numatom;

}

void scanmass(double mass[MAXNUMATOM],char *filename, int numatom)
{
  int i,j;
  double d;
  FILE *file;

  if ((file=fopen(filename,"r"))== NULL )
    exit(1);

  for (i=0;i<numatom;++i)
  {
     fscanf(file,"%lf",&d);
     mass[i]=d;
  }

  fclose(file);

}

int scantraj(double coord_tag[MAXNUMATOM][3],int numatom, FILE *inputfile)
{
  int i,j;
  double d;

  for (i=0;i<numatom;++i)
  {
    for (j=0;j<3;++j)
    {
      fscanf(inputfile,"%lf",&d);
      coord_tag[i][j]=d;
      if (coord_tag[i][j] == 111111.0)
	return 0;
    }
  }

  return 1;
}

void test_velo(double coord_tag[MAXNUMATOM][3],double coord_tag_old[MAXNUMATOM][3],double velo_tag[MAXNUMATOM][3],int numatom, double deltat)
{
  int i,j;

  for (i=0;i<numatom;++i)
  {
    for (j=0;j<3;++j)
    {
      velo_tag[i][j] = (coord_tag[i][j]-coord_tag_old[i][j])/deltat;
    }
  }
}
