//////////////////////////////////////////////
//     Helical Configuration of Chain       //
//////////////////////////////////////////////

#include <stdio.h>

int main(int argc, char argv[])
{
  double r[3];
  double pi[3];
  double psi[3];

  double matA[3][3],matAnum[3][3][3];
  double vecB[3],vectBnum[3][3];

  double costheta,d,rou;
  double coord[100][3];

  struct protein prot;
  char *outputfilename;

  for (i=0;i<3;++i)
  {
    fomRmat(cos(pi[i]),sin(pi[i]),matAnum[i]);
    fomTmat(cos(psi[i]),sin(psi[i]),vectBnum[i]);
  }

  multMM(matAnum[2],matAnum[0],matAdummy);
  multMM(matAdummy,matAnum[1],matA);
  multMV(matAdummy,vectBnum[1],vectBdummy1);
  multMV(matAnum[2],vectBnum[0],vectBdummy2);

  for (i=0;i<3;++i)
  {
    vectB[i]=vectBdummy1[i]+vectBdummy2[i]+vectBnum[2][i]
  }

  costheta = (matA[0][0]+matA[1][1]+matA[2][2]-1)/2.0;
  d = sqrt(pow(vectB[0]*(matA[0][2]+matA[2][0])+vectB[1]*(matA[1][2]+matA[2][1])+vectB[2]*(matA[2][2]-matA[0][0]-matA[1][1]),2)/(3.0-matA[2][2]-matA[0][0]-matA[1][1])*(matA[2][2]-matA[0][0]-matA[1][1]+1));
  rou = sqrt((vectB[0]*vectB[0]+vectB[1]*vectB[1]+vectB[2]*vectB[2]-d*d)/(3.0-matA[2][2]-matA[1][1]-matA[0][0]));

  build(costheta,d,rou,coord);

  writepdb(outfilename,prot);

  return 0;
}

void build(double costheta, double d, double rou, double coord[12][3])
{
  int i;

  

}
