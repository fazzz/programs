#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bestfit.h"
#include "quaternion.h"

#include "f2c.h"
#include "clapack.h"

double Coff0=1.0,Coff1=1.0,Coff2=1.0,rmsd2=0.0;

void mmult(double coordA[MAXNUMATOM][3], double coordB[MAXNUMATOM][3], double matrix[3][3], int numatom,double mass[MAXNUMATOM]);
void fomKmat(double Kmat[4][4], double mat[3][3]);
void fomPolynominal(/*double C0, double C1, double C2,*/ double mat[3][3], double Kmat[4][4]);
double Newton_Rapson(/*double C0, double C1, double C2,*/ double det);
double CalcG(double coordA[MAXNUMATOM][3],int numatom);
void transCentroid(double coordA[MAXNUMATOM][3], double coordB[MAXNUMATOM][3], double mass[MAXNUMATOM], int numatom);
void transMotion(double velo[MAXNUMATOM][3], double mass[MAXNUMATOM], int numatom);

double bestfit(double coord_ref[MAXNUMATOM][3]
	       ,double coord_tag[MAXNUMATOM][3],double velo_tag[MAXNUMATOM][3]
	       ,double mass[MAXNUMATOM],  int numatom)
{
  int i,j;
  double coord_bestfit[4],coord_tag_dummy[4],velo_bestfit[4],velo_tag_dummy[4];
  double mat[3][3];
  double K[4][4],Kdummy[16];
//  double Coff0=1.0,Coff1=1.0,Coff2=1.0;
  double det,lambda=0.0,GA=0.0,GB=0.0,rmsd=0.0;
  double q[4],eigenvalue[4],work[16];
  long int info=0,lwork=16,n=4;

    transCentroid(coord_ref,coord_tag,mass,numatom);

// for (i=0;i<numatom;++i)
//  printf("%12.8lf %12.8lf %12.8lf\n",coord_ref[i][0],coord_ref[i][1],coord_ref[i][2]);


    GA=CalcG(coord_ref,numatom);
    GB=CalcG(coord_tag,numatom);

    det=0.00001;

    mmult(coord_ref,coord_tag,mat,numatom,mass);
    fomKmat(K,mat);
 
    for (i=0;i<4;++i)
      eigenvalue[i]=0.0;

    for (i=0;i<16;++i)
      work[i]=0.0;

    Kdummy[0]=K[0][0];  Kdummy[4]=K[0][1];  Kdummy[8]=K[0][2];  Kdummy[12]=K[0][3];
    Kdummy[1]=K[1][0];  Kdummy[5]=K[1][1];  Kdummy[9]=K[1][2];  Kdummy[13]=K[1][3];
    Kdummy[2]=K[2][0];  Kdummy[6]=K[2][1];  Kdummy[10]=K[2][2]; Kdummy[14]=K[2][3];
    Kdummy[3]=K[3][0];  Kdummy[7]=K[3][1];  Kdummy[11]=K[3][2]; Kdummy[15]=K[3][3];

    dsyev_("V", "U", &n, Kdummy, &n, eigenvalue, work, &lwork, &info);
    for (i=0;i<4;++i)
    {
      q[i]=Kdummy[12+i];
    }
    for (i=0;i<numatom;++i)
    {
      coord_tag_dummy[0]=0.0;
      velo_tag_dummy[0]=0.0;
      for (j=1;j<4;++j)
      {
	coord_tag_dummy[j]=coord_tag[i][j-1];
	velo_tag_dummy[j]=velo_tag[i][j-1];
      }
      quaternion_rotation(coord_tag_dummy,q,coord_bestfit);
      quaternion_rotation(velo_tag_dummy,q,velo_bestfit);
      for (j=0;j<3;++j)
      {
	coord_tag[i][j]=coord_bestfit[j+1];
	velo_tag[i][j]=velo_bestfit[j+1];
      }
   }


    transMotion(velo_tag,mass,numatom);



  
//  printf("%lf %lf %lf\n",Coff0,Coff1,Coff2);
//    fomPolynominal(/*Coff0,Coff1,Coff2,*/mat,K);
//  printf("%lf %lf %lf\n",Coff0,Coff1,Coff2);

//    lambda=Newton_Rapson(/*Coff0,Coff1,Coff2,*/det);

    //    rmsd2=(abs(GA+GB-2.0*lambda))/numatom;
    //    rmsd2=sqrt(rmsd2);

    for (i=0;i<numatom;++i)
    {
      for (j=0;j<3;++j)
      {
	rmsd += (coord_tag[i][j]-coord_ref[i][j])*(coord_tag[i][j]-coord_ref[i][j]);
      }
    }
      
    rmsd = sqrt(abs(rmsd)/numatom);

    return rmsd;
}


void mmult(double coordA[MAXNUMATOM][3], double coordB[MAXNUMATOM][3], double matrix[3][3], int numatom,double mass[MAXNUMATOM])
{
	int i,j,k;

	for (i=0;i<3;++i)
		for (j=0;j<3;++j)
			matrix[i][j] = 0.0;

	for (i=0;i<3;++i)
	{
		for (j=0;j<3;++j)
		{
			for (k=0;k<numatom;++k)
			{
				matrix[i][j] += coordA[k][i]*coordB[k][j]*mass[k];
			}
		}
	}
}

void fomKmat(double Kmat[4][4], double mat[3][3])
{
	Kmat[0][0] =  mat[0][0]+mat[1][1]+mat[2][2];
	Kmat[1][1] =  mat[0][0]-mat[1][1]-mat[2][2];
	Kmat[2][2] = -mat[0][0]+mat[1][1]-mat[2][2];
	Kmat[3][3] = -mat[0][0]-mat[1][1]+mat[2][2];

	Kmat[0][1] = mat[1][2]-mat[2][1];
	Kmat[0][2] = mat[2][0]-mat[0][2];
	Kmat[0][3] = mat[0][1]-mat[1][0];

	Kmat[1][2] = mat[0][1]+mat[1][0];
	Kmat[2][3] = mat[1][2]+mat[2][1];
	Kmat[1][3] = mat[2][0]+mat[0][2];

	Kmat[1][0] = Kmat[0][1];
	Kmat[2][0] = Kmat[0][2];
	Kmat[2][1] = Kmat[1][2];

	Kmat[3][0] = Kmat[0][3];
	Kmat[3][1] = Kmat[1][3];
	Kmat[3][2] = Kmat[2][3];

}

void fomPolynominal(/*double C0, double C1, double C2,*/double mat[3][3],double Kmat[4][4])
{
  double D=(mat[0][1]*mat[0][1]+mat[0][2]*mat[0][2]-mat[1][0]*mat[1][0]-mat[2][0]*mat[2][0])
          *(mat[0][1]*mat[0][1]+mat[0][2]*mat[0][2]-mat[1][0]*mat[1][0]-mat[2][0]*mat[2][0]);
  double E,F,G,H,I;

  double Sxzp=mat[0][2]+mat[2][0];
  double Sxyp=mat[0][1]+mat[1][0];
  double Syzp=mat[1][2]+mat[2][1];
  double Sxzm=mat[0][2]-mat[2][0];
  double Sxym=mat[0][1]-mat[1][0];
  double Syzm=mat[1][2]-mat[2][1];

  double Syyzz=mat[1][1]*mat[2][2]-mat[1][2]*mat[2][1];

  double Spmm=mat[0][0]-mat[1][1]-mat[2][2];
  double Sppp=mat[0][0]+mat[1][1]+mat[2][2];
  double Spmp=mat[0][0]-mat[1][1]+mat[2][2];
  double Sppm=mat[0][0]+mat[1][1]-mat[2][2];

  double Coff0_2;

  E= (-mat[0][0]*mat[0][0]+mat[1][1]*mat[1][1]+mat[2][2]*mat[2][2]+mat[1][2]*mat[1][2]+mat[2][1]*mat[2][1]
      -2.0*Syyzz)
    *(-mat[0][0]*mat[0][0]+mat[1][1]*mat[1][1]+mat[2][2]*mat[2][2]+mat[1][2]*mat[1][2]+mat[2][1]*mat[2][1]
      +2.0*Syyzz);

  F= (-Sxzp*Syzm+Sxym*Spmm)*(-Sxzm*Syzp+Sxym*Spmp);

  G= (-Sxzp*Syzp-Sxyp*Sppm)*(-Sxzm*Syzm-Sxyp*Sppp);

  H= (Sxyp*Syzp+Sxzp*Spmp)*(-Sxzm*Syzm+Sxzp*Sppp);

  I= (Sxyp*Syzm+Sxzm*Spmm)*(-Sxym*Syzp+Sxzm*Sppm);
  
  Coff2=-2.0*(mat[0][0]*mat[0][0]+mat[0][1]*mat[0][1]+mat[0][2]*mat[0][2]
	     +mat[1][0]*mat[1][0]+mat[1][1]*mat[1][1]+mat[1][2]*mat[1][2]
	     +mat[2][0]*mat[2][0]+mat[2][1]*mat[2][1]+mat[2][2]*mat[2][2]);

  Coff1=8.0*(mat[0][0]*mat[1][2]*mat[2][1]+mat[1][1]*mat[2][0]*mat[0][2]+mat[2][2]*mat[0][1]*mat[1][0])
       -8.0*(mat[0][0]*mat[1][1]*mat[2][2]+mat[1][2]*mat[2][0]*mat[0][1]+mat[2][1]*mat[1][0]*mat[0][2]);
 
  Coff0=D+E+F+G+H+I;
  Coff0      =(Kmat[0][0]*Kmat[1][1]-Kmat[0][1]*Kmat[0][1])*(Kmat[2][2]*Kmat[3][3]-Kmat[2][3]*Kmat[2][3])
             +(Kmat[0][1]*Kmat[0][2]-Kmat[0][0]*Kmat[1][2])*(Kmat[1][2]*Kmat[3][3]-Kmat[2][3]*Kmat[0][3])
             +(Kmat[0][0]*Kmat[1][3]-Kmat[0][1]*Kmat[0][3])*(Kmat[1][2]*Kmat[2][3]-Kmat[2][2]*Kmat[1][3])
             +(Kmat[0][1]*Kmat[1][2]-Kmat[1][1]*Kmat[0][2])*(Kmat[0][2]*Kmat[3][3]-Kmat[2][3]*Kmat[1][3])
             +(Kmat[1][1]*Kmat[0][3]-Kmat[0][1]*Kmat[1][3])*(Kmat[0][2]*Kmat[2][3]-Kmat[2][2]*Kmat[0][3])
             +(Kmat[0][2]*Kmat[1][3]-Kmat[1][2]*Kmat[0][3])*(Kmat[0][2]*Kmat[1][3]-Kmat[1][2]*Kmat[0][3]);
  //  printf("%lf %lf %lf\n",Coff0,Coff1,Coff2);
}

double Newton_Rapson(/*double C0, double C1, double C2, */double det)
{
  int i;
  double lambda=1.0, lambda_old=0.0;
  //  printf("%d\n",abs(lambda-lambda_old));
  i=abs(lambda-lambda_old);
  while (abs(lambda-lambda_old)>det)
  {
    lambda_old=lambda;
    lambda=lambda-(pow(lambda,4)+Coff2*pow(lambda,2)+Coff1*lambda+Coff0)/(4.0*pow(lambda,3)+2.0*Coff2*lambda+Coff1);
  }
  //  printf("%lf\n",lambda);
  return lambda;
}

double CalcG(double coordA[MAXNUMATOM][3],int numatom)
{
  int i,j;
  double G=0.0;

  for (i=0;i<numatom;++i)
  {
    for (j=0;j<3;++j)
    {
      G +=  coordA[i][j]*coordA[i][j];
    }
  }

  return G;

}

void transCentroid(double coordA[MAXNUMATOM][3], double coordB[MAXNUMATOM][3], double mass[MAXNUMATOM], int numatom)
{
  int i,j;
  double cm_A[3],cm_B[3];
  double summass;

  summass=0.0;

    for (j=0;j<3;++j)
    {
      cm_A[j] = 0.0;
      cm_B[j] = 0.0; 
    }

  for (i=0;i<numatom;++i)
  {
    for (j=0;j<3;++j)
    {
      cm_A[j] += coordA[i][j]*mass[i];
      cm_B[j] += coordB[i][j]*mass[i]; 
    }
    summass += mass[i];
  }

 for (i=0;i<3;++i)
 {
   cm_A[i] = cm_A[i]/(summass);
   cm_B[i] = cm_B[i]/(summass); 
 }

  for (i=0;i<numatom;++i)
  {
    for (j=0;j<3;++j)
    {
      coordA[i][j] -= cm_A[j]; 
      coordB[i][j] -= cm_B[j]; 
    }
  }

}

void transMotion(double velo[MAXNUMATOM][3], double mass[MAXNUMATOM], int numatom)
{
  int i,j;
  double v_COM[3];
  double summass;

  summass=0.0;

    for (j=0;j<3;++j)
    {
      v_COM[j] = 0.0;
    }

  for (i=0;i<numatom;++i)
  {
    for (j=0;j<3;++j)
    {
      v_COM[j] += velo[i][j]*mass[i];
    }
    summass += mass[i];
  }

 for (i=0;i<3;++i)
 {
   v_COM[i] = v_COM[i]/(summass);
 }

  for (i=0;i<numatom;++i)
  {
    for (j=0;j<3;++j)
    {
      velo[i][j] -= v_COM[j]; 
    }
  }

}

