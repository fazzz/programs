
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "rmsd.h"
#include "EF.h"
#include "PT.h"
#include "f2c.h"
#include "clapack.h"

double Coff0=1.0,Coff1=1.0,Coff2=1.0,rmsd2=0.0;

double rmsd_qcp(double *coordA,double *coordB,int numatom, int MODE) {
  int i,j,k,l;
  int numatomp;
  double mat[3][3];
  double K[4][4];
  double det,lambda=0.0,GA=0.0,GB=0.0,rmsd=0.0,lambda_ini;
  double *crdA,*crdB;
  double *mass;

  crdA=(double *)gcemalloc(sizeof(double)*numatom*3);
  crdB=(double *)gcemalloc(sizeof(double)*numatom*3);
  mass=(double *)gcemalloc(sizeof(double)*numatom);


  l=0;
  for (j=0;j<numatom;++j) {
    if (MODE==CA) {
      if (strncmp(AP.IGRAPH[j],"CA",2)==0) {
	mass[l]=AP.AMASS[j];
	for (k=0;k<3;++k) {
	  crdA[l*3+k]=coordA[j*3+k];
	  crdB[l*3+k]=coordB[j*3+k];
	}
	++l;
      }
    }
    else if (MODE==HV) {
      if (strncmp(AP.IGRAPH[j],"H",1)!=0) {
	mass[l]=AP.AMASS[j];
	for (k=0;k<3;++k) {
	  crdA[l*3+k]=coordA[j*3+k];
	  crdB[l*3+k]=coordB[j*3+k];
	}
	++l;
      }
    }
    else {
      mass[l]=AP.AMASS[j];
      for (k=0;k<3;++k) {
	crdA[l*3+k]=coordA[j*3+k];
	crdB[l*3+k]=coordB[j*3+k];
      }
      ++l;
    }
  }
  numatomp=l;

  transCentroid(crdA,crdB,numatomp,mass);

  rmsd=0.0;
  for (i=0;i<numatomp;++i)
    for (j=0;j<3;++j)
      rmsd+=sqrt((crdA[i*3+j]-crdB[i*3+j])*(crdA[i*3+j]-crdB[i*3+j]));
  rmsd=rmsd/numatomp;

  GA=CalcG(crdA,numatomp);
  GB=CalcG(crdB,numatomp);
  //  det=0.0001;
  det=1.0e-7;

  mmult(crdA,crdB,mat,numatomp,mass);
  fomKmat(K,mat);
 
  lambda_ini=(GA+GB)/2.0;
  fomPolynominal(mat,K);  
  lambda=Newton_Rapson(lambda_ini,det);
  
  rmsd=sqrt((GA+GB-2.0*lambda)/numatomp);

  return rmsd;
}


void mmult(double *coordA, double *coordB, double matrix[3][3], int numatom, double *mass) {
  int i,j,k;
  
  for (i=0;i<3;++i) for (j=0;j<3;++j)   matrix[i][j] = 0.0;

  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      for (k=0;k<numatom;++k)
	matrix[i][j] += coordA[k*3+i]*coordB[k*3+j];
	  /***//*AP.AMASS*//*mass[k];*/
}

void fomKmat(double Kmat[4][4], double mat[3][3]){
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

void fomPolynominal(/*double C0, double C1, double C2,*/double mat[3][3],double Kmat[4][4]) {
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
}

double Newton_Rapson(/*double C0, double C1, double C2, */double lambda_ini,  double det) {
  int i,j=0;
  double lambda=0.0/*=1.0*/, lambda_old=0.0;
  i=abs(lambda-lambda_old);
  lambda=lambda_ini;
  do 
  {
    lambda_old=lambda;
    lambda=lambda-(pow(lambda,4)+Coff2*pow(lambda,2)+Coff1*lambda+Coff0)/(4.0*pow(lambda,3)+2.0*Coff2*lambda+Coff1);
    ++j;
  } while (abs(lambda-lambda_old)<det);
  
return lambda;
}

double CalcG(double *coordA,int numatom) {
  int i,j;
  double G=0.0;

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) G +=  coordA[i*3+j]*coordA[i*3+j];
  
  return G;

}

void transCentroid(double *coordA, double *coordB, int numatom, double *mass) {
  int i,j;
  double cm_A[3],cm_B[3];
  double summass;

  summass=0.0;

  for (j=0;j<3;++j)  {
    cm_A[j] = 0.0;
    cm_B[j] = 0.0; 
  }
  
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      cm_A[j] += coordA[i*3+j]*/*AP.AMASS*/mass[i];
      cm_B[j] += coordB[i*3+j]*/*AP.AMASS*/mass[i]; 
    }
    summass += /*AP.AMASS*/mass[i];
  }
  
 for (i=0;i<3;++i)  {
   cm_A[i] = cm_A[i]/(summass);
   cm_B[i] = cm_B[i]/(summass); 
 }
 
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j)  {
      coordA[i*3+j] -= cm_A[j]; 
      coordB[i*3+j] -= cm_B[j]; 
    }
  }
  
}


double InnerProduct(double *A, double **coords1, double **coords2, const int len, const double *weight)
{
    double          x1, x2, y1, y2, z1, z2;
    int             i;
    const double   *fx1 = coords1[0], *fy1 = coords1[1], *fz1 = coords1[2];
    const double   *fx2 = coords2[0], *fy2 = coords2[1], *fz2 = coords2[2];
    double          G1 = 0.0, G2 = 0.0;

    A[0] = A[1] = A[2] = A[3] = A[4] = A[5] = A[6] = A[7] = A[8] = 0.0;

    if (weight != NULL)
    {
        for (i = 0; i < len; ++i)
        {
            x1 = weight[i] * fx1[i];
            y1 = weight[i] * fy1[i];
            z1 = weight[i] * fz1[i];

            G1 += x1 * fx1[i] + y1 * fy1[i] + z1 * fz1[i];

            x2 = fx2[i];
            y2 = fy2[i];
            z2 = fz2[i];

            G2 += weight[i] * (x2 * x2 + y2 * y2 + z2 * z2);

            A[0] +=  (x1 * x2);
            A[1] +=  (x1 * y2);
            A[2] +=  (x1 * z2);

            A[3] +=  (y1 * x2);
            A[4] +=  (y1 * y2);
            A[5] +=  (y1 * z2);

            A[6] +=  (z1 * x2);
            A[7] +=  (z1 * y2);
            A[8] +=  (z1 * z2);   
        }
    }
    else
    {
        for (i = 0; i < len; ++i)
        {
            x1 = fx1[i];
            y1 = fy1[i];
            z1 = fz1[i];

            G1 += x1 * x1 + y1 * y1 + z1 * z1;

            x2 = fx2[i];
            y2 = fy2[i];
            z2 = fz2[i];

            G2 += (x2 * x2 + y2 * y2 + z2 * z2);

            A[0] +=  (x1 * x2);
            A[1] +=  (x1 * y2);
            A[2] +=  (x1 * z2);

            A[3] +=  (y1 * x2);
            A[4] +=  (y1 * y2);
            A[5] +=  (y1 * z2);

            A[6] +=  (z1 * x2);
            A[7] +=  (z1 * y2);
            A[8] +=  (z1 * z2);  
        }
    }

    return (G1 + G2) * 0.5;
}


int FastCalcRMSDAndRotation(double *rot, double *A, double *rmsd, double E0, int len, double minScore)
{
    double Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz;
    double Szz2, Syy2, Sxx2, Sxy2, Syz2, Sxz2, Syx2, Szy2, Szx2,
           SyzSzymSyySzz2, Sxx2Syy2Szz2Syz2Szy2, Sxy2Sxz2Syx2Szx2,
           SxzpSzx, SyzpSzy, SxypSyx, SyzmSzy,
           SxzmSzx, SxymSyx, SxxpSyy, SxxmSyy;
    double C[4];
    int i;
    double mxEigenV; 
    double oldg = 0.0;
    double b, a, delta, rms, qsqr;
    double q1, q2, q3, q4, normq;
    double a11, a12, a13, a14, a21, a22, a23, a24;
    double a31, a32, a33, a34, a41, a42, a43, a44;
    double a2, x2, y2, z2; 
    double xy, az, zx, ay, yz, ax; 
    double a3344_4334, a3244_4234, a3243_4233, a3143_4133,a3144_4134, a3142_4132; 
    double evecprec = 1e-6;
    double evalprec = 1e-9;

    Sxx = A[0]; Sxy = A[1]; Sxz = A[2];
    Syx = A[3]; Syy = A[4]; Syz = A[5];
    Szx = A[6]; Szy = A[7]; Szz = A[8];

    Sxx2 = Sxx * Sxx;
    Syy2 = Syy * Syy;
    Szz2 = Szz * Szz;

    Sxy2 = Sxy * Sxy;
    Syz2 = Syz * Syz;
    Sxz2 = Sxz * Sxz;

    Syx2 = Syx * Syx;
    Szy2 = Szy * Szy;
    Szx2 = Szx * Szx;

    SyzSzymSyySzz2 = 2.0*(Syz*Szy - Syy*Szz);
    Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;

    C[2] = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
    C[1] = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz);

    SxzpSzx = Sxz + Szx;
    SyzpSzy = Syz + Szy;
    SxypSyx = Sxy + Syx;
    SyzmSzy = Syz - Szy;
    SxzmSzx = Sxz - Szx;
    SxymSyx = Sxy - Syx;
    SxxpSyy = Sxx + Syy;
    SxxmSyy = Sxx - Syy;
    Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

    C[0] = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
         + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
         + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz))
         + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz))
         + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz))
         + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz));

    mxEigenV = E0;
    for (i = 0; i < 50; ++i)
    {
        oldg = mxEigenV;
        x2 = mxEigenV*mxEigenV;
        b = (x2 + C[2])*mxEigenV;
        a = b + C[1];
        delta = ((a*mxEigenV + C[0])/(2.0*x2*mxEigenV + b + a));
        mxEigenV -= delta;
        if (fabs(mxEigenV - oldg) < fabs((evalprec)*mxEigenV))
            break;
    }

    if (i == 50) 
       fprintf(stderr,"\nMore than %d iterations needed!\n", i);

    rms = sqrt(2.0 * (E0 - mxEigenV)/len);
    (*rmsd) = rms;

    if (minScore > 0) 
        if (rms < minScore)
            return (-1); // Don't bother with rotation. 

    a11 = SxxpSyy + Szz-mxEigenV; a12 = SyzmSzy; a13 = - SxzmSzx; a14 = SxymSyx;
    a21 = SyzmSzy; a22 = SxxmSyy - Szz-mxEigenV; a23 = SxypSyx; a24= SxzpSzx;
    a31 = a13; a32 = a23; a33 = Syy-Sxx-Szz - mxEigenV; a34 = SyzpSzy;
    a41 = a14; a42 = a24; a43 = a34; a44 = Szz - SxxpSyy - mxEigenV;
    a3344_4334 = a33 * a44 - a43 * a34; a3244_4234 = a32 * a44-a42*a34;
    a3243_4233 = a32 * a43 - a42 * a33; a3143_4133 = a31 * a43-a41*a33;
    a3144_4134 = a31 * a44 - a41 * a34; a3142_4132 = a31 * a42-a41*a32;
    q1 =  a22*a3344_4334-a23*a3244_4234+a24*a3243_4233;
    q2 = -a21*a3344_4334+a23*a3144_4134-a24*a3143_4133;
    q3 =  a21*a3244_4234-a22*a3144_4134+a24*a3142_4132;
    q4 = -a21*a3243_4233+a22*a3143_4133-a23*a3142_4132;

    qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

/* The following code tries to calculate another column in the adjoint matrix when the norm of the 
   current column is too small.
   Usually this commented block will never be activated.  To be absolutely safe this should be
   uncommented, but it is most likely unnecessary.  
*/
    if (qsqr < evecprec)
    {
        q1 =  a12*a3344_4334 - a13*a3244_4234 + a14*a3243_4233;
        q2 = -a11*a3344_4334 + a13*a3144_4134 - a14*a3143_4133;
        q3 =  a11*a3244_4234 - a12*a3144_4134 + a14*a3142_4132;
        q4 = -a11*a3243_4233 + a12*a3143_4133 - a13*a3142_4132;
        qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4;

        if (qsqr < evecprec)
        {
            double a1324_1423 = a13 * a24 - a14 * a23, a1224_1422 = a12 * a24 - a14 * a22;
            double a1223_1322 = a12 * a23 - a13 * a22, a1124_1421 = a11 * a24 - a14 * a21;
            double a1123_1321 = a11 * a23 - a13 * a21, a1122_1221 = a11 * a22 - a12 * a21;

            q1 =  a42 * a1324_1423 - a43 * a1224_1422 + a44 * a1223_1322;
            q2 = -a41 * a1324_1423 + a43 * a1124_1421 - a44 * a1123_1321;
            q3 =  a41 * a1224_1422 - a42 * a1124_1421 + a44 * a1122_1221;
            q4 = -a41 * a1223_1322 + a42 * a1123_1321 - a43 * a1122_1221;
            qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4;

            if (qsqr < evecprec)
            {
                q1 =  a32 * a1324_1423 - a33 * a1224_1422 + a34 * a1223_1322;
                q2 = -a31 * a1324_1423 + a33 * a1124_1421 - a34 * a1123_1321;
                q3 =  a31 * a1224_1422 - a32 * a1124_1421 + a34 * a1122_1221;
                q4 = -a31 * a1223_1322 + a32 * a1123_1321 - a33 * a1122_1221;
                qsqr = q1*q1 + q2 *q2 + q3*q3 + q4*q4;
                
                if (qsqr < evecprec)
                {
                    /* if qsqr is still too small, return the identity matrix. */
                    rot[0] = rot[4] = rot[8] = 1.0;
                    rot[1] = rot[2] = rot[3] = rot[5] = rot[6] = rot[7] = 0.0;

                    return(0);
                }
            }
        }
    }

    normq = sqrt(qsqr);
    q1 /= normq;
    q2 /= normq;
    q3 /= normq;
    q4 /= normq;

    a2 = q1 * q1;
    x2 = q2 * q2;
    y2 = q3 * q3;
    z2 = q4 * q4;

    xy = q2 * q3;
    az = q1 * q4;
    zx = q4 * q2;
    ay = q1 * q3;
    yz = q3 * q4;
    ax = q1 * q2;

    rot[0] = a2 + x2 - y2 - z2;
    rot[1] = 2 * (xy + az);
    rot[2] = 2 * (zx - ay);
    rot[3] = 2 * (xy - az);
    rot[4] = a2 - x2 + y2 - z2;
    rot[5] = 2 * (yz + ax);
    rot[6] = 2 * (zx + ay);
    rot[7] = 2 * (yz - ax);
    rot[8] = a2 - x2 - y2 + z2;

    return (1);
}


void CenterCoords(double **coords, const int len, const double *weight)
{
    int             i;
    double          xsum, ysum, zsum, wsum;
    double         *x = coords[0], *y = coords[1], *z = coords[2];

    xsum = ysum = zsum = 0.0;

    if (weight != NULL)
    {
        wsum = 0.0;
        for (i = 0; i < len; ++i)
        {
            xsum += weight[i] * x[i];
            ysum += weight[i] * y[i];
            zsum += weight[i] * z[i];
            
            wsum += weight[i];
        }

        xsum /= wsum;
        ysum /= wsum;
        zsum /= wsum;
    }
    else
    {
        for (i = 0; i < len; ++i)
        {
            xsum += x[i];
            ysum += y[i];
            zsum += z[i];
        }

        xsum /= len;
        ysum /= len;
        zsum /= len;
    }

    for (i = 0; i < len; ++i)
    {
        x[i] -= xsum;
        y[i] -= ysum;
        z[i] -= zsum;
    }
}


double CalcRMSDRotationalMatrix(double **coords1, double **coords2, const int len, double *rot, const double *weight)
{
    double A[9], rmsd;
    /* center the structures */
    CenterCoords(coords1, len, weight);
    CenterCoords(coords2, len, weight);

    /* calculate the (weighted) inner product of two structures */
    double E0 = InnerProduct(A, coords1, coords2, len, weight);

    /* calculate the RMSD & rotational matrix */
    FastCalcRMSDAndRotation(rot, A, &rmsd, E0, len, -1);

    return rmsd;
}

