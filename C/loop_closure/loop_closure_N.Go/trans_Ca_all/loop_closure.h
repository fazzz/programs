#include <stdio.h>
#include <string.h>
#include <math.h>

#define pi 3.14159265
#define NUM 20

double alpha=9.0*pi/180.0;
double beta=70.5*pi/180.0;
double rou=3.519,eta=1.436;
double p1[3];
double p2[3];
double p3[3];

int loop_closure(double s[3], double u[3], double v[3], double psi1[4], double psi2[4], double psi3[4], double phi1[4], double phi2[4], double phi3[4]);
double calc_gfunc(double u[3], double cospsi1, double sinpsi1, double cospsi2, double sinpsi2, double cosphi1, double sinphi1, double cosphi2, double sinphi2, double alpha, double beta);
void calc_r(double alpha, double beta, double psi1, double s[3], double rou, double eta, double r[3]);
void mmult(double m1[3][3], double m2[3][3], double m1m2[3][3]);
void mvult(double m1[3][3], double v2[3], double m1v2[3]);
void fomTmat(double cos, double sin, double mat[3][3]);
void fomRmat(double cos, double sin, double mat[3][3]);
int calc_phi1(double r[3], double w, double rou, double eta, double cosphi1[2], double sinphi1[2],double qz[2][3]);
int calc_phi2(double x,double w,double rou,double eta,double alpha,double beta,double cosphi2[2],double sinphi2[2]);
void calc_psi2(double w,double rou,double eta,double alpha,double beta,double qz[2][3],double cosphi2[2],double sinphi2[2], double cospsi2[2][2], double sinpsi2[2][2]);
double calc_psi3(double psi1,double psi2,double phi1,double phi2,double u[3]);
double calc_phi3(double psi1,double psi2,double psi3,double phi1,double phi2,double v[3]);
void fomITmat(double cos, double sin, double mat[3][3]);
void fomIRmat(double cos, double sin, double mat[3][3]);

int loop_closure(double s[3], double u[3],  double v[3], double psi1[4], double psi2[4], double psi3[4], double phi1[4], double phi2[4], double phi3[4])
{
	char outname[100];
	int i,j,k,l,m,n,test1,test2,len,numans;
	double x,y,z;
	double f[10000][4];
	double r[3],w;
	double cospsi2[2][2],sinpsi2[2][2],cosphi1[2],sinphi1[2],cosphi2[2],sinphi2[2];
	double qz[2][3];
	double test3[10000];
	double psi1dummy,psi3dummy,phi3dummy;
    double psi1ans[10],psi2ans[10],psi3ans[10],phi1ans[10],phi2ans[10],phi3ans[10];
	double pep1[4][3],pep2[4][3];

	double g[10000][4];
    FILE *in;

	for (i=0;i<360*NUM;++i)
	{
		test3[i] = 0;
		for (j=0;j<4;++j)
			g[i][j] = 0.0;
	}

	m = 0;
	n = 0;
	for (i=0;i<360*NUM;++i)
	{
		psi1dummy = (double)i*pi/(180.0*(double)NUM);

		// equ 21
		calc_r(alpha, beta, psi1dummy, s, rou, eta, r);
		// equ 30
		w = (r[0]*r[0]+r[1]*r[1]+r[2]*r[2]-2.0*rou*r[0])/(2.0*eta);

		// equ 29,31
		test1 = calc_phi1(r,w,rou,eta,cosphi1,sinphi1,qz);
		// equ 35
		test2 = calc_phi2(r[0],w,rou,eta,alpha,beta,cosphi2,sinphi2);
		if (test1 == 1 && test2 == 1)
		{
			test3[i] = 1;
			// equ 36
			calc_psi2(w,rou,eta,alpha,beta,qz,cosphi2,sinphi2,cospsi2,sinpsi2);

			l=0;
			for (j=0;j<2;++j)
			{
				for (k=0;k<2;++k)
				{
					g[i][l] = calc_gfunc(u,cos(psi1dummy),sin(psi1dummy),cospsi2[j][k],sinpsi2[j][k],cosphi1[j],sinphi1[j],cosphi2[k],sinphi2[k],alpha,beta);
					if (-0.0003 <= g[i][l] && 0.0003 >= g[i][l])
					{
						psi1ans[m] = psi1dummy;
						psi2ans[m] = acos(cospsi2[j][k]);
						if (sinpsi2[j][k] < 0.0)
						{
							psi2ans[m] = -psi2ans[m]+2.0*pi;
						}
						phi1ans[m] = acos(cosphi1[j]);
						if (sinphi1[j] < 0.0)
						{
							phi1ans[m] = -phi1ans[m]+2.0*pi;
						}
						phi2ans[m] = acos(cosphi2[k]);
						if (sinphi2[k] < 0.0)
						{
							phi2ans[m] = -phi2ans[m]+2.0*pi;
						}
                        if ((psi3ans[m]=calc_psi3(psi1ans[m],psi2ans[m],phi1ans[m],phi2ans[m],u))<0.0)
                            psi3ans[m]+=2.0*pi;
                        if ((phi3ans[m]=calc_phi3(psi1ans[m],psi2ans[m],psi3ans[m],phi1ans[m],phi2ans[m],v))<0.0)
                            phi3ans[m]+=2.0*pi;
						++m;
					}
					++l;
				}
			}
		}
	}
	numans = m;


	printf("\n");
	printf("num  psi1  psi2  psi3  phi1  phi2  phi3     \n");
	for (i=0;i<numans;++i)
	{
        psi1[i]=psi1ans[i]*180.0/pi;
        psi2[i]=psi2ans[i]*180.0/pi;
        psi3[i]=psi3ans[i]*180.0/pi;
        phi1[i]=phi1ans[i]*180.0/pi;
        phi2[i]=phi2ans[i]*180.0/pi;
        phi3[i]=phi3ans[i]*180.0/pi;
		printf("%-d: %5.1lf %5.1lf %5.1lf %5.1lf %5.1lf %5.1lf \n",i,psi1[i],psi2[i],psi3[i],phi1[i]-180.0,phi2[i]-180.0,phi3[i]-180.0);
	}

	return numans;
}

double calc_gfunc(double u[3], double cospsi1, double sinpsi1, double cospsi2, double sinpsi2, double cosphi1, double sinphi1, double cosphi2, double sinphi2, double alpha, double beta)
{
	double g;
	double Talpha[3][3],Tbeta[3][3],Rpsi1[3][3],Rphi1[3][3],Rpsi2[3][3],Rphi2[3][3];
	double TaRs1[3][3],TaRs1Tb[3][3],TaRs1TbRh1[3][3],TaRs1TbRh1Ta[3][3],TaRs1TbRh1TaRs2[3][3],TaRs1TbRh1TaRs2Tb[3][3],TaRs1TbRh1TaRs2TbRh2[3][3],TaRs1TbRh1TaRs2TbRh2Ta[3][3];
	double cosalpha,sinalpha,cosbeta,sinbeta;

	cosalpha = cos(alpha);
	sinalpha = sin(alpha);
	cosbeta = cos(beta);
	sinbeta = sin(beta);

	fomTmat(cosalpha, sinalpha, Talpha);
	fomTmat(cosbeta, sinbeta, Tbeta);
	fomRmat(cospsi1, sinpsi1, Rpsi1);
	fomRmat(cospsi2, sinpsi2, Rpsi2);
	fomRmat(cosphi1, sinphi1, Rphi1);
	fomRmat(cosphi2, sinphi2, Rphi2);

	mmult(Talpha,Rpsi1,TaRs1);
	mmult(TaRs1,Tbeta,TaRs1Tb);
	mmult(TaRs1Tb,Rphi1,TaRs1TbRh1);
	mmult(TaRs1TbRh1,Talpha,TaRs1TbRh1Ta);
	mmult(TaRs1TbRh1Ta,Rpsi2,TaRs1TbRh1TaRs2);
	mmult(TaRs1TbRh1TaRs2,Tbeta,TaRs1TbRh1TaRs2Tb);
	mmult(TaRs1TbRh1TaRs2Tb,Rphi2,TaRs1TbRh1TaRs2TbRh2);
	mmult(TaRs1TbRh1TaRs2TbRh2,Talpha,TaRs1TbRh1TaRs2TbRh2Ta);

	g = u[0]*TaRs1TbRh1TaRs2TbRh2Ta[0][0]+u[1]*TaRs1TbRh1TaRs2TbRh2Ta[1][0]+u[2]*TaRs1TbRh1TaRs2TbRh2Ta[2][0]-cos(beta);

	return g;
}

void calc_r(double alpha, double beta, double psi1, double s[3], double rou, double eta, double r[3])
{
	int i,j;
	double TalphaInv[3][3],TbetaInv[3][3],Rpsi1Inv[3][3],TbR1[3][3],TbR1Ta[3][3];
	double r_dummy[3];

	for (i=0;i<3;++i)
	{
		r[i] = 0.0;
		for (j=0;j<3;++j)
		{
			TalphaInv[i][j] = 0.0;
			TbetaInv[i][j] = 0.0;
			Rpsi1Inv[i][j] = 0.0;
		}
	}

	TalphaInv[0][0] = cos(alpha);
	TalphaInv[0][1] = sin(alpha);
	TalphaInv[1][0] = -sin(alpha);
	TalphaInv[1][1] = cos(alpha);
	TalphaInv[2][2] = 1.0;

	TbetaInv[0][0] = cos(beta);
	TbetaInv[0][1] = sin(beta);
	TbetaInv[1][0] = -sin(beta);
	TbetaInv[1][1] = cos(beta);
	TbetaInv[2][2] = 1.0;

	Rpsi1Inv[0][0] = 1.0;
	Rpsi1Inv[1][1] = cos(psi1);
	Rpsi1Inv[1][2] = sin(psi1);
	Rpsi1Inv[2][1] = -sin(psi1);
	Rpsi1Inv[2][2] = cos(psi1);

	mmult(TbetaInv,Rpsi1Inv,TbR1);
	mmult(TbR1,TalphaInv,TbR1Ta);

	r_dummy[0] = s[0]-rou;
	r_dummy[1] = s[1]-eta;
	r_dummy[2] = s[2];

	for (i=0;i<3;++i)
	{
		for (j=0;j<3;++j)
		{
			r[i] += TbR1Ta[i][j]*r_dummy[j];
		}
	}
}

void mmult(double m1[3][3], double m2[3][3], double m1m2[3][3])
{
	int i,j,k;

	for (i=0;i<3;++i)
		for (j=0;j<3;++j)
			m1m2[i][j] = 0.0;

	for (i=0;i<3;++i)
	{
		for (j=0;j<3;++j)
		{
			for (k=0;k<3;++k)
			{
				m1m2[i][j] += m1[i][k]*m2[k][j];
			}
		}
	}
}

void mvult(double m1[3][3], double v2[3], double m1v2[3])
{
    int i,j;

    for (i=0;i<3;++i)
    {
        m1v2[i]=0.0;
    }
    
    for (i=0;i<3;++i)
    {
        for (j=0;j<3;++j)
        {
            m1v2[i]+=m1[i][j]*v2[j];
        }
    }
}
    
    
void fomTmat(double cos, double sin, double mat[3][3])
{
	int i,j;

	for (i=0;i<3;++i)
	{
		for (j=0;j<3;++j)
		{
			mat[i][j] = 0.0;
		}
	}

	
	mat[0][0] = cos;
	mat[0][1] = -sin;
	mat[1][0] = sin;
	mat[1][1] = cos;
	mat[2][2] = 1.0;

}

void fomRmat(double cos, double sin, double mat[3][3])
{
	int i,j;

	for (i=0;i<3;++i)
	{
		for (j=0;j<3;++j)
		{
			mat[i][j] = 0.0;
		}
	}
	mat[0][0] = 1.0;
	mat[1][1] = cos;
	mat[1][2] = -sin;
	mat[2][1] = sin;
	mat[2][2] = cos;

}

void fomITmat(double cos, double sin, double mat[3][3])
{
	int i,j;

	for (i=0;i<3;++i)
	{
		for (j=0;j<3;++j)
		{
			mat[i][j] = 0.0;
		}
	}

	
	mat[0][0] = cos;
	mat[0][1] = sin;
	mat[1][0] = -sin;
	mat[1][1] = cos;
	mat[2][2] = 1.0;

}

void fomIRmat(double cos, double sin, double mat[3][3])
{
	int i,j;

	for (i=0;i<3;++i)
	{
		for (j=0;j<3;++j)
		{
			mat[i][j] = 0.0;
		}
	}
	mat[0][0] = 1.0;
	mat[1][1] = cos;
	mat[1][2] = sin;
	mat[2][1] = -sin;
	mat[2][2] = cos;

}

int calc_phi1(double r[3], double w, double rou, double eta, double cosphi1[2], double sinphi1[2],double qz[2][3])
{
	double D;

	D = r[1]*r[1]+r[2]*r[2]-w*w;
	if (D <= 0.0000 && D >= -0.00001)
	{
		D = 0.0;
	}

	if (D >= 0.0)
	{
		cosphi1[0] = (r[1]*w+r[2]*sqrt(D))/(r[1]*r[1]+r[2]*r[2]);
		sinphi1[0] = (r[2]*w-r[1]*sqrt(D))/(r[1]*r[1]+r[2]*r[2]);
		cosphi1[1] = (r[1]*w-r[2]*sqrt(D))/(r[1]*r[1]+r[2]*r[2]);
		sinphi1[1] = (r[2]*w+r[1]*sqrt(D))/(r[1]*r[1]+r[2]*r[2]);
		qz[0][0] = r[0]-rou;
		qz[0][1] = w-eta;
		qz[0][2] = sqrt(D);
		qz[1][0] = r[0]-rou;
		qz[1][1] = w-eta;
		qz[1][2] = -sqrt(D);
	}
	else
	{
		return 0;
	}

	return 1;
}

int calc_phi2(double x,double w,double rou,double eta,double alpha,double beta,double cosphi2[2],double sinphi2[2])
{
	cosphi2[0] = (rou*cos(beta)-(x-rou)*cos(alpha)-(w-eta)*sin(alpha))/(eta*sin(beta));
	cosphi2[1] = cosphi2[0];

	if (-1.00001 <= cosphi2[0] && cosphi2[0] <= 1.00001)
	{
		if(-1.00001 <= cosphi2[0] && cosphi2[0] <= -1.0)
		{
			cosphi2[0] = -1.0;
			cosphi2[1] = cosphi2[0];
		}
		else if(1.0 <= cosphi2[0] && cosphi2[0] <= 1.00001)
		{
			cosphi2[0] = 1.0;
			cosphi2[1] = cosphi2[0];
		}

		sinphi2[0] =  sqrt(1.0-cosphi2[0]*cosphi2[0]);
		sinphi2[1] = -sqrt(1.0-cosphi2[0]*cosphi2[0]);
	}
	else
	{
		return 0;
	}

	return 1;
}

void calc_psi2(double w,double rou,double eta,double alpha,double beta,double qz[2][3],double cosphi2[2],double sinphi2[2], double cospsi2[2][2], double sinpsi2[2][2])
{
	int i,j;
	double fact;

	for (i=0;i<2;++i)
	{
		for (j=0;j<2;++j)
		{
			fact = (-qz[i][0]*sin(alpha)+qz[i][1]*cos(alpha))*(-qz[i][0]*sin(alpha)+qz[i][1]*cos(alpha))+qz[i][2]*qz[i][2];
			cospsi2[i][j] = (
				             (rou*sin(beta)+eta*cosphi2[j]*cos(beta))*(-qz[i][0]*sin(alpha)+qz[i][1]*cos(alpha))
				            +eta*qz[i][2]*sinphi2[j]
				            )
				            /fact;

			sinpsi2[i][j] = (
				             (rou*sin(beta)+eta*cosphi2[j]*cos(beta))*qz[i][2]
				        -eta*(-qz[i][0]*sin(alpha)+qz[i][1]*cos(alpha))*sinphi2[j]
				            )
							/fact;

			if (cospsi2[i][j] > 1.0)
			{
				cospsi2[i][j] = 1.0;
			}
			else if (cospsi2[i][j] < -1.0)
			{
				cospsi2[i][j] = -1.0;
			}
			if (sinpsi2[i][j] > 1.0)
			{
				sinpsi2[i][j] = 1.0;
			}
			else if (sinpsi2[i][j] < -1.0)
			{
				sinpsi2[i][j] = -1.0;
			}
		}
	}
}


double calc_psi3(double psi1,double psi2,double phi1,double phi2,double u[3])
{
    double psi3;
	double ITalpha[3][3],ITbeta[3][3],IRpsi1[3][3],IRphi1[3][3],IRpsi2[3][3],IRphi2[3][3];
	double cosalpha,sinalpha,cosbeta,sinbeta,cospsi3,sinpsi3;
    double vect[3];
    double ITaRh2[3][3],ITaRh2Tb[3][3],ITaRh2TbRs2[3][3],ITaRh2TbRs2Ta[3][3],ITaRh2TbRs2TaRh1[3][3],ITaRh2TbRs2TaRh1Tb[3][3],ITaRh2TbRs2TaRh1TbRs1[3][3],ITaRh2TbRs2TaRh1TbRs1Ta[3][3];

	fomITmat(cos(alpha), sin(alpha), ITalpha);
	fomITmat(cos(beta), sin(beta), ITbeta);
	fomIRmat(cos(psi1), sin(psi1), IRpsi1);
	fomIRmat(cos(psi2), sin(psi2), IRpsi2);
	fomIRmat(cos(phi1), sin(phi1), IRphi1);
	fomIRmat(cos(phi2), sin(phi2), IRphi2);

	mmult(ITalpha,IRphi2,ITaRh2);
	mmult(ITaRh2,ITbeta,ITaRh2Tb);
	mmult(ITaRh2Tb,IRpsi2,ITaRh2TbRs2);
	mmult(ITaRh2TbRs2,ITalpha,ITaRh2TbRs2Ta);
	mmult(ITaRh2TbRs2Ta,IRphi1,ITaRh2TbRs2TaRh1);
	mmult(ITaRh2TbRs2TaRh1,ITbeta,ITaRh2TbRs2TaRh1Tb);
	mmult(ITaRh2TbRs2TaRh1Tb,IRpsi1,ITaRh2TbRs2TaRh1TbRs1);
	mmult(ITaRh2TbRs2TaRh1TbRs1,ITalpha,ITaRh2TbRs2TaRh1TbRs1Ta);

    mvult(ITaRh2TbRs2TaRh1TbRs1Ta,u,vect);
    cospsi3=vect[1]/sin(beta);
    sinpsi3=vect[2]/sin(beta);

    if (sinpsi3 > 0.0)
        psi3=acos(cospsi3);
    else
        psi3=-acos(cospsi3);

    return psi3;
}


double calc_phi3(double psi1,double psi2, double psi3,double phi1,double phi2,double v[3])
{
    double phi3;
	double ITalpha[3][3],ITbeta[3][3],IRpsi1[3][3],IRphi1[3][3],IRpsi2[3][3],IRphi2[3][3],IRpsi3[3][3];
	double cosalpha,sinalpha,cosbeta,sinbeta,cosphi3,sinphi3;
    double vect[3];
    double ITbRs3[3][3],ITbRs3Ta[3][3],ITbRs3TaRh2[3][3],ITbRs3TaRh2Tb[3][3],ITbRs3TaRh2TbRs2[3][3],ITbRs3TaRh2TbRs2Ta[3][3],ITbRs3TaRh2TbRs2TaRh1[3][3],ITbRs3TaRh2TbRs2TaRh1Tb[3][3],ITbRs3TaRh2TbRs2TaRh1TbRs1[3][3],ITbRs3TaRh2TbRs2TaRh1TbRs1Ta[3][3];

	fomITmat(cos(alpha), sin(alpha), ITalpha);
	fomITmat(cos(beta), sin(beta), ITbeta);
	fomIRmat(cos(psi1), sin(psi1), IRpsi1);
	fomIRmat(cos(psi2), sin(psi2), IRpsi2);
	fomIRmat(cos(phi1), sin(phi1), IRphi1);
	fomIRmat(cos(phi2), sin(phi2), IRphi2);
	fomIRmat(cos(psi3), sin(psi3), IRpsi3);

	mmult(ITbeta,IRpsi3,ITbRs3);
	mmult(ITbRs3,ITalpha,ITbRs3Ta);
	mmult(ITbRs3Ta,IRphi2,ITbRs3TaRh2);
	mmult(ITbRs3TaRh2,ITbeta,ITbRs3TaRh2Tb);
	mmult(ITbRs3TaRh2Tb,IRpsi2,ITbRs3TaRh2TbRs2);
	mmult(ITbRs3TaRh2TbRs2,ITalpha,ITbRs3TaRh2TbRs2Ta);
	mmult(ITbRs3TaRh2TbRs2Ta,IRphi1,ITbRs3TaRh2TbRs2TaRh1);
	mmult(ITbRs3TaRh2TbRs2TaRh1,ITbeta,ITbRs3TaRh2TbRs2TaRh1Tb);
	mmult(ITbRs3TaRh2TbRs2TaRh1Tb,IRpsi1,ITbRs3TaRh2TbRs2TaRh1TbRs1);
	mmult(ITbRs3TaRh2TbRs2TaRh1TbRs1,ITalpha,ITbRs3TaRh2TbRs2TaRh1TbRs1Ta);

    mvult(ITbRs3TaRh2TbRs2TaRh1TbRs1Ta,v,vect);
    cosphi3=vect[1];
    sinphi3=vect[2];
    
    if (sinphi3 > 0.0)
        phi3=acos(cosphi3);
    else
        phi3=-acos(cosphi3);

    return phi3;
}


