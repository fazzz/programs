#include <stdio.h>
#include <math.h>
#define pi 3.14159265

double calc_u(double u[3], double cospsi, double sinpsi, double cosphi, double sinphi,double cospsi2, double sinpsi2, double cosphi2, double sinphi2,double cospsi3, double sinpsi3, double cosalpha, double sinalpha, double cosbeta, double sinbeta);
double calc_s(double s[3], double cospsi, double sinpsi, double cosphi, double sinphi,double cospsi2, double sinpsi2, double cosphi2, double sinphi2, double cosalpha, double sinalpha, double cosbeta, double sinbeta, double q[3]);



void mmult(double m1[3][3], double m2[3][3], double m1m2[3][3]);
void mvult(double m1[3][3], double v2[3], double m1v2[3]);
void fomTmat(double cos, double sin, double mat[3][3]);
void fomRmat(double cos, double sin, double mat[3][3]);

int main(int argc, char argv[])
{
	double g;
	double psi,phi,psi2,phi2,psi3;
	double alpha,beta;
	double u[3],s[3],q[3];
	FILE *in, *out;

	alpha = 9.0*pi/180.0;
	beta = 70.5*pi/180.0;

	q[0] = 3.519;
	q[1] = 1.436;
	q[2] = 0.00;

/*	if (argc <= 1)
	{
		psi = -57.0*pi/180.0;
		phi = -48.0*pi/180.0+pi;
	}
	else
	{*/
		if ((in = fopen("inmkip.txt","r")) == NULL)
		{
			printf("cannot open \n");
			exit(0);
		}
		fscanf(in,"%lf %lf %lf %lf %lf", &psi, &phi, &psi2, &phi2, &psi3);
		psi = psi*pi/180.0;
		phi = phi*pi/180.0+pi;
		psi2 = psi2*pi/180.0;
		phi2 = phi2*pi/180.0+pi;
		psi3 = psi3*pi/180.0;
		fclose(in);
//	}

	g = calc_u(u, cos(psi), sin(psi), cos(phi), sin(phi),cos(psi2), sin(psi2), cos(phi2), sin(phi2), cos(psi3), sin(psi3), cos(alpha) ,sin(alpha), cos(beta) ,sin(beta));
	calc_s(s, cos(psi), sin(psi), cos(phi), sin(phi),cos(psi2), sin(psi2), cos(phi2), sin(phi2), cos(alpha) ,sin(alpha), cos(beta) ,sin(beta) ,q);
	printf("%lf %lf %lf\n",u[0],u[1],u[2]);
	printf("%lf %lf %lf\n",s[0],s[1],s[2]);
	if((out = fopen("inah.txt","w")) == NULL)
	{
		exit(1);
	}
	fprintf(out,"%lf %lf %lf\n",u[0],u[1],u[2]);
	fprintf(out,"%lf %lf %lf\n",s[0],s[1],s[2]);
	fclose(out);
}

double calc_u(double u[3], double cospsi, double sinpsi, double cosphi, double sinphi,double cospsi2, double sinpsi2, double cosphi2, double sinphi2,double cospsi3, double sinpsi3, double cosalpha, double sinalpha, double cosbeta, double sinbeta)
{
	int i,j;
	double g;
	double Talpha[3][3],Tbeta[3][3],Rpsi1[3][3],Rphi1[3][3],Rpsi2[3][3],Rphi2[3][3];
	double TaRs1[3][3],TaRs1Tb[3][3],TaRs1TbRh1[3][3],TaRs1TbRh1Ta[3][3],TaRs1TbRh1TaRs2[3][3],TaRs1TbRh1TaRs2Tb[3][3],TaRs1TbRh1TaRs2TbRh2[3][3],TaRs1TbRh1TaRs2TbRh2Ta[3][3];
	double vect[3];

	fomTmat(cosalpha, sinalpha, Talpha);
	fomTmat(cosbeta, sinbeta, Tbeta);
	fomRmat(cospsi, sinpsi, Rpsi1);
	fomRmat(cospsi2, sinpsi2, Rpsi2);
	fomRmat(cosphi, sinphi, Rphi1);
	fomRmat(cosphi2, sinphi2, Rphi2);

	mmult(Talpha,Rpsi1,TaRs1);
	mmult(TaRs1,Tbeta,TaRs1Tb);
	mmult(TaRs1Tb,Rphi1,TaRs1TbRh1);
	mmult(TaRs1TbRh1,Talpha,TaRs1TbRh1Ta);
	mmult(TaRs1TbRh1Ta,Rpsi2,TaRs1TbRh1TaRs2);
	mmult(TaRs1TbRh1TaRs2,Tbeta,TaRs1TbRh1TaRs2Tb);
	mmult(TaRs1TbRh1TaRs2Tb,Rphi2,TaRs1TbRh1TaRs2TbRh2);
	mmult(TaRs1TbRh1TaRs2TbRh2,Talpha,TaRs1TbRh1TaRs2TbRh2Ta);

	vect[0] = cosbeta;
	vect[1] = sinbeta*cospsi3;
	vect[2] = sinbeta*sinpsi3;

	mvult(TaRs1TbRh1TaRs2TbRh2Ta,vect,u);
	g = u[0]*TaRs1TbRh1TaRs2TbRh2Ta[0][0]+u[1]*TaRs1TbRh1TaRs2TbRh2Ta[1][0]+u[2]*TaRs1TbRh1TaRs2TbRh2Ta[2][0]-cosbeta;
	return g;
}

double calc_s(double s[3], double cospsi, double sinpsi, double cosphi, double sinphi,double cospsi2, double sinpsi2, double cosphi2, double sinphi2, double cosalpha, double sinalpha, double cosbeta, double sinbeta, double q[3])
{
	int i,j;
	double Talpha[3][3],Tbeta[3][3],Rpsi1[3][3],Rphi1[3][3],Rpsi2[3][3],Rphi2[3][3];
	double TaRs1[3][3],TaRs1Tb[3][3],TaRs1TbRh1[3][3],TaRs1TbRh1Ta[3][3],TaRs1TbRh1TaRs2[3][3],TaRs1TbRh1TaRs2Tb[3][3],TaRs1TbRh1TaRs2TbRh2[3][3],TaRs1TbRh1TaRs2TbRh2Ta[3][3];
	double s1[3],s2[3];

	fomTmat(cosalpha, sinalpha, Talpha);
	fomTmat(cosbeta, sinbeta, Tbeta);
	fomRmat(cospsi, sinpsi, Rpsi1);
	fomRmat(cospsi2, sinpsi2, Rpsi2);
	fomRmat(cosphi, sinphi, Rphi1);
	fomRmat(cosphi2, sinphi2, Rphi2);

	mmult(Talpha,Rpsi1,TaRs1);
	mmult(TaRs1,Tbeta,TaRs1Tb);
	mmult(TaRs1Tb,Rphi1,TaRs1TbRh1);
	mmult(TaRs1TbRh1,Talpha,TaRs1TbRh1Ta);
	mmult(TaRs1TbRh1Ta,Rpsi2,TaRs1TbRh1TaRs2);
	mmult(TaRs1TbRh1TaRs2,Tbeta,TaRs1TbRh1TaRs2Tb);
	mmult(TaRs1TbRh1TaRs2Tb,Rphi2,TaRs1TbRh1TaRs2TbRh2);
	mmult(TaRs1TbRh1TaRs2TbRh2,Talpha,TaRs1TbRh1TaRs2TbRh2Ta);

	mvult(TaRs1TbRh1,q,s1);
	mvult(TaRs1TbRh1TaRs2TbRh2,q,s2);

	for (i=0;i<3;++i)
	{
		s[i] = q[i]+s1[i]+s2[i];
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
		m1v2[i] = 0.0;
	}

	for (i=0;i<3;++i)
	{
		for (j=0;j<3;++j)
		{
			m1v2[i] += m1[i][j]*v2[j];
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


