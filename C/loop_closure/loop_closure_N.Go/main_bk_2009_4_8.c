#include <stdio.h>
#include <string.h>
#include <math.h>

#define pi 3.14159265

double calc_gfunc(double u[3], double cospsi1, double sinpsi1, double cospsi2, double sinpsi2, double cosphi1, double sinphi1, double cosphi2, double sinphi2, double alpha, double beta);
void calc_r(double alpha, double beta, double psi1, double s[3], double rou, double eta, double r[3]);
void mmult(double m1[3][3], double m2[3][3], double m1m2[3][3]);
void fomTmat(double cos, double sin, double mat[3][3]);
void fomRmat(double cos, double sin, double mat[3][3]);
int calc_phi1(double r[3], double w, double rou, double eta, double cosphi1[2], double sinphi1[2],double qz[2][3]);
int calc_phi2(double x,double w,double rou,double eta,double alpha,double beta,double cosphi2[2],double sinphi2[2]);
void calc_psi2(double w,double rou,double eta,double alpha,double beta,double qz[2][3],double cosphi2[2],double sinphi2[2], double cospsi2[2][2], double sinpsi2[2][2]);

int main(int argc, char *argv[])
{
	char outname[100];
	int i,j,k,l,m,n,test1,test2,len,numans;
	double x,y,z;
	double alpha,beta,rou,eta;
	double u[3],s[3];
	double psi1,r[3],w;
	double cospsi2[2][2],sinpsi2[2][2],cosphi1[2],sinphi1[2],cosphi2[2],sinphi2[2];
	double qz[2][3];
	double test3[4000];
	double psi1ans[4000],psi2ans[4000],phi1ans[4000],phi2ans[4000];
	double rans[4000],xans[4000];

	double g[4000][4];

	FILE *in,*out;
//  caluculation of alpha , beta
//	get_initial_values(pep1, pep2);

	alpha = 9.0*pi/180.0;
	beta = 70.5*pi/180.0;

	rou = 3.519;
//	eta = 1.436;
	eta = 1.242;
	if (argc == 0)
	{
		exit(1);
	}

	if ((in = fopen(argv[1],"r")) == NULL)
	{
		printf("cannot open %s\n",argv[1]);
		exit(0);
	}

	fscanf(in,"%lf %lf %lf",&x,&y,&z);
	u[0] = x;
	u[1] = y;
	u[2] = z;
	fscanf(in,"%lf %lf %lf",&x,&y,&z);
	s[0] = x;
	s[1] = y;
	s[2] = z;

	fclose(in);

//	u[0] = 0.059;
//	u[1] =-0.852;
//	u[2] =-0.521;
//	s[0] = 8.821;
//	s[1] = 1.842;
//	s[2] = -3.101;

	for (i=0;i<3600;++i)
	{
		test3[i] = 0;
		for (j=0;j<4;++j)
			g[i][j] = 0.0;
	}

	m = 0;
	n = 0;
	for (i=0;i<3600;++i)
	{
		psi1 = (double)i*pi/1800.0;

		// equ 21
		calc_r(alpha, beta, psi1, s, rou, eta, r);
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

			rans[n] = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
			xans[n] = r[0];
			++n;
			l=0;
			for (j=0;j<2;++j)
			{
				for (k=0;k<2;++k)
				{
					g[i][l] = calc_gfunc(u,cos(psi1),sin(psi1),cospsi2[j][k],sinpsi2[j][k],cosphi1[j],sinphi1[j],cosphi2[k],sinphi2[k],alpha,beta);
					++l;
					if (-0.001 <= g[i][l] && 0.001 >= g[i][l])
					{
						psi1ans[m] = psi1;
						psi2ans[m] = acos(cospsi2[j][k]);
						phi1ans[m] = acos(cosphi1[j]);
						phi2ans[m] = acos(cosphi2[k]);
					}
				}
			}
		}
	}
	numans = m;

	for (i=0;i<3600;++i)
	{
		if (test3[i] == 1)
			printf("%-5d %5.3lf %5.3lf %5.3lf %5.3lf\n",i,g[i][0],g[i][1],g[i][2],g[i][3]);
	}

	len = strlen(argv[1])+6;
	sprintf(outname,"g%s.txt",argv[1]);
	outname[len] = '\0';


	if ((out = fopen(outname,"w")) == NULL)
	{
		printf("cannot open g.txt\n");
		exit(0);
	}

	fprintf(out,"  psi1      g1      g2     g3     g4\n",i,g[i][0],g[i][1],g[i][2],g[i][3]);
	for (i=0;i<3600;++i)
	{
		psi1 = (double)i*pi/1800.0;
		if (test3[i] == 1)
			fprintf(out,"%-5.3lf %5.3lf %5.3lf %5.3lf %5.3lf\n",psi1*180.0/pi,g[i][0],g[i][1],g[i][2],g[i][3]);
	}

	fclose(out);

	sprintf(outname,"a%s.txt",argv[1]);
	outname[len] = '\0';

	if ((out = fopen(outname,"w")) == NULL)
	{
		printf("cannot open g.txt\n");
		exit(0);
	}

	for (i=0;i<numans;++i)
	{
		fprintf(out,"%-d: %5.1lf %5.1lf %5.1lf %5,1lf\n",i,psi1ans[i]*180.0/pi,psi2ans[i]*180.0/pi,phi1ans[i]*180.0/pi,phi2ans[i]*180.0/pi);
	}

	fclose(out);

	sprintf(outname,"x%s.txt",argv[1]);
	outname[len] = '\0';

	if ((out = fopen(outname,"w")) == NULL)
	{
		printf("cannot open g.txt\n");
		exit(0);
	}

	for (i=0;i<n;++i)
	{
		fprintf(out,"%5.4lf %5.4lf\n",xans[i],rans[i]);
	}

	fclose(out);

	return 0;
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
	mmult(TaRs1TbRh1TaRs2Tb,Talpha,TaRs1TbRh1TaRs2TbRh2Ta);

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

int calc_phi1(double r[3], double w, double rou, double eta, double cosphi1[2], double sinphi1[2],double qz[2][3])
{
	double D;

	D = r[1]*r[1]+r[2]*r[2]-w*w;
	if (D >= 0.0)
	{
		cosphi1[0] = (r[1]*w+r[2]*sqrt(D))/(r[1]*r[1]+r[2]*r[2]);
		sinphi1[0] = (r[2]*w-r[1]*sqrt(D))/(r[1]*r[1]+r[2]*r[2]);
		cosphi1[1] = (r[1]*w-r[2]*sqrt(D))/(r[1]*r[1]+r[2]*r[2]);
		sinphi1[1] = (r[2]*w+r[1]*sqrt(D))/(r[1]*r[1]+r[2]*r[2]);
		qz[0][0] = r[0]-rou;
		qz[0][1] = r[0]-eta;
		qz[0][2] = D;
		qz[1][0] = r[0]-rou;
		qz[1][1] = r[0]-eta;
		qz[1][2] = -D;
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
	cosphi2[1] = (rou*cos(beta)-(x-rou)*cos(alpha)-(w-eta)*sin(alpha))/(eta*sin(beta));
	if (-1.0 <= cosphi2[0] && cosphi2[0] <= 1.0)
	{
		sinphi2[0] =  sqrt(1-cosphi2[0]*cosphi2[0]);
		sinphi2[1] = -sqrt(1-cosphi2[0]*cosphi2[0]);
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

	for (i=0;i<2;++i)
	{
		for (j=0;j<2;++j)
		{
			cospsi2[i][j] = (
				             (rou*sin(beta)+eta*cosphi2[j]*cos(beta))*(-qz[i][0]*sin(alpha)+qz[i][1]*cos(alpha))
				            +eta*qz[i][2]*sinphi2[j]
				            )
				            /(
				            (-qz[i][0]*sin(alpha)+qz[i][1]*cos(alpha))*(-qz[i][0]*sin(alpha)+qz[i][1]*cos(alpha))+qz[i][2]*qz[i][2]
				             );
			sinpsi2[i][j] = (
				             (rou*sin(beta)+eta*cosphi2[j]*cos(beta))*qz[i][2]
				            +eta*(-qz[i][0]*sin(alpha)+qz[i][1]*cos(alpha))*sinphi2[j]
				            )
				            /(
				            (-qz[i][0]*sin(alpha)+qz[i][1]*cos(alpha))*(-qz[i][0]*sin(alpha)+qz[i][1]*cos(alpha))+qz[i][2]*qz[i][2]
				             );
		}
	}
}
