#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "physics.h"

// 0-1 の範囲の一様分布
#define rnd() (1.0 / (RAND_MAX+1.0))*rand()
#define POOLSIZE 97

static double pool[POOLSIZE];

void nrmd_set_initial_velo(int nNumClut);
double randw(void);
void init_randw(void);

// 二面角速度初期値の設定を行う関数
void set_initial_velo(void)
{
	int nNumClut;

	for(nNumClut=1;nNumClut<prot.DOF;++nNumClut)
	{
		nrmd_set_initial_velo(nNumClut);
	}
}

// 二面角速度の代入を行う関数
void nrmd_set_initial_velo(int nNumClut)
{
	double t_x, u_x;

	double x,y,z;

	double momentum_clust;

	double kT;

	momentum_clust = clust[nNumClut].Inertia_clust[2][2]*(1.660e-27)*1.0e-20;
//	momentum_clust = clust[nNumClut].Inertia_clust_total*(1.660e-27)*1.0e-20;

	kT = T_Kelvin_Initial*k_B_J;

	x = rnd();

	y = log(1-rnd());

	z = -2 * kT /momentum_clust * y;

	t_x = sqrt(z);
	u_x = 2 * PI * rnd();

	clust[nNumClut].ddihedang[0]=t_x*cos(u_x)*1.0e-12/10;

}

// 一様分布の発生を行う関数_1
void init_randw(void)
{
	int i;

	for(i = 0;i < POOLSIZE; i++)
	{
		pool[i] = rnd();
	}
}

// 一様分布の発生を行う関数_2
double randw(void)
{
	static int i = POOLSIZE -1;
	double r;

	i = (int)(POOLSIZE * pool[i]);
	r = pool[i];
	pool[i] = rnd();

	if (r > 0.5)
		return 1.0;
	else
		return -1.0;
}
