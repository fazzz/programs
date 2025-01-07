#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include "genrand.h"

// 一様乱数の生成を種種の方法で行う関数郡

/////////////////////////////////////////////////
/////////////////////////////////////////////////

// 線形合同法を用いて一様乱数の生成を行う関数郡

// ”種”の初期化
void init_rnd_lcm(unsigned long x)
{
  seed = x;
}

// 0 以上 ULONG_MAX 以下の整数の
// 一様乱数の生成を行う関数
unsigned long irnd_lcm(void)
{
  seed = seed * 1566083941UL + 1;
  return seed;
}

// 0 以上 1 以下の実数の
// 一様乱数の生成を行う関数
double rnd_lcm(void)
{
  return (1.0/(ULONG_MAX+1.0))*irnd_lcm();
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////

// rand() を用いて
// 0 以上 1 以下の実数の
// 一様乱数の生成を行う関数
double rnd(void)
{
	double U;

    U = (1.0/(RAND_MAX + 1.0)) * (rand() + 0.5);

	return U;
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
