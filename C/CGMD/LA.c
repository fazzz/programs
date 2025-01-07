#include <stdio.h>
#include <stdlib.h>

#include "LA.h"

void mmult(double m1[MAXN][MAXN], double m2[MAXN][MAXN], double m1m2[MAXN][MAXN], int n) {
  int i,j,k;

  for (i=0;i<n;++i)
    for (j=0;j<n;++j)
      m1m2[i][j] = 0.0;

  for (i=0;i<n;++i) {
    for (j=0;j<n;++j) {
      for (k=0;k<n;++k)	{
	m1m2[i][j] += m1[i][k]*m2[k][j];
      }
    }
  }
}
