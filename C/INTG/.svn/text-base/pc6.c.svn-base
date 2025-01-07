#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "INTG.h"
#include "LA.h"

void intg_set(void) {
  int i,j;

  GearsConst5[0] = 3.0/16.0;
  GearsConst5[1] = 251.0/360.0;
  GearsConst5[2] = 1.0;
  GearsConst5[3] = 11.0/18.0;
  GearsConst5[4] = 1.0/6.0;
  GearsConst5[5] = 1.0/60.0;
  
  GearsConst2[0] = 3.0/20.0;
  GearsConst2[1] = 251.0/360.0;
  GearsConst2[2] = 1.0;
  GearsConst2[3] = 11.0/18.0;
  GearsConst2[4] = 1.0/6.0;
  GearsConst2[5] = 1.0/60.0;
  
  GearsConst4[0] = 251.0/720.0;
  GearsConst4[1] = 1.0;
  GearsConst4[2] = 11.0/12.0;
  GearsConst4[3] = 1.0/3.0;
  GearsConst4[4] = 1.0/24.0;
	
  for (i=0;i<6;++i)
    for (j=0;j<6;++j)
      if (i != j)
	Telar_Matrix[i][j] = 0.0;
      else
	Telar_Matrix[i][i] = 1.0;
	
  Telar_Matrix[0][1] = 1.0;
  Telar_Matrix[0][2] = 1.0;
  Telar_Matrix[0][3] = 1.0;
  Telar_Matrix[0][4] = 1.0;
  Telar_Matrix[0][5] = 1.0;
  Telar_Matrix[1][2] = 2.0;
  Telar_Matrix[1][3] = 3.0;
  Telar_Matrix[1][4] = 4.0;
  Telar_Matrix[1][5] = 5.0;
  Telar_Matrix[2][3] = 3.0;
  Telar_Matrix[2][4] = 6.0;
  Telar_Matrix[2][5] = 10.0;
  Telar_Matrix[3][4] = 4.0;
  Telar_Matrix[3][5] = 10.0;
  Telar_Matrix[4][5] = 5.0;

}

void intg_pc5(int flag,double prev[6],double corv[6], double dt,double acc) {
  int i,j;

  if(flag==1){  
    for (i=0;i<6;++i)
      prev[i] = 0.0;
    for (i=0;i<6;++i)
      for (j=0;j<6;++j)
	prev[i] += Telar_Matrix[i][j]*corv[j];
  }
  else {
    for (i=0;i<6;++i)
      corv[i]=0.0;
    for (i=0;i<6;++i)
      corv[i]=prev[i]+GearsConst5[i]*(0.5*dt*dt*acc-prev[2]);
  }
}

void intg_pc4(int flag,double prev[5],double corv[5], double dt,double vel) {
  int i,j;

  if(flag==1){  
    for (i=0;i<5;++i)
      prev[i] = 0.0;
    for (i=0;i<5;++i)
      for (j=0;j<5;++j)
	prev[i] += Telar_Matrix[i][j]*corv[j];
  }
  else {
    for (i=0;i<5;++i)
      corv[i]=0.0;
    for (i=0;i<5;++i)
      corv[i]=prev[i]+GearsConst4[i]*(dt*vel-prev[1]);
  }
}
