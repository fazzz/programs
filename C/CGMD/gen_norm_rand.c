#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gener.h"
#include "genrand.h"

// ���� ave �A���U vari �̈ꎟ�����K���z�̐������s���֐�
double nrmd_Box_Muller_method(int num, double ave, double vari)
{
  static int sw = 0;
  static double t, u;

  double rnd_1;
  double rnd_2;

  // ���`�����@���g��
  if(num == 1)
  {
    rnd_1 = rnd_lcm();
    rnd_2 = rnd_lcm();
  }
  // �W���̊֐� rand() ���g��
  else
  {
    rnd_1 = rnd();
    rnd_2 = rnd();
  }

  if (sw == 0)
  {
    sw = 1;
    t = sqrt(-2.0* vari* log(1.0-rnd_1));
    u = 2.0*PI*rnd_2;
    return ave + t * cos(u);
  }
  else
  {
    sw = 0;
    return ave + t* sin(u);
  }
}
