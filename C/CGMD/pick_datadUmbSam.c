#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "gener.h"
#include "ABA.h" // #include "ABA_multimer.h"
#include "physics.h"
#include "force.h"
#include "MD.h"
#include "BD.h"
#include "UmsSan.h"

// setudouno�̎擾���s���֐�
void cndUmbSanscan(void)
{
  int i;
  double x;
  int y;
 
  FILE* input;

  if ((input=fopen("UmbSan.in","r")) == NULL)
  {
    V_K_US =0.0;

    dihed_ref_US=0.0;
  }
  else
  {
    // ��ʊp���̌��q�̔ԍ��̎擾
    for (i=0;i<5;++i)
    {
      fscanf(input,"%d",&y);
      atom_pair_US[i] = y;
    }


    fscanf(input,"%lf",&x);
    V_K_US = x;

    fscanf(input,"%lf",&x);
    dihed_ref_US = x;

    fclose(input);
  }

}
