#include <stdio.h>

#include "Vis_MD.h"

// �Ԉ�����C���v�b�g�Ɍx�����o���֐�
void usage (char *progname)
{
  printf ("usage: input the sequence\n");
  printf ("example: %s -h -s TYR GLY GLY PHE MET\n");
  printf ("example: %s -l -s TYR GLY GLY PHE MET\n");
  printf ("example: %s -h -i MetEnk.clt -f MetEnk.seq\n");
  printf ("example: %s -l -i MetEnk.clt -f MetEnk.seq\n");
}
