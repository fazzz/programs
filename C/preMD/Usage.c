#include <stdio.h>

#include "Vis_MD.h"

// 間違ったインプットに警告を出す関数
void usage (char *progname)
{
  printf ("usage: input the sequence\n");
  printf ("example: %s -h -s TYR GLY GLY PHE MET\n");
  printf ("example: %s -l -s TYR GLY GLY PHE MET\n");
  printf ("example: %s -h -i MetEnk.clt -f MetEnk.seq\n");
  printf ("example: %s -l -i MetEnk.clt -f MetEnk.seq\n");
}
