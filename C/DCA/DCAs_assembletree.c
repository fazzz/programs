
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "DCA.h"
#include "EF.h"
#include "LA.h"

AST *DCAs_assemble_tree_testcase(int numclut,int *numtree){
  AST *assemble_tree;

  *numtree=13;

  assemble_tree=(AST *)gcemalloc(sizeof(AST)*(*numtree));

  assemble_tree[0].leafflag=0;
  assemble_tree[0].num=1;
  assemble_tree[0].left=2;
  assemble_tree[0].right=3;
  assemble_tree[0].refright=5;
  assemble_tree[0].refleft=1;

  assemble_tree[1].leafflag=0;
  assemble_tree[1].num=2;
  assemble_tree[1].left=4;
  assemble_tree[1].right=5;
  assemble_tree[1].refright=3;
  assemble_tree[1].refleft=1;

  assemble_tree[2].leafflag=0;
  assemble_tree[2].num=3;
  assemble_tree[2].left=6;
  assemble_tree[2].right=13;
  assemble_tree[2].refright=7;
  assemble_tree[2].refleft=5;

  assemble_tree[3].leafflag=0;
  assemble_tree[3].num=4;
  assemble_tree[3].left=7;
  assemble_tree[3].right=8;
  assemble_tree[3].refright=2;
  assemble_tree[3].refleft=1;

  assemble_tree[4].leafflag=0;
  assemble_tree[4].num=5;
  assemble_tree[4].left=9;
  assemble_tree[4].right=10;
  assemble_tree[4].refright=4;
  assemble_tree[4].refleft=3;

  assemble_tree[5].leafflag=0;
  assemble_tree[5].num=6;
  assemble_tree[5].left=11;
  assemble_tree[5].right=12;
  assemble_tree[5].refright=6;
  assemble_tree[5].refleft=5;

  assemble_tree[6].leafflag=1;
  assemble_tree[6].right=1;
  assemble_tree[6].num=7;

  assemble_tree[7].leafflag=1;
  assemble_tree[7].right=2;
  assemble_tree[7].num=8;

  assemble_tree[8].leafflag=1;
  assemble_tree[8].right=3;
  assemble_tree[8].num=9;

  assemble_tree[9].leafflag=1;
  assemble_tree[9].right=4;
  assemble_tree[9].num=10;

  assemble_tree[10].leafflag=1;
  assemble_tree[10].right=5;
  assemble_tree[10].num=11;

  assemble_tree[11].leafflag=1;
  assemble_tree[11].right=6;
  assemble_tree[11].num=12;

  assemble_tree[12].leafflag=1;
  assemble_tree[12].right=7;
  assemble_tree[12].num=13;

  return assemble_tree;
}

AST *DCAs_assemble_tree_testcase_2(int numclut,int *numtree){
  AST *assemble_tree;

  *numtree=8;

  assemble_tree=(AST *)gcemalloc(sizeof(AST)*(*numtree));

  assemble_tree[0].leafflag=0;
  assemble_tree[0].num=1;
  assemble_tree[0].left=2;
  assemble_tree[0].right=8;
  assemble_tree[0].refright=5;
  assemble_tree[0].refleft=1;

  assemble_tree[1].leafflag=0;
  assemble_tree[1].num=2;
  assemble_tree[1].left=5;
  assemble_tree[1].right=6;
  assemble_tree[1].refright=2;
  assemble_tree[1].refleft=1;

  assemble_tree[2].leafflag=0;
  assemble_tree[2].num=3;
  assemble_tree[2].left=4;
  assemble_tree[2].right=8;
  assemble_tree[2].refright=5;
  assemble_tree[2].refleft=3;

  assemble_tree[3].leafflag=0;
  assemble_tree[3].num=4;
  assemble_tree[3].left=7;
  assemble_tree[3].right=8;
  assemble_tree[3].refright=4;
  assemble_tree[3].refleft=3;

  assemble_tree[4].leafflag=1;
  assemble_tree[4].num=5;
  assemble_tree[4].right=1;

  assemble_tree[5].leafflag=1;
  assemble_tree[5].num=6;
  assemble_tree[5].right=2;

  assemble_tree[6].leafflag=1;
  assemble_tree[6].num=7;
  assemble_tree[6].right=3;

  assemble_tree[7].leafflag=1;
  assemble_tree[7].num=8;
  assemble_tree[7].right=4;

  assemble_tree[7].leafflag=1;
  assemble_tree[7].num=9;
  assemble_tree[7].num=5;

  return assemble_tree;
}
