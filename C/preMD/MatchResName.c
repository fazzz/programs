#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Vis_MD.h"

// 標準残基データ中にインプットした残基名はあるか
int MatchResName(char *name) {
  int h,i;
  char name1[20],name2[20];
  
  h = hash(name);
  elementoftable *element;
  
  for (element = table[h]; element != NULL; element = element->next){
    for(i=0;element->nameofthisresidue[i]!='\0';++i)
      name1[i]=element->nameofthisresidue[i];
    for(i=0;*name!='\0';++i){
      name2[i]=*(name);
      ++name;
    }
    name2[i]='\0';
    if (name2[0] == element->nameofthisresidue[0] && name2[1] == element->nameofthisresidue[1] && name2[2] == element->nameofthisresidue[2]) {
      return 1;
    }
    if (strcmp(name1,name2) == 0) {
      return 1;
    }
  }

  return -1;
}
