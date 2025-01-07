#include <stdio.h>
#include <stdlib.h>

#include "Vis_MD.h"

elementoftable *LURData(char nameo[10],int create,residuedata dataofthisresdiue) {
  int h,i;

  h = hash(nameo);
  elementoftable *thiselement;

  for (thiselement = table[h]; thiselement != NULL; thiselement = thiselement->next) {
    if ((i=strcmp(nameo, thiselement->nameofthisresidue)) == 0) {
      return thiselement;
    }
  }
  if (create == 0){
    thiselement = (elementoftable *)emalloc(sizeof(elementoftable));
    for (i=0;nameo[i]!='\0';++i)
      thiselement->nameofthisresidue[i] = nameo[i];
    thiselement->nameofthisresidue[i]='\0';
    thiselement->dataofthisresdiue = dataofthisresdiue;
    thiselement->next = table[h];
    table[h] = thiselement;
  }
  return thiselement;
}

