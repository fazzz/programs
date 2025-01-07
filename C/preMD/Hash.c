#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


#include "Vis_MD.h"

#define MULTIPLIER /*39*//*37*/31
//#define MULTIPLIER 1000

unsigned int hash(char * str) {
  unsigned int h;
  unsigned char *p;

  h = 0;
  for (p = (unsigned char *)str; *p != '\0'; p++)
    h = MULTIPLIER * h + *p;
  h = h%NHASH;
  //printf("h=%d hh=%d %s\n",h, hh,str);
  /****************************/
  /* int i,h=0;		      */
  /* int L=strlen(str);	      */
  /* 			      */
  /* for (i=0;i<L;++i) {      */
  /*   h=MULTIPLIER*h+str[i]; */
  /* }			      */
  /* h = h % NHASH;	      */
  /****************************/
  return h;
  //  printf("h=%d  %s\n",h, str);
}
