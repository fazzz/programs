#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ASMD.h"
#include "efunc.h"

int pickMDdata(FILE *mddata, FILE *out){
  int i;
  char name;

  for (i=0;i<numstep;){
    fscanf(mddata,"%s",&name);
    if (strncmp(name,"steps",5)) {
      fscanf(mddata,"%d",&i);
    }
    if (strncmp(name,"total",5)) {
      fscanf(mddata,"%d",&temp);
      fscanf(mddata,"%d",i);
    }
    if (strncmp(name,"T_kelvin",8)) {
      ;
    }
    if (strncmp(name,"toal_energy",11)) {
      ;
    }
    if (strncmp(name,"toal_energy_vertial",19)) {
      ;
    }
    if (strncmp(name,"kinetic_energy8",15)) {
      ;
    }
    if (strncmp(name,"kinetic_energy_vertial",22)) {
      ;
    }
    if (strncmp(name,"potential_energy",15)) {
      ;
    }
    if (strncmp(name,"potential_energy_vertial",23)) {
      ;
    }
    else { 
      ;
    }
  }
}
