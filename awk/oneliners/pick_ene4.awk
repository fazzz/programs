#!/bin/gawk -f

BEGIN{
  i=0
#  printf("step ene enev kine kinev pote potev")
  printf("step ene enev kine kinev pote potev\n")
}

$1 ~ /(toal_energy)/{
  ++i
#    printf("%d %lf ",i,$3) #> "energy_profile.txt"
  if ((i % intv)==0) {
    printf("%d %lf ",i,$3) #> "energy_profile.txt"
#    printf("%lf ",$3) #> "energy_profile.txt"
  }
}

$1 ~ /(toal_vertial_energy)/{
    if ((i % intv)==0) {
	printf("%lf ",$3) #> "energy_profile.txt"
    }
}

##################################################
# $1 ~ /(kinetic_energy6)/{			 #
#       printf("%lf ",$3) > "energy_profile.txt" #
# }						 #
##################################################

$1 ~ /(kinetic_energy8)/{
    if ((i % intv)==0) {
	printf("%lf ",$3) #> "energy_profile.txt"
    }
}

$1 ~ /(kinetic_energy_vertial)/{
    if ((i % intv)==0) {
	printf("%lf ",$3) #> "energy_profile.txt"
    }
}

$1 ~ /potential_energy_real/{
    if ((i % intv)==0) {
      printf("%lf ",$3) #> "energy_profile.txt"
    }
}

$1 ~ /(potential_energy_vertial)/{
    if ((i % intv)==0) {
	printf("%lf\n",$3) #> "energy_profile.txt"
    }
}

