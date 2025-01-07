#!/bin/gawk -f

BEGIN{
  i=0
  printf("step ene enev kine kinev pote potev dihede\n")
}

$1 ~ /(toal_energy)/{
  ++i
    printf("%d %lf ",i,$3) #> "energy_profile.txt"
}

$1 ~ /(toal_vertial_energy)/{
      printf("%lf ",$3) #> "energy_profile.txt"
}

##################################################
# $1 ~ /(kinetic_energy6)/{			 #
#       printf("%lf ",$3) > "energy_profile.txt" #
# }						 #
##################################################

$1 ~ /(kinetic_energy8)/{
      printf("%lf ",$3) #> "energy_profile.txt"
}


$1 ~ /(kinetic_energy_vertial)/{
      printf("%lf ",$3) #> "energy_profile.txt"
}

$1 ~ /potential_energy_real/{
      printf("%lf ",$3) #> "energy_profile.txt"
}

$1 ~ /(potential_energy_vertial)/{
    printf("%lf ",$3) #> "energy_profile.txt"
}

$1 ~ /(dihedral_energy)/{
    printf("%lf\n",$3) #> "energy_profile.txt"
}
