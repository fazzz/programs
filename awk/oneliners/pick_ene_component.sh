#!/bin/gawk -f

BEGIN{
  i=0
#  printf("step ene dihede es LJ 14es 14LJ\n")
}

$1 ~ /(potential_energy_real)/{
    ++i
    printf("%d %lf ",i,$3) #> "energy_profile.txt"
}

$1 ~ /dihedral_energy/{
      printf("%lf ",$3) #> "energy_profile.txt"
}

$1 ~ /elect_energy/{
    printf("%lf ",$3) #> "energy_profile.txt"
}

$1 ~ /VDW_energy/{
    printf("%lf ",$3) #> "energy_profile.txt"
}

$1 ~ /1_4_elect_energy/{
    printf("%lf ",$3) #> "energy_profile.txt"
}

$1 ~ /1_4_VDW_energy/{
    printf("%lf\n",$3) #> "energy_profile.txt"
}
