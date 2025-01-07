#!/bin/gawk -f

BEGIN{
  i=0
  printf("step pote(dihed) pote(e1-4) pote(LJ1-4) pote(e1-5) pote(LJ1-5) pote(d+1-4+1-5) \n")
}


$1 ~ /(dihedral_energy)/{
    ++i
    printf("%d %lf ",i,$3) #> "energy_profile.txt"
}

$1 ~ /(elect_energy)/{
    printf("%lf ",$3) #> "energy_profile.txt"
}

$1 ~ /(VDW_energy)/{
    printf("%lf ",$3) #> "energy_profile.txt"
}

$1 ~ /(1_4_elect_energy)/{
    printf("%lf ",$3) #> "energy_profile.txt"
}

$1 ~ /(1_4_VDW_energy)/{
    printf("%lf\n",$3) #> "energy_profile.txt"
}


