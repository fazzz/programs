#!/bin/gawk -f

BEGIN{
  i=0
  printf("step enet dihed natt repult dihedt \n")
}

$1 ~ /(potential_energy)/{
  ++i
  printf("\n%d %lf ",i,$3) #> "energy_profile.txt"
}

$1 ~ /(native_contact_energy)/{
  printf(" %lf ",$3) #> "energy_profile.txt"
}

$1 ~ /(non_native_repul_energy)/{
  printf(" %lf ",$3) #> "energy_profile.txt"
}

$1 ~ /(dihedral_angle_energy)/{
  printf(" %lf ",$3) #> "energy_profile.txt"
}
