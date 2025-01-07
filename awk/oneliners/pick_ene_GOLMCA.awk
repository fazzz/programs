#!/bin/gawk -f

BEGIN{
  i=0
  printf("step p_tot p_nat p_rep p_dih p_ang p_bon ")
}

$1 ~ /(p_tot)/{
  ++i
  printf("\n%d %lf ",i,$3) #> "energy_profile.txt"
}

$1 ~ /(p_nat)/{
  printf(" %lf ",$3) #> "energy_profile.txt"
}

$1 ~ /(p_rep)/{
  printf(" %lf ",$3) #> "energy_profile.txt"
}

$1 ~ /(p_dih)/{
  printf(" %lf ",$3) #> "energy_profile.txt"
}

$1 ~ /(p_ang)/{
  printf(" %lf ",$3) #> "energy_profile.txt"
}

$1 ~ /(p_bon)/{
  printf(" %lf",$3) #> "energy_profile.txt"
}
