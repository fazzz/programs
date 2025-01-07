#!/bin/gawk -f

BEGIN{
  i=0
  printf("step ene dihede es LJ 14es 14LJ\n")
}

$1 ~ /(E_t)/{
    ++i
    printf("%d %lf ",i,$3) 
}

$1 ~ /p_tot/{
    printf("%lf ",$3) 
}

$1 ~ /p_nat/{
    printf("%lf ",$3) 
}

$1 ~ /p_rep/{
    printf("%lf ",$3) 
}

$1 ~ /p_14LJ/{
    p_14L=$3;
    printf("%lf ",$3) 
}

$1 ~ /p_dih/{
    p_dih=$3;
    printf("%lf ",$3) 
}

$1 ~ /p_ang/{
    p_ang=$3;
    printf("%lf ",$3) 
}

$1 ~ /p_bon/{
    p_bon=$3;
    printf("%lf %lf\n",$3,p_14L+p_dih+P_ang+p_bon) 
}
