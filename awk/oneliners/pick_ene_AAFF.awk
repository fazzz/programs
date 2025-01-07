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

$1 ~ /p_es/{
    printf("%lf ",$3) 
}

$1 ~ /p_LJ/{
    printf("%lf ",$3) 
}

$1 ~ /p_14es/{
    printf("%lf ",$3) 
}

$1 ~ /p_14LJ/{
    printf("%lf ",$3) 
}

$1 ~ /p_dih/{
    printf("%lf ",$3) 
}

$1 ~ /p_ang/{
    printf("%lf ",$3) 
}

$1 ~ /p_bon/{
    printf("%lf \n",$3) 
}
