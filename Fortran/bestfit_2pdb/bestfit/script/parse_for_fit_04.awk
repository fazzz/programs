BEGIN{FIELDWIDTHS = "6 5 1 4 1 3 1 1 4 1 3 8 8 8 6 6 6 4"
    nr=0;
    nnr=0
    flag=0;
}

$1~/ATOM/ && $4~/CA/{if($9!=nr){nr=$9;nnr=nnr+1;printf("%s %4d %4d\n",$4,NR,nnr)}}

#$1~/ATOM/ && $4~/CA/{printf("%s %4d %4d\n",$4,NR,nnr)} #$1~/ATOM/ && $4~/CA/{printf("%s %4d %4d\n",$4,$2,nnr)}
