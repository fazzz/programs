BEGIN{
    printf("[ MATCH ] \n")
    n=0

}

{
    printf("%4d ",$3);
    n=n+1;
}

n%15==0{
    n=0;
    printf("\n");
}
