#BEGIN{printf("match=(name \"CA\" and ( ")}
#
#{if(NR==1){printf("resnr %d ",$1)}else{printf(" or resnr %d ",$1)}}
#
#END{printf("));\nmatch;\n")}

{printf("%d \n",$1)}

