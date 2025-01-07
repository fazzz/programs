BEGIN{n=0;flag=0;} #BEGIN{n=1;flag=0;}

{if($1==":"){flag=1}else{flag=0}}

$2!="-"{n=n+1;if(flag==1){print n}}
