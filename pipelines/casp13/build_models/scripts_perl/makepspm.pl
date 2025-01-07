my @ret;
my $pssm = $ARGV[0];
my $outname = $ARGV[1];
if(-f $outname){
	#die;
}

my @aapt = split(//,"ARNDCQEGHILKMFPSTWYV");
open(OUT,">".$outname)or die;
open(IN,$pssm) or die;
while(my $ss = <IN>){
	if($ss =~ /Last/){
		last;
	}
}
my $head = <IN>;

print OUT "\nLast position-specific scoring matrix computed.\n";

print OUT "            A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V \n";


$head =~ s/^[\s]+//;
$head =~ s/[\s]+$//;
my @sp = split(/ +/,$head);
my @similarity;
while(my $ss = <IN>){
	$ss =~ s/[\r\n]//g;
	if($ss =~ /Lambda/){
		last;
	}
	if($ss =~ /^([\s]*[^\s]+[\s]*[^\s]+[\s]*)([^\s].+)/){
		my $h = $1;
		my $p = $2;
		my $aacode = "";
		
		if($h =~ /^[\s]*([^\s]+)[\s]*([^\s]+)/){
			$aacode = $2;
			if($aacode !~ /[A-Z]/){
				print $aacode."\n";
				die;
			}
			
		}else{
			die;
		}	
		my @pt = split(/[\s]+/,$p);
		my @rt;
		my $allzero = 1;
#			for(my $ii = 0; $ii < 20;$ii++){
		for(my $ii = 20; $ii < 40;$ii++){
			if(!defined $pt[$ii] || $pt[$ii] =~ /[a-zA-Z]/|| $pt[$ii] !~ /[0-9\.]/){
				print $pssm."\n";
				print $ss."\n";
				die;
			}
			if($pt[$ii] != 0){
				$allzero = 0;
			}
			push(@rt,$pt[$ii]);
		}
		if($allzero == 1){
			@rt = ();
			for(my $ii = 0; $ii < 20;$ii++){
				if($aapt[$ii] eq $aacode){
					push(@rt,100);
				}else{
					push(@rt,0);
				}
			}
		}
		
		print OUT $h." ".join(" ",@rt);
		print OUT "\n";
	}
}
print OUT "\n";
close(IN);
close(OUT);