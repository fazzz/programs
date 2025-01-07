package fastaloader;


use strict;
use warnings;






sub load_seq{
	my $filename = $_[0];
	if(!-f $filename){
		print "nofile named ".$filename."\n";
		die;
		
		
	}
	open(IN,$filename);
	my @ret;
	my $name = "dummy";
	my $seq ="";
	while(my $ss = <IN>){
		if($ss =~ /^[\s]*>/){
			if(length($seq) > 0){
				push(@ret,$name.$seq);
				$seq = "";
			}
			$name = $ss;
		}else{
			if($ss =~ />/){
				print "illigal line.\n";
				print $ss;
				die;
			}
			if($ss =~ /^[\r\n]*$/){
			}else{
				$seq .= $ss;
			}
		}
		
		
	}
	close(IN);
	
	if(length($seq) > 0){
		push(@ret,$name.$seq);
		$seq = "";
	}
	return \@ret;
}

sub load_seq_hash{
	my $filename = $_[0];
	my @array = @{&load_seq($filename)};
	my %ret;
	foreach my $aa(@array){
		
		if($aa =~ /^[\s]*>([^\r\n]+)[\r\n](.+)$/s){
			my $name = $1;
			my $seq = $2;
			my @nname = split(/>/,$name);
			foreach my $nn(@nname){
				my $gn = $nn;
				$gn =~ s/[\s]//g;
				if(length($gn) > 0){
					$nn =~ /^[\s]*([^\s]+)/;
					my $kname = $1;
					my $kdesc = "";
					if($nn =~ /^[\s]*[^\s]+[\s]+(.+)/){
						$kdesc = $1;
					}
					if(defined $ret{$kname}){
						print "Duplicated sequence name in $filename\n".$kname."\n".$seq."\n";
					}
					#$ret{$kname} = ">".$kname."  ".$kdesc."\n".$seq;
					$ret{$kname} = ">".$name."\n".$seq;
					last;#please comment out if there are multiple names 
				}
			}
		}else{
			print "Abnormal sequence found.\n".$aa;
			
		}
	}
	return \%ret;
}


sub get_seq_fasta{
	my $h = toHash($_[0]);
	return ${$h}{"seq"};
}

sub get_desc_fasta{
	my $h = toHash($_[0]);
	return ${$h}{"desc"};
}

sub get_name_fasta{
	my $h = toHash($_[0]);
	return ${$h}{"name"};
}

sub toHash{
	my $seq = $_[0];
	my @lin = split(/[\r\n]+/,$seq);
	my %ret;
	if($lin[0] =~ /^[\s]*>[\s]*([^\s]+)/){
		$ret{"name"} = $1;
	}
	if($lin[0] =~ /^[\s]*>[\s]*[^\s]+[\s]+([^\s][^\r\n]+)/){
		$ret{"desc"} = $1;
	}else{
		$ret{"desc"} = "";
	}
	my $s = "";
	for(my $ii = 1;$ii <= $#lin;$ii++){
		$s .= $lin[$ii];
	}
	$ret{"seq"} = $s;
	return \%ret;
}

sub multi_fasta_to_phylip{
	my @array = @{$_[0]};
	my @name;
	my @seq;
	my $maxlen = 0;;
	foreach my $aa(@array){
		my $h = toHash($aa);
		my $n = ${$h}{"name"};
		my $s = ${$h}{"seq"};
		$s =~ s/[\s]//g;
		my $le = length($s);
		if($maxlen < $le){
			$maxlen = $le;
		}
		push(@name,$n);
		push(@seq,$s);
	}
	my @ret;
	for(my $ii = 0;$ii <= $#name;$ii++){
		my $n = $name[$ii];
		my $s = $seq[$ii];
		my $l = length($s);
		for(;$l < $maxlen;$l++){
			$s .= "-";
		}
		push(@ret,$n." ".$s);
		
	}
	my $sn = $#name+1;
	unshift(@ret,$sn." ".$maxlen);
	return join("\n",@ret);
	
}





1;