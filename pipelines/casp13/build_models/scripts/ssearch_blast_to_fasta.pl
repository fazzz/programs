use strict;
use warnings;



require "blastparser.pl";

if(!defined $ARGV[1]){
	#print "usage:perl ssearch_blast_to_fasta.pl <infilename> <outfiledir> <threshold(optional>\n";
	#exit(0);
}
my $fname = $ARGV[0];
my $outfiledir = $ARGV[1];
my $threshold_ = 0.001;


if(!-d $outfiledir){
	mkdir($outfiledir);
}
$outfiledir .= "/";
$outfiledir =~ s/\/+$/\//;

if(defined $ARGV[2]){
	$threshold_ = $ARGV[2];
}
my $buff = "";
my @ret;
open(B_IN,$fname)or die;
my $flag = 0;
while(my $ss = <B_IN>){
	if($ss =~ /^Query=/){
		if($ss =~ /^Query=[\s]*$/){
			$ss =~ s/Query=/Query=dummylabel/g;
		}
		$buff = $ss;
		last;
	}
}
while(my $ss = <B_IN>){
	if($ss =~ /^Query=/){
		if(length($buff) > 0){
			if($buff =~ /Query[\s]*=[\s]*([^\s]+)/){
				my $sname = $1;
				# my $fname = $sname;
				# $fname =~ s/[^A-Za-z0-9\.]/_/g;
				
				# my $outfilename = $outfiledir.$fname.".blastout.fas";
				# my $outfilename = $outfiledir.$fname.".fas";
				my $outfilename = $fname.".fas";
				print $outfilename;
				
				my $res = blastparser::get_result_array($buff);
				save_fasta($outfilename,$sname,$res,$threshold_);
				
			}
			$buff = "";
		}
		
		if($ss =~ /^Query=[\s]+$/){
			$ss =~ s/^Query=/Query=dummylabel/;
		}
	#	if(rand() > 0.999){
	#		last;
	#	}
	}
	$buff .= $ss;
	
	
}
if($buff =~ /Query=/){
	if($buff =~ /Query[\s]*=[\s]*([^\s]+)/){
		my $sname = $1;
		# my $fname = $sname;
		# $fname =~ s/[^A-Za-z0-9\.]/_/g;
		
		# my $outfilename = $outfiledir.$fname.".blastout.fas";
		# my $outfilename = $outfiledir.$fname.".fas";
		my $outfilename = $fname.".fas";
		print $outfilename;
		
		my $res = blastparser::get_result_array($buff);
		save_fasta($outfilename,$sname,$res,$threshold_);
		
	}
}
close(B_IN);



sub save_fasta{
	my $outfilename = $_[0];
	my $queryname = $_[1];
	my $res = $_[2];
	my $threshold = $_[3];
	
	if(length($outfilename) > 600){
		print substr($outfilename,0,100)."\n is too long.";
		die;
	}
	my $rlen = @{$res};
	my %printed;
	my $maxquery = "";
	my @pline;
	for(my $ii = 0;$ii < $rlen;$ii++){
		my @rresarray = @{blastparser::parse_blast_result(${$res}[$ii])};
		foreach my $ra(@rresarray){
			my %rres = %{$ra};
			if($rres{"Expect"} <= $threshold){
				my $name = $rres{"name"};
				my $seq = $rres{"s_seq"};
				my $qseq = $rres{"q_seq"};
				$qseq =~ s/[^A-Za-z]//g;
				if(length($maxquery) < length($qseq)){
					$maxquery = $qseq;
				}
				if(!defined $printed{$name}){
					$printed{$name} = 1;
				}else{
					$name .= "_".$printed{$name};
					$printed{$name}++;
				}
				push(@pline,">".$name."  ".$rres{"Expect"}."\n");
				push(@pline,$seq."\n");
			}
		}
	}

	open(POUT,"> ".$outfilename);
	print POUT ">".$queryname."\n";
	print POUT $maxquery;
	print POUT "\n";
	foreach my $ll(@pline){
		print POUT $ll;
	}
	
	close(POUT);
	
}



