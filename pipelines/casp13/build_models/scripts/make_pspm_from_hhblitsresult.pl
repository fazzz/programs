#!/bin/perl

require "/home/tsukasa/casp_database/forte/scripts/a3mprocess.pl";


my $infile = $ARGV[0];
my $outa3m = $infile.".a3m";
my $outhhm = $infile.".hhm";
my $outpspm = $infile.".pspmhh";
my $outascii_psi = $infile.".asciihhpsi";
my $outpspm_psi = $infile.".pspmhhpsi";
my $db = "/home/tsukasa/casp_database/hh/uniclust30_2017_10/uniclust30_2017_10";


# if(-f $outa3m || -f $outhhm || -f $outpspm){
# 	print "Outputfile already exists.\n";
# 	die;
# }

#system("/home/tsukasa/lcl/bin/hhblits -n 3 -cpu 1 -all -i $infile -o $infile.hhr -oa3m $outa3m -ohhm $outhhm -d $db");
system("/home/tsukasa/lcl/bin/hhblits -n 3 -cpu 20 -i $infile -o $infile.hhr -oa3m $outa3m -ohhm $outhhm -d $db");



my $outfile = $outa3m.".fas$$";
my $prffile = $outpspm;
my $infile = $outa3m;
my $cou = 0;




open(IN,$infile);
while(my $ss = <IN>){
	if($ss =~ />/){
		$cou++;
	}
}
close(IN);

if($cou < 6000){
	# a3mprocess::a3m_to_fas($infile,$outfile);
	system("python /home/tsukasa/casp_database/forte/scripts/a3m_to_fas.py $infile $outfile");
	system("java -jar /home/tsukasa/casp_database/forte/scripts/FastaToProfile.jar -i $outfile -o $prffile");
	
	
	
#	unlink($outfile);
}else{
	open(OUT,"> $outfile.cddin");
	open(IN,$infile);
	my $snum = 0;
	while(my $ss = <IN>){
		if($ss =~ />([^\s]+)/){
			print OUT ">".$1."___$snum\n";
			$snum++;
		}else{
			$ss =~ s/[^A-Za-z]//g;
			print OUT uc $ss."\n";
		}
	}
	close(IN);
	close(OUT);
	system("/home/tsukasa/casp13_preparation/cdhit-4.6.8/cd-hit -i $outfile.cddin -o $outfile.cddout -c 0.9");
	my %store;
	open(IN,$outfile.".cddout");
	while(my $ss = <IN>){
		if($ss =~ />([^\s]+)/){
			my $name = $1;
			$name =~ s/___[0-9]+[\s]*$//;
			$store{$name} = 100;
			#print $name."\n";
		}
	}
	close(IN);
	open(OUT,"> $outfile.profin");
	open(IN,$infile);
	$snum = 0;
	my $flag = 0;
	my $head = <IN>;
	print OUT $head;
	$head = <IN>;
	print OUT $head;
	while(my $ss = <IN>){
		if($ss =~ />([^\s]+)/){
			my $name = $1;
			if(defined $store{$name}){
				$flag = 1;
			}else{
				$flag = 0;
			}
		}
		if($flag == 1){
			print OUT $ss;
		}
	}
	close(IN);
	close(OUT);
	
	# a3mprocess::a3m_to_fas("$outfile.profin",$outfile);
	system("python /home/tsukasa/casp_database/forte/scripts/a3m_to_fas.py $outfile.profin $outfile");
	system("java -Xmx512g -jar /home/tsukasa/casp_database/forte/scripts/FastaToProfile.jar -i $outfile -o $prffile");
	#system("java -jar /home/tsukasa/casp_database/forte/scripts/FastaToProfile.jar -i $outfile -o $prffile");
#	unlink($outfile);
	unlink($outfile.".profin");
	unlink($outfile.".cddin");
	unlink($outfile.".cddout");
	unlink($outfile.".cddout.clstr");
}


my $psiblast = "/home/tsukasa/psiblastexb250 ";
my $nrdb = "/home/tsukasa/casp_database/blast_nr/nr";

my $seqcoung = 0;
open(IN,$outfile);
while(my $ss = <IN>){
	if($ss =~ />/){
		$seqcount++;
	}
}
close(IN);
if($seqcount > 0){
	system("$psiblast -num_iterations 1 -num_threads 1 -num_descriptions 100000 -db $nrdb -in_msa $outfile -out_ascii_pssm $outascii_psi > /dev/null");
}else{
#	system("$psiblast -num_iterations 3 -num_descriptions 100000 -num_alignments 0 -db $nrdb -query $infile -out_ascii_pssm $outpspm_psi > /dev/null");

}
unlink($outfile);
makepspm($outascii_psi,$outpspm_psi);

sub makepspm{
	my @ret;
	my $pssm = $_[0];
	my $outname = $_[1];
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
}
