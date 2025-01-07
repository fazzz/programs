use strict;
use warnings;

#require "../scripts/fastaloader.pl";

require "/home/tsukasa/casp_database/forte/scripts/blastparser.pl";
require "/home/tsukasa/casp_database/forte/scripts/alignment_recover.pl";



#PSIBLAST ‚ÌŒ‹‰Ê‚©‚ç MSA ‚ðì‚é
#CDHIT,MAFFT,CDD ‚ÌFASTA ‚ª•K—v

my $psiblast = "/home/tsukasa/psiblastexb250";
my $nrdb = "/home/tsukasa/casp_database/blast_nr/nr";


my $infilename = $ARGV[0];
my $outfilename = $infilename.".pspmpb";
my $asciifilename = $infilename.".asciipb";

# if(-f $outfilename || -f $asciifilename){
# 	print "Outputfile already exists.\n";
# 	die;
# }
#waitUntil();#‘¼‚Ìì‹Æ‚ªI‚í‚é‚Ü‚Å‘Ò‚Â

my @allseq = @{fastaloader::load_seq($infilename)};

my $nname = fastaloader::get_name_fasta($allseq[0])."_".$$;
$nname =~ s/[^a-zA-Z0-9\.\-]/_/g;
my $blastoutname =  "$nname.blastout";
my $tmpname = $infilename;


system("$psiblast -num_descriptions 100000  -num_iterations 5 -save_pssm_after_last_round -query $tmpname -db $nrdb -out $blastoutname -out_ascii_pssm $asciifilename ");
makepspm($asciifilename,$outfilename);
unlink($blastoutname);


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
