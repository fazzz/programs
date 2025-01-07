package a3mprocess;


use strict;
use warnings;


#require "c:\\dir_important\\fastaloader.pl";
#require "/home/tsukasa/casp_database/forte/scripts/fastaloader.pl";
require "/home/yamamori/work/casp13/scripts_forte/fastaloader.pl";



sub a3m_to_fas{
	my $inname = $_[0];
	my $outname = $_[1];
	
	my @inseq = @{fastaloader::load_seq($inname)};
	
	if(-f $outname){
		print "$outname already exists. Please check.";
		die;
	}
	my %insertx_before;
	my %dupname;
	my %processed;
	for(my $ii = 1;$ii <= $#inseq;$ii++){
		my $ss = $inseq[$ii];
		my $seq = fastaloader::get_seq_fasta($ss);
		
		$seq =~ s/[\s]//g;
		my @aa = split(//,$seq);
		
		$processed{$ii} = \@aa;
		my @buff;
		my $cou = 0;
		my $icou = 0;#ƒMƒƒƒbƒv‚Ì”
		for(my $ii = 0;$ii <= $#aa;$ii++){
			if($aa[$ii] =~ /[a-z]/){
				$icou++;
				if(!defined $insertx_before{$cou}){
					$insertx_before{$cou} = 0;
				}
				if($insertx_before{$cou} < $icou){#
					$insertx_before{$cou} = $icou;
				}
			}else{
				$cou++;
				$icou = 0;
			}
		}
	}
	
	
	my @plines;
	for(my $ii = 0;$ii <= $#inseq;$ii++){
		my $ss = $inseq[$ii];
		my $seq = fastaloader::get_seq_fasta($ss);
		my $name = fastaloader::get_name_fasta($ss);
		my $desc = fastaloader::get_desc_fasta($ss);
		
		my $aa;
		if(!defined $processed{$ii}){
			$seq =~ s/[\s]//g;
			my @caa = split(//,$seq);
			$aa = \@caa;
		}else{
			$aa =$processed{$ii};
		}
		my @buff;
		my $cou = 0;
		my $icou = 0;
		my $alen = $#{$aa};
		for(my $ii = 0;$ii <= $alen;$ii++){
			if(${$aa}[$ii] =~ /[a-z]/){
				$icou++;
				push(@buff,${$aa}[$ii]);
			}else{
				if(!defined $insertx_before{$cou}){
					$insertx_before{$cou} = 0;
				}
				while($icou < $insertx_before{$cou}){
					push(@buff,"-");
					$icou++;
				}
				push(@buff,${$aa}[$ii]);
				$icou = 0;
				$cou++;
			}
			
			if($ii == $alen){
				if(!defined $insertx_before{$cou}){
					$insertx_before{$cou} = 0;
				}
				while($icou < $insertx_before{$cou}){
					push(@buff,"-");
					$icou++;
				}
			}
		}
		
		if(defined $dupname{$name}){
			push(@plines, ">$name.($dupname{$name})  $desc\n");
			push(@plines, join("",@buff));
			push(@plines, "\n\n");
		}else{
			push(@plines,  ">$name  $desc\n");
			push(@plines,  join("",@buff));
			push(@plines,  "\n\n");
			
		}
		$dupname{$name}++;
	}
	
	open(OUT,">".$outname);
	print OUT join("",@plines);
	close(OUT);
}

sub fas_to_a3m{
	my $inname = $_[0];
	my $outname = $_[1];
	
	my @inseq_tmp = @{fastaloader::load_seq($inname)};
	my @inseq;
	for(my $ii = 0;$ii <= $#inseq_tmp;$ii++){
		my $seq = fastaloader::get_seq_fasta($inseq_tmp[$ii]);
		$seq =~ s/[^A-Za-z]//g;
		if(length($seq) > 0){
			push(@inseq,$inseq_tmp[$ii]);
		}
	}
	@inseq_tmp = ();
	
	open(OUT,">".$outname);
	my %mask;
	for(my $ii = 0;$ii == 0;$ii++){
		my $ss = $inseq[$ii];
		my $name = fastaloader::get_name_fasta($ss);
		my $seq = fastaloader::get_seq_fasta($ss);
		$seq =~ s/[\s]//g;
		my @aa = split(//,$seq);
		for(my $ii = 0;$ii <= $#aa;$ii++){
			if($aa[$ii] =~ /[A-Za-z]/){#‰½‚ç‚©‚ÌƒAƒ~ƒmŽ_‚ª‘¶Ý‚·‚éê‡‚Í 1AƒMƒƒƒbƒv‚Æ‚È‚Á‚Ä‚¢‚éê‡‚Í -1
				$mask{$ii} = 1;
			}else{
				$mask{$ii} = -1;
			}
		}
	}

	for(my $ii = 0;$ii <= $#inseq;$ii++){
		my $ss = $inseq[$ii];
		my $seq = fastaloader::get_seq_fasta($ss);
		my $name = fastaloader::get_name_fasta($ss);
		my $desc = fastaloader::get_desc_fasta($ss);
		
		$seq =~ s/[\s]//g;
		
		
		my @aa = split(//,$seq);
		my @buff;
		my $cou = 0;
		for(my $ii = 0;$ii <= $#aa;$ii++){
			if($aa[$ii] =~ /[A-Za-z]/){
				if($mask{$ii} == 1){
					push(@buff,uc $aa[$ii]);
				}else{#æ“ª”z—ñ‚ªƒMƒƒƒbƒv‚Ìê‡¬•¶Žš
					push(@buff,lc $aa[$ii]);
				}
			}else{
				if($mask{$ii} == 1){
					push(@buff,$aa[$ii]);
					#•’ÊƒMƒƒƒbƒv‚¾‚¯
				}else{
					#æ“ª”z—ñ‚à‚±‚Ì”z—ñ‚àƒMƒƒƒbƒv‚Å‚ ‚éê‡‚Í‰½‚à‘}“ü‚µ‚È‚¢
				}
			}
		}
		print OUT ">$name  $desc\n";
		print OUT join("",@buff);
		print OUT "\n\n";
		
		
	}
	close(OUT);
}



1;
