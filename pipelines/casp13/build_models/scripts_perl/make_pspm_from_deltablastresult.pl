use strict;
use warnings;

#require "../scripts/fastaloader.pl";

#require "/home/tsukasa/casp_database/forte/scripts/blastparser.pl";
#require "/home/tsukasa/casp_database/forte/scripts/alignment_recover.pl";

require "/home/yamamori/work/casp13/scripts_forte/blastparser.pl";
require "/home/yamamori/work/casp13/scripts_forte/alignment_recover.pl";

#DELTABLAST の結果から MSA を作る
#CDHIT,MAFFT,CDD のFASTA が必要


my $cdhit = "/home/tsukasa/casp13_preparation/cdhit-4.6.8/cd-hit";
my $javapath = "java";
#my $fastatoprofile = "/home/tsukasa/casp_database/forte/scripts/FastaToProfile.jar";
my $fastatoprofile = "/home/yamamori/work/casp13/scripts_forte/FastaToProfile.jar";
#my $deltablast = "/home/tsukasa/lcl/bin/deltablast";
my $deltablast = "/home/yamamori/opt2/ncbi-blast-2.8.1+/bin/deltablast";

# my $pseudodb = "/home/tsukasa/casp_database/forte/data/newdb/hoge";
my $pseudodb = "/home/tsukasa/casp_database/blast_nr/nr";
my $rpsdb = "/home/tsukasa/casp_database/blast_cdd_delta/cdd_delta";
my $cddfasdir = "/home/tsukasa/casp_database/cdd_fasta/";#cdd の FASTA




my $infilename = $ARGV[0];
my $outfilename = $infilename.".pspmdb";
my $asciifilename = $infilename.".asciidb";
my $outmsaname = $infilename.".msa".$$;


# if(-f $outfilename || -f $asciifilename || -f $outmsaname){
# 	print "Outputfile already exists.\n";
# 	die;
# }
#waitUntil();#他の作業が終わるまで待つ

my @allseq = @{fastaloader::load_seq($infilename)};

my $nname = fastaloader::get_name_fasta($allseq[0])."_".$$;
$nname =~ s/[^a-zA-Z0-9\.\-]/_/g;
my $blastoutname =  "$nname.blastout";
my $tmpname = $infilename;


system("$deltablast -show_domain_hits -num_alignments 500 -num_iterations 1 -query $tmpname -db $pseudodb -rpsdb $rpsdb -out $blastoutname -out_ascii_pssm $asciifilename");
my $queryname = fastaloader::get_name_fasta($allseq[0]);
my $queryseq = fastaloader::get_seq_fasta($allseq[0]);

add_one_seq_in_msa_with_deltablastresult($queryname,$queryseq,$blastoutname,$cddfasdir,$outmsaname,0.05001);
unlink($blastoutname);


system("$javapath -Xmx24g -jar $fastatoprofile -i $outmsaname -o $outfilename -pseudocount blosum62");
unlink($outmsaname);






sub add_one_seq_in_msa_with_deltablastresult{
	my $queryname = $_[0];
	my $queryseq = $_[1];
	my $blastresult = $_[2];
	my $msafasdir =  $_[3];#あらかじめアラインされたファイル
	
	my $outfilename = $_[4];
	my $ethreshold = $_[5];
	
	
	my @blast = @{blastparser::get_multi_result_array($blastresult)};
	my @hash = @{blastparser::parse_blast_result($blast[0])};
	
	@blast = ();
	
	
	
	
	
	my %pairs;
	
	
	foreach my $bb(@hash){
		my %bres = %{$bb};
		
		if($bres{"Expect"} > $ethreshold){
			next;
		}
		my %hs;
		$hs{$queryname} = $bres{"q_seq"};
		$hs{"lcl|consensus"} = $bres{"s_seq"};
		$hs{"indexoffset>"} = $bres{"q_start"};
		
		if(!defined $bres{"desc"}){
			foreach my $kk(keys %bres){
				print $kk."\t".$bres{$kk}."\n";
			}
			die;
		}
		#print $bres{"desc"}."\n";
		if($bres{"desc"} =~ /^([^\s,;]+)/){
			my $filename = $msafasdir.$1.".FASTA";#ファイルが無いかもしれない？？
			if(!-f $filename){
				print STDERR "no file $filename\n";
				die;
			}
			if(!defined $pairs{$filename}){
				my @tmp;
				$pairs{$filename} = \@tmp;
			}
			push(@{$pairs{$filename}},\%hs);
		}
	}
	
	my $msacou = 0;;;
	my @all_msa;
	foreach my $kkp(keys %pairs){
		my @fas = @{fastaloader::load_seq($kkp)};
		$msacou += ($#fas+1)*($#{$pairs{$kkp}}+1);
		push(@all_msa,@fas);
	}
	my %storename;
	
	if($msacou > 5999){
		my @pseq = @{reduce_alignment(\@all_msa)};
		foreach my $pp(@pseq){
			my $pn = fastaloader::get_name_fasta($pp);
			$storename{$pn} = 100;
		}
	}else{
		
		foreach my $pp(@all_msa){
			my $pn = fastaloader::get_name_fasta($pp);
			$storename{$pn} = 100;
		}
	}
	
	
	my $cou = 0;;
	
	my @msas;
	foreach my $kkp(keys %pairs){
	#	print $kkp."\n";
	}
	foreach my $kkp(keys %pairs){
		my @al = @{$pairs{$kkp}};
		foreach my  $aa(@al){
			my %msa;
			my %desc;
			my @names;
			my @fas = @{fastaloader::load_seq($kkp)};
			foreach my $ff(@fas){
				my $name = fastaloader::get_name_fasta($ff);
				my $seq = fastaloader::get_seq_fasta($ff);
				my $desc = fastaloader::get_desc_fasta($ff);
				if(defined $storename{$name} || $name =~ /consensus/ ){
					$msa{$name} = $seq;
					$desc{$name} = $desc;
					push(@names,$name);
				}
			}
			@fas = ();
			@hash = ();
			#print STDERR "$kkp###\n";
			my @tmppairs;
			push(@tmppairs,$aa);#二番目以降も入れる
			#my %msares = %{alignment_recover::add_one_seq_in_msa_($queryname,$queryseq,$pairs{$kkp},\%msa,1)};
			my %msares = %{alignment_recover::add_one_seq_in_msa_($queryname,$queryseq,\@tmppairs,\%msa,1)};
			my @ret;
			push(@ret,">".$queryname."\n".$msares{$queryname});
			foreach my $nn(@names){
				push(@ret,">".$nn." ".$desc{$nn}."\n".$msares{$nn});
			}
			push(@msas,\@ret);
			$cou++;
		}
	}
	my $stock_name;
	my $stock_desc;
	my $stock_aa;
	my $stock;
	if($#msas == -1){
		my @tmp = (">".$queryname."\n".$queryseq);
		$stock = \@tmp;
		
	}else{
		my @allmsa;
		foreach my $mm(@msas){
			my ($n,$d,$s) = fastaProcess($mm);
			my @g;
			push(@g,$n);
			push(@g,$d);
			push(@g,$s);
			
			push(@allmsa,\@g);
		}
		$stock = fastaMerge_2(\@allmsa);
	}
	
	
	
	
	if(defined $outfilename){
		open(OUT,">".$outfilename);
		
		print OUT join("\n",@{$stock});
		close(OUT);
	}
	return $stock;
}




sub fastaMerge_2{
	
	
	my @allsets = @{$_[0]};
	my @alldesc;
	my @allname;
	my @allseq;
	my %refseq;
	foreach my $aset(@allsets){
		my ($resa1,$resa3,$resa_aa) = (${$aset}[0],${$aset}[1],${$aset}[2]);
		
		my @fa_aa = @{$resa_aa};
		my @faname = @{$resa1};
		my @fadesc = @{$resa3};
		
		
		my $anum = $#fa_aa+1;
		
		
		my $s = join("",@{$fa_aa[0]});
		$s =~ s/[^A-Za-z]//g;
		my $reflen = length($s);
		
		for(my $ii = 0;$ii < $anum;$ii++){
			my @tmp;
			
			for(my $jj = 0;$jj < $reflen;$jj++){
				push(@tmp,"");
				push(@tmp,"");
			}
			push(@tmp,"");
			push(@allseq,\@tmp);
			if($ii == 0){
				$refseq{$#allseq} = 100;
			}
		}
		
		push(@allname,@faname);
		push(@alldesc,@fadesc);
	}
	my $offset = 0;
	
	my %maxlen;
	foreach my $aset(@allsets){
		my ($resa1,$resa3,$resa_aa) = (${$aset}[0],${$aset}[1],${$aset}[2]);
		my @fa_aa = @{$resa_aa};
		my @faname = @{$resa1};
		my @fadesc = @{$resa3};
		
		
		my $anum = $#fa_aa+1;
		
		
		my $reflen = $#{$fa_aa[0]}+1;
		my $columnpos = 0;
		for(my $jj = 0;$jj < $reflen;$jj++){
			my $matched = 0;
			if(${${$resa_aa}[0]}[$jj] =~ /[a-zA-Z]/){
				$columnpos++;
				$matched = 1;
			}
			for(my $ii = 0;$ii < $anum;$ii++){
				if($matched == 1){
					${$allseq[$offset+$ii]}[$columnpos*2-1] .= ${${$resa_aa}[$ii]}[$jj];
					$maxlen{$columnpos*2-1} = 1;
				}else{
					${$allseq[$offset+$ii]}[$columnpos*2] .= ${${$resa_aa}[$ii]}[$jj];
					if(!defined $maxlen{$columnpos*2}){
						$maxlen{$columnpos*2} = length(${$allseq[$offset+$ii]}[$columnpos*2]);
					}elsif($maxlen{$columnpos*2} < length(${$allseq[$offset+$ii]}[$columnpos*2])){
						$maxlen{$columnpos*2} = length(${$allseq[$offset+$ii]}[$columnpos*2]);
					}
				}
			}
		}
		$offset += $anum;
	}
	
	my @ret;
	for(my $ii = 0;$ii <= $#allseq;$ii++){
		if($ii != 0 && defined $refseq{$ii}){
			next;
		}
		my @seq;
		my @aseq = @{$allseq[$ii]};
		for(my $jj = 0; $jj <= $#aseq;$jj++){
			my $res = $aseq[$jj];
			push(@seq,$res);
			if(!defined $maxlen{$jj}){
				$maxlen{$jj} = 0;
			}
			if($maxlen{$jj} > length($res)){
				my $pl = $maxlen{$jj} - length($res);
				for(my $nn = 0;$nn < $pl;$nn++){
					push(@seq,"-");
				}
			}
		}
		 
		my $jseq = join("",@seq);
		my $name = $allname[$ii];
		my $desc = $alldesc[$ii];
		if($jseq =~ /[A-Za-z]/){
			push(@ret,">".$name."  ".$desc."\n".$jseq."\n");
		}
		
	}
	return \@ret;
}

sub fastaMerge_{

	my ($resa1,$resa3,$resa2) = ($_[0],$_[1],$_[2]);
	my ($resb1,$resb3,$resb2) = ($_[3],$_[4],$_[5]);
	
	my @fa_aa = @{$resa2};
	my @fb_aa = @{$resb2};
	my @faname = @{$resa1};
	my @fbname = @{$resb1};
	
	my @fadesc = @{$resa3};
	my @fbdesc = @{$resb3};
	
	#最初の配列の何番目のアミノ酸に A 組の何カラム目が当たるか
	my %seqa_map;
	my %gapa;#最初の配列の何番目のアミノ酸の前にギャップがあり、カラムが当てはまっているか
	
	
	if($faname[0] ne $fbname[0]){
		my $fs = $faname[0]."\n".join("",@{$fa_aa[0]});
		my $fb = $fbname[0]."\n".join("",@{$fb_aa[0]});
		print STDERR "First sequences are  different ".$fs."\n".$fb."\n";
		die;
	}
	
	my @refa = @{$fa_aa[0]};
	my $refseq_A = uc join("",@refa);
	$refseq_A =~ s/[\s\-\.]//g;
	if($refseq_A =~ /[^A-Za-z]/){
		print $refseq_A."  $& cannot be processed\n";
	}
	
	my $cou = 0;
	my $gapina = 0;
	for(my $ii = 0;$ii <= $#refa;$ii++){
		if($refa[$ii] =~ /[a-zA-Z]/){
			$seqa_map{$cou} = $ii;
			$cou++;
		}else{
			if(!defined $gapa{$cou}){
				my @tmp;
				$gapa{$cou} = \@tmp;
			}
			push(@{$gapa{$cou}},$ii);
			$gapina++;
		}
	}
	
	
	my %seqb_map;
	my @refb = @{$fb_aa[0]};
	my %gapb;#最初の配列の何番目のアミノ酸の前にギャップがあり、カラムが当てはまっているか
	my $refseq_B = uc join("",@refb);
	$refseq_B =~ s/[\s\-\.]//g;
	if($refseq_B =~ /[^A-Za-z]/){
		print STDERR $refseq_B."  $& cannot be processed\n";
		die;
	}
	if($refseq_A ne $refseq_B){
		print STDERR $refseq_A."\n\n".$refseq_B."  sequence_difference.\n";
		die;
	}
	$cou = 0;
	my $gapinb = 0;
	for(my $ii = 0;$ii <= $#refb;$ii++){
		if($refb[$ii] =~ /[a-zA-Z]/){
			$seqb_map{$cou} = $ii;
			$cou++;
		}else{
			if(!defined $gapb{$cou}){
				my @tmp;
				$gapb{$cou} = \@tmp;
			}
			push(@{$gapb{$cou}},$ii);
			$gapinb++;
		}
	}
	
	my $msalength = $cou+$gapina+$gapinb;
	my @resfas;
	for(my $ii= 0; $ii <= $#fa_aa+$#fb_aa;$ii++){
		my @tmp;
		push(@resfas,\@tmp);
	}
	
	
	if(1==0){#ギャップ部分は全部ずらす
		for(my $ii = 0;$ii <= $cou;$ii++){
			if(defined $gapa{$ii}){
				foreach my $pos(@{$gapa{$ii}}){
					my $row = 0;
					for(my $aa = 0;$aa <= $#fa_aa;$aa++){
						push(@{$resfas[$row]},${$fa_aa[$aa]}[$pos]);
						$row++;
					}
					for(my $bb = 0;$bb <= $#fb_aa;$bb++){
						push(@{$resfas[$row]},"-");
						$row++;
					}
				}
			}
			if(defined $gapb{$ii}){
				foreach my $pos(@{$gapb{$ii}}){
					my $row = 0;
					for(my $aa = 0;$aa <= $#fa_aa;$aa++){
						push(@{$resfas[$row]},"-");
						$row++;
					}
					for(my $bb = 0;$bb <= $#fb_aa;$bb++){
						push(@{$resfas[$row]},${$fb_aa[$bb]}[$pos]);
						$row++;
					}
				}
			}
			if($ii != $cou){
				my $row = 0;
				for(my $aa = 0;$aa <= $#fa_aa;$aa++){
					push(@{$resfas[$row]},${$fa_aa[$aa]}[$seqa_map{$ii}]);
					$row++;
				}
				for(my $bb = 0;$bb <= $#fb_aa;$bb++){
					push(@{$resfas[$row]},${$fb_aa[$bb]}[$seqb_map{$ii}]);
					$row++;
				}
			}
			
		}
	}else{#ギャップ部分のアラインメントはバラバラ
	
		for(my $ii = 0;$ii <= $cou;$ii++){
			if(defined  $gapa{$ii} && defined  $gapb{$ii}){
				my $gapsizea = @{$gapa{$ii}};
				my $gapsizeb = @{$gapb{$ii}};
				
				my $gapsize = $gapsizea;
				if($gapsize < $gapsizeb){
					$gapsize = $gapsizeb;				
				}
				
				for(my $pp = 0;$pp < $gapsize;$pp++){
					my $apos = -1;
					my $bpos = -1;
					if($pp < $gapsizea){
						$apos = ${$gapa{$ii}}[$pp];
					}
					
					if($pp < $gapsizeb){
						$bpos = ${$gapb{$ii}}[$pp];
					}
					
					
					my $row = 0;
					for(my $aa = 0;$aa <= $#fa_aa;$aa++){
						if($apos == -1){
							push(@{$resfas[$row]},"-");
						}else{
							push(@{$resfas[$row]},${$fa_aa[$aa]}[$apos]);
						}
						$row++;
					}
					for(my $bb = 0;$bb <= $#fb_aa;$bb++){
						if($bpos == -1){
							push(@{$resfas[$row]},"-");
						}else{
							push(@{$resfas[$row]},${$fb_aa[$bb]}[$bpos]);
						}
						$row++;
					}
				}
				
			}else{
				if(defined $gapa{$ii}){
					foreach my $pos(@{$gapa{$ii}}){
						my $row = 0;
						for(my $aa = 0;$aa <= $#fa_aa;$aa++){
							push(@{$resfas[$row]},${$fa_aa[$aa]}[$pos]);
							$row++;
						}
						for(my $bb = 0;$bb <= $#fb_aa;$bb++){
							push(@{$resfas[$row]},"-");
							$row++;
						}
					}
				}
				if(defined $gapb{$ii}){
					foreach my $pos(@{$gapb{$ii}}){
						my $row = 0;
						for(my $aa = 0;$aa <= $#fa_aa;$aa++){
							push(@{$resfas[$row]},"-");
							$row++;
						}
						for(my $bb = 0;$bb <= $#fb_aa;$bb++){
							push(@{$resfas[$row]},${$fb_aa[$bb]}[$pos]);
							$row++;
						}
					}
				}
			}
			if($ii != $cou){
				my $row = 0;
				for(my $aa = 0;$aa <= $#fa_aa;$aa++){
					push(@{$resfas[$row]},${$fa_aa[$aa]}[$seqa_map{$ii}]);
					$row++;
				}
				for(my $bb = 0;$bb <= $#fb_aa;$bb++){
					push(@{$resfas[$row]},${$fb_aa[$bb]}[$seqb_map{$ii}]);
					$row++;
				}
			}
			
		}
	
	}
	
	
	
	my %remseq;
	if(1 == 0){
		for(my $bb = 1;$bb <= $#fb_aa;$bb++){
			my $brow = $#fa_aa+1+$bb;
			for(my $aa = 1;$aa <= $#fa_aa;$aa++){
				my $arow = $aa;
				my ($id,$cov) = seq_ident_($resfas[$arow],$resfas[$brow]);
				if($id == -1){
					$remseq{$bb} = 100;
					last;
				}
				if($id == 1.0 && $cov > 0.8){
					my $r = aa_merge_($resfas[$arow],$resfas[$brow]);
					$resfas[$arow] = $r;
					$remseq{$bb} = 100;
					print "merged.\n";
					last;
				}
				#print $id."\t $cov\n";
			}
		}
	}
	
	
	my @ret_name;
	my @ret_desc;
	my @ret_aa;
	my $row = 0;
	my %used;
	
	for(my $aa = 0;$aa <= $#fa_aa;$aa++){
	
		my $dname = $faname[$aa];
		$dname =~ s/\.[0-9]+$//;
		
		my $sname = $faname[$aa];
		
		if(defined $used{$sname}){
			$sname = $dname.".".$used{$sname};
			if(defined $used{$dname}){
				$sname = $dname.".".$used{$dname};
				while(-f $sname){
					$sname = $dname.".".$used{$dname};
					$used{$dname}++;
				}
			}
		}
		$used{$sname}++;
		$used{$dname}++;
		
		push(@ret_name,$sname);
		push(@ret_desc,$fadesc[$aa]);
		push(@ret_aa,$resfas[$row]);
		$row++;
	}
	for(my $bb = 0;$bb <= $#fb_aa;$bb++){
		if(defined $remseq{$bb}){
		}elsif($bb == 0){#最初の配列をデバッグのために入れる
		#	push(@ret,">".$fbname[$bb]."  ".$fbdesc[$bb]."\n".join("",@{$resfas[$row]})."\n");
		}else{
			
			my $dname = $fbname[$bb];
			$dname =~ s/\.[0-9]+$//;
			
			my $sname = $fbname[$bb];
			
			if(defined $used{$sname}){
				$sname = $dname.".".$used{$sname};
				if(defined $used{$dname}){
					$sname = $dname.".".$used{$dname};
					while(-f $sname){
						$sname = $dname.".".$used{$dname};
						$used{$dname}++;
					}
				}
			}
			$used{$sname}++;
			$used{$dname}++;
				
			push(@ret_name,$sname);
			push(@ret_desc,$fbdesc[$bb]);
			push(@ret_aa,$resfas[$row]);
		}
		$row++;
	}
	
	return \@ret_name,\@ret_desc,\@ret_aa;
}





sub fastaMerge{
	my $fa = $_[0];
	my $fb = $_[1];
	
	my ($resa1,$resa3,$resa2) = fastaProcess($fa);
	my ($resb1,$resb3,$resb2) = fastaProcess($fb);
	
	my @fa_aa = @{$resa2};
	my @fb_aa = @{$resb2};
	my @faname = @{$resa1};
	my @fbname = @{$resb1};
	
	my @fadesc = @{$resa3};
	my @fbdesc = @{$resb3};
	
	#最初の配列の何番目のアミノ酸に A 組の何カラム目が当たるか
	my %seqa_map;
	my %gapa;#最初の配列の何番目のアミノ酸の前にギャップがあり、カラムが当てはまっているか
	
	
	if($faname[0] ne $fbname[0]){
		print STDERR "First sequences are  different ".${$fa}[0]."\n".${$fb}[0]."\n";
		die;
	}
	
	my @refa = @{$fa_aa[0]};
	my $refseq_A = uc join("",@refa);
	$refseq_A =~ s/[\s\-\.]//g;
	if($refseq_A =~ /[^A-Za-z]/){
		print $refseq_A."  $& cannot be processed\n";
	}
	
	my $cou = 0;
	my $gapina = 0;
	for(my $ii = 0;$ii <= $#refa;$ii++){
		if($refa[$ii] =~ /[a-zA-Z]/){
			$seqa_map{$cou} = $ii;
			$cou++;
		}else{
			if(!defined $gapa{$cou}){
				my @tmp;
				$gapa{$cou} = \@tmp;
			}
			push(@{$gapa{$cou}},$ii);
			$gapina++;
		}
	}
	
	
	my %seqb_map;
	my @refb = @{$fb_aa[0]};
	my %gapb;#最初の配列の何番目のアミノ酸の前にギャップがあり、カラムが当てはまっているか
	my $refseq_B = uc join("",@refb);
	$refseq_B =~ s/[\s\-\.]//g;
	if($refseq_B =~ /[^A-Za-z]/){
		print STDERR $refseq_B."  $& cannot be processed\n";
		die;
	}
	if($refseq_A ne $refseq_B){
		print STDERR $refseq_A."\n\n".$refseq_B."  sequence_difference.\n";
		die;
	}
	$cou = 0;
	my $gapinb = 0;
	for(my $ii = 0;$ii <= $#refb;$ii++){
		if($refb[$ii] =~ /[a-zA-Z]/){
			$seqb_map{$cou} = $ii;
			$cou++;
		}else{
			if(!defined $gapb{$cou}){
				my @tmp;
				$gapb{$cou} = \@tmp;
			}
			push(@{$gapb{$cou}},$ii);
			$gapinb++;
		}
	}
	
	my $msalength = $cou+$gapina+$gapinb;
	my @resfas;
	for(my $ii= 0; $ii <= $#fa_aa+$#fb_aa;$ii++){
		my @tmp;
		push(@resfas,\@tmp);
	}
	
	
	for(my $ii = 0;$ii <= $cou;$ii++){
		if(defined $gapa{$ii}){
			foreach my $pos(@{$gapa{$ii}}){
				my $row = 0;
				for(my $aa = 0;$aa <= $#fa_aa;$aa++){
					push(@{$resfas[$row]},${$fa_aa[$aa]}[$pos]);
					$row++;
				}
				for(my $bb = 0;$bb <= $#fb_aa;$bb++){
					push(@{$resfas[$row]},"-");
					$row++;
				}
			}
		}
		if(defined $gapb{$ii}){
			foreach my $pos(@{$gapb{$ii}}){
				my $row = 0;
				for(my $aa = 0;$aa <= $#fa_aa;$aa++){
					push(@{$resfas[$row]},"-");
					$row++;
				}
				for(my $bb = 0;$bb <= $#fb_aa;$bb++){
					push(@{$resfas[$row]},${$fb_aa[$bb]}[$pos]);
					$row++;
				}
			}
		}
		if($ii != $cou){
			my $row = 0;
			for(my $aa = 0;$aa <= $#fa_aa;$aa++){
				push(@{$resfas[$row]},${$fa_aa[$aa]}[$seqa_map{$ii}]);
				$row++;
			}
			for(my $bb = 0;$bb <= $#fb_aa;$bb++){
				push(@{$resfas[$row]},${$fb_aa[$bb]}[$seqb_map{$ii}]);
				$row++;
			}
		}
		
	}
	
	
	
	my %remseq;
	for(my $bb = 1;$bb <= $#fb_aa;$bb++){
		my $brow = $#fa_aa+1+$bb;
		for(my $aa = 1;$aa <= $#fa_aa;$aa++){
			my $arow = $aa;
			my ($id,$cov) = seq_ident_($resfas[$arow],$resfas[$brow]);
			if($id == 1.0 && $cov > 0.8){
				my $r = aa_merge_($resfas[$arow],$resfas[$brow]);
				$resfas[$arow] = $r;
				$remseq{$bb} = 100;
				print "merged.\n";
				last;
			}
			#print $id."\t $cov\n";
		}
	}
	
	
	
	if(1== 0){
	
		my $snum_threshold = 3000;
		my @entropy_base;
		my $eff=1000;
		my $sf = $resfas[0];
		for(my $uu = 0;$uu <= $#{$sf};$uu++){
			if(${$sf}[$uu] =~ /[a-zA-Z]/){
				my %tmp;
				$entropy_base[$uu] = \%tmp;
			}else{
			}
		}
		my $sb = 0;
		for(my $aa = 0;$aa <= $#fa_aa;$aa++){
			my $sb = $resfas[$aa];
			for(my $uu = 0;$uu <= $#{$sf};$uu++){
				if(${$sb}[$uu] =~ /[A-Za-z]/){
					if(defined $entropy_base[$uu]){
						$entropy_base[$uu]{${$sb}[$uu]}++;
					}
				}
			}
		}
		my @ao = split(//,"ACDEFGHIKLMNPQRSTVWY");
		############途中
		
		
		for(my $bb = 1;$bb <= $#fb_aa;$bb++){
			my $brow = $#fa_aa+1+$bb;
			my $sb = $resfas[$bb];
			for(my $uu = 0;$uu <= $#{$sb};$uu++){
				if(${$sb}[$uu] =~ /[A-Za-z]/){
					if($entropy_base[$uu] > -1 && $entropy_base[$uu] < $eff){
						$snum_threshold++;
					}
				}
			}
		}
	}
	
	
	my @ret;
	my $row = 0;
	my %used;
	for(my $aa = 0;$aa <= $#fa_aa;$aa++){
		
		my $dname = $faname[$aa];
		$dname =~ s/\.[0-9]+$//;
		
		my $sname = $faname[$aa];
		
		if(defined $used{$sname}){
			$sname = $dname.".".$used{$sname};
			if(defined $used{$dname}){
				$sname = $dname.".".$used{$dname};
				while(-f $sname){
					$sname = $dname.".".$used{$dname};
					$used{$dname}++;
				}
			}
		}
		$used{$sname}++;
		$used{$dname}++;
		
		push(@ret,">".$sname."  ".$fadesc[$aa]."\n".join("",@{$resfas[$row]})."\n");
		$row++;
	}
	for(my $bb = 0;$bb <= $#fb_aa;$bb++){
		if(defined $remseq{$bb}){
		}elsif($bb == 0){#最初の配列をデバッグのために入れる
		#	push(@ret,">".$fbname[$bb]."  ".$fbdesc[$bb]."\n".join("",@{$resfas[$row]})."\n");
		}else{
		
		
			my $dname = $fbname[$bb];
			$dname =~ s/\.[0-9]+$//;
			
			my $sname = $fbname[$bb];
			
			if(defined $used{$sname}){
				$sname = $dname.".".$used{$sname};
				if(defined $used{$dname}){
					$sname = $dname.".".$used{$dname};
					while(-f $sname){
						$sname = $dname.".".$used{$dname};
						$used{$dname}++;
					}
				}
			}
			$used{$sname}++;
			$used{$dname}++;
			
			
			push(@ret,">".$sname."  ".$fbdesc[$bb]."\n".join("",@{$resfas[$row]})."\n");
		}
		$row++;
	}
	
	return \@ret;
}




sub aa_merge_{
	my $aa_a = $_[0];
	my $aa_b = $_[1];
	
	
	if($#{$aa_a} != $#{$aa_b}){
		die;
	}
	
	my @ret;
	for(my $ii = 0;$ii <= $#{$aa_a};$ii++){
		my $af = 0;
		my $bf = 0;
		if(${$aa_a}[$ii] =~ /[A-Za-z]/){
			push(@ret,${$aa_a}[$ii]);
		}elsif(${$aa_b}[$ii] =~ /[A-Za-z]/){
			push(@ret,${$aa_b}[$ii]);
		}else{
			push(@ret,${$aa_a}[$ii]);
		}
		
	}
	return \@ret;
	
}





#配列の identity, coverage を返す
sub seq_ident_{
	my $aa_a = $_[0];
	my $aa_b = $_[1];
	
	if($#{$aa_a} != $#{$aa_b}){
		return 0,0;
	}
	my $ide = 0;
	my $count = 0;
	my $acou = 0;
	my $bcou = 0;
	for(my $ii = 0;$ii <= $#{$aa_a};$ii++){
		my $af = 0;
		my $bf = 0;
		if(${$aa_a}[$ii] =~ /[A-Za-z]/){
			$af = 1;
			$acou++;
		}
		if(${$aa_b}[$ii] =~ /[A-Za-z]/){
			$bf = 1;
			$bcou++;
		}
		if($af*$bf == 1){
			$count++;
			if(${$aa_a}[$ii] eq ${$aa_b}[$ii]){
				$ide++;
			}
		}
		
	}
	my $cov = $count/$acou;
	if($acou > $bcou){
		if($bcou == 0){
			return -1,-1;
		}
		$cov = $count/$bcou;
	}
	if($count == 0){
		return 0,0;
	}
	return $ide/$count,$cov;
}

#配列の identity, coverage を返す
sub seq_ident{
	my $seq_a = $_[0];
	my $seq_b = $_[1];
	
	
	$seq_a =~ s/[\s]//g;
	$seq_b =~ s/[\s]//g;
	
	
	if(length($seq_a) != length($seq_b)){
		return 0.0,0.0;
	}
	
	
	my @aa_a = split(//,$seq_a);
	my @aa_b = split(//,$seq_b);
	
	
	my $ide = 0;
	my $count = 0;
	my $acou = 0;
	my $bcou = 0;
	for(my $ii = 0;$ii <= $#aa_a;$ii++){
		my $af = 0;
		my $bf = 0;
		if($aa_a[$ii] =~ /[A-Za-z]/){
			$af = 1;
			$acou++;
		}
		if($aa_b[$ii] =~ /[A-Za-z]/){
			my $bf = 1;
			$bcou++;
		}
		if($af*$bf == 1){
			$count++;
			if($aa_a[$ii] eq $aa_b[$ii]){
				$ide++;
			}
		}
		
	}
	my $cov = $count/$acou;
	if($acou > $bcou){
		$cov = $count/$bcou;
	}
	return $ide,$cov;
}


sub waitUntil{
	sleep(60);
	while(1==1){
		sleep(100);
		my @dd = `date`;
		print $dd[0];
		my @line = `qstat|wc`;
		my $flag = 0;
		foreach my $ll(@line){
			if($ll =~ /^[\s]+([0-9])[\s]+/){
				my $cou = $1;
				if($cou == 0){
					$flag = 1;
				}
			}
			print $ll;
		}
		if($flag == 1){
			last;
		}
	}
	print "Let's move to the next.\n";
}





sub reduce_alignment{
	my @allfas = @{$_[0]};
	my $dummyfile = "dummy_dummy.".$$.rand().".fas";
	
	my %used;
	my %num_name_map;
	open(OUT,">".$dummyfile);
	for(my $ii = 0;$ii <= $#allfas;$ii++){
		my $ff = $allfas[$ii];
		my $name = fastaloader::get_name_fasta($ff);
		my $seq = fastaloader::get_seq_fasta($ff);
		my $aname = $name;
		my $seqcode = "seq.".$ii;
		$num_name_map{$seqcode} = $name;
		print OUT ">".$seqcode."\n";
		$seq =~ s/[^A-Za-z]//g;
		print OUT $seq."\n";
	}
	close(OUT);
	
	system($cdhit." -i $dummyfile -o $dummyfile.cdout1 -c 0.9");
	
	my $cou = 0;
	my %s_90;
	open(IN,"$dummyfile.cdout1");
	while(my $ss = <IN>){
		if($ss =~ />[\s]*([^\s]+)/){
			$s_90{$num_name_map{$1}} = 100;
			$cou++;
		}
	}
	close(IN);
	if($cou > 6000){
		unlink("$dummyfile.cdout1");
		unlink("$dummyfile.cdout1.clstr");
		system($cdhit." -i $dummyfile -o $dummyfile.cdout1 -c 0.8");
		%s_90 = ();
		open(IN,"$dummyfile.cdout1");
		while(my $ss = <IN>){
			if($ss =~ />[\s]*([^\s]+)/){
				$s_90{$num_name_map{$1}} = 100;
			}
		}
		close(IN);
	}
	
	
	unlink($dummyfile);
	unlink("$dummyfile.cdout1");
	unlink("$dummyfile.cdout1.clstr");
	my @ret;
	
	for(my $ii = 0;$ii <= $#allfas;$ii++){
		my $ff = $allfas[$ii];
		my $name = fastaloader::get_name_fasta($ff);
		if(defined $s_90{$name}){
			push(@ret,$ff);
		}
	}
	return \@ret;
	
}









sub fastaProcess{
	my $fa = $_[0];
	
	my @faname;
	my @fa_aa;
	my @fadesc;
	foreach my $ff(@{$fa}){
		my @pt = split(/[\r\n]/,$ff);
		my $name = "";
		my $desc = "";
		my $seq = "";
		for(my $ii = 0;$ii <= $#pt;$ii++){
			if($pt[$ii] =~ /^[\s]*>[\s]*([^\s]+)/){
				my $pname = $1;
				if(length($name) > 0){
					die;
					
				}
				$name = $pname;
				if(	$pt[$ii] =~ /^[\s]*>[\s]*[^\s]+[\s]+([^\r\n]+)/){
					my $pdesc = $1;
					$desc = $pdesc;
				}
				
			}else{
				$pt[$ii] =~ s/[\s]//g;
				if(length($name) > 0){
					$seq .= $pt[$ii];
				}else{
					if(length($pt[$ii]) != 0){
						print STDERR $pt[$ii]." was not processed.\n";
					}
				}
			}
		}
		my @aa = split(//,$seq);
		push(@fa_aa,\@aa);
		push(@faname,$name);
		push(@fadesc,$desc);
	}
	return \@faname,\@fadesc,\@fa_aa;
}
