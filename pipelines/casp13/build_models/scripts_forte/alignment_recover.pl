package alignment_recover;
use strict;
use warnings;

#require "c:\\dir_important\\fastaloader.pl";
#require "c:\\dir_important\\blastparser.pl";


#require "/home/tsukasa/casp_database/forte/scripts/fastaloader.pl";
require "/home/yamamori/work/casp13/scripts_forte/fastaloader.pl";
#require "/home/tsukasa/casp_database/forte/scripts/blastparser.pl";
require "/home/yamamori/work/casp13/scripts_forte/blastparser.pl";




#test();



sub test{
	my $targetname = "query";
	my $targetseq = "ZZZABCGHIJKLMNZZZ";
	
	my %hs1 = (
		"query" => "ABC--GHIJKL",
		"seq_1" => "---EFGHIJK-"
	);
	my %hs1b = (
		"query" => "ABC--GHIJKL",
		"seq_1" => "ABC--------"
	);
	my %hs2 = (
		"query" => "ABC--GHIJKLMN",
		"seq_2" => "---DEGHIJKL-N"
	);
	my @pair = (\%hs1,\%hs2,\%hs1b);
	my %dest = (
	"seq_1" =>  "ABC-EFGHIJK--",
	"seq_2" =>  "---DE-GHIJKLN"
	);
	
	my %aligned = %{add_one_seq_in_msa_($targetname,$targetseq,\@pair,\%dest)};
	foreach my $kk(keys %aligned){
		print $kk." ".$aligned{$kk}."\n";
	}
}


sub add_one_seq_in_msa_with_blastresult{
	my $queryname = $_[0];
	my $queryseq = $_[1];
	my $blastresult = $_[2];
	my $msafas =  $_[3];
	
	my $outfilename = $_[4];
	
	my $remove_runover = 0;#MSA は、NC 末にマッチが無い部分を出力しない	
	if(defined $_[5]){
		$remove_runover = $_[5];
	}
	my @fas = @{fastaloader::load_seq($msafas)};
	
	my @blast = @{blastparser::get_result_array_file($blastresult)};
	my @hash = @{blastparser::parse_blast_result($blast[0])};
	@blast = ();
	
	
	my @pair;
	foreach my $bb(@hash){
		my %bres = %{$bb};
		my %hs;
		$hs{$queryname} = $bres{"q_seq"};
		$hs{$bres{"name"}} = $bres{"s_seq"};
		push(@pair,\%hs);
	}
	my %msa;
	my %desc;
	my @names;
	foreach my $ff(@fas){
		my $name = fastaloader::get_name_fasta($ff);
		my $seq = fastaloader::get_seq_fasta($ff);
		my $desc = fastaloader::get_desc_fasta($ff);
		$msa{$name} = $seq;
		$desc{$name} = $desc;
		push(@names,$name);
	}
	@fas = ();
	@hash = ();
	
	my %msares = add_one_seq_in_msa($queryname,$queryseq,\@pair,\%msa,$remove_runover);
	if(defined $outfilename){
		open(OUT,">".$outfilename);
		print OUT ">".$queryname."\n".$msares{$queryname}."\n";
		foreach my $nn(@names){
			print OUT ">".$nn." ".$desc{$nn}."\n";
			print OUT $msares{$nn}."\n";
		}
		close(OUT);
	}else{
		my @ret;
		push(@ret,">".$queryname."\n".$msares{$queryname});
		foreach my $nn(@names){
			push(@ret,">".$nn." ".$desc{$nn}."\n".$msares{$nn});
		}
		return \@ret;
	}
}





#msa を壊さずに一本の配列を追加する

#並べたい配列の名前
#並べたいアミノ酸配列
#ペアワイズアラインメント、スコアの高い順#名前＝＞アミノ酸配列 というハッシュエントリが2つずつ入っている
#貼り付けたい MSA #配列の名前=>アミノ酸配列
#の順で変数を渡すと、配列の名前=>並べられたアミノ酸配列
#というハッシュが返ってくる
sub add_one_seq_in_msa_{
	my $targetname = $_[0];#並べたい配列の名前
	
	my $targetseq = uc $_[1];#並べたいアミノ酸配列
	
	my @pair = @{$_[2]};#ペアワイズアラインメント、スコアの高い順
	                    #名前＝＞アミノ酸配列 というハッシュエントリが2つずつ入っている
						#"indexoffset>"＝＞開始点 が、同一アミノ酸配列を持つドメインが複数ある場合のためにある
	
	my %dest_alignment = %{$_[3]};#貼り付けたい MSA 
	                              #配列の名前=>アミノ酸配列
	
	
	my $remove_runover = 0;#MSA は、NC 末にマッチが無い部分を出力しない	
	if(defined $_[4]){
		$remove_runover = $_[4];
	}
	my $dummyrun = 0;#dummy アラインメントが出来た時に削除しない
	if(defined $_[5]){
		$dummyrun = $_[5];
	}
	
	
	
	$targetseq =~ s/[^A-Z]//g;
	my @aas = split(//,$targetseq);
	
	my @seq_ali_map;#MSA 上で何番目のカラムに来るか
	foreach my $aa(@aas){
		push(@seq_ali_map,-1);
	}
	
	
	
	my $dummyname = "dummy_dummy_dummy".rand().$$;
	my %mappedpos;#既にマップされたカラム
	foreach my $hs(@pair){
		my $subseq = "";
		my $subname = "";
		my $quseq = "";
		foreach my $kk(keys %{$hs}){
			if($targetname eq $kk){
				$quseq = uc ${$hs}{$kk};
			}else{
				$subseq = uc ${$hs}{$kk};
				$subname = $kk;
			}
		}
		
		if(length($subname) == 0 && length($quseq) == 0){
			next;
		}
		if(length($subname) == 0){
			print STDERR "couldn't find subjct seq in pairwise alignment.\n";
			die;
		}
		if(length($quseq) == 0){
			print STDERR "couldn't find target seq in pairwise alignment.\n";
			die;
		}
		
		if(!defined $dest_alignment{$subname}){
			print STDERR "couldn't find $subname in MSA. skipped\n";
			next;
		}
		
		my $tmpm = uc $dest_alignment{$subname};
		$tmpm =~ s/[^A-Z]//g;
		
		my $tmps = $subseq;
		$tmps =~ s/[^A-Z]//g;
		
		
		my $tmpqs = uc $quseq;
		$tmpqs =~ s/[^A-Z]//g;
		
		
		my $offsets = index($tmpm,$tmps);
		
		if($offsets == -1){
			print STDERR "couldn't find $subname. $tmpm\n$tmps\n\ntrying mafft.\n";
			#die;
			#next;
			#print STDERR "different sequence ware found $targetname. $targetseq\n$tmpqs \ndummy sequence will be attached. \n";
			
			
			my $mafftpath = "mafft";
			my $mafftout  = "tmpmout".$$.rand();
			my $mafftin  = "tmpmin".$$.rand();
			open(PPOUT,">".$mafftin);
			print PPOUT ">$dummyname  \n".$tmps."\n";
			print PPOUT ">$subname  \n".$tmpm."\n";
			close(PPOUT);
			
			system($mafftpath." --auto $mafftin > $mafftout 2> /dev/null");
			my %ppfas = %{fastaloader::load_seq_hash($mafftout)};
			my $dseq = fastaloader::get_seq_fasta($ppfas{"$dummyname"});
			my $zseq = fastaloader::get_seq_fasta($ppfas{"$subname"});
			
			
			
			
			$dseq =~ s/[\s]//g;
			$zseq =~ s/[\s]//g;
			
			my @chk1 = split(//,$dseq);
			my @chk2 = split(//,$zseq);
			my $diffcou = 0;
			for(my $l = 0;$l <= $#chk1;$l++){
				if($chk2[$l] =~ /[A-Z]/ && $chk1[$l] =~ /[A-Z]/){
					if($chk2[$l] ne $chk1[$l]){
						$diffcou++;
					}
					
				}
			}
			if($diffcou > 3){
			
				print STDERR "the sequence is extremely different\n$dseq.\n$zseq\n$diffcou\n";
				#die;
			}
			
			
			unlink($mafftin);
			unlink($mafftout);
			my @pll;
			my %pllpair;
			$pllpair{$dummyname} = $dseq;
			$pllpair{$subname} = $zseq;
			push(@pll,\%pllpair);
			
			
			
			my %ppres = %{add_one_seq_in_msa_($dummyname,$tmps,\@pll,\%dest_alignment,0,1)};
			%dest_alignment = %ppres;
			$subname = $dummyname;
			$subseq = $tmps;
			$quseq = $tmpqs;
			#$offsets = 0;
			#next;
		}
		
		
		my $offsetq = 0;
		if(defined ${$hs}{"indexoffset>"} && ${$hs}{"indexoffset>"} > 1){
			$offsetq = index($targetseq,$tmpqs,${$hs}{"indexoffset>"}-1);
		}else{
			$offsetq = index($targetseq,$tmpqs);
		}
		if($offsetq < 0){
		
			print STDERR "different sequence ware found $targetname. $targetseq\n$tmpqs \ndummy sequence will be attached. \n";
			next;
		
		}
		$subseq =~ s/[\s]//g;
		$quseq =~ s/[\s]//g;
		
		
		my @ss = split(//,$subseq);
		
		
		my $ssmsa = $dest_alignment{$subname};
		$ssmsa =~ s/[\s]//g;
		my @ssmsa = split(//,$ssmsa);
		
		
		my @sspair = split(//,$subseq);
		my %submap;#ss 上の、どの残基がどのカラムに来るか。
		
		
		my $soff = 0;
		my $mscou = 0;
		for(my $ii = 0;$ii <= $#ssmsa;$ii++){
			if($ssmsa[$ii]  =~ /[A-Za-z]/){
				my $m = $&;
				if($mscou >= $offsets){
					while($soff <= $#sspair){
						if($sspair[$soff] =~ /[a-zA-Z]/){
							my $a = $&;
							if(uc $m ne uc $a){
								if($a !~ /[Xx]/){#blast がマスクしているとき
									print "sequence mismatch?? error in code $m $a\n$subseq\n$ssmsa";
									die;
								}
							}
							$submap{$soff} = $ii;#pairwisealignment 上で $soff に来る残基は $ii 番目のカラムに来る。
							$soff++;
							last;
						}
						$soff++;
					}
				}
				$mscou++;
			}
		}
		
		
		my @qq = split(//,$quseq);
		my $qcou = 0;
		for(my $ii = 0;$ii <= $#qq;$ii++){
			if($qq[$ii] =~ /[A-Za-z]/){
				if(defined $submap{$ii}){
					if($seq_ali_map[$qcou+$offsetq] == -1){
						my $pos = $qcou+$offsetq;
						my $cpos = $submap{$ii};
						if(!defined $mappedpos{$cpos}){#他の残基が既に MSA 上のそのカラムにある場合
							
							$mappedpos{$cpos} = 100;
							
							for(my $jj = $pos-1;$jj >= 0;$jj--){
								if($seq_ali_map[$jj]  > $cpos){#前に並べられたアラインメントと齟齬が無いか調べる
									$cpos = -1;
									last;
								}elsif($seq_ali_map[$jj] > -1){
									last;
								}
							}
							for(my $jj = $pos+1;$jj <= $#seq_ali_map;$jj++){
								if($seq_ali_map[$jj] > -1 && $seq_ali_map[$jj]  < $cpos){#前に並べられたアラインメントと齟齬が無いか調べる
									$cpos = -1;
									last;
								}elsif($seq_ali_map[$jj] > -1){
									last;
								}
							}
							if($cpos != -1){
								$seq_ali_map[$pos] = $submap{$ii};
							}
						}
					}
				}
				if(uc $qq[$ii] ne uc $aas[$qcou+$offsetq]){
					if($qq[$ii] !~ /[Xx]/){
						print "sequence mismatch?? error in code $qq[$ii] $aas[$qcou+$offsetq]\n$quseq\n$targetseq";
						die;
					}
				}
				$qcou++;
			}
		}
		
		
		my $flag = 1;
		foreach my $ss(@seq_ali_map){
			if($ss == -1){
				$flag = -1;
				last;
			}
		}
		if($flag == 1){
			last;
		}
	}
	my %insertion;
	#MSA上にギャップを入れる場合。ギャップ場所=>そこに入るアミノ酸
	my %alignedpos_aa;
	my $gcou = "";
	my $headremove = -2;#前後のミスマッチは削除する場合
	my $tailremove = -1;#前後のミスマッチは削除する場合
	for(my $ii = 0;$ii <= $#seq_ali_map;$ii++){
		if($seq_ali_map[$ii] != -1){
			if($headremove == -2){
				$headremove = $seq_ali_map[$ii]-1;
			}
			$tailremove = $seq_ali_map[$ii]+1;
			$alignedpos_aa{$seq_ali_map[$ii]} = $gcou.$aas[$ii];
			if(length($gcou) > 0){
				$gcou =~ s/./-/g;
				$insertion{$seq_ali_map[$ii]} = $gcou;
			}
			
			$gcou = "";
		}else{
			$gcou .= $aas[$ii];
		}
	}
	my $tail = $gcou;#末尾にギャップを入れる場合
	my $tailgap = $gcou;
	$tailgap =~ s/./-/g;
	my %restmp;
	my $destlen = 0;
	my $destref = "";
	foreach my $kk(keys %dest_alignment){
		if($destlen == 0){
			$destlen = length($dest_alignment{$kk});
			$destref = $dest_alignment{$kk};
		}elsif($destlen != length($dest_alignment{$kk})){
			print STDERR "MSA length was different please check.\n$targetname\n\n$dummyrun\n$destref\n".$dest_alignment{$kk}."\n";
			die;
		}
		my $ps = $dest_alignment{$kk};
		$ps =~ s/[\s]//g;
		my @paa = split(//,$ps);
		
		for(my $ii = 0;$ii < $destlen;$ii++){
			if($remove_runover == 1 &&($headremove >= $ii || $tailremove <= $ii)){
				$paa[$ii] = "";
			}
			if(defined $insertion{$ii}){
				$paa[$ii] = $insertion{$ii}.$paa[$ii];
			}
		}
		if(length($tail) > 0){
			$paa[$destlen] = $tailgap;
		}
		$restmp{$kk} = \@paa;
		
	}
	my @que_aligned;
	for(my $ii = 0;$ii < $destlen;$ii++){
		if($remove_runover == 1 &&($headremove >= $ii || $tailremove <= $ii)){
			next;
		}
		if(defined $alignedpos_aa{$ii}){
			push(@que_aligned,$alignedpos_aa{$ii});
		}else{
			push(@que_aligned,"-");
		}
	}
	if(length($tail) > 0){
		push(@que_aligned,$tail);
	}
	$restmp{$targetname} = \@que_aligned;
	
	my %res;
	foreach my $kk(keys %restmp){
		if($kk eq $dummyname && $dummyrun==0){
			next;
		}
		if(!@{$restmp{$kk}}){
			print STDERR $kk." ignored.\n";
			print STDERR join("",@{$restmp{$kk}});
			next;
		}
		
		$res{$kk} = join("",@{$restmp{$kk}});
	}
	return \%res;
}






1;
