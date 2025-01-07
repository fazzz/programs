package blastparser;

use strict;
use warnings;

#test();
sub test{
#	require "c:\\dir_important\\fastaloader.pl";
#	my @mres = @{blastparser::get_multi_result_array("test.blastout")};

#	foreach my $mm(@mres){
#		if($mm =~ /Query[\s]*=[\s]*([^\s]+)/){
#			print $1."\n";
#			
#		}
#		my @res = @{blastparser::get_result_array($mm)};
#		foreach my $rr(@res){
#			
#			foreach my $pres(@{blastparser::parse_blast_result($rr)}){
#				my %cr = %{$pres};
#				print $cr{"name"}."\t".$cr{"Expect"}."\n";
#			}
#		}
#		blastalign("test.blastout","test.fas",0.0001);
#	}
	
	
	
	
	my @qarray = (
	"ABCDEFGHIJKL",
	"BCDEFGHIJKL",
	"CDEFGHIJKL",
	"ABCDEFGHIJ",
	
	"---CDEFGHIJK",
	"--BC--DEFGHIJK",
	"ABCDEFGHIJ--K---"
	);
	
	my @sarray = (
	"ABCDEFGHIJ--",
	"BCD--GHIJKL",
	"CD-FGHI-KL",
	"A-CDEFGHIJ",
	
	"XYZCDEFGHIJK",
	"XYBCZZDEFGHIJK",
	"ABCDEFGHIJZZKZZZ"
	);
	my @starts =(
	1,
	2,
	3,
	1,
	3,
	2,
	1
	);
	
	my @res = @{max_redundant_alignment2(\@qarray,\@starts,\@sarray)};
	foreach my $rr(@res){
		print $rr."\n";
	}
	

}


#blastalign('C:\dummy\work\request\toxoplasma\ip3\TgME49-233220_477_766.blastout','C:\dummy\work\request\toxoplasma\ip3\TgME49-233220_477_766.blast_aligned',0.001);
#require "a3mprocess.pl";
#a3mprocess::fas_to_a3m('C:\dummy\work\request\toxoplasma\ip3\TgME49-233220_477_766.blast_aligned','C:\dummy\work\request\toxoplasma\ip3\TgME49-233220_477_766.blast_a3m');



sub blastalign{#blast の結果ファイルからアラインメントを生成する
	my $infilename = $_[0];
	my $outfilename = $_[1];
	my $threshold = $_[2];
	my $sourceseq = $_[3];#もし全長が欲しい場合
	my @mres = @{get_multi_result_array($infilename)};
	my $bbcou = 0;
	for(my $m = 0;$m <= $#mres;$m++){
		
		if($bbcou == 0){
			open(BBOUT,">".$outfilename);
		}else{
			open(BBOUT,">".$outfilename."-".$bbcou);
		
		}
		$bbcou++;
		if($mres[$m] =~ /Query[\s]*=[\s]*([^\s]+)/){#これが無い場合がある。。。
			my $queryname = $1;
			my @res = @{get_result_array($mres[$m])};
			my @hashres;
			my @names;
			my @descs;
			if(defined $sourceseq){
				my %bhash;
				$bhash{"name"} = $queryname;
				$bhash{"Expect"} = 0.0;
				$bhash{"Identities"} = "100/100";
				$bhash{"Positives"} = "100/100";
				
				
				$bhash{"score"} = 10000000;
				$bhash{"s_seq"} = $sourceseq;
				$bhash{"q_seq"} = $sourceseq;
				$bhash{"s_start"} = 1;
				$bhash{"s_end"} = length($sourceseq);
				$bhash{"q_start"} = 1;
				$bhash{"q_end"} = length($sourceseq);
				
				
				
				push(@names,$bhash{"name"});
				push(@descs,"Expect=".$bhash{"Expect"}.", Identities=".$bhash{"Identities"}.",Positives=".$bhash{"Positives"}.", ");
				push(@hashres,\%bhash);
				
				
				
			}
			
			my $cou = 0;
			foreach my $rr(@res){
				my @bres = @{parse_blast_result($rr)};
				foreach my $bb(@bres){
					my %bhash = %{$bb};
					if($bhash{"Expect"} < $threshold){
						
						#my %khash = %{pairalign_postprocess(\%bhash)};#しない方が良い
						#push(@names,"query");
						#push(@names,$bhash{"name"});
						#push(@hashres,\%khash);
						
						
						
						if(defined $sourceseq){
							my $cch = uc $bhash{"q_seq"};
							$cch  =~ s/[^A-Z]//g;
							if(uc $cch eq uc $sourceseq){
								next;
							}
						}
						
						push(@names,$bhash{"name"});
						push(@descs,"Expect=".$bhash{"Expect"}.", Identities=".$bhash{"Identities"}.",Positives=".$bhash{"Positives"}.", ");
						push(@hashres,\%bhash);
						
						
						
					}
				}
				if($cou++ > 100){
				#	last;
				}
			}
			
			
			if($#hashres == -1 || ($#hashres == 0 && defined $sourceseq)){
				if(defined $sourceseq){
					print BBOUT ">".$queryname."\n$sourceseq\n";
				}
			}else{
				my @msa = @{result_to_msa(\@hashres)};
				@hashres = ();
				$mres[$m] = "";
				
				my %nameused;
				
				print BBOUT ">".$queryname." query\n".$msa[0]."\n";
				
				$nameused{$queryname} = 2;
				unshift(@names,"");
				unshift(@descs,"");
				for(my $ii = 1;$ii <= $#msa;$ii++){
					if(defined $nameused{$names[$ii]}){
						print BBOUT ">".$names[$ii]."_".$nameused{$names[$ii]}." ".$descs[$ii]."\n".$msa[$ii]."\n";
						$nameused{$names[$ii]}++;
					}else{
						print BBOUT ">".$names[$ii]." ".$descs[$ii]."\n".$msa[$ii]."\n";
						$nameused{$names[$ii]} = 2;
					}
				}
			}
		}
		close(BBOUT);
	}
}	




sub pairalign_postprocess{#ローカルアラインメントの質を高めるために何らかの処理をする
	my %hash = %{$_[0]};
	my @qaa = split(//,$hash{"q_seq"});
	my @saa = split(//,$hash{"s_seq"});
	my @qres;
	my @sres;
	
	my $flag = 0;
	my @bstart;
	my @bend;
	for(my $ii = 0;$ii <= $#qaa;$ii++){
		if($qaa[$ii] eq "-" || $saa[$ii] eq "-"){
			if($flag < 0){
				for(my $jj = $ii-1;$jj >= 0;$jj--){
					if($qaa[$jj] eq $saa[$jj]){
						push(@bstart,$jj+1);
						push(@bend,$ii-1);
						last;
					}
				}
			}
			$flag = $ii;
		}else{
			if($flag > -1){
				my $gstart = $ii;
				for(my $jj = $ii;$jj <= $#qaa;$jj++){
					if($qaa[$jj] eq $saa[$jj]){
						if($gstart != $jj){
							push(@bstart,$gstart);
							push(@bend,$jj-1);
						}
						last;
					}
				}
				$flag = -1;
			}
		}
	}
	
	for(my $jj = $#qaa;$jj >= 0;$jj--){
		if($qaa[$jj] eq $saa[$jj]){
			push(@bstart,$jj+1);
			push(@bend, $#qaa);
			last;
		}
	}
	for(my $bb = 0;$bb <= $#bstart;$bb++){
		if($bstart[$bb] > $bend[$bb]){
		
		}else{
			my @ext;
			for(my $jj = $bstart[$bb];$jj <= $bend[$bb];$jj++){
				 push(@ext,"-");
			}
			my $ps = join("",@ext);
			$qaa[$bend[$bb]] .= $ps;
			$saa[$bstart[$bb]] = $ps.$saa[$bstart[$bb]];
			
		}
	}
	my @qs = split(//,join("",@qaa));
	my @ss = split(//,join("",@saa));
	my @sr;
	my @qr;
	
	for(my $ii = 0;$ii <= $#qs;$ii++){
		if($qs[$ii] eq "-" && $ss[$ii] eq "-"){
			
		}else{
			push(@sr,$ss[$ii]);
			push(@qr,$qs[$ii]);
		}
	}
	$hash{"q_seq"} = join("",@qr);
	$hash{"s_seq"} = join("",@sr);
	return \%hash;
	
}

sub get_forte_result_array{
	my $fname = $_[0];
	my $buff = "";
	open(B_IN,$fname)or die;
	my @ret;
	my $name = "";
	my $flag = 0;
	while(my $ss = <B_IN>){
		$buff .= $ss;
		if($ss =~ /^[\s]*[0-9]+[\s]+([^\s][^\r\n]+)$/){
			my $li = $1;
			if(length($buff) > 0){
				push(@ret,$buff);
			}
			$buff ="";
		}
	}
	close(B_IN);
	return \@ret;
}


sub parse_forte_result{
	my $region =$_[0];
	my @lines = split(/[\r\n]+/,$region);
	my %current;
	my $query="";
	my $target="";
	my @ret;
	my $qline = 0;
	for(my $ii = 0;$ii <= $#lines;){
		$lines[$ii] =~ s/[\r\n]//g;
		my $li = $lines[$ii];
		if($li =~ /^[\s]*[0-9]+[\s]+([^\s][^\r\n]+)$/){
			if($li =~ /([^\s]+)[\s]+([^\s]+)$/){
				my $name = $1;
				my $score = $2;
				$current{"name"} = $name;
				$current{"score"} = $score;
				$current{"s_seq"} = $target;
				$current{"q_seq"} = $query;
				$current{"s_start"} = 1;
				$current{"s_end"} = length($target);
				$current{"q_start"} = 1;
				$current{"q_end"} = length($target);
				my %tmp = %current;
				push(@ret,\%tmp);
				%current = ();
				$query = "";
				$target = "";
			}
			$ii++;
		}else{
			if($li =~ /^[^\s]+[\s]+([A-Za-z\-]+)$/){
				$query .= $1;
				$ii++;
				if($lines[$ii] =~ /^concor/){
					$ii++;
				}
				
				$lines[$ii] =~ s/[\r\n]//g;
				$li = $lines[$ii];
				if($li =~ /^[^\s]+[\s]+([A-Za-z\-]+)$/){
					$target .= $1;
				}else{
					print "formaterror \n".$region."\n";
					die;
				}
			}
			$ii++;
		}
	}
	
	return \@ret;
}

sub result_to_msa{
	my @res = @{$_[0]};#parse_blast_resut で作成した Hit のハッシュ。
	
	my @qarray;
	my @sarray;
	my @qstart;
	foreach my $rr(@res){
		my %hash = %{$rr};
		push(@qarray,$hash{"q_seq"});
		push(@sarray,$hash{"s_seq"});
		push(@qstart,$hash{"q_start"});
	}
	
	#return max_redundant_alignment(\@qarray,\@qstart,\@sarray);
	return max_redundant_alignment2(\@qarray,\@qstart,\@sarray);
}


sub max_redundant_alignment2{
	my @qarray = @{$_[0]};#ペアワイズでの query の配列
	my @startarray = @{$_[1]};#query の開始点
	my @sarray = @{$_[2]};#subject の配列
	my %gaps;
	my @ret;
	if(!defined $startarray[0]){
		print $sarray[0]."\n";
		die;
	}
	my $startmin = $startarray[0]-1;
	my $endmax = 0;
	my @que_fulllength;
	
	
	for(my $ii = 0;$ii <= $#qarray;$ii++){
		my $que = $qarray[$ii];
		$que =~ s/[\s]//g;
		my @q = split(//,$que);
		my $cou = $startarray[$ii]-1;
		my $tmpstart = $cou; 
		if($cou < $startmin){
			$startmin = $cou;
		}
		$que =~ s/[^A-Za-z]//g;
		my $end = length($que)+$cou;
		if($end > $endmax){
			$endmax = $end;
		}
	}
	my @maxseq;
	for(my $ii = 0;$ii <= $#sarray;$ii++){
		my @tmp;
		for(my $jj = 0;$jj < $endmax;$jj++){
			push(@tmp,"");#偶数番目にはギャップが入り、奇数番目にはアラインされた残基が入る
			push(@tmp,"");
		}
		push(@tmp,"");
		push(@maxseq,\@tmp);
	}
	
	
	my @queryaa;
	for(my $jj = 0;$jj < $endmax;$jj++){
		push(@queryaa,"");
		push(@queryaa,"");
	}
	push(@queryaa,"");
	push(@maxseq,\@queryaa);
	for(my $ii = 0;$ii <= $#qarray;$ii++){
		my $que = $qarray[$ii];
		$que =~ s/[\s]//g;
		my @q = split(//,$que);
		my $sbb = $sarray[$ii];
		$sbb =~ s/[\s]//g;
		my @s = split(//,$sbb);
		my $start = -1;
		my $end = -1;
		my $buff = "";
		my $cou = $startarray[$ii]-1;
		my $target_ali = "";
		my $que_ali = "";
		
		
		my $lastpos = $startarray[$ii]-1;
		for(my $qi = 0;$qi <= $#q;$qi++){
			if($q[$qi] =~ /[A-Za-z]/){
				$lastpos++;
			}
		}
		
		for(my $qi = 0;$qi <= $#q;$qi++){
			if($q[$qi] =~ /[A-Za-z]/){
				$cou++;
				${$maxseq[$ii]}[$cou*2-1] = $s[$qi];
				if($queryaa[$cou*2-1] ne "" && $queryaa[$cou*2-1] ne $q[$qi]){
					
					die $queryaa[$cou*2-1]."\t".$q[$qi];
				}
				$queryaa[$cou*2-1] = $q[$qi];
			}else{
				if($cou == $startarray[$ii]-1){
					${$maxseq[$ii]}[0] .= $s[$qi]; 
				}elsif($cou == $lastpos){
					${$maxseq[$ii]}[$endmax*2] .= $s[$qi]; 
				}else{
					${$maxseq[$ii]}[$cou*2] .= $s[$qi]; 
				}
			}
		}
		
	}
	
	my @longest;
	for(my $ii = 0;$ii <= $endmax*2;$ii++){
		my $ll = 0;
		for(my $jj = 0;$jj <= $#maxseq;$jj++){
			if(length(${$maxseq[$jj]}[$ii]) > $ll){
				$ll = length(${$maxseq[$jj]}[$ii])
			}
			
		}
		
		for(my $jj = 0;$jj <= $#maxseq;$jj++){
			while($ll > length(${$maxseq[$jj]}[$ii])){
				${$maxseq[$jj]}[$ii] .= "-";
			}
		}
	}
	
	
	push(@ret,join("",@{$maxseq[$#maxseq]}));
	for(my $jj = 0;$jj <= $#maxseq-1;$jj++){
		push(@ret,join("",@{$maxseq[$jj]}));
	}
	return \@ret;
	
}
sub max_redundant_alignment{
	my @qarray = @{$_[0]};#ペアワイズでの query の配列
	my @startarray = @{$_[1]};#query の開始点
	my @sarray = @{$_[2]};#subject の配列
	my %gaps;
	my @ret;
	my $startmin = $startarray[0]-1;
	my $endmax = 0;
	for(my $ii = 0;$ii <= $#qarray;$ii++){
		my $que = $qarray[$ii];
		$que =~ s/[\s]//g;
		my @q = split(//,$que);
		my $start = -1;
		my $end = -1;
		my $buff = "";
		my $cou = $startarray[$ii]-1;
		my $tmpstart = $cou; 
		
		for(my $jj = 0;$jj <= $#q;$jj++){
			if($q[$jj] =~ /[\-\.]/){
				$buff .= $q[$jj];
				if($cou  == $startarray[$ii]-1){
					$tmpstart--;
				}
			}else{
				if(length($buff) > 0){
					if(!defined $gaps{$cou}){
						$gaps{$cou} = "";
					}
					if(length($gaps{$cou}) < length($buff)){
						$gaps{$cou} = $buff;
					}
					$buff = "";
				}
				$cou++;
			}
		}
		
		if($tmpstart < $startmin){
			$startmin = $tmpstart;
		}
		if(length($buff) > 0){
			if(!defined $gaps{$cou}){
				$gaps{$cou} = "";
			}
			if(length($gaps{$cou}) < length($buff)){
				$gaps{$cou} = $buff;
			}
			$buff = "";
		}
	}
	
	
	for(my $ii = 0;$ii <= $#qarray;$ii++){
		my $que = $qarray[$ii];
		$que =~ s/[\s]//g;
		my @q = split(//,$que);
		my $sbb = $sarray[$ii];
		$sbb =~ s/[\s]//g;
		my @s = split(//,$sbb);
		my $start = -1;
		my $end = -1;
		my $buff = "";
		my $cou = $startarray[$ii]-1;
		my $target_ali = "";
		my $que_ali = "";
		for(my $jj = 0;$jj <= $#q;$jj++){
			if($q[$jj] =~ /[\-\.]/){
				$buff .= $q[$jj];
				
			}else{
				if($cou  == $startarray[$ii]-1){
					my $pst = "";
					for(my $kk = $startmin;$kk < $cou-length($target_ali);$kk++){
						$pst .= "-";
						if(defined $gaps{$kk}){
							$pst .= $gaps{$kk};
							
						}
					}
					$target_ali = $pst.$target_ali;
					$que_ali = $pst.$que_ali;
				}
			
				if(defined $gaps{$cou}){
					if(length($gaps{$cou}) > length($buff)){
						for(my $kk = length($buff);$kk < length($gaps{$cou});$kk++){
							$target_ali .= "-";
							$que_ali .= "-";
						}
					}
				}
				
				$buff = "";
				$cou++;
			}
			if(!defined $s[$jj]){
				print join("",@qarray)."\n";
				print join("",@sarray)."\n";
				die;
			}
			$target_ali .= $s[$jj];
			$que_ali .= $q[$jj];
		}
		
		if(defined $gaps{$cou}){
			if(length($gaps{$cou}) > length($buff)){
				for(my $kk = length($buff);$kk < length($gaps{$cou});$kk++){
					$target_ali .= "-";
					$que_ali .= "-";
				}
			}
		}
		push(@ret,$que_ali);
		push(@ret,$target_ali);
	}
	my $mlen = length($ret[0]);
	
	foreach my $rr(@ret){
		if($mlen < length($rr)){
			$mlen = length($rr);
		}
	}
	
	
	for(my $ii = 0;$ii <= $#ret;$ii++){
		while($mlen > length($ret[$ii])){
			$ret[$ii] .= "-";
		}
	}
	
	
	
	return \@ret;
	
}

sub map_seq{
	my $seq = $_[0];
	my $target = $_[1];
	my $target_st = $_[2];
	my $target_en = $_[3];
	$seq =~ s/[\s]//g;
	$target =~ s/[\s]//g;
	my @ss = split(//,$seq);
	my @tar = split(//,$target);
	
	my @part;
	for(my $ii = $target_st;$ii <= $target_en;$ii++){
		push(@part,$ss[$ii-1]);
	}
	my @ret;
	my $cou = 0;
	foreach my $tt(@tar){
		if($tt =~ /[A-Za-z]/){
			if(!defined $part[$cou]){
				print $seq."\n".$target_st."\t".$target_en."\n".$target."\n";
			}else{
				push(@ret,$part[$cou]);
			}
			$cou++;
		}else{
			push(@ret,"-");
		}
	}
	return join("",@ret);
}

sub get_multi_result_array{
	my $fname = $_[0];
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
				push(@ret,$buff);
				$buff = "";
			}
			
			if($ss =~ /^Query=[\s]+$/){
				$ss =~ s/^Query=/Query=dummylabel/;
			}
		}
		$buff .= $ss;
		
		
	}
	if($buff =~ /Query=/){
		push(@ret,$buff);
	}
	close(B_IN);
	if(!defined $ret[0]){
		print $fname." is not a blast result file?\n";
		die;
	}
	return \@ret;


}


sub get_result_array_file{
	my $fname = $_[0];
	my $buff = "";
	open(B_IN,$fname)or die;
	my @ret;
	my $name = "";
	my $flag = 0;
	while(my $ss = <B_IN>){
		if($ss =~ /^Query=[\s]*([^\s]+)[\s]*/){
			last;
		}
	}
	while(my $ss = <B_IN>){
		if($ss =~ /^>[\s]*([^\s]+)/){
			my $n = $1;
			
			if(length($buff) > 0 && length($name) > 0){
				push(@ret,$buff);
			}
			$buff = "";
			$name = $n;
		}
		
		if($ss =~ /^Lambda      K  /){
			last;
		}
		$buff .= $ss;
		
	}
	if(length($buff) > 0 && length($name) > 0){
		push(@ret,$buff);
		$buff = "";
	}
	close(B_IN);
	return \@ret;
}


sub get_result_array{
	my $fbuff = $_[0];
	my $buff = "";
	my @ret;
	my $name = "";
	my $flag = 0;
	
	$fbuff =~ s/\r\n/\n/g;
	$fbuff =~ s/\r/\n/g;
	my @pt = split(/\n/,$fbuff);
	foreach my $ss (@pt){
		if($ss =~ /^>[\s]*([^\s]+)/){
			my $n = $1;
			
			if(length($buff) > 0 && length($name) > 0){
				push(@ret,$buff);
			}
			$buff = "";
			$name = $n;
		}
		
		if($ss =~ /^Lambda      K  /){
			last;
		}
		$buff .= $ss."\n";;
		
	}
	if(length($buff) > 0 && length($name) > 0){
		push(@ret,$buff);
		$buff = "";
	}
	close(B_IN);
	return \@ret;
}





sub resulthash_to_fasta{
	my %hash = %{$_[0]};
	my @ret;
	push(@ret,">query $hash{'q_start'}-$hash{'q_end'}\n".$hash{"q_seq"});
	push(@ret,">".$hash{"name"}." $hash{'s_start'}-$hash{'s_end'}\n".$hash{"s_seq"});
	
	
	
	return \@ret;
	
	
}
#
# Score =   188 bits (479),  Expect = 4e-060, Method: Composition-based stats.
# Identities = 161/164 (98%), Positives = 161/164 (98%), Gaps = 0/164 (0%)
#
sub parse_blast_result{
	my $region =$_[0];
	$region =~ s/\r\n/\n/g;
	$region =~ s/\r/\n/g;
	my @lines = split(/\n/,$region);
	my %current;
	my @ret;
	my $qline = 0;
	for(my $ii = 0;$ii <= $#lines;){
		$lines[$ii] =~ s/[\r\n]//g;
		if($lines[$ii] =~ /^>[\s]*([^\s]+)/){
			my $hitname = $1;
			my $hitdesc = "";#一行しか取れない
			if($lines[$ii] =~ /^>[\s]*[^\s]+[\s]+([^\r\n]+)/){
				$hitdesc = $1;
			}
			for(;$ii <= $#lines;$ii++){
				if($lines[$ii] =~  /^Length[\s]*=/){
					
					if($lines[$ii+1] =~ /^[\s]*$/){
						last;
					}
				}
			}
			if($#lines <= $ii){
				last;
			}
			while(1==1){
				$current{"name"} = $hitname;;
				$current{"desc"} = $hitdesc;
				$current{"full"} =$region;
				if(!defined $hitname){
					die;
				}
				for(;$ii <= $#lines;$ii++){
					if($lines[$ii] =~ /^ Score =/){
						$ii--;
						last;
					}
				}
				
				
				if($#lines <= $ii){
					last;
				}
				
				for(my $jj = $ii+1;$jj <= $#lines;$jj++){
					if($lines[$jj] =~ /^Query /){
						$qline = $jj;
						last;
					}
					my @pt = split(/,/,$lines[$jj]);
					foreach my $pp(@pt){
						if($pp =~ /[\s]*([^\s].+)[\s]*=[\s]*(.+)/){
							my $lab = $1;
							my $val = $2;
							$lab =~ s/[\s]//g;
							if(defined $current{$lab}){
								print $lab." is already exists. ".$current{$lab}." ".$val."\n";
							}else{
								$current{$lab} = $val;
							}
							if($val =~ /\(([^\)]+)\)/){
								$current{$lab."%"} = $1;
							}
						}elsif($pp =~ /[\s]*([^\s].+)[\s]*\:[\s]*(.+)/){
							my $lab = $1;
							my $val = $2;
							$lab =~ s/[\s]//g;
							$current{$lab} = $val;
						}
					}
				}
				
				
				$current{"s_seq"} = "";
				$current{"q_seq"} = "";
				$current{"s_start"} = -1;
				$current{"s_end"} = -1;
				$current{"q_start"} = -1;
				$current{"q_end"} = -1;
				my $jj;
				my $multiflag = 0;
				for($jj = $qline;$jj <= $#lines;$jj++){
					if($lines[$jj] =~ /^ Score /){
						$multiflag = 1;
						last;
					}
					if($lines[$jj] =~ /^>/){
						$multiflag = 0;
						last;
					}
					$lines[$jj] =~ s/[\r\n]//g;
					if($lines[$jj] =~ /^Query[\s]+([0-9]+)[\s]*([^0-9\s]+)[\s]+([0-9]+)/){
						my $st = $1;
						my $en = $3;
						my $seq = $2;
						$current{"q_seq"} .= $seq;
						if($current{"q_start"} < 0){
							$current{"q_start"} = $st;
							
						}
						if($current{"q_end"} > $en || $current{"q_start"} > $st){
							print "format error \n".$region;
							die;
						}
						$current{"q_end"} = $en;
						
					}elsif($lines[$jj] =~ /^Query[\s]+(-+)[\s]+/){
						$current{"q_seq"} .= $1;
					
					}

					
					if($lines[$jj] =~ /^Sbjct[\s]+([0-9]+)[\s]*([^0-9\s]+)[\s]+([0-9]+)/){
						my $st = $1;
						my $en = $3;
						my $seq = $2;
						$current{"s_seq"} .= $seq;
						if($current{"s_start"} < 0){
							$current{"s_start"} = $st;
							
						}
						if($current{"s_end"} > $en || $current{"s_start"} > $st){
							print "format error \n".$region;
							die;
						}
						$current{"s_end"} = $en;
					}elsif($lines[$jj] =~ /^Sbjct[\s]+(-+)[\s]+/){
						$current{"s_seq"} .= $1;
					
					}
				}
				my %tmp = %current;
				push(@ret,\%tmp);
				%current = ();
				if($multiflag == 1){
					$ii = $jj-1;
				}else{
					$ii = $jj;
					last;
				}
			}
			
		}else{
			$ii++;
		}
		
	}
	
	return \@ret;
}


1;
