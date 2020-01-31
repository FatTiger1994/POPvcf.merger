#!/usr/bin/perl
#if ($ARGV[2]){
#	$full_print=1;
#}else{
#	$full_print=0;
#}
$filter=$ARGV[3];

@chrom = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22);
foreach $chr (@chrom){
	#get array ABCD
	#my @infile;
	my @arraya;
	my @arrayb;
	my @arrayc;
	my @arrayd;
	
	open(INFILE,"$ARGV[0]/$ARGV[1]_chr$chr.txt");
	while(<INFILE>){
		push @arraya,(split /\t/)[1];
		push @arrayb,(split /\t/)[2];
	}
	close INFILE;
	open(INFILE,"$ARGV[0]/$ARGV[2]_chr$chr.txt");
	while(<INFILE>){
		push @arrayc,(split /\t/)[1];
		push @arrayd,(split /\t/)[2];
	}
	close INFILE;
	
	#analysis
	for($i=0;$i<@arraya;$i=$i+1){
		$cal=0;
		for($j=0;$j<@arrayc;$j=$j+1){
			if(@arrayb[$i] < @arrayc[$j]){
				#ABCD
				next;
			}elsif(@arraya[$i] > @arrayd[$j]){
				#CDAB
				next;
			}elsif(@arraya[$i] <= @arrayc[$j] and @arrayc[$j] <= @arrayb[$i] and @arrayb[$i] <= @arrayd[$j]){
				#ACBD
				#$AC=@arrayc[$j] - @arraya[$i]+1;
				$CB=@arrayb[$i] - @arrayc[$j]+1;
				#$BD=@arrayd[$j] - @arrayb[$i]+1;
				$AB=@arrayb[$i] - @arraya[$i]+1;
				$CD=@arrayd[$j] - @arrayc[$j]+1;
				if($CB > $filter * $AB and $CB > $filter * $CD){
					$cal=1;
					last;
				}
			}elsif(@arrayc[$j] <= @arraya[$i] and @arraya[$i] <= @arrayd[$j] and @arrayd[$j] <= @arrayb[$i]){
				#CADB
				#$CA=@arraya[$i] - @arrayc[$j]+1;
				$AD=@arrayd[$j] - @arraya[$i]+1;
				#$DB=@arrayb[$i] - @arrayd[$j]+1;
				$AB=@arrayb[$i] - @arraya[$i]+1;
				$CD=@arrayd[$j] - @arrayc[$j]+1;
				if($AD > $filter * $AB and $AD > $filter * $CD){
					$cal=1;
					last;
				}
			}elsif(@arraya[$i] <= @arrayc[$j] and @arrayc[$j] <= @arrayd[$j] and @arrayd[$j] <= @arrayb[$i]){
				#ACDB
				#$AC=@arrayc[$j] - @arraya[$i]+1;
				$CD=@arrayd[$j] - @arrayc[$j]+1;
				#$DB=@arrayb[$i] - @arrayd[$j]+1;
				$AB=@arrayb[$i] - @arraya[$i]+1;
				if($CD > $filter * $AB){
					$cal=1;
					last;
				}
			}elsif(@arrayc[$j] <= @arraya[$i] and @arraya[$i] <= @arrayb[$i] and @arrayb[$i] <= @arrayd[$j]){
				#CABD
				#$CA=@arraya[$i] - @arrayc[$j]+1;
				$AB=@arrayb[$i] - @arraya[$i]+1;
				#$BD=@arrayd[$j] - @arrayb[$i]+1;
				$CD=@arrayd[$j] - @arrayc[$j]+1;
				if($AB > $filter * $CD){
					$cal=1;
					last;
				}
			}
		}
		open(OUT,">>$ARGV[4]");
		#if($full_print==0){
		#chomp $refR;
		if($cal==1){
			print OUT "1\n";
		}else{
			print OUT "0\n";
		}
		#}else{
		#	if($cal==1){
		#		$abcov=sprintf "%.2f",$abcov;
		#		$cdcov=sprintf "%.2f",$cdcov;
		#		print OUT @infile[$i]."\t50%_overlap\t".$length."\t".$abcov."\t".$cdcov."\t".$meandis."\n";
		#	}else{
		#		print OUT @infile[$i]."\t.\t.\t.\t.\t.\n";
		#	}
		#}
	}
	close (OUT);
}
