#!/usr/bin/perl

my $overlap_filter=0.5;

my @array_chr;
my @bigblock_start;
my @bigblock_end;
my @be_start;
my @be_end;
my @be_cn;
#read bigblock
open(INFILE,"$ARGV[0]");
while(<INFILE>){
	push @array_chr,(split /\t/)[0];
	push @bigblock_start,(split /\t/)[1];
	push @bigblock_end,(split /\t/)[2];
}
chomp @bigblock_end;
close INFILE;

my $chr=@array_chr[0];
#read original sample data
open(INFILE,"$ARGV[1]");
while(<INFILE>){
	push @be_start,(split /\t/)[2];
	push @be_end,(split /\t/)[3];
	push @be_cn,(split /\t/)[4];
}
close INFILE;
chomp @be_cn;
#按bigblock依次进行遍历，对每个bigblock里所有的SV进行remerge
#print ("@bigblock_start\n");
open(OUT,">>$ARGV[0]_remerge");
for($bigblock_index=0;$bigblock_index<@bigblock_start;$bigblock_index++){
	#open(OUT,">>blocks/block_$block_index");
	my @block_start;
	my @block_end;
	my @block_value;
	my @block_cn;
	#获取该bigblock内的所有SV，并对所有SV加权为1
	for($be_index=0;$be_index<@be_start;$be_index++){
		if(@be_start[$be_index] >= @bigblock_start[$bigblock_index] and @be_end[$be_index] <= @bigblock_end[$bigblock_index]){
			#print OUT @be_id[$be_index]."\t".@be_cn[$be_index]."\t".@be_start[$be_index]."\t".@be_end[$be_index]."\n";
			push @block_start, @be_start[$be_index];
			push @block_end, @be_end[$be_index];
			push @block_cn, @be_cn[$be_index];
			push @block_value, 1;
		}
	}
	#通过两两比对，找出与其他SV的overlap总值最高的region，将其设置为core，并将其他region按权重向上靠拢，以加权平均值计算边界
	my $cal=1;
	while($cal){
		undef $tmp_merge_core;
		undef @tmp_merge_list;
		$out_sum_overlap_val=0;
		#open(OUT,">>block141_test.txt");
		for($i=0;$i<@block_start;$i++){
			#$samp=$i;
			$tmp_sum_overlap_val=0;
			for($j=0;$j<@block_start;$j++){
				#print @block_start[$i]."\t".@block_end[$i]."\t".@block_start[$j]."\t".@block_end[$j]."\n";
				if($i == $j){
					next;
				}
				if(@block_end[$i] <= @block_start[$j] or @block_end[$j] <= @block_start[$i] or @block_cn[$i] != @block_cn[$j]){
					next;
				}
				if(@block_start[$i] <= @block_start[$j] and @block_end[$i] <= @block_end[$j]){
					$overlap_val = (@block_end[$i] - @block_start[$j] + 1) / (@block_end[$j] - @block_start[$i] + 1);
					if($overlap_val >= $overlap_filter){
						$tmp_sum_overlap_val+=$overlap_val;
					}
				}elsif(@block_start[$j] <= @block_start[$i] and @block_end[$j] <= @block_end[$i]){
					$overlap_val = (@block_end[$j] - @block_start[$i] + 1) / (@block_end[$i] - @block_start[$j] + 1);
					if($overlap_val >= $overlap_filter){
						$tmp_sum_overlap_val+=$overlap_val;
					}
				}elsif(@block_start[$i] <= @block_start[$j] and @block_end[$j] <= @block_end[$i]){
					$overlap_val = (@block_end[$j] - @block_start[$j] + 1) / (@block_end[$i] - @block_start[$i] + 1);
					if($overlap_val >= $overlap_filter){
						$tmp_sum_overlap_val+=$overlap_val;
					}
				}elsif(@block_start[$j] <= @block_start[$i] and @block_end[$i] <= @block_end[$j]){
					$overlap_val = (@block_end[$i] - @block_start[$i] + 1) / (@block_end[$j] - @block_start[$j] + 1);
					if($overlap_val >= $overlap_filter){
						$tmp_sum_overlap_val+=$overlap_val;
					}
				}else{
					print "FAIL:cannot calculate overlap!\n";
					exit (1);
				}
			}
			#print "******************************\n".$tmp_sum_overlap_val."\t".$out_sum_overlap_val."\n";
			if ($tmp_sum_overlap_val > $out_sum_overlap_val){
				$out_sum_overlap_val = $tmp_sum_overlap_val;
				$tmp_merge_core = $i;
				#print $tmp_sum_overlap_val."\n";
				#print $tmp_sum_overlap_val."\t".$out_sum_overlap_val."\n";
				#print "CORE:\t".$tmp_merge_core."\n";
			}
		}
		#print $tmp_merge_core."\n";
		if($out_sum_overlap_val == 0){
			last;
		}

		push @tmp_merge_list, $tmp_merge_core;
		for($j=0;$j<@block_start;$j++){
			if($tmp_merge_core == $j){
				next;
			}
			if(@block_end[$tmp_merge_core] <= @block_start[$j] or @block_end[$j] <= @block_start[$tmp_merge_core] or @block_cn[$tmp_merge_core] != @block_cn[$j]){
				next;
			}
			if(@block_start[$tmp_merge_core] <= @block_start[$j] and @block_end[$tmp_merge_core] <= @block_end[$j]){
				$overlap1 = (@block_end[$tmp_merge_core] - @block_start[$j] + 1) / (@block_end[$tmp_merge_core] - @block_start[$tmp_merge_core] + 1);
				$overlap2 = (@block_end[$tmp_merge_core] - @block_start[$j] + 1) / (@block_end[$j] - @block_start[$j] + 1);
				if($overlap1 >= $overlap_filter and $overlap2 >= $overlap_filter){
					push @tmp_merge_list, $j;
				}
			}elsif(@block_start[$j] <= @block_start[$tmp_merge_core] and @block_end[$j] <= @block_end[$tmp_merge_core]){
				$overlap1 = (@block_end[$j] - @block_start[$tmp_merge_core] + 1) / (@block_end[$tmp_merge_core] - @block_start[$tmp_merge_core] + 1);
				$overlap2 = (@block_end[$j] - @block_start[$tmp_merge_core] + 1) / (@block_end[$j] - @block_start[$j] + 1);
				if($overlap1 >= $overlap_filter and $overlap2 >= $overlap_filter){
					push @tmp_merge_list, $j;
				}
			}elsif(@block_start[$tmp_merge_core] <= @block_start[$j] and @block_end[$j] <= @block_end[$tmp_merge_core]){
				$overlap1 = (@block_end[$j] - @block_start[$j] + 1) / (@block_end[$tmp_merge_core] - @block_start[$tmp_merge_core] + 1);
				if($overlap1 >= $overlap_filter){
					push @tmp_merge_list, $j;
				}
			}elsif(@block_start[$j] <= @block_start[$tmp_merge_core] and @block_end[$tmp_merge_core] <= @block_end[$j]){
				$overlap1 = (@block_end[$tmp_merge_core] - @block_start[$tmp_merge_core] + 1) / (@block_end[$j] - @block_start[$j] + 1);
				if($overlap1 >= $overlap_filter){
					push @tmp_merge_list, $j;
				}
			}else{
				print "FAIL:cannot calculate overlap!\n";
				exit (1);
			}
		}
		#print ("@tmp_merge_list\n")."\n";
		my $block_out_start=0;
		my $block_out_end=0;
		my $block_out_value=0;
		my $block_out_cn=@block_cn[$tmp_merge_core];
		#print $count."\n";
		@sorted_tmp_merge_list = sort { $a <=> $b } @tmp_merge_list;
		#print ("@sorted_tmp_merge_list\n")."\n";
		#print "BLOCK::\n";
		#for($t=0;$t<@block_start;$t++){
			#print @block_start[$t]."\t".@block_end[$t]."\t".@block_value[$t]."\t".@block_cn[$t]."\n";
		#}
		for($i=0;$i<@tmp_merge_list;$i++){
			$block_out_start+=@block_start[@sorted_tmp_merge_list[$i]-$i] * @block_value[@sorted_tmp_merge_list[$i]-$i];
			#print @block_start[@sorted_tmp_merge_list[$i]-$i]."\t".@block_end[@sorted_tmp_merge_list[$i]-$i]."\t".@block_value[@sorted_tmp_merge_list[$i]-$i]."\t".@block_cn[@sorted_tmp_merge_list[$i]-$i]."\n";
			$block_out_end+=@block_end[@sorted_tmp_merge_list[$i]-$i] * @block_value[@sorted_tmp_merge_list[$i]-$i];
			$block_out_value+=@block_value[@sorted_tmp_merge_list[$i]-$i];
			splice @block_start, @sorted_tmp_merge_list[$i]-$i, 1;
			splice @block_end, @sorted_tmp_merge_list[$i]-$i, 1;
			splice @block_value, @sorted_tmp_merge_list[$i]-$i, 1;
			splice @block_cn, @sorted_tmp_merge_list[$i]-$i, 1;
		}
		$block_out_start = $block_out_start / $block_out_value;
		$block_out_end = $block_out_end / $block_out_value;
		push @block_start, round($block_out_start);
		push @block_end, round($block_out_end);
		push @block_value, $block_out_value;
		push @block_cn, $block_out_cn;
	#print ("@block_start\n")."\n";
	}
	my $block_count = @block_start;		#added 2019.2.19
	for($i=0;$i<@block_start;$i++){
		print OUT $chr."\t".@block_start[$i]."\t".@block_end[$i]."\t".@block_cn[$i]."\t".$block_count."\n";
	}
}
close (OUT);

#system("sort -n -k 1 -k 2 -k 3 $ARGV[0]_remerge > $ARGV[0]_remerge_sorted");
sub round()
{
	($l_number) = @_;

	$l_number = $l_number+0.5;
	$l_number = int($l_number);

	return $l_number;
}
