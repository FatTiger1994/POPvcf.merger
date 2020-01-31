#!/usr/bin/perl

$filter=$ARGV[3];

#open(INFILE,"$ARGV[2]");
#@genotype_list=<INFILE>;
#chomp @genotype_list;
#close INFILE;
#my $counter=0;
open(OUT,">>$ARGV[2]_genotype_list");
open(OUT2,">>$ARGV[2]_genotype_list_index");

@chrom = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24);
foreach $chr (@chrom){
	#get array ABCD
	
	my @infile_start;
	my @infile_end;
	my @infile_cn;
	my @ref_start;
	my @ref_end;
	my @ref_cn;
	#my @genotype_list;
	#infile=merged list
	open(INFILE,"$ARGV[0]_chr$chr.txt");
	while(<INFILE>){
		push @infile_start,(split /\t/)[1];
		push @infile_end,(split /\t/)[2];
		push @infile_cn,(split /\t/)[3];
	}
	chomp @infile_cn;
	close INFILE;
	#ref=sample's region
	open(INFILE,"$ARGV[1]_chr$chr.txt");
	while(<INFILE>){
		push @ref_start,(split /\t/)[1];
		push @ref_end,(split /\t/)[2];
		push @ref_cn,(split /\t/)[3];
	}
	close INFILE;
	chomp @ref_cn;

	#open(INFILE,"$ARGV[2]");
	#@genotype_list=<INFILE>;
	#chomp @genotype_list;
	#close INFILE;
	
	#analysis

	for($i=0;$i<@infile_start;$i=$i+1){
		my $cal=0;
		my $catch=0;
		#my $sv_dif=0;
		my $tmp_start=@infile_start[$i];
		my $tmp_end=@infile_end[$i];
		for($j=0;$j<@ref_start;$j=$j+1){
			if(@infile_end[$i] < @ref_start[$j]){
				#ABCD
				next;
			}elsif(@infile_start[$i] > @ref_end[$j]){
				#CDAB
				next;
			}elsif(@infile_start[$i] <= @ref_start[$j] and @ref_start[$j] <= @infile_end[$i] and @infile_end[$i] <= @ref_end[$j]){
				#ACBD
				$catch=1;
				if(@infile_cn[$i] ne @ref_cn[$j]){
					next;
				}
				$CB=@infile_end[$i] - @ref_start[$j]+1;
				$AB=@infile_end[$i] - @infile_start[$i]+1;
				$CD=@ref_end[$j] - @ref_start[$j]+1;
				if($CB > $filter * $AB and $CB > $filter * $CD){
					$cal=1;
					$tmp_start=@ref_start[$j];
					$tmp_end=@ref_end[$j];
					last;
				}
			}elsif(@ref_start[$j] <= @infile_start[$i] and @infile_start[$i] <= @ref_end[$j] and @ref_end[$j] <= @infile_end[$i]){
				#CADB
				$catch=1;
				if(@infile_cn[$i] ne @ref_cn[$j]){
					next;
				}
				$AD=@ref_end[$j] - @infile_start[$i]+1;
				$AB=@infile_end[$i] - @infile_start[$i]+1;
				$CD=@ref_end[$j] - @ref_start[$j]+1;
				if($AD > $filter * $AB and $AD > $filter * $CD){
					$cal=1;
					$tmp_start=@ref_start[$j];
					$tmp_end=@ref_end[$j];
					last;
				}
			}elsif(@infile_start[$i] <= @ref_start[$j] and @ref_start[$j] <= @ref_end[$j] and @ref_end[$j] <= @infile_end[$i]){
				#ACDB
				$catch=1;
				if(@infile_cn[$i] ne @ref_cn[$j]){
					next;
				}
				$CD=@ref_end[$j] - @ref_start[$j]+1;
				$AB=@infile_end[$i] - @infile_start[$i]+1;
				if($CD > $filter * $AB){
					$cal=1;
					$tmp_start=@ref_start[$j];
					$tmp_end=@ref_end[$j];
					last;
				}
			}elsif(@ref_start[$j] <= @infile_start[$i] and @infile_start[$i] <= @infile_end[$i] and @infile_end[$i] <= @ref_end[$j]){
				#CABD
				$catch=1;
				if(@infile_cn[$i] ne @ref_cn[$j]){
					next;
				}
				$CD=@ref_end[$j] - @ref_start[$j]+1;
				$AB=@infile_end[$i] - @infile_start[$i]+1;
				if($AB > $filter * $CD){
					$cal=1;
					$tmp_start=@ref_start[$j];
					$tmp_end=@ref_end[$j];
					last;
				}
			}
		}
		if($chr == 23){
			$chr="X";
		}elsif($chr == 24){
			$chr="Y";
		}
		if($cal == 1){
			print OUT "chr".$chr.":".$tmp_start."-".$tmp_end."\n";
			print OUT2 "match\t".@infile_cn[$i]."\n";
		}elsif($catch == 1){
			print OUT "chr".$chr.":".$tmp_start."-".$tmp_end."\n";
			print OUT2 "mismatch\t".@infile_cn[$i]."\n";
		}else{
			print OUT "chr".$chr.":".$tmp_start."-".$tmp_end."\n";
			print OUT2 "unmatch\t".@infile_cn[$i]."\n";
		}
		#$counter=$counter+1;
	}
}
print OUT "exit\nEOF\n";
close (OUT);
close (OUT2);
