#!/usr/bin/perl
use strict;

my @BBstart;
my @BBend;
my @BBsvtype;
my @SSid;
my @SSstart;
my @SSend;
my @SSsvtype;
my %hash;
my $filter=0.5

open(INFILE,"$ARGV[0]"); #bigblock
while(<INFILE>){
	#push @BBchr,(split /\t/)[0];
	push @BBstart,(split /\t/)[1];
	push @BBend,(split /\t/)[2];
	push @BBsvtype,(split /\t/)[3];
}
chomp @BBsvtype;
close INFILE;

open(INFILE,"$ARGV[1]"); #sample sorted
while(<INFILE>){
	push @SSid,(split /\t/)[0];
	push @SSstart,(split /\t/)[2];
	push @SSend,(split /\t/)[3];
	push @SSsvtype,(split /\t/)[4];
}
chomp @SSsvtype;
close INFILE;

open(INFILE,"$ARGV[2]"); #sample list
while(<INFILE>){
	push @sampleid,(split /\t/)[0];
}
chomp @sampleid;
close INFILE;

open(OUT,">>test.out.txt");
#analysis
for(my $j=0;$j<@BBstart;$j=$j+1){
	for (my $index=0;$index<@sampleid;$index++){
		$hash{"@sampleid[$index]"}=0;
	}
	for(my $i=0;$i<@SSstart;$i=$i+1){
		if(@SSend[$i] < @BBstart[$j]){
			#ABCD
			next;
		}elsif(@SSstart[$i] > @BBend[$j]){
			#CDAB
			last;
		}elsif(@SSstart[$i] <= @BBstart[$j] and @BBstart[$j] <= @SSend[$i] and @SSend[$i] <= @BBend[$j] and @SSsvtype[$I] eq @BBsvtype[$j]){
			#ACBD
			$AC=@BBstart[$j] - @SSstart[$i]+1;
			$CB=@SSend[$i] - @BBstart[$j]+1;
			$BD=@BBend[$j] - @SSend[$i]+1;
			$AB=@SSend[$i] - @SSstart[$i]+1;
			$CD=@BBend[$j] - @BBstart[$j]+1;
			if($CB > $filter * $AB and $CB > $filter * $CD){
				$hash{"@SSid[$i]"}=1;
			}
		}elsif(@BBstart[$j] <= @SSstart[$i] and @SSstart[$i] <= @BBend[$j] and @BBend[$j] <= @SSend[$i] and @SSsvtype[$I] eq @BBsvtype[$j]){
			#CADB
			$CA=@SSstart[$i] - @BBstart[$j]+1;
			$AD=@BBend[$j] - @SSstart[$i]+1;
			$DB=@SSend[$i] - @BBend[$j]+1;
			$AB=@SSend[$i] - @SSstart[$i]+1;
			$CD=@BBend[$j] - @BBstart[$j]+1;
			if($AD > $filter * $AB and $AD > $filter * $CD){
				$hash{"@SSid[$i]"}=1;
			}
		}elsif(@SSstart[$i] <= @BBstart[$j] and @BBstart[$j] <= @BBend[$j] and @BBend[$j] <= @SSend[$i] and @SSsvtype[$I] eq @BBsvtype[$j]){
			#ACDB
			$AC=@BBstart[$j] - @SSstart[$i]+1;
			$CD=@BBend[$j] - @BBstart[$j]+1;
			$DB=@SSend[$i] - @BBend[$j]+1;
			$AB=@SSend[$i] - @SSstart[$i]+1;
			if($CD > $filter * $AB){
				$hash{"@SSid[$i]"}=1;
			}
		}elsif(@BBstart[$j] <= @SSstart[$i] and @SSstart[$i] <= @SSend[$i] and @SSend[$i] <= @BBend[$j] and @SSsvtype[$I] eq @BBsvtype[$j]){
			#CABD
			$CA=@SSstart[$i] - @BBstart[$j]+1;
			$AB=@SSend[$i] - @SSstart[$i]+1;
			$BD=@BBend[$j] - @SSend[$i]+1;
			$CD=@BBend[$j] - @BBstart[$j]+1;
			if($AB > $filter * $CD){
				$hash{"@SSid[$i]"}=1;
			}
		}
	}
}
close (OUT);