#!/usr/bin/perl
use strict;

my @BBchr;
my @BBstart;
my @BBend;
my @BBsvtype;
my @SSid;
my @SSstart;
my @SSend;
my @SSsvtype;
#my @SSGT;
my @sampleid;
my %hash;
my $filter=0.5;

open(INFILE,"$ARGV[0]"); #bigblock
while(<INFILE>){
	push @BBchr,(split /\t/)[0];
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
	#push @SSGT,(split /\t/)[5];
}
chomp @SSsvtype;
close INFILE;

open(INFILE,"$ARGV[2]"); #sample list
while(<INFILE>){
	push @sampleid,(split /\t/)[0];
}
chomp @sampleid;
close INFILE;

open(OUT,">>out.vcf");
#analysis
for(my $j=0;$j<@BBstart;$j=$j+1){
	for (my $index=0;$index<@sampleid;$index++){
		$hash{"@sampleid[$index]"}="0";
	}
	for(my $i=0;$i<@SSstart;$i=$i+1){
		if(@SSsvtype[$i] ne @BBsvtype[$j]){
			next;
		}else{
			if(@SSend[$i] < @BBstart[$j]){
				#ABCD
				next;
			}elsif(@SSstart[$i] > @BBend[$j]){
				#CDAB
				last;
			}elsif(@SSstart[$i] <= @BBstart[$j] and @BBstart[$j] <= @SSend[$i] and @SSend[$i] <= @BBend[$j]){
				#ACBD
				my $AC=@BBstart[$j] - @SSstart[$i]+1;
				my $CB=@SSend[$i] - @BBstart[$j]+1;
				my $BD=@BBend[$j] - @SSend[$i]+1;
				my $AB=@SSend[$i] - @SSstart[$i]+1;
				my $CD=@BBend[$j] - @BBstart[$j]+1;
				if($CB > $filter * $AB and $CB > $filter * $CD){
					$hash{"@SSid[$i]"}=1;
				}
			}elsif(@BBstart[$j] <= @SSstart[$i] and @SSstart[$i] <= @BBend[$j] and @BBend[$j] <= @SSend[$i]){
				#CADB
				my $CA=@SSstart[$i] - @BBstart[$j]+1;
				my $AD=@BBend[$j] - @SSstart[$i]+1;
				my $DB=@SSend[$i] - @BBend[$j]+1;
				my $AB=@SSend[$i] - @SSstart[$i]+1;
				my $CD=@BBend[$j] - @BBstart[$j]+1;
				if($AD > $filter * $AB and $AD > $filter * $CD){
					$hash{"@SSid[$i]"}=1;
				}
			}elsif(@SSstart[$i] <= @BBstart[$j] and @BBstart[$j] <= @BBend[$j] and @BBend[$j] <= @SSend[$i]){
				#ACDB
				my $AC=@BBstart[$j] - @SSstart[$i]+1;
				my $CD=@BBend[$j] - @BBstart[$j]+1;
				my $DB=@SSend[$i] - @BBend[$j]+1;
				my $AB=@SSend[$i] - @SSstart[$i]+1;
				if($CD > $filter * $AB){
					$hash{"@SSid[$i]"}=1;
				}
			}elsif(@BBstart[$j] <= @SSstart[$i] and @SSstart[$i] <= @SSend[$i] and @SSend[$i] <= @BBend[$j]){
				#CABD
				my $CA=@SSstart[$i] - @BBstart[$j]+1;
				my $AB=@SSend[$i] - @SSstart[$i]+1;
				my $BD=@BBend[$j] - @SSend[$i]+1;
				my $CD=@BBend[$j] - @BBstart[$j]+1;
				if($AB > $filter * $CD){
					$hash{"@SSid[$i]"}=1;
				}
			}
		}
	}
	print OUT @BBchr[$j]."\t".@BBstart[$j]."\t".@BBend[$j]."\t".@BBsvtype[$j];
	for (my $index=0;$index<@sampleid;$index++){
		print OUT "\t".$hash{"@sampleid[$index]"};
	}
	print OUT "\n";
}
close (OUT);