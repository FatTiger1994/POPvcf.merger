#!/usr/bin/perl
	
my @genotype;
my @infile_end;
my @ref_start;
my @ref_end;
#my @genotype_list;
#infile=merged list
open(INFILE,"$ARGV[0]");
@genotype=<INFILE>;
chomp @genotype;
close INFILE;
#ref=sample's region
open(INFILE,"$ARGV[1]");
while(<INFILE>){
	push @index,(split /\t/)[0];
	push @index_cn,(split /\t/)[1];
}
chomp @index_cn;
close INFILE;
#open(INFILE,"$ARGV[2]");
#@genotype_list=<INFILE>;
#chomp @genotype_list;
#close INFILE;

#analysis
open(OUT,">>$ARGV[2]_rechecked");
for($i=0;$i<@genotype;$i=$i+1){
	if($index[$i] eq "match"){
		if(@index_cn[$i] eq "-1"){
			if(@genotype[$i] < 1.5){
				print OUT @genotype[$i]."\n";
			}else{
				print OUT "NA\n";
			}
		}elsif(@index_cn[$i] eq "1"){
			if(@genotype[$i] > 2.5){
				print OUT @genotype[$i]."\n";
			}else{
				print OUT "NA\n";
			}
		}else{
			print OUT @genotype[$i]."\n";
		}
	}elsif($index[$i] eq "mismatch"){
		print OUT "2\n";
	}elsif($index[$i] eq "unmatch"){
		if(@genotype[$i] < 2.5 and @genotype[$i] > 1.5){
			print OUT @genotype[$i]."\n";
		}else{
			print OUT "NA\n";
		}
	}else{
		print "ERROR!!\n";
		exit (1);
	}
}
close (OUT);
