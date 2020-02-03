#!/usr/bin/perl

my @samples;
my @vcf_list;

open(INDEL,"$ARGV[0]");
while(<INDEL>){
	push @chr,(split /\t/)[1];
	push @start,(split /\t/)[2];
	push @end,(split /\t/)[3];
	push @svtype,(split /\t/)[4];
}
close INDEL;
chomp @svtype;
my $chrom=@chr[1];
my @block_start;
my @block_end;

my $tmp_start=@start[0];
my $tmp_end=@end[0];
my $tmp_svtype=@svtype[0];
for($i=1;$i<@start;$i=$i+1){
	if(@start[$i] <= $tmp_end and @svtype[$i] eq $tmp_svtype){
		if(@end[$i] > $tmp_end){
			$tmp_end=@end[$i];
		}
	}else{
		push (@block_start,$tmp_start);
		push (@block_end,$tmp_end);
		push (@block_svtype,$tmp_svtype);
		$tmp_start=@start[$i];
		$tmp_end=@end[$i];
		$tmp_svtype=@svtype[$i];
	}
}
push (@block_start,$tmp_start);
push (@block_end,$tmp_end);
push (@block_svtype,$tmp_svtype);

open(OUT,">>$ARGV[1]/cnvr_bigblock_chr$chrom.txt");
for($j=0;$j<@block_start;$j=$j+1){
	print OUT $chrom."\t".@block_start[$j]."\t".@block_end[$j]."\t".@block_svtype."\n";
}
close (OUT);
