#!/usr/bin/perl

$vcf_file=${ARGV[0]};
$sample_file=${ARGV[1]};
$gap_loc=${ARGV[2]};
$out_loc=${ARGV[3]};
$chr=${ARGV[4]};
$filter=0;

my @vcf_list;
my @sample_list;
my @arraya;
my @arraycn;
my @arrayb;
my @infile;
#读取vcf list
open(INFILE,"${vcf_file}");
@vcf_list=<INFILE>;
chomp @vcf_list;
close INFILE;
#读取sample list
open(INFILE,"${sample_file}");
@sample_list=<INFILE>;
chomp @sample_list;
close INFILE;

#遍历samples
for($vcf_index=0;$vcf_index<@vcf_list;$vcf_index++){
	my @array_chr;
	my @array_start;
	my @array_svtype;
	my @array_end;
	my @vcf_file;

	#open(VCFFILE,@vcf_list[$vcf_index]);
	open(VCFFILE,"gzip -dc @vcf_list[$vcf_index] |") or die ("can not open file!\n");
	@vcf_file=<VCFFILE>;
	close VCFFILE;
	chomp @vcf_file;
	for($i=0;$i<@vcf_file;$i++){
		if (@vcf_file[$i] =~ /[#]{1,2}/){
			next;
		}elsif(@vcf_file[$i] !~ "SVTYPE="){
			next;			
		}elsif(@vcf_file[$i] !~ "END="){
			next;
		}else{
			@splitted_vcf_file=split("\t",@vcf_file[$i]);
			push (@array_chr,@splitted_vcf_file[0]);
			push (@array_start,@splitted_vcf_file[1]);
			push (@array_svtype,@splitted_vcf_file[7]);
			push (@array_end,@splitted_vcf_file[7]);
		}
	}
	foreach $string (@array_chr){
		$string=~ s/(chr)([0-9xXyY]{1,2})/$2/;
	}
	foreach $string (@array_svtype){
		$string=~ s/$/;/;
		$string=~ s/(.*;SVTYPE=)([A-Z]{3,5})(;.*)/$2/;
	}
	foreach $string (@array_end){
		$string=~ s/^/;/;
		$string=~ s/(.*;END=)([0-9]{3,9})(;.*)/$2/;
	}
	#foreach $string (@array_svtype){
	#	$string=~ s/DEL/-1/;
	#	$string=~ s/DUP/1/;
	#	$string=~ s/CNV/+-1/;
	#	$string=~ s/INS/i1/;
	#	$string=~ s/INV/0/;
	#}
	#遍历samples_chr
	for($array_index=0;$array_index<@array_chr;$array_index++){
		if(@array_chr[$array_index] ne ${chr}){
			next;
#		}elsif(@array_chr[$array_index] == ${chr}){
		}else{
			if(@array_svtype[$array_index] eq "INS"){
				@array_start[$array_index] = @array_start[$array_index]-50;
				@array_end[$array_index] = @array_end[$array_index]+50;
			}
			$tmp_singlesv=@sample_list[$vcf_index]."\t".$chr."\t".@array_start[$array_index]."\t".@array_end[$array_index]."\t".@array_svtype[$array_index];
			push (@infile,$tmp_singlesv);
			push (@arraya,@array_start[$array_index]);
			push (@arrayb,@array_end[$array_index]);
			push (@arraycn,@array_svtype[$array_index]);
#		}else{
#			print "Cannot recognize ".@array_chr[$array_index]."!!!\n";
#			exit (1);
		}
	}
}
close (OUT);

#获得对应chr的sampleID,start,end
#过滤gap

my @arrayc;
my @arrayd;
open(GAPFILE,"${gap_loc}_chr${chr}");
while(<GAPFILE>){
	push @arrayc,(split /\t/)[1];
	push @arrayd,(split /\t/)[2];
}
close GAPFILE;

open(OUT,">>${out_loc}/cnvr_nogap_chr${chr}.txt");
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
			$AC=@arrayc[$j] - @arraya[$i]+1;
			$CB=@arrayb[$i] - @arrayc[$j]+1;
			$BD=@arrayd[$j] - @arrayb[$i]+1;
			$AB=@arrayb[$i] - @arraya[$i]+1;
			$CD=@arrayd[$j] - @arrayc[$j]+1;
			if($CB > $filter * $AB){
				$cal=1;
				last;
			}
		}elsif(@arrayc[$j] <= @arraya[$i] and @arraya[$i] <= @arrayd[$j] and @arrayd[$j] <= @arrayb[$i]){
			#CADB
			$CA=@arraya[$i] - @arrayc[$j]+1;
			$AD=@arrayd[$j] - @arraya[$i]+1;
			$DB=@arrayb[$i] - @arrayd[$j]+1;
			$AB=@arrayb[$i] - @arraya[$i]+1;
			$CD=@arrayd[$j] - @arrayc[$j]+1;
			if($AD > $filter * $AB){
				$cal=1;
				last;
			}
		}elsif(@arraya[$i] <= @arrayc[$j] and @arrayc[$j] <= @arrayd[$j] and @arrayd[$j] <= @arrayb[$i]){
			#ACDB
			$AC=@arrayc[$j] - @arraya[$i]+1;
			$CD=@arrayd[$j] - @arrayc[$j]+1;
			$DB=@arrayb[$i] - @arrayd[$j]+1;
			$AB=@arrayb[$i] - @arraya[$i]+1;
			if($CD > $filter * $AB){
				$cal=1;
				last;
			}
		}elsif(@arrayc[$j] <= @arraya[$i] and @arraya[$i] <= @arrayb[$i] and @arrayb[$i] <= @arrayd[$j]){
			#CABD
			$cal=1;
			last;
		}
	}
	if($cal==0){
		print OUT @infile[$i]."\n";
	}
}
close (OUT);
