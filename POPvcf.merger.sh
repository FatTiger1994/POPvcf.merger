w#!/bin/bash

WD=`pwd`
bin="/picb/humpopg-bigdata5/wangyimin/scripts/POPvcf.merger"
gap_loc="/picb/humpopg-bigdata5/wangyimin/scripts/PGG.SVpipeline.filter.annotate/gap_filterer/hg19.Gap.bed"

function BE_block_merge()
{
	${bin}/POPvcf.bigblock.pl ${WD}/merge/cnvr_nogap_chr${1}_sorted.txt ${WD}/merge
	${bin}/POPvcf.bigblock.remerge.pl ${WD}/merge/cnvr_bigblock_chr${1}.txt ${WD}/merge/cnvr_nogap_chr${1}_sorted.txt
}

FLAG_vcf_list=0
#FLAG_sample_index=0
FLAG_sample_list=0
FLAG_root_list=0
#FLAG_filter=0
#overlapping filter of sample & merged list (default 0.5)
filter=0.5

if [ x$1 != x ]
then
	while getopts "v:s:g" arg #选项后面的冒号表示该选项需要参数
	do
		case $arg in
		v)
			vcf_list=$OPTARG
			FLAG_vcf_list=1
			;;
		s)
			sample_list=$OPTARG
			FLAG_sample_list=1
			;;
		g)
			genotype=1
			echo "Input files have genotype infomation..."
			;;
		?)  #当有不认识的选项的时候arg为?
			echo "FAIL: Unknown argument. Please check again."
			exit 1
			;;
		esac
	done
	
	until [ $FLAG_vcf_list -eq 1 -a $FLAG_sample_list -eq 1 ]
	do
		echo "FAIL: Variables uncomplete. Please check again."
		exit 1
	done

	CHRs=$(echo {1..22} X Y)
	
	mkdir ${WD}/merge 2>/dev/null
	cd ${WD}/merge
	echo "Merging CNVR without gap regions..."

	for i in ${CHRs}
	do
		# generate cnvr_nogap_chr*.txt
		${bin}/POPvcf.CNVR.pl ${vcf_list} ${sample_list} ${gap_loc} ${WD}/merge ${i} &
	done

	wait
	echo "done!"

	for i in ${CHRs}
	do
		echo "merging chr "${i}
		# cnvr_nogap_chr*.txt ----> cnvr_nogap_chr*_sorted.txt
		cat ${WD}/merge/cnvr_nogap_chr${i}.txt | sort -k 5,5 -k 2n -k 3n -k 4n > ${WD}/merge/cnvr_nogap_chr${i}_sorted.txt
		# cnvr_nogap_chr*_sorted.txt ----> cnvr_bigblock_chr*.txt
		# cnvr_bigblock_chr*.txt ----> cnvr_bigblock_chr*.txt_remerge
		BE_block_merge $i &
	done

	wait
	echo "done!"
	cd ${WD}
	echo "#CHR	start	end	svtype" > out.vcf.tmp
	cat ${sample_list} >> out.vcf.tmp
	cat out.vcf.tmp | sed ":a;N;s/\n/\t/g;ta" >> out.vcf
	rm out.vcf.tmp

	for i in ${CHRs}
	do
		sort -n -k 1 -k 2 -k 3 ${WD}/merge/cnvr_bigblock_chr${i}.txt_remerge > ${WD}/merge/cnvr.bigblock.chr${i}.remerge.sorted
		if [ $genotype -eq 1 ]
		then
			${bin}/POPvcf.summary.withGT.pl ${WD}/merge/cnvr.bigblock.chr${i}.remerge.sorted ${WD}/merge/cnvr_nogap_chr${i}_sorted.txt ${sample_list}
		else
			${bin}/POPvcf.summary.pl ${WD}/merge/cnvr.bigblock.chr${i}.remerge.sorted ${WD}/merge/cnvr_nogap_chr${i}_sorted.txt ${sample_list}
		fi
	done

	#cat ${WD}/merged_genotype_list | cut -f 1,2,3 | sed 's/^23/X/g' | sed 's/^24/Y/g' > ${WD}/output/CNVR.txt
	#echo "#chr start end "$SAMPLEs | sed 's/ /\t/g' > ${WD}/output/out.txt
	#paste ${WD}/output/CNVR.txt ${WD}/genotype/genotype_rechecked.txt >> ${WD}/output/out.txt
	##added 2019.2.19
	#echo "#chr	start	end	Block_Complexity" > ${WD}/output/Block_Complexity.txt
	#cat ${WD}/merged_genotype_list | cut -f 1,2,3,5 >> ${WD}/output/Block_Complexity.txt

	echo "All tasks complete!"
else
	echo "	[USAGE]
	POPvcf.merger.sh [-v ] [-s ]
	-v dir of vcf_list (one column, list of vcf files, should be named vcf_list.txt)
	-s dir of sample_list (one column, sampleid)
	-g use -g if genotype is offered
	example: POPvcf.merger.sh -v /picb/humpopg-bigdata5/wangyimin/NGS/svmap/KGP_HighCov/20190628_WashU_Manta/vcf.list -s /picb/humpopg-bigdata5/wangyimin/NGS/svmap/KGP_HighCov/20190628_WashU_Manta/sample.list -g"
fi
