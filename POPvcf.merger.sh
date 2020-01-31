#!/bin/bash

WD=`pwd`
bin="/picb/humpopg-bigdata5/wangyimin/scripts/POPvcf.merger"
gap_loc="/picb/humpopg-bigdata5/wangyimin/scripts/PGG.SVpipeline.filter.annotate/gap_filterer/hg19.Gap.bed"

function BE_block_merge()
{
	${bin}/get_bigblock.pl ${WD}/merge/cnvr_nogap_chr${1}_sorted.txt ${WD}/merge
	${bin}/get_bigblock_remerge.pl ${WD}/merge/cnvr_bigblock_chr${1}.txt ${WD}/merge/cnvr_nogap_chr${1}_sorted.txt
}

function check_time()
{
	local time_subroutine=$(getconf CLK_TCK)
	local start_time=$(awk '{print $22}' /proc/$1/stat 2>/dev/null)
	if [ ! $start_time ]
	then
		return 0
	fi
	local sys_uptime=$(awk '{print $1}' /proc/uptime)
	local pid_uptime=$((${sys_uptime%.*} - ${start_time}/${time_subroutine}))
	#echo "wait_time="$pid_uptime
	if [ ${pid_uptime} -ge ${time_limit} ]
	then
		kill -9 $1 
		return 2
	else
		return 1
	fi
}

function run_cnvnator()
{
	rootfile=`grep $1 ${root_list}`
	if [ ! -f $rootfile ]
	then
		echo $rootfile" doesn't exist! Please check again."
		exit
	fi
	cnvnator -root $rootfile -genotype 100 <${WD}/samplecheck/${1}_genotype_list >${WD}/genotype/${1}_genotype.txt 2>/dev/null&
}

function run_genotype()
{
	local tmp_done=0
	run_cnvnator $1
	local pid=$!
	sleep ${time_start}
	check_time $pid
	local return_value=$?
	while [ $tmp_done -eq 0 ]
	do
		while [ $return_value -eq 1 ]
		do
			sleep 10
			check_time $pid
			return_value=$?
		done
		if [ $return_value -eq 2 ]
		then
			echo "Task $1 had beed killed: Time out. Recalculating..."
			rm ${WD}/genotype/$1_genotype.txt
			tmp_done=1
			return 1
		else
			tmp_done=1
		fi
	done
	grep -v 'male' ${WD}/genotype/${1}_genotype.txt >${WD}/genotype/${1}_genotype_filtered.txt
	cut -d ' ' -f 4 ${WD}/genotype/${1}_genotype_filtered.txt >${WD}/genotype/${1}
	#return 1
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
	while getopts "v:s:r:" arg #选项后面的冒号表示该选项需要参数
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
		#r)
		#	root_list=$OPTARG
		#	FLAG_root_list=1
		#	;;
		#f)
		#	filter=$OPTARG
		#	#FLAG_filter=1
		#	;;
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
	
	if [ ! -f "${WD}/merged_genotype_list" ]
	then
		mkdir ${WD}/merge 2>/dev/null
		cd ${WD}/merge
		echo "1. Converting metaSV to CNVR without gap regions..."
		for i in {1..24}
		do
			# generate cnvr_nogap_chr*.txt
			${bin}/Convert_metaSV2nogap_CNVR.pl ${vcf_list} ${sample_list} ${gap_loc} ${WD}/merge ${i} &
		done
		wait
		echo "done!"

		CHRs=$(echo {1..22} X Y)
		for i in $CHRs
		do
			echo "merging chr "${i}
			# cnvr_nogap_chr*.txt ----> cnvr_nogap_chr*_sorted.txt
			sort -n -k 2 -k 3 -k 4 ${WD}/merge/cnvr_nogap_chr${i}.txt > ${WD}/merge/cnvr_nogap_chr${i}_sorted.txt
			# cnvr_nogap_chr*_sorted.txt ----> cnvr_bigblock_chr*.txt
			# cnvr_bigblock_chr*.txt ----> cnvr_bigblock_chr*.txt_remerge
			BE_block_merge $i &
		done
		wait
		echo "done!"

		cat ${WD}/merge/cnvr_bigblock_chr*.txt_remerge > ${WD}/merge/cnvr.bigblock.remerge.unsort
		sort -n -k 1 -k 2 -k 3 ${WD}/merge/cnvr.bigblock.remerge.unsort > ${WD}/merge/cnvr.bigblock.remerge.sorted
		#mergeContinuousBlock.pl vcf_list.BE4mat.merge.uniBlock_sorted
		#cd ${WD}
		cat ${WD}/merge/cnvr.bigblock.remerge.sorted > ${WD}/merged_genotype_list
		#cat genotype_list | sed 's/^/chr/g' | sed 's/chr23/chrX/g' | sed 's/chr24/chrY/g' | sed 's/\t/ /g' > genotype_list.txt
		#echo $'exit\nEOF' >> genotype_list.txt
	else
		echo "1. MetaSV samples merged."
	fi

	CNVR_count=`cat ${WD}/merge/cnvr.bigblock.remerge.sorted | wc -l`
	SAMPLEs=`cat ${sample_list}`
	time_limit=`expr $CNVR_count / 150`
	#time_limit=864000
	time_start=`expr $CNVR_count / 500 + 100`

	if [ ! -f ${WD}/samplecheck/sample_transformat.log ]
	then
		mkdir ${WD}/samplecheck 2>/dev/null
		cd ${WD}/samplecheck
		echo "2. MetaSV data transformating..."
		for i in {1..24}
		do
			grep "^${i}	" ${WD}/merge/cnvr.bigblock.remerge.sorted > ${WD}/samplecheck/genotype_list_chr${i}.txt&
			for samp in $SAMPLEs
			do
				grep "$samp" ${WD}/merge/cnvr_nogap_chr${i}.txt | cut -f 2,3,4,5 > ${WD}/samplecheck/${samp}_transformated_chr${i}.txt&
			done
			wait
		done
		timenow=`date`
		echo "Sample transformated!	$timenow" > ${WD}/samplecheck/sample_transformat.log
	else
		echo "2. MetaSV data transformated."
	fi
		
	if [ ! -f ${WD}/samplecheck/genotype_lists.log ]
	then
		cd ${WD}/samplecheck
		echo "3. Making genotype lists..."
		unfinish=1
		while [ $unfinish -eq 1 ]
		do
			unfinish=0
			tasks=0
			for samp in $SAMPLEs
			do	
				genotype_count_tmp=`cat ${WD}/samplecheck/${samp}_genotype_list_index 2>/dev/null | wc -l`
				if [ $genotype_count_tmp -eq $CNVR_count ]
				then
					continue
				else
					unfinish=1
					#rm ${WD}/samplecheck/${samp}_transformated_chr* 2>/dev/null
					rm ${WD}/samplecheck/${samp}_genotype_list* 2>/dev/null
					#for i in {1..24}
					#do
					#	grep "^${i}	" ${WD}/samplecheck/${samp}_transformated.txt > ${WD}/samplecheck/${samp}_transformated_chr${i}.txt
					#done
					${bin}/genotype.list.maker.pl ${WD}/samplecheck/genotype_list ${WD}/samplecheck/${samp}_transformated ${samp} ${filter}&
					tasks=`expr $tasks + 1`
					if [ $tasks -ge 100 ]
					then
						wait
						tasks=`expr $tasks - 100`
					fi
				fi
			done
			wait
			echo "Rechecking..."
		done
		timenow=`date`
		echo "${bin}/genotype.list.maker.pl ${WD}/samplecheck/genotype_list ${WD}/samplecheck/${samp}_transformated ${samp} ${filter}		$timenow" > ${WD}/samplecheck/genotype_lists.log
	else
		echo "3. Genotype lists exist."
	fi

	if [ ! -f "${WD}/genotype/genotype.txt" ]
	then
		mkdir ${WD}/genotype 2>/dev/null
		cd ${WD}/genotype
		echo "4. Running genotype..."
		tasks=0
		unfinish=1
		while [ $unfinish -eq 1 ]
		do
			unfinish=0
			for i in $SAMPLEs
			do
				sample_count=`cat ${i} 2>/dev/null | wc -l` 
				if [ $sample_count -eq $CNVR_count ]
				then
					continue
				else
					unfinish=1
					rm ${WD}/genotype/${i}* 2>/dev/null
					run_genotype $i &
					echo "---->"$i" "`date`
					tasks=`expr $tasks + 1`
					if [ $tasks -ge 10 ]
					then
						wait
						tasks=`expr $tasks - 10`
					fi
				fi
			done
			wait
			echo "Rechecking..."
		done
		echo "done!"

		echo "Merging all samples..."
		paste_SAMPLEs=(`cat ${sample_list}`)
		paste ${paste_SAMPLEs[@]} >> ${WD}/genotype/genotype.txt
		echo "done!"
	else
		echo "4. Genotype results exists."
	fi

	if [ ! -f "${WD}/genotype/genotype_rechecked.txt" ]
	then
		cd ${WD}/genotype
		echo "5. Rechecking genotype..."
		tasks=0
		unfinish=1
		while [ $unfinish -eq 1 ]
		do
			unfinish=0
			for i in $SAMPLEs
			do
				sample_count=`cat ${i}_rechecked 2>/dev/null | wc -l` 
				if [ $sample_count -eq $CNVR_count ]
				then
					continue
				else
					unfinish=1
					rm ${WD}/genotype/${i}_rechecked 2>/dev/null
					${bin}/genotype.rechecker.pl ${WD}/genotype/${i} ${WD}/samplecheck/${i}_genotype_list_index ${i}&
					tasks=`expr $tasks + 1`
					if [ $tasks -ge 30 ]
					then
						wait
						tasks=`expr $tasks - 30`
					fi
				fi
			done
			wait
			echo "Rechecking..."
		done
		j=0
		for i in $SAMPLEs
		do
			paste_SAMPLEs[$j]=$i"_rechecked"
			j=`expr $j + 1`
		done
		mkdir ${WD}/output 2>/dev/null
		paste ${paste_SAMPLEs[@]} >> ${WD}/genotype/genotype_rechecked.txt
		cp ${WD}/genotype/genotype_rechecked.txt ${WD}/output/genotype_reddchecked.txt
		echo "done!"
	else
		echo "5. Genotype rechecked."
	fi

	cat ${WD}/merged_genotype_list | cut -f 1,2,3 | sed 's/^23/X/g' | sed 's/^24/Y/g' > ${WD}/output/CNVR.txt
	echo "#chr start end "$SAMPLEs | sed 's/ /\t/g' > ${WD}/output/out.txt
	paste ${WD}/output/CNVR.txt ${WD}/genotype/genotype_rechecked.txt >> ${WD}/output/out.txt
	#added 2019.2.19
	echo "#chr	start	end	Block_Complexity" > ${WD}/output/Block_Complexity.txt
	cat ${WD}/merged_genotype_list | cut -f 1,2,3,5 >> ${WD}/output/Block_Complexity.txt

	echo "All tasks complete!"
else
	echo "	[USAGE]
	PGG.SVpipeline.merge.genotype.sh [-v ] [-s ] [-r ]
	-v dir of vcf_list (one column, list of vcf files, should be named vcf_list.txt)
	-s dir of sample_list (one column, sampleid)
	-r dir of root_list (one column, list of root files)
	example: PGG.SVpipeline.merge.genotype.sh -v /picb/humpopg-bigdata5/wangyimin/NGS/svmap/test/vcf_list.txt -s /picb/humpopg-bigdata5/wangyimin/NGS/svmap/test/sample.list -r /picb/humpopg-bigdata5/wangyimin/NGS/svmap/test/root_list.txt"
fi
