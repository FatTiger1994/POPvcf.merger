#!/bin/bash
WD=`pwd`
TAG_SUF=".vcf.gz"
FLAG_LOC=0

if [ x$1 != x ]
then
	while getopts "i:s:" arg #选项后面的冒号表示该选项需要参数
	do
		case $arg in
		i)
			echo "	input file's location:
			$OPTARG" #参数存在$OPTARG中
			TAG_LOC=$OPTARG
			FLAG_LOC=1
			;;
		s)
			echo "	suffix:
			$OPTARG"
			TAG_SUF=$OPTARG
			;;
		?)  #当有不认识的选项的时候arg为?
			echo "FAIL: Unknown argument. Please check again."
			exit 1
			;;
		esac
	done

	until [ $FLAG_LOC -eq 1 ]
	do
		echo "FAIL: No location input. Please check again."
		exit 1
	done

	rm vcf.list sample.list 2>/dev/null
	
	cd ${TAG_LOC}
	LISTLOC=`pwd`
	LIST1=`ls ${LISTLOC}/*${TAG_SUF}`
	LIST2=`ls *${TAG_SUF}`
	
	for i in ${LIST1[@]}
	do
		echo $i >> ${WD}/vcf.list
	done
	for j in ${LIST2[@]}
	do
		echo $j | cut -f 1 -d '.' | cut -f 1 -d '_' >> ${WD}/sample.list
	done

else
	echo "USAGE
	POPvcf.list.generater.sh [-i ] [-s ]
	-i location of vcf files
	-s common suffex of vcf files
	example:POPvcf.list.generater.sh -i /picb/humpopg-bigdata/archive/20200102_1000G_2504_high_coverage_SV/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage_SV/working/20190628_WashU_Manta_Smoove/Manta/ -s .vcf.gz"
fi
