#!/bin/bash
WD=`pwd`
LOC=$1
rm vcf.list sample.list 2>/dev/null

cd ${LOC}
LISTLOC=`pwd`
LIST1=`ls ${LISTLOC}/*.vcf.gz`
LIST2=`ls *.vcf.gz`

for i in ${LIST1[@]}
do
	echo $i >> ${WD}/vcf.list
done
for j in ${LIST2[@]}
do
	echo $j | cut -f 1 -d '.' >> ${WD}/sample.list
done
