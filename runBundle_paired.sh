#!/bin/bash

#echo "start"
for i in $@
do	
	EXP=${i%%_1.*fa*}
	FA1=$i
	FA2=$EXP"_2.fastq.trimmed.fastq"	
	if [ ! -d $(dirname $i)/rnaEditor ]; then
  		mkdir $(dirname $i)/rnaEditor
	fi
	echo "run RNAEditor on $EXP with $FA1 and $FA2"
	python RnaEdit.py -i $i $FA2 -c configuration_paired.txt &

done
