#!/bin/bash

#echo "start"
for i in $@
do	
	EXP=$(sed 's/\.fa.*//' <<< $(basename $i));
	if [ ! -d $(dirname $i)/rnaEditor ]; then
  		mkdir $(dirname $i)/rnaEditor
	fi
	echo "run RNAEditor on" $EXP 
	python RnaEdit.py -i $i -c configuration_single.txt	&

done
