#!/bin/bash

#echo "start"
for i in $@
do	
	EXP=$(sed 's/\.fa.*//' <<< $(basename $i));
	if [ ! -d $(dirname $i)/rnaEditor ]; then
  		mkdir $(dirname $i)/rnaEditor
	fi
	echo "run RNAEditor on" $EXP 
	python RnaEdit.py -i $i -r /media/Storage/databases/rnaEditor_annotations/human/human_g1k_v37.fasta -s /media/Storage/databases/rnaEditor_annotations/human/dbsnp_135.b37.vcf -m /media/Storage/databases/rnaEditor_annotations/human/hapmap_3.3.b37.sites.vcf -g /media/Storage/databases/rnaEditor_annotations/human/1000G_omni2.5.b37.sites.vcf -e /media/Storage/databases/rnaEditor_annotations/human/NHLBI_Exome_Sequencing_Project_6500SI.vcf -a /media/Storage/databases/rnaEditor_annotations/human/Alu_repeats_noCHR.bed -G /media/Storage/databases/rnaEditor_annotations/human/genes.gtf -o $(dirname $i)/rnaEditor/$EXP -d /usr/local/bin/	&

done