'''
Created on 05.06.2014

@author: david
'''

#!/usr/bin/env python

from VariantSet import VariantSet

var = VariantSet("/media/ATLAS_NGS_storage/Till/fastq/demultiplexed/rnaEditor/Icm4.vcf")
bedFile="/media/Storage/databases/rnaEditor_annotations/human/Alu_repeats_noCHR.bed"

alu,nonAlu=var.splitByBed(bedFile)

savedAlu=VariantSet("/media/ATLAS_NGS_storage/Till/fastq/demultiplexed/rnaEditor/Icm4.alu.vcf")
savedNonAlu=VariantSet("/media/ATLAS_NGS_storage/Till/fastq/demultiplexed/rnaEditor/Icm4.nonAlu.vcf")

print "original vars", len(var)
print "alu:",len(alu), "nonAlu:", len(nonAlu)
print "savedAlu:",len(savedAlu), "savedNonAlu:", len(savedNonAlu)