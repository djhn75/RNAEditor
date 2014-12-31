'''
Created on 05.06.2014

@author: david
'''
from Genome import Genome
from Helper import Helper
from collections import defaultdict
from VariantSet import Variant
from VariantSet import VariantSet
from CallEditingSites import CallEditingSites
import gtfHandler
from itertools import izip
from copy import copy, deepcopy
import sys

#vcfFile=sys.argv[1]
#vcfFile="/media/Storage/bio-data/David/Kostas/scrambleN/scrambleN.nonAlu_Y.vcf"

"""
ces = CallEditingSites(bamFile="/media/Storage/bio-data/David/Kostas/scrambleN/scrambleN.realigned.marked.recalibrated.bam",
                        refGenome="/media/Storage/databases/rnaEditor_annotations/human/human_g1k_v37.fasta", 
                        dbsnp="/media/Storage/databases/rnaEditor_annotations/human/dbsnp_135.b37_Y.vcf", 
                        hapmap="/media/Storage/databases/rnaEditor_annotations/human/hapmap_3.3.b37.sites.vcf",
                        omni="/media/Storage/databases/rnaEditor_annotations/human/1000G_omni2.5.b37.sites.vcf",
                        esp="/media/Storage/databases/rnaEditor_annotations/human/NHLBI_Exome_Sequencing_Project_6500SI.vcf", 
                        aluRegions="/media/Storage/databases/rnaEditor_annotations/human/Alu_repeats_noCHR.bed", 
                        gtfFile="/media/Storage/databases/rnaEditor_annotations/human/genes_Y.gtf", 
                        outfilePrefix=vcfFile[:vcfFile.rfind(".")], 
                        sourceDir="/usr/local/bin/")
"""
g=Genome("/media/Storage/databases/rnaEditor_annotations/human/genes_Y.gtf")
variants= VariantSet("/media/Storage/bio-data/David/Kostas/scrambleN/scrambleN_Y.vcf")


#variants.deleteOverlappsFromVcf("/media/Storage/databases/rnaEditor_annotations/human/dbsnp_135.b37_Y.vcf")
variants.annotateVariantDict(g)
variants.printGeneList(g, "dink.gvf", True)
gbc=g.getGenesByChromosome()

for gene in gbc["Y"]:
    print gene.geneId
#print len(gbc["Y"])

#ces.removeHomopolymers(variants, "/media/Storage/bio-data/David/Kostas/scrambleN/scrambleN.nonAlu_Y.vcfs", 4)

#ces.start()
