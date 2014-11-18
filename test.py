'''
Created on 05.06.2014

@author: david
'''
from Transcriptome import Transcriptome
from Helper import Helper
from collections import defaultdict
from vcfHandler import Variant
import vcfHandler
from CallEditingSites import CallEditingSites



gtfFile = "/media/Storage/databases/rnaEditor_annotations/human/genes_small.gtf"
transcriptome = Transcriptome()
genesByChromosome = transcriptome.createTranscriptomeFromFile(gtfFile)

variantsB = "/media/Storage/bio-data/David/Kostas/scrambleN/scrambleN.vcf"
variants = vcfHandler.parseVcfFile_variantSetByChromosome(variantsB)
vcfHandler.sortVariantDict(variants)

for key in genesByChromosome.keys():
    for gene in genesByChromosome[key]:
        print gene


genes = transcriptome.getGenesByGeneID()
gene = genes["ENSG00000163131"]
print gene
#print information about CTSS
print str(gene.codingExons)
gene.printInfo()

print "%s %d $d $s %s CDS:%s Exons:%s" % (gene.chromosome,gene.start,gene.end,gene.names[0],gene.geneType,str(gene.codingExons),str(gene.exons))

#for geneType in a.geneByTypes:
#genesByChromosome = transcriptome.getGenesByChromosome()
   

startTime = Helper.getTime()
vcfFile="/media/Storage/databases/rnaEditor_annotations/human/dbsnp_135.b37.vcf"
variantsB = "/media/Storage/bio-data/David/Kostas/scrambleN/scrambleN.vcf"

variants = vcfHandler.parseVcfFile_variantSetByChromosome(variantsB)

vcfHandler.sortVariantDict(variants)

""""a= Helper.removeVariantsAFromVariantsB(vcfFile, variantsB)
i=0
for chr in a.keys():
    for variant in a[chr]:
        i+=1
duration = Helper.getTime() - startTime
print "\t[DONE]" + " Duration [" + str(duration) + "] " + str(i)
print "finished"
"""


#scrambleN=vcfHandler.parseVcfFile("/media/Storage/bio-data/David/Kostas/scrambleN/scrambleN.vcf")


edit = CallEditingSites(bamFile="/media/Storage/bio-data/David/Kostas/scrambleN/scrambleN.realigned.marked.recalibrated.bam",
                        refGenome="/media/Storage/databases/rnaEditor_annotations/human/human_g1k_v37.fasta", 
                        dbsnp="/media/Storage/databases/rnaEditor_annotations/human/dbsnp_135.b37_chr8.vcf", 
                        hapmap="/media/Storage/databases/rnaEditor_annotations/human/hapmap_3.3.b37.sites.vcf",
                        omni="/media/Storage/databases/rnaEditor_annotations/human/1000G_omni2.5.b37.sites.vcf",
                        esp="/media/Storage/databases/rnaEditor_annotations/human/NHLBI_Exome_Sequencing_Project_6500SI.vcf", 
                        aluRegions="/media/Storage/databases/rnaEditor_annotations/human/Alu_repeats_noCHR.bed", 
                        geneAnnotationFile="/media/Storage/databases/rnaEditor_annotations/human/genes_small.gtf", 
                        outfilePrefix="/media/Storage/bio-data/David/Kostas/scrambleN/scrambleN", 
                        sourceDir="/usr/local/bin/", threads=12, standCall=0, standEmit=0,
                        edgeDistance=5, keepTemp=True, overwrite=False)




A = edit.start()

#print len(A)
#vcfHandler.printVariantDict(A,"dink.vcf")

print "finished"