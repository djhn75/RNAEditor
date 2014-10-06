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
"""
gtfFile = "/media/Storage/databases/rnaEditor_annotations/human/genes_small.gtf"

#features = Genes.construct_gene_list_from_file(gtfFile, log)

a = Transcriptome()
a.createTranscriptomeFromFile(gtfFile)

#for geneType in a.geneByTypes:
genesByChromosome = a.getGenesByChromosome()
for gene in genesByChromosome["1"]:
    
    #for transcript in gene.transcripts:
    #print str(transcript.names) + "\t" + str(transcript.start) + "\t" + str(transcript.stop)
    print gene
    
    

startTime = Helper.getTime()
vcfFile="/media/Storage/databases/rnaEditor_annotations/human/dbsnp_135.b37.vcf"
variantsB = "/media/Storage/bio-data/David/Kostas/scrambleN/scrambleN.vcf"

Helper.readVcfFile(vcfFile)


a= Helper.removeVariantsAFromVariantsB(vcfFile, variantsB)
i=0
for chr in a.keys():
    for variant in a[chr]:
        i+=1
duration = Helper.getTime() - startTime
print "\t[DONE]" + " Duration [" + str(duration) + "] " + str(i)
print "finished"
"""


#scrambleN=vcfHandler.parseVcfFile("/media/Storage/bio-data/David/Kostas/scrambleN/scrambleN.vcf")

scrambleN = open("/media/Storage/bio-data/David/Kostas/scrambleN/scrambleN.vcf")
dbSNP = open("/media/Storage/databases/rnaEditor_annotations/human/dbsnp_135.b37.vcf")


"""for line in scrambleN:
    if line.startswith("#"):
        continue
    A.add(getVariantTuble(line))

print len(A)
"""
"""
cmd= ['/usr/local/bin/samtools', 'view', '/media/Storage/bio-data/David/Kostas/scrambleN/scrambleN.realigned.marked.recalibrated.bam', '4:81000939-81000939']
a = Helper.getCommandOutput(cmd)
print a
"""

edit = CallEditingSites(bamFile="/media/Storage/bio-data/David/Kostas/scrambleN/scrambleN.realigned.marked.recalibrated.bam",
                        refGenome="/media/Storage/databases/rnaEditor_annotations/human/human_g1k_v37.fasta", 
                        dbsnp="/media/Storage/databases/rnaEditor_annotations/human/dbsnp_135.b37_chr8.vcf", 
                        hapmap="/media/Storage/databases/rnaEditor_annotations/human/hapmap_3.3.b37.sites.vcf",
                        omni="/media/Storage/databases/rnaEditor_annotations/human/1000G_omni2.5.b37.sites.vcf",
                        esp="/media/Storage/databases/rnaEditor_annotations/human/NHLBI_Exome_Sequencing_Project_6500SI.vcf", 
                        aluRegions="/media/Storage/databases/rnaEditor_annotations/human/Alu_repeats_noCHR.bed", 
                        geneAnnotationFile="/media/Storage/databases/rnaEditor_annotations/human/genes.gtf", 
                        outfilePrefix="/media/Storage/bio-data/David/Kostas/scrambleN/scrambleN", 
                        sourceDir="/usr/local/bin/", threads=12, standCall=0, standEmit=0,
                        edgeDistance=5, keepTemp=True, overwrite=False)




A = edit.start()

"""
for line in dbSNP:
    if line.startswith("#"):
        continue
    Btupple = getVariantTuble(line)
    if Btupple in A:
        #A.discard(Btupple)
        A.remove(Btupple)
        del Adict[Btupple]
        
"""
#print len(A)
#vcfHandler.printVariantDict(A,"dink.vcf")

print "finished"