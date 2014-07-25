'''
Created on 05.06.2014

@author: david
'''
from Transcriptome import Transcriptome
from Helper import Helper


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
    
    
"""
startTime = Helper.getTime()
vcfFile="/media/Storage/databases/rnaEditor_annotations/human/dbsnp_135.b37.vcf"
variantsB = "/media/Storage/bio-data/David/Kostas/scrambleN/scrambleN.vcf"

Helper.readVcfFile(vcfFile)

"""
a= Helper.removeVariantsAFromVariantsB(vcfFile, variantsB)
i=0
for chr in a.keys():
    for variant in a[chr]:
        i+=1
duration = Helper.getTime() - startTime
print "\t[DONE]" + " Duration [" + str(duration) + "] " + str(i)
print "finished"
"""