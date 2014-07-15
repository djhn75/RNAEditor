'''
Created on 05.06.2014

@author: david
'''
from Transcriptome import Transcriptome



gtfFile = "/media/Storage/databases/rnaEditor_annotations/human/genes_small.gtf"

#features = Genes.construct_gene_list_from_file(gtfFile, log)

a = Transcriptome()
a.createTranscriptomeFromFile(gtfFile)

#for geneType in a.geneByTypes:
for gene in a.geneByTypes["protein_coding"]:
    
    #for transcript in gene.transcripts:
    #print str(transcript.names) + "\t" + str(transcript.start) + "\t" + str(transcript.stop)
    print gene