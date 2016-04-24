'''
Created on 03.07.2014

@author: david
'''
from Transcript import Transcript


class Gene(object):
    '''
    represents a Gene with certain transcripts
    '''
    
    def __init__(self, geneId, chromosome, strand, geneType, names, exons, codingExons):
        self.geneId = geneId        #usually the ENSG-Number
        self.chromosome = chromosome
        self.strand = strand        # minus = false; plus = True
        self.geneType = geneType    #the source of the transcript
        self.names = names          #list of names
        self.exons = sorted(exons, reverse = not strand) 
        self.codingExons = sorted(codingExons, reverse = not strand)
        self.start = min(e[0] for e in self.exons)
        self.end = max(e[1] for e in self.exons)
        
        self.transcripts = []
    __slots__ =["geneId","geneType","chromosome","start","end","strand","names","exons","codingExons", "transcripts"]
    
    
    def addTranscript(self, transcript):
        '''
        adds a trinscript Object to the transcript list
        '''
        if type(transcript) is Transcript:
            self.transcripts.append(transcript)

    def __str__(self, *args, **kwargs):
        return "\n".join([self.geneId,str(self.names), "%s:%s-%s" % (str(self.chromosome),str(self.start), str(self.end))])
    
    def printInfo(self):
        '''
        prints all the important information like Gene name and all the exons
        '''
        print("GeneName: %s" % str(tuple(self.names)))
        print("GeneId: %s" % self.geneId)
        print("Position: %s:%s-%s %s" % (str(self.chromosome),str(self.start), str(self.end),"+" if self.strand else "-"))
        print("CDS: %s" % str(self.codingExons))
        print("Exons: %s" % str(self.exons))
        print("Biotype: %s" % self.geneType)