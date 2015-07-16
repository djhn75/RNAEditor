'''
Created on 03.07.2014

@author: david
'''

class Transcript(object):
    '''
    represents a transcript (ENST-Number)
    '''

    def __init__(self, gene, transcriptId, names, protId, exonIndices, codingExonIndices,
                    codingFrames, startCodons, stopCodons,
                    startAnalysis = None, stop = None, codingStart = None, codingStop = None):


        self.gene = gene
        self.transcriptId = transcriptId
        self.names = names
        self.protId                = protId
        self.exonIndices           = exonIndices
        self.codingExonIndices    = codingExonIndices
        self.codingFrames          = codingFrames
        self.startCodons           = startCodons
        self.stopCodons            = stopCodons
        self.startAnalysis                    = min(self.gene.exons[e][0] for e in self.exonIndices) if startAnalysis == None else startAnalysis
        self.stop                    = max(self.gene.exons[e][1] for e in self.exonIndices) if stop == None else stop
        if codingStart != None:
            self.codingStart = codingStart
        else:
            self.codingStart      = min(self.gene.codingExons[e][0] for e in self.codingExonIndices) if len(self.codingExonIndices) else None
        if codingStop != None:
            self.codingStop = codingStop
        else:
            self.codingStop      = max(self.gene.codingExons[e][1] for e in self.codingExonIndices) if len(self.codingExonIndices) else None
    
    
       
        __slots__ = ["gene","transcriptId","protId", "startAnalysis","stop","codingStart",
             "codingStop","names","exonIndices","codingExonIndices",
             "codingFrames","virtualExonIndices","virtualCodingExonIndices",
             "startCodons","stopCodons"]
