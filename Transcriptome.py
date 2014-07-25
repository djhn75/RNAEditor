'''
Created on 11.06.2014

@author: david
'''
from collections import defaultdict
import gtfHandler
import gzip
from Helper import Helper
from Gene import Gene
from itertools import izip
from array import array
from Transcript import Transcript

'''
classdocs
'''
class Transcriptome(object):
    
    def __init__(self):
        '''
        Constructor
        '''
        
        #Dict / List for all Genes
        self.geneList = []
        self.genesByChromosome = defaultdict(list)
        
        # list of all features
        self.featureTypes = set() #set of all possible featuresTypes
        self.geneTypes  = set() #set of all possible geneTypes (sources)


        #
        #   gene_id -> beg/end/cdna_id/exon_index
        #   gene_id -> (gene_name,contig,strand)
        #
        self.uniqGeneSet = set()
        self.uniqGene_to_source = dict()
        self.uniqGene_to_names          = defaultdict(set)
        self.uniqGene_to_transcriptIds       = defaultdict(set)
        

        #
        #   transcriptId
        #
        self.transcriptId_to_names        = defaultdict(set)
        self.transcriptId_to_protId      = dict()
        self.transcriptId_to_startCodons = defaultdict(tuple)
        self.transcriptId_to_stopCodons  = defaultdict(tuple)
        self.transcriptId_to_exons        = defaultdict(set)
        self.transcriptId_to_cds = defaultdict(set)
        self.transcriptId_to_codingFrames= defaultdict(set)

        self.transcriptIds        = defaultdict(set)
        
        
        
    def parseGtf(self,gtfFile):
        """
            fill the dictionarys with the correspending attributes
        """
        try:
            for lineNumber, f in enumerate(gtfHandler.iterator(gtfFile)):
                uniqGene = f.geneId, f.chr, f.strand #(ENSG,chr,strand)
                interval = (f.start, f.end)
                exonNumber = int(f.attributes["exon_number"]) - 1
                
                self.featureTypes.add(f.featureType)
                if "transcript_name" in f.attributes:
                    self.transcriptId_to_names[f.transcriptId].add(f.attributes["transcript_name"])
                if "gene_name" in f.attributes:
                    self.uniqGene_to_names[uniqGene].add(f.attributes["gene_name"])
        
                self.uniqGeneSet.add(uniqGene)
                self.uniqGene_to_source[uniqGene] = f.source
                self.geneTypes.add(f.source)
                self.uniqGene_to_transcriptIds[uniqGene].add(f.transcriptId)
                
                
                #react for different feature types
                if f.featureType == "exon":
                    self.transcriptId_to_exons[f.transcriptId].add(interval)
                elif f.featureType == "CDS":
                    proteinId = f.attributes["protein_id"]
                    if (uniqGene in self.transcriptId_to_protId and self.transcriptId_to_protId[uniqGene] != proteinId):
                        Helper.warning("Transcript [%s] has many Protein IDs" % uniqGene)
                    self.transcriptId_to_protId[f.transcriptId] = proteinId
                    self.transcriptId_to_codingFrames[f.transcriptId].add((exonNumber, f.frame))
                    self.transcriptId_to_cds[f.transcriptId].add(interval)
                
                elif f.featureType == "start_codon":
                    self.transcriptId_to_startCodons[f.transcriptId] += (exonNumber, interval),
                elif f.featureType == "stop_codon":
                    self.transcriptId_to_stopCodons[f.transcriptId] += (exonNumber, interval),
                    
        except:
            print lineNumber, "    " + f.geneId
            raise
        
    
    def assembleTranscriptome(self):
        transcriptsByType = defaultdict(list)
        
        #construct Genes
        for uniqGene in self.uniqGeneSet:
            geneId, chr, strand = uniqGene
            
            geneNames = list(self.uniqGene_to_names[uniqGene])
            geneType = self.uniqGene_to_source[uniqGene]
            geneExons = set()
            geneCds = set()
            
            #get exons from all transcripts
            for transcriptId in self.uniqGene_to_transcriptIds[uniqGene]:
                geneExons |= self.transcriptId_to_exons[transcriptId]
                geneCds  |= self.transcriptId_to_cds[transcriptId]
            
            #usually gtf Files are sorted, but this can't be assumed    
            geneExons = sorted(geneExons, reverse = not strand)    
            geneCds = sorted(geneCds, reverse = not strand)
            
            gene = Gene(geneId, chr, strand, geneType, geneNames, geneExons, geneCds)
            
            
            geneExons = dict(izip(geneExons,xrange(1000000)))
            geneCds = dict(izip(geneCds,xrange(1000000))) 
            
            
            self.geneList.append(gene)
            
           
            #construct transcripts
            for transcriptId in self.uniqGene_to_transcriptIds[uniqGene]:
                transcriptNames = self.transcriptId_to_names[transcriptId]
            
                protId = self.transcriptId_to_protId[transcriptId] if transcriptId in self.transcriptId_to_protId else None
                exons = sorted(self.transcriptId_to_exons[transcriptId])
                codingExons = sorted(self.transcriptId_to_cds[transcriptId])
                exonIndices = array('H', [geneExons[e] for e in exons])
                codingExonIndices = array('H', [geneCds[e] for e in codingExons])
                codingFrames = array('H', [int(frame) for exonNumber, frame in sorted(self.transcriptId_to_codingFrames[transcriptId])])
                startCodon = tuple([interval for exonNumber, interval in sorted(self.transcriptId_to_startCodons[transcriptId])])
                stopCodon = tuple([interval for exonNumber, interval in sorted(self.transcriptId_to_stopCodons[transcriptId])])
                
                if len(codingExons) != len(codingFrames):
                    raise Exception("Number of coding Exons and Frames differ for %s %s" % geneId, transcriptId)
                
                
                transcript = Transcript(gene, transcriptId,list(transcriptNames), protId, exonIndices, codingExonIndices, codingFrames, startCodon, stopCodon)
                
                gene.addTranscript(transcript)
            
        self.genesByChromosome = self.getGenesByChromosome()
            
    def createTranscriptomeFromFile(self,gtfFilePath):
        """
        Construct Transcriptome from GTF File
            Saves all the information in dictionarys
        """
        startTime = Helper.getTime()
        Helper.info(" [%s] Parsing Gene Data from %s" % (startTime.strftime("%c"),gtfFilePath))
        
        #parse GTF file
        if gtfFilePath.endswith(".gz"):
            gtfFile = gzip.open(gtfFilePath)
        else:
            gtfFile = open(gtfFilePath)
        
        #parse GTF file
        self.parseGtf(gtfFile)
        self.assembleTranscriptome()
        
        duration = Helper.getTime() -startTime
        Helper.info(" Finished parsing in %s" % (str(duration)))
        
        #delete unneccesarry variables
        del self.uniqGeneSet
        del self.uniqGene_to_source
        del self.uniqGene_to_names
        del self.uniqGene_to_transcriptIds
        

        #
        #   transcriptId
        #
        del self.transcriptId_to_names 
        del self.transcriptId_to_protId
        del self.transcriptId_to_startCodons
        del self.transcriptId_to_stopCodons
        del self.transcriptId_to_exons
        del self.transcriptId_to_cds
        del self.transcriptId_to_codingFrames

        del self.transcriptIds
    
    """
    Finds the gene wich are overlapping the given region
    """
    def findOverlappingGenes(self, chromosome, start, stop):
       
        if len(self.genesByChromosome):
            raise Exception("GeneByChromosome Dictionary is empty")
        
        overlappingGenes = []
        
        try:
            for gene in self.genesByChromosome[chromosome]:
                if gene.start < start and gene.end > stop:
                    overlappingGenes.append(gene)
        except KeyError:
            Helper.warning("chromosome %s not found in geneByChromosome" % chromosome)
        return overlappingGenes
    
    """
    Returns a dictionary with chromosomes as key and all the genes on the chromosome as values
    {"1":[Gene1,Gene2....]}  
    """
    def getGenesByChromosome(self):
        genesByChr = defaultdict(list)
        if len(self.geneList):
            raise Exception("Gene List is empty")
        else:
            for gene in self.geneList:
                genesByChr[gene.chromosome].append(gene)
            return genesByChr    
        
        