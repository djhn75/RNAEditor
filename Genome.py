'''
Created on 11.06.2014

@author: david
'''
'''
classdocs
'''

#from __builtin__ import SyntaxError
from array import array
from collections import defaultdict
import gzip
#from itertools import izip
from operator import attrgetter

from Gene import Gene
from Helper import Helper
from Transcript import Transcript
import gtfHandler


class Genome(object):
    
    def __init__(self,gtfFile,logFile=None,textField=0):
        '''
        Constructor
        '''
        
        self.logFile=logFile
        self.textField =textField
        
        startTime = Helper.getTime()
        Helper.info(" [%s] Assembling Genome from %s" % (startTime.strftime("%c"),gtfFile),self.logFile,self.textField)
        
        
        #Dict / List for all Genes
        self.geneList = []
        self.genesByChromosome = defaultdict(list)
        
        # list of all features
        self.featureTypes = set() #set of all possible featuresTypes
        self.geneTypes  = set() #set of all possible geneTypes (sources)


        #   gene_id -> (gene_name,contig,strand)
        self.uniqGeneSet = set()
        self.uniqGene_to_source = dict()
        self.uniqGene_to_names          = defaultdict(set)
        self.uniqGene_to_transcriptIds       = defaultdict(set)
        
        #   transcriptId
        self.transcriptId_to_names        = defaultdict(set)
        self.transcriptId_to_protId      = dict()
        self.transcriptId_to_startCodons = defaultdict(tuple)
        self.transcriptId_to_stopCodons  = defaultdict(tuple)
        self.transcriptId_to_exons        = defaultdict(set)
        self.transcriptId_to_cds = defaultdict(set)
        self.transcriptId_to_codingFrames= defaultdict(set)

        self.transcriptIds        = defaultdict(set)
        
        #this fills also the dictionary self.geneByChromosome
        self.createTranscriptomeFromFile(gtfFile)
        self.genesByChromosome = self.getGenesByChromosome()
        
        Helper.printTimeDiff(startTime,self.logFile,self.textField)
        
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
            raise SyntaxError("In line %s with Gene %s! Check your gtf file" % lineNumber,f.geneId)
        
    def assembleTranscriptome(self):
        '''
        Loop over uniqueGeneSet, in which the ENSG-IDs are saved, 
        and assemble all the transcripts and exons for this gene and save it as a Gene object.
        This gene obeject ist then added to the geneList of this Genome object        
        '''
        #transcriptsByType = defaultdict(list)
        
        #construct Genes
        for uniqGene in self.uniqGeneSet:
            geneId, chromosome, strand = uniqGene
            
            geneNames = list(self.uniqGene_to_names[uniqGene])
            geneType = self.uniqGene_to_source[uniqGene]
            geneExons = set()
            geneCds = set()
            
            #get exons from all transcripts
            for transcriptId in self.uniqGene_to_transcriptIds[uniqGene]:
                geneExons |= self.transcriptId_to_exons[transcriptId] # add new exon tuples to geneExons set
                geneCds  |= self.transcriptId_to_cds[transcriptId]
            
            #usually gtf Files are sorted, but this can't be assumed    
            geneExons = sorted(geneExons, reverse = not strand)    
            geneCds = sorted(geneCds, reverse = not strand)
            
            gene = Gene(geneId, chromosome, strand, geneType, geneNames, geneExons, geneCds)
            
            geneExons = dict(zip(geneExons,xrange(1000000)))
            geneCds = dict(zip(geneCds,xrange(1000000))) 
            
            
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
            
    def createTranscriptomeFromFile(self,gtfFilePath):
        """
        Construct Genome from GTF File
        Saves all the information in dictionarys
        
        This function calls internally:
            -parseGTF
            -assembleTranscriptome
        Returns a list with all the genes for each chromosome in a dictionary
        return: genesByChromosome
        """
        #startTime = Helper.getTime()
        #Helper.info(" [%s] Parsing Gene Data from %s" % (startTime.strftime("%c"),gtfFilePath))
        
        #check fileType of Genome File
        if gtfFilePath.endswith(".gz"):
            gtfFile = gzip.open(gtfFilePath)
        else:
            gtfFile = open(gtfFilePath)
        
        #parse GTF file
        self.parseGtf(gtfFile)
        self.assembleTranscriptome()
        
        #duration = Helper.getTime() -startTime
        #Helper.info(" Finished parsing in %s" % (str(duration)))
        
        del self.featureTypes
        del self.geneTypes
        
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
     
    
    def getGenesByChromosome(self):
        """
        Returns a dictionary with chromosomes as key and all the genes on the chromosome as values
        The genes are also sorted
        {"1":[Gene1,Gene2....]}  
        """
        
        genesByChr = defaultdict(list)
        if len(self.geneList) == 0:
            raise Exception("Gene List is empty")
        else:
            for gene in self.geneList:
                genesByChr[gene.chromosome].append(gene)
            
            #sort genes by by start and end
            for key in genesByChr.keys():
                genesByChr[key] = sorted(genesByChr[key], key=attrgetter('start','end'))
        
        return genesByChr    
        
    def getGenesByGeneID(self):
        """
        Returns a dictionary with geneID (ENSG000001) as key and the gene object as value
        {"ENSG000001":GeneObject;"ENSG000002":GeneObject2}  
        """
        genesByGeneID = defaultdict()
        for gene in self.geneList:
            genesByGeneID[gene.geneId]=gene
        return genesByGeneID
    
    def annotateRegion(self,chromosome,start,stop):
        """
            returns information for the given region like (3'UTR,Exon,Intron,5'UTR)
        
        Gene ......|---------gene1-----------|...................|-----gene2------|...
        pos1 ...|----pos1----|......................................
        pos2 ................|---pos2----|..........................
        pos3 ..........................|-------pos3----|............
        pos4 ..|------------------------pos4---------------|........
        pos5 |pos5|..................................................
        pos6...|----------------------------------------pos6-------------------------|.....
        """
       
       
       
        for gene in self.genesByChromosome[chromosome]:
            segment=set()
            
            if start < gene.start < stop and stop < gene.end:
                #case 1
                pass 
        
        #TODO: write this function like annotate position
        
    def annotatePosition(self,chromosome, position):
        '''
        returns the gene and information for the given position like (3'UTR,Exon,Intron,5'UTR)
        :param chromosome: String
        :param position: Int
        :return list of Tuples List[(gene,segment1;segment2..)...]
        '''
        result=[]

        #Loop over the genes of the chromosome
        for gene in self.genesByChromosome[chromosome]:
            #check if the position is in the current gene
            
            #geneName=None
            segment = set()
            
            if gene.start < position < gene.end:
                #geneName=gene.names[0]    
                
                if len(gene.codingExons)>0:
                    #check if position is in front of first coding exon
                    if (gene.strand and position < gene.codingExons[0][0]) or (not gene.strand and position > gene.codingExons[0][1]):
                        for exon in gene.exons:
                            if exon[0]<position<exon[1]:
                                segment.add("5'UTR")
                                #continue
                    #check if position is behind the last coding exon
                    elif (gene.strand and position > gene.codingExons[-1][1]) or not gene.strand and position < gene.codingExons[-1][0] :
                        for exon in gene.exons:
                            if exon[0]<position<exon[1]:
                                segment.add("3'UTR")
                                #continue
                    
                    for cds in gene.codingExons:
                        if cds[0] < position < cds[1]:
                            segment.add("coding-exon")
                            #continue
                else:
                    for exon in gene.exons:
                        if exon[0] < position < exon[1]:
                            segment.add("noncoding-exon")
                            #continue
                
                if len(segment)==0:
                    segment.add("intron")
                    
                result.append((gene,tuple(segment)))
                #if len(segment)>=2:
                #   print(chromosome,gene.geneId,int(position),segment)
                    
        if result == []:
            result.append(("-",tuple(["intergenic"]))) 
        #return geneName, segment    
        return result
        
  