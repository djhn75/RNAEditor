'''
Created on 05.06.2014

@author: david
'''
from Helper import Helper
from collections import defaultdict, OrderedDict
import os
import operator
from copy import copy
#from exceptions import KeyError
from Genome import Genome
import collections
import numpy as np
from random import shuffle
import sys
from pysam import Samfile
from Gene import Gene

class Variant:
    '''
    reflects a Variant
    '''

    def __init__(self, chromosome, position, id, ref, alt, qual, filter, info):
        self.chromosome = chromosome
        self.position = position
        self.id = id
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filter
        self.attributes = info
        
    def __iter__(self):
        return self.var

class VariantSet(object):
    '''
    handles a vcfFile and stores the results internally as a Dictionary with Tuple of (chromosome,pos,ref,alt) as keys and a the VariantObject as value
    '''
    def __init__(self,vcfFile,logFile=None,textField=0):
        self.logFile=logFile
        self.textField=textField
        self.variantDict = self.parseVcf(vcfFile)
        #return self.parseVcfFile_variantSetByChromosome(vcfFile)
    
    def __iter__(self):
        return iter((self.variantDict.values()))
    
    def __add__(self,other):
        newVariantSet = copy(self)
        newDict = {}
        newDict.update(self.variantDict)
        newDict.update(other.variantDict)
        newVariantSet.variantDict=newDict
        return newVariantSet
    
    def readline(self,line):
        '''
        process one line of the vcf file
        <chromosome> <position> <identifier> <reference_base> <alternative_base> <quality> <filter> {attributes}
        '''
        vcfList = line.rstrip().split("\t")
        
        try:
            vcfList[1] = int(vcfList[1]) #position of SNP
            vcfList[5] = float(vcfList[5]) if vcfList[5] !="." else 0.0

        except ValueError:
            raise ValueError("Error in line '%s'" % " ".join(line))

        #parse info
        info = vcfList[7]
        #trim comments
        info=info[:info.find("#")].rstrip()
        
        values = map(lambda x: x.strip(), info.split(";"))
        #[:-1])
        
        attributes={}
        for info in values:
            info = map( lambda x: x.strip(), info.split("="))
            if len(info)>1:
                name, value=info[0], info[1]
                try:
                    value=float(value)
                    value=int(value)
                except ValueError:
                        pass
                except TypeError:
                        pass    
               
                if name == "BaseCounts":
                    value=value.replace("'","")
                    value=value.replace("[","")
                    value=value.replace("]","")
                    value = value.split(",")
                if name == "GI":
                    a=[]
                    for anno in value.split(","):
                        #TODO: Delete the next line later, this is because of a tailing comma which was removed
                        if anno=="": continue
                        gene,segments=anno.split(":")
                        a.append((gene,set(segments.split("|"))))
                    value=a
                attributes[name]=value
        
        vcfList[7]=attributes    
        return vcfList
    
    def checkVariantType(self,variants):
            '''
            Checks if the type of the argument is a str or a file
            returns a dictionary of all the variants
            '''
            if type(variants) == dict:
                return variants
            elif type(variants) == file or type(variants) == str:
                variants = self.parseVcf(variants)
                return variants
                
            else: 
                raise TypeError("variants has wrong type, need variantDict, str or file, %s found" % type(variants))
                  
    def iterator(self,infile):
        while True:
            line = infile.readline()
            if not line: raise StopIteration
            if line.startswith("#"): continue #skip comments
            vcfList=self.readline(line)
            variant = Variant(vcfList[0],vcfList[1],vcfList[2],vcfList[3],vcfList[4],vcfList[5],vcfList[6],vcfList[7])
            yield variant
    
    def getVarPosListByChromosome(self):
        '''
        return: all the variant positions by chromosome
        {"1":[4,6,8,45,67],"2":[6,9,67,69].....}
        This is only needed for the cluster algorithm later on
        '''
        varPosList=defaultdict(list)
        for v in self.variantDict.values():
            varPosList[v.chromosome].append(v.position)
            
        #make numpy array out of the lists
        for chromosome in varPosList.keys():
            varPosList[chromosome]=np.asarray(varPosList[chromosome])
        return varPosList
    
    def getVariantListByChromosome(self):
        '''
        @return: variants as Dictionary with chromosome as key and a list of VariantObjects as values
        {"1":[VariantObject1,VariantObject2....],"2":[VariantObject1,VariantObject2....]}
        '''
        variantsByChromosome = defaultdict(list)
        for v in self.variantDict.values():
            variantsByChromosome[v.chromosome].append(v)
        
        #Helper.printTimeDiff(startTime)
        return variantsByChromosome
    
    def parseVcf(self,vcfFile):
        '''
        Imports a given Variant File and returns the variants as Dictionary with Tuple of (chromosome,pos,ref,alt) as key and a the VariantObject as value
        {(1,45435,"A","G"):VariantObject1,(1,45435,"A","G"):VariantObject1,.....}
        
        '''
        startTime = Helper.getTime()
        Helper.info(" [%s] Parsing Variant Data from %s" % (startTime.strftime("%c"),vcfFile),self.logFile,self.textField)
        
        #check correct Type
        if type(vcfFile) == str:
            if os.path.getsize(vcfFile) == 0: #getsize raises OSError if file is not existing
                raise IOError("%s File is empty" % vcfFile)
            vcfFile = open(vcfFile,"r")
        elif type(vcfFile) != file:
            raise TypeError("Invalid type in 'parseVcfFile' (need string or file, %s found)" % type(vcfFile)) 
            
        variantDict = OrderedDict()
        for v in self.iterator(vcfFile):
            variantDict[(v.chromosome,v.position,v.ref,v.alt)]=v
            #variantDict[(v.chromosome,v.position)]=v
        
        Helper.printTimeDiff(startTime,self.logFile,self.textField)
        return variantDict
    
    def printVariantDict(self,outfile):
        '''
        print the variants from the dictionary to the outfile if defined
        '''
        if type(outfile) == str:
            try:
                outfile=open(outfile,"w")
            except IOError:
                Helper.warning("Could not open %s to write Variant" % outfile ,self.logFile,self.textField)
        if type(outfile) != file:   
            raise AttributeError("Invalid outfile type in 'printVariantDict' (need string or file, %s found)" % type(outfile))
        
        startTime=Helper.getTime()
        Helper.info("[%s] Print Variants to %s" %  (startTime.strftime("%c"),outfile.name),self.logFile,self.textField)
            
        outfile.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "\n"]))
        for v in self.variantDict.values():
            attributeString=""
            for key in v.attributes.keys():
                if key=="BaseCounts":
                    attributeString+= "BaseCounts=" + ",".join(v.attributes["BaseCounts"]) + ";"
                    continue                
                elif key =="GI":
                    a=""
                    for anno in v.attributes["GI"]:
                        gene,segment = anno
                        if gene == "-":
                            a += gene+":"+"|".join(segment)  
                        else:
                            if type(gene)==str: #when variantDict was not annotated yet
                                a+=gene +":"+"|".join(segment)+","
                            else:     
                                a+=gene.names[0]+":"+"|".join(segment)+","
                            
                    attributeString+=key+"="+a[:-1]+";"
                    continue
                attributeString+= key+"="+str(v.attributes[key])+";"
            outfile.write("\t".join([v.chromosome,str(v.position),v.id,v.ref,v.alt,str(v.qual),v.filter, attributeString+"\n"]))    


    def topGenes(self,sumDict, fileName,number=20,value=4):
        if number > len(sumDict):
            if len(sumDict)<1:
                Helper.warning("no edited genes found", self.logFile, self.textField)
                return
            Helper.warning("The number of top genes you wanted is bigger than the number of edited genes", self.logFile, self.textField)
            number=len(sumDict)
        if value > 4:
            Helper.error("sumDict only hold four values", self.logFile, self.textField)
        
        
        counts=collections.OrderedDict(sorted(sumDict.items(), key=lambda t: t[1][value],reverse=True)[:number])
        barNameTuple=()
        valueMatrix=[[]]
        for array in counts.values():
            valueMatrix[0].append(array[value])
        for gene in counts.keys():
            barNameTuple+=(gene.names[0],)

        if value==0:
            barName="3'-UTR"
        elif value==1:
            barName="5'-UTR"
        elif value==2:
            barName="Exonic"
        elif value==3:
            barName="Intronic"
        elif value==4:
            barName="Total"
        
        yLim=max(max(i) for i in valueMatrix)+1
        Helper.createBarplot(valueMatrix, fileName, barNameTuple, [barName], width=0.35, title="Highly Edited Genes",yLim=yLim,barText=False,yText="Editing Counts")

    def printGeneList(self,genome,outfile,printSummary=True):
        '''
        print List of genes with all the variants
        Gene-Variation-File
        "Gene_ID","gene_Name","SEGMENT","#CHROM","GENE_START","GENE_STOP","VAR_POS","REF","ALT","QUAL","BaseCount(A,C,T,G)"
        
        Gene Summary File
        "Gene_ID",Gene_Name,#3'UTR,#5'UTR,#EXON,'INTRON,#TOTAL
        :param genome:  object of class Genome
        :param outfile: 
        :param printSummary: boolean wether to print summary-file
        '''

        sumDict={}
        
        if type(genome) != Genome:
            raise AttributeError("Type of genome is %s, but has to be an object of Genome" % type(genome))
        
        if type(outfile) == str:
            try:
                outfile=open(outfile,"w")
                
            except IOError:
                Helper.warning("Could not open %s to write Variant" % outfile ,self.logFile,self.textField)
        if type(outfile) != file:   
            raise AttributeError("Invalid outfile type in 'printVariantDict' (need string or file, %s found)" % type(outfile))
        
        startTime=Helper.getTime()
        Helper.info("[%s] Print Genes and Variants to %s" %  (startTime.strftime("%c"),outfile.name),self.logFile,self.textField)
        
        sumFile=open(outfile.name[:outfile.name.rfind(".")]+".summary","w")        
        
        outfile.write("\t".join(["#Gene_ID","Name","SEGMENT","#CHROM","GENE_START","GENE_STOP","VAR_ID","VAR_POS",
                                 "REF","ALT","QUAL","#A","#C","#G","#T","Reads_Total","Edited_Reads","Editing_Ratio","\n"]))
        
        for v in self.variantDict.values():
            anno = v.attributes["GI"]
            for a in anno:
                gene,segments = a
                totalReads=str(int(sum(map(int,v.attributes["BaseCounts"]))))
                if v.ref =="A" and v.alt == "G":
                    editedReads=str(v.attributes["BaseCounts"][2])
                    ratio=str(round(float(editedReads)/float(totalReads),2))
                elif (v.ref=="T" and v.alt=="C"):
                    editedReads=str(v.attributes["BaseCounts"][1])
                    ratio=str(round(float(editedReads)/float(totalReads),2))
                else:
                    editedReads="0"
                    ratio="0"
                    
                
                if gene == "-":
                    out=["-", "-",",".join(segments),v.chromosome,"-","-",v.id,str(v.position),
                                             v.ref,v.alt,str(v.qual),"\t".join(v.attributes["BaseCounts"]),totalReads,editedReads,ratio,"\n"]
                    outfile.write("\t".join(out))
                else:
                    out=[gene.geneId, gene.names[0],",".join(segments),v.chromosome,str(gene.start),str(gene.end),v.id,str(v.position),
                                             v.ref,v.alt,str(v.qual),"\t".join(v.attributes["BaseCounts"]),totalReads,editedReads,ratio,"\n"]
                    outfile.write("\t".join(out))
                
                #count variations per gene
                if gene not in sumDict:
                    sumDict[gene]= [0,0,0,0,0]
                
                for seg in segments:
                    if seg == "3'UTR":
                        sumDict[gene][0]+=1
                    elif seg == "5'UTR":
                        sumDict[gene][1]+=1
                    elif seg in ("coding-exon","noncoding-exon"):
                        sumDict[gene][2]+=1
                    elif seg == "intron":
                        sumDict[gene][3]+=1
                    sumDict[gene][4]+=1
        
        

                     
        #print number of variants per gene
        if printSummary:

            
            sumDictGeneIds=set()
            sumFile.write("\t".join(["#Gene_ID","Name","#3'UTR","#5'UTR","#EXON","INTRON","#TOTAL","\n"]))
            for gene in sumDict.keys():
                numbers=map(str,sumDict[gene])
                if gene=="-":
                    sumFile.write("\t".join(["intergenic","-"]+["-","-","-","-",numbers[4]]+["\n"]))
                else:
                    sumFile.write("\t".join([gene.geneId,gene.names[0]]+numbers+["\n"]))
                    sumDictGeneIds.add(gene.geneId)        
            #print non effected Genes
            #this was added to have the whole set og genes in the summary file
            #so that it is easier to compare results in Excel
            genesByGeneId=genome.getGenesByGeneID()
            a=set(genesByGeneId.keys())
            b=sumDictGeneIds
            nonEffectedGenes = a-b
            for geneId in nonEffectedGenes:
                gene=genesByGeneId[geneId]
                sumFile.write("\t".join([gene.geneId,gene.names[0]]+["0","0","0","0","0",]+["\n"]))
            
            ################################################################
            ############    Draw Barplots with high edited Genes ###########
            ################################################################
            '''
            outdir = outfile.name[:outfile.name.rfind("/")+1]
            tmp=outfile.name[outfile.name.rfind("/")+1:]
            sampleName=tmp[:tmp.find(".")
                           ]
            fileName=outdir+"html/"+sampleName+".editedGenes(3UTR).png"
            self.topGenes(sumDict,fileName, 20, 0)
            
            fileName=outdir+"html/"+sampleName+".editedGenes(5UTR).png"
            self.topGenes(sumDict,fileName, 20, 1)
            
            fileName=outdir+"html/"+sampleName+".editedGenes(Exon).png"
            self.topGenes(sumDict,fileName, 20, 2)
            
            fileName=outdir+"html/"+sampleName+".editedGenes(Intron).png"
            self.topGenes(sumDict,fileName, 20, 3)
            
            fileName=outdir+"html/"+sampleName+".editedGenes(Total).png"
            
            if "-" in sumDict.keys():
                del sumDict["-"] #delete intergenics, because we only we only want to show highly edited Genes!!!
            self.topGenes(sumDict,fileName, 20, 4)
            '''
                            
    def printClusters(self, outFile):
        
        if type(outFile) == str:
            try:
                outFile=open(outFile,"w")
                
            except IOError:
                Helper.warning("Could not open %s to write Variant" % outFile ,self.logFile,self.textField)
        if type(outFile) != file:   
            raise AttributeError("Invalid outfile type in 'printVariantDict' (need string or file, %s found)" % type(outFile))
        
        startTime=Helper.getTime()
        Helper.info("[%s] Print Clusters to %s" %  (startTime.strftime("%c"),outFile.name),self.logFile,self.textField)
        
        
        outFile.write("\t".join(["#Chr","Start","Stop","IslandID","GeneID","Gene Symbol","Cluster Length","Number of Editing_sites","Editing_rate","\n"]))
        
        for cluster in self.clusterDict.keys():
            end = max(v.position for v in self.clusterDict[cluster])
            start = min(v.position for v in self.clusterDict[cluster])
            
            length = end - start
            editingRate=float(len(self.clusterDict[cluster]))/float(length)
            geneIdSet=set()
            geneNameSet=set()
            for v in self.clusterDict[cluster]:
                try: 
                    gene = v.attributes['GI'][0][0]
                    if type(gene) == Gene:
                        geneIdSet.add(gene.geneId)
                        geneNameSet |= set(gene.names)
                        #geneList.append(v.attributes['GI'][0][0])
                    else:
                        geneIdSet.add("Intergenic")
                        geneNameSet.add("Intergenic")
                except KeyError:
                    geneIdSet.add("N/A") #when variant has no attribute GI
            
            outFile.write("\t".join([v.chromosome,str(start),str(end),"Island"+str(cluster), #Chr","Start","Stop","Cluster Name",
                                     ",".join(map(str,geneIdSet)),",".join(map(str,geneNameSet)), #"GeneID","Gene Symbol"
                                     str(length),str(len(self.clusterDict[cluster])),'%1.2f'%float(editingRate),"\n"]))
            
    def getVariantTuble(self,line):
        '''
        returns a tuple of (chromosome, position, alt, ref) from a line of a vcfFile
        '''
        line=line.split("\t")
        try:
            for alt in line[4].split(","):
                tuple = (line[0],int(line[1]),line[3],alt)
                yield tuple
        except IndexError:
            raise ValueError("Error in line '%s'" % " ".join(line))
    
    def getVariantByGene(self):
        '''
        Returns a dictionary with geneId as key and all the variants on the gene as values
        The genes are also sorted
        {"1":[Gene1,Gene2....]}  
        '''
        variantByGene=defaultdict(set)
        
        try:
            for v in self.variantDict.values():
                for anno in v.attributes["GI"]:
                    gene,segment = anno
                    variantByGene[gene].add(v)
        except KeyError:
            raise KeyError("Variant has no attribute GI. Try to run 'annotateVariantDict' before to get GeneInfo")
        return variantByGene
    
    def deleteOverlapsFromVcf(self,variants):
        '''
        delete the variants from 'variantsA' which also are in 'variantsB'
        '''

        variantSetA = set(self.variantDict.keys())
        
        #detrmine type of variantB
        if type(variants) == str:
            variantsB = open(variants)
        elif type(variants) != file:
            raise TypeError("variantB has wrong type, need str or file, %s found" % type(variantsB))
        #TODO: variants could also be another object of VariantsSet
        
        #get Start time
        startTime = Helper.getTime()
        Helper.info(" [%s] Delete overlapps from %s" % (startTime.strftime("%c"),variantsB.name),self.logFile,self.textField)

        for line in variantsB:
            if line.startswith("#"):
                continue
            for varTuple in self.getVariantTuble(line):
                if varTuple in variantSetA:
                #A.discard(varTuple)
                    variantSetA.remove(varTuple)
                    del self.variantDict[varTuple]
        
        #calculate duration 
        Helper.printTimeDiff(startTime,self.logFile,self.textField)
    
    def getOverlapsFromBed(self,bedFile,getNonOverlaps=False):
        '''
        returns overlaps from bed file features
        :param bedFile: as string or file
        :param getNonOverlaps: boolean
        :return new variantSet of overlaps 
        '''
        
        if type(bedFile) == str:
            bedFile = open(bedFile)
        elif type(bedFile) != file:
            raise TypeError("bedFile has wrong type, need str or file, %s found" % type(bedFile))
        
        startTime=Helper.getTime()
        Helper.info("[%s] Delete overlaps from %s" %  (startTime.strftime("%c"),bedFile.name) ,self.logFile,self.textField)
        
        variantsByChromosome = self.getVariantListByChromosome() 
        overlapps = set()
        for line in bedFile:
            try:
                sl = line.split("\t") 
                #if "\t" in line else line.split(" ")
                chromosome,start,stop = sl[:3]
                start,stop=(int(start),int(stop))
            except ValueError:
                raise ValueError("Error in line '%s'" % line)
            
            for v in variantsByChromosome[chromosome]:
                if start < v.position < stop:
                    overlapps.add((v.chromosome,v.position,v.ref,v.alt))
                     
        if getNonOverlaps:
            overlapps = set(self.variantDict.keys()) - overlapps #delete all accept the ones which are overlapping
        
        newSet={}
        for variantTuple in overlapps:
            #del self.variantDict[variantTuple]
            newSet[variantTuple]=self.variantDict[variantTuple]
        
        Helper.printTimeDiff(startTime, self.logFile,self.textField)
        return newSet
    
    def splitByBed(self,bedFile):
        '''
        returns overlaps and nonOverlaps from bed file features
        :param bedFile: as string or file
        :param getNonOverlaps: boolean
        '''
        
        if type(bedFile) == str:
            bedFile = open(bedFile)
        elif type(bedFile) != file:
            raise TypeError("bedFile has wrong type, need str or file, %s found" % type(bedFile))
        
        startTime=Helper.getTime()
        Helper.info("[%s] Split Variants by Bed File %s" %  (startTime.strftime("%c"),bedFile.name) ,self.logFile,self.textField)
        
        variantsByChromosome = self.getVariantListByChromosome() 
        overlapSet = set()
        i=0
        for line in bedFile:
            
            try:
                sl = line.split("\t") 
                #if "\t" in line else line.split(" ")
                chromosome,start,stop = sl[:3]
                start,stop=(int(start),int(stop))
            except ValueError:
                raise ValueError("Error in line '%s'" % line)
            
            for v in variantsByChromosome[chromosome]:
                if start < v.position < stop:
                    overlapSet.add((v.chromosome,v.position,v.ref,v.alt))
            i+=1
            if i %100000==0:
                Helper.status("%s Bed Feautes parsed" % i, self.logFile,self.textField,"grey")
        
        
        Helper.info("finished parsing Bed file", self.logFile,self.textField)
        Helper.printTimeDiff(startTime, self.logFile,self.textField)
               
        #nonOverlapSet = set(self.variantDict.keys()) - overlapSet #delete all accept the ones which are overlapping
        
        
        overlaps = {key: self.variantDict[key] for key in self.variantDict if key in overlapSet}
        
        Helper.info("finished creating overlaps", self.logFile,self.textField)
        Helper.printTimeDiff(startTime, self.logFile,self.textField)
        
        nonOverlaps = {key: self.variantDict[key] for key in self.variantDict if key not in overlapSet}
        
        """
        overlaps={}
        for variantTuple in overlapSet:
            #del self.variantDict[variantTuple]
            overlaps[variantTuple]=self.variantDict[variantTuple]
        
        nonOverlaps={}
        for variantTuple in nonOverlapSet:
            nonOverlaps[variantTuple]=self.variantDict
        """
        
        Helper.printTimeDiff(startTime, self.logFile,self.textField)
        return overlaps, nonOverlaps

    def sortVariantDict(self,variantDict):
        '''
        Sorts a VariantDictionary by the variant position
        :param variantDict:
        '''
        #if type(variantDict) != list:
        #    raise TypeError("variants has wrong type, need variantDict, %s found" % type(variantDict))
        for key in variantDict.keys():
            variantDict[key] = sorted(variantDict[key], key=operator.attrgetter('position'))

    def annotateVariantDict(self,genome):
        '''
        adds the corresponding Gene and the exact segment wehre the SNP appears
        :param genome: Genome
        '''
        startTime = Helper.getTime()
        Helper.info(" [%s] Annotating Variants" % (startTime.strftime("%c")),self.logFile,self.textField)
        for v in self.variantDict.values():
            anno = genome.annotatePosition(v.chromosome,v.position) #[(gene1,segment1;segment2;..)..]
            GI=[]
            for a in anno:
                GI.append(a)
            v.attributes["GI"]=GI
        
        Helper.printTimeDiff(startTime,self.logFile,self.textField)
            
    def createClusters(self,eps=50,minSamples=5):
        
        islandCounter=0
        eps=int(eps)
        minSamples=int(minSamples)
        variantsByChromosome = self.getVariantListByChromosome()
        self.clusterDict=defaultdict(list)
        for chr in variantsByChromosome.keys():
            posList = [v.position for v in variantsByChromosome[chr]] #position of all variants from that chromosome
            
            labels = self.getLabels(posList,eps,minSamples) #actually doing db clustering
            n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
            
            
            if n_clusters_ > 0:
                #loop over labels and variants
                tmpDict=defaultdict(list)
                for var,label in zip(variantsByChromosome[chr],labels):
                    #clusterdict{1:[var1,var2],2:[var5,var8]}
                    if label==-1:
                        continue
                    
                    tmpDict[label].append(var)
                    
                #set new label for clusterdict, to avoid overwriting
                for label in tmpDict.keys():
                    self.clusterDict[islandCounter]=tmpDict.pop(label)
                    islandCounter+=1
                                    
    def getLabels(self,positionList,eps=10, minSamples=5):
        """Perform DBSCAN clustering from vector array.
    
        Parameters
        ----------
        X: array [int1,int1]
            Array of Samples. 
            In this case it should be the positions of the variations in the genome per chromosome
    
        eps: float, optional
            The maximum distance between two samples for them to be considered
            as in the same neighborhood.
    
        minSamples: int, optional
            The number of samples in a neighborhood for a point to be considered
            as a core point.
    
    
        Returns
        -------
        core_samples: array [n_core_samples]
            Indices of core samples.
    
        labels : array [n_samples]
            Cluster labels for each point.  Noisy samples are given the label -1.
    
        """
        if not eps > 0.0:
            raise ValueError("eps must be positive.")
    
        X = np.asarray(positionList)
        n = X.shape[0] #get number of elements (not sure) 
    
        index_order=range(n)
        shuffle(index_order)
        
        
        distanceMatrix = self.calculate1dDistanceMatrix(X,eps)
    
        # Calculate neighborhood for all samples. This leaves the original point
        # in, which needs to be considered later (i.e. point i is the
        # neighborhood of point i. While True, its useless information)
    
        #distanceMatrix = [np.where(x <= eps)[0] for x in distanceMatrix]
        
        # Initially, all samples are noise.
        labels = -np.ones(n, dtype=np.int)
    
        # A list of all core samples found.
        core_samples = []
    
        # label_num is the label given to the new cluster
        label_num = 0
    
        # Look at all samples and determine if they are core.
        # If they are then build a new cluster from them.
        for index in index_order:
            # Already classified
            if labels[index] != -1:
                continue
    
            # get neighbors from distanceMatrix or ballTree
            index_neighborhood = []
    
            index_neighborhood = distanceMatrix[index]
            
    
            # Too few samples to be core
            if len(index_neighborhood) < minSamples:
                continue
    
            core_samples.append(index)
            labels[index] = label_num
            # candidates for new core samples in the cluster.
            candidates = [index]
    
            while len(candidates) > 0:
                new_candidates = []
                # A candidate is a core point in the current cluster that has
                # not yet been used to expand the current cluster.
                for c in candidates:
                    c_neighborhood = []
                    
                    c_neighborhood = distanceMatrix[c]
                    
                    noise = np.where(labels[c_neighborhood] == -1)[0] #indexes of candidate neigbours which do not belong to a cluster yet
                    noise = c_neighborhood[noise]
                    labels[noise] = label_num
                    for neighbor in noise:
                        n_neighborhood = []
                        
                        n_neighborhood = distanceMatrix[neighbor]
                        
                        # check if its a core point as well
                        if len(n_neighborhood) >= minSamples:
                            # is new core point
                            new_candidates.append(neighbor)
                            core_samples.append(neighbor)
                # Update candidates for next round of cluster expansion.
                candidates = new_candidates
            # Current cluster finished.
            # Next core point found will start a new cluster.
            label_num += 1
        #return core_samples, labels
        
        return labels
    
    def calculate1dDistanceMatrix(self,lst,eps):
        '''
        creates a distance matrix for the given vector
        :param lst: vector of samples       
        :return: np.array(diffMatrix)
        '''
        if not isinstance(lst, (list, tuple, np.ndarray)):
            raise TypeError("Paramer has to be eithe a List or a Tuple found %s" % type(lst))
        if not all(isinstance(item, (int,float)) for item in lst):
            raise TypeError("List should only contain numbers")
        lst = np.asarray(lst)
        diffMatrix=[]

        for l1 in lst:
            diffList=[]
            
            diffList= abs(lst-l1)
            diffList = np.where(diffList<=eps)[0]
            diffMatrix.append(diffList)

        return np.asarray(diffMatrix)

    def deleteNonEditingBases(self):
        startTime=Helper.getTime()
        Helper.info("Delete non Editing Bases (keep only T->C and A->G)",self.logFile,self.textField)
        
        for varTuple in self.variantDict.keys():
            chr,pos,ref,alt = varTuple
            if (ref =="A" and alt == "G") or (ref=="T" and alt=="C"):
                pass
            else:
                del self.variantDict[varTuple]

    def __len__(self):
        return len(self.variantDict)
    
    def removeEdgeMismatches(self,bamFile,minDistance, minBaseQual):
        startTime=Helper.getTime()
        minDistance=int(minDistance)
        counter=0;j=0  
        num_lines = len(self.variantDict)
        Helper.info(" [%s] remove Missmatches from the first %s bp from read edges" % (startTime.strftime("%c"),str(minDistance)),self.logFile,self.textField)
        
        bamFile = Samfile(bamFile, "rb")
        
        for varKey in self.variantDict.keys():
            variant = self.variantDict[varKey]
            
            counter+=1
            if counter%10000==0:
                Helper.status('%s mm parsed ' % counter ,self.logFile, self.textField,"grey")
            
            keepSNP=False
            varPos=variant.position-1
            iter = bamFile.pileup(variant.chromosome, variant.position-1, variant.position)
            #walks up the region wich overlap this position
            for x in iter:
                if x.pos == varPos:
                    for pileupread in x.pileups: #walk through the single reads
                        if not pileupread.is_del and not pileupread.is_refskip:
                            distance=abs(pileupread.alignment.alen-pileupread.query_position) if pileupread.alignment.is_reverse else pileupread.query_position
                            if distance >= minDistance:
                                #check readBase and Base Quality
                                if pileupread.alignment.query_sequence[pileupread.query_position] == variant.alt and pileupread.alignment.query_qualities[pileupread.query_position]>=minBaseQual:
                                #if pileupread.alignment.query_sequence[pileupread.query_position] == variant.alt:
                                    keepSNP=True
                                    
            if keepSNP==False:
                j+=1
                del self.variantDict[varKey]
        
        Helper.status('%s of %svariants were deleted' % (j,num_lines), self.logFile, self.textField,"black") 
        Helper.printTimeDiff(startTime, self.logFile, self.textField)
        bamFile.close()