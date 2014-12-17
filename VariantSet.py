'''
Created on 05.06.2014

@author: david
'''
from Helper import Helper
from collections import defaultdict
import os
import operator
from copy import copy


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
        


class VariantSet(object):
    '''
    handles a vcfFile and stores the results internally as a Dictionary with Tuple of (chromosome,pos,ref,alt) as keys and a the VariantObject as value
    '''
    def __init__(self,vcfFile):
        self.variantDict = self.parseVcf(vcfFile)
        #return self.parseVcfFile_variantSetByChromosome(vcfFile)
    
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
        
        values = map(lambda x: x.strip(), info.split(";")[:-1])
        
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
                    value = value.split(",")
                if name == "GI":
                    a=[]
                    for anno in value.split(","):
                        gene,segments=anno.split(":")
                        a.append((gene,set(segments.split("|"))))
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
    
    def getVariantSetByChromosome(self):
        '''
        returns the variants as Dictionary with chromosome as key and a list of VariantObjects as values
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
        Helper.info(" [%s] Parsing Variant Data from %s" % (startTime.strftime("%c"),vcfFile))
        
        #check correct Type
        if type(vcfFile) == str:
            if os.path.getsize(vcfFile) == 0: #getsize raises OSError if file is not existing
                raise IOError("%s File is empty" % vcfFile)
            vcfFile = open(vcfFile,"r")
        elif type(vcfFile) != file:
            raise TypeError("Invalid type in 'parseVcfFile' (need string or file, %s found)" % type(vcfFile)) 
            
        variantDict = {}
        for v in self.iterator(vcfFile):
            variantDict[(v.chromosome,v.position,v.ref,v.alt)]=v
            #variantDict[(v.chromosome,v.position)]=v
        
        Helper.printTimeDiff(startTime)
        return variantDict
    
    def printVariantDict(self,outfile):
        '''
        print the variants from the dictionary to the outfile if defined
        '''
        if type(outfile) == str:
            try:
                outfile=open(outfile,"w")
            except IOError:
                Helper.warning("Could not open %s to write Variant" % outfile )
        if type(outfile) != file:   
            raise AttributeError("Invalid outfile type in 'printVariantDict' (need string or file, %s found)" % type(outfile))
            
        outfile.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "\n"]))
        for v in self.variantDict.values():
            attributeString=""
            for key in v.attributes.keys():
                if key =="GI":
                    a=""
                    for anno in v.attributes["GI"]:
                        gene,segment = anno
                        if gene == "-":
                            a += gene+":"+"|".join(segment)  
                        else: 
                            gene.geneId+":"+"|".join(segment)+","  
                    attributeString+=key+"="+a+";"
                    continue
                attributeString+= key+"="+str(v.attributes[key])+";"
            outfile.write("\t".join([v.chromosome,str(v.position),v.id,v.ref,v.alt,str(v.qual),v.filter, attributeString+"\n"]))    

    def printGeneList(self,outfile,printSummary=True):
        '''print List of genes with all the variants
        Gene-Variation-File
        "Gene_ID","SEGMENT","#CHROM","GENE_START","GENE_STOP","VAR_POS","REF","ALT","QUAL","BaseCount(A,C,T,G)"
        
        Gene Summary File
        "Gene_ID",#3'UTR,#5'UTR,#EXON,'INTRON,#TOTAL
        
        '''
        sumDict={}
        
        if type(outfile) == str:
            sumFile=outfile[:outfile.rfind(".")]+".summary"
            try:
                outfile=open(outfile,"w")
                
            except IOError:
                Helper.warning("Could not open %s to write Variant" % outfile )
        if type(outfile) != file:   
            raise AttributeError("Invalid outfile type in 'printVariantDict' (need string or file, %s found)" % type(outfile))
        
        sumFile=open(outfile.name[:outfile.name.rfind(".")]+".summary","w")
        outfile.write("\t".join(["#Gene_ID","SEGMENT","#CHROM","GENE_START","GENE_STOP","VAR_POS","REF","ALT","QUAL","BaseCount(A,C,T,G)"]))
        for v in self.variantDict.values():
            anno = v.attributes["GI"]
            for a in anno:
                #TODO: change return type of annotate postion to gene-object
                gene,segments = a
                if gene == "-":
                    outfile.write("\t".join(["-", "-",",".join(segments),v.chromosome,"-","-",v.id,str(v.position),v.ref,v.alt,str(v.qual),",".join(v.attributes["BaseCounts"]),"\n"]))
                else:
                    outfile.write("\t".join([gene.geneId, gene.names[0],",".join(segments),v.chromosome,str(gene.start),str(gene.end),v.id,str(v.position),v.ref,v.alt,str(v.qual),",".join(v.attributes["BaseCounts"]),"\n"]))
                
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
            sumFile.write("\t".join(["#Gene_ID","#3'UTR","#5'UTR","#EXON","INTRON","#TOTAL","\n"]))
            for gene in sumDict.keys():
                numbers=map(str,sumDict[gene])
                if gene=="-":
                    sumFile.write("\t".join(["intergenic","-"]+["-","-","-","-",numbers[4]]+["\n"]))
                else:
                    sumFile.write("\t".join([gene.geneId,gene.names[0]]+numbers+["\n"]))
                
    def getVariantTuble(self,line):
        '''
        returns a tuple of (chromosome, position, alt, ref) from a line of a vcfFile
        '''
        line=line.split("\t")
        try:
            tuple = (line[0],int(line[1]),line[3],line[4])
            #tuple = (line[0],int(line[1]))
        except IndexError:
            raise ValueError("Error in line '%s'" % " ".join(line))
        return tuple
    
    def deleteOverlappsFromVcf(self,variants):
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
        Helper.info(" [%s] Delete overlapps from %s" % (startTime.strftime("%c"),variantsB.name))

        for line in variantsB:
            if line.startswith("#"):
                continue
            varTuple = self.getVariantTuble(line)
            if varTuple in variantSetA:
                #A.discard(varTuple)
                variantSetA.remove(varTuple)
                del self.variantDict[varTuple]
        
        #calculate duration 
        Helper.printTimeDiff(startTime)
    
    def getOverlappsFromBed(self,bedFile,getNonOverlapps=False):
        startTime=Helper.getTime()
        
        
        if type(bedFile) == str:
            bedFile = open(bedFile)
        elif type(bedFile) != file:
            raise TypeError("bedFile has wrong type, need str or file, %s found" % type(bedFile))
        
        startTime=Helper.getTime()
        Helper.info(("Delete overlaps with %s" % bedFile.name))
        
        variantsByChromosome = self.getVariantSetByChromosome() 
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
                     
        if getNonOverlapps:
            overlapps = set(self.variantDict.keys()) - overlapps #delete all accept the ones which are overlapping
        
        newSet={}
        for variantTuple in overlapps:
            #del self.variantDict[variantTuple]
            newSet[variantTuple]=self.variantDict[variantTuple]
        
        Helper.printTimeDiff(startTime)
        return newSet

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
        Helper.info(" [%s] Annotating Variants" % (startTime.strftime("%c")))
        for v in self.variantDict.values():
            anno = genome.annotatePosition(v.chromosome,v.position) #[(gene1,segment1;segment2;..)..]
            GI=[]
            for a in anno:
                GI.append(a)
            v.attributes["GI"]=GI
        
        Helper.printTimeDiff(startTime)
        #TODO: Test this stupid shit
            
