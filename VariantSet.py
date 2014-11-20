'''
Created on 05.06.2014

@author: david
'''
from Helper import Helper
from collections import defaultdict
import os
import operator


class Variant:
    '''
    reflects a Variant
    '''
    
    def __init__(self, chromosome, position, id, ref, alt, qual, filter, info, baseCounts ):
        self.chromosome = chromosome
        self.position = position
        self.id = id
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filter
        self.info = info
        if type(baseCounts) == list and len(baseCounts) == 4 or baseCounts == None:
            self.baseCounts = baseCounts #[A,C,G,T] number of reads which contain the appropriate base 
        else:
            raise TypeError("Parameter Basecount is malformed: %s " % str(baseCounts))


class VariantSet(object):
    
    def __init__(self,vcfFile):
        return self.parseVcfFile_variantSetByChromosome(vcfFile)
    
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
        
        baseCountPos = info.find("BaseCounts")
        if baseCountPos != -1:
            #extract BaseCounts from list and save it in a list
            vcfList[8] = info[baseCountPos:info.find(";", baseCountPos)].split("=")[1].split(",") 
        else:
            vcfList[8] = None
        return vcfList
    
    def checkVariantType(self,variants):
            '''
            Checks if the type of the argument is a str or a file
            returns a dictionary of all the variants
            '''
            if type(variants) == dict:
                return variants
            elif type(variants) == file or type(variants) == str:
                variants = self.parseVcfFile_variantsDict(variants)
                return variants
                
            else: 
                raise TypeError("variants has wrong type, need variantDict, str or file, %s found" % type(variants))
         
             
    def iterator(self,infile):
    
        while 1:
            line = infile.readline()
            if not line: raise StopIteration
            if line.startswith("#"): continue #skip comments
            vcfList=self.readline(line)
            variant = Variant(vcfList[0],vcfList[1],vcfList[2],vcfList[3],vcfList[4],vcfList[5],vcfList[6],vcfList[7],vcfList[8])
            yield variant
    
    def parseVcfFile_variantSetByChromosome(self,vcfFile):
        '''
        Imports a given Variant File and returns the variants as Dictionary with chromosome as key and a list of VariantObjects as values
        {"1":[VariantObject1,VariantObject2....],"2":[VariantObject1,VariantObject2....]}
        
        '''
        startTime = Helper.getTime()
        Helper.info(" [%s] Parsing Variant Data from %s" % (startTime.strftime("%c"),vcfFile))
        
        #check correct Type
        if type(vcfFile) == str:
            vcfFile = open(vcfFile)
        elif type(vcfFile) != file:
            raise TypeError("Invalid type in 'parseVcfFile' (need string or file, %s found)" % type(vcfFile)) 
            
        variantsByChromosome = defaultdict(list)
        for v in self.iterator(vcfFile):
            variantsByChromosome[v.chromosome].append(v)
        
        Helper.printTimeDiff(startTime)
        return variantsByChromosome
    
    def parseVcfFile_variantsDict(self,vcfFile):
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
    
    
    def printVariantDict(self,variantDict,outfile=None):
        '''
        print the variants from the dictionary to the outfile if defined
        '''
        if type(outfile) == str:
            try:
                outfile=open(outfile,"w")
            except IOError:
                Helper.warning("Could not open %s to write Variant" % outfile )
        elif type(outfile) == file:   
            outfile.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "\n"]))
            for v in variantDict.values():
                
                outfile.write("\t".join([v.chromosome,str(v.position),v.id,v.ref,v.alt,str(v.qual),v.filter, v.info, "\n"]))    
        else:
            raise AttributeError("Invalid outfile type in 'parseVcfFile' (need string or file, %s found)" % type(outfile))
            
    def printVariantSet(self,variantSet,variantDict,outfile=None):
        '''
        prints variants from the set to the outfile if defined
        ''' 
        if outfile != None:
            try:
                outfile=open(outfile,"w")
            except IOError:
                Helper.warning("Could not open %s to write Variant" % outfile )
            
            for key in variantDict.values():
                v = variantDict[key] 
                outfile.write("\t".join([v.chromosome,str(v.position),v.ref,v.alt,str(v.qual),v.filter,v.info,"\n"]))    
        else:
            for key in variantSet:
                v = variantDict[key]
                print "\t".join([v.chromosome,str(v.position),v.ref,v.alt,str(v.qual),v.filter,v.info,"\n"])
            
    
    
    def getVariantTuble(self,line):
        '''
        returns a tuple of (chromosome, position, alt, ref) from a vcfFile
        '''
        line=line.split("\t")
        try:
            tuple = (line[0],int(line[1]),line[3],line[4])
            #tuple = (line[0],int(line[1]))
        except IndexError:
            raise ValueError("Error in line '%s'" % " ".join(line))
        return tuple
    
    def deleteOverlappsFromA(self,variantsA,variantsB):
        '''
        delete the variants from 'variantsA' which also are in 'variantsB'
        '''
    
        
        #determine type of variantsA 
        if type(variantsA) == dict:
            variantSetA = set(variantsA.keys())
        elif type(variantsA) == file or type(variantsA) == str:
            variantsA = self.parseVcfFile_variantsDict(variantsA)
            variantSetA = set(variantsA.keys())
        else: 
            raise TypeError("variantsA has wrong type, need variantDict, str or file, %s found" % type(variantsA))
        
        #detrmine type of variantB
        if type(variantsB) == str:
            variantsB = open(variantsB)
        elif type(variantsB) != file:
            raise TypeError("variantB has wrong type, need str or file, %s found" % type(variantsB))
        
        #get Start time
        startTime = Helper.getTime()
        Helper.info(" [%s] Delete overlapps from %s" % (startTime.strftime("%c"),variantsB.name))
    
        
        
        for line in variantsB:
            if line.startswith("#"):
                continue
            Btupple = self.getVariantTuble(line)
            if Btupple in variantSetA:
                #A.discard(Btupple)
                variantSetA.remove(Btupple)
                del variantsA[Btupple]
        
        #calculate duration 
        Helper.printTimeDiff(startTime)
        return  variantsA    
    
    def sortVariantDict(self,variantDict):
        '''
        Sorts a VariantDictionary by the variant position
        :param variantDict:
        '''
        #if type(variantDict) != list:
        #    raise TypeError("variants has wrong type, need variantDict, %s found" % type(variantDict))
        for key in variantDict.keys():
            variantDict[key] = sorted(variantDict[key], key=operator.attrgetter('position'))
            