'''
Created on 05.06.2014

@author: david
'''
from Helper import Helper
from collections import defaultdict
import os
import operator


class Variant:
    """
    handle one feature (line) of a GTF file
    """
    
    def __init__(self):
        self.chromosome = ""
        self.position = 0
        self.id = ""
        self.ref = ""
        self.alt = ""
        self.qual = 0
        self.filter = ""
        self.info = ""
        self.baseCounts = [0,0,0,0] #[A,C,G,T] number of reads which contain the appropriate base 
        
    
    def readline(self,line):
        """
        process one line of the vcf file
        <chromosome> <position> <identifier> <reference_base> <alternative_base> <quality> <filter> {attributes}
        """
        line = line.rstrip().split("\t")
        
        try:
            self.chromosome = line[0]
            self.position = int(line[1])
            self.id = line[2]
            self.ref = line[3]
            self.alt = line[4]
            self.qual = float(line[5]) if line[5] !="." else 0.0
            self.filter = line[6]
            self.info = line[7]
        except ValueError:
            raise ValueError("Error in line '%s'" % " ".join(line))

        
        #parse info
        info = line[7]
        #trim comments
        info=info[:info.find("#")].rstrip()
        
        baseCountPos = info.find("BaseCounts")
        if baseCountPos != -1:
            #extract BaseCounts from list and save it in a list
            self.baseCounts = info[baseCountPos:info.find(";", baseCountPos)].split("=")[1].split(",") 

def checkVariantType(variants):
        """
        Checks if the type of the argument is a str or a file
        returns a dictionary of all the variants
        """
        if type(variants) == dict:
            return variants
        elif type(variants) == file or type(variants) == str:
            variants = parseVcfFile_variantsDict(variants)
            return variants
            
        else: 
            raise TypeError("variants has wrong type, need variantDict, str or file, %s found" % type(variants))
     
         
def iterator(infile):

    while 1:
        line = infile.readline()
        if not line: raise StopIteration
        if line.startswith("#"): continue #skip comments
        gtf = Variant()
        gtf.readline(line)
        yield gtf

def parseVcfFile_variantSetByChromosome(vcfFile):
    """
    Imports a given Variant File and returns the variants as Dictionary with chromosome as key and a list of VariantObjects as values
    {"1":[VariantObject1,VariantObject2....],"2":[VariantObject1,VariantObject2....]}
    
    """
    startTime = Helper.getTime()
    Helper.info(" [%s] Parsing Variant Data from %s" % (startTime.strftime("%c"),vcfFile))
    
    #check correct Type
    if type(vcfFile) == str:
        vcfFile = open(vcfFile)
    elif type(vcfFile) != file:
        raise TypeError("Invalid type in 'parseVcfFile' (need string or file, %s found)" % type(vcfFile)) 
        
    variantsByChromosome = defaultdict(list)
    for v in iterator(vcfFile):
        variantsByChromosome[v.chromosome].append(v)
    
    Helper.printTimeDiff(startTime)
    return variantsByChromosome

def parseVcfFile_variantsDict(vcfFile):
    """
    Imports a given Variant File and returns the variants as Dictionary with Tuple of (chromosome,pos,ref,alt) as key and a the VariantObject as value
    {(1,45435,"A","G"):VariantObject1,(1,45435,"A","G"):VariantObject1,.....}
    
    """
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
    for v in iterator(vcfFile):
        variantDict[(v.chromosome,v.position,v.ref,v.alt)]=v
        #variantDict[(v.chromosome,v.position)]=v
    
    Helper.printTimeDiff(startTime)
    return variantDict


def printVariantDict(variantDict,outfile=None):
    """
    print the variants from the dictionary to the outfile if defined
    """
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
        
def printVariantSet(variantSet,variantDict,outfile=None):
    """
    prints variants from the set to the outfile if defined
    """ 
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
        


def getVariantTuble(line):
    """
    returns a tuple of (chromosome, position, alt, ref) from a vcfFile
    """
    line=line.split("\t")
    try:
        tuple = (line[0],int(line[1]),line[3],line[4])
        #tuple = (line[0],int(line[1]))
    except IndexError:
        raise ValueError("Error in line '%s'" % " ".join(line))
    return tuple

def deleteOverlappsFromA(variantsA,variantsB):
    """
    delete the variants from 'variantsA' which also are in 'variantsB'
    """

    
    #determine type of variantsA 
    if type(variantsA) == dict:
        variantSetA = set(variantsA.keys())
    elif type(variantsA) == file or type(variantsA) == str:
        variantsA = parseVcfFile_variantsDict(variantsA)
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
        Btupple = getVariantTuble(line)
        if Btupple in variantSetA:
            #A.discard(Btupple)
            variantSetA.remove(Btupple)
            del variantsA[Btupple]
    
    #calculate duration 
    Helper.printTimeDiff(startTime)
    return  variantsA    

def sortVariantDict(variantDict):
    '''
    Sorts a VariantDictionary by the variant position
    :param variantDict:
    '''
    #if type(variantDict) != list:
    #    raise TypeError("variants has wrong type, need variantDict, %s found" % type(variantDict))
    for key in variantDict.keys():
        variantDict[key] = sorted(variantDict[key], key=operator.attrgetter('position'))
        