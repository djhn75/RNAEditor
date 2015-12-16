'''
Created on 05.06.2014

@author: david
'''

#!/usr/bin/env python

import pysam
from Helper import Helper
import os
from VariantSet import VariantSet


class Dink():
    
    def __init__(self):
        class RnaEdit():
            def __init__(self):
                self.logFile=None
                self.textField=0
        self.rnaEdit = RnaEdit()        
        
        
    '''remove missmatches from homopolymers'''
    def removeHomopolymers(self,variants,outFile,distance):
        startTime=Helper.getTime()
        Helper.info(" [%s] remove Missmatches from homopolymers " % (startTime.strftime("%c")),self.rnaEdit.logFile,self.rnaEdit.textField)
        
        tempBedFile = open(outFile+"_tmp.bed","w+")
        tempSeqFile = outFile + "_tmp.tsv"
        
        refGenome = "/media/Storage/databases/rnaEditor_annotations/human/human_g1k_v37.fasta"
        fastaFile = pysam.FastaFile(self.rnaEdit.refGenome)
        mmNumberTotal = len(variants.variantDict)
        #print temporary BedFile
        numberPassed=0
        for key in variants.variantDict.keys():
            chr,position,ref,alt = key
            startPos = position - distance if position >= distance else 0
            endpos = position + distance
            sequence = fastaFile.fetch(chr,startPos,endpos)
            pattern = ref*distance
            
            """ !!!Test if this gives better results
                !!!ONLY DELETE IF MM IS AT THE END OF A HOMOPOLYMER NUKLEOTIDES
            if sequence.startswith(pattern):
                del mmDict[site] 
            elif sequence.endswith(pattern):
                del mmDict[site]
            """
            if pattern in sequence:
                try:
                    del variants.variantDict[key]
                except KeyError:
                    pass
            else:
                numberPassed+=1
                    
        #output statistics
        Helper.info("\t\t %d out of %d passed the Homopolymer-Filter" % (numberPassed, mmNumberTotal),self.rnaEdit.logFile,self.rnaEdit.textField)
        Helper.printTimeDiff(startTime,self.rnaEdit.logFile,self.rnaEdit.textField)

refGenome = "/media/Storage/databases/rnaEditor_annotations/human/human_g1k_v37.fasta"
vcfFile = "/media/Storage/bio-data/test/rnaEditor/scrambleN/scrambleN.editingSites.vcf"
output= "/media/Storage/bio-data/test/rnaEditor/scrambleN/scrambleN"

variants = VariantSet(vcfFile)
d=Dink()

d.removeHomopolymers(variants, output, 5)