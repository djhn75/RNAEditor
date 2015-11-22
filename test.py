'''
Created on 05.06.2014

@author: david
'''

#!/usr/bin/env python

from VariantSet import VariantSet
from Helper import Helper
import re
import pysam




def removeEdgeMissmatches(variants,bamFile,minDistance, minBaseQual):
        '''delete variants from Bam file which appear near read edges'''
        #Loop through vcf-File
            #call overlapping reads with samtools view
            #loop over reads
                #discard variants wich appear ONLY near edges
                #write the rest to the output file
        startTime=Helper.getTime()
        minDistance=int(minDistance)
        
        counter=0    
        
        num_lines = len(variants.variantDict)
        Helper.info(" [%s] remove Missmatches from the first %s bp from read edges" % (startTime.strftime("%c"),str(minDistance)))
        
        for varKey in variants.variantDict.keys():
            variant = variants.variantDict[varKey]
            snpPos = variant.position
            position=variant.chromosome+":" + str(snpPos) + "-" + str(snpPos) 
            #line[1]+"-"+line[1]
            keepSNP=False
            
            command = ["/usr/local/bin/"+"samtools", "view", bamFile, position] 
            samout = Helper.getCommandOutput(command).splitlines() #get the reads wich are overlapping the snp region
            for samLine in samout: #loop over reads
                samfields=samLine.split()
                try:
                    flag,startPos,mapQual,cigar,sequence,seqQual = samfields[1],int(samfields[3]),samfields[4],samfields[5],samfields[9],samfields[10]
                except ValueError:
                    raise ValueError("Error in line '%s'" % " ".join(samfields))

                readPos=0
                mmReadPos=0
                cigarNums=re.split("[MIDNSHP]", cigar)[:-1]
                cigarLetters=re.split("[0-9]+",cigar)[1:]
                
                #loop over read cigar (check insertions,deletions and skipped regions) 
                for i in range(len(cigarLetters)): #parse over single read
                    if cigarLetters[i] in {"I","S","H"}: #Insertion, Soft Clipping and Hard Clipping
                        readPos = readPos + int(cigarNums[i])
                    elif cigarLetters[i] in {"D","N"}: #Deletions and skipped Regions
                        startPos = startPos + int(cigarNums[i])
                    elif cigarLetters[i] in {"M"}: #Matches
                        for j in range(int(cigarNums[i])):
                            if startPos == snpPos:
                                mmReadPos = readPos
                            readPos += 1
                            startPos += 1         
                if mmReadPos != 0: #happens when snp is in non matching regions (like Insertion, Soft Clipping and Hard Clipping, Deletions ans skipped regions)    
                    edgeDistance = snpPos - startPos
                
                    #only remove the snps from first minDistance bases
                    revStrand = int(flag) & 16 #check if 5th bit is set (results: 0 = +strand; 16= -strand)
                    if (revStrand == 0 and mmReadPos > minDistance) or (revStrand == 16 and mmReadPos < readPos - minDistance):
                        mmBaseQual= ord(seqQual[mmReadPos])
                        mmReadBase= sequence[mmReadPos]
                        if(mmBaseQual >= minBaseQual + 33) and (mmReadBase == variant.alt): #check for quality of the base and the read contains the missmatch
                            keepSNP=True

            if not keepSNP:
                del variants.variantDict[varKey]    
            counter+=1
            if counter % 10000 == 0: #print out current status
                Helper.status(str(counter) + " of " + str(num_lines) + " missmatches finished",)




            
bamFile="/media/Storage/bio-data/David/test/Icm4.realigned.marked.recalibrated.bam"

var = VariantSet("/media/Storage/bio-data/David/test/Icm4.1.vcf")
removeEdgeMissmatches2(var, bamFile, 26, 20)

#removeEdgeMissmatches(var, bamFile, 26, 20)

"""
bamFile="/media/ATLAS_NGS_storage/Till/fastq/demultiplexed/rnaEditor/Icm4.realigned.marked.recalibrated.bam"
bamFile = pysam.AlignmentFile(bamFile, "rb")
iter = bamFile.fetch("1", 14906, 14907)
for x in iter:
    print ("coverage at base %s = %s" % (x.pos))

bamFile.close()

"""
"""
var.deleteNonEditingBases()

bedFile="/media/Storage/databases/rnaEditor_annotations/human/Alu_repeats_noCHR.bed"
#bedFile="/media/Storage/databases/rnaEditor_annotations/human/Alu_repeats_noCHR_Y.bed"


alu,nonAlu=var.splitByBed(bedFile)

savedAlu=VariantSet("/media/ATLAS_NGS_storage/Till/fastq/demultiplexed/rnaEditor/Icm4.alu.vcf")
savedNonAlu=VariantSet("/media/ATLAS_NGS_storage/Till/fastq/demultiplexed/rnaEditor/Icm4.nonAlu.vcf")

savedAlu.deleteNonEditingBases()
savedNonAlu.deleteNonEditingBases()

print "original vars", len(var)
#print "alu:",len(alu), "nonAlu:", len(nonAlu)
print "savedAlu:",len(savedAlu), "savedNonAlu:", len(savedNonAlu)

"""


