#!/usr/bin/python
'''
Created on Nov 13, 2013
Script to add Annotation to a VCF file 
@author: David John
'''
import argparse, sys
from string import split
from Helper import Helper

parser = argparse.ArgumentParser(description="Remap the missmatch reads to confirm unique mapping")
parser.add_argument("-v","--variants",dest="variantFile", help='Variant file in VCF-format', required=True)
parser.add_argument("-b","--bam",dest="bamFile", help="bam-file with reads", required=True)
parser.add_argument("-r","--refGenome",dest="refGenome", help="Reference Genome", required=True)
parser.add_argument("-o","--outFile",help="File to print remaining missmatches", required=True)

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
    

args= parser.parse_args()

#set Files
variantFile = open(args.variantFile)
bamFile=args.bamFile
tempFasta = args.variantFile + "_tmp.fa"
print tempFasta

geneHash = {}
for line in variantFile:
    line=line.split("\t")
    chromosome,snpPos=line[0],line[1]
    position=line[0]+":"+snpPos+"-"+snpPos
    
    samout = Helper.getCommandOutput("samtools view -F 1024 " + bamFile + " " + position).splitlines() #-F 1024 to filter out duplicate reads
    for samLine in samout:
        samfields=samLine.split()
        flag,startPos,mapQual,cigar,sequence,seqQual = samfields[1],samfields[3],samfields[4],samfields[5],samfields[9],samfields[10]
                       
        #print fasta file
        print "> " + position + "\n" + sequence
        #print >>