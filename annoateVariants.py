#!/usr/bin/python
'''
This script takes several vcf files with editing sites and dertemines editing islands
Created on 24.09.17

@author: david
'''

import argparse, os, sys
from string import split
from VariantSet import VariantSet
from Genome import Genome



parser = argparse.ArgumentParser(description='reanalyze editing islands.')
parser.add_argument("-f", "--files", metavar="N",type=str, nargs="+", help="list of variant files (vcf)", required=True)
parser.add_argument("-g", "--genome", metavar="N",type=str, help="Genome file (GTF", required=True)
parser.add_argument("-o", "--out", metavar = "N", type = str, help = "outputDir", default = "./")

args = parser.parse_args()

genome = Genome(args.genome)
for file in args.files:  # loop through all files
    samplename=file[file.rfind('/')+1:file.rfind('.vcf')]
    variants = VariantSet(file)

    variants.annotateVariantDict(genome)
    variants.printVariantDict(args.out + samplename +'.annotated.vcf')