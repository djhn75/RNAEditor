#!/usr/bin/python

'''
    Created on May 17, 2013
    Main Class to detect RNA-editing in a given FastQ file
    @author: david
'''

from Helper import Helper
from MapFastq import MapFastq
from CallEditingSites import CallEditingSites
import multiprocessing, argparse


class RnaEdit(object):
    
    
    

    def __init__(self, fastqFiles, refGenome, dbsnp, 
                 hapmap, omni, esp, 
                 aluRegions, geneAnnotation, outfilePrefix="default", 
                 sourceDir="/usr/local/bin", threads=multiprocessing.cpu_count()-1,maxDiff=0.04, 
                 seedDiff=2, paired=False, standCall=0, 
                 standEmit=0, edgeDistance=3, keepTemp=False, 
                 overwrite=False):
              
        self.mapFastQ=MapFastq(fastqFiles, refGenome, dbsnp, outfilePrefix, sourceDir, threads, maxDiff, seedDiff, paired, keepTemp, overwrite)
        mapResultFile=self.mapFastQ.start()
        
        #print mapResultFile + " was created \t Mapping Process finished"
        del self.mapFastQ
        
        self.callEditSites=CallEditingSites(mapResultFile, refGenome, dbsnp, 
                                            hapmap, omni, esp, 
                                            aluRegions, geneAnnotation, outfilePrefix, 
                                            sourceDir, threads, standCall,
                                            standEmit, edgeDistance, keepTemp, 
                                            overwrite)
        self.callEditSites.start()
     
    def __del__(self):
        #del self.mapFastQ
        del self.callEditSites
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='map FastQ Files to the given genome and realigns the reads for SNP-calling.',)
    parser.add_argument('-i', '--input', metavar='Fastq-Files',nargs='+', type=str, help='Input fastq files (maximum two for paire-end-sequencing)', required=True)
    parser.add_argument('-r', '--RefGenome', metavar='Fasta-File', help='File that contains the reference sequences', type=str, default='/media/media/databases/human/human_g1k_v37.fa')
    parser.add_argument('-s', '--dbsnp', help=' SNP database (dbSNP) in VCF format (downloaded from the GATK homepage)', type=str, default='/media/databases/human/dbsnp_135.b37.vcf')
    parser.add_argument('-m', '--hapmap', help='hapmap database in vcf format (see GATK homepage)', type=str, default='/media/databases/human/hapmap_3.3.b37.sites.vcf')
    parser.add_argument('-g', '--omni', help='1000 Genome variants in vcf format (see GATK homepage)', type=str, default='/media/databases/human/1000G_omni2.5.b37.sites.vcf')
    parser.add_argument('-e', '--esp', help='Exome Sequencing Project variants', type=str, default='/media/databases/human/NHLBI_Exome_Sequencing_Project_6500SI.vcf')
    parser.add_argument('-a', '--aluRegions', help='Alu-Regions downloaded fron the UCSC table browser', type=str, default='/media/databases/human/hg19/rna-editing/Alu_repeats_noChr.bed')
    parser.add_argument('-G', '--geneAnnotation', help='Gene annotation File in bed format', type=argparse.FileType('r'), default='/media/databases/human/hg19/UCSC_Genes_noChr.ucsc')
    parser.add_argument('-o', '--output', metavar='output-prefix', type=str,help='prefix that is written in front of the output files', default="default")
    parser.add_argument('-d', '--sourceDir', help='- Directory to all the tools [default: bin/]', default='bin/', type=Helper.readable_dir)
    parser.add_argument('-t', '--threads', help='number of threads', type=int, default=multiprocessing.cpu_count()-1)
    parser.add_argument('-n', '--maxDiff', help=' maximum Number of mismatches in the reads (int) or error rate in percentage (float)[0.04]', type=float, default=0.04)
    parser.add_argument('--seedDiff', help='maximum Number of mismatches in the seed sequence (int)[2]', type=int, default=2)
    parser.add_argument('-p', '--paired', help="Use this paramater if you have paired end reads [false]", action='store_true', default=False)
    parser.add_argument('-sc', '--standCall', help='-The minimum phred-scaled confidence threshold at which variants should be considered as true (int) [0]', type=int, default=0)
    parser.add_argument('-se', '--standEmit', help=' The minimum phred-scaled confidence threshold at which variants should be emitted (int)[0]', type=int, default=0)
    parser.add_argument('-ed', '--edgeDistance', help='The minimum edge distance of the SNPs', type=int, default=3)
    parser.add_argument('--keepTemp', help='keep the intermediate Files [False]', action='store_true', default=False)
    parser.add_argument('--overwrite', help='overwrite existing Files [False]', action='store_true', default=False)
    
    args = parser.parse_args()
    
    edit=RnaEdit(args.input, args.RefGenome, args.dbsnp,
                 args.hapmap, args.omni, args.esp, 
                 args.aluRegions, args.geneAnnotation, args.output, 
                 args.sourceDir, args.threads, args.maxDiff, 
                 args.seedDiff, args.paired, args.standCall,
                 args.standEmit,args.edgeDistance, args.keepTemp, 
                 args.overwrite)
    del edit
    
    
    
else:
    pass    
    
    