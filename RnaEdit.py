#!/usr/bin/python

'''
    Created on May 17, 2013
    Main Class to detect RNA-editing in a given FastQ file
    @author: david
'''

from Helper import Helper, Parameters
from MapFastq import MapFastq
from CallEditingSites import CallEditingSites
import multiprocessing, argparse, os


class RnaEdit(object):

    def __init__(self, fastqFiles, refGenome, dbsnp, 
                 hapmap, omni, esp, 
                 aluRegions, geneAnnotation, outfilePrefix="default", 
                 sourceDir="/usr/local/bin", threads=multiprocessing.cpu_count()-1,maxDiff=0.04, 
                 seedDiff=2, paired=False, standCall=0, 
                 standEmit=0, edgeDistance=3, keepTemp=False, 
                 overwrite=False,runNumber=0):
        
        
        #check if the input Files are there
        self.checkDependencies(fastqFiles, refGenome, dbsnp, 
                 hapmap, omni, esp, 
                 aluRegions, geneAnnotation, outfilePrefix, 
                 sourceDir, threads,maxDiff, 
                 seedDiff, paired, standCall, 
                 standEmit, edgeDistance, keepTemp, 
                 overwrite)
        
        """
        START MAPPING
        """
        self.mapFastQ=MapFastq(fastqFiles, refGenome, dbsnp, outfilePrefix, sourceDir, threads, maxDiff, seedDiff, paired, keepTemp, overwrite,runNumber)
        mapResultFile=self.mapFastQ.start()
        
        #print mapResultFile + " was created \t Mapping Process finished"
        del self.mapFastQ
        
        """
        START CALLING EDITING SITES
        """
        self.callEditSites=CallEditingSites(mapResultFile, refGenome, dbsnp, 
                                            hapmap, omni, esp, 
                                            aluRegions, geneAnnotation, outfilePrefix, 
                                            sourceDir, threads, standCall,
                                            standEmit, edgeDistance, keepTemp, 
                                            overwrite,runNumber)
        
        #self.callEditSites.start()
        del self.callEditSites
        
        Helper.status("rnaEditor Finished with %s" % outfilePrefix)
     
    def __del__(self):
        #del self.mapFastQ
        pass
        

    def checkDependencies(self,fastqFiles, refGenome, dbsnp, 
                 hapmap, omni, esp, 
                 aluRegions, geneAnnotation, outfilePrefix, 
                 sourceDir, threads,maxDiff, 
                 seedDiff, paired, standCall, 
                 standEmit, edgeDistance, keepTemp, 
                 overwrite):
        '''
        Checks the existence of the necessary packages and tools
        :param sourceDir: folder which contains all the software
        '''
        Helper.newline(1)
        Helper.info("CHECK DEPENDENCIES")
        
        #check if all tools are there
        if not os.path.isfile(sourceDir+"bwa"):
            Helper.error("BWA not found in %s" % sourceDir)
        if not os.path.isfile(sourceDir+"picard-tools/SortSam.jar"):
            Helper.error("SortSam.jar not found in %s" % sourceDir+"picard-tools")
        if not os.path.isfile(sourceDir+"picard-tools/MarkDuplicates.jar"):
            Helper.error("MarkDuplicates.jar not found in %s" % sourceDir+"picard-tools")
        if not os.path.isfile(sourceDir+"GATK/GenomeAnalysisTK.jar"):
            Helper.error("GenomeAnalysisTK.jar not found in %s" % sourceDir+"GATK/")
        if not os.path.isfile(sourceDir+"bedtools/fastaFromBed"):
            Helper.error("fastaFromBed not found in %s" % sourceDir+"bedtools/")
        if not os.path.isfile(sourceDir+"blat"):
            Helper.error("blat not found in %s" % sourceDir)
        if not os.path.isfile(sourceDir+"samtools"):
            Helper.error("samtools not found in %s" % sourceDir)
        if not os.system("java -version")==0:
            Helper.error("Java could not be found, Please install java")
        
        #check if all files are there
        if not os.path.isfile(refGenome):
            Helper.error("Could not find Reference Genome in %s: " % refGenome)
        # Files for BWA
        if not os.path.isfile(refGenome+".amb"):
            Helper.error("Could not find %s.amb" % refGenome)
            Helper.error("run: 'bwa index %s' to create it" % refGenome)
        if not os.path.isfile(refGenome+".ann"):
            Helper.error("Could not find %s.ann" % refGenome)
            Helper.error("run: 'bwa index %s' to create it" % refGenome)
        if not os.path.isfile(refGenome+".bwt"):
            Helper.error("Could not find %s.bwt" % refGenome)
            Helper.error("run: 'bwa index %s' to create it" % refGenome)
        if not os.path.isfile(refGenome+".pac"):
            Helper.error("Could not find %s.pac" % refGenome)
            Helper.error("run: 'bwa index %s' to create it" % refGenome)
        if not os.path.isfile(refGenome+".sa"):
            Helper.error("Could not find %s.sa" % refGenome)
            Helper.error("run: 'bwa index %s' to create it" % refGenome)
        
        #Files for GATK
        
        if not os.path.isfile(refGenome.replace(".fastq",".dict")):
            Helper.error("Could not find %s" % refGenome.replace(".fastq",".dict"))
            Helper.error("run: 'java -jar %s/picard-tools/CreateSequenceDictionary.jar R=%s  O= %s.dict' to create it" % (sourceDir,refGenome,refGenome))
        if not os.path.isfile(refGenome+".fai"):
            Helper.error("Could not find %s.sai" % refGenome)
            Helper.error("run: 'samtools faidx %s' to create it" % refGenome)
    
        #SNP databases
        if not os.path.isfile(dbsnp):
            Helper.error("Could not find %s: " % dbsnp)
        if not os.path.isfile(hapmap):
            Helper.error("Could not find %s: " % hapmap)
        if not os.path.isfile(omni):
            Helper.error("Could not find %s: " % omni)
        if not os.path.isfile(esp):
            Helper.error("Could not find %s: " % esp)
        
        #region Files
        if not os.path.isfile(aluRegions):
            Helper.error("Could not find %s: " % aluRegions)
        if not os.path.isfile(geneAnnotation):
            Helper.error("Could not find %s: " % geneAnnotation)

if __name__ == '__main__':
    Parameters.readDefaults()
    
    parser = argparse.ArgumentParser(description='map FastQ Files to the given genome and realigns the reads for SNP-calling.',)
    parser.add_argument('-i', '--input', metavar='Fastq-Files',nargs='+', type=str, help='Input fastq files (maximum two for paire-end-sequencing)', required=True)
    parser.add_argument('-r', '--RefGenome', metavar='Fasta-File', help='File that contains the reference sequences', type=str, default=Parameters.refGenome)
    parser.add_argument('-s', '--dbsnp', help=' SNP database (dbSNP) in VCF format (downloaded from the GATK homepage)', type=str, default=Parameters.dbSNP)
    parser.add_argument('-m', '--hapmap', help='hapmap database in vcf format (see GATK homepage)', type=str, default=Parameters.hapmap)
    parser.add_argument('-g', '--omni', help='1000 Genome variants in vcf format (see GATK homepage)', type=str, default=Parameters.omni)
    parser.add_argument('-e', '--esp', help='Exome Sequencing Project variants', type=str, default=Parameters.esp)
    parser.add_argument('-a', '--aluRegions', help='Alu-Regions downloaded fron the UCSC table browser', type=str, default=Parameters.aluRegions)
    parser.add_argument('-G', '--geneAnnotation', help='Gene annotation File in bed format', type=str, default=Parameters.gtfFile)
    parser.add_argument('-o', '--output', metavar='output-prefix', type=str,help='prefix that is written in front of the output files', default=Parameters.output)
    parser.add_argument('-d', '--sourceDir', help='- Directory to all the tools [default: bin/]', default=Parameters.binary, type=str)
    parser.add_argument('-t', '--threads', help='number of threads', type=int, default=multiprocessing.cpu_count()-1)
    parser.add_argument('-n', '--maxDiff', help=' maximum Number of mismatches in the reads (int) or error rate in percentage (float)[0.04]', type=float, default=Parameters.maxDiff)
    parser.add_argument('--seedDiff', help='maximum Number of mismatches in the seed sequence (int)[2]', type=int, default=Parameters.seedDiff)
    parser.add_argument('-p', '--paired', help="Use this paramater if you have paired end reads [false]", action='store_true', default=Parameters.paired)
    parser.add_argument('-sc', '--standCall', help='-The minimum phred-scaled confidence threshold at which variants should be considered as true (int) [0]', type=int, default=Parameters.standCall)
    parser.add_argument('-se', '--standEmit', help=' The minimum phred-scaled confidence threshold at which variants should be emitted (int)[0]', type=int, default=Parameters.standEmit)
    parser.add_argument('-ed', '--edgeDistance', help='The minimum edge distance of the SNPs', type=int, default=Parameters.edgeDistance)
    parser.add_argument('--keepTemp', help='keep the intermediate Files [False]', action='store_true', default=Parameters.keepTemp)
    parser.add_argument('--overwrite', help='overwrite existing Files [False]', action='store_true', default=Parameters.overwrite)
    parser.add_argument('--index', type=int, help='Tab on which to output (int)', default=0)
    
    args = parser.parse_args()
    
    edit=RnaEdit(args.input, args.RefGenome, args.dbsnp,
                 args.hapmap, args.omni, args.esp, 
                 args.aluRegions, args.geneAnnotation, args.output, 
                 args.sourceDir, args.threads, args.maxDiff, 
                 args.seedDiff, args.paired, args.standCall,
                 args.standEmit,args.edgeDistance, args.keepTemp, 
                 args.overwrite,args.index)
    
    del edit
    
    
    
else:
    pass    
    
    