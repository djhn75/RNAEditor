#!/usr/bin/python
'''
Created on May 22, 2013

@author: david
'''

import argparse, os, multiprocessing
from Helper import Helper
from string import rfind


class MapFastq(object):
    '''
    Maps a fastQ file to the given genome
    '''
    
    def __init__(self,fastqFile,refGenome,dbsnp,outfilePrefix="default",sourceDir="bin/",
                 threads=multiprocessing.cpu_count()-1,maxDiff=0.04,seedDiff=2,
                 paired=False,keepTemp=False,overwrite=False):
        '''
        Constructor
        '''
        self.debug=True
        
        self.fastqFile=fastqFile
        self.refGenome=refGenome
        self.dbsnp=dbsnp
        if outfilePrefix=="default":
            self.outfilePrefix=self.fastqFile[0:self.fastqFile.rfind(".")]
        else:
            self.outfilePrefix=outfilePrefix
        self.sourceDir=sourceDir
        self.threads=str(threads)
        self.maxDiff=str(maxDiff)
        self.seedDiff=str(seedDiff)
        self.paired=paired
        self.keepTemp=keepTemp
        self.overwrite=overwrite
        
        
        self.logFile=open(self.outfilePrefix + ".log","w+")
        if self.debug==True:
            self.printAttributes()
        
        self.checkDependencies()
    
    
    def printAttributes(self):
        print
        print "*** MAP READS WITH FOLLOWING ATTRIBUTES ***"
        print "\t FastQ-File: " + self.fastqFile
        print "\t outfilePrefix:" + self.outfilePrefix
        print "\t refGenome:" + self.refGenome
        print "\t dbsnp:" + self.dbsnp
        print "\t sourceDir:" + self.sourceDir
        print "\t threads:" + self.threads
        print "\t maxDiff:" + self.maxDiff
        print "\t seedDiff:" + self.seedDiff
        print "\t paired:" + str(self.paired)
        print "\t keepTemp:" + str(self.keepTemp)
        print "\t overwrite:" + str(self.overwrite)
        print

    '''check if all the needed files are there'''
    def checkDependencies(self):
        if not os.path.exists(self.dbsnp):
            Exception("dbSNP File: Not found!!!")
        if not os.path.exists(self.refGenome):
            Exception("reference Genome File: Not found!!!")

    def start(self):
        
        #Align Fastq Reads to the Genome
        saiFile=self.outfilePrefix+".sai"
        cmd = [self.sourceDir+"bwa", "aln" , "-t",self.threads, "-n", self.maxDiff , "-k", self.seedDiff, self.refGenome, self.fastqFile]
        Helper.proceedCommand("Align Reads with BWA", cmd, self.fastqFile, saiFile, self.logFile, self.overwrite)
        
        #convert sai to sam
        samFile=self.outfilePrefix+".sam"
        #TODO:check for paired
        cmd = [self.sourceDir + "bwa", "samse", "-r", "@RG\tID:A\tLB:A\tSM:A\tPL:ILLUMINA\tPU:HiSEQ2000", self.refGenome, saiFile, self.fastqFile]
        Helper.proceedCommand("convert sai to sam", cmd, saiFile, samFile, self.logFile, self.overwrite)
    
        #convert sam to bam
        bamFile=self.outfilePrefix+".bam"
        cmd=["java", "-Xmx4G", "-jar", self.sourceDir + "picard-tools/SortSam.jar", "INPUT=" + samFile, "OUTPUT=" + bamFile, "SO=coordinate", "VALIDATION_STRINGENCY=LENIENT", "CREATE_INDEX=true"]
        Helper.proceedCommand("convert sam to bam", cmd, saiFile, bamFile, self.logFile, self.overwrite)
        
        
        #run Alignement with tophat
        """
        bamFile=self.outfilePrefix+"/accepted_hits.bam"
        cmd=[self.sourceDir + "tophat/tophat2", "--no-coverage-search","--keep-fasta-order", "-p", "12", "--rg-id", "A","--rg-sample","A","--rg-library","illumina","--rg-platform-unit","HiSeq", "-o", self.outfilePrefix, self.refGenome, self.fastqFile ]
        print cmd
        Helper.proceedCommand("Map reads with tophat", cmd, self.fastqFile, bamFile, self.logFile, self.overwrite)
        """
        
        #sort bam
        #sortBamFile=self.outfilePrefix+".bam"
        #cmd=["java", "-Xmx4G", "-jar", self.sourceDir + "picard-tools/SortSam.jar", "INPUT=" + bamFile, "OUTPUT=" + sortBamFile, "SO=coordinate", "VALIDATION_STRINGENCY=LENIENT", "CREATE_INDEX=true"]
        #Helper.proceedCommand("sortq bam", cmd, bamFile, sortBamFile, self.logFile, self.overwrite)
        
        #Add read group ONLY NEEDED WHEN MAPPED WITH TOPHAT
        #rgFile=self.outfilePrefix+".bam"
        #cmd=["java", "-Xmx4G", "-jar", self.sourceDir + "picard-tools/AddOrReplaceReadGroups.jar", "INPUT=" + bamFile, "OUTPUT=" + rgFile, "SO=coordinate", "VALIDATION_STRINGENCY=LENIENT", "CREATE_INDEX=true", "ID=A", "LB=A", "SM=A", "PL=illumina", "PU=HiSeq2000", "SM=A"]
        #Helper.proceedCommand("Add read Groups", cmd, bamFile, rgFile, self.logFile, self.overwrite)
        
        
        #Identify Target Regions for realignment
        intervalFile=self.outfilePrefix+".indels.intervals"
        cmd=["java","-Xmx8G","-jar",self.sourceDir + "GATK/GenomeAnalysisTK.jar", "-nt",self.threads, "-T", "RealignerTargetCreator", "-R", self.refGenome, "-I", bamFile, "-o", intervalFile,"-l", "ERROR"]
        Helper.proceedCommand("Identify Target Regions for realignment", cmd, bamFile, intervalFile, self.logFile, self.overwrite)
        
        #Proceed Realignement
        realignedFile=self.outfilePrefix+".realigned.bam"
        cmd=["java","-Xmx8G","-jar",self.sourceDir + "GATK/GenomeAnalysisTK.jar", "-T", "IndelRealigner", "-R", self.refGenome, "-I", bamFile, "-l", "ERROR", "-targetIntervals", intervalFile, "-o", realignedFile]
        Helper.proceedCommand("Proceed Realignement", cmd, intervalFile, realignedFile, self.logFile, self.overwrite)
        
        #mark PCR duplicates
        markedFile=self.outfilePrefix+".realigned.marked.bam"
        cmd=["java","-Xmx8G","-jar",self.sourceDir + "picard-tools/MarkDuplicates.jar","INPUT=" + realignedFile, "OUTPUT=" + markedFile, "METRICS_FILE="+self.outfilePrefix+".pcr.metrics", "VALIDATION_STRINGENCY=LENIENT", "CREATE_INDEX=true"]
        Helper.proceedCommand("mark PCR duplicates", cmd, realignedFile, markedFile, self.logFile, self.overwrite)
        
        #Find Quality Score recalibration spots
        recalFile=self.outfilePrefix+".recalSpots.grp"
        cmd=["java","-Xmx8G","-jar",self.sourceDir + "GATK/GenomeAnalysisTK.jar", "-T", "BaseRecalibrator", "-l", "ERROR", "-R", self.refGenome, "-knownSites", self.dbsnp, "-I", markedFile, "-cov", "CycleCovariate", "-cov", "ContextCovariate", "-o", recalFile]
        Helper.proceedCommand("Find Quality Score recalibration spots", cmd, markedFile, recalFile, self.logFile, self.overwrite)
        
        #proceed Quality Score recalibration
        recaledBamFile=self.outfilePrefix+".realigned.marked.recalibrated.bam"
        cmd=["java","-Xmx8G","-jar",self.sourceDir + "GATK/GenomeAnalysisTK.jar", "-T", "PrintReads","-l", "ERROR", "-R", self.refGenome, "-I", markedFile, "-BQSR", recalFile, "-o",recaledBamFile]
        Helper.proceedCommand("Proceed Quality Score recalibration", cmd, recalFile, recaledBamFile, self.logFile, self.overwrite)
        
        return recaledBamFile
        
    def __del__(self):
        if self.keepTemp==False:
            os.remove(self.outfilePrefix+".sai")
            os.remove(self.outfilePrefix+".sam")
            #os.remove(self.outfilePrefix+".bam")
            os.remove(self.outfilePrefix+".indels.intervals")
            os.remove(self.outfilePrefix+".realigned.bam")
            os.remove(self.outfilePrefix+".realigned.bai")
            os.remove(self.outfilePrefix+".realigned.marked.bam")
            os.remove(self.outfilePrefix+".realigned.marked.bai")
            os.remove(self.outfilePrefix+".recalSpots.grp")
            #os.remove(self.outfilePrefix+".realigned.marked.recalibrated.bam")
    

'''
when the script is called directly
''' 
if __name__ == '__main__':
    #parse command line arguments and set defaults
    parser = argparse.ArgumentParser(description='map FastQ Files to the given genome and realigns the reads for SNP-calling.')
    parser.add_argument('-i', '--input', metavar='Fastq-File', type=argparse.FileType('r'), help='Input fastq file that should be mapped to the genome', required=True)
    parser.add_argument("-r", "--RefGenome", metavar='Fasta-File', help="File that contains the reference sequences", type=argparse.FileType('r'), default='/media/databases/human/human_g1k_v37.fa')
    parser.add_argument('-s', '--dbsnp', help=' SNP database (dbSNP) in VCF format (downloaded from the GATK homepage)', type=argparse.FileType('r'), default='/media/databases/human/dbsnp_135.b37.vcf')
    parser.add_argument('-o', '--output', metavar='output-prefix', type=str,help='prefix that is written in front of the output files', default="default")
    parser.add_argument('-d', '--sourceDir', help='- Directory to all the tools [default: /bin/]', default='bin/', type=Helper.readable_dir)
    parser.add_argument('-t', '--threads', help='number of threads', type=int, default=multiprocessing.cpu_count()-1)
    parser.add_argument('-n', '--maxDiff', help=' maximum Number of mismatches in the reads (int) or error rate in percentage (float)[0.04]', type=float, default=0.04)
    parser.add_argument('--seedDiff', help='maximum Number of mismatches in the seed sequence (int)[2]', type=int, default=2)
    parser.add_argument('-p', '--paired', help="Use this paramater if you have paired end reads [false]", action='store_true', default=False)
    parser.add_argument('--keepTemp', help='keep the intermediate Files [False]', action='store_true', default=False)
    parser.add_argument('--overwrite', help='overwrite existing Files [False]', action='store_true', default=False)
    
    args = parser.parse_args()
    
    mapFastQ=MapFastq(args.input.name, args.RefGenome.name, args.dbsnp.name,args.output, args.sourceDir, args.threads, args.maxDiff, args.seedDiff, args.paired, args.keepTemp, args.overwrite)
    mapFastQ.start()
    del mapFastQ