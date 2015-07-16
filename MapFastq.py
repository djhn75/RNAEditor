#!/usr/bin/python
'''
Created on May 22, 2013

@author: david
'''

import argparse, os, multiprocessing
from Helper import Helper
#from RnaEdit import RnaEdit



class MapFastq(object):
    '''
    Maps a fastQ file to the given genome
    '''
    def __init__(self,rnaEdit):
        '''
        Constructor
        '''
        
        self.rnaEdit=rnaEdit
        
        #check read Quality encoding and convert to phred33 quality if necessary
        for i in range(len(self.rnaEdit.fastqFiles)):
            if Helper.isPhred33Encoding(self.rnaEdit.fastqFiles[i], 100, self.rnaEdit.logFile, self.rnaEdit.textField) == False:
                self.rnaEdit.fastqFiles[i]=Helper.convertPhred64toPhred33(self,self.rnaEdit.fastqFiles[i],self.rnaEdit.params.output+ "_" + str(i+1) + "_phred33.fastq",self.rnaEdit.logFile,self.rnaEdit.textField)
                
        
        
        #set fastQ files
        if self.rnaEdit.params.paired ==True:
            self.fastqFile1=self.rnaEdit.fastqFiles[0]
            self.fastqFile2=self.rnaEdit.fastqFiles[1]
        elif self.rnaEdit.params.paired==False:
            self.fastqFile = self.rnaEdit.fastqFiles[0]


        self.printAttributes()
        
        #self.checkDependencies()
    
    def printAttributes(self):
        print
        Helper.info("*** MAP READS WITH FOLLOWING ATTRIBUTES ***", self.rnaEdit.logFile,self.rnaEdit.textField) 
        if self.rnaEdit.params.paired:
            Helper.info("\t FastQ-File_1: " + self.fastqFile1,self.rnaEdit.logFile,self.rnaEdit.textField)
            Helper.info("\t FastQ-File_2: " + self.fastqFile2,self.rnaEdit.logFile,self.rnaEdit.textField)
        else:
            Helper.info("\t FastQ-File: " + self.fastqFile,self.rnaEdit.logFile,self.rnaEdit.textField)
        Helper.info("\t outfilePrefix:" + self.rnaEdit.params.output,self.rnaEdit.logFile,self.rnaEdit.textField)
        Helper.info("\t refGenome:" + self.rnaEdit.params.refGenome,self.rnaEdit.logFile,self.rnaEdit.textField)
        Helper.info("\t dbsnp:" + self.rnaEdit.params.dbsnp,self.rnaEdit.logFile,self.rnaEdit.textField)
        Helper.info("\t sourceDir:" + self.rnaEdit.params.sourceDir,self.rnaEdit.logFile,self.rnaEdit.textField)
        Helper.info("\t threads:" + self.rnaEdit.params.threads,self.rnaEdit.logFile,self.rnaEdit.textField)
        Helper.info("\t maxDiff:" + self.rnaEdit.params.maxDiff,self.rnaEdit.logFile,self.rnaEdit.textField)
        Helper.info("\t seedDiff:" + self.rnaEdit.params.seedDiff,self.rnaEdit.logFile,self.rnaEdit.textField)
        Helper.info("\t paired:" + str(self.rnaEdit.params.paired),self.rnaEdit.logFile,self.rnaEdit.textField)
        Helper.info("\t keepTemp:" + str(self.rnaEdit.params.keepTemp),self.rnaEdit.logFile,self.rnaEdit.textField)
        Helper.info("\t overwrite:" + str(self.rnaEdit.params.overwrite),self.rnaEdit.logFile,self.rnaEdit.textField)
        Helper.info("",self.rnaEdit.logFile,self.rnaEdit.textField)
    
        
    def startAnalysis(self):   
        recaledBamFile=self.rnaEdit.params.output+".realigned.marked.recalibrated.bam"
        if os.path.isfile(recaledBamFile):
            Helper.info("* * * [Skipping] Mapping result File already exists * * *",self.rnaEdit.logFile,self.rnaEdit.textField)
            self.rnaEdit.logFile.flush()
            return recaledBamFile
        
        
        if self.rnaEdit.params.paired == True:  #For paired end sequencing
            #Align first Fastq Reads to the Genome
            saiFile1=self.rnaEdit.params.output+"_1.sai"
            cmd = [self.rnaEdit.params.sourceDir+"bwa", "aln" , "-t",self.rnaEdit.params.threads, "-n", self.rnaEdit.params.maxDiff , "-k", self.rnaEdit.params.seedDiff, self.rnaEdit.params.refGenome, self.fastqFile1]
            Helper.proceedCommand("Align first Reads with BWA", cmd, self.fastqFile1, saiFile1, self.rnaEdit)
            
            #Align second Fastq Reads to the Genome
            saiFile2=self.rnaEdit.params.output+"_2.sai"
            cmd = [self.rnaEdit.params.sourceDir+"bwa", "aln" , "-t",self.rnaEdit.params.threads, "-n", self.rnaEdit.params.maxDiff , "-k", self.rnaEdit.params.seedDiff, self.rnaEdit.params.refGenome, self.fastqFile2]
            Helper.proceedCommand("Align second Reads with BWA", cmd, self.fastqFile2, saiFile2, self.rnaEdit)
        
            #convert sai to sam
            samFile=self.rnaEdit.params.output+".sam"
            cmd = [self.rnaEdit.params.sourceDir + "bwa", "sampe", "-r", "@RG\tID:A\tLB:A\tSM:A\tPL:ILLUMINA\tPU:HiSEQ2000", self.rnaEdit.params.refGenome, saiFile1, saiFile2, self.fastqFile1, self.fastqFile2]
            Helper.proceedCommand("convert sai to sam", cmd, saiFile1, samFile, self.rnaEdit)
        elif self.rnaEdit.params.paired == False:  #For single end sequencing
            #Align Fastq Reads to the Genome
            saiFile=self.rnaEdit.params.output+".sai"
            cmd = [self.rnaEdit.params.sourceDir+"bwa", "aln" , "-t",self.rnaEdit.params.threads, "-n", self.rnaEdit.params.maxDiff , "-k", self.rnaEdit.params.seedDiff, self.rnaEdit.params.refGenome, self.fastqFile]
            Helper.proceedCommand("Align Reads with BWA", cmd, self.fastqFile, saiFile, self.rnaEdit)
            
            #convert sai to sam
            samFile=self.rnaEdit.params.output+".sam"

            cmd = [self.rnaEdit.params.sourceDir + "bwa", "samse", "-r", "@RG\tID:A\tLB:A\tSM:A\tPL:ILLUMINA\tPU:HiSEQ2000", self.rnaEdit.params.refGenome, saiFile, self.fastqFile]
            Helper.proceedCommand("convert sai to sam", cmd, saiFile, samFile, self.rnaEdit)
        
        #convert sam to bam
        bamFile=self.rnaEdit.params.output+".bam"
        cmd=["java", "-Xmx4G", "-jar", self.rnaEdit.params.sourceDir + "picard-tools/SortSam.jar", "INPUT=" + samFile, "OUTPUT=" + bamFile, "SO=coordinate", "VALIDATION_STRINGENCY=LENIENT", "CREATE_INDEX=true"]
        Helper.proceedCommand("convert sam to bam", cmd, samFile, bamFile, self.rnaEdit)
        
        #return bamFile
        
        
        #run Alignement with tophat
        """
        bamFile=self.rnaEdit.params.output+"/accepted_hits.bam"
        cmd=[self.rnaEdit.params.sourceDir + "tophat/tophat2", "--no-coverage-search","--keep-fasta-order", "-p", "12", "--rg-id", "A","--rg-sample","A","--rg-library","illumina","--rg-platform-unit","HiSeq", "-o", self.rnaEdit.params.output, self.rnaEdit.params.refGenome, self.fastqFile ]
        print cmd
        Helper.proceedCommand("Map reads with tophat", cmd, self.rnaEdit.params.fastqFile, bamFile, self.rnaEdit.)
        """
        
        #sort bam
        #sortBamFile=self.rnaEdit.params.output+".bam"
        #cmd=["java", "-Xmx4G", "-jar", self.rnaEdit.params.sourceDir + "picard-tools/SortSam.jar", "INPUT=" + bamFile, "OUTPUT=" + sortBamFile, "SO=coordinate", "VALIDATION_STRINGENCY=LENIENT", "CREATE_INDEX=true"]
        #Helper.proceedCommand("sortq bam", cmd, bamFile, sortBamFile, self.rnaEdit)
        
        #Add read group ONLY NEEDED WHEN MAPPED WITH TOPHAT
        #rgFile=self.rnaEdit.params.output+".bam"
        #cmd=["java", "-Xmx4G", "-jar", self.rnaEdit.params.sourceDir + "picard-tools/AddOrReplaceReadGroups.jar", "INPUT=" + bamFile, "OUTPUT=" + rgFile, "SO=coordinate", "VALIDATION_STRINGENCY=LENIENT", "CREATE_INDEX=true", "ID=A", "LB=A", "SM=A", "PL=illumina", "PU=HiSeq2000", "SM=A"]
        #Helper.proceedCommand("Add read Groups", cmd, bamFile, rgFile, self.rnaEdit)
        
        
        #Identify Target Regions for realignment
        intervalFile=self.rnaEdit.params.output+".indels.intervals"
        cmd=["java","-Xmx16G","-jar",self.rnaEdit.params.sourceDir + "GATK/GenomeAnalysisTK.jar", "-nt",self.rnaEdit.params.threads, "-T", "RealignerTargetCreator", "-R", self.rnaEdit.params.refGenome, "-I", bamFile, "-o", intervalFile,"-l", "ERROR"]
        Helper.proceedCommand("Identify Target Regions for realignment", cmd, bamFile, intervalFile, self.rnaEdit)
        
        #Proceed Realignement
        realignedFile=self.rnaEdit.params.output+".realigned.bam"
        cmd=["java","-Xmx16G","-jar",self.rnaEdit.params.sourceDir + "GATK/GenomeAnalysisTK.jar", "-T", "IndelRealigner", "-R", self.rnaEdit.params.refGenome, "-I", bamFile, "-l", "ERROR", "-targetIntervals", intervalFile, "-o", realignedFile]
        Helper.proceedCommand("Proceed Realignement", cmd, intervalFile, realignedFile, self.rnaEdit)
        
        #mark PCR duplicates
        markedFile=self.rnaEdit.params.output+".realigned.marked.bam"
        cmd=["java","-Xmx16G","-jar",self.rnaEdit.params.sourceDir + "picard-tools/MarkDuplicates.jar","INPUT=" + realignedFile, "OUTPUT=" + markedFile, "METRICS_FILE="+self.rnaEdit.params.output+".pcr.metrics", "VALIDATION_STRINGENCY=LENIENT", "CREATE_INDEX=true"]
        Helper.proceedCommand("mark PCR duplicates", cmd, realignedFile, markedFile, self.rnaEdit)
        
        #Find Quality Score recalibration spots
        recalFile=self.rnaEdit.params.output+".recalSpots.grp"
        cmd=["java","-Xmx16G","-jar",self.rnaEdit.params.sourceDir + "GATK/GenomeAnalysisTK.jar", "-T", "BaseRecalibrator", "-l", "ERROR", "-R", self.rnaEdit.params.refGenome, "-knownSites", self.rnaEdit.params.dbsnp, "-I", markedFile, "-cov", "CycleCovariate", "-cov", "ContextCovariate", "-o", recalFile]
        Helper.proceedCommand("Find Quality Score recalibration spots", cmd, markedFile, recalFile, self.rnaEdit)
        
        #proceed Quality Score recalibration
        cmd=["java","-Xmx16G","-jar",self.rnaEdit.params.sourceDir + "GATK/GenomeAnalysisTK.jar", "-T", "PrintReads","-l", "ERROR", "-R", self.rnaEdit.params.refGenome, "-I", markedFile, "-BQSR", recalFile, "-o",recaledBamFile]
        Helper.proceedCommand("Proceed Quality Score recalibration", cmd, recalFile, recaledBamFile, self.rnaEdit)
        
        return recaledBamFile
        
    def __del__(self):
        if self.rnaEdit.params.keepTemp==False:
            if os.path.isfile(self.rnaEdit.params.output+".sai"):
                os.remove(self.rnaEdit.params.output+".sai")
            if os.path.isfile(self.rnaEdit.params.output+".sam"):
                os.remove(self.rnaEdit.params.output+".sam")
            if os.path.isfile(self.rnaEdit.params.output+".bam"):
                os.remove(self.rnaEdit.params.output+".bam")
            if os.path.isfile(self.rnaEdit.params.output+".bai"):
                os.remove(self.rnaEdit.params.output+".bai")
            if os.path.isfile(self.rnaEdit.params.output+".indels.intervals"):
                os.remove(self.rnaEdit.params.output+".indels.intervals")
            if os.path.isfile(self.rnaEdit.params.output+".realigned.bam"):
                os.remove(self.rnaEdit.params.output+".realigned.bam")
            if os.path.isfile(self.rnaEdit.params.output+".realigned.bai"):
                os.remove(self.rnaEdit.params.output+".realigned.bai")
            if os.path.isfile(self.rnaEdit.params.output+".realigned.marked.bam"):
                os.remove(self.rnaEdit.params.output+".realigned.marked.bam")
            if os.path.isfile(self.rnaEdit.params.output+".realigned.marked.bai"):
                os.remove(self.rnaEdit.params.output+".realigned.marked.bai")
            if os.path.isfile(self.rnaEdit.params.output+".recalSpots.grp"):
                os.remove(self.rnaEdit.params.output+".recalSpots.grp")
            #os.remove(self.outfilePrefix+".realigned.marked.recalibrated.bam")
    


def checkDependencies(args):
    '''
    Checks the existence of the necessary packages and tools
    :param sourceDir: folder which contains all the software
    '''
    Helper.newline(1)
    Helper.info("CHECK DEPENDENCIES")
    
    #check if all tools are there
    if not os.path.isfile(args.sourceDir+"bwa"):
        Helper.error("BWA not found in %s" % args.sourceDir)
    if not os.path.isfile(args.sourceDir+"picard-tools/SortSam.jar"):
        Helper.error("SortSam.jar not found in %s" % args.sourceDir+"picard-tools")
    if not os.path.isfile(args.sourceDir+"picard-tools/MarkDuplicates.jar"):
        Helper.error("MarkDuplicates.jar not found in %s" % args.sourceDir+"picard-tools")
    if not os.path.isfile(args.sourceDir+"GATK/GenomeAnalysisTK.jar"):
        Helper.error("GenomeAnalysisTK.jar not found in %s" % args.sourceDir+"GATK/")
    if not os.path.isfile(args.sourceDir+"bedtools/fastaFromBed"):
        Helper.error("fastaFromBed not found in %s" % args.sourceDir+"bedtools/")
    if not os.path.isfile(args.sourceDir+"samtools"):
        Helper.error("samtools not found in %s" % args.sourceDir)
    if not os.system("java -version")==0:
        Helper.error("Java could not be found, Please install java")
    
    #check if all files are there
    if not os.path.isfile(args.RefGenome):
        Helper.error("Could not find Reference Genome in %s: " % args.RefGenome)
    # Files for BWA
    if not os.path.isfile(args.RefGenome+".amb"):
        Helper.error("Could not find %s.amb" % args.RefGenome)
        Helper.error("run: 'bwa index %s' to create it" % args.RefGenome)
    if not os.path.isfile(args.RefGenome+".ann"):
        Helper.error("Could not find %s.ann" % args.RefGenome)
        Helper.error("run: 'bwa index %s' to create it" % args.RefGenome)
    if not os.path.isfile(args.RefGenome+".bwt"):
        Helper.error("Could not find %s.bwt" % args.RefGenome)
        Helper.error("run: 'bwa index %s' to create it" % args.RefGenome)
    if not os.path.isfile(args.RefGenome+".pac"):
        Helper.error("Could not find %s.pac" % args.RefGenome)
        Helper.error("run: 'bwa index %s' to create it" % args.RefGenome)
    if not os.path.isfile(args.RefGenome+".sa"):
        Helper.error("Could not find %s.sa" % args.RefGenome)
        Helper.error("run: 'bwa index %s' to create it" % args.RefGenome)
    
    #Files for GATK
    if not os.path.isfile(args.RefGenome.replace(".fastq",".dict")):
        Helper.error("Could not find %s" % args.RefGenome.replace(".fastq",".dict"))
        Helper.error("run: 'java -jar %s/picard-tools/CreateSequenceDictionary.jar R=%s  O= %s.dict' to create it" % (args.sourceDir,args.RefGenome,args.RefGenome.replace(".fastq",".dict")))
    if not os.path.isfile(args.RefGenome+".fai"):
        Helper.error("Could not find %s.fai" % args.RefGenome)
        Helper.error("run: 'samtools faidx %s' to create it" % args.RefGenome)

    #SNP databases
    if not os.path.isfile(args.dbsnp):
        Helper.error("Could not find %s: " % args.dbsnp)


'''
when the script is called directly
''' 
if __name__ == '__main__':
    #parse command line arguments and set defaults
    parser = argparse.ArgumentParser(description='map FastQ Files to the given genome and realigns the reads for SNP-calling.')
    parser.add_argument('-i', '--input', metavar='Fastq-File', type='+', help='Input fastq files (maximum two for paire-end-sequencing)', required=True)
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
    checkDependencies(args)
    
    mapFastQ=MapFastq(args.input, args.RefGenome.name, args.dbsnp.name,args.output, args.sourceDir, args.threads, args.maxDiff, args.seedDiff, args.paired, args.keepTemp, args.overwrite)
    mapFastQ.startAnalysis()
    del mapFastQ