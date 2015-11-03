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
import traceback
from PyQt4 import QtGui, QtCore


import gc


class RnaEdit(QtCore.QThread):

    def __init__(self, fastqFiles, params, textField):
        QtCore.QThread.__init__(self)
        if isinstance(params, Parameters):
            self.params = params
        else:
            Helper.error("Params has to be Instance of Parameters")
        if isinstance(textField, QtGui.QTextEdit) or textField==0:
            self.textField=textField
        else:
            Helper.error("textField has to be Instance of QtGui.QTextEdit or 0")
        

        self.fastqFiles=fastqFiles
        
        #hold the running Popen object
        self.runningCommand=False
        self.isTerminated = False
        #check if the input Files are there
        
        #hold basic statistic values of the run
        basicStatDict={}
            
        
        #set directory where the outputFiles should be written to
        if self.params.output=="default":
            if self.fastqFiles[0].endswith("realigned.marked.recalibrated.bam"):
                self.sampleName=fastqFiles[0][fastqFiles[0].rfind("/")+1:fastqFiles[0].rfind(".realigned.marked.recalibrated.bam")]
                self.outdir=fastqFiles[0][0:fastqFiles[0].rfind("/")+1]
            else:
                self.sampleName=fastqFiles[0][fastqFiles[0].rfind("/")+1:fastqFiles[0].rfind(".")]
                # outdir = /path/to/output/rnaEditor/samplename/
                self.outdir=fastqFiles[0][0:fastqFiles[0].rfind("/")+1]+"rnaEditor/"+self.sampleName+"/"
            
            #output=/path/to/output/rnaEditor/samplename/samplename
            self.params.output=self.outdir+self.sampleName
            if not os.path.exists(self.outdir):
                os.makedirs(self.outdir)
            
            #create folder for html output
            if not os.path.exists(self.outdir+"/html"):
                os.makedirs(self.outdir+"/html")
        
        
        self.checkDependencies()
        
        #check if the input Files are there
        self.printParameters()
        
    def run(self):
        try:
            self.startAnalysis()
        except Exception:
            Helper.error("RnaEditor Failed",self.rnaEdit.logFile,self.rnaEdit.textField)
        
        
    def startAnalysis(self):
        """
        START MAPPING
        """
        if self.fastqFiles[0].endswith(".bam"):
            Helper.info("Bam File given. Skip mapping", self.logFile, self.textField)
            self.mapFastQ=None
            mapResultFile=self.fastqFiles[0]
        else:
            self.mapFastQ=MapFastq(self)
            mapResultFile=self.mapFastQ.startAnalysis()

        """
        START CALLING EDITING SITES
        """
        self.callEditSites=CallEditingSites(mapResultFile,self)
        result = self.callEditSites.startAnalysis()
        
        #finished
        self.isTerminated=True
        
        self.emit(QtCore.SIGNAL("taskDone"), self.params.output)
        
        Helper.status("rnaEditor Finished with %s" % self.params.output, self.logFile, self.textField)
        self.cleanUp()
    
    def stopSafely(self):
        self.quit()
        Helper.info("Analysis was stopped by User", self.logFile, self.textField)
    
    def stopImmediately(self):
        if hasattr(self, 'callEditSites'):
            self.callEditSites.cleanUp()
        self.isTerminated=True
        
        if self.runningCommand != False:
            self.runningCommand.kill()
        else:
            self.terminate()
            self.wait()
        Helper.error("Analysis was terminated by User", self.logFile, self.textField)
         
    def cleanUp(self):
        #print "deleteAssay " + str(self)
        if self.runningCommand != False:
            self.runningCommand.kill()

        try:
            if self.mapFastQ!=None:
                self.mapFastQ.cleanUp()
            del self.mapFastQ
        except AttributeError:
            Helper.error("could not delete MapFastQ instance", self.logFile, self.textField)
        try:
            self.callEditSites.cleanUp()
            del self.callEditSites
        except AttributeError:
            Helper.error("could not delete RnaEdit instance", self.logFile, self.textField)
        
    def checkDependencies(self):
        """checks if all files are there
        if all programs are installed properly and if the output directory is writable"""
        try:
            self.logFile=open(self.params.output + ".log","w+")
        except IOError:
            Helper.error("Cannot open Log File", textField=self.textField)

        if type(self.fastqFiles) == list:
            self.fastqFiles=self.fastqFiles
        elif type(self.fastqFile) == str:
            self.fastqFiles=[self.fastqFiles]
        else:
            Helper.error("FastQ File has wrong variable type",self.logFile,self.textField)
        
        for file in self.fastqFiles:
            if not os.path.isfile(file):
                Helper.error("Could not find: %s" %file,self.logFile,self.textField)
            
        '''
        Checks the existence of the necessary packages and tools
        :param sourceDir: folder which contains all the software
        '''
        Helper.newline(1)
        Helper.info("CHECK DEPENDENCIES",self.logFile,self.textField)
        
        #check if all tools are there
        if not os.path.isfile(self.params.sourceDir+"bwa"):
            Helper.error("BWA not found in %s" % self.params.sourceDir,self.logFile,self.textField)
        if not os.path.isfile(self.params.sourceDir+"picard-tools/SortSam.jar"):
            Helper.error("SortSam.jar not found in %s" % self.params.sourceDir+"picard-tools",self.logFile,self.textField)
        if not os.path.isfile(self.params.sourceDir+"picard-tools/MarkDuplicates.jar"):
            Helper.error("MarkDuplicates.jar not found in %s" % self.params.sourceDir+"picard-tools",self.logFile,self.textField)
        if not os.path.isfile(self.params.sourceDir+"GATK/GenomeAnalysisTK.jar"):
            Helper.error("GenomeAnalysisTK.jar not found in %s" % self.params.sourceDir+"GATK/",self.logFile,self.textField)
        if not os.path.isfile(self.params.sourceDir+"bedtools/fastaFromBed"):
            Helper.error("fastaFromBed not found in %s" % self.params.sourceDir+"bedtools/",self.logFile,self.textField)
        if not os.path.isfile(self.params.sourceDir+"blat"):
            Helper.error("blat not found in %s" % self.params.sourceDir,self.logFile,self.textField)
        if not os.path.isfile(self.params.sourceDir+"samtools"):
            Helper.error("samtools not found in %s" % self.params.sourceDir,self.logFile,self.textField)
        if not os.system("java -version")==0:
            Helper.error("Java could not be found, Please install java",self.logFile,self.textField)
        
        
        
        #check if all files are there
        if not os.path.isfile(self.params.refGenome):
            Helper.error("Could not find Reference Genome in %s: " % self.params.refGenome,self.logFile,self.textField)
        
        # Files for BWA
        if not os.path.isfile(self.params.refGenome+".amb"):
            Helper.error("Could not find %s.amb" % self.params.refGenome,self.logFile,self.textField)
            Helper.error("run: 'bwa index %s' to create it" % self.params.refGenome,self.logFile,self.textField)
        if not os.path.isfile(self.params.refGenome+".ann"):
            Helper.error("Could not find %s.ann" % self.params.refGenome,self.logFile,self.textField)
            Helper.error("run: 'bwa index %s' to create it" % self.params.refGenome,self.logFile,self.textField)
        if not os.path.isfile(self.params.refGenome+".bwt"):
            Helper.error("Could not find %s.bwt" % self.params.refGenome,self.logFile,self.textField)
            Helper.error("run: 'bwa index %s' to create it" % self.params.refGenome,self.logFile,self.textField)
        if not os.path.isfile(self.params.refGenome+".pac"):
            Helper.error("Could not find %s.pac" % self.params.refGenome,self.logFile,self.textField)
            Helper.error("run: 'bwa index %s' to create it" % self.params.refGenome,self.logFile,self.textField)
        if not os.path.isfile(self.params.refGenome+".sa"):
            Helper.error("Could not find %s.sa" % self.params.refGenome,self.logFile,self.textField)
            Helper.error("run: 'bwa index %s' to create it" % self.params.refGenome,self.logFile,self.textField)

        
        #Files for GATK
        
        if not os.path.isfile(self.params.refGenome.replace(".fastq",".dict")):
            Helper.error("Could not find %s" % self.params.refGenome.replace(".fastq",".dict"),self.logFile,self.textField)
            Helper.error("run: 'java -jar %s/picard-tools/CreateSequenceDictionary.jar R=%s  O= %s.dict' to create it" % (self.params.sourceDir,self.params.refGenome,self.params.refGenome),self.logFile,self.textField)
        if not os.path.isfile(self.params.refGenome+".fai"):
            Helper.error("Could not find %s.sai" % self.params.refGenome,self.logFile,self.textField)
            Helper.error("run: 'samtools faidx %s' to create it" % self.params.refGenome,self.logFile,self.textField)
    
        #SNP databases
        if not os.path.isfile(self.params.dbsnp):
            Helper.error("Could not find %s: " % self.params.dbsnp,self.logFile,self.textField)
        if not os.path.isfile(self.params.hapmap):
            Helper.error("Could not find %s: " % self.params.hapmap,self.logFile,self.textField)
        if not os.path.isfile(self.params.omni):
            Helper.error("Could not find %s: " % self.params.omni,self.logFile,self.textField)
        if not os.path.isfile(self.params.esp):
            Helper.error("Could not find %s: " % self.params.esp,self.logFile,self.textField)
            
        #region Files
        if not os.path.isfile(self.params.aluRegions):
            Helper.error("Could not find %s: " % self.params.aluRegions,self.logFile,self.textField)
            
        if not os.path.isfile(self.params.gtfFile):
            Helper.error("Could not find %s: " % self.params.gtfFile,self.logFile,self.textField)

        
        Helper.info("Dependencies satisfied", self.logFile, self.textField)
        
        #check if display can activate 

    def printParameters(self):

        Helper.info("*** Start RnaEditor with: ***", self.logFile,self.textField) 
        if self.params.paired:
            Helper.info("\t FastQ-File_1: " + self.fastqFiles[0],self.logFile,self.textField)
            Helper.info("\t FastQ-File_2: " + self.fastqFiles[1],self.logFile,self.textField)
        else:
            Helper.info("\t FastQ-File: " + self.fastqFiles[0],self.logFile,self.textField)
        Helper.info("\t outfilePrefix:" + self.params.output,self.logFile,self.textField)
        Helper.info("\t refGenome:" + self.params.refGenome,self.logFile,self.textField)
        Helper.info("\t dbsnp:" + self.params.dbsnp,self.logFile,self.textField)
        Helper.info("\t sourceDir:" + self.params.sourceDir,self.logFile,self.textField)
        Helper.info("\t threads:" + self.params.threads,self.logFile,self.textField)
        Helper.info("\t maxDiff:" + self.params.maxDiff,self.logFile,self.textField)
        Helper.info("\t seedDiff:" + self.params.seedDiff,self.logFile,self.textField)
        Helper.info("\t paired:" + str(self.params.paired),self.logFile,self.textField)
        Helper.info("\t keepTemp:" + str(self.params.keepTemp),self.logFile,self.textField)
        Helper.info("\t overwrite:" + str(self.params.overwrite),self.logFile,self.textField)
        Helper.info("",self.logFile,self.textField)

        def createDiagrams(self):
            pass



if __name__ == '__main__':
    
    
    parser = argparse.ArgumentParser(description='map FastQ Files to the given genome and realigns the reads for SNP-calling.',)
    parser.add_argument('-i', '--input', metavar='Fastq-Files',nargs='+', type=str, help='Input fastq files (maximum two for paire-end-sequencing)', required=True)
    parser.add_argument('-c', '--conf', metavar='Configuration File', type=str, help='Configuration File used to read Parameters for RnaEditor', required=True)
    
    
    args = parser.parse_args()
    
    
    parameters = Parameters(args.conf) 
    edit=RnaEdit(args.input,parameters,0)
    
    
    edit.start()
    edit.wait()
    Helper.createDiagramms(edit.params.output)
    Helper.printResultHtml(edit.params.output)
    del edit
    
    
    
else:
    pass    
    
    