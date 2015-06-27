'''
Created on May 22, 2013

@author: david
'''

from datetime import datetime, date, time
import argparse, sys, os, subprocess, errno
from collections import defaultdict
from multiprocessing import Queue
from time import sleep


class Parameters():
    '''
    Reads and saves the default values from the configuration File 'configurtion.txt'
    '''
    refGenome=""
    dbSNP=""
    hapmap=""
    omni=""
    esp=""
    aluRegions=""
    gtfFile = ""
    output = ""
    binary = ""
    maxDiff = 0.04
    seedDiff = 2
    paired = False
    standCall = 0
    standEmit = 0
    edgeDistance = 3
    keepTemp = False
    overwrite = False
    threads = 5
    
    
    @staticmethod
    def readDefaults(file="configuration.txt"):
        confFile = open(file)
        for line in confFile:
            
            if line.startswith("#"):
                continue
            if line == "\n":
                continue
            
            line=line.rstrip()
            id, value = line.split("=")
            id=id.strip()
            value=value.strip()

            if id=="refGenome":
                Parameters.refGenome=value
            elif id=="dbSNP":
                Parameters.dbSNP=value
            elif id=="hapmap":
                Parameters.hapmap=value
            elif id=="omni":
                Parameters.omni=value
            elif id=="esp":
                Parameters.esp=value
            elif id == "aluRegions":
                Parameters.aluRegions=value
            elif id == "gtfFile":
                Parameters.gtfFile=value
            elif id == "output":
                Parameters.output=value
            elif id == "binary":
                Parameters.binary=value
            elif id == "maxDiff":
                Parameters.maxDiff=float(value)
            elif id == "seedDiff":
                Parameters.seedDiff=int(value)
            elif id == "paired":
                #Parameters.paired=float(value)
                if str(value).lower() in ("yes", "y", "true",  "t", "1"): Parameters.paired = True
                if str(value).lower() in ("no",  "n", "false", "f", "0", "0.0", "", "none", "[]", "{}"): Parameters.paired=False
            elif id == "standCall":
                Parameters.standCall=int(value)
            elif id == "standEmit":
                Parameters.standEmit=int(value)    
            elif id == "edgeDistance":
                Parameters.edgeDistance=int(value)
            elif id == "threads":
                Parameters.threads=int(value)
            elif id == "keepTemp":
                #Parameters.paired=float(value)
                if str(value).lower() in ("yes", "y", "true",  "t", "1"): Parameters.keepTemp = True
                if str(value).lower() in ("no",  "n", "false", "f", "0", "0.0", "", "none", "[]", "{}"): Parameters.keepTemp=False
            elif id == "overwrite":
                #Parameters.paired=float(value)
                if str(value).lower() in ("yes", "y", "true",  "t", "1"): Parameters.overwrite = True
                if str(value).lower() in ("no",  "n", "false", "f", "0", "0.0", "", "none", "[]", "{}"): Parameters.overwrite=False
           
class Helper():
    '''
    Helpfunctions
    '''
    
    '''
    check if given directory is a readable directory and give the right data type 
    '''
    
    prefix = "*** "
    praefix = " ***"
    
    #dummy element is added to the array to avoid 0/1 problem from the Tab array and these arrays
    #otherwise i had to add -1 every time i want to access the following arrays
    runningAssaysThreads=["dummy"]
    assays=["dummy"] #RNAeditor objects
    
    assayIsRunning=['dummy'] #stores TRUE if the assay at that position was started successfully
    runningAssaysTabs=["dummy"] #tabs from the user interface
    runningCommand=["dummy"] #saves subprocess.POPEN objects to kill them later of False if no process is running
    index=0
    assayCount=0
    
    @staticmethod
    def getSampleName(fq):
                #get name from input File
        if fq.endswith(".fastq"):
            sampleName = fq[fq.rfind("/")+1:fq.rfind(".fastq")]
        elif fq.endswith(".fq"):
            sampleName = fq[fq.rfind("/")+1:fq.rfind(".fq")]
        elif fq.endswith(".bam"):
            sampleName = fq[fq.rfind("/")+1:fq.rfind(".bam")]
        else:
            return None
        return sampleName
    
    @staticmethod
    def readable_dir(prospective_dir):
        if not os.path.isdir(prospective_dir):
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(prospective_dir))
        if os.access(prospective_dir, os.R_OK):
            return prospective_dir
        else:
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))
    
    '''
    return current time
    '''
    @staticmethod
    def getTime():
        curr_time = datetime.now()
        #return "["+curr_time.strftime("%c")+"]"
        return curr_time
    
    """
    converts the inputFile to phred33 Quality and writes it into the ourdir
    """
    @staticmethod
    def convertPhred64toPhred33(self,fastqFile,outFile,logFile,runNumber):
        startTime=Helper.getTime()
        Helper.info("[" + startTime.strftime("%c") + "] * * * convert Quality encoding: " + fastqFile[fastqFile.rfind("/")+1:]   + " * * *",logFile,runNumber)
        
        
        if os.path.exists(outFile):
            Helper.info("* * * [Skipping] Result File already exists * * *",logFile,runNumber)
            return outFile
        
        outFile = open(outFile,"w")
        fastqFile=open(fastqFile,"r")
        
        lineNumber = 0 
        for line in fastqFile:
            lineNumber+=1
            if lineNumber%4==0:
                a=[]
                for char in line.rstrip():
                    phredQual=ord(char)-64
                    phredChar=chr(phredQual+33)
                    a.append(phredChar)
                outFile.write("".join(a) + "\n")
            else:
                outFile.write(line)
        outFile.close()
        return outFile.name
    
    """
    chech in the first lines if the quality encoding is phred33
    """
    @staticmethod
    def isPhred33Encoding(inFastqFile,lines,logFile, runNumber):
        fastqFile=open(inFastqFile,"r")
        lineNumber=0
        lines=lines*4
        for line in fastqFile:
            lineNumber+=1
            if lineNumber%4==0:
                for char in line.rstrip():
                    if ord(char)>74 and ord(char)<105:
                        #print line.rstrip()
                        fastqFile.close()
                        return False
                    if ord(char)>105:
                        Helper.error("%s has no valid quality encoding. \n\t Please use a valid FastQ file??" % fastqFile.name,logFile, runNumber)
            if lineNumber > lines:
                fastqFile.close()
                return True
            
        Helper.error("%s has less than %i Sequences. \n These are not enough reads for editing detection!!" % (fastqFile.name,lines),logFile, runNumber)
    
    '''
    run a specific NGS-processing-step on the system
    '''
    @staticmethod
    def proceedCommand(description,cmd,infile,outfile,logFile,overwrite=False,runNumber=0):
        startTime=Helper.getTime()
        Helper.info("[" + startTime.strftime("%c") + "] * * * " + description + " * * *",logFile,runNumber)
        
        
        #check if infile exists
        if not os.path.isfile(infile):
            Helper.error(infile + "does not exist, Error in previous Step",logFile,runNumber)
            #Exception(infile + "does not exist, Error in previous Step")
            #exit(1)
        
        #check if outfile already exists
        
        if not os.path.isfile(outfile) or overwrite==True:
            if outfile == "None":
                resultFile=None
            else:
                resultFile=open(outfile,"w+")
            try:    
                
                #print " ".join(cmd),resultFile,logFile
                
                #retcode = subprocess.call(cmd, stdout=resultFile, stderr=logFile)
                Helper.runningCommand[runNumber] = subprocess.Popen(cmd, stdout=resultFile, stderr=logFile)
                retcode=Helper.runningCommand[runNumber].wait()
                """while retcode==None:
                    #print "check if process is still running"
                    sleep(10)
                    retcode=Helper.runningCommand[runNumber].wait()
                """
                print retcode
                
                #del Helper.runningCommand[runNumber]
                Helper.runningCommand[runNumber]=False
                if retcode != 0:
                    if retcode == -9:
                        Helper.error(description+ " canceled by User!!!",logFile,runNumber)
                    else:
                        Helper.error(description+ " failed!!!",logFile,runNumber)
                    
                    if resultFile!=None:
                        os.remove(resultFile.name)
                    #exit(1)
            except OSError, o:
                if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
                    Helper.error(cmd[0] + " Command not found on this system",logFile,runNumber)
                    if resultFile!=None:
                        os.remove(resultFile.name)
                    #exit(1)
                else:
                    Helper.error(cmd[0] + o.strerror,logFile,runNumber)
                    if resultFile!=None:
                        os.remove(resultFile.name)
                    #exit(1)
            Helper.printTimeDiff(startTime, logFile, runNumber)
        else:
            print "\t [SKIP] File already exist",logFile,runNumber

    """
    return a dictionary whith chromosome as keys and a set of variants as values
    variantDict={chromosome:(variantPos1,variantPos2,....)}
    """
    @staticmethod
    def getPositionDictFromVcfFile(vcfFile,runNumber):
        variantFile=open(vcfFile)
        variantDict=defaultdict(set)
        Helper.info("reading Variants from %s" % vcfFile,runNumber)
        for line in variantFile:
            #skip comments
            if line.startswith("#"): continue
            line=line.split("\t")
            chromosome,position,ref,alt = line[0],line[1],line[3],line[4]
            variantDict[chromosome].add(position)
        return variantDict
    
    @staticmethod
    def removeVariantsAFromVariantsB(variantsDictA,variantsDictB):
        
        
        if type(variantsDictA) is str:
            variantsDictA = Helper.returnVariantDictFromVcfFile(variantsDictA)
        if type(variantsDictB) is str:
            variantsDictB = Helper.returnVariantDictFromVcfFile(variantsDictB)
            
        resultDict = defaultdict(set)
        for variant in variantsDictB.iterkeys():
            if variant not in variantsDictA[chr]:
                resultDict[chr].append(variant)
        return resultDict
       
    @staticmethod
    def returnVariantDictFromVcfFile(vcfFile):
        """
        returns the vcfFile as a two instance dictionary with chromosome as first key and a Tuple of (position,ref,alt) as second key  and  the rest of the vcfLine as a list
        {chr1: {(position1, 'A', 'G'): [dbSNP_id, ' quality', 'filter', 'attributes'], (4, 324, 'dsgdf', 'dsfsd'): [42, 243, 324]}})
        """
        vcfFile=open(vcfFile)
        vcfDict = defaultdict(dict)
        for line in vcfFile:
            #skip comments
            if line.startswith("#"): continue
            line=line.rstrip().split("\t")
            chromosome,position,ref, alt = line[0],line[1], line[3], line[4]
            vcfDict[chromosome][position,ref,alt]=([line[2]]+line[5:])
        return vcfDict



    @staticmethod
    def getCommandOutput(command):
        #print command
        #print os.path.dirname(os.path.abspath(__file__))
        #print os.getcwd()
        return subprocess.check_output(command)
    
    
    
    @staticmethod
    def printTimeDiff(startTime,logFile=None,textField=0):
        duration = Helper.getTime() - startTime
        if textField!=0:
            #currentAssay = Helper.runningAssays[textField] 
            textField.append("\t" + Helper.prefix + "[DONE] Duration [" + str(duration) + "]"  + Helper.praefix)
        if logFile!=None:
            logFile.write("\t" + Helper.prefix + "[DONE] Duration [" + str(duration) + "]"  + Helper.praefix + "\n")
        
        sys.stderr.write("\t" + Helper.prefix + "[DONE] Duration [" + str(duration) + "]"  + Helper.praefix + "\n")
    @staticmethod
    def newline (quantity=1,logFile=None,textField=0):
        if textField!=0:
            #currentAssay = Helper.runningAssays[runNumber] 
            textField.append("\n"*quantity)
        if logFile!=None:
            logFile.write("\n"*quantity)
        sys.stderr.write("\n"*quantity)
    @staticmethod
    def info (message,logFile=None,textField=0):
        if textField!=0:
            #currentAssay = Helper.runningAssays[runNumber] 
            textField.append(Helper.prefix + "INFO:    "  + message + Helper.praefix)
        if logFile!=None:
            logFile.write(Helper.prefix + "INFO:    "  + message + Helper.praefix + "\n")
        sys.stderr.write(Helper.prefix + "INFO:    "  + message + Helper.praefix + "\n")
    @staticmethod
    def warning (message,logFile=None,textField=0):
        if textField!=0:
            textField.append("\n\n" + Helper.prefix + "WARNING:    " + message + Helper.praefix + "\n\n")
        if logFile!=None:
            logFile.write(Helper.prefix + "WARNING:    "  + message + Helper.praefix + "\n")
        sys.stderr.write("\n\n" + Helper.prefix + "WARNING:    " + message + Helper.praefix + "\n\n")
    @staticmethod
    def error (message,logFile=None,textField=0):
        if textField!=0:
            textField.append("\n\n" + Helper.prefix + "ERROR:    "  + message + Helper.praefix + "\n\n")
        if logFile!=None:
            logFile.write(Helper.prefix + "ERROR:    "  + message + Helper.praefix + "\n")
        #sys.stderr.write("\n\n" + Helper.prefix + "ERROR:    " + message + Helper.praefix + "\n\n")
        raise Exception("\n\n" + Helper.prefix + "ERROR:    " + message + Helper.praefix + "\n\n")
    @staticmethod
    def debug (message,logFile=None,textField=0):
        if textField!=0:
            textField.append(Helper.prefix + "DEBUG:    "  + message + Helper.praefix)
        if logFile!=None:
            logFile.write(Helper.prefix + "DEBUG:    "  + message + Helper.praefix + "\n")
        sys.stderr.write(Helper.prefix + message + Helper.praefix + "\n")
    @staticmethod
    def status(message,logFile=None,textField=0):
        if textField!=0:
            #currentAssay = Helper.runningAssays[runNumber] 
            textField.append(Helper.prefix + "STATUS:    "  + message + Helper.praefix)
        if logFile!=None:
            logFile.write(Helper.prefix + "STATUS:    "  + message + Helper.praefix + "\n")
        sys.stdout.write("\r" + Helper.prefix + "STATUS:    "  + message + Helper.praefix + "\n")
        sys.stdout.flush()