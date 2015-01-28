'''
Created on May 22, 2013

@author: david
'''

from datetime import datetime, date, time
import argparse, sys, os, subprocess, errno
from collections import defaultdict




class Helper():
    '''
    Helpfunctions
    '''
    
    '''
    check if given directory is a readable directory and give the right data type 
    '''
    
    prefix = "*** "
    praefix = " ***"
    
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
    def convertPhred64toPhred33(self,fastqFile,outFile,logFile):
        startTime=Helper.getTime()
        logFile.write("[" + startTime.strftime("%c") + "] * * * convert Quality encoding: " + fastqFile[fastqFile.rfind("/")+1:]   + " * * *")
        logFile.flush()
        print "[" + startTime.strftime("%c") + "] * * * convert Quality encoding for " + fastqFile[fastqFile.rfind("/")+1:]   + " * * *"
        
        if os.path.exists(outFile):
            print >> self.logFile, "* * * [Skipping] Result File already exists * * *"
            self.logFile.flush()
            print "* * * [Skipping] Result File already exists * * *"
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
    def isPhred33Encoding(inFastqFile,lines):
        fastqFile=open(inFastqFile,"r")
        lineNumber=0
        lines=lines*4
        for line in fastqFile:
            lineNumber+=1
            if lineNumber%4==0:
                for char in line.rstrip():
                    if ord(char)>74:
                        #print line.rstrip()
                        fastqFile.close()
                        return False
                
            if lineNumber > lines:
                fastqFile.close()
                return True
    
    '''
    run a specific NGS-processing-step on the system
    '''
    @staticmethod
    def proceedCommand(description,cmd,infile,outfile,logFile,overwrite=False):
        startTime=Helper.getTime()
        Helper.info("[" + startTime.strftime("%c") + "] * * * " + description + " * * *")
        
        
        #check if infile exists
        if not os.path.isfile(infile):
            
            print infile + "does not exist, Error in previous Step"
            #Exception(infile + "does not exist, Error in previous Step")
            exit(1)
        
        #check if outfile already exists
        
        if not os.path.isfile(outfile) or overwrite==True:
            if outfile == "None":
                resultFile=None
            else:
                resultFile=open(outfile,"w+")
            try:    
                #os.popen(cmd)
                #retcode = subprocess.call(cmd,shell=True)
                #print cmd,resultFile,logFile
                
                retcode = subprocess.call(cmd, stdout=resultFile, stderr=logFile)
                if retcode != 0:
                    print >> sys.stderr, "Error: " + description + " failed"
                    if resultFile!=None:
                        os.remove(resultFile.name)
                    exit(1)
            except OSError, o:
                if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
                    print >> sys.stderr, "Error: " + cmd[0] + " Command not found on this system"
                    if resultFile!=None:
                        os.remove(resultFile.name)
                    exit(1)
                else:
                    print >> sys.stderr, "Error: " + cmd[0] + o.strerror
                    if resultFile!=None:
                        os.remove(resultFile.name)
                    exit(1)
            duration=Helper.getTime()-startTime
            print >> logFile, "\t[DONE]" + " Duration [" + str(duration) + "]"
            logFile.flush()
            print "\t[DONE]" + " Duration [" + str(duration) + "]"
        else:
            print "\t [SKIP] File already exist"

    """
    return a dictionary whith chromosome as keys and a set of variants as values
    variantDict={chromosome:(variantPos1,variantPos2,....)}
    """
    @staticmethod
    def getPositionDictFromVcfFile(vcfFile):
        variantFile=open(vcfFile)
        variantDict=defaultdict(set)
        Helper.info("reading Variants from %s" % vcfFile)
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
    def printTimeDiff(startTime,logFile=None):
        duration = Helper.getTime() - startTime
        if logFile!=None:
            logFile.write("\t" + Helper.prefix + "[DONE] Duration [" + str(duration) + "]"  + Helper.praefix + "\n")
        
        sys.stderr.write("\t" + Helper.prefix + "[DONE] Duration [" + str(duration) + "]"  + Helper.praefix + "\n")
    @staticmethod
    def newline (quantity=1,logFile=None):
        if logFile!=None:
            logFile.write("\n"*quantity)
        sys.stderr.write("\n"*quantity)
    @staticmethod
    def info (message,logFile=None):
        if logFile!=None:
            logFile.write(Helper.prefix + "INFO:    "  + message + Helper.praefix + "\n")
        sys.stderr.write(Helper.prefix + "INFO:    "  + message + Helper.praefix + "\n")
    @staticmethod
    def warning (message,logFile=None):
        if logFile!=None:
            logFile.write(Helper.prefix + "WARNING:    "  + message + Helper.praefix + "\n")
        sys.stderr.write("\n\n" + Helper.prefix + "WARNING:    " + message + Helper.praefix + "\n\n")
    @staticmethod
    def error (message,logFile=None):
        if logFile!=None:
            logFile.write(Helper.prefix + "ERROR:    "  + message + Helper.praefix + "\n")
        #sys.stderr.write("\n\n" + Helper.prefix + "ERROR:    " + message + Helper.praefix + "\n\n")
        raise Exception("\n\n" + Helper.prefix + "ERROR:    " + message + Helper.praefix + "\n\n")
    @staticmethod
    def debug (message,logFile=None):
        if logFile!=None:
            logFile.write(Helper.prefix + "DEBUG:    "  + message + Helper.praefix + "\n")
        sys.stderr.write(Helper.prefix + message + Helper.praefix + "\n")
    @staticmethod
    def status(message,logFile=None):
        if logFile!=None:
            logFile.write(Helper.prefix + "STATUS:    "  + message + Helper.praefix + "\n")
        sys.stdout.write("\r" + Helper.prefix + "STATUS:    "  + message + Helper.praefix + "\n")
        sys.stdout.flush()