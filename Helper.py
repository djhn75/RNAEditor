'''
Created on May 22, 2013

@author: david
'''


from datetime import datetime, date, time
import argparse, sys, os, subprocess, errno


class Helper():
    '''
    Helpfunctions
    '''
    
    '''
    check if given directory is a readable directory and give the right data type 
    '''
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
    
    '''
    run a specific NGS-processing-step on the system
    '''
    @staticmethod
    def proceedCommand(description,cmd,infile,outfile,logFile,overwrite=False):
        startTime=Helper.getTime()
        print >> logFile, "[" + startTime.strftime("%c") + "] * * * " + description + " * * *"
        logFile.flush()
        print "[" + startTime.strftime("%c") + "] * * * " + description + " * * *"
        
        #check if infile exists
        if not os.path.exists(infile):
            
            print infile + "does not exist, Error in previous Step"
            #Exception(infile + "does not exist, Error in previous Step")
            exit(1)
        
        #check if outfile already exists
        
        if not os.path.exists(outfile) or overwrite==True:
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
                    os.remove(resultFile.name)
                    exit(1)
            except OSError, o:
                if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
                    print >> sys.stderr, "Error: " + cmd[0] + " Command not found on this system"
                    os.remove(resultFile.name)
                    exit(1)
                else:
                    print >> sys.stderr, "Error: " + cmd[0] + o.strerror
                    os.remove(resultFile.name)
                    exit(1) 
            duration=Helper.getTime()-startTime
            print >> logFile, "\t[DONE]" + " Duration [" + str(duration) + "]"
            logFile.flush()
            print "\t[DONE]" + " Duration [" + str(duration) + "]"
        else:
            print "\t [SKIP] File already exist"

    @staticmethod
    def getCommandOutput(command):
        return subprocess.check_output(command, shell=True)