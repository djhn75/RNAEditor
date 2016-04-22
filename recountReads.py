'''
Created on 02.03.2015
This script recalculates the number of read per base and joins the .
If there is an 

@author: david
'''

import argparse, os, re, sys

import pysam

from Helper import Helper


parser = argparse.ArgumentParser(description='Merges the GVF Files and recalculates the base Counts after RnaEditor is finished.')
parser.add_argument('-f', '--files', metavar='N', type=str, nargs='+', help='the list of files')
parser.add_argument('-b', '--bams', metavar='N', type=str, nargs='+', help='the list of bam files')
parser.add_argument('-t', '--top', metavar='N', type=str, nargs="+", help='list of header names (space separated)')
parser.add_argument('-o', '--outFile', metavar='output File', type=str,help='Output File', default="baseCounts_combined.txt")
parser.add_argument('-c', '--columns', metavar='N', type=int, nargs='+', help='columns to keep (space separated)',default=[2],)
parser.add_argument('-k', '--keys', metavar='N', nargs='+', type=int, help='columnnumber on which to join',default=[1])
parser.add_argument('-d', '--delimiter', metavar='N', type=str, help='delimiter', default="\t")
parser.add_argument('-e', '--empty', metavar='N', type=str, help='Sign for empty Values', default="--")
args = parser.parse_args()

startTime = Helper.getTime()


def fillDicts(files,columns,keys):
    '''
        creates the table and fills the set of keys
    '''
    fileNumber=len(files)
    fileCounter=0
    keySet=()
    fileCounter=0
    for file in files: #loop through all files
        i=0
        Helper.info("Get information from %s" % file)
        file = open(file)
        
        for line in file: #loop through current file
            line = line.split()
            keyTuple=()
            for k in keys:
                keyTuple=keyTuple+(line[k-1],)
            
            value=[]
            for column in columns: #get the needed values
                try:
                    value.append(line[column-1])
                except IndexError:
                    raise ValueError("Not enough rows in line: %s in file %s" % (" ".join(line),file.name))
            
            if keyTuple in keySet:
                #currentDefaultList=idDict[keyTuple]
                #currentDefaultList[fileCounter]=value
                #idDict[keyTuple]=currentDefaultList
                idDict[keyTuple][fileCounter]=value #replace filecounter List with values from current File
            else:
                currentDefaultList=[["--"]*len(columns)]*len(files) #create default list, with all values empty
                currentDefaultList[fileCounter]=value
                idDict[keyTuple]=currentDefaultList
                keySet=keySet+(keyTuple,)
            
            i+=1
            if i % 1000 == 0:
                Helper.status("%s lines parsed" % i)
        fileCounter+=1
    return idDict,keySet

def getBaseCount(reads, varPos):
    '''
    
    :param reads: 
    :param varPos:
    '''
    '''
        returns the baseCount for the 
    '''
    baseCount = {'A':0,'C':0,'G':0,'T':0}
    for read in reads:
        
        readPos=0
        mmReadPos=0
        startPos = read.pos
        try:
            cigarNums=re.split("[MIDNSHP]", read.cigarstring)[:-1]
            cigarLetters=re.split("[0-9]+",read.cigarstring)[1:]
        except (TypeError):
            continue    #for unmapped reads the cigarstring is empty
                        #to avoid a query for unmapped reads all the 
                        #time the error is catched and the read will be skipped
            #raise TypeError("Invalid Cigar String %s" % read.cigarstring)
        
        for i in range(len(cigarLetters)): #parse over single read
            if cigarLetters[i] in {"I","S","H"}: #Insertion, Soft Clipping and Hard Clipping
                readPos = readPos + int(cigarNums[i])
            elif cigarLetters[i] in {"D","N"}: #Deletions and skipped Regions
                startPos = startPos + int(cigarNums[i])
            elif cigarLetters[i] in {"M"}: #Matches
                for j in range(int(cigarNums[i])):
                    if startPos == varPos:
                        mmReadPos = readPos
                        mmReadBase= read.seq[mmReadPos]
                        try:
                            baseCount[mmReadBase]+=1 #increase number for the base at the mm pos
                        except (KeyError):
                            sys.stderr.write("unknown Base %s \n" % mmReadBase)
                            
                    readPos += 1
                    startPos += 1

    return map(str,[baseCount['A'],baseCount['C'],baseCount['G'],baseCount['T']])



'''check right order of bam and gtf files'''
if len(args.bams) != len(args.files):
    raise ValueError("Number of GVF and Bam Files should be the same")
counter=0    
for bam in args.bams:
    bam=os.path.basename(bam)
    gtf=os.path.basename(args.files[counter])
    
    bamBase = bam.split(".")[0]
    gtfBase = gtf.split(".")[0]
     
    if bamBase != gtfBase:
        raise ValueError("BAM and GVF Files have to be in the same order")
    counter+=1


'''Fill the header'''
idDict = {}
if args.top==None:
    header=[] 
    for file in args.files:
        header.append(os.path.basename(file))   
else:
    header = args.top

'''fill table'''
idDict,keySet = fillDicts(args.files, args.columns,args.keys)

'''recount Reads'''
fileCounter=0
defaultList= ["--"]*len(args.columns)
for bamFile in args.bams:
    i=0
    #Helper.status("recounting Reads for %s" % bamFile)    
    Helper.info("recounting Reads from %s" % bamFile)
    samfile = pysam.AlignmentFile(bamFile, "rb")
    for keyTuple in keySet[1:]:
        i+=1
        '''check if basecount is unset for current condition''' 
        if idDict[keyTuple][fileCounter] == defaultList: 
            chr,startAnalysis = keyTuple[3],int(keyTuple[7])-1 #pysam is zero based        
            reads=samfile.fetch(chr, startAnalysis, startAnalysis+1)
            baseCount = getBaseCount(reads,startAnalysis)
            idDict[keyTuple][fileCounter] = baseCount
    if counter % 1000 == 0:
        Helper.status("%s out of %s editing sites finished" % (i,len(keySet)))
    
    fileCounter+=1
    
        
'''write the results to the output file'''
outFile = open(args.outFile,"w")     
deli="\t"*len(args.columns)
outFile.write("\t"*len(args.keys)+deli.join(header)+"\n")
for keyTuple in keySet:
    output=list(keyTuple)
    for v in idDict[keyTuple]:
        output=output+v
    outFile.write("\t".join(output)+"\n")
