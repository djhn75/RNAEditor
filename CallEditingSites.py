'''
Created on May 23, 2013

@author: david
'''

import argparse, multiprocessing, os, sys, re
from Helper import Helper
from genericpath import exists
from fileinput import close



class CallEditingSites(object):
    '''
    classdocs
    '''

    def printAttributes(self):
        print
        print "*** CALL VARIANTS WITH FOLLOWING ATTRIBUTES ***"
        print "\t Bam-File: " + self.bamFile
        print "\t outfilePrefix:" + self.outfilePrefix
        print "\t refGenome:" + self.refGenome
        print "\t dbsnp:" + self.dbsnp
        print "\t HapMap:" + self.hapmap
        print "\t 1000G Omni:" + self.omni
        print "\t Alu-Regions:" + self.aluRegions
        
        print "\t sourceDir:" + self.sourceDir
        print "\t threads:" + self.threads
        print "\t StandCall:" + self.standCall
        print "\t standEmit:" + self.standEmit
        print "\t keepTemp:" + str(self.keepTemp)
        print "\t overwrite:" + str(self.overwrite)
        print

    def __init__(self, bamFile, refGenome, dbsnp,
                 hapmap, omni, esp, 
                 aluRegions, geneAnnotationFile, outfilePrefix="default",
                 sourceDir="/usr/local/bin/", threads=multiprocessing.cpu_count()-1,standCall=0,
                 standEmit=0, edgeDistance=6, keepTemp=False, 
                 overwrite=False):
        '''
        Constructor
        set all the class Arguments
        '''
        self.debug=True
        
        self.bamFile=bamFile
        self.refGenome=refGenome
        self.dbsnp=dbsnp
        self.hapmap=hapmap
        self.omni=omni
        self.esp=esp
        self.aluRegions=aluRegions
        self.geneAnnotationFile = geneAnnotationFile
        if outfilePrefix=="default":
            self.outfilePrefix=self.bamFile[0:self.bamFile.rfind(".realigned")]
        else:
            self.outfilePrefix=outfilePrefix
        self.sourceDir=sourceDir
        self.threads=str(threads)
        self.standCall=str(standCall)
        self.standEmit=str(standEmit)
        self.edgeDistance=edgeDistance
        self.keepTemp=keepTemp
        self.overwrite=overwrite
        
        self.features = Helper.readGeneFetaures(self.geneAnnotationFile)
        
        
        self.logFile=open(self.outfilePrefix + ".log","a")
        if self.debug==True:
            self.printAttributes()
        
        self.checkDependencies()
    
    
    
    '''check if all the needed files are there'''
    def checkDependencies(self):
        #TODO: check for index Files
        if not os.path.exists(self.dbsnp):
            Exception("dbSNP File: Not found!!!")
        if not os.path.exists(self.refGenome):
            Exception("reference Genome File: Not found!!!")
        if not os.path.exists(self.hapmap):
            Exception("reference Genome File: Not found!!!")
        if not os.path.exists(self.refGenome):
            Exception("reference Genome File: Not found!!!")
        
    '''delete variants from Bam file which appear near read edges'''
    def removeEdgeMissmatches(self,vcfFile,bamFile,minDistance, minBaseQual, outFile):
        #pass
        #Loop through vcf-File
            #call overlapping reads
            #loop over reads
                #discard variants wich appear ONLY near edges
                #write the rest to the output file
        startTime=Helper.getTime()
        description = "remove Missmatches from the first " + str(minDistance) + "bp of the reads"
        print >> self.logFile, "[" + startTime.strftime("%c") + "] * * * " + description + " * * *"
        self.logFile.flush()
        print "[" + startTime.strftime("%c") + "] * * * " + description + " * * *"
        
        num_lines = str(sum(1 for line in open(vcfFile)))
        
        vcfFile=open(vcfFile,"r")
        counter=0
        if not exists(outFile):
            outFile=open(outFile,"w")
        else:
            print "\t [SKIP] File already exist"
            return
        for line in vcfFile:
            line=line.split("\t")
            snpPos=int(line[1])
            mmBase = line[4]
            position=line[0]+":"+line[1]+"-"+line[1]
            keepSNP=False
            
            #print position, str(minDistance)
            #samout= os.system("samtools view " + bamFile + " " + position)
            #command="samtools view " + bamFile + " " + position
            command = ["samtools", "view", bamFile, position]
            samout = Helper.getCommandOutput(command).splitlines()
            for samLine in samout:
                samfields=samLine.split()
                flag,startPos,mapQual,cigar,sequence,seqQual = samfields[1],int(samfields[3]),samfields[4],samfields[5],samfields[9],samfields[10]
                readPos=0
                mmReadPos=0
                cigarNums=re.split("[MIDNSHP]", cigar)[:-1]
                cigarLetters=re.split("[0-9]+",cigar)[1:]
                
                
                for i in range(len(cigarLetters)): #parse over single read
                    if cigarLetters[i] in {"I","S","H"}: #Insertion, Soft Clipping and Hard Clipping
                        readPos = readPos + int(cigarNums[i])
                    elif cigarLetters[i] in {"D","N"}: #Deletions and skipped Regions
                        startPos = startPos + int(cigarNums[i])
                    elif cigarLetters[i] in {"M"}: #Matches
                        for j in range(int(cigarNums[i])):
                            if startPos == snpPos:
                                mmReadPos = readPos 
                            readPos += 1
                            startPos += 1
                
                if mmReadPos != 0:
                                   
                    edgeDistance = int(snpPos) - int(startPos)
                
                    #only remove the snps from first 6 bases
                    revStrand = int(flag) & 16
                    if (revStrand == 0 and mmReadPos > minDistance) or (revStrand == 16 and mmReadPos < readPos - minDistance):
                        mmBaseQual= ord(seqQual[mmReadPos])
                        mmReadBase= sequence[mmReadPos]
                        if(mmBaseQual >= minBaseQual + 33) and (mmReadBase == mmBase): #check for quality of the base and the read contains the missmatch
                            keepSNP=True
                    #print "   ".join([str(revStrand),str(keepSNP),str(edgeDistance),str(len(sequence)),flag,startPos,mapQual,cigar,sequence,seqQual])
                #print distance
            
            
            if keepSNP:
                outFile.write("\t".join(line))  #print SNP
            counter+=1
            if counter % 1000 == 0:
                sys.stdout.write("\r" + str(counter) + " of " + num_lines + " missmatches finished")
                sys.stdout.flush()
        sys.stdout.write("\r" + num_lines + " of " + num_lines + " missmatches finished")
        sys.stdout.flush()
        duration=Helper.getTime()-startTime
        print >> self.logFile, "\t[DONE]" + " Duration [" + str(duration) + "]"
        self.logFile.flush()
        print "\t[DONE]" + " Duration [" + str(duration) + "]"    
            #sys.exit(0)
    

    
    
    ''' remove variant near splice junctions'''
    def removeIntronSpliceJunction(self,vcfFile,annotationFile,outFile):
        startTime=Helper.getTime()
        description = "remove Missmatches from intronic splice junctions"
        print >> self.logFile, "[" + startTime.strftime("%c") + "] * * * " + description + " * * *"
        self.logFile.flush()
        print "[" + startTime.strftime("%c") + "] * * * " + description + " * * *"
        
        annotationFile=open(annotationFile)
        if not exists(outFile):
            outFile=open(outFile,"w")
        else:
            print >> self.logFile, "\t [SKIP] File already exist"
            print "\t [SKIP] File already exist"
            return
        
        geneDict={} #save all genes according to the Chromosome in an Dictionary
        for line in annotationFile:
            line=line.split()
            chr=line[1]
            if chr in geneDict:
                geneDict[chr] = geneDict[chr] + [line]
            elif chr not in geneDict:
                geneDict[chr] = [line]
        annotationFile.close()
        
        vcfFile = open(vcfFile)
        for line in vcfFile:
            keepSnp=True
            line=line.split()
            mmChr,mmPos = line[0],int(line[1])
            if mmChr not in geneDict.keys():
                continue 
            for gene in geneDict[mmChr]:
                
                refChr,refStart,refStop = gene[1],int(gene[3]),int(gene[4])
                if mmPos > refStart and mmPos < refStop: #check if is inside of gene location 
                    numberExons = int(gene[7])
                    exonStarts = gene[8].split(",")
                    exonEnds = gene[9].split(",")
                    for i in range(numberExons): #check if variant lies in Exon
                        if int(exonStarts[i])-4 < mmPos and int(exonStarts[i])+1 > mmPos: #read is in front of exon
                            keepSnp=False
                        elif int(exonEnds[i]) < mmPos and  int(exonEnds[i])+4 > mmPos:  #read is at the end of exon
                            keepSnp=False
            if keepSnp == True:
                outFile.write("\t".join(line) + "\n")
        
        duration=Helper.getTime()-startTime
        print >> self.logFile, "\t[DONE]" + " Duration [" + str(duration) + "]"
        self.logFile.flush()
        print "\t[DONE]" + " Duration [" + str(duration) + "]"
    
    
    '''remove missmatches from homopolymers'''
    def removeHomopolymers(self,vcfFile,outFile,distance):
        startTime=Helper.getTime()
        description = "remove Missmatches from homopolymers"
        print >> self.logFile, "[" + startTime.strftime("%c") + "] * * * " + description + " * * *"
        self.logFile.flush()
        print "[" + startTime.strftime("%c") + "] * * * " + description + " * * *"
        
        tempBedFile = open(vcfFile+"_tmp.bed","w+")
        tempSeqFile = vcfFile + "_tmp.tsv"
        
        vcfFile=open(vcfFile,"r")
        mmDict = {}
        #print BedFile
        for line in vcfFile:
            lineSplit = line.split()
            chr,position,referenceNucleotide = lineSplit[0],lineSplit[1],lineSplit[3]
            siteNuc = chr+":"+position+"-" + referenceNucleotide
            mmDict[siteNuc] = line
            
            startPos = int(position) - distance if position >= distance else 0
            endpos = int(position) + distance
            
            tempBedFile.write("\t".join([chr,str(startPos),str(endpos),siteNuc])+"\n")
        
        tempBedFile.close()
        #run fastaFromBed
        cmd=[self.sourceDir+"bedtools/fastaFromBed", "-name", "-tab", "-fi", self.refGenome, "-bed", tempBedFile.name, "-fo", tempSeqFile]
        Helper.proceedCommand("catch surrounding sequences of Missmatches", cmd, tempBedFile.name, tempSeqFile, self.logFile, self.overwrite)
        
        mmNumberTotal = len(mmDict)
        
        #read sequence file
        tempSeqFile= open(tempSeqFile)
        for line in tempSeqFile:
            siteNuc,sequence = line.split()
            site,referenceNucleotide = siteNuc.split("-")
            #check if mm sorounding sequence are homopolymer nukleotides
            
            pattern = referenceNucleotide*distance
            
            """ !!!Test if this gives better results
                !!!ONLY DELETE IF MM IS AT THE END OF A HOMOPOLYMER NUKLEOTIDES
            if sequence.startswith(pattern):
                del mmDict[site] 
            elif sequence.endswith(pattern):
                del mmDict[site]
            """
            if pattern in sequence:
                del mmDict[siteNuc]
            
        #output the surviving Missmatches
        outFile = open(outFile,"w+")        
        for site in mmDict.keys():
            outFile.write(mmDict[site])
        
        #output statistics
        print >> self.logFile, "\t\t %d out of %d passed the Homopolymer-Filter" % (len(mmDict), mmNumberTotal)
        print "\t\t %d out of %d passed the Homopolymer-Filter" % (len(mmDict), mmNumberTotal)
        
        if self.keepTemp == False:
            os.remove(tempBedFile.name)
            os.remove(tempSeqFile.name)   
            
        
    '''do blat search (delete variants from reads that are not uniquely mapped)'''
    def blatSearch(self,vcfFile, outFile, minBaseQual, minMissmatch):
        startTime=Helper.getTime()
        description = "look for non uniquely mapped reads by blat"
        print >> self.logFile, "[" + startTime.strftime("%c") + "] * * * " + description + " * * *"
        self.logFile.flush()
        print "[" + startTime.strftime("%c") + "] * * * " + description + " * * *"
       
        counter=0
        geneHash = {}
        variantFile=open(vcfFile,"r")
        tempFasta = vcfFile + "_tmp.fa"
        if not os.path.isfile(tempFasta) or not os.path.getsize(tempFasta) > 0: #check if temFast exists and is not empty. If it exist it will not be created again
            tempFastaFile=open(tempFasta,"w+")
            num_lines = str(sum(1 for line in open(vcfFile)))
            
            #write missmatch read to fasta file
            for line in variantFile:
                line=line.split("\t")
                chromosome,snpPos,mmBase = line[0], int(line[1]), line[4]
                position=line[0]+":"+ str(snpPos)+"-"+str(snpPos)
                missmatchReadCount=1
    
                samout = Helper.getCommandOutput([self.sourceDir+"samtools", "view", "-F", "1024", self.bamFile, position]).splitlines() #-F 1024 to filter out duplicate reads
                for samLine in samout:
                    samfields=samLine.split()
                    flag,startPos,mapQual,cigar,sequence,seqQual = samfields[1],int(samfields[3]),samfields[4],samfields[5],samfields[9],samfields[10]
                    readPos=0
                    mmReadPos=0
                    keepRead=False
                    cigarNums=re.split("[MIDNSHP]", cigar)[:-1]
                    cigarLetters=re.split("[0-9]+",cigar)[1:]
                    
                    for i in range(len(cigarLetters)): #parse over single read
                        if cigarLetters[i] in {"I","S","H"}: #Insertion, Soft Clipping and Hard Clipping
                            readPos = readPos + int(cigarNums[i])
                        elif cigarLetters[i] in {"D","N"}: #Deletions and skipped Regions
                            startPos = startPos + int(cigarNums[i])
                        elif cigarLetters[i] in {"M"}: #Matches
                            for j in range(int(cigarNums[i])):
                                if startPos == snpPos:
                                    mmReadPos = readPos
                                    mmBaseQual= ord(seqQual[mmReadPos])
                                    mmReadBase= sequence[mmReadPos]
                                    if(mmBaseQual >= minBaseQual + 33) and (mmReadBase == mmBase): #check for quality of the base and the read contains the missmatch
                                        keepRead=True
                                readPos += 1
                                startPos += 1
                    if keepRead == True: #if read contains the missmatch 
                        tempFastaFile.write("> " + chromosome + "-" + str(snpPos) + "-" + str(missmatchReadCount) + "\n" + sequence + "\n")
                        missmatchReadCount += 1
    
                counter += 1
                if counter % 1000 == 0:
                    sys.stdout.write("\r" + str(counter) + " of " + num_lines + " missmatched reads written")
                    sys.stdout.flush()
        
            
            print "\n created fasta file " + tempFasta
            tempFastaFile.close()
        variantFile.close()    
        
        
        #do blat search
        pslFile=outFile+".psl"
        if not os.path.isfile(pslFile) or not os.path.getsize(pslFile) > 0:
            cmd = [self.sourceDir+"blat","-stepSize=5","-repMatch=2253", "-minScore=20","-minIdentity=0","-noHead", self.refGenome, tempFasta, pslFile]
            print cmd
            Helper.proceedCommand("do blat search for unique reads",cmd,tempFasta, "None", self.logFile, self.overwrite)
            
        
        startTime=Helper.getTime()
        description = "check for variants that pass the blat criteria"
        print >> self.logFile, "[" + startTime.strftime("%c") + "] * * * " + description + " * * *"
        self.logFile.flush()
        print "[" + startTime.strftime("%c") + "] * * * " + description + " * * *"
        
        if not os.path.isfile(outFile):
            #open psl file
            pslFile=open(pslFile)
            blatDict={}
            for line in pslFile: #summarize the blat hits
                pslFields = line.split()
                name = pslFields[9]
                blatScore = [pslFields[0], pslFields[13], pslFields[17], pslFields[18], pslFields[20]] # #of Matches, targetName, blockCount, blockSize, targetStarts 
                if name in blatDict:
                    blatDict[name] = blatDict[name] + [blatScore]
                else:
                    blatDict[name] = [blatScore]
    
    
            siteDict = {}
            discardDict = {}
            
        
            #loop over blat Hits
            for pslKey in blatDict.keys():      #Loop over all blat hits of mmReads to observe the number of Alignements   
                keepSNP=False
                chr,pos=pslKey.split("-")[0:2]
                
                site = ":".join([chr,pos])
                pslLine = blatDict[pslKey]
                largestScore=0
                largestScoreLine=pslLine[0]
                scoreArray=[]
                
                #look for largest blatScore and save the largest line too
                for blatHit in pslLine: 
                    lineScore=int(blatHit[0])
                    scoreArray.append(lineScore)
                    if lineScore > largestScore:
                        largestScore = lineScore
                        largestScoreLine=blatHit
                pos=int(pos)
                scoreArray.sort(reverse=True)
                if len(scoreArray) < 2:   #test if more than one blat Hit exists
                    scoreArray.append(0)
                if chr == largestScoreLine[1] and scoreArray[1] < scoreArray[0]*0.95: #check if same chromosome and hit is lower the 95 perchen of first hit
                    blockCount,blockSizes,blockStarts = int(largestScoreLine[2]),largestScoreLine[3].split(",")[:-1],largestScoreLine[4].split(",")[:-1]
                    for i in range(blockCount):
                        startPos = int(blockStarts[i])+1
                        endPos = startPos + int(blockSizes[i])
                        if pos >= startPos and pos < endPos: #check if alignement overlaps missmatch
                            keepSNP = True
                
                if keepSNP == True:
                    if site in siteDict:
                        siteDict[site]+=1
                    else:
                        siteDict[site]=1
                elif keepSNP == False: #when read not passes the blat criteria
                    if site in discardDict:
                        discardDict[site]+=1
                    else:
                        discardDict[site]=1
            pslFile.close() 
            
                    
            #loop through variant file again and check what passes the blat criteria
            variantFile=open(vcfFile)
            resultFile=open(outFile,"w+")
            
            mmNumberTotal=0
            mmNumberTooSmall=0
            mmReadsSmallerDiscardReads=0 
            for line in variantFile:
                line = line.split()
                name=":".join(line[0:2])
                numberBlatReads=0
                numberDiscardReads=0
                if name in siteDict:
                    numberBlatReads = siteDict[name]
                if name in discardDict:
                    numberDiscardReads = discardDict[name]
                
                if numberBlatReads >= minMissmatch and numberBlatReads >= numberDiscardReads:
                    resultFile.write("\t".join(line)+"\n")
                
                
                #count statistics
                if numberBlatReads < minMissmatch:
                    mmNumberTooSmall+=1
                elif numberBlatReads < numberDiscardReads: #check if more readd fit the blat criteria than not
                    mmReadsSmallerDiscardReads+=1    
                mmNumberTotal+=1
            variantFile.close()
            
            if self.keepTemp == False:
                os.remove(tempFasta)
                os.remove(pslFile.name)
            
            #output statisticsttkkg
            
            duration=Helper.getTime()-startTime
            mmPassedNumber=mmNumberTotal-(mmNumberTooSmall+mmReadsSmallerDiscardReads)
            print >> self.logFile, "\t\t %d out of %d passed blat criteria" % (mmPassedNumber, mmNumberTotal)
            print >> self.logFile, "\t\t %d Missmatches had fewer than %d missmatching-Reads." % (mmNumberTooSmall, minMissmatch)
            print >> self.logFile, "\t\t %d Missmatches had more missaligned reads than correct ones." % (mmReadsSmallerDiscardReads)
            self.logFile.flush()
            print "\t\t %d out of %d passed blat criteria (%d percent)" % (mmPassedNumber, mmNumberTotal,(mmPassedNumber/mmNumberTotal)*100)
            print "\t\t %d Missmatches had fewer than %d missmatching-Reads." % (mmNumberTooSmall, minMissmatch)
            print "\t\t %d Missmatches had more missaligned reads than correct ones." % (mmReadsSmallerDiscardReads)
            
        duration=Helper.getTime()-startTime
        print >> self.logFile, "\t[DONE]" + " Duration [" + str(duration) + "]"
        self.logFile.flush()
        print "\t[DONE]" + " Duration [" + str(duration) + "]"

    #write two variantFiles together
    
    '''combine two variant File'''
    def combineVariants(self,aluSites,nonAuluSites,outFile):
        aluSites=open(aluSites,"r")
        nonAuluSites = open(nonAuluSites,"r")
        outFile = open(outFile,"w+")
        for line in aluSites:
            outFile.write(line)
        for line in nonAuluSites:
            outFile.write(line)
        
        
    def annotateVariants(self,annotationFile,variantFile,outFile):
        startTime=Helper.getTime()
        description = "Annotate Missmatches"
        print >> self.logFile, "[" + startTime.strftime("%c") + "] * * * " + description + " * * *"
        self.logFile.flush()
        print "[" + startTime.strftime("%c") + "] * * * " + description + " * * *"
        
        #check if file exitsts
        if not exists(outFile):
            outFile=open(outFile,"w")
        else:
            print >> self.logFile, "\t [SKIP] File already exist"
            print "\t [SKIP] File already exist"
            return
        
        #write header
        outFile.write("\t".join(["#CHROM","POS","ID","REF","ALT","QUAL","GENE","INFO"])+"\n")
        
        """annotationFile=open(annotationFile,"r")
        geneDict={} #save all genes according to the Chromosome in an Dictionary
        for line in annotationFile:
            line=line.split()
            chr=line[1]
            if chr in geneDict:
                geneDict[chr] = geneDict[chr] + [line]
            elif chr not in geneDict:
                geneDict[chr] = [line]
        annotationFile.close()
"""        
        
        
        
        geneCountDict = {}
        vcfFile = open(variantFile) 
        for line in vcfFile:
            keepSnp=True
            line=line.split()
            mmChr, mmPos, mmRef, mmAlt = line[0],int(line[1]), line[3], line[4]
            
            
            if mmChr not in geneDict.keys(): #if chromosome not in GeneDict
                continue 
            
            if not ((mmRef=="T" and  mmAlt=="C") or (mmRef=="A" and mmAlt=="G")): #only proceed with A->G or T->C events
                continue
            
            location = "intergenic"
            geneName = "-"
            for gene in geneDict[mmChr]: #loop through all chromosome Genes
                
                geneId,refChr,strand,refStart,refStop = gene[0],gene[1],gene[2],int(gene[3]),int(gene[4])
                if mmPos > refStart and mmPos < refStop: #check if is inside of gene location 
                    numberExons = int(gene[7])
                    exonStarts = gene[8].split(",")
                    exonEnds = gene[9].split(",")
                    geneName=geneId
                    if geneId  in geneCountDict:
                        countList = geneCountDict[geneId]
                    else:
                        #countList = [chr,refStart,refStop,5'utr,3'utr,intronic,exonic,totalNumber]
                        countList = [refChr,refStart,refStop,0,0,0,0,0]
                    if mmPos < int(exonStarts[0]):
                        if strand == "+":
                            location = "5'UTR"
                            countList[3]+=1
                        elif strand == "-":
                            location = "3'UTR"
                            countList[4]+=1
                    elif mmPos > int(exonEnds[-2]):
                        if strand == "+":
                            location = "3'UTR"
                            countList[4]+=1
                        elif strand == "-":
                            location = "5'UTR"
                            countList[3]+=1
                    else:
                        for i in range(numberExons): #check if variant lies in Exon
                            if int(exonStarts[i]) < mmPos and int(exonEnds[i]) > mmPos: #read is in front of exon
                                location = "exonic"
                                countList[5]+=1
                        if location != "exonic":
                            location = "intronic"
                            countList[6]+=1
                    
                    countList[7]+=1 #count total number of missmatches in 
                    geneCountDict[geneId]=countList
            
            #write ouputFile 
            outFile.write("\t".join([line[0],line[1],line[2],line[3],line[4],line[5],geneName,location]) + "\n")
            
        #print count Table        
        geneCountsFile = open(self.outfilePrefix + ".editedGenes.counts","w+")
        geneCountsFile.write("\t".join(["Gene", "start", "stop","#5'utr","#3'utr","#intronic","#exonic","#total"])+"\n")
        for geneName in geneCountDict.keys():
            temp=map(str,geneCountDict[geneName]) #convert list to list of str
            geneCountsFile.write(geneName + "\t" + "\t".join(temp)+"\n")               
            
            
            
            
            
        
        duration=Helper.getTime()-startTime
        print >> self.logFile, "\t[DONE]" + " Duration [" + str(duration) + "]"
        self.logFile.flush()
        print "\t[DONE]" + " Duration [" + str(duration) + "]"
            
    def __del__(self):
        if self.keepTemp==False:
            os.remove(self.outfilePrefix+".vcf")
            os.remove(self.outfilePrefix+".no_dbsnp.vcf")
           
            os.remove(self.outfilePrefix+".no_dbsnp.no_1000genome.vcf")
            os.remove(self.outfilePrefix+".no_dbsnp.no_1000genome.no_esp.vcf")
            os.remove(self.outfilePrefix+".no_dbsnp.no_1000genome.no_esp.noStartMM.vcf")
            os.remove(self.outfilePrefix+".nonAlu.vcf")
            os.remove(self.outfilePrefix+".nonAlu.noSpliceSites.vcf")
            os.remove(self.outfilePrefix+".nonAlu.noSpliceSites.noHomo.vcf")
            
    
    def start(self):
        #Rough variant calling with GATK
        vcfFile=self.outfilePrefix+".vcf"
        cmd = ["java","-Xmx6G","-jar",self.sourceDir + "GATK/GenomeAnalysisTK.jar", 
               "-T","UnifiedGenotyper","-R", self.refGenome, "-glm", "SNP","-I", self.bamFile, 
               "-D", self.dbsnp, "-o", vcfFile, "-metrics", self.outfilePrefix+".snp.metrics", "-nt", self.threads, "-l","ERROR",
               "-stand_call_conf", self.standCall, "-stand_emit_conf", self.standEmit,"-A", "Coverage", "-A", "AlleleBalance","-A", "BaseCounts"]
        #print cmd
        Helper.proceedCommand("Call variants", cmd, self.bamFile, vcfFile, self.logFile, self.overwrite)
    
        #delete SNPs from dbSNP
        noDbsnp=self.outfilePrefix+".no_dbsnp.vcf"
        cmd = [self.sourceDir+"bedtools/intersectBed","-v","-a",vcfFile,"-b",self.dbsnp]
        #print cmd
        Helper.proceedCommand("delete SNPs from dbSNP",cmd, vcfFile, noDbsnp, self.logFile, self.overwrite)
        
        #delete variants from 1000 Genome Project
        no1000G = self.outfilePrefix + ".no_dbsnp.no_1000genome.vcf"
        cmd =  [self.sourceDir+"bedtools/intersectBed","-v","-a",noDbsnp,"-b",self.omni]
        Helper.proceedCommand("delete SNPs from 1000 Genome Omni database",cmd, noDbsnp, no1000G, self.logFile, self.overwrite)
        
        #delete variants from UW exome calls
        noEsp = self.outfilePrefix + ".no_dbsnp.no_1000genome.no_esp.vcf"
        cmd =  [self.sourceDir+"bedtools/intersectBed","-v","-a",no1000G,"-b",self.esp]
        Helper.proceedCommand("delete SNPs from Exome Sequencing Project",cmd, noDbsnp, noEsp, self.logFile, self.overwrite)
        
        #erase artificial missmatches from read-starts
        noStartMissmatches= self.outfilePrefix + ".no_dbsnp.no_1000genome.no_esp.noStartMM.vcf"
        self.removeEdgeMissmatches(noEsp, self.bamFile, self.edgeDistance, 25, noStartMissmatches)
        
        #split non-Alu and Alu regions
        nonAlu = self.outfilePrefix + ".nonAlu.vcf"
        alu = self.outfilePrefix + ".alu.vcf"
        cmd =  [self.sourceDir+"bedtools/intersectBed","-v","-a",noStartMissmatches,"-b",self.aluRegions]
        Helper.proceedCommand("write variants from non-alu regions",cmd, noStartMissmatches, nonAlu, self.logFile, self.overwrite)  # write nonAlu-Regions
        cmd =  [self.sourceDir+"bedtools/intersectBed","-a",noStartMissmatches,"-b",self.aluRegions]
        Helper.proceedCommand("write variants from alu regions",cmd, noStartMissmatches, alu, self.logFile, self.overwrite) #write alu-regions
        
        #proceed with non-Alu reads only!!!
        #erase variants from intronic splice junctions
        noSpliceSites = self.outfilePrefix + ".nonAlu.noSpliceSites.vcf"
        self.removeIntronSpliceJunction(nonAlu, self.geneAnnotationFile, noSpliceSites)
        
        #remove sites from simple repeats
        
        #cmd =  [self.sourceDir+"bedtools/intersectBed","-v","-a",noStartMissmatches,"-b",self.simpleRepeats]
        
        #erase variants from homopolymer runs
        noHomo = self.outfilePrefix + ".nonAlu.noSpliceSites.noHomo.vcf"
        self.removeHomopolymers(noSpliceSites, noHomo, 4)
        
        #do blat search
        blatOutfile = self.outfilePrefix + ".nonAlu.noSpliceSites.noHomo.blat.vcf"
        self.blatSearch(noHomo, blatOutfile, 25, 1)
        
        #combine alu and non Alu sites
        combinedFile = self.outfilePrefix + ".editingSites.vcf"
        self.combineVariants(alu, blatOutfile, combinedFile)
        
        annotatedFile = self.outfilePrefix + ".annotated.vcf"
        self.annotateVariants(self.geneAnnotationFile, combinedFile, annotatedFile)
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='output vatiants from a given .bam file.')
    parser.add_argument('-i', '--input', metavar='bam-File', type=argparse.FileType('r'), help='Input bam file from which variants should be called', required=True)
    parser.add_argument("-r", "--RefGenome", metavar='Fasta-File', help="File that contains the reference sequences", type=argparse.FileType('r'), default='/media/media/databases/human/human_g1k_v37.fa')
    parser.add_argument('-s', '--dbsnp', help='SNP database (dbSNP) in VCF format (downloaded from the GATK homepage)', type=argparse.FileType('r'), default='/media/media/databases/human/dbsnp_135.b37.vcf')
    parser.add_argument('-m', '--hapmap', help='hapmap database in vcf format (see GATK homepage)', type=argparse.FileType('r'), default='/media/media/databases/human/hapmap_3.3.b37.sites.vcf')
    parser.add_argument('-g', '--omni', help='1000 Genome variants in vcf format (see GATK homepage)', type=argparse.FileType('r'), default='/media/media/databases/human/1000G_omni2.5.b37.sites.vcf')
    parser.add_argument('-e', '--esp', help='Exome Sequencing Project variants', type=argparse.FileType('r'), default='/media/media/media/databases/human/NHLBI_Exome_Sequencing_Project_6500SI.vcf')
    parser.add_argument('-a', '--aluRegions', help='Alu-Regions downloaded fron the UCSC table browser', type=argparse.FileType('r'), default='/media/media/media/databases/human/hg19/rna-editing/Alu_repeats_noChr.bed')
    parser.add_argument('-G', '--geneAnnotation', help='Gene annotation File in bed format', type=argparse.FileType('r'), default='/media/media/databases/human/hg19/UCSC_Genes_noChr.ucsc')
    parser.add_argument('-o', '--output', metavar='output-prefix', type=str,help='prefix that is written in front of the output files', default="default")
    parser.add_argument('-d', '--sourceDir', help='- Directory to all the tools [default: bin/]', default='bin/', type=Helper.readable_dir)
    parser.add_argument('-t', '--threads', help='number of threads', type=int, default=multiprocessing.cpu_count()-1)
    parser.add_argument('-sc', '--standCall', help='-The minimum phred-scaled confidence threshold at which variants should be considered as true (int) [0]', type=int, default=0)
    parser.add_argument('-se', '--standEmit', help=' The minimum phred-scaled confidence threshold at which variants should be emitted (int)[0]', type=int, default=0)
    parser.add_argument('-ed', '--edgeDistance', help='The minimum edge distance of the SNPs', type=int, default=6)
    parser.add_argument('--keepTemp', help='keep the intermediate Files [False]', action='store_true', default=False)
    parser.add_argument('--overwrite', help='overwrite existing Files [False]', action='store_true', default=False)
    
    args = parser.parse_args()

    call=CallEditingSites(args.input.name, args.RefGenome.name, args.dbsnp.name, 
                          args.hapmap.name, args.omni.name, args.esp, 
                          args.aluRegions, args.geneAnnotation, args.output, 
                          args.sourceDir, args.threads, args.standCall, 
                          args.standEmit, args.edgeDistance, args.keepTemp, 
                          args.overwrite)