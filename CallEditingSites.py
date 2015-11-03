'''
Created on May 23, 2013

@author: david
'''

import os, sys, re
from Helper import Helper
from VariantSet import VariantSet
from Genome import Genome
from copy import copy
import gc



class CallEditingSites(object):
    '''
    classdocs
    '''
    def __init__(self, bamFile, rnaEdit):
        '''
        Constructor
        set all the class Arguments
        '''
        self.debug=True
        
        self.bamFile=bamFile
        self.rnaEdit=rnaEdit
        


    def printAttributes(self):
        
        Helper.info("*** CALL VARIANTS WITH FOLLOWING ATTRIBUTES ***", self.rnaEdit.logFile, self.rnaEdit.textField) 
        Helper.info( "\t Bam-File: " + self.bamFile, self.rnaEdit.logFile, self.rnaEdit.textField) 
        Helper.info( "\t outfilePrefix:" + self.rnaEdit.params.output, self.rnaEdit.logFile, self.rnaEdit.textField) 
        Helper.info( "\t refGenome:" + self.rnaEdit.params.refGenome, self.rnaEdit.logFile, self.rnaEdit.textField) 
        Helper.info( "\t dbsnp:" + self.rnaEdit.params.dbsnp, self.rnaEdit.logFile, self.rnaEdit.textField) 
        Helper.info( "\t HapMap:" + self.rnaEdit.params.hapmap, self.rnaEdit.logFile, self.rnaEdit.textField) 
        Helper.info( "\t 1000G Omni:" + self.rnaEdit.params.omni, self.rnaEdit.logFile, self.rnaEdit.textField) 
        Helper.info( "\t Alu-Regions:" + self.rnaEdit.params.aluRegions, self.rnaEdit.logFile,self.rnaEdit.textField) 
        
        Helper.info( "\t sourceDir:" + self.rnaEdit.params.sourceDir, self.rnaEdit.logFile, self.rnaEdit.textField) 
        Helper.info( "\t threads:" + self.rnaEdit.params.threads, self.rnaEdit.logFile, self.rnaEdit.textField) 
        Helper.info( "\t StandCall:" + self.rnaEdit.params.standCall, self.rnaEdit.logFile, self.rnaEdit.textField) 
        Helper.info( "\t standEmit:" + self.rnaEdit.params.standEmit, self.rnaEdit.logFile, self.rnaEdit.textField) 
        Helper.info( "\t keepTemp:" + str(self.rnaEdit.params.keepTemp), self.rnaEdit.logFile, self.rnaEdit.textField) 
        Helper.info( "\t overwrite:" + str(self.rnaEdit.params.overwrite), self.rnaEdit.logFile, self.rnaEdit.textField) 
        


          
    
    def removeEdgeMissmatches(self,variants,bamFile,minDistance, minBaseQual):
        '''delete variants from Bam file which appear near read edges'''
        #Loop through vcf-File
            #call overlapping reads with samtools view
            #loop over reads
                #discard variants wich appear ONLY near edges
                #write the rest to the output file
        startTime=Helper.getTime()
        minDistance=int(minDistance)
        
        counter=0    
        
        num_lines = len(variants.variantDict)
        Helper.info(" [%s] remove Missmatches from the first %s bp from read edges" % (startTime.strftime("%c"),str(minDistance)),self.rnaEdit.logFile,self.rnaEdit.textField)
        
        for varKey in variants.variantDict.keys():
            variant = variants.variantDict[varKey]
            snpPos = variant.position
            position=variant.chromosome+":" + str(snpPos) + "-" + str(snpPos) 
            #line[1]+"-"+line[1]
            keepSNP=False
            
            command = [self.rnaEdit.params.sourceDir+"samtools", "view", bamFile, position] 
            samout = Helper.getCommandOutput(command).splitlines() #get the reads wich are overlapping the snp region
            for samLine in samout: #loop over reads
                samfields=samLine.split()
                try:
                    flag,startPos,mapQual,cigar,sequence,seqQual = samfields[1],int(samfields[3]),samfields[4],samfields[5],samfields[9],samfields[10]
                except ValueError:
                    raise ValueError("Error in line '%s'" % " ".join(samfields))

                readPos=0
                mmReadPos=0
                cigarNums=re.split("[MIDNSHP]", cigar)[:-1]
                cigarLetters=re.split("[0-9]+",cigar)[1:]
                
                #loop over read cigar (check insertions,deletions and skipped regions) 
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
                if mmReadPos != 0: #happens when snp is in non matching regions (like Insertion, Soft Clipping and Hard Clipping, Deletions ans skipped regions)    
                    edgeDistance = snpPos - startPos
                
                    #only remove the snps from first minDistance bases
                    revStrand = int(flag) & 16 #check if 5th bit is set (results: 0 = +strand; 16= -strand)
                    if (revStrand == 0 and mmReadPos > minDistance) or (revStrand == 16 and mmReadPos < readPos - minDistance):
                        mmBaseQual= ord(seqQual[mmReadPos])
                        mmReadBase= sequence[mmReadPos]
                        if(mmBaseQual >= minBaseQual + 33) and (mmReadBase == variant.alt): #check for quality of the base and the read contains the missmatch
                            keepSNP=True

            if not keepSNP:
                del variants.variantDict[varKey]    
            counter+=1
            if counter % 10000 == 0: #print out current status
                Helper.status(str(counter) + " of " + str(num_lines) + " missmatches finished",self.rnaEdit.logFile,self.rnaEdit.textField)
          
    def removeIntronicSpliceJunctions(self,variants,genome,distance=4): 
        '''
        remove variant near splice junctions and returns the other variants
        :param variants: VariantSet
        :param genome: object of the class Genome
        '''
        startTime=Helper.getTime()
        
        Helper.info(" [%s] remove Missmatches from the intronic splice junctions " % (startTime.strftime("%c")),self.rnaEdit.logFile,self.rnaEdit.textField)
        #TODO Finish this fucking fuction
        
        geneDict = genome.getGenesByChromosome()

        for key in variants.variantDict.keys():
            delVar=False
            chromosome,position,ref,alt = key
            for gene in geneDict[chromosome]:
                if gene.start < position < gene.end:#check if is inside of gene location
                    for exon in gene.codingExons:
                        if (exon[0]-distance < position < exon[0]) or (exon[1] < position < exon[1]+distance):
                            #print(key)
                            delVar=True
            if delVar:
                del variants.variantDict[key]
                            
        Helper.printTimeDiff(startTime,self.rnaEdit.logFile,self.rnaEdit.textField)
        
    '''remove missmatches from homopolymers'''
    def removeHomopolymers(self,variants,outFile,distance):
        startTime=Helper.getTime()
        Helper.info(" [%s] remove Missmatches from homopolymers " % (startTime.strftime("%c")),self.rnaEdit.logFile,self.rnaEdit.textField)
        
        tempBedFile = open(outFile+"_tmp.bed","w+")
        tempSeqFile = outFile + "_tmp.tsv"
        
        #print temporary BedFile
        for key in variants.variantDict.keys():
            chr,position,ref,alt = key
            siteNuc = ",".join([chr,str(position),ref,alt])
            startPos = position - distance if position >= distance else 0
            endpos = position + distance
            
            tempBedFile.write("\t".join([chr,str(startPos),str(endpos),siteNuc])+"\n")
        
        tempBedFile.close()
        #run fastaFromBed
        cmd=[self.rnaEdit.params.sourceDir+"bedtools/fastaFromBed", "-name", "-tab", "-fi", self.rnaEdit.params.refGenome, "-bed", tempBedFile.name, "-fo", tempSeqFile]
        Helper.proceedCommand("catch surrounding sequences of Missmatches", cmd, tempBedFile.name, tempSeqFile,self.rnaEdit)
        
        mmNumberTotal = len(variants.variantDict)
        
        #read sequence file
        tempSeqFile= open(tempSeqFile)
        for line in tempSeqFile:
            siteNuc,sequence = line.split()
            try:
                chr,position,ref,alt = siteNuc.split(",",3)
            except (ValueError):
                raise ValueError("Failed to read line: %s" % line)
            #check if mm sorounding sequence are homopolymer nukleotides
            
            pattern = ref*distance
            
            """ !!!Test if this gives better results
                !!!ONLY DELETE IF MM IS AT THE END OF A HOMOPOLYMER NUKLEOTIDES
            if sequence.startswith(pattern):
                del mmDict[site] 
            elif sequence.endswith(pattern):
                del mmDict[site]
            """
            if pattern in sequence:
                try:
                    del variants.variantDict[(chr,int(position),ref,alt)]
                except KeyError:
                    pass
                
        #output statistics
        Helper.info("\t\t %d out of %d passed the Homopolymer-Filter" % (mmNumberTotal, mmNumberTotal),self.rnaEdit.logFile,self.rnaEdit.textField)
        Helper.printTimeDiff(startTime,self.rnaEdit.logFile,self.rnaEdit.textField)
        
        tempSeqFile.close()
        
        if self.rnaEdit.params.keepTemp == False:
            os.remove(tempBedFile.name)
            os.remove(tempSeqFile.name)   
                
    '''do blat search (delete variants from reads that are not uniquely mapped)'''
    def blatSearch(self,variants, outFile, minBaseQual, minMissmatch):
        startTime=Helper.getTime()
        Helper.info(" [%s] Search non uniquely mapped reads" % (startTime.strftime("%c")),self.rnaEdit.logFile,self.rnaEdit.textField)
        
        counter=0
        geneHash = {}
        tempFasta = outFile + "_tmp.fa"
        if not os.path.isfile(tempFasta) or not os.path.getsize(tempFasta) > 0: #check if temFast exists and is not empty. If it exist it will not be created again
            tempFastaFile=open(tempFasta,"w+")
            mmNumberTotal = len(variants.variantDict)
            
            #write missmatch read to fasta file
            for key in variants.variantDict.keys():
                chromosome,position,ref,alt = key
                samPos=chromosome+":"+ str(position)+"-"+str(position)
                missmatchReadCount=1
    
                samout = Helper.getCommandOutput([self.rnaEdit.params.sourceDir+"samtools", "view", "-F", "1024", self.bamFile, samPos]).splitlines() #-F 1024 to filter out duplicate reads
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
                                if startPos == position:
                                    mmReadPos = readPos
                                    mmBaseQual= ord(seqQual[mmReadPos])
                                    mmReadBase= sequence[mmReadPos]
                                    if(mmBaseQual >= minBaseQual + 33) and (mmReadBase == alt): #check for quality of the base and the read contains the missmatch
                                        keepRead=True
                                readPos += 1
                                startPos += 1
                    if keepRead == True: #if read contains the missmatch 
                        tempFastaFile.write("> "+chromosome+"-"+str(position)+"-"+ref+"-"+alt+"-"+str(missmatchReadCount)+"\n"+sequence+"\n")
                        missmatchReadCount += 1
    
                counter += 1
                if counter % 1000 == 0:
                    sys.stdout.write("\r" + str(counter) + " of " + str(mmNumberTotal) + " variants done")
                    sys.stdout.flush()
        
            Helper.info("\n created fasta file " + tempFasta,self.rnaEdit.logFile,self.rnaEdit.textField)
            Helper.printTimeDiff(startTime,self.rnaEdit.logFile,self.rnaEdit.textField)
            tempFastaFile.close()
                
        
        #do blat search
        pslFile=outFile+".psl"
        if not os.path.isfile(pslFile) or not os.path.getsize(pslFile) > 0:
            cmd = [self.rnaEdit.params.sourceDir+"blat","-stepSize=5","-repMatch=2253", "-minScore=20","-minIdentity=0","-noHead", self.rnaEdit.params.refGenome, tempFasta, pslFile]
            #print cmd
            Helper.proceedCommand("do blat search for unique reads",cmd,tempFasta, "None", self.rnaEdit)
        
        Helper.info(" [%s] look for non uniquely mapped reads by blat" % (startTime.strftime("%c")),self.rnaEdit.logFile,self.rnaEdit.textField)    
        
        if not os.path.isfile(outFile):
            #open psl file
            pslFile=open(pslFile)
            blatDict={}
            for line in pslFile: #summarize the blat hits
                pslFields = line.split()
                chr,pos,ref,alt,mmReadCount = pslFields[9].split("-")
                varTuple=(chr,int(pos),ref,alt)
                blatScore = [pslFields[0], pslFields[13], pslFields[17], pslFields[18], pslFields[20]] # #of Matches, targetName, blockCount, blockSize, targetStarts 
                if varTuple in blatDict:
                    blatDict[varTuple] = blatDict[varTuple] + [blatScore]
                else:
                    blatDict[varTuple] = [blatScore]

            siteDict = {}
            discardDict = {}
            
            #loop over blat Hits
            for varTuple in blatDict.keys():      #Loop over all blat hits of mmReads to observe the number of Alignements   
                keepSNP=False
                chr,pos,ref,alt=varTuple             
                pslLine = blatDict[varTuple]
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
                    if varTuple in siteDict:
                        siteDict[varTuple]+=1
                    else:
                        siteDict[varTuple]=1
                elif keepSNP == False: #when read not passes the blat criteria
                    if varTuple in discardDict:
                        discardDict[varTuple]+=1
                    else:
                        discardDict[varTuple]=1
            pslFile.close() 
            
                    
            #loop through variants again and check what passes the blat criteria
            
            
            mmNumberTotal=0
            mmNumberTooSmall=0
            mmReadsSmallerDiscardReads=0 
            for key in variants.variantDict.keys():
                numberBlatReads=0
                numberDiscardReads=0
                if key in siteDict:
                    numberBlatReads = siteDict[key]
                if key in discardDict:
                    numberDiscardReads = discardDict[key]
                
                if  numberBlatReads <= minMissmatch and numberBlatReads <= numberDiscardReads:
                    del variants.variantDict[key]
                
                
                #count statistics
                if numberBlatReads < minMissmatch:
                    mmNumberTooSmall+=1
                elif numberBlatReads < numberDiscardReads: #check if more readd fit the blat criteria than not
                    mmReadsSmallerDiscardReads+=1    
                mmNumberTotal+=1
            
            if self.rnaEdit.params.keepTemp == False:
                os.remove(tempFasta)
                os.remove(pslFile.name)
            
            #output statisticsttkkg
            mmPassedNumber=mmNumberTotal-(mmNumberTooSmall+mmReadsSmallerDiscardReads)
            
            Helper.info("\t\t %d out of %d passed blat criteria" % (mmPassedNumber, mmNumberTotal),self.rnaEdit.logFile,self.rnaEdit.textField)
            Helper.info("\t\t %d Missmatches had fewer than %d missmatching-Reads." % (mmNumberTooSmall, minMissmatch),self.rnaEdit.logFile,self.rnaEdit.textField)
            Helper.info("\t\t %d Missmatches had more missaligned reads than correct ones." % (mmReadsSmallerDiscardReads),self.rnaEdit.logFile,self.rnaEdit.textField)
            
        Helper.printTimeDiff(startTime,self.rnaEdit.logFile,self.rnaEdit.textField)

    
   
    
            
    def cleanUp(self):
        #print [x for x in gc.get_objects()]
        #print str(self) + " cleaned up"
        self.genome = None
        
        if self.rnaEdit.params.keepTemp==False:
            #os.remove(self.rnaEdit.params.output+".vcf")
            os.remove(self.rnaEdit.params.output+"_tmp.bep")
            os.remove(self.rnaEdit.params.output+"_tmp.tsv")
            os.remove(self.rnaEdit.params.output+"_blat.psl")
            os.remove(self.rnaEdit.params.output+"_blat_tmp.fa")
            #os.remove(self.rnaEdit.params.output+".nonAlu.vcf")
            #os.remove(self.rnaEdit.params.output+".nonAlu.noSpliceSites.vcf")
            #os.remove(self.rnaEdit.params.output+".nonAlu.noSpliceSites.noHomo.vcf")
        
            
    def deleteNonEditingBases(self,variants):
        startTime=Helper.getTime()
        Helper.info("Delete non Editing Bases (keep only T->C and A->G)",self.rnaEdit.logFile,self.rnaEdit.textField)
        
        for varTuple in variants.variantDict.keys():
            chr,pos,ref,alt = varTuple
            if (ref =="A" and alt == "G") or (ref=="T" and alt=="C"):
                pass
            else:
                del variants.variantDict[varTuple]
    
    def startAnalysis(self):
        '''Proceeds all the steps to detect editing Sites from a bam File
        
        @return: 0 on success and 1 if analysis was canceled by user
        '''
        
        '''check if result file already exists''' 
        if os.path.isfile(self.rnaEdit.params.output+".editingSites.clusters") and self.rnaEdit.params.overwrite==False:
            print "\t [SKIP] Final result file already exist",self.rnaEdit.logFile,self.rnaEdit.textField
            return 1
        
        
        #Rough variant calling with GATK
        self.printAttributes()
        
        #create transcriptome from GTF-File
        #startTime = Helper.getTime()
        #Helper.info(" [%s] Parsing Gene Data from %s" % (startTime.strftime("%c"),self.rnaEdit.params.gtfFile),self.rnaEdit.logFile,self.rnaEdit.textField)
        
        #duration = Helper.getTime() -startTime
        #Helper.info(" Finished parsing in %s" % (str(duration)),self.rnaEdit.logFile,self.rnaEdit.textField)
    
                
        vcfFile=self.rnaEdit.params.output+".vcf"
        cmd = ["java","-Xmx6G","-jar",self.rnaEdit.params.sourceDir + "GATK/GenomeAnalysisTK.jar", 
               "-T","UnifiedGenotyper","-R", self.rnaEdit.params.refGenome, "-glm", "SNP","-I", self.bamFile, 
               "-D", self.rnaEdit.params.dbsnp, "-o", vcfFile, "-metrics", self.rnaEdit.params.output+".snp.metrics", "-nt", self.rnaEdit.params.threads, "-l","ERROR",
               "-stand_call_conf", self.rnaEdit.params.standCall, "-stand_emit_conf", self.rnaEdit.params.standEmit,"-A", "Coverage", "-A", "AlleleBalance","-A", "BaseCounts"]
        #print cmd
        Helper.proceedCommand("Call variants", cmd, self.bamFile, vcfFile, self.rnaEdit)
        
        #check if file already exists
        if not os.path.isfile(self.rnaEdit.params.output+"noSNPs.vcf") or self.rnaEdit.params.overwrite==True:
            #read in initial SNPs
            variants = VariantSet(vcfFile,self.rnaEdit.logFile,self.rnaEdit.textField)

    
            '''delete SNPs from dbSNP'''
            variants.deleteOverlappsFromVcf(self.rnaEdit.params.dbsnp)
            
            '''delete variants from 1000 Genome Project'''
            variants.deleteOverlappsFromVcf(self.rnaEdit.params.omni)
            
            '''delete variants from UW exome calls'''
            variants.deleteOverlappsFromVcf(self.rnaEdit.params.esp)
            
            '''annotate all Variants'''
            #variants.annotateVariantDict(self.genome)
            
            '''save variants if something goes wrong'''
            variants.printVariantDict(self.rnaEdit.params.output+"noSNPs.vcf")
        else:
            if not os.path.isfile(self.rnaEdit.params.output+"noReadEdges.vcf"):
                variants = VariantSet(self.rnaEdit.params.output+"noSNPs.vcf",self.rnaEdit.logFile,self.rnaEdit.textField)
                
        
        
        if not os.path.isfile(self.rnaEdit.params.output+"noReadEdges.vcf") or self.rnaEdit.params.overwrite==True:
            '''erase artificial missmatches at read-edges from variants'''
            #self.removeEdgeMissmatches(variants, self.bamFile, self.rnaEdit.params.edgeDistance, 25)
            
            '''save variants if something goes wrong'''
            variants.printVariantDict(self.rnaEdit.params.output+"noReadEdges.vcf")
        else:
            variants = VariantSet(self.rnaEdit.params.output+"noReadEdges.vcf",self.rnaEdit.logFile,self.rnaEdit.textField)
            
            
        '''get non-Alu Variants'''
        nonAluVariants=copy(variants)
        nonAluVariants.variantDict=variants.getOverlappsFromBed(self.rnaEdit.params.aluRegions,getNonOverlapps=True)
        
        '''get Alu Variants'''
        aluVariants=copy(variants)
        aluVariants.variantDict=variants.getOverlappsFromBed(self.rnaEdit.params.aluRegions,getNonOverlapps=False)
        
        
        '''Read Genome'''
        self.genome = Genome(self.rnaEdit.params.gtfFile)
        
        
        #print out variants from Alu regions
        aluVariants.annotateVariantDict(self.genome)
        aluVariants.printVariantDict(self.rnaEdit.params.output+".alu.vcf")
        aluVariants.printGeneList(self.genome,self.rnaEdit.params.output+".alu.gvf", printSummary=True)
        
        ##############################################
        ###   proceed with non-Alu reads only!!!    ##
        ##############################################
        
        #erase variants from intronic splice junctions
        self.removeIntronicSpliceJunctions(nonAluVariants, self.genome)
        
        #erase variants from homopolymer runs
        self.removeHomopolymers(nonAluVariants,self.rnaEdit.params.output, 4)
        
        #do blat search
        blatOutfile = self.rnaEdit.params.output + "_blat"
        self.blatSearch(nonAluVariants, blatOutfile, 25, 1)
        
        #print nonAlu variants
        nonAluVariants.printVariantDict(self.rnaEdit.params.output+".nonAlu.vcf")
        nonAluVariants.printGeneList(self.genome,self.rnaEdit.params.output+".nonAlu.gvf", printSummary=True)
        
        #print nonAlu editing Sites
        self.deleteNonEditingBases(nonAluVariants)
        nonAluVariants.printVariantDict(self.rnaEdit.params.output+".editingSites.nonAlu.vcf")
        nonAluVariants.printGeneList(self.genome,self.rnaEdit.params.output+".editingSites.nonAlu.gvf",printSummary=True)
        nonAluVariants.createClusters(eps=50,minSamples=5)
        nonAluVariants.printVariantDict(self.rnaEdit.params.output+".editingSites.nonAlu.clusters")
        #print Alu editing Sites
        self.deleteNonEditingBases(aluVariants)
        aluVariants.printVariantDict(self.rnaEdit.params.output+".editingSites.alu.vcf")
        aluVariants.printGeneList(self.genome,self.rnaEdit.params.output+".editingSites.alu.gvf",printSummary=True)
        aluVariants.createClusters(eps=50,minSamples=5)
        aluVariants.printVariantDict(self.rnaEdit.params.output+".editingSites.alu.clusters")
        
        #combine alu and non Alu sites
        variants=aluVariants+nonAluVariants
        self.deleteNonEditingBases(variants)
        
        #print Final tables
        variants.printVariantDict(self.rnaEdit.params.output+".editingSites.vcf")
        variants.printGeneList(self.genome,self.rnaEdit.params.output+".editingSites.gvf",printSummary=True)
        variants.createClusters(eps=50,minSamples=5)
        variants.printClusters(self.rnaEdit.params.output+".editingSites.clusters")
        
        return 1


def checkDependencies(args):
    '''
    Checks the existence of the necessary packages and tools
    :param sourceDir: folder which contains all the software
    '''
    Helper.newline(1)
    Helper.info("CHECK DEPENDENCIES")
    
    #check if all tools are there
    if not os.path.isfile(args.sourceDir+"GATK/GenomeAnalysisTK.jar"):
        Helper.error("GenomeAnalysisTK.jar not found in %s" % args.sourceDir+"GATK/")
    if not os.path.isfile(args.sourceDir+"bedtools/fastaFromBed"):
        Helper.error("fastaFromBed not found in %s" % args.sourceDir+"bedtools/")
    if not os.path.isfile(args.sourceDir+"blat"):
        Helper.error("blat not found in %s" % args.sourceDir)
    if not os.path.isfile(args.sourceDir+"samtools"):
        Helper.error("samtools not found in %s" % args.sourceDir)
    if not os.system("java -version")==0:
        Helper.error("Java could not be found, Please install java")
    
    #check if all files are there
    if not os.path.isfile(args.RefGenome):
        Helper.error("Could not find Reference Genome in %s: " % args.RefGenome)
    
    #Files for GATK
    if not os.path.isfile(args.RefGenome.replace(".fastq",".dict")):
        Helper.error("Could not find %s" % args.RefGenome.replace(".fastq",".dict"))
        Helper.error("run: 'java -jar %s/picard-tools/CreateSequenceDictionary.jar R=%s  O= %s' to create it" % (args.sourceDir,args.RefGenome,args.RefGenome.replace(".fastq",".dict")))
    if not os.path.isfile(args.RefGenome+".fai"):
        Helper.error("Could not find %s.fai" % args.RefGenome)
        Helper.error("run: 'samtools faidx %s' to create it" % args.RefGenome)

    #SNP databases
    if not os.path.isfile(args.dbsnp):
        Helper.error("Could not find %s: " % args.dbsnp)
    if not os.path.isfile(args.hapmap):
        Helper.error("Could not find %s: " % args.hapmap)
    if not os.path.isfile(args.omni):
        Helper.error("Could not find %s: " % args.omni)
    if not os.path.isfile(args.esp):
        Helper.error("Could not find %s: " % args.esp)
    
    #region Files
    if not os.path.isfile(args.aluRegions):
        Helper.error("Could not find %s: " % args.aluRegions)
    if not os.path.isfile(args.gtfFile):
        Helper.error("Could not find %s: " % args.gtfFile)        
        
