'''
Created on May 23, 2013

@author: david
'''

import argparse, multiprocessing, os, sys, re
from Helper import Helper
from genericpath import exists

from VariantSet import VariantSet

from Genome import Genome
from copy import copy



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
                 aluRegions, gtfFile, outfilePrefix="default",
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
        self.genome = gtfFile
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
        
        #self.features = Helper.readGeneFeatures(self.genome)
        
        
        self.logFile=open(self.outfilePrefix + ".log","a")
        if self.debug==True:
            self.printAttributes()
        
        self.checkDependencies()
        
        #create transcriptome from GTF-File
        self.genome = Genome(gtfFile)
        
    
    
    
    '''check if all the needed files are there'''
    def checkDependencies(self):
        #TODO: check for index Files
        if not os.path.exists(self.dbsnp):
            Exception("dbSNP File: Not found!!!")
        if not os.path.exists(self.bamFile):
            Exception(".bam File: Not found!!!")
        if not os.path.exists(self.hapmap):
            Exception("HapMap variant File File: Not found!!!")
        if not os.path.exists(self.esp):
            Exception("ESP variant File File: Not found!!!")
        if not os.path.exists(self.aluRegions):
            Exception("AluRegion File File: Not found!!!")
        if not os.path.exists(self.refGenome):
            Exception("reference Genome File: Not found!!!")
        
    '''delete variants from Bam file which appear near read edges'''
    def removeEdgeMissmatches(self,variants,bamFile,minDistance, minBaseQual):
        #Loop through vcf-File
            #call overlapping reads with samtools view
            #loop over reads
                #discard variants wich appear ONLY near edges
                #write the rest to the output file
        startTime=Helper.getTime()
        
        counter=0    
        
        num_lines = len(variants.variantDict)
        Helper.info(" [%s] remove Missmatches from the first %s bp from read edges" % (startTime.strftime("%c"),str(minDistance)))
        
        for varKey in variants.variantDict.keys():
            variant = variants.variantDict[varKey]
            snpPos = variant.position
            position=variant.chromosome+":" + str(snpPos) + "-" + str(snpPos) 
            #line[1]+"-"+line[1]
            keepSNP=False
            
            command = [self.sourceDir+"samtools", "view", bamFile, position] 
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
            if counter % 100 == 0: #print out current status
                Helper.status(str(counter) + " of " + str(num_lines) + " missmatches finished")
        
    
    def removeIntronicSpliceJunctions(self,variants,genome,distance=4): 
        '''
        remove variant near splice junctions and returns the other variants
        :param variants: VariantSet
        :param genome: object of the class Genome
        '''
        startTime=Helper.getTime()
        Helper.info(" [%s] remove Missmatches from the intronic splice junctions " % (startTime.strftime("%c")))
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
                            
        Helper.printTimeDiff(startTime)
        
    '''remove missmatches from homopolymers'''
    def removeHomopolymers(self,variants,outFile,distance):
        startTime=Helper.getTime()
        Helper.info(" [%s] remove Missmatches from homopolymers " % (startTime.strftime("%c")))
        
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
        cmd=[self.sourceDir+"bedtools/fastaFromBed", "-name", "-tab", "-fi", self.refGenome, "-bed", tempBedFile.name, "-fo", tempSeqFile]
        Helper.proceedCommand("catch surrounding sequences of Missmatches", cmd, tempBedFile.name, tempSeqFile, self.logFile, self.overwrite)
        
        mmNumberTotal = len(variants.variantDict)
        
        #read sequence file
        tempSeqFile= open()
        for line in tempSeqFile:
            siteNuc,sequence = line.split()
            try:
                chr,position,ref,alt = siteNuc.split(",")
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
                del variants.variantDict[(chr,int(position),ref,alt)]
                
        #output statistics
        Helper.info("\t\t %d out of %d passed the Homopolymer-Filter" % (mmNumberTotal, mmNumberTotal))
        Helper.printTimeDiff(startTime)
        
        if self.keepTemp == False:
            os.remove(tempBedFile.name)
            os.remove(tempSeqFile.name)   
                
    '''do blat search (delete variants from reads that are not uniquely mapped)'''
    def blatSearch(self,variants, outFile, minBaseQual, minMissmatch):
        startTime=Helper.getTime()
        Helper.info(" [%s] Search non uniquely mapped reads" % (startTime.strftime("%c")))
        
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
    
                samout = Helper.getCommandOutput([self.sourceDir+"samtools", "view", "-F", "1024", self.bamFile, samPos]).splitlines() #-F 1024 to filter out duplicate reads
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
        
            Helper.info("\n created fasta file " + tempFasta)
            Helper.printTimeDiff(startTime)
            tempFastaFile.close()
                
        
        #do blat search
        pslFile=outFile+".psl"
        if not os.path.isfile(pslFile) or not os.path.getsize(pslFile) > 0:
            cmd = [self.sourceDir+"blat","-stepSize=5","-repMatch=2253", "-minScore=20","-minIdentity=0","-noHead", self.refGenome, tempFasta, pslFile]
            #print cmd
            Helper.proceedCommand("do blat search for unique reads",cmd,tempFasta, "None", self.logFile, self.overwrite)
        
        Helper.info(" [%s] look for non uniquely mapped reads by blat" % (startTime.strftime("%c")))    
        
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
            
            if self.keepTemp == False:
                os.remove(tempFasta)
                os.remove(pslFile.name)
            
            #output statisticsttkkg
            mmPassedNumber=mmNumberTotal-(mmNumberTooSmall+mmReadsSmallerDiscardReads)
            
            Helper.info("\t\t %d out of %d passed blat criteria" % (mmPassedNumber, mmNumberTotal))
            Helper.info("\t\t %d Missmatches had fewer than %d missmatching-Reads." % (mmNumberTooSmall, minMissmatch))
            Helper.info("\t\t %d Missmatches had more missaligned reads than correct ones." % (mmReadsSmallerDiscardReads))
            
        Helper.printTimeDiff(startTime)

            
    def __del__(self):
        if self.keepTemp==False:
            pass
            #os.remove(self.outfilePrefix+".vcf")
            #os.remove(self.outfilePrefix+".no_dbsnp.vcf")
            #os.remove(self.outfilePrefix+".no_dbsnp.no_1000genome.vcf")
            #os.remove(self.outfilePrefix+".no_dbsnp.no_1000genome.no_esp.vcf")
            #os.remove(self.outfilePrefix+".no_dbsnp.no_1000genome.no_esp.noStartMM.vcf")
            #os.remove(self.outfilePrefix+".nonAlu.vcf")
            #os.remove(self.outfilePrefix+".nonAlu.noSpliceSites.vcf")
            #os.remove(self.outfilePrefix+".nonAlu.noSpliceSites.noHomo.vcf")
            
    def deleteNonEditingBases(self,variants):
        startTime=Helper.getTime()
        Helper.info("Delete non Editing Bases (keep only T->C and A->G)")
        
        for varTuple in variants.variantDict.keys():
            chr,pos,ref,alt = varTuple
            if (ref =="A" and alt == "G") or (ref=="T" and alt=="C"):
                pass
            else:
                del variants.variantDict[varTuple]
    
    def start(self):
        #Rough variant calling with GATK
        vcfFile=self.outfilePrefix+".vcf"
        cmd = ["java","-Xmx6G","-jar",self.sourceDir + "GATK/GenomeAnalysisTK.jar", 
               "-T","UnifiedGenotyper","-R", self.refGenome, "-glm", "SNP","-I", self.bamFile, 
               "-D", self.dbsnp, "-o", vcfFile, "-metrics", self.outfilePrefix+".snp.metrics", "-nt", self.threads, "-l","ERROR",
               "-stand_call_conf", self.standCall, "-stand_emit_conf", self.standEmit,"-A", "Coverage", "-A", "AlleleBalance","-A", "BaseCounts"]
        #print cmd
        Helper.proceedCommand("Call variants", cmd, self.bamFile, vcfFile, self.logFile, self.overwrite)
        
        #read in initial SNPs
        variants = VariantSet(vcfFile)
        
        
        #annotate all Variants
        variants.annotateVariantDict(self.genome)
        #print len(rawSnps)
        '''
        #delete SNPs from dbSNP
        variants.deleteOverlappsFromVcf(self.dbsnp)
        #print len(noDbsnp)
       
        #delete variants from 1000 Genome Project
        variants.deleteOverlappsFromVcf(self.omni)
        #print len(noOmni)
        
        #delete variants from UW exome calls
        variants.deleteOverlappsFromVcf(self.esp)
        #print len(noEsp)
        
        #erase artificial missmatches at read-edges from variants
        self.removeEdgeMissmatches(variants, self.bamFile, self.edgeDistance, 25)
        '''
        nonAluVariants=copy(variants)
        nonAluVariants.variantDict=variants.getOverlappsFromBed(self.aluRegions,getNonOverlapps=True)
        
        aluVariants=copy(variants)
        aluVariants.variantDict=variants.getOverlappsFromBed(self.aluRegions,getNonOverlapps=False)
        
        
        
        #print out variants from Alu regions
        
        aluVariants.printVariantDict(self.outfilePrefix+".alu.vcf")
        aluVariants.printGeneList(self.genome,self.outfilePrefix+".alu.gvf", printSummary=True)
        
        #proceed with non-Alu reads only!!!
        #erase variants from intronic splice junctions
        self.removeIntronicSpliceJunctions(nonAluVariants, self.genome)
        
        
        #erase variants from homopolymer runs
        self.removeHomopolymers(nonAluVariants,self.outfilePrefix, 4)
        
        #do blat search
        blatOutfile = self.outfilePrefix + "_blat"
        self.blatSearch(nonAluVariants, blatOutfile, 25, 1)
        
        #print nonAlu variants
        nonAluVariants.printVariantDict(self.outfilePrefix+".nonAlu.vcf")
        nonAluVariants.printGeneList(self.genome,self.outfilePrefix+".nonAlu.gvf", printSummary=True)
        
        variants=aluVariants+nonAluVariants
        self.deleteNonEditingBases(variants)
        
        
        variants.printVariantDict(self.outfilePrefix+".editingSites.vcf")
        variants.printGeneList(self.genome,self.outfilePrefix+".editingSites.gvf", printSummary=True)
        #combine alu and non Alu sites
        
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='output vatiants from a given .bam file.')
    parser.add_argument('-i', '--input', metavar='bam-File', type=argparse.FileType('r'), help='Input bam file from which variants should be called', required=True)
    parser.add_argument("-r", "--RefGenome", metavar='Fasta-File', help="File that contains the reference sequences", type=argparse.FileType('r'), default='/media/media/databases/human/human_g1k_v37.fa')
    parser.add_argument('-s', '--dbsnp', help='SNP database (dbSNP) in VCF format (downloaded from the GATK homepage)', type=argparse.FileType('r'), default='/media/media/databases/human/dbsnp_135.b37.vcf')
    parser.add_argument('-m', '--hapmap', help='hapmap database in vcf format (see GATK homepage)', type=argparse.FileType('r'), default='/media/media/databases/human/hapmap_3.3.b37.sites.vcf')
    parser.add_argument('-g', '--omni', help='1000 Genome variants in vcf format (see GATK homepage)', type=argparse.FileType('r'), default='/media/media/databases/human/1000G_omni2.5.b37.sites.vcf')
    parser.add_argument('-e', '--esp', help='Exome Sequencing Project variants', type=argparse.FileType('r'), default='/media/media/media/databases/human/NHLBI_Exome_Sequencing_Project_6500SI.vcf')
    parser.add_argument('-a', '--aluRegions', help='Alu-Regions downloaded fron the UCSC table browser', type=argparse.FileType('r'), default='/media/media/media/databases/human/hg19/rna-editing/Alu_repeats_noChr.bed')
    parser.add_argument('-gtf', '--gtfFile', help='Gene annotation File in GTF format', type=argparse.FileType('r'))
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
                          args.aluRegions, args.gtfFile, args.output, 
                          args.sourceDir, args.threads, args.standCall, 
                          args.standEmit, args.edgeDistance, args.keepTemp, 
                          args.overwrite)