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
from pysam import Samfile
from pysam import Fastafile

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
    
    
  
    def removeIntronicSpliceJunctions(self,variants,genome,distance=4): 
        '''
        remove variant near splice junctions and returns the other variants
        :param variants: VariantSet
        :param genome: object of the class Genome
        '''
        startTime=Helper.getTime()
        
        Helper.info(" [%s] remove Missmatches from the intronic splice junctions " % (startTime.strftime("%c")),self.rnaEdit.logFile,self.rnaEdit.textField)

        
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
        
        refGenome = "/media/Storage/databases/rnaEditor_annotations/human/human_g1k_v37.fasta"
        fastaFile = Fastafile(self.rnaEdit.params.refGenome)
        mmNumberTotal = len(variants.variantDict)
        #print temporary BedFile
        numberPassed=0
        for key in variants.variantDict.keys():
            chr,position,ref,alt = key
            startPos = position - distance if position >= distance else 0
            endpos = position + distance
            sequence = fastaFile.fetch(chr,startPos,endpos)
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
                    del variants.variantDict[key]
                except KeyError:
                    pass
            else:
                numberPassed+=1
                    
        #output statistics
        Helper.info("\t\t %d out of %d passed the Homopolymer-Filter" % (numberPassed, mmNumberTotal),self.rnaEdit.logFile,self.rnaEdit.textField)
        Helper.printTimeDiff(startTime,self.rnaEdit.logFile,self.rnaEdit.textField)
                
    '''do blat search (delete variants from reads that are not uniquely mapped)'''
    def blatSearch(self,variants, outFile, minBaseQual, minMissmatch):
        startTime=Helper.getTime()
        Helper.info(" [%s] Search non uniquely mapped reads" % (startTime.strftime("%c")),self.rnaEdit.logFile,self.rnaEdit.textField)
        
        bamFile= Samfile(self.bamFile,"rb")
        #create Fasta file for blat to remap the variant overlapping reads
        tempFasta = outFile + "_tmp.fa"
        if not os.path.isfile(tempFasta) or not os.path.getsize(tempFasta) > 0: #check if temFast exists and is not empty. If it exist it will not be created again
            tempFastaFile=open(tempFasta,"w+")
            mmNumberTotal = len(variants.variantDict)
            
            #############################################
            #########    CREATE FASTA FILE        #######
            #############################################
            Helper.info(" [%s] Create fasta file for blat " % (startTime.strftime("%c")),self.rnaEdit.logFile,self.rnaEdit.textField)
            counter=1
            
            if len(variants.variantDict.keys()) == 0:
                Helper.error("No Variants left" ,self.rnaEdit.logFile,self.rnaEdit.textField)
            
            for varKey in variants.variantDict.keys(): 
                variant=variants.variantDict[varKey]
                varPos=variant.position-1
                iter = bamFile.pileup(variant.chromosome, variant.position-1, variant.position)               
                alignements=[]
                for x in iter:
                    if x.pos == varPos:
                        #loop over reads of that position
                        for pileupread in x.pileups:
                            if not pileupread.is_del and not pileupread.is_refskip:
                                if pileupread.alignment.query_sequence[pileupread.query_position] == variant.alt and pileupread.alignment.query_qualities[pileupread.query_position]>=minBaseQual:
                                #if pileupread.alignment.query_sequence[pileupread.query_position] == variant.alt:
                                    alignements.append(pileupread.alignment.seq)
                
                if len(alignements)>=minMissmatch:
                    missmatchReadCount=0
                    for sequence in alignements:
                        tempFastaFile.write("> "+variant.chromosome+"-"+str(variant.position)+"-"+variant.ref+"-"+variant.alt+"-"+str(missmatchReadCount)+"\n"+sequence+"\n")
                        missmatchReadCount += 1
                        
                counter += 1
                if counter % 1000 == 0:
                    sys.stdout.write("\r" + str(counter) + " of " + str(mmNumberTotal) + " variants done")
                    Helper.info(str(counter) + " of " + str(mmNumberTotal) + " variants done", self.rnaEdit.logFile,self.rnaEdit.textField)
                    sys.stdout.flush()
        
            Helper.info("\n created fasta file " + tempFasta,self.rnaEdit.logFile,self.rnaEdit.textField)
            Helper.printTimeDiff(startTime,self.rnaEdit.logFile,self.rnaEdit.textField)
            tempFastaFile.close()
                
        #############################
        #####   do blat search  #####
        #############################
        pslFile=outFile+".psl"
        if not os.path.isfile(pslFile) or not os.path.getsize(pslFile) > 0:
            cmd = [self.rnaEdit.params.sourceDir+"blat","-stepSize=5","-repMatch=2253", "-minScore=20","-minIdentity=0","-noHead", self.rnaEdit.params.refGenome, tempFasta, pslFile]
            #print cmd
            Helper.proceedCommand("do blat search for unique reads",cmd,tempFasta, "None", self.rnaEdit)
        Helper.info(" [%s] Blat finished" % (startTime.strftime("%c")),self.rnaEdit.logFile,self.rnaEdit.textField)
        Helper.info(" [%s] Parse Blat output to look for non uniquely mapped reads" % (startTime.strftime("%c")),self.rnaEdit.logFile,self.rnaEdit.textField)    
        
        if not os.path.isfile(outFile):
            #open psl file
            pslFile=open(pslFile,"r")
            blatDict={}
            
            for line in pslFile: #summarize the blat hits
                pslFields = line.split()
                chr,pos,ref,alt,mmReadCount = pslFields[9].split("-")
                varTuple=(chr,int(pos),ref,alt)
                try:
                    blatScore = [pslFields[0], pslFields[13], pslFields[17], pslFields[18], pslFields[20]] # #of Matches, targetName, blockCount, blockSize, targetStarts 
                except IndexError:
                    Helper.warning("Not enough Values in '%s' (Skip)" % line, self.rnaEdit.logFile,self.rnaEdit.textField)
                    continue
                if varTuple in blatDict:
                    blatDict[varTuple] = blatDict[varTuple] + [blatScore]
                else:
                    blatDict[varTuple] = [blatScore]

            siteDict = {}
            discardDict = {}
            Helper.info(" [%s] Analyse Blat hits (Slow)" % (startTime.strftime("%c")),self.rnaEdit.logFile,self.rnaEdit.textField)    
        
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
                if chr == largestScoreLine[1] and scoreArray[1] < scoreArray[0]*0.95: #check if same chromosome and hit is lower the 95 percent of first hit
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
            
            ##############################################################################        
            #####        loop through variants and delete invalid variants          ######
            ##############################################################################
            Helper.info(" [%s] Deleting invalid variants" % (startTime.strftime("%c")),self.rnaEdit.logFile,self.rnaEdit.textField)    
        
            
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
                elif numberBlatReads < numberDiscardReads: #check if more reads fit the blat criteria than not
                    mmReadsSmallerDiscardReads+=1    
                mmNumberTotal+=1
            
            if self.rnaEdit.params.keepTemp == False:
                os.remove(tempFasta)
                os.remove(pslFile.name)
            
            #output statistics
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
            if os.path.isfile(self.rnaEdit.params.output+"_tmp.bed"):
                os.remove(self.rnaEdit.params.output+"_tmp.bed")
            if os.path.isfile(self.rnaEdit.params.output+"_tmp.tsv"):
                os.remove(self.rnaEdit.params.output+"_tmp.tsv")
            if os.path.isfile(self.rnaEdit.params.output+".noBlat.vcf.psl"):
                os.remove(self.rnaEdit.params.output+".noBlat.vcf.psl")
            if os.path.isfile(self.rnaEdit.params.output+"noBlat.vcf_tmp.fa"):
                os.remove(self.rnaEdit.params.output+"noBlat.vcf_tmp.fa")
            if os.path.isfile(self.rnaEdit.params.output+".noSNPs.vcf"):
                os.remove(self.rnaEdit.params.output+".noSNPs.vcf")
            if os.path.isfile(self.rnaEdit.params.output+".noSpliceJunction.vcf"):
                os.remove(self.rnaEdit.params.output+".noSpliceJunction.vcf")
            if os.path.isfile(self.rnaEdit.params.output+".noHomo.vcf"):
                os.remove(self.rnaEdit.params.output+".noHomo.vcf")
            if os.path.isfile(self.rnaEdit.params.output+".noReadEdges.vcf"):
                os.remove(self.rnaEdit.params.output+".noReadEdges.vcf")
            if os.path.isfile(self.rnaEdit.params.output+".nonAlu.vcf"):
                os.remove(self.rnaEdit.params.output+".nonAlu.vcf")
            if os.path.isfile(self.rnaEdit.params.output+".nonAlu.noSpliceSites.vcf"):
                os.remove(self.rnaEdit.params.output+".nonAlu.noSpliceSites.vcf")
            if os.path.isfile(self.rnaEdit.params.output+".nonAlu.noSpliceSites.noHomo.vcf"):
                os.remove(self.rnaEdit.params.output+".nonAlu.noSpliceSites.noHomo.vcf")
            
    def startAnalysis(self):
        '''Proceeds all the steps to detect editing Sites from a bam File
        
        @return: 0 on success and 1 if analysis was canceled by user
        '''
        
        '''check if result file already exists''' 
        if os.path.isfile(self.rnaEdit.params.output+".editingSites.clusters") and self.rnaEdit.params.overwrite==False:
            Helper.status("\t [SKIP] Final result file already exist",self.rnaEdit.logFile,self.rnaEdit.textField,"green")
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
        
        
        #################################
        ###   Delete known SNPs!!!    ###
        #################################
        #check if file already exists
        if not os.path.isfile(self.rnaEdit.params.output+".noSNPs.vcf") or self.rnaEdit.params.overwrite==True:
            #read in initial SNPs
            variants = VariantSet(vcfFile,self.rnaEdit.logFile,self.rnaEdit.textField)

    
            '''delete SNPs from dbSNP'''
            variants.deleteOverlapsFromVcf(self.rnaEdit.params.dbsnp)
            
            '''delete variants from 1000 Genome Project'''
            if self.rnaEdit.params.omni != "None":
                variants.deleteOverlapsFromVcf(self.rnaEdit.params.omni)
            
            '''delete variants from UW exome calls'''
            if self.rnaEdit.params.esp != "None":
                variants.deleteOverlapsFromVcf(self.rnaEdit.params.esp)
            
            '''annotate all Variants'''
            #variants.annotateVariantDict(self.genome)
            
            '''save variants if something goes wrong'''
            variants.printVariantDict(self.rnaEdit.params.output+".noSNPs.vcf")
        else:
            if not os.path.isfile(self.rnaEdit.params.output+".noReadEdges.vcf"):
                variants = VariantSet(self.rnaEdit.params.output+".noSNPs.vcf",self.rnaEdit.logFile,self.rnaEdit.textField)
                
        
        ###############################################
        ###   Delete variants from read edges!!!    ###
        ###############################################
        if not os.path.isfile(self.rnaEdit.params.output+".noReadEdges.vcf") or self.rnaEdit.params.overwrite==True:
            '''erase artificial missmatches at read-edges from variants'''
            variants.removeEdgeMismatches(self.bamFile, self.rnaEdit.params.edgeDistance, 25)
            #self.removeEdgeMissmatches(variants, self.bamFile, self.rnaEdit.params.edgeDistance, 25)
            
            '''save variants if something goes wrong'''
            variants.printVariantDict(self.rnaEdit.params.output+".noReadEdges.vcf")
        else:
            if not os.path.isfile(self.rnaEdit.params.output+".alu.vcf") or not os.path.isfile(self.rnaEdit.params.output+".nonAlu.vcf"):
                variants = VariantSet(self.rnaEdit.params.output+".noReadEdges.vcf",self.rnaEdit.logFile,self.rnaEdit.textField)
            
        
        ###############################################
        ###   split Alu- and non-Alu Variants!!!    ###
        ###############################################  
        
        if (not os.path.isfile(self.rnaEdit.params.output+".alu.vcf") or not os.path.isfile(self.rnaEdit.params.output+".nonAlu.vcf")) or self.rnaEdit.params.overwrite==True:
            '''get non-Alu Variants'''
            nonAluVariants=copy(variants)
            #nonAluVariants.variantDict=variants.getOverlapsFromBed(self.rnaEdit.params.aluRegions,getNonOverlaps=True)
            
            '''get Alu Variants'''
            aluVariants=copy(variants)
            #aluVariants.variantDict=variants.getOverlapsFromBed(self.rnaEdit.params.aluRegions,getNonOverlaps=False)
            aluVariants.variantDict,nonAluVariants.variantDict = variants.splitByBed(self.rnaEdit.params.aluRegions)
            aluVariants.printVariantDict(self.rnaEdit.params.output+".alu.vcf")
            nonAluVariants.printVariantDict(self.rnaEdit.params.output+".nonAlu.vcf")
        else:     
            aluVariants = VariantSet(self.rnaEdit.params.output+".alu.vcf",self.rnaEdit.logFile,self.rnaEdit.textField)
            if not os.path.isfile(self.rnaEdit.params.output+".noSpliceJunction.vcf"):
                nonAluVariants = VariantSet(self.rnaEdit.params.output+".nonAlu.vcf",self.rnaEdit.logFile,self.rnaEdit.textField)
        
        
        #print out variants from Alu regions
        #
        
        ##############################################
        ###   proceed with non-Alu reads only!!!    ##
        ##############################################
        
        ##############################################
        ###   Remove intronic Splice junction!!!    ##
        ##############################################
        self.genome = Genome(self.rnaEdit.params.gtfFile,self.rnaEdit.logFile,self.rnaEdit.textField)
        #erase variants from intronic splice junctions
        if not os.path.isfile(self.rnaEdit.params.output+".noSpliceJunction.vcf") or self.rnaEdit.params.overwrite==True:
            self.removeIntronicSpliceJunctions(nonAluVariants, self.genome)
            nonAluVariants.printVariantDict(self.rnaEdit.params.output+".noSpliceJunction.vcf")
        else:
            if not os.path.isfile(self.rnaEdit.params.output+".noHomo.vcf"):
                nonAluVariants = VariantSet(self.rnaEdit.params.output+".noSpliceJunction.vcf",self.rnaEdit.logFile,self.rnaEdit.textField)
        
        
        ##############################################
        ### erase variants from homopolymers!!! ##
        ##############################################
        if not os.path.isfile(self.rnaEdit.params.output+".noHomo.vcf") or self.rnaEdit.params.overwrite==True:
            self.removeHomopolymers(nonAluVariants,self.rnaEdit.params.output, 4)
            nonAluVariants.printVariantDict(self.rnaEdit.params.output+".noHomo.vcf")
        else:
            if not os.path.isfile(self.rnaEdit.params.output+".noBlat.vcf"):
                nonAluVariants = VariantSet(self.rnaEdit.params.output+".noHomo.vcf",self.rnaEdit.logFile,self.rnaEdit.textField)
            
        
        ##############################################
        ###     erase duplicate mapped reads!!!     ##
        ##############################################
        if not os.path.isfile(self.rnaEdit.params.output+".noBlat.vcf") or self.rnaEdit.params.overwrite==True:
            blatOutfile = self.rnaEdit.params.output + ".noBlat.vcf"
            self.blatSearch(nonAluVariants, blatOutfile, 25, 2)
            
            #print nonAlu variants
            nonAluVariants.printVariantDict(self.rnaEdit.params.output+".noBlat.vcf")
        else:
            if not os.path.isfile(self.rnaEdit.params.output+".editingSites.nonAlu.vcf"):
                nonAluVariants = VariantSet(self.rnaEdit.params.output+".noBlat.vcf",self.rnaEdit.logFile,self.rnaEdit.textField)
                #nonAluVariants.deleteNonEditingBases()
                #nonAluVariants.printVariantDict(self.rnaEdit.params.output+".editingSites.nonAlu.vcf")
            else:
                nonAluVariants = VariantSet(self.rnaEdit.params.output+".editingSites.nonAlu.vcf",self.rnaEdit.logFile,self.rnaEdit.textField)
        #nonAluVariants.printGeneList(self.genome,self.rnaEdit.params.output+".nonAlu.gvf", printSummary=True)
        
        #print nonAlu editing Sites
        nonAluVariants.deleteNonEditingBases()
        nonAluVariants.annotateVariantDict(self.genome)
        nonAluVariants.printVariantDict(self.rnaEdit.params.output+".editingIslands.bed")
        nonAluVariants.printGeneList(self.genome,self.rnaEdit.params.output+".editingSites.nonAlu.gvf",printSummary=True)
        nonAluVariants.createClusters(eps=50,minSamples=5)
        nonAluVariants.printClusters(self.rnaEdit.params.output+".editingIslands.bed")
        #print Alu editing Sites
        aluVariants.deleteNonEditingBases()
        aluVariants.annotateVariantDict(self.genome)
        aluVariants.printVariantDict(self.rnaEdit.params.output+".editingSites.alu.vcf")
        aluVariants.printGeneList(self.genome,self.rnaEdit.params.output+".editingSites.alu.gvf",printSummary=True)
        aluVariants.createClusters(eps=50,minSamples=5)
        aluVariants.printClusters(self.rnaEdit.params.output+".editingIslands.bed")
        
        #combine alu and non Alu sites
        variants=aluVariants+nonAluVariants
        variants.deleteNonEditingBases()
        
        #print Final tables
        '''Read Genome'''
        
        variants.annotateVariantDict(self.genome)
        
        variants.printVariantDict(self.rnaEdit.params.output+".editingSites.vcf")
        variants.printGeneList(self.genome,self.rnaEdit.params.output+".editingSites.gvf",printSummary=True)
        variants.createClusters(eps=50,minSamples=5)
        variants.printClusters(self.rnaEdit.params.output+".editingIslands.bed")
        
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

    
    #region Files
    if not os.path.isfile(args.aluRegions):
        Helper.error("Could not find %s: " % args.aluRegions)
    if not os.path.isfile(args.gtfFile):
        Helper.error("Could not find %s: " % args.gtfFile)        
        
