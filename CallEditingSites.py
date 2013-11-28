'''
Created on May 23, 2013

@author: david
'''

import argparse, multiprocessing, os, sys, re
from Helper import Helper
from genericpath import exists

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
        
        print "\t sourceDir:" + self.sourceDir
        print "\t threads:" + self.threads
        print "\t StandCall:" + self.standCall
        print "\t standEmit:" + self.standEmit
        print "\t keepTemp:" + str(self.keepTemp)
        print "\t overwrite:" + str(self.overwrite)
        print

    def __init__(self,bamFile,refGenome,dbsnp,hapmap,omni,esp, aluRegions, outfilePrefix="default",sourceDir="/usr/local/bin/",
                 threads=multiprocessing.cpu_count()-1,standCall=0,standEmit=0, edgeDistance=6,
                 keepTemp=False,overwrite=False):
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
            snpPos=line[1]
            position=line[0]+":"+snpPos+"-"+snpPos
            keepSNP=False
            
            #print position, str(minDistance)
            #samout= os.system("samtools view " + bamFile + " " + position)
            samout = Helper.getCommandOutput("samtools view " + bamFile + " " + position).splitlines()
            for samLine in samout:
                samfields=samLine.split()
                flag,startPos,mapQual,cigar,sequence,seqQual = samfields[1],samfields[3],samfields[4],samfields[5],samfields[9],samfields[10]
                readPos=1
                mmReadPos=0
                cigarNums=re.split("[MIDNSHP]", cigar)[:-1]
                cigarLetters=re.split("[0-9]+",cigar)[1:]
                
                
                for i in range(len(cigarLetters)): #parse over single read
                    if cigarLetters[i] in {"I","S","H"}: #Insertion, Soft Clipping and Hard Clipping
                        readPos = readPos + int(cigarNums[i])
                    elif cigarLetters[i] in {"D","N"}: #Deletions and skipped Regions
                        startPos = startPos + int(cigarNums[i])
                    elif cigarLetters[i] in {"M"}: #Matches
                        for j in range(cigarNums[i]):
                            if startPos == snpPos:
                                mmReadPos = readPos 
                            readPos += 1
                            startPos += 1
                
                if mmReadPos != 0:
                                   
                    edgeDistance = int(snpPos) - int(startPos)
                
                    #only remove the snps from first 6 bases
                    revStrand = int(flag) & 16
                    if (revStrand == 0 and mmReadPos > minDistance) or (revStrand == 16 and mmReadPos < readPos - minDistance):
                        if(ord(seqQual[mmReadPos]) > minBaseQual + 33): #check for quality of the base
                            keepSNP=True
                    #print "   ".join([str(revStrand),str(keepSNP),str(edgeDistance),str(len(sequence)),flag,startPos,mapQual,cigar,sequence,seqQual])
                #print distance
            
            
            if keepSNP:
                outFile.write("\t".join(line))  #print SNP
            counter+=1
            if counter % 1000 == 0:
                sys.stdout.write("\r" + str(counter) + " of " + num_lines + " missmatches finished")
                sys.stdout.flush()
        duration=Helper.getTime()-startTime
        print >> self.logFile, "\t[DONE]" + " Duration [" + str(duration) + "]"
        self.logFile.flush()
        print "\t[DONE]" + " Duration [" + str(duration) + "]"    
            #sys.exit(0)
                
    def __del__(self):
        pass
    
    def start(self):
        #Rough variant calling with GATK
        vcfFile=self.outfilePrefix+".vcf"
        cmd = ["java","-Xmx4G","-jar",self.sourceDir + "GATK/GenomeAnalysisTK.jar", 
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
        nonAlu = self.outfilePrefix + ".nonAluEditingSites.vcf"
        alu = self.outfilePrefix + ".aluEditingSites.vcf"
        cmd =  [self.sourceDir+"bedtools/intersectBed","-v","-a",noStartMissmatches,"-b",self.aluRegions]
        Helper.proceedCommand("write variants from non-alu regions",cmd, noStartMissmatches, nonAlu, self.logFile, self.overwrite)  # write nonAlu-Regions
        cmd =  [self.sourceDir+"bedtools/intersectBed","-a",noStartMissmatches,"-b",self.aluRegions]
        Helper.proceedCommand("write variants from alu regions",cmd, noStartMissmatches, alu, self.logFile, self.overwrite) #write alu-regions
        
        #erase variants from intronic splice junctions
        
        
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='output vatiants from a given .bam file.')
    parser.add_argument('-i', '--input', metavar='bam-File', type=argparse.FileType('r'), help='Input bam file from which variants should be called', required=True)
    parser.add_argument("-r", "--RefGenome", metavar='Fasta-File', help="File that contains the reference sequences", type=argparse.FileType('r'), default='/media/media/databases/human/human_g1k_v37.fasta')
    parser.add_argument('-s', '--dbsnp', help='SNP database (dbSNP) in VCF format (downloaded from the GATK homepage)', type=argparse.FileType('r'), default='/media/media/databases/human/dbsnp_135.b37.vcf')
    parser.add_argument('-m', '--hapmap', help='hapmap database in vcf format (see GATK homepage)', type=argparse.FileType('r'), default='/media/media/databases/human/hapmap_3.3.b37.sites.vcf')
    parser.add_argument('-g', '--omni', help='1000 Genome variants in vcf format (see GATK homepage)', type=argparse.FileType('r'), default='/media/media/databases/human/1000G_omni2.5.b37.sites.vcf')
    parser.add_argument('-e', '--esp', help='Exome Sequencing Project variants', type=argparse.FileType('r'), default='/media/media/media/databases/human/NHLBI_Exome_Sequencing_Project_6500SI.vcf')
    parser.add_argument('-a', '--aluRegions', help='Alu-Regions downloaded fron the UCSC table browser', type=argparse.FileType('r'), default='/media/media/media/databases/human/hg19/rna-editing/Alu_repeats.bed')
    parser.add_argument('-o', '--output', metavar='output-prefix', type=str,help='prefix that is written in front of the output files', default="default")
    parser.add_argument('-d', '--sourceDir', help='- Directory to all the tools [default: bin/]', default='bin/', type=Helper.readable_dir)
    parser.add_argument('-t', '--threads', help='number of threads', type=int, default=multiprocessing.cpu_count()-1)
    parser.add_argument('-sc', '--standCall', help='-The minimum phred-scaled confidence threshold at which variants should be considered as true (int) [0]', type=int, default=0)
    parser.add_argument('-se', '--standEmit', help=' The minimum phred-scaled confidence threshold at which variants should be emitted (int)[0]', type=int, default=0)
    parser.add_argument('-ed', '--edgeDistance', help='The minimum edge distance of the SNPs', type=int, default=6)
    parser.add_argument('--keepTemp', help='keep the intermediate Files [False]', action='store_true', default=False)
    parser.add_argument('--overwrite', help='overwrite existing Files [False]', action='store_true', default=False)
    
    args = parser.parse_args()

    call=CallEditingSites(args.input.name, args.RefGenome.name, args.dbsnp.name, args.hapmap.name, args.omni.name, args.esp, args.output, args.sourceDir, args.threads, args.standCall, args.standEmit, args.edgeDistance, args.keepTemp, args.overwrite)