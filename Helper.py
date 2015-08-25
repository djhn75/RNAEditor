'''
Created on May 22, 2013

@author: david
'''

from datetime import datetime
import argparse, sys, os, subprocess, errno
from collections import defaultdict, OrderedDict
import traceback
import ui
import numpy as np
import matplotlib.pyplot as plt




class Parameters():
    '''
    Reads and saves the default values from the configuration File 'configurtion.txt'
    '''

    
    def __init__(self,source="configuration.txt"):
        '''
        creates an parameter object which contains all the parameters for RnaEditor
        :param source: either a QWidget or a textFile from where the parameters are read from
        '''
        if type(source)==str:
            self.readDefaultsFromFile(source)
        elif isinstance(source, ui.InputTab.InputTab):
            self.getParametersFromInputTab(source)
        else:
            Helper.error("Parameter source has wrong Type [str or QWidget]")
    
    
    def getParametersFromInputTab(self,inputTab):
        '''
        get the Parameters and update the default Parameters from the Default class 
        '''
        self.refGenome = str(inputTab.refGenomeTextBox.text())
        self.gtfFile = str(inputTab.gtfFileTextBox.text())
        self.dbsnp = str(inputTab.dbsnpTextBox.text())
        self.hapmap = str(inputTab.hapmapTextBox.text())
        self.omni = str(inputTab.omniTextBox.text())
        self.esp = str(inputTab.espTextBox.text())
        self.aluRegions = str(inputTab.aluRegionsTextBox.text())
        self.output = str(inputTab.outputTextBox.text())
        self.sourceDir = str(inputTab.sourceDirTextBox.text())
        
        self.threads = str(inputTab.threadsSpinBox.value())
        self.maxDiff = str(inputTab.maxDiffSpinBox.value())
        self.seedDiff = str(inputTab.seedSpinBox.value())
        self.standCall = str(inputTab.standCallSpinBox.value())
        self.standEmit = str(inputTab.standEmitSpinBox.value())
        self.edgeDistance = str(inputTab.edgeDistanceSpinBox.value())
        self.paired = inputTab.pairedCheckBox.isChecked()
        self.overwrite = inputTab.overwriteCheckBox.isChecked()
        self.keepTemp = inputTab.keepTempCheckBox.isChecked()
        
            
    def readDefaultsFromFile(self,file):
        try:
            confFile = open(file)
        except IOError:
            Helper.error("Unable to open configuration file")
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
                self.refGenome=value
            elif id=="dbSNP":
                self.dbsnp=value
            elif id=="hapmap":
                self.hapmap=value
            elif id=="omni":
                self.omni=value
            elif id=="esp":
                self.esp=value
            elif id == "aluRegions":
                self.aluRegions=value
            elif id == "gtfFile":
                self.gtfFile=value
            elif id == "output":
                self.output=value
            elif id == "sourceDir":
                self.sourceDir=value
            elif id == "maxDiff":
                self.maxDiff=str(value)
            elif id == "seedDiff":
                self.seedDiff=str(value)
            elif id == "paired":
                #Parameters.paired=float(value)
                if str(value).lower() in ("yes", "y", "true",  "t", "1"): self.paired = True
                if str(value).lower() in ("no",  "n", "false", "f", "0", "0.0", "", "none", "[]", "{}"): self.paired=False
            elif id == "standCall":
                self.standCall=str(value)
            elif id == "standEmit":
                self.standEmit=str(value)    
            elif id == "edgeDistance":
                self.edgeDistance=str(value)
            elif id == "threads":
                self.threads=str(value)
            elif id == "keepTemp":
                #Parameters.paired=float(value)
                if str(value).lower() in ("yes", "y", "true",  "t", "1"): self.keepTemp = True
                if str(value).lower() in ("no",  "n", "false", "f", "0", "0.0", "", "none", "[]", "{}"): self.keepTemp=False
            elif id == "overwrite":
                #Parameters.paired=float(value)
                if str(value).lower() in ("yes", "y", "true",  "t", "1"): self.overwrite = True
                if str(value).lower() in ("no",  "n", "false", "f", "0", "0.0", "", "none", "[]", "{}"): self.overwrite=False
                

           
class Helper():
    '''
    Helpfunctions
    '''

    
    prefix = "*** "
    praefix = " ***"
    
    #dummy element is added to the array to avoid 0/1 problem from the Tab array and these arrays
    #otherwise i had to add -1 every time i want to access the following arrays
    runningThreads=["dummy"]

    
    
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
    

    @staticmethod
    def getTime():
        '''
        return current time
        '''
        curr_time = datetime.now()
        #return "["+curr_time.strftime("%c")+"]"
        return curr_time
    
    
    @staticmethod
    def convertPhred64toPhred33(fastqFile,outFile,logFile,textField):
        """
        converts the inputFile to phred33 Quality and writes it into the ourdir
        """
        startTime=Helper.getTime()
        Helper.info("[" + startTime.strftime("%c") + "] * * * convert Quality encoding: " + fastqFile[fastqFile.rfind("/")+1:]   + " * * *",logFile,textField)
        
        
        if os.path.exists(outFile):
            Helper.info("* * * [Skipping] Result File already exists * * *",logFile,textField)
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
    
    @staticmethod
    def isPhred33Encoding(inFastqFile,lines,logFile, runNumber):
        """
        check in the first lines if the quality encoding is phred33
        """
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
    
    @staticmethod
    def proceedCommand(description,cmd,infile,outfile,rnaEdit):
        '''
        run a specific NGS-processing-step on the system
        '''
        logFile=rnaEdit.logFile
        textField=rnaEdit.textField
        overwrite=rnaEdit.params.overwrite
        
        startTime=Helper.getTime()
        Helper.info("[" + startTime.strftime("%c") + "] * * * " + description + " * * *",logFile,textField)
        
        
        #check if infile exists
        if not os.path.isfile(infile):
            Helper.error(infile + "does not exist, Error in previous Step",logFile,textField)
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
                rnaEdit.runningCommand = subprocess.Popen(cmd, stdout=resultFile, stderr=logFile)
                retcode = rnaEdit.runningCommand.wait()
                """while retcode==None:
                    #print "check if process is still running"
                    sleep(10)
                    retcode=Helper.runningCommand[runNumber].wait()
                """
                #print retcode
                
                #del Helper.runningCommand[runNumber]
                rnaEdit.runningCommand=False
                if retcode != 0:
                    if retcode == -9:
                        Helper.error(description+ " canceled by User!!!",logFile,textField)
                    else:
                        Helper.error(description+ " failed!!!",logFile,textField)
                    
                    if resultFile!=None:
                        os.remove(resultFile.name)
                    #exit(1)
            except OSError, o:
                if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
                    Helper.error(cmd[0] + " Command not found on this system",logFile,textField)
                    if resultFile!=None:
                        os.remove(resultFile.name)
                    #exit(1)
                else:
                    Helper.error(cmd[0] + o.strerror,logFile,textField)
                    if resultFile!=None:
                        os.remove(resultFile.name)
                    #exit(1)
            Helper.printTimeDiff(startTime, logFile, textField)
        else:
            print "\t [SKIP] File already exist",logFile,textField

    @staticmethod
    def getPositionDictFromVcfFile(vcfFile,runNumber):
        """
        return a dictionary whith chromosome as keys and a set of variants as values
        variantDict={chromosome:(variantPos1,variantPos2,....)}
        """
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
    def countOccurrences(inFile,column=0,logFile=None,textField=0):
        '''
        Counts how often a value appears in the given column
        :param file: 
        :param column: hold the data wich should be counted
        '''
        if type(column)!=int:
            column=int(column)
        if type(inFile) == str:
            try:
                inFile=open(inFile)
            except IOError:
                Helper.warning("Could not open %s to write Variant" % file ,logFile,textField)
        if type(inFile) != file:   
            raise AttributeError("Invalid file type in 'countOccurrences' (need string or file, %s found)" % type(inFile))
        
        keySet=()
        countDict={}
        
        for line in inFile:
            if line.startswith("#"):
                continue
            value=line.split()[column]
            if value in keySet:
                countDict[value]+=1
            else:
                keySet+=(value,)
                countDict[value]=1
        
        return countDict
    
    @staticmethod
    def createBarplot(valueMatrix,fileName,barNamesTuple,legendTuple,width=0.25,title="",yLim=None, barText=True, yText=""):
        '''
        
        :param valueMatrix: [[ValuesBar1][ValuesBar2][ValuesBar3]]
        :param fileName: 
        :param barNamesTuple: Name of Groups (has to be equal to Number of Values per Group (len(ValueArray[0])) )
        :param legendTuple: Name of each bar (equal to len(valueArray))
        :param width:
        :param title:
        '''
        valueLen=len(valueMatrix[0])
        
        
        
        for values in valueMatrix:
            assert valueLen==len(values), "ValueMatrix has to be Symmetric"
         
        assert len(legendTuple) == len(valueMatrix), "legendTuple has to have the same length as ValueArray "
        assert len(barNamesTuple) == len(valueMatrix[0]), "barNamesTuple has to have the same length as ValueArray[0] "
        
        ind = np.arange(len(valueMatrix[0]))  # the x locations for the groups
        fig, ax = plt.subplots()
        plt.subplots_adjust(bottom=0.24)
        ax.set_title(title)
        ax.set_ylabel(yText)
        ax.set_xticks(ind+width)
        ax.set_xticklabels( barNamesTuple, rotation='vertical' )
        ax.set_xlim(min(ind)-0.1,max(ind)+(width*len(valueMatrix))+0.1)
        if yLim!=None:
            ax.set_ylim(0,yLim)
        
        color=['b','g','r','y']*len(valueMatrix)
        
        
        rects=()
        
        def autolabel(rects):# attach some text labels
            for rect in rects:
                height = rect.get_height()
                ax.text(rect.get_x()+rect.get_width()/2., height/2  , '%1.0f'%float(height), ha='center', va='bottom', fontsize=8, rotation='vertical', )
        
        i=0
        for values,c in zip(valueMatrix,color):
            rect=ax.bar((ind+width*i), values, width, color=c, )
            rects+=(rect,)
            
            if barText==True:
                autolabel(rect)
            i+=1
        ax.legend( rects, legendTuple )
        
        fig.savefig(fileName)
        
    @staticmethod    
    def createDiagramms(output,logFile=None,textField=0):
        '''
        writes all the diagrams wich aree then showd in the resultTab
        :param output: output variable of Params.putput
        '''
        outdir = output[0:output.rfind("/")+1]
        sampleName=output[output.rfind("/")+1:]
        #print outdir, sampleName
        #################################################
        ####               Basecount Plot            ####
        #################################################
        counts1=Helper.getMMBaseCounts(output+".alu.vcf")
        counts2=Helper.getMMBaseCounts(output+".nonAlu.vcf")
        
        fileName=outdir+"html/"+sampleName+"_baseCounts.png"
        valueMatrix=[counts1.values(),counts2.values()]
        Helper.createBarplot(valueMatrix, fileName, counts1.keys(), ("Alu","non-Alu"),width=0.4,title="Variants per Base",yText="Number")
        
    
        #################################################
        ####       Editing per Position Plot         ####
        #################################################
        fileName=outdir+"html/"+sampleName+"_EditingPositions.png"
        counts1=Helper.countOccurrences(output+".editingSites.alu.gvf", 2, logFile, textField)
        counts2=Helper.countOccurrences(output+".editingSites.nonAlu.gvf", 2, logFile, textField) 
        
        valueMatrix=[Helper.getPercentage(counts1.values()),Helper.getPercentage(counts2.values())]
        Helper.createBarplot(valueMatrix, fileName, counts1.keys(), ("Alu","non-Alu"),width=0.4,title="Editing Sites per Position",yLim=100,yText="Precentage")
        
    @staticmethod
    def printResultHtml(output,logFile=None,textField=0):
        '''
        print the HTML file wich is shown when rnaEditor is finished
        :param output: output prefix from rnaEdit object
        '''
        outdir = output[0:output.rfind("/")+1]
        sampleName=output[output.rfind("/")+1:]
        htmlOutPrefix=outdir+"html/"+sampleName
        
        outDict={"title":"Result Page for "+ sampleName,
                 "sampleName":sampleName,
              "baseCounts":htmlOutPrefix+"_baseCounts.png",
              "editingPositions":htmlOutPrefix+"_EditingPositions.png",
              "3UTR":htmlOutPrefix+".editedGenes(3UTR).png",
              "5UTR":htmlOutPrefix+".editedGenes(5UTR).png",
              "Exon":htmlOutPrefix+".editedGenes(Exon).png",
              "Intron":htmlOutPrefix+".editedGenes(Intron).png",
              "Total":htmlOutPrefix+".editedGenes(Total).png"}
        

        outfile=open(output+".html","w+")
        
        outfile.write("""<!DOCTYPE html PUBLIC '-//W3C//DTD XHTML 1.0 Strict//EN'
    'http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd'>
<html xmlns='http://www.w3.org/1999/xhtml' xml:lang='en' lang='en'>
<head>
    <title>sampleName Results</title>
    <meta http-equiv='content-type' content='text/html;charset=utf-8' />
    <meta name='RnaEditor' content='Results from RnaEditor' />
<style type='text/css'>
    body {background-color: #F6F6F6;}
    html, body {height: 100%;margin:0px;padding: 0px;font-family: sans-serif;}
    h1,h2,h3,h4,h5,h6{margin-top: 1em;margin-bottom: 0.6em;}
    h2{    font-size: 100%;font-weight: bold;}
    h3{    font-size: 90%;    font-weight: bold;}
    #left{width: 11em;position: absolute;}
    #toc{margin-top: 40px;margin-left: 5px;}
    #toc ul{padding-left: 10px;list-style-type:none;list-style-image:none;font-size: 90%;line-height: 1.9;}
    #page, #header{z-index: 0;position: relative;padding: 0.1em 1.5em;margin-left: 11em;margin-right: 1em;border-width: 1px 1px 1px 1px;border-style: solid;border-color: #A7D7F9;}
    #page{border-width: 0px 1px 1px 1px;background: #fff;}
    #header{background-position: left bottom;background-color: #fff;}
    imgl {height: 70%;width: 90%;max-width: 850px;max-height: 550px;}
    figure {padding: 5px;float: left;border: 1px solid #cccccc;border-radius: 5px;
}
    figure img {border-radius: 3px 3px 0 0;}
    figure figcaption {font-size: 70%; padding: 2px 4px 2px 4px;background-color: #636363;color: #cccccc;font-style: italic;border-radius: 0 0 3px 3px;text-align:center}
</style>
</head>""")
        
        htmlString="""
<body>
    <div id='left'>
        <div id='logo' class='logo'>
            <img src='../ui/icons/rnaEditor_512x512.png' height=140 width=140>
        </div>
    
            <!-- Table of Content-->
        <div id='toc' class='toc'>
        <div id='toctitle'><h2>Content</h2></div>
        <ul>
            <li class='toclevel-1 tocsection-1'><a href='#basicStats'><span class='tocnumber'>1</span> <span class='toctext'>Basic Statistic</span></a></li>
            <li class='toclevel-1 tocsection-2'><a href='#NucleotideChanges'><span class='tocnumber'>2</span> <span class='toctext'>Nucleotide Changes</span></a></li>
            <li class='toclevel-1 tocsection-3'><a href='#EditingPerPosition'><span class='tocnumber'>3</span> <span class='toctext'>Editing Sites per Position</span></a></li>
            <li class='toclevel-1 tocsection-4'><a href='#EditedGenes'><span class='tocnumber'>4</span> <span class='toctext'>Highly Edited Genes</span></a></li>
            <li class='toclevel-1 tocsection-5'><a href='#EditedIslands'><span class='tocnumber'>5</span> <span class='toctext'>Editing Islands</span></a></li>
        </ul>
        </div>
    </div>    
    
    <div id='header'>
        <h1>Results for %(sampleName)s</h1>
    </div>
    
    
    <div id='page' class='page'>

        <!-- Content -->
        <div id='content' class='content'>
            
            <h2><span class='mw-headline' id='basicStats'>Basic Statistic</span></h2>
                <p>
                    Basic statistic:
                </p>
            
            <h2><span class='mw-headline' id='NucleotideChanges'>Nucleotide Changes</span></h2>
                <p>Nucleotide changes after all the filters have been applied. 
                   A high amount of A->G and T-C missmatches is indicative for a high editing ratio.</p>
                <img src='%(baseCounts)s'  alt='Base Counts' >
                
                
            <h2><span class='mw-headline' id='EditingPerPosition'>Editing Sites per Position</span></h2>
                <p>Positions where editing takes place.</p>
                <figure>
                    <img src='%(editingPositions)s'  alt='Editing Positions' >
                    <figcaption>Number of editing sites in different gene regions</figcaption>
                </figure>
            <h2><span class='mw-headline' id='EditedGenes'>Highly Edited Genes</span></h2>
                    <p>This paragraph shows highly edited genes for each segment of the genes.</p>
                    <h3>3' UTR</h3>
                        <img src='%(3UTR)s' alt='highly edited genes in 3'UTR'>
                    <h3>5' UTR</h3>
                        <img src='%(5UTR)s' alt='highly edited genes in 5'UTR'>
                    <h3>Exons</h3>
                        <img src='%(Exon)s' alt='highly edited genes in exons'>
                    <h3>Introns</h3>
                        <img src='%(Intron)s' alt='highly edited genes in introns'>
                    <h3>Total</h3>
                        <img src='%(Total)s' alt='highly edited genes (total)'>
                
                <ul>
                    <li><a rel='nofollow' class='external text' href='http://www.academia-net.de/alias/Profil/1033482'>Stefanie Dimmeler</a> in der Datenbank renommierter Wissenschaftlerinnen <a href='https://de.wikipedia.org/wiki/AcademiaNet' title='AcademiaNet'>AcademiaNet</a></li>
                </ul>
                
            <h2><span class='mw-headline' id='EditedIslands'>Editing Islands</span></h2>
                <p>Following Regions are had editing Islands</p>
        </div>
    </div>
</body>
</html>
        
"""%outDict
        outfile.write(htmlString)
        

    @staticmethod
    def getPercentage(list):
        '''
        returns the percentage as a list
        :param list:
        '''
        array=[]
        summe=float(sum(list))
        for value in list:
            array.append(round((float(value)/summe)*100.0,2))
        #print array    
        return array   

    @staticmethod
    def getMMBaseCounts(vcfFile):
        '''
        Count the number of base mismatches and return the Dictionary with the numbers
        :param vcf File:
        '''
        if type(vcfFile) == str:
            if os.path.getsize(vcfFile) == 0: #getsize raises OSError if file is not existing
                raise IOError("%s File is empty" % vcfFile)
            vcfFile = open(vcfFile,"r")
        elif type(vcfFile) != file:
            raise TypeError("Invalid type in 'parseVcfFile' (need string or file, %s found)" % type(vcfFile)) 
        
        mmBaseCounts=OrderedDict([("A->C",0),("A->G",0),("A->T",0),("C->A",0),("C->G",0),("C->T",0),("G->A",0),("G->C",0),("G->T",0),("T->A",0),("T->C",0),("T->G",0)])
        
        for line in vcfFile:
            if line.startswith('#'): continue
            line=line.split()
            if line[3]=="A" and line[4]=="C": mmBaseCounts["A->C"]+=1
            elif line[3]=="A" and line[4]=="G": mmBaseCounts["A->G"]+=1
            elif line[3]=="A" and line[4]=="T": mmBaseCounts["A->T"]+=1
            elif line[3]=="C" and line[4]=="A": mmBaseCounts["C->A"]+=1
            elif line[3]=="C" and line[4]=="G": mmBaseCounts["C->G"]+=1
            elif line[3]=="C" and line[4]=="T": mmBaseCounts["C->T"]+=1
            elif line[3]=="G" and line[4]=="A": mmBaseCounts["G->A"]+=1
            elif line[3]=="G" and line[4]=="C": mmBaseCounts["G->C"]+=1
            elif line[3]=="G" and line[4]=="T": mmBaseCounts["G->T"]+=1
            elif line[3]=="T" and line[4]=="A": mmBaseCounts["T->A"]+=1
            elif line[3]=="T" and line[4]=="C": mmBaseCounts["T->C"]+=1
            elif line[3]=="T" and line[4]=="G": mmBaseCounts["T->G"]+=1
            
        return mmBaseCounts
        
    
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
        print(traceback.format_exc())
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