'''
Created on May 22, 2013

@author: david
'''

from datetime import datetime
import argparse, sys, os, subprocess, errno
from collections import defaultdict, OrderedDict
import traceback
import ui
from numpy import arange
from PyQt5 import QtCore
from PyQt5.QtCore import pyqtSignal, QObject
#from matplotlib.backends import qt_compat
#from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.pyplot import subplots_adjust, subplots
from shutil import copyfile
from io import IOBase

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
        self.intronDistance = str(inputTab.intronDistanceSpinBox.value())
        self.minPts = str(inputTab.minPtsSpinBox.value())
        self.eps = str(inputTab.epsSpinBox.value())
        self.paired = inputTab.pairedCheckBox.isChecked()
        self.overwrite = inputTab.overwriteCheckBox.isChecked()
        self.keepTemp = inputTab.keepTempCheckBox.isChecked()
                  
    def readDefaultsFromFile(self,file="configuration.txt"):
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
            elif id == "intronDistance":
                self.intronDistance=str(value)
            elif id== "minPts":
                self.minPts=str(value)
            elif id== "eps":
                self.eps=str(value)
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
    colors = ["red","green","blue","olive","grey"]
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
        :rtype: String to converted FastQ file
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
                    if phredQual<0: #correct unvalid values
                        phredQual=0
                    elif phredQual>40:
                        phredQual=40
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
            
        Helper.error("%s has less than %i Sxequences. \n These are not enough reads for editing detection!!" % (fastqFile.name,lines),logFile, runNumber)
    
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
            Helper.error(infile + " does not exist, Error in previous Step",logFile,textField)
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
            except OSError as o:
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
            Helper.printTimeDiff(startTime, logFile, textField, "green")
        else:
            Helper.info("\t [SKIP] File already exist",logFile,textField, "green") 

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
            line=line.rstrip().split()
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
        if not (isinstance(inFile, IOBase)):   
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
        
        ind = arange(len(valueMatrix[0]))  # the x locations for the groups
        fig, ax = subplots()
        
        subplots_adjust(bottom=0.24)
        ax.set_title(title)
        ax.set_ylabel(yText)
        ax.set_xticks(ind+width)
        ax.set_xticklabels( barNamesTuple, rotation='vertical' )
        ax.set_xlim(min(ind)-0.1,max(ind)+(width*len(valueMatrix))*1.1)
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
    def printResultHtml(stats,logFile=None,textField=0):
        '''
        print the HTML file wich is shown when rnaEditor is finished
        :param output: output prefix from rnaEdit object
        '''
        
        htmlOutPrefix="html/"+stats.sampleName
        
        #copy rnaEditor logo to htmlOutPrefix
        fileDir = os.path.dirname(os.path.realpath(__file__))
        copyfile(fileDir + '/ui/icons/rnaEditor_512x512.png',stats.outdir+"html/rnaEditor_512x512.png")
        
        outDict={"title":"Result Page for "+ stats.sampleName,
                 "sampleName":stats.sampleName,
                 "icon":stats.outdir+"html/rnaEditor_512x512.png",
              "baseCounts":htmlOutPrefix+"_baseCounts.png",
              "editingPositions":htmlOutPrefix+"_EditingPositions.png",
              "3UTR":htmlOutPrefix+".editedGenes(3UTR).png",
              "5UTR":htmlOutPrefix+".editedGenes(5UTR).png",
              "Exon":htmlOutPrefix+".editedGenes(Exon).png",
              "Intron":htmlOutPrefix+".editedGenes(Intron).png",
              "Total":htmlOutPrefix+".editedGenes(Total).png",
              "currentTime":Helper.getTime().strftime("%d.%m.%Y %H:%M"),
              "totalAluNumber": str(stats.totalAluNumber),
              "totalNonAluNumber" : str(stats.totalNonAluNumber),
              "totalNumber": str(stats.totalNumber),
              "percentageEditing": str(stats.percentageEditing)+"%",
              "baseCoutHtmlTable": str(stats.baseCountHTMLTable),
              "editingPositionHTMLTable": str(stats.editingPositionHTMLTable),
              "utr3HtmlTable":  str(stats.utr3HtmlTable),
              "utr5HtmlTable": str(stats.utr5HtmlTable),
              "exonHtmlTable": str(stats.exonHtmlTable),
              "intronHtmlTable": str(stats.intronHtmlTable),
              "totalHtmlTable": str(stats.totalHtmlTable)
              }
        

        outfile=open(stats.output+".html","w+")
        
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
    h3{margin-left:1em;font-size: 90%;    font-weight: bold;}
    .contentBlock{margin-top: 1em;border-width: 1px 1px 3px 1px;border-style: solid;border-color: #A7D7F9; max-width: 950px}
    #left{width: 11em;position: absolute;}
    #toc{margin-top: 40px;margin-left: 5px;}
    #toc ul{padding-left: 10px;list-style-type:none;list-style-image:none;font-size: 90%;line-height: 1.9;}
    #page, #header{z-index: 0;position: relative;padding: 0.1em 1.5em;margin-left: 11em;margin-right: 1em;border-width: 1px 1px 1px 1px;border-style: solid;border-color: #A7D7F9;}
    #page{border-width: 0px 1px 1px 1px;background: #fff;}
    #header{background-position: left bottom;background-color: #fff;}
    img {}
    table {font-size: 80%;background-color: #DDD;width:70%;max-width: 850px;margin-left:auto;margin-right:auto; text-align: center;}
    th {text-align: center;background-color: #0000FF;color: #FFF;padding: 0.4em;}      
    td {text-align: center;background-color: #FFF;color: #000;padding: 0.4em;}
    .geneTable{width: 500px;}

    figure {padding: 5px;border: 1px solid #cccccc;border-radius: 5px;max-width: 900px}
    figure img {border-radius: 3px 3px 0 0;height: 70%;width: 90%;max-width: 850px;max-height: 550px;}
    figure figcaption {font-size: 70%; padding: 2px 4px 2px 4px;background-color: #636363;color: #cccccc;font-style: italic;border-radius: 0 0 3px 3px;text-align:center}
</style>
</head>""")
        
        htmlString="""
<body>
    <div id='left'>
        <div id='logo' class='logo'>
            <img src='%(icon)s' height=140 width=140>
        </div>
    
            <!-- Table of Content-->
        <div id='toc' class='toc'>
        <div id='toctitle'><h2>Content</h2></div>
        <ul>
            <li class='toclevel-1 tocsection-1'><a href='#basicStats'><span class='tocnumber'>1</span> <span class='toctext'>Basic Statistic</span></a></li>
            <li class='toclevel-1 tocsection-2'><a href='#NucleotideChanges'><span class='tocnumber'>2</span> <span class='toctext'>Nucleotide Changes</span></a></li>
            <li class='toclevel-1 tocsection-3'><a href='#EditingPerPosition'><span class='tocnumber'>3</span> <span class='toctext'>Editing Sites per Position</span></a></li>
            <ul>
                <li class='toclevel-2 tocsection-3.1'><a href='#3utr'><span class='tocnumber'>3.1</span> <span class='toctext'>3'UTR</span></a></li>
                <li class='toclevel-2 tocsection-3.2'><a href='#5utr'><span class='tocnumber'>3.2</span> <span class='toctext'>5'UTR</span></a></li>
                <li class='toclevel-2 tocsection-3.3'><a href='#exon'><span class='tocnumber'>3.3</span> <span class='toctext'>Exons</span></a></li>
                <li class='toclevel-2 tocsection-3.4'><a href='#intron'><span class='tocnumber'>3.4</span> <span class='toctext'>Introns</span></a></li>
                <li class='toclevel-2 tocsection-3.4'><a href='#total'><span class='tocnumber'>3.5</span> <span class='toctext'>Total</span></a></li>
            </ul>
            <li class='toclevel-1 tocsection-4'><a href='#EditedGenes'><span class='tocnumber'>4</span> <span class='toctext'>Highly Edited Genes</span></a></li>
            <li class='toclevel-1 tocsection-5'><a href='#EditedIslands'><span class='tocnumber'>5</span> <span class='toctext'>Editing Islands</span></a></li>
        </ul>
        </div>
    </div>    
    
    <div id='header'>
        <h1 >Results for %(sampleName)s</h1>
        <p style="text-align: right;">%(currentTime)s</p>
    </div>
    
    
    <div id='page' class='page'>

        <!-- Content -->
        <div id='content' class='content'>
            <div class='contentBlock'>
            <h2><span class='mw-headline' id='basicStats'>Basic Statistic</span></h2>
                <p>
                    <table>
                        <tbody>
                            <tr>
                                <th>Measure</th>
                                <th>Value</th>        
                            </tr>
                            <tr>
                                <td>Total Number of editing sites:</td>
                                <td>%(totalNumber)s</td>
                            </tr>
                            <tr>
                                <td>in Alu Regions:</td>
                                <td>%(totalAluNumber)s</td>
                            </tr>
                            <tr>
                                <td>in non-Alu Regions:</td>
                                <td>%(totalNonAluNumber)s</td>
                            </tr>
                            <tr>
                                <td>Percentage of edited Genes:</td>
                                <td>%(percentageEditing)s</td>
                            </tr>
                        </tbody>
                    </table>
                </p>
            </div>
            <div class='contentBlock'>
            <h2><span class='mw-headline' id='NucleotideChanges'>Nucleotide Changes</span></h2>
                <p>Nucleotide changes after all the filters have been applied. 
                   A high amount of A->G and T->C missmatches is indicative for a high editing rate.</p>
                <figure>
                    <img src='%(baseCounts)s'  alt='Base Counts' >
                    <figcaption>Number of mismatsch types</figcaption>
                </figure>
                %(baseCoutHtmlTable)s
            </div>
            <div class='contentBlock'>
            <h2><span class='mw-headline' id='EditingPerPosition'>Editing Sites per Position</span></h2>
                <p>Positions where editing takes place.</p>
                <figure>
                    <img src='%(editingPositions)s'  alt='Editing Positions' >
                    <figcaption>Number of editing sites in different gene regions</figcaption>
                </figure>
                %(editingPositionHTMLTable)s
            </div>
            <div class='contentBlock'>
            <h2><span class='mw-headline' id='EditedGenes'>Highly Edited Genes</span></h2>
                    <p>This paragraph shows highly edited genes for each segment of the genes.</p>
                    <h3 id='3utr'>3' UTR</h3>
                        <figure>
                            <img src='%(3UTR)s' alt='highly edited genes in 3'UTR'>
                            <figcaption>Highly edited Genes in 3'UTR Regions</figcaption>
                        </figure>
                        %(utr3HtmlTable)s
                    <h3 id='5utr'>5' UTR</h3>
                        <figure>
                            <img src='%(5UTR)s' alt='highly edited genes in 5'UTR'>
                            <figcaption>Highly edited Genes in 5'UTR Regions</figcaption>
                        </figure>
                        %(utr5HtmlTable)s
                    <h3 id='exon'>Exons</h3>
                        <figure>
                            <img src='%(Exon)s' alt='highly edited genes in exons'>
                            <figcaption>Highly edited Genes in Exon Regions</figcaption>
                        </figure>
                        %(exonHtmlTable)s
                    <h3 id='intron'>Introns</h3>
                        <figure>
                            <img src='%(Intron)s' alt='highly edited genes in introns'>
                            <figcaption>Highly edited Genes in Intron Regions</figcaption>
                        </figure>
                        %(intronHtmlTable)s
                    <h3 id='total'>Total</h3>
                        <figure>
                            <img src='%(Total)s' alt='highly edited genes (total)'>
                            <figcaption>Highly edited Genes in 5'UTR Regions</figcaption>
                        </figure>
                        %(totalHtmlTable)s
            </div>   
            <div class='contentBlock'>
            <h2><span class='mw-headline' id='EditedIslands'>Editing Islands</span></h2>
                <p> Detected Editing Islands</p>
            </div>
                
            
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
        elif not (isinstance(vcfFile, IOBase)):
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
    def printTimeDiff(startTime,logFile=None,textField=0,color="green"):
        duration = Helper.getTime() - startTime
        if textField!=0:
            #currentAssay = Helper.runningAssays[textField] 
            if color in Helper.colors:
                textField.append(str("<font color=\""+color+"\">{}</font>").format("[DONE] Duration [" + str(duration) + "]"  + Helper.praefix + "\n"))
            else:
                textField.append(Helper.prefix + "[DONE] Duration [" + str(duration) + "]"  + Helper.praefix + "\n")
        if logFile!=None:
            logFile.write("\t" + Helper.prefix + "[DONE] Duration [" + str(duration) + "]"  + Helper.praefix + "\n")
        
        sys.stderr.write("\t" + Helper.prefix + "[DONE] Duration [" + str(duration) + "]"  + Helper.praefix + "\n")

    @staticmethod
    def newline(quantity=1,
                 logFile=None,
                 textField=0):
        if textField!=0:
            # currentAssay = Helper.runningAssays[runNumber]
            textField.append("\n"*quantity)
        if logFile != None:
            logFile.write("\n"*quantity)
            logFile.flush()
        sys.stderr.write("\n"*quantity)
    @staticmethod
    def info (message,logFile=None,textField=0,color="olive"):
        if textField!=0:
            if color in Helper.colors:
                textField.append(str("<font color=\""+color+"\">{}</font>").format(Helper.prefix + "STATUS:    "  + message + Helper.praefix))
            else:
                textField.append(Helper.prefix + "INFO:    "  + message + Helper.praefix)
        if logFile!=None:
            logFile.write(Helper.prefix + "INFO:    "  + message + Helper.praefix + "\n")
            logFile.flush()
        sys.stderr.write(Helper.prefix + "INFO:    "  + message + Helper.praefix + "\n")
    @staticmethod
    def warning (message,logFile=None,textField=0):
        if textField!=0:
            textField.append("\n\n" + Helper.prefix + "WARNING:    " + message + Helper.praefix + "\n\n")
        if logFile!=None:
            logFile.write(Helper.prefix + "WARNING:    "  + message + Helper.praefix + "\n")
            logFile.flush()
        sys.stderr.write("\n\n" + Helper.prefix + "WARNING:    " + message + Helper.praefix + "\n\n")
    @staticmethod
    def error (message,logFile=None,textField=0,color="red"):
        if textField!=0:
            if color in Helper.colors:
                textField.append(str("<font color=\""+color+"\">{}</font>").format(Helper.prefix + "STATUS:    "  + message + Helper.praefix))
            else:
                textField.append(Helper.prefix + "ERROR:    "  + message + Helper.praefix)
        if logFile!=None:
            logFile.write(Helper.prefix + "ERROR:    "  + message + Helper.praefix + "\n")
            logFile.flush()
        print(traceback.format_exc())
        #sys.stderr.write("\n\n" + Helper.prefix + "ERROR:    " + message + Helper.praefix + "\n\n")
        raise Exception("\n\n" + Helper.prefix + "ERROR:    " + message + Helper.praefix + "\n\n")
    @staticmethod
    def debug (message,logFile=None,textField=0):
        if textField!=0:
            textField.append(Helper.prefix + "DEBUG:    "  + message + Helper.praefix)
        if logFile!=None:
            logFile.write(Helper.prefix + "DEBUG:    "  + message + Helper.praefix + "\n")
            logFile.flush()
        sys.stderr.write(Helper.prefix + message + Helper.praefix + "\n")
    @staticmethod
    def status(message,logFile=None,textField=0,color=None,bold=False):
        if textField!=0:
            #currentAssay = Helper.runningAssays[runNumber] 
            if color in Helper.colors:
                if bold==True:
                    textField.append(str("<font color=\""+color+"\"><b>{}</b></font>").format(Helper.prefix + "STATUS:    "  + message + Helper.praefix))
                else:
                    textField.append(str("<font color=\""+color+"\">{}</font>").format(Helper.prefix + "STATUS:    "  + message + Helper.praefix))
            else:
                textField.append(Helper.prefix + "STATUS:    "  + message + Helper.praefix)
        if logFile!=None:
            logFile.write(Helper.prefix + "STATUS:    "  + message + Helper.praefix + "\n")
            logFile.flush()
        sys.stdout.write("\r" + Helper.prefix + "STATUS:    "  + message + Helper.praefix + "\n")
        sys.stdout.flush()
        
        
        
        
class Communicate(QObject):

    taskDone = pyqtSignal(str)