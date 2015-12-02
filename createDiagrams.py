'''
Created on Nov 30, 2015

@author: david
'''
import argparse 
from Helper import Helper
import os
from collections import OrderedDict



def parseSummaryFile(sumFile,logFile=None,textField=0):
    '''
    Parses a .summary file from an rnaEditor output directory and returns it as an ordered dict
    Note: unedited Genes will be skipped
    :param sumFile: .summary file of rnaEditor
    :param logFile:
    :param textField:
    :return: OrderedDict {GeneName1:[GeneId1,3'UTR,5'UTR,EXON,Intron,Total]}
    '''
    
    
    if type(sumFile)==str:
        try:
            sumFile=open(sumFile,"r")
        except IOError:
            Helper.warning("Could not open %s to write Variant" % sumFile ,logFile,textField)
    elif type(sumFile)==file:
        pass
    else:
        raise TypeError("Summary file hat to be path or file object", logFile, textField)
    
    dict=OrderedDict()
    for line in sumFile:
        if line.startswith("#"): continue #skip comments
        line = line.rstrip().split()
        if line[6]<1:continue #skip unedited genes
        try:
            v=map(int,line[2:7])
        except ValueError:
            v=line[2:7]
        dict[line[0]]=[line[1]]+v
        
    return dict
        
def topGenes(sumDict, fileName,number=20,value=5, logFile=None,textField=0):
        if number > len(sumDict):
            if len(sumDict)<1:
                Helper.warning("no edited genes found", logFile, textField)
                return
            Helper.warning("The given gene number is bigger than the number of total edited genes", logFile, textField)
            number=len(sumDict)
        if value not in (1,2,3,4,5):
            Helper.error("sumDict hast to be between 1 an 5", logFile, textField)
        
        
        counts=OrderedDict(sorted(sumDict.items(), key=lambda t: t[1][value],reverse=True)[:number])
        barNameTuple=()
        valueMatrix=[[]]
        for array in counts.values():
            valueMatrix[0].append(array[value])
        for gene in counts.keys():
            barNameTuple+=(counts[gene][0],)

        if value==1:
            barName="3'-UTR"
        elif value==2:
            barName="5'-UTR"
        elif value==3:
            barName="Exonic"
        elif value==4:
            barName="Intronic"
        elif value==5:
            barName="Total"
        
        yLim=max(max(i) for i in valueMatrix)+1
        Helper.createBarplot(valueMatrix, fileName, barNameTuple, [barName], width=0.35, title="Highly Edited Genes",yLim=yLim,barText=False,yText="Editing Counts")


def createDiagramms(output, geneNumber=20,logFile=None,textField=0):
        '''
        writes all the diagrams wich aree then showd in the resultTab
        :param output: output variable of Params.output
        '''
        Helper.info("Creating Diagrams for %s" % output, logFile, textField)
        
        outdir = output[0:output.rfind("/")+1]
        sampleName=output[output.rfind("/")+1:]
        
        if not os.path.exists(outdir+"html/"):
            os.makedirs(outdir+"html/")
        
        #print outdir, sampleName
        #################################################
        ####               Basecount Plot            ####
        #################################################
        counts1=Helper.getMMBaseCounts(output+".alu.vcf")
        counts2=Helper.getMMBaseCounts(output+".nonAlu.vcf")
        
        file=open(outdir+"html/"+sampleName+"_baseCounts.txt","w")
        file.write("MM    alu    nonAlu\n")
        for key in counts1.keys():
            file.write("\t".join([str(key),str(counts1[key]),str(counts2[key]),"\n"]))
        file.close()
        
        fileName=outdir+"html/"+sampleName+"_baseCounts.png"
        
        valueMatrix=[counts1.values(),counts2.values()]
        Helper.createBarplot(valueMatrix, fileName, counts1.keys(), ("Alu","non-Alu"),width=0.4,title="Variants per Base",yText="Number")
        
    
        #################################################
        ####       Editing per Position Plot         ####
        #################################################
        fileName=outdir+"html/"+sampleName+"_EditingPositions.png"
        fileNamePercentage=outdir+"html/"+sampleName+"_EditingPositions(Percentage).png"
        counts1=Helper.countOccurrences(output+".editingSites.alu.gvf", 2, logFile, textField)
        counts2=Helper.countOccurrences(output+".editingSites.nonAlu.gvf", 2, logFile, textField) 
        
        file=open(outdir+"html/"+sampleName+"_editingSites.txt","w")
        file.write("Position    alu    nonAlu\n")
        for key in counts1.keys():
            file.write("\t".join([str(key),str(counts1[key]),str(counts2[key]),"\n"]))
        file.close()
        
        #set values to 0 if they dont exist in the opposite file
        
        a=[counts1["3'UTR"],counts1["5'UTR"],counts1["coding-exon"],counts1["noncoding-exon"],counts1["intron"],counts1["intergenic"]]
        b=[counts2["3'UTR"],counts2["5'UTR"],counts2["coding-exon"],counts2["noncoding-exon"],counts2["intron"],counts2["intergenic"]]
        valueMatrix=[a,b]
        Helper.createBarplot(valueMatrix, fileName, counts1.keys(), ("Alu","non-Alu"),width=0.4,title="Editing Sites per Position",yText="Total Counts")
        
        valueMatrix=[Helper.getPercentage(a),Helper.getPercentage(b)]
        Helper.createBarplot(valueMatrix, fileNamePercentage, counts1.keys(), ("Alu","non-Alu"),width=0.4,title="Editing Sites per Position",yLim=100,yText="Precentage")
        
        
        #################################################
        ####           Edited Genes Plot            ####
        #################################################
        sumDict=parseSummaryFile(output+".editingSites.summary", logFile, textField)
        del sumDict["intergenic"]
        fileName=outdir+"html/"+sampleName+".editedGenes(3UTR).png"
        topGenes(sumDict,fileName, geneNumber, 1)
            
        fileName=outdir+"html/"+sampleName+".editedGenes(5UTR).png"
        topGenes(sumDict,fileName, geneNumber, 2)
            
        fileName=outdir+"html/"+sampleName+".editedGenes(Exon).png"
        topGenes(sumDict,fileName, geneNumber, 3)
            
        fileName=outdir+"html/"+sampleName+".editedGenes(Intron).png"
        topGenes(sumDict,fileName, geneNumber, 4)
            
        fileName=outdir+"html/"+sampleName+".editedGenes(Total).png"
            
        if "intergenic" in sumDict.keys():
            del sumDict["-"] #delete intergenics, because we only we only want to show highly edited Genes!!!
        topGenes(sumDict,fileName, geneNumber, 5)
        
        Helper.printResultHtml(output, logFile, textField)
        
parser = argparse.ArgumentParser(description='create Diagrams for the output of RNA Editor results')
parser.add_argument('-o', '--output', metavar='str', type=str, help="output directory and Samplename (Example: '/home/rnaEditor/sample/sample')", required=True)
parser.add_argument('-n', '--geneNumber', metavar='int', type=int, help='Number of genes shown in the figure', default=20)
args = parser.parse_args()


if __name__ == '__main__':
    createDiagramms(args.output,args.geneNumber)