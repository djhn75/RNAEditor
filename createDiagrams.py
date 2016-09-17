'''
Created on Nov 30, 2015

@author: david
'''
import argparse 
from Helper import Helper
import os
from collections import OrderedDict


class Stats():
    def __init__(self, output):
        self.output = output
        


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
    totalGenes=0
    for line in sumFile:
        if line.startswith("#"): continue #skip comments
        line = line.rstrip().split()
        totalGenes+=1
        if int(line[6])<1: continue #skip unedited genes
        try:
            v=map(int,line[2:7])
        except ValueError:
            v=line[2:7]
        dict[line[0]]=[line[1]]+v
        
    return dict,totalGenes
        
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
        elif value  ==2:
            barName="5'-UTR"
        elif value==3:
            barName="Exonic"
        elif value==4:
            barName="Intronic"
        elif value==5:
            barName="Total"
        
        yLim=max(max(i) for i in valueMatrix)+1
        Helper.createBarplot(valueMatrix, fileName, barNameTuple, [barName], width=0.35, title="Highly Edited Genes",yLim=yLim,barText=False,yText="Editing Counts")
        
        
        file= open(fileName.replace("png","txt"),"w")
        file.write("\t".join(["Gene_Symbol","Number_of_editing_sites"])+"\n")
        htmlStr="<table class='geneTable'><tr><th>GeneName</th><th>Number of editing sites</th></tr>"
        for gene in counts:
            htmlStr+="<tr><td>%s</td><td>%s</td></tr>"%(counts[gene][0],counts[gene][value])
            geneName=counts[gene][0]
            numbers=str(counts[gene][value])
            file.write("\t".join([geneName,numbers]) +"\n")
        htmlStr+="</table>"
        return htmlStr

        
        

def createDiagramms(output, geneNumber=20,logFile=None,textField=0):
        '''
        writes all the diagrams wich aree then showd in the resultTab
        All the values are stored in an instance of the Class Stats
        
        :param output: output variable of Params.output
        '''
        Helper.info("Creating Diagrams for %s" % output, logFile, textField)
        
        stats = Stats(output)
        
        stats.output = output
        stats.outdir = output[0:output.rfind("/")+1]
        stats.sampleName=output[output.rfind("/")+1:]
        
        if not os.path.exists(stats.outdir+"html/"):
            os.makedirs(stats.outdir+"html/")
        
        #print outdir, sampleName
        #################################################
        ####               Basecount Plot            ####
        #################################################
        counts1=Helper.getMMBaseCounts(output+".alu.vcf")
        counts2=Helper.getMMBaseCounts(output+".noBlat.vcf") #use the var file after all the filters have been applied for nonAlu regions
        
        
        stats.totalAluNumber = counts1["A->G"] + counts1["T->C"]
        stats.totalNonAluNumber = counts2["A->G"] + counts2["T->C"]
        stats.totalNumber = stats.totalAluNumber + stats.totalNonAluNumber
        
        
        #write the baseCounts to a file to open with Excel later
        stats.baseCountHTMLTable="<table><tr><th>Missmatch Type</th><th>Alu</th><th>Non Alu</th></tr>"
        file=open(stats.outdir+"html/"+stats.sampleName+"_baseCounts.txt","w")
        file.write("\t".join(["MM","alu","nonAlu"])+"\n")
        for keyAlu,keyNonAlu in zip(counts1.keys(),counts2.keys()):
            file.write("\t".join([str(keyAlu),str(counts1[keyAlu]),str(counts2[keyAlu])])+"\n")
            stats.baseCountHTMLTable+="<tr><td>%s</td><td>%s</td><td>%s</td></tr>"%(keyAlu,str(counts1[keyAlu]),str(counts2[keyAlu]))
        file.close()
        stats.baseCountHTMLTable+="</table>"
        fileName=stats.outdir+"html/"+stats.sampleName+"_baseCounts.png"
        
        valueMatrix=[counts1.values(),counts2.values()]
        Helper.createBarplot(valueMatrix, fileName, counts1.keys(), ("Alu","non-Alu"),width=0.4,title="Variants per Base", barText=False, yText="Number",)
        
    
        #################################################
        ####       Editing per Position Plot         ####
        #################################################
        fileName=stats.outdir+"html/"+stats.sampleName+"_EditingPositions.png"
        fileNamePercentage=stats.outdir+"html/"+stats.sampleName+"_EditingPositions(Percentage).png"
        counts1=Helper.countOccurrences(output+".editingSites.alu.gvf", 2, logFile, textField)
        counts2=Helper.countOccurrences(output+".editingSites.nonAlu.gvf", 2, logFile, textField) 
        
        file=open(stats.outdir+"html/"+stats.sampleName+"_editingSites.txt","w")
        file.write("\t".join(["Position","alu","nonAlu"])+"\n")
        
        for key in counts1.keys():
            if key in counts2.keys():
                file.write("\t".join([str(key),str(counts1[key]),str(counts2[key])])+"\n")
            else:
                file.write("\t".join([str(key), str(counts1[key]), "--"]) + "\n")
        for key in counts2.keys():
            if key not in counts1.keys():
                file.write("\t".join([str(key), "--", str(counts2[key])]) + "\n")
        file.close()
        
        
        #set values to 0 if they dont exist in the opposite file
        orderList = ["3'UTR","5'UTR","coding-exon","noncoding-exon","intron","intergenic"]
        aluPositions,nonAluPositions = [],[]
        for key in orderList:
            aluPositions.append(counts1[key]) if key in counts1.keys() else aluPositions.append(0.000000001)
            nonAluPositions.append(counts2[key]) if key in counts2.keys() else nonAluPositions.append(0.00000001)
            """if key in counts1.keys():
                aluPositions.append(counts1[key])
            else:
                aluPositions.append(0)
        for key in orderList:
            nonAluPositions.append(counts2[key]) 
            """
        sumAlu,sumNonAlu = sum(aluPositions),sum(nonAluPositions)
        #aluPositions=[counts1["3'UTR"],counts1["5'UTR"],counts1["coding-exon"],counts1["noncoding-exon"],counts1["intron"],counts1["intergenic"]]
        #nonAluPositions=[counts2["3'UTR"],counts2["5'UTR"],counts2["coding-exon"],counts2["noncoding-exon"],counts2["intron"],counts2["intergenic"]]
        barNames=["3'UTR","5'UTR","coding-exon","noncoding-exon","intron","intergenic"]
        valueMatrix=[aluPositions,nonAluPositions]
        
        Helper.createBarplot(valueMatrix, fileName, barNames, ("Alu","non-Alu"),width=0.4,title="Editing Sites per Position", barText=False, yText="Total Counts")
        
        valueMatrixPercentage=[Helper.getPercentage(aluPositions),Helper.getPercentage(nonAluPositions)]
        Helper.createBarplot(valueMatrixPercentage, fileNamePercentage, barNames, ("Alu","non-Alu"),width=0.4,title="Editing Sites per Position",yLim=100,yText="Precentage")
        
        #make String for the HTML Table and write to a table for Ecxel
        file=open(stats.outdir+"html/"+stats.sampleName+"_editingSites.txt","w")
        file.write("\t".join(["Position","alu","nonAlu"])+"\n")
        stats.editingPositionHTMLTable="<table><tr><th>Editing Position</th><th>Total Alu</th><th>Alu Percentage</th><th>Total Non Alu</th><th>Non Alu Percentage</th></tr>"
        for key in orderList:
            alu = counts1[key] if key in counts1.keys() else 0
            nonAlu = counts2[key] if key in counts2.keys() else 0
            aluNumber=str(counts1[key]) if key in counts1.keys() else "0"
            nonAluNumber=str(counts2[key]) if key in counts2.keys() else "0"
            aluPercentage= str(round(float(alu)/sumAlu,3)*100)+" %"
            nonAluPercentage= str(round(float(nonAlu)/sumNonAlu,3)*100)+" %"
            file.write("\t".join([str(key),aluNumber,nonAluNumber])+"\n")
            stats.editingPositionHTMLTable+="<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>"%(key,aluNumber, aluPercentage,nonAluNumber, nonAluPercentage)
        file.close()
        stats.editingPositionHTMLTable+="</table>"
        
       
            
        
        #################################################
        ####           Edited Genes Plot            ####
        #################################################
        sumDict,totalGenes=parseSummaryFile(output+".editingSites.summary", logFile, textField)
        stats.percentageEditing = round(float(len(sumDict))/float(totalGenes)*100.0, 2)
        if "intergenic" in sumDict.keys():
            del sumDict["intergenic"] 
        fileName=stats.outdir+"html/"+stats.sampleName+".editedGenes(3UTR).png"
        stats.utr3HtmlTable = topGenesDict = topGenes(sumDict,fileName, geneNumber, 1)
           
        fileName=stats.outdir+"html/"+stats.sampleName+".editedGenes(5UTR).png"
        stats.utr5HtmlTable = topGenes(sumDict,fileName, geneNumber, 2)
            
        fileName=stats.outdir+"html/"+stats.sampleName+".editedGenes(Exon).png"
        stats.exonHtmlTable = topGenes(sumDict,fileName, geneNumber, 3)
            
        fileName=stats.outdir+"html/"+stats.sampleName+".editedGenes(Intron).png"
        stats.intronHtmlTable = topGenes(sumDict,fileName, geneNumber, 4)
            
        if "intergenic" in sumDict.keys():
            del sumDict["-"] #delete intergenics, because we only we only want to show highly edited Genes!!!
        fileName=stats.outdir+"html/"+stats.sampleName+".editedGenes(Total).png"
        stats.totalHtmlTable = topGenes(sumDict,fileName, geneNumber, 5)
        
        Helper.printResultHtml(stats, logFile, textField)
        
parser = argparse.ArgumentParser(description='create Diagrams for the output of RNA Editor results')
parser.add_argument('-o', '--output', metavar='str', type=str, help="output directory and Samplename (Example: '/home/rnaEditor/sample/sample')", required=True)
parser.add_argument('-n', '--geneNumber', metavar='int', type=int, help='Number of genes shown in the figure', default=20)
args = parser.parse_args()


if __name__ == '__main__':
    createDiagramms(args.output,args.geneNumber)