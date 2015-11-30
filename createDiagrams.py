'''
Created on Nov 30, 2015

@author: david
'''
import argparse 
from Helper import Helper
import os




def createDiagramms(output,logFile=None,textField=0):
        '''
        writes all the diagrams wich aree then showd in the resultTab
        :param output: output variable of Params.output
        '''
        Helper.status("Creating Diagrams for %s" % output, logFile, textField)
        
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
        
parser = argparse.ArgumentParser(description='create Diagrams for the output of RNA Editor results')
parser.add_argument('-o', '--output', metavar='str', type=str, help="output directory and Samplename (Example: '/home/rnaEditor/sample/sample')", required=True)
args = parser.parse_args()


if __name__ == '__main__':
    createDiagramms(args.output)