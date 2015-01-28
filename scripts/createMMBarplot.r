#!/usr/bin/R

Args <- commandArgs(trailingOnly=TRUE)
inFile <- Args[1]
plot <- paste(Args[1],".mmCount.pdf",sep="")
outFile <- paste(Args[1],".mmCount.txt",sep="")


table_ScrambleN <- read.table(inFile,sep="\t", quote="",dec=".", skip=0)
#table_ScrambleN <- read.table("/home//uchida/workspace/rnaEditor/testdata/adar1.nonAlu.noSpliceSites.noHomo.vcf",sep="\t", quote="",dec=".", skip=0)

class(table_ScrambleN$V4)
class(table_ScrambleN$V5)

vector_ScrambleN <- c("A->C"=0,"A->G"=0,"A->T"=0,"C->A"=0,"C->G"=0,"C->T"=0,"G->A"=0,"G->C"=0,"G->T"=0,"T->A"=0,"T->C"=0,"T->G"=0)

total_ScrambleN<-nrow(table_ScrambleN)

vector_ScrambleN["A->C"] = nrow(subset(table_ScrambleN,V4 == "A" & V5 == "C"))
vector_ScrambleN["A->G"] = nrow(subset(table_ScrambleN,V4 == "A" & V5 == "G"))
vector_ScrambleN["A->T"] = nrow(subset(table_ScrambleN,V4 == "A" & V5 == "T"))
vector_ScrambleN["C->A"] = nrow(subset(table_ScrambleN,V4 == "C" & V5 == "A"))
vector_ScrambleN["C->G"] = nrow(subset(table_ScrambleN,V4 == "C" & V5 == "G"))
vector_ScrambleN["C->T"] = nrow(subset(table_ScrambleN,V4 == "C" & V5 == "T"))
vector_ScrambleN["G->A"] = nrow(subset(table_ScrambleN,V4 == "G" & V5 == "A"))
vector_ScrambleN["G->C"] = nrow(subset(table_ScrambleN,V4 == "G" & V5 == "C"))
vector_ScrambleN["G->T"] = nrow(subset(table_ScrambleN,V4 == "G" & V5 == "T"))
vector_ScrambleN["T->A"] = nrow(subset(table_ScrambleN,V4 == "T" & V5 == "A"))
vector_ScrambleN["T->C"] = nrow(subset(table_ScrambleN,V4 == "T" & V5 == "C"))
vector_ScrambleN["T->G"] = nrow(subset(table_ScrambleN,V4 == "T" & V5 == "G"))

pdf(file=plot)
#counts<-rbind(vector_Adar1,vector_Adar2,vector_ScrambleG,vector_ScrambleN)
colors<-c("blue","green","yellow","red")
#rownames(counts)<-c("Adar1","Adar2","ScrambleG","ScrambleN")
barplot(vector_ScrambleN,xlab="Mismatch type",ylab="Number of mismatches",main="",col=colors,beside=TRUE,names.arg=rownames(vector_ScrambleN))
#legend("top", legend=rownames(vector_ScrambleN), fill = colors)
#write.table(counts,file="/media/media/bio-data/Konstantinos-Stellos/mRNA_Kat_121005/fastq/pooled/mapped/numberOfMissmatchesInAlu.txt")
dev.off()