'''
Created on 02.03.2015

@author: david
'''
from Helper import Helper

for i in range(1,1000):
    line=[]
    command = ["samtools", "view", "/media/Storage/bio-data/David/Kostas/scrambleN/scrambleN.realigned.marked.recalibrated.bam", "1:150703594-150703594"] 
    samout = Helper.getCommandOutput(command).splitlines() #get the reads wich are overlapping the snp region
    if i % 10 == 0:
        print i
    for samline in samout:
        line.append(samline)

