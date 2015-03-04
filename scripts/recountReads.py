'''
Created on 02.03.2015
This script recalculates the number of read per base after they were joined.
If there is an 

@author: david
'''

import pysam

for i in range(1,1000):
    samfile = pysam.AlignmentFile("/media/Storage/bio-data/David/Kostas/scrambleN/scrambleN.realigned.marked.recalibrated.bam", "rb")
    s=samfile.fetch('1', 150703594,150703594)
    list=[]
    if i % 10 ==0:
        print i
    for read in s:
         #print read
         pass
         list.append(read)