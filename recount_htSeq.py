'''
Created on 02.03.2015

@author: david
'''

import HTSeq

samfile = HTSeq.BAM_Reader("/media/Storage/bio-data/David/Kostas/scrambleN/scrambleN.realigned.marked.recalibrated.bam")

for i in range(1,1000):

    s=samfile.fetch( region = "1:150703594-150703594" ) #fetching reads in a region
    list=[]
    if i % 10 ==0:
        print i
    for read in s:
         #print read
         pass
         list.append(read)