'''
Created on 02.03.2015
This script recalculates the number of read per base after they were joined.
If there is an 

@author: david
'''

import pysam
from Helper import Helper
import argparse, os

parser = argparse.ArgumentParser(description='Merge tables.')
parser.add_argument('-f', '--files', metavar='N', type=str, nargs='+', help='the list of files')
parser.add_argument('-t', '--top', metavar='N', type=str, nargs="+", help='list of header names (space separated)')
parser.add_argument('-c', '--columns', metavar='N', type=int, nargs='+', help='columns to keep (space separated)',default=[2],)
parser.add_argument('-k', '--keys', metavar='N', nargs='+', type=int, help='columnnumber on which to join',default=[1])
parser.add_argument('-d', '--delimiter', metavar='N', type=str, help='delimiter', default="\t")
parser.add_argument('-e', '--empty', metavar='N', type=str, help='Sign for empty Values', default="--")
args = parser.parse_args()

startTime = Helper.getTime()
samfile = pysam.AlignmentFile("/media/Storage/bio-data/David/Kostas/scrambleN/scrambleN.realigned.marked.recalibrated.bam", "rb")
for i in range(1,1000):
    
    s=samfile.fetch('1', 100703594,150703999)
    list=[]
    if i % 10 ==0:
        print i
    for read in s:
         #print read
         pass
         list.append(read)
         
Helper.printTimeDiff(startTime)