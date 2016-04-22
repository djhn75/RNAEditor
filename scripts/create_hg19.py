#!/usrbin/python
#download all the necessary Files to the given Folder

#@author: david


import argparse, sys
from string import split, replace, ascii_letters
from subprocess import STDOUT
import subprocess
import urllib2
import zipfile


parser = argparse.ArgumentParser(description='Download all the files for RNAeditor')
parser.add_argument('-o', '--outDir', metavar='N', type=str, help='Directory where RNAeditor stores the Annotation Files', required=True)
args = parser.parse_args()


url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz"

outDir=args.outDir
file_name = url.split('/')[-1]

cmd=["curl",url]
subprocess.call(cmd,stdout=open(outDir+file_name,"w"))

"""
"u = urllib2.urlopen(url)
f = open(outDir+file_name, 'wb')
meta = u.info()
file_size = int(meta.getheaders("Content-Length")[0])
print "Downloading: %s Bytes: %s" % (file_name, file_size)

file_size_dl = 0
block_sz = 8192
while True:
    buffer = u.read(block_sz)
    if not buffer:
        break

    file_size_dl += len(buffer)
    f.write(buffer)
    status = r"%10d  [%3.2f%%]" % (file_size_dl, file_size_dl * 100. / file_size)
    status = status + chr(8)*(len(status)+1)
    print status,

f.close()
"""