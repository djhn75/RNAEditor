'''
Created on 05.06.2014

@author: david
'''

#!/usr/bin/env python
from VariantSet import VariantSet

v = VariantSet("/media/Storage/bio-data/Eva/rnaEditor/1024_Endothelzellen_LNA92a_R1.trimmed/1024_Endothelzellen_LNA92a_R1.trimmed.noSpliceJunction.vcf")

bamFile="/media/Storage/bio-data/Eva/rnaEditor/1024_Endothelzellen_LNA92a_R1.trimmed/1024_Endothelzellen_LNA92a_R1.trimmed.noDup.realigned.recalibrated.bam"
noEdges = v.bl