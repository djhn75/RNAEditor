'''
Created on 05.06.2014

@author: david
'''

#!/usr/bin/env python
# a bar plot with errorbars
import numpy as np
import matplotlib.pyplot as plt
import os
from collections import OrderedDict
from Helper import Helper
from VariantSet import VariantSet
from Genome import Genome


output="/media/ATLAS_NGS_storage/David/Kostas/rnaEditor/adar1/adar1"
Helper.printResultHtml(output)



      
"""def deleteNonEditingBases(variants):
    startTime=Helper.getTime()
    Helper.info("Delete non Editing Bases (keep only T->C and A->G)")
    
    for varTuple in variants.variantDict.keys():
        chr,pos,ref,alt = varTuple
        if (ref =="A" and alt == "G") or (ref=="T" and alt=="C"):
            pass
        else:
            del variants.variantDict[varTuple]


aluVariants = VariantSet("/media/ATLAS_NGS_storage/David/Kostas/rnaEditor/adar1/adar1.alu.vcf")
output="/media/ATLAS_NGS_storage/David/Kostas/rnaEditor/adar1/adar1"
genome=Genome("/media/Storage/databases/rnaEditor_annotations/human/genes.gtf")

#print Alu editing Sites
aluVariants.annotateVariantDict(genome)
#deleteNonEditingBases(aluVariants)
aluVariants.printVariantDict(output+".editingSites.alu.vcf")
aluVariants.printGeneList(genome,output+".editingSites.alu.gvf",printSummary=True)


Helper.createDiagramms(output)
"""

"""
variants= VariantSet("/media/Storage/bio-data/David/Kostas/scrambleN/scrambleN_1.vcf")
Yclust = dbCluster()

varPosList = []


for v in variants:
    varPosList.append(v.position)
varPosList = np.asarray(varPosList)

for eps in range(211,211):
    for min_samples in range(3,50):
        #print('EPS: %s, minSamples: %s' % (eps,min_samples))
        Yclust.dbscan(varPosList, eps=eps, min_samples=min_samples)
        core_sample_indices, labels = Yclust.coreSamples, Yclust.labels
        
        # Number of clusters in labels, ignoring noise if present.w
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
        labelNames = set(labels)
        X = []
        for el in varPosList:
            X.append([el,0])
        X = np.array(X)

        if n_clusters_ > 0:
            print('EPS: %s, minSamples: %s, #Clusters: %d, SC: %0.3f' % (eps, min_samples, n_clusters_, metrics.silhouette_score(X, labels)))
            #print('#Clusters: %d, SC: %0.3f' % (n_clusters_,metrics.silhouette_score(X, labels)))
        else:
            continue
            #print('EPS: %s, minSamples: %s, #Clusters: %d, SC: %0.3f' % (eps, min_samples, n_clusters_, 0))
            
        #print('Estimated number of clusters: %d' % n_clusters_)
        #print("Silhouette Coefficient: %0.3f" % metrics.silhouette_score(X, labels))
        #print("Silhouette Coefficient: %0.3f" % metrics.silhouette_score(X, labels))

Yclust.dbscan(varPosList, eps=2, min_samples=5)

print "%d number of variants" % len(varPosList)

core_sample_indices, labels = Yclust.coreSamples, Yclust.labels
# Number of clusters in labels, ignoring noise if present.w
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

labelNames = set(labels)
#labels = list(labels)


for label in labelNames:
    if label !=-1:
        #get indices of poinst which belong to the current dbCluster
        labelIndices = np.where(labels == label)[0]
        clusterPoints = varPosList[labelIndices]
        #print "%d" % label
        #print "%s" % ",".join(map(str,clusterPoints))
        mini=min(clusterPoints)
        maxi= max(clusterPoints)
        mean= std(clusterPoints)
        dens= float(len(clusterPoints))/ float((maxi-mini))
        print "dbCluster %d: [%s] ,min: %d, max: %d, mean: %d, Density: %f" % (int(label),",".join(map(str,clusterPoints)),mini, maxi, mean,dens)

X = []
for el in varPosList:
    X.append([el,0])
X = np.array(X)


print('Estimated number of clusters: %d' % n_clusters_)
print("Silhouette Coefficient: %0.3f" % metrics.silhouette_score(X, labels))


"""
