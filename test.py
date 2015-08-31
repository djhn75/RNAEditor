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
from dbCluster import dbCluster
from sklearn import metrics
from _collections import defaultdict


output="/media/Storage/bio-data/David/Kostas/rnaEditor/adar1/adar1"
#output="/media/ATLAS_NGS_storage/David/Kostas/rnaEditor/adar1/adar1"
#Helper.printResultHtml(output)
#Helper.createDiagramms(output)

"""aluVariants = VariantSet(output+".alu.vcf")
genome=Genome("/media/Storage/databases/rnaEditor_annotations/human/genes.gtf")
aluVariants.annotateVariantDict(genome)
aluVariants.printGeneList(genome,output+".testgvf",printSummary=True)
"""
#Helper.printResultHtml(output)


variants= VariantSet("/media/Storage/bio-data/David/Kostas/rnaEditor/adar1/adar1.editingSites.alu.vcf")
Yclust = dbCluster()

varPosListByChromosome = variants.getVarPosListByChromosome()

"""
for chromosome in varPosListByChromosome.keys():
    for eps in range(210,211):
        for min_samples in range(3,50):
            #print('EPS: %s, minSamples: %s' % (eps,min_samples))
            Yclust.dbscan(varPosListByChromosome[chromosome], eps=eps, min_samples=min_samples)
            core_sample_indices, labels = Yclust.coreSamples, Yclust.labels
            
            # Number of clusters in labels, ignoring noise if present.w
            n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
            labelNames = set(labels)
            X = []
            for el in varPosListByChromosome[chromosome]:
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
"""
eps=50,
min_samples=4
islandCounter=0
for chromosome in varPosListByChromosome.keys():
    Yclust.dbscan(varPosListByChromosome[chromosome], eps, min_samples)
    core_sample_indices, labels = Yclust.coreSamples, Yclust.labels
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    labelNames = set(labels)
    X = []
    clusterDict=defaultdict(list)
    for varPos,label in zip(varPosListByChromosome[chromosome],labels):
        clusterDict[label].append(varPos)
        X.append([varPos,0])
    X = np.array(X)

    if n_clusters_ > 0:
        #print('EPS: %s, minSamples: %s, #Clusters: %d, SC: %0.3f' % (eps, min_samples, n_clusters_, metrics.silhouette_score(X, labels)))
        for key in clusterDict.keys():
            if key==-1:
                continue
            
            biggest=max(clusterDict[key])
            smallest=min(clusterDict[key])
            length=biggest-smallest
            editingRate=float(len(clusterDict[key]))/float(length)
            print "\t".join(map(str,[chromosome,biggest,smallest,"editingIsland"+str(islandCounter),length,editingRate,len(clusterDict[key])]))
            islandCounter+=1
        #print('#Clusters: %d, SC: %0.3f' % (n_clusters_,metrics.silhouette_score(X, labels)))

#print "%d number of variants" % len(varPosListByChromosome)



core_sample_indices, labels = Yclust.coreSamples, Yclust.labels
# Number of clusters in labels, ignoring noise if present.w
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

"""
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

