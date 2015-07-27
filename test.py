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





N = 12
ind = np.arange(N)  # the x locations for the groups
width = 0.25       # the width of the bars
fig, ax = plt.subplots()

aluBaseCounts=Helper.getMMBaseCounts("/media/ATLAS_NGS_storage/David/Kostas/rnaEditor/adar1/adar1.alu.vcf")
aluMeans = aluBaseCounts.values()
rects1 = ax.bar(ind, aluMeans, width, color='r', )

nonAluBaseCounts=Helper.getMMBaseCounts("/media/ATLAS_NGS_storage/David/Kostas/rnaEditor/adar1/adar1.nonAlu.vcf")
nonAluMeans = nonAluBaseCounts.values()
rects2 = ax.bar(ind+width, nonAluMeans, width, color='b', )

# add some text for labels, title and axes ticks
ax.set_ylabel('Number')
ax.set_title('Variants per Base')
ax.set_xticks(ind+width)
ax.set_xticklabels( ("A->C","A->G","A->T","C->A","C->G","C->T","G->A","G->C","G->T","T->A","T->C","T->G") )

ax.legend( (rects1[0], rects2[0]), ('Alu', 'nonAlu') )

def autolabel(rects):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
                ha='center', va='bottom')

autolabel(rects1)
autolabel(rects2)

#fig.show()

fig.savefig("/media/Storage/bio-data/David/mmBarplot.pdf")



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
