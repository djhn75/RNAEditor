'''
Created on 05.06.2014

@author: david
'''

from VariantSet import VariantSet
from dbCluster import dbCluster
import numpy as np
from sklearn import metrics
from numpy.ma.core import mean, std

variants= VariantSet("/media/Storage/bio-data/David/Kostas/scrambleN/scrambleN_1.vcf")
Yclust = dbCluster()

varPosList = []


for v in variants:
    varPosList.append(v.position)
varPosList = np.asarray(varPosList)

for eps in range(2,100):
    for min_samples in range(3,50):
        print('EPS: %s, minSamples: %s' % (eps,min_samples))
        Yclust.dbscan(varPosList, eps=eps, min_samples=min_samples)
        core_sample_indices, labels = Yclust.coreSamples, Yclust.labels
        
        # Number of clusters in labels, ignoring noise if present.w
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
        labelNames = set(labels)
        X = []
        for el in varPosList:
            X.append([el,0])
        X = np.array(X)


        #print('#Clusters: %d, SC: %0.3f' % (n_clusters_,metrics.silhouette_score(X, labels)))
        print('Estimated number of clusters: %d' % n_clusters_)
        print("Silhouette Coefficient: %0.3f" % metrics.silhouette_score(X, labels))
        #print("Silhouette Coefficient: %0.3f" % metrics.silhouette_score(X, labels))

"""Yclust.dbscan(varPosList, eps=20, min_samples=5)

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

ces = CallEditingSites(bamFile="/media/Storage/bio-data/David/Kostas/scrambleN/scrambleN.realigned.marked.recalibrated.bam",
                        refGenome="/media/Storage/databases/rnaEditor_annotations/human/human_g1k_v37.fasta", 
                        dbsnp="/media/Storage/databases/rnaEditor_annotations/human/dbsnp_135.b37_Y.vcf", 
                        hapmap="/media/Storage/databases/rnaEditor_annotations/human/hapmap_3.3.b37.sites.vcf",
                        omni="/media/Storage/databases/rnaEditor_annotations/human/1000G_omni2.5.b37.sites.vcf",
                        esp="/media/Storage/databases/rnaEditor_annotations/human/NHLBI_Exome_Sequencing_Project_6500SI.vcf", 
                        aluRegions="/media/Storage/databases/rnaEditor_annotations/human/Alu_repeats_noCHR.bed", 
                        gtfFile="/media/Storage/databases/rnaEditor_annotations/human/genes_Y.gtf", 
                        outfilePrefix=vcfFile[:vcfFile.rfind(".")], 
                        sourceDir="/usr/local/bin/")

g=Genome("/media/Storage/databases/rnaEditor_annotations/human/genes_Y.gtf")"""




#variants.deleteOverlappsFromVcf("/media/Storage/databases/rnaEditor_annotations/human/dbsnp_135.b37_Y.vcf")
"""variants.annotateVariantDict(g)
variants.printGeneList(g, "dink.gvf", True)
gbc=g.getGenesByChromosome()

for gene in gbc["Y"]:
    print gene.geneId"""
#print len(gbc["Y"])

#ces.removeHomopolymers(variants, "/media/Storage/bio-data/David/Kostas/scrambleN/scrambleN.nonAlu_Y.vcfs", 4)

#ces.start()

