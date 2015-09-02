'''
Created on 13.02.2015

@author: david
'''
import numpy as np
from random import shuffle
from sklearn.metrics.pairwise import pairwise_distances 
from VariantSet import VariantSet
from _collections import defaultdict


class DbCluster():
    def __init__(self,variants, eps=10, minSamples=5):
        self.eps=eps
        self.minSamples=minSamples
        return self.createClusters(variants)

